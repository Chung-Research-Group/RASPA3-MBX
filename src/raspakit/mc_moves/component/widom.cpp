module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

module mc_moves_widom;

#ifdef USE_STD_IMPORT
import std;
#endif

import component;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import randomnumbers;
import system;
import energy_factor;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import running_energy;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_mbx;
import mc_moves_move_types;

double MC_Moves::WidomMove(RandomNumber& random, System& system, std::size_t selectedComponent)
{
  std::size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
  MoveTypes move = MoveTypes::Widom;
  Component& component = system.components[selectedComponent];
  std::chrono::system_clock::time_point t1, t2;

  // Update move statistics for Widom insertion move.
  component.mc_moves_statistics.addTrial(move);

  double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb = system.forceField.cutOffCoulomb;
  Component::GrowType growType = component.growType;

  // Attempt to grow a new molecule using Configurational Bias Monte Carlo (CBMC) insertion.
  t1 = std::chrono::system_clock::now();
  std::optional<ChainGrowData> growData = CBMC::growMoleculeSwapInsertion(
      random, component, system.hasExternalField, system.forceField, system.simulationBox, 
      system.interpolationGrids, system.externalFieldInterpolationGrid,
      system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), system.beta, growType,
      cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedMolecule, 1.0, false, false);
  t2 = std::chrono::system_clock::now();

  component.mc_moves_cputime[move]["NonEwald"] += (t2 - t1);
  system.mc_moves_cputime[move]["NonEwald"] += (t2 - t1);

  // If molecule growth failed, terminate the move.
  if (!growData) return 0.0;

  [[maybe_unused]] std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

  // Check if the new molecule is inside blocked pockets; if so, abort the move.
  if (system.insideBlockedPockets(component, newMolecule))
  {
    return 0.0;
  }

  // Update statistics for successfully constructed molecules.
  component.mc_moves_statistics.addConstructed(move);

  // Compute the energy difference in Ewald Fourier space due to the new molecule.
  t1 = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newMolecule, {});
  t2 = std::chrono::system_clock::now();

  component.mc_moves_cputime[move]["Ewald"] += (t2 - t1);
  system.mc_moves_cputime[move]["Ewald"] += (t2 - t1);

  // Compute the tail corrections for the energy due to the new molecule.
  t1 = std::chrono::system_clock::now();
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), newMolecule, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newMolecule, {});
  t2 = std::chrono::system_clock::now();

  component.mc_moves_cputime[move]["Tail"] += (t2 - t1);
  system.mc_moves_cputime[move]["Tail"] += (t2 - t1);

  // Compute the correction factor from Ewald and tail energy differences.
  double correctionFactorEwald =
      // std::exp(-system.beta * (energyFourierDifference.potentialEnergy()));
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy()));

  double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);

  // MBX Calculator
  if (system.useMBX)
  {
    // Compute the total energy difference from FF. Now it just calculates everything all again, no matter it has been calculated before or not. 
    // We can optimize this later by reusing the calculated energy difference from the CBMC growth and retrace steps. But it's not taking much time anyway, so we can leave it for now.
    // Compute tail energy difference due to long-range corrections
    RunningEnergy tailEnergyDifferenceInterMolecule = Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                system.spanOfMoleculeAtoms(), newMolecule, {});
    RunningEnergy tailEnergyDifferenceFrameworkMolecule = Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), newMolecule, {});
    RunningEnergy tailEnergyDifference = tailEnergyDifferenceInterMolecule + tailEnergyDifferenceFrameworkMolecule;

    std::optional<RunningEnergy> frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), newMolecule, {});
    std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), newMolecule, {});

    RunningEnergy energyDifferenceFF = frameworkMolecule.value() + interMolecule.value() + energyFourierDifference + tailEnergyDifference;  

    t1 = std::chrono::system_clock::now();
    // Energy of the system after the insertion of new trial molecule.
    // MBX will crash if the newly inserted atoms overlap the exisiting atoms. We have not added the check for that
    // as the check has already been placed in interMolecule and frameworkMolecule FF based calculation.
    RunningEnergy newTotalEnergy = Interactions::computeMBXEnergy(system, system.components, system.simulationBox, system.framework,
                                                       selectedComponent, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
                                                       newMolecule, true);
    t2 = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["MBX"] += (t2 - t1);
    system.mc_moves_cputime[move]["MBX"] += (t2 - t1);

    // Energy of the system before the insertion of trial molecule
    RunningEnergy oldTotalEnergy = system.runningEnergies;
    
    // MBX energy difference before and after the insertion move old and new configuration 
    RunningEnergy energyDifferenceMBX{};
    energyDifferenceMBX.mbxEnergy = newTotalEnergy.mbxEnergy - oldTotalEnergy.mbxEnergy;
    
    // The energyDifference for frameworkMoleculeVDW contribution as obtained from forceField.
    energyDifferenceMBX.frameworkMoleculeVDW = frameworkMolecule.frameworkMoleculeVDW;
    energyDifferenceMBX.tail = tailEnergyDifferenceFrameworkMolecule.tail;

    // Logging
    std::cerr << "MBX_E " 
              << energyDifferenceMBX.potentialEnergy() 
              << " "
              << "FF_E " 
              << energyDifferenceFF.potentialEnergy() 
              << " "
              << "XYZ_1 " 
              << newMolecule[0].position.x << " " << newMolecule[0].position.y << " " << newMolecule[0].position.z << " "
              << "XYZ_2 " 
              << newMolecule[1].position.x << " " << newMolecule[1].position.y << " " << newMolecule[1].position.z << " "
              << "XYZ_3 " 
              << newMolecule[2].position.x << " " << newMolecule[2].position.y << " " << newMolecule[2].position.z << " "
              << "\n";

    double correctionFactorMBX = std::exp(-system.beta * (energyDifferenceMBX.potentialEnergy() - energyDifferenceFF.potentialEnergy()))

    return correctionFactorEwald * correctionFactorMBX * growData->RosenbluthWeight / idealGasRosenbluthWeight;
  }
  
  return correctionFactorEwald * growData->RosenbluthWeight / idealGasRosenbluthWeight;
}
