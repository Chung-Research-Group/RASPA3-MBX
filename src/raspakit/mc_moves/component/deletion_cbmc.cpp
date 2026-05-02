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

module mc_moves_deletion_cbmc;

#ifdef USE_STD_IMPORT
import std;
#endif

import double3;
import double3x3;
import simd_quatd;
import component;
import atom;
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
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import interactions_mbx;
import units;
import mc_moves_move_types;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::deletionMoveCBMC(RandomNumber& random, System& system,
                                                                            std::size_t selectedComponent,
                                                                            std::size_t selectedMolecule)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::SwapCBMC;
  Component& component = system.components[selectedComponent];

  // Increment the count of swap deletion moves for the selected component
  component.mc_moves_statistics.addTrial(move, 1);

  // Proceed only if there is at least one molecule of the selected component
  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
  {
    // Get a reference to the molecule being deleted
    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);
    std::copy(system.electricField.begin(), system.electricField.end(), system.electricFieldNew.begin());
    // std::span<double3> electricFieldMoleculeNew = system.spanElectricFieldNew(selectedComponent, selectedMolecule);

    // Retrieve cutoff distances from the force field
    double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
    double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
    double cutOffCoulomb = system.forceField.cutOffCoulomb;
    Component::GrowType growType = component.growType;

    // Retrace the molecule for the swap deletion using CBMC algorithm
    time_begin = std::chrono::system_clock::now();
    ChainRetraceData retraceData = CBMC::retraceMoleculeSwapDeletion(
        random, system.components[selectedComponent], system.hasExternalField, system.forceField, system.simulationBox,
        system.interpolationGrids, system.externalFieldInterpolationGrid, system.framework,
        system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), system.beta, growType, cutOffFrameworkVDW,
        cutOffMoleculeVDW, cutOffCoulomb, molecule);
    time_end = std::chrono::system_clock::now();

    // Update the CPU time statistics for the non-Ewald part of the move
    component.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);

    // Compute the energy difference in Fourier space due to the deletion
    time_begin = std::chrono::system_clock::now();
    RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, {}, molecule);
    time_end = std::chrono::system_clock::now();
    // Update the CPU time statistics for the Ewald part of the move
    component.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

    // Compute the tail energy difference due to the deletion
    time_begin = std::chrono::system_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                system.spanOfMoleculeAtoms(), {}, molecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                   system.spanOfFrameworkAtoms(), {}, molecule);
    time_end = std::chrono::system_clock::now();
    // Update the CPU time statistics for the tail corrections
    component.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

    // Update the constructed count for the move statistics
    component.mc_moves_statistics.addConstructed(move, 1);

    RunningEnergy polarizationDifference;
    if (system.forceField.computePolarization)
    {
      std::vector<Atom> old_molecule = std::vector(molecule.begin(), molecule.end());
      std::vector<double3> old_electric_field = std::vector<double3>(old_molecule.size());

      Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                    system.spanOfFrameworkAtoms(), {},
                                                                    old_electric_field, {}, old_molecule);

      Interactions::computeEwaldFourierElectricFieldDifference(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
          system.totalEik, system.forceField, system.simulationBox, {}, old_electric_field, {}, old_molecule);

      // Compute polarization energy difference
      polarizationDifference = Interactions::computePolarizationEnergyDifference(system.forceField, {},
                                                                                 old_electric_field, {}, old_molecule);
    }

    // Calculate the correction factor for Ewald summation
    double correctionFactorEwald =
        std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy() +
                                 polarizationDifference.potentialEnergy()));

    // Compute acceptance probability factors
    double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (system.beta * fugacity * system.simulationBox.volume);

    double Pacc = preFactor * idealGasRosenbluthWeight / retraceData.RosenbluthWeight;
    std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
    double biasTransitionMatrix = system.getTMMCBiasFactor(selectedComponent, false);

    // Check if the new macrostate is within the allowed TMMC range
    if (system.doTMMC)
    {
      std::size_t newN = oldN - 1;
      std::pair<std::size_t, std::size_t> minmax = system.getTMMCMinMax(selectedComponent);
      if (newN < minmax.first)
      {
        return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
      }
    }

    // Energy of the system before the insertion of trial molecule
    RunningEnergy oldTotalEnergy = system.runningEnergies;

    RunningEnergy energyDifferenceFF =
        -retraceData.energies + energyFourierDifference + tailEnergyDifference + polarizationDifference;
    RunningEnergy energyDifferenceMBX;
    RunningEnergy newTotalEnergy = RunningEnergy();
    std::vector<double> mbxEnergyLog(7, 0.0);
    if (!system.useMBX)
    {
      // Energy logging
      // std::cerr << "deletion_cbmc" << "," << selectedComponent << ","
      //           << (system.numberOfIntegerMoleculesPerComponent[selectedComponent] - 1) << ","
      //           << (oldTotalEnergy.potentialEnergy() + energyDifferenceFF.potentialEnergy()) << ","
      //           << (oldTotalEnergy.frameworkMoleculeVDW + energyDifferenceFF.frameworkMoleculeVDW) << ","
      //           << (oldTotalEnergy.moleculeMoleculeVDW + energyDifferenceFF.moleculeMoleculeVDW) << ","
      //           << (oldTotalEnergy.tail + energyDifferenceFF.tail) << ","
      //           << (oldTotalEnergy.frameworkMoleculeCharge + energyDifferenceFF.frameworkMoleculeCharge) << ","
      //           << (oldTotalEnergy.moleculeMoleculeCharge + energyDifferenceFF.moleculeMoleculeCharge) << ","
      //           << ((oldTotalEnergy.ewald_fourier + energyDifferenceFF.ewald_fourier) +
      //               (oldTotalEnergy.ewald_self + energyDifferenceFF.ewald_self) +
      //               (oldTotalEnergy.ewald_exclusion + energyDifferenceFF.ewald_exclusion))
      //           << "," << energyDifferenceFF.potentialEnergy() << "," << Pacc << "\n";
    }
    else
    {
      // Compute the total energy difference from FF. Now it just calculates everything all again, no matter it has been
      // calculated before or not. We can optimize this later by reusing the calculated energy difference from the CBMC
      // growth and retrace steps. But it's not taking much time anyway, so we can leave it for now.

      // Now we calculate the MBX energy difference.
      // We calculate the system energy difference before and after the CMBC Reinsertion
      time_begin = std::chrono::system_clock::now();
      newTotalEnergy = Interactions::computeMBXEnergy(
          system, system.components, system.simulationBox, system.framework, selectedComponent,
          system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), molecule, false, &mbxEnergyLog);

      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["MBX"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["MBX"] += (time_end - time_begin);

      // MBX energy difference before and after the insertion move old and new configuration
      energyDifferenceMBX.mbxEnergy = newTotalEnergy.mbxEnergy - oldTotalEnergy.mbxEnergy;

      // The energyDifference for frameworkMoleculeVDW contribution as obtained from forceField
      energyDifferenceMBX.frameworkMoleculeVDW = -retraceData.energies.frameworkMoleculeVDW;

      // Compute tail energy difference due to long-range corrections
      RunningEnergy tailEnergyDifferenceFrameworkMolecule = Interactions::computeFrameworkMoleculeTailEnergyDifference(
          system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), {}, molecule);
      energyDifferenceMBX.tail = tailEnergyDifferenceFrameworkMolecule.tail;

      // Add the correction factor, exp(-beta*DeltaDeltaE)
      Pacc *= std::exp(-system.beta * (energyDifferenceMBX.potentialEnergy() - energyDifferenceFF.potentialEnergy()));
    }

    // Apply acceptance/rejection rule
    if (random.uniform() < biasTransitionMatrix * Pacc)
    {
      component.mc_moves_statistics.addAccepted(move, 1);

      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      system.deleteMolecule(selectedComponent, selectedMolecule, molecule);

      if (system.useMBX)
      {
        // Energy logging
        std::cerr << "deletion_cbmc" << "," << selectedComponent << ","
                << (system.numberOfIntegerMoleculesPerComponent[selectedComponent]) << ","
                << (oldTotalEnergy.potentialEnergy() + energyDifferenceMBX.potentialEnergy()) << ","
                << (oldTotalEnergy.frameworkMoleculeVDW + energyDifferenceMBX.frameworkMoleculeVDW) << ","
                << (oldTotalEnergy.tail + energyDifferenceMBX.tail) << "," << newTotalEnergy.mbxEnergy << ","
                << (mbxEnergyLog[1] /= Units::EnergyToKCalPerMol) << ","  // e2b
                << (mbxEnergyLog[2] /= Units::EnergyToKCalPerMol) << ","  // e3b
                << (mbxEnergyLog[3] /= Units::EnergyToKCalPerMol) << ","  // e4b
                << (mbxEnergyLog[4] /= Units::EnergyToKCalPerMol) << ","  // edisp
                << (mbxEnergyLog[5] /= Units::EnergyToKCalPerMol) << ","  // eelec_perm
                << (mbxEnergyLog[6] /= Units::EnergyToKCalPerMol) << ","  // eelec_ind
                << energyDifferenceMBX.potentialEnergy() << "," << Pacc << "\n";

        return {-energyDifferenceMBX, double3(Pacc, 1.0 - Pacc, 0.0)};
      }

      // return the -energyDifferenceFF here.
      return {-energyDifferenceFF, double3(Pacc, 1.0 - Pacc, 0.0)};
      // Below is the original return.
      // return {retraceData.energies - energyFourierDifference - tailEnergyDifference - polarizationDifference,
      //         double3(Pacc, 1.0 - Pacc, 0.0)};
    };
    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  // Return default values if no molecules are available for deletion
  return {std::nullopt, double3(0.0, 1.0, 0.0)};
}
