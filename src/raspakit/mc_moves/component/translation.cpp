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

module mc_moves_translation;

#ifdef USE_STD_IMPORT
import std;
#endif

import component;
import molecule;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import energy_factor;
import energy_status;
import energy_status_inter;
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import interactions_mbx;
import units;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::translationMove(RandomNumber &random, System &system,
                                                       std::size_t selectedComponent, std::size_t selectedMolecule,
                                                       const std::vector<Component> &components, Molecule &molecule,
                                                       std::span<Atom> molecule_atoms)
{
  double3 displacement{};
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::Translation;
  Component &component = system.components[selectedComponent];

  // Randomly select a direction (0 for x, 1 for y, 2 for z)
  std::size_t selectedDirection = std::size_t(3.0 * random.uniform());

  // Get the maximum displacement allowed for the selected component
  double maxDisplacement = component.mc_moves_statistics.getMaxChange(move, selectedDirection);

  // Compute the displacement along the selected direction
  displacement[selectedDirection] = maxDisplacement * 2.0 * (random.uniform() - 0.5);

  // Update move statistics for the selected direction
  component.mc_moves_statistics.addTrial(move, selectedDirection);

  // Construct the trial molecule with the new displacement
  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      components[selectedComponent].translate(molecule, molecule_atoms, displacement);

  // Check if the trial molecule is inside blocked pockets
  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    return std::nullopt;
  }

  // Update move construction statistics
  component.mc_moves_statistics.addConstructed(move, selectedDirection);
  
  std::vector<double3> electricFieldMoleculeNew(molecule_atoms.size());
  std::vector<double3> electricFieldMoleculeOld(molecule_atoms.size());

  // Energy of the system before the insertion of trial molecule
  RunningEnergy oldTotalEnergy = system.runningEnergies;
  RunningEnergy energyDifference;

  if (system.useMBX)
  {
    std::vector<double> mbxEnergyLog(7, 0);         // Vector to store energylog values

    time_begin = std::chrono::system_clock::now();
    RunningEnergy newTotalEnergy = Interactions::computeMBXEnergy(system, components, system.simulationBox, system.framework,
                                                       selectedComponent, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
                                                       trialMolecule.second, true, &mbxEnergyLog);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["MBX"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["MBX"] += (time_end - time_begin);
    
    // MBX energy difference before and after the insertion move old and new configuration 
    energyDifference.mbxEnergy = newTotalEnergy.mbxEnergy - oldTotalEnergy.mbxEnergy;

    // Compute framework-molecule energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkMolecule;

    frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), trialMolecule.second, molecule_atoms);
    
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Framework-Molecule"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Framework-Molecule"] += (time_end - time_begin);
    if (!frameworkMolecule.has_value()) return std::nullopt;

    energyDifference.frameworkMoleculeVDW = frameworkMolecule.value().frameworkMoleculeVDW; 
    
    // compute framework-molecule tail energy contribution
    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifferenceFrameworkMolecule = Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfFrameworkAtoms(), trialMolecule.second,
                                                              molecule_atoms);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

    energyDifference.tail = tailEnergyDifferenceFrameworkMolecule.tail;
    
    // Here you can add logging commands
    std::cerr << "translation" << ","
              << selectedComponent << ","
              << system.numberOfIntegerMoleculesPerComponent[selectedComponent] << ","
              << (oldTotalEnergy.potentialEnergy() + energyDifference.potentialEnergy()) << ","
              << (oldTotalEnergy.frameworkMoleculeVDW + energyDifference.frameworkMoleculeVDW) << ","
              << (oldTotalEnergy.tail + energyDifference.tail) << ","
              << newTotalEnergy.mbxEnergy << ","
              << (mbxEnergyLog[1] /= Units::EnergyToKCalPerMol) << ","  // e2b
              << (mbxEnergyLog[2] /= Units::EnergyToKCalPerMol) << ","  // e3b
              << (mbxEnergyLog[3] /= Units::EnergyToKCalPerMol) << ","  // e4b
              << (mbxEnergyLog[4] /= Units::EnergyToKCalPerMol) << ","  // edisp
              << (mbxEnergyLog[5] /= Units::EnergyToKCalPerMol) << ","  // eelec_perm
              << (mbxEnergyLog[6] /= Units::EnergyToKCalPerMol) << ","  // eelec_ind
              << energyDifference.potentialEnergy() << ","
              << (std::exp(-system.beta * energyDifference.potentialEnergy())) << "\n";

  }
  else
  {

    // Compute external field energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, 
        system.externalFieldInterpolationGrid, trialMolecule.second, molecule_atoms);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["ExternalField-Molecule"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["ExternalField-Molecule"] += (time_end - time_begin);
    if (!externalFieldMolecule.has_value()) return std::nullopt;

    // Compute framework-molecule energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> frameworkMolecule;
    if (system.forceField.computePolarization)
    {
      frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second,
          molecule_atoms);
    }
    else
    {
      frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), trialMolecule.second, molecule_atoms);
    }
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Framework-Molecule"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Framework-Molecule"] += (time_end - time_begin);
    if (!frameworkMolecule.has_value()) return std::nullopt;

    // Compute molecule-molecule energy contribution
    time_begin = std::chrono::system_clock::now();
    std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms);
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Molecule-Molecule"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Molecule-Molecule"] += (time_end - time_begin);
    if (!interMolecule.has_value()) return std::nullopt;

    // Compute Ewald energy contribution
    time_begin = std::chrono::system_clock::now();
    RunningEnergy ewaldFourierEnergy;
    if (system.forceField.computePolarization)
    {
      ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
          system.totalEik, system.forceField, system.simulationBox, electricFieldMoleculeNew, electricFieldMoleculeOld,
          trialMolecule.second, molecule_atoms);
    }
    else
    {
      ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
          system.simulationBox, trialMolecule.second, molecule_atoms);
    }
    time_end = std::chrono::system_clock::now();
    component.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

    RunningEnergy polarizationDifference;
    if (system.forceField.computePolarization)
    {
      // Compute polarization energy difference
      polarizationDifference = Interactions::computePolarizationEnergyDifference(
          system.forceField, electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second, molecule_atoms);
    }

    // Calculate the total energy difference
    energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                    ewaldFourierEnergy + polarizationDifference; 

    // Here you can add logging commands
    std::cerr << "translation" << ","
              << selectedComponent << ","
              << system.numberOfIntegerMoleculesPerComponent[selectedComponent] << ","
              << (oldTotalEnergy.potentialEnergy() + energyDifference.potentialEnergy()) << ","
              << (oldTotalEnergy.frameworkMoleculeVDW + energyDifference.frameworkMoleculeVDW) << ","
              << (oldTotalEnergy.moleculeMoleculeVDW + energyDifference.moleculeMoleculeVDW) << ","
              << (oldTotalEnergy.tail + energyDifference.tail) << ","
              << (oldTotalEnergy.frameworkMoleculeCharge + energyDifference.frameworkMoleculeCharge) << ","
              << (oldTotalEnergy.moleculeMoleculeCharge + energyDifference.moleculeMoleculeCharge) << ","
              << ((oldTotalEnergy.ewald_fourier + energyDifference.ewald_fourier) +
                 (oldTotalEnergy.ewald_self + energyDifference.ewald_self) +
                 (oldTotalEnergy.ewald_exclusion + energyDifference.ewald_exclusion)) << ","
              << energyDifference.potentialEnergy() << ","
              << (std::exp(-system.beta * energyDifference.potentialEnergy())) << "\n";
  }
  
  // Apply acceptance/rejection rule based on Metropolis criterion
  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    // Update acceptance statistics
    component.mc_moves_statistics.addAccepted(move, selectedDirection);

    // Accept the Ewald move
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

    // Update the molecule atoms with the trial positions
    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), molecule_atoms.begin());
    molecule = trialMolecule.first;

    // Update the electric field if polarization is computed
    if (system.forceField.computePolarization && !system.useMBX)
    {
      std::span<double3> electricFieldMolecule = system.spanElectricFieldNew(selectedComponent, selectedMolecule);
      std::copy(electricFieldMoleculeNew.begin(), electricFieldMoleculeNew.end(), electricFieldMolecule.begin());
    }

    return energyDifference;
  };
  return std::nullopt;
}
