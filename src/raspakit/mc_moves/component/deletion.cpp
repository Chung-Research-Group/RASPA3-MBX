module;

module mc_moves_deletion;

import std;

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

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::deletionMove(RandomNumber& random, System& system,
                                                                        std::size_t selectedComponent,
                                                                        std::size_t selectedMolecule)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::Swap;
  Component& component = system.components[selectedComponent];

  // Increment swap deletion move counts for the selected component
  component.mc_moves_statistics.addTrial(move, 1);

  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
  {
    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

    // Increment constructed swap deletion move counts
    component.mc_moves_statistics.addConstructed(move, 1);

    double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double preFactor = double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (system.beta * fugacity * system.simulationBox.volume);

    RunningEnergy oldTotalEnergy = system.runningEnergies;
    RunningEnergy energyDifference;

    if (system.useMBX)
    {
      time_begin = std::chrono::system_clock::now();
      // Energy of the system after the insertion of new trial molecule.
      std::vector<double> mbxEnergyLog(7, 0);  // Vector to store energylog values
      // MBX will crash if the newly inserted atoms overlap the exisiting atoms. We have not added the check for that
      // as the check has already been placed in interMolecule and frameworkMolecule FF based calculation.
      RunningEnergy newTotalEnergy = Interactions::computeMBXEnergy(
          system, system.components, system.simulationBox, system.framework, selectedComponent,
          system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), molecule, false, &mbxEnergyLog);
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["MBX"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["MBX"] += (time_end - time_begin);

      energyDifference.mbxEnergy = newTotalEnergy.mbxEnergy - oldTotalEnergy.mbxEnergy;

      // Compute framework-molecule energy contribution
      time_begin = std::chrono::system_clock::now();

      std::optional<RunningEnergy> frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), {}, molecule);
      if (!frameworkMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Framework-Molecule"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Framework-Molecule"] += (time_end - time_begin);

      energyDifference.frameworkMoleculeVDW = frameworkMolecule.value().frameworkMoleculeVDW;

      // compute framework-molecule tail energy contribution
      time_begin = std::chrono::system_clock::now();
      std::optional<RunningEnergy> tailEnergyDifferenceFrameworkMolecule =
          Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                     system.spanOfFrameworkAtoms(), {}, molecule);
      if (!frameworkMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};
      time_end = std::chrono::system_clock::now();
      component.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

      energyDifference.tail = tailEnergyDifferenceFrameworkMolecule.value().tail;

      // Energy logging
#if DEBUG
      std::cerr << "deletion" << "," << selectedComponent << ","
                << (system.numberOfIntegerMoleculesPerComponent[selectedComponent] - 1) << ","
                << (oldTotalEnergy.potentialEnergy() + energyDifference.potentialEnergy()) << ","
                << (oldTotalEnergy.frameworkMoleculeVDW + energyDifference.frameworkMoleculeVDW) << ","
                << (oldTotalEnergy.tail + energyDifference.tail) << "," << newTotalEnergy.mbxEnergy << ","
                << (mbxEnergyLog[1] /= Units::EnergyToKCalPerMol) << ","  // e2b
                << (mbxEnergyLog[2] /= Units::EnergyToKCalPerMol) << ","  // e3b
                << (mbxEnergyLog[3] /= Units::EnergyToKCalPerMol) << ","  // e4b
                << (mbxEnergyLog[4] /= Units::EnergyToKCalPerMol) << ","  // edisp
                << (mbxEnergyLog[5] /= Units::EnergyToKCalPerMol) << ","  // eelec_perm
                << (mbxEnergyLog[6] /= Units::EnergyToKCalPerMol) << ","  // eelec_ind
                << energyDifference.potentialEnergy() << ","
                << preFactor * (std::exp(-system.beta * energyDifference.potentialEnergy())) << "\n";
#endif
    }
    else
    {
      // Copy the current electric field if polarization is computed
      std::vector<double3> electricFieldMoleculeOld(molecule.size());

      // Compute external field energy contribution
      std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
          system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid, {},
          molecule);
      if (!externalFieldMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

      // Compute framework-molecule energy contribution
      std::optional<RunningEnergy> frameworkMolecule;
      if (system.forceField.computePolarization)
      {
        frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
            system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
            system.spanOfFrameworkAtoms(), {}, electricFieldMoleculeOld, {}, molecule);
      }
      else
      {
        frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
            system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
            system.spanOfFrameworkAtoms(), {}, molecule);
      }
      if (!frameworkMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

      // Compute molecule-molecule energy contribution
      std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
          system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), {}, molecule);
      if (!interMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

      // Compute Ewald Fourier energy difference
      time_begin = std::chrono::system_clock::now();
      RunningEnergy energyFourierDifference;
      if (system.forceField.computePolarization)
      {
        energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
            system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
            system.totalEik, system.forceField, system.simulationBox, {}, electricFieldMoleculeOld, {}, molecule);
      }
      else
      {
        energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
            system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik,
            system.forceField, system.simulationBox, {}, molecule);
      }
      time_end = std::chrono::system_clock::now();

      // Update CPU time statistics for Ewald calculations
      component.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

      // Compute tail correction energy difference
      time_begin = std::chrono::system_clock::now();
      [[maybe_unused]] RunningEnergy tailEnergyDifference =
          Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfMoleculeAtoms(), {}, molecule) +
          Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                     system.spanOfFrameworkAtoms(), {}, molecule);
      time_end = std::chrono::system_clock::now();

      // Update CPU time statistics for tail corrections
      component.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);
      system.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

      RunningEnergy polarizationDifference;
      if (system.forceField.computePolarization)
      {
        // Compute polarization energy difference
        polarizationDifference = Interactions::computePolarizationEnergyDifference(
            system.forceField, {}, electricFieldMoleculeOld, {}, molecule);
      }

      // Get the total difference in energy
      energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                         energyFourierDifference + tailEnergyDifference + polarizationDifference;

      // Energy logging
#if DEBUG
      std::cerr << "deletion" << "," << selectedComponent << ","
                << (system.numberOfIntegerMoleculesPerComponent[selectedComponent] - 1) << ","
                << (oldTotalEnergy.potentialEnergy() + energyDifference.potentialEnergy()) << ","
                << (oldTotalEnergy.frameworkMoleculeVDW + energyDifference.frameworkMoleculeVDW) << ","
                << (oldTotalEnergy.moleculeMoleculeVDW + energyDifference.moleculeMoleculeVDW) << ","
                << (oldTotalEnergy.tail + energyDifference.tail) << ","
                << (oldTotalEnergy.frameworkMoleculeCharge + energyDifference.frameworkMoleculeCharge) << ","
                << (oldTotalEnergy.moleculeMoleculeCharge + energyDifference.moleculeMoleculeCharge) << ","
                << ((oldTotalEnergy.ewald_fourier + energyDifference.ewald_fourier) +
                    (oldTotalEnergy.ewald_self + energyDifference.ewald_self) +
                    (oldTotalEnergy.ewald_exclusion + energyDifference.ewald_exclusion))
                << "," << energyDifference.potentialEnergy() << ","
                << (preFactor * (std::exp(-system.beta * energyDifference.potentialEnergy()))) << "\n";
#endif
    }

    // Calculate the acceptance probability
    double Pacc = preFactor * std::exp(-system.beta * energyDifference.potentialEnergy());
    std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
    double biasTransitionMatrix = system.getTMMCBiasFactor(selectedComponent, false);

    // Check if TMMC is enabled and if new state is below minimum macrostate
    if (system.doTMMC)
    {
      std::size_t newN = oldN - 1;
      std::pair<std::size_t, std::size_t> minmax = system.getTMMCMinMax(selectedComponent);

      if (newN < minmax.first)
      {
        return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
      }
    }

    // Apply acceptance/rejection rule
    if (random.uniform() < biasTransitionMatrix * Pacc)
    {
      // Move accepted; update acceptance statistics
      component.mc_moves_statistics.addAccepted(move, 1);

      // Accept Ewald move and delete molecule from system
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      system.deleteMolecule(selectedComponent, selectedMolecule, molecule);

      return {-energyDifference, double3(Pacc, 1.0 - Pacc, 0.0)};
    };
    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  // No molecules to delete; return default values
  return {std::nullopt, double3(0.0, 1.0, 0.0)};
}
