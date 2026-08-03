module;

module mc_moves_volume;

import std;

import component;
import atom;
import molecule;
import int3;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import energy_status;
import energy_status_inter;
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
#ifdef BUILD_MBX
import interactions_mbx;
#endif
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::volumeMove(RandomNumber &random, System &system)
{
  if (system.useMBX && system.framework.has_value())
  {
    throw std::runtime_error(
        "MBX isotropic volume moves are only supported for systems without a framework");
  }
  if (system.useMBX && system.hasExternalField)
  {
    throw std::runtime_error(
        "MBX isotropic volume moves are not supported with an external field");
  }

  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::VolumeChange;

  // Update volume move counts
  system.mc_moves_statistics.addTrial(move);

  RunningEnergy oldTotalEnergy = system.runningEnergies;
#ifdef BUILD_MBX
  if (system.useMBX)
  {
    time_begin = std::chrono::steady_clock::now();
    oldTotalEnergy = Interactions::computeMBXEnergySystem(
        system, system.components, system.simulationBox, system.framework, {}, system.spanOfMoleculeAtoms());
    time_end = std::chrono::steady_clock::now();
    system.mc_moves_cputime[move][Move::Timing::MBX] += (time_end - time_begin);
  }
#else
  if (system.useMBX)
  {
    throw std::runtime_error(
        "MBX volume energy was requested, but this RASPA binary was built without BUILD_MBX support");
  }
#endif
  // Calculate the total number of molecules
  double numberOfMolecules = static_cast<double>(std::accumulate(system.numberOfIntegerMoleculesPerComponent.begin(),
                                                                 system.numberOfIntegerMoleculesPerComponent.end(), 0));
  double oldVolume = system.simulationBox.volume;
  double maxVolumeChange = system.mc_moves_statistics.getMaxChange(move);

  // Propose a new volume change
  double newVolume = std::exp(std::log(oldVolume) + maxVolumeChange * (2.0 * random.uniform() - 1.0));

  // Compute scaling factor for box dimensions
  double scale = std::pow(newVolume / oldVolume, 1.0 / 3.0);

  SimulationBox newBox = system.simulationBox.scaled(scale);
  std::pair<std::vector<Molecule>, std::vector<Atom>> newPositions =
      system.scaledCenterOfMassPositions(system.simulationBox, newBox);

  double cutOffFrameworkVDW_stored = system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW_stored = system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb_stored = system.forceField.cutOffCoulomb;
  double ewald_alpha_stored = system.forceField.EwaldAlpha;
  int3 ewald_k_stored = system.forceField.numberOfWaveVectors;
  std::size_t reciprocalCutOffSquared_stored = system.forceField.reciprocalIntegerCutOffSquared;

  const auto restoreCutOffs = [&system, cutOffFrameworkVDW_stored, cutOffMoleculeVDW_stored,
                               cutOffCoulomb_stored, ewald_alpha_stored, ewald_k_stored,
                               reciprocalCutOffSquared_stored]()
  {
    system.forceField.cutOffFrameworkVDW = cutOffFrameworkVDW_stored;
    system.forceField.cutOffMoleculeVDW = cutOffMoleculeVDW_stored;
    system.forceField.cutOffCoulomb = cutOffCoulomb_stored;
    system.forceField.EwaldAlpha = ewald_alpha_stored;
    system.forceField.numberOfWaveVectors = ewald_k_stored;
    system.forceField.reciprocalIntegerCutOffSquared = reciprocalCutOffSquared_stored;
  };

  system.forceField.initializeAutomaticCutOff(newBox);

  time_begin = std::chrono::steady_clock::now();
  // Compute new intermolecular energy
  RunningEnergy newTotalInterEnergy =
      Interactions::computeInterMolecularEnergy(system.forceField, newBox, newPositions.second);
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  // Do not pass a severely compressed/overlapping trial to MBX. The total-energy kernel does not return
  // an overlap sentinel, so inspect its VDW contribution before evaluating the many-body model.
  if (system.useMBX &&
      (!std::isfinite(newTotalInterEnergy.moleculeMoleculeVDW) ||
       newTotalInterEnergy.moleculeMoleculeVDW > system.forceField.energyOverlapCriteria))
  {
    restoreCutOffs();
    return std::nullopt;
  }

  time_begin = std::chrono::steady_clock::now();
  // Compute new tail corrections. The tail energy is position-independent, so a volume change only rescales it
  // through the box volume; reuse the maintained effective pseudo-atom-type counts (O(nType^2)).
  RunningEnergy newTotalTailEnergy = Interactions::computeInterMolecularTailEnergyAggregated(
      system.forceField, newBox, system.effectiveNumberOfPseudoAtomsVDW, system.fractionalPseudoAtomCountsPerGroup);
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

  time_begin = std::chrono::steady_clock::now();
  // Compute new Ewald Fourier energy
  RunningEnergy newTotalEwaldEnergy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.trialEik,
      system.forceField, newBox, system.components, system.numberOfMoleculesPerComponent, newPositions.second,
      system.netChargeFramework);
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Sum up all energy contributions
  RunningEnergy newTotalEnergy = newTotalInterEnergy + newTotalTailEnergy + newTotalEwaldEnergy;

  // The intra-molecular energies have not changed by the com-scaling
  newTotalEnergy.bond = oldTotalEnergy.bond;
  newTotalEnergy.ureyBradley = oldTotalEnergy.ureyBradley;
  newTotalEnergy.bend = oldTotalEnergy.bend;
  newTotalEnergy.inversionBend = oldTotalEnergy.inversionBend;
  newTotalEnergy.outOfPlaneBend = oldTotalEnergy.outOfPlaneBend;
  newTotalEnergy.torsion = oldTotalEnergy.torsion;
  newTotalEnergy.improperTorsion = oldTotalEnergy.improperTorsion;
  newTotalEnergy.bondBond = oldTotalEnergy.bondBond;
  newTotalEnergy.bondBend = oldTotalEnergy.bondBend;
  newTotalEnergy.bondTorsion = oldTotalEnergy.bondTorsion;
  newTotalEnergy.bendBend = oldTotalEnergy.bendBend;
  newTotalEnergy.bendTorsion = oldTotalEnergy.bendTorsion;
  newTotalEnergy.intraVDW = oldTotalEnergy.intraVDW;
  newTotalEnergy.intraCoul = oldTotalEnergy.intraCoul;

  std::array<double, 7> mbxEnergyTerms{};
  if (system.useMBX)
  {
#ifdef BUILD_MBX
    time_begin = std::chrono::steady_clock::now();
    newTotalEnergy = Interactions::computeMBXEnergySystem(
        system, system.components, newBox, system.framework, {}, newPositions.second, mbxEnergyTerms);
    time_end = std::chrono::steady_clock::now();
    system.mc_moves_cputime[move][Move::Timing::MBX] += (time_end - time_begin);
#else
    throw std::runtime_error(
        "MBX volume energy was requested, but this RASPA binary was built without BUILD_MBX support");
#endif
  }

  // Update constructed move counts
  system.mc_moves_statistics.addConstructed(move);

  const double energyDifference = newTotalEnergy.potentialEnergy() - oldTotalEnergy.potentialEnergy();
  const double acceptanceProbability =
      std::exp((numberOfMolecules + 1.0) * std::log(newVolume / oldVolume) -
               (system.pressure * (newVolume - oldVolume) + energyDifference) * system.beta);

  // Apply acceptance/rejection rule
  if (random.uniform() < acceptanceProbability)
  {
    // Move accepted: update system state
    system.mc_moves_statistics.addAccepted(move);

    system.simulationBox = newBox;
    std::copy(newPositions.first.begin(), newPositions.first.end(), system.moleculeData.begin());
    // scaledCenterOfMassPositions returns molecule atoms only; framework atoms are untouched.
    std::copy(newPositions.second.begin(), newPositions.second.end(), system.spanOfMoleculeAtoms().begin());

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    system.writeAcceptedEnergyLog("volume", std::numeric_limits<std::size_t>::max(), newTotalEnergy,
                                  mbxEnergyTerms, energyDifference, acceptanceProbability);

    return newTotalEnergy - oldTotalEnergy;
  }

  restoreCutOffs();

  return std::nullopt;
}

std::optional<RunningEnergy> MC_Moves::anisotropicVolumeMove(RandomNumber& random, System& system)
{
  if (system.useMBX)
  {
    throw std::runtime_error("Anisotropic volume moves are not supported with MBX");
  }

  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::AnisotropicVolumeChange;

  system.mc_moves_statistics.addTrial(move);

  RunningEnergy oldTotalEnergy = system.runningEnergies;
  const double numberOfMolecules = static_cast<double>(std::accumulate(
      system.numberOfIntegerMoleculesPerComponent.begin(), system.numberOfIntegerMoleculesPerComponent.end(), 0));
  const SimulationBox& oldBox = system.simulationBox;
  const double oldVolume = oldBox.volume;

  const double3 maxVolumeChange(
      system.mc_moves_statistics.getMaxChange(move, 0),
      system.mc_moves_statistics.getMaxChange(move, 1),
      system.mc_moves_statistics.getMaxChange(move, 2));

  const double3 scale(std::exp(maxVolumeChange.x * (2.0 * random.uniform() - 1.0)),
                      std::exp(maxVolumeChange.y * (2.0 * random.uniform() - 1.0)),
                      std::exp(maxVolumeChange.z * (2.0 * random.uniform() - 1.0)));

  SimulationBox newBox = oldBox.scaled(scale);

  // Scale molecules by their center of mass: only the center-of-mass fractional position is kept fixed while
  // each molecule is translated rigidly. This preserves the internal geometry of both rigid and flexible
  // molecules, so the intramolecular energies remain unchanged (see below).
  std::pair<std::vector<Molecule>, std::vector<Atom>> newPositions =
      system.scaledCenterOfMassPositions(oldBox, newBox);

  const double cutOffFrameworkVDW_stored = system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW_stored = system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb_stored = system.forceField.cutOffCoulomb;
  const double ewald_alpha_stored = system.forceField.EwaldAlpha;
  const int3 ewald_k_stored = system.forceField.numberOfWaveVectors;
  const std::size_t reciprocalCutOffSquared_stored = system.forceField.reciprocalIntegerCutOffSquared;

  system.forceField.initializeAutomaticCutOff(newBox);

  time_begin = std::chrono::steady_clock::now();
  RunningEnergy newTotalInterEnergy =
      Interactions::computeInterMolecularEnergy(system.forceField, newBox, newPositions.second);
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  time_begin = std::chrono::steady_clock::now();
  // Tail energy is position-independent; reuse the maintained effective pseudo-atom-type counts (O(nType^2)).
  RunningEnergy newTotalTailEnergy = Interactions::computeInterMolecularTailEnergyAggregated(
      system.forceField, newBox, system.effectiveNumberOfPseudoAtomsVDW, system.fractionalPseudoAtomCountsPerGroup);
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

  time_begin = std::chrono::steady_clock::now();
  RunningEnergy newTotalEwaldEnergy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.trialEik,
      system.forceField, newBox, system.components, system.numberOfMoleculesPerComponent, newPositions.second,
      system.netChargeFramework);
  time_end = std::chrono::steady_clock::now();
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  RunningEnergy newTotalEnergy = newTotalInterEnergy + newTotalTailEnergy + newTotalEwaldEnergy;

  newTotalEnergy.bond = oldTotalEnergy.bond;
  newTotalEnergy.ureyBradley = oldTotalEnergy.ureyBradley;
  newTotalEnergy.bend = oldTotalEnergy.bend;
  newTotalEnergy.inversionBend = oldTotalEnergy.inversionBend;
  newTotalEnergy.outOfPlaneBend = oldTotalEnergy.outOfPlaneBend;
  newTotalEnergy.torsion = oldTotalEnergy.torsion;
  newTotalEnergy.improperTorsion = oldTotalEnergy.improperTorsion;
  newTotalEnergy.bondBond = oldTotalEnergy.bondBond;
  newTotalEnergy.bondBend = oldTotalEnergy.bondBend;
  newTotalEnergy.bondTorsion = oldTotalEnergy.bondTorsion;
  newTotalEnergy.bendBend = oldTotalEnergy.bendBend;
  newTotalEnergy.bendTorsion = oldTotalEnergy.bendTorsion;
  newTotalEnergy.intraVDW = oldTotalEnergy.intraVDW;
  newTotalEnergy.intraCoul = oldTotalEnergy.intraCoul;

  system.mc_moves_statistics.addConstructed(move);

  const double volumeRatio = scale.x * scale.y * scale.z;
  const double pressureWork =
      oldVolume * (system.pressureTensorDiagonal.x * (scale.x - 1.0) +
                   system.pressureTensorDiagonal.y * (scale.y - 1.0) +
                   system.pressureTensorDiagonal.z * (scale.z - 1.0));

  const double energyDifference = newTotalEnergy.potentialEnergy() - oldTotalEnergy.potentialEnergy();
  const double acceptanceProbability =
      std::exp((numberOfMolecules + 1.0) * std::log(volumeRatio) -
               (pressureWork + energyDifference) * system.beta);

  if (random.uniform() < acceptanceProbability)
  {
    system.mc_moves_statistics.addAccepted(move);

    system.simulationBox = newBox;
    std::copy(newPositions.first.begin(), newPositions.first.end(), system.moleculeData.begin());
    // scaledCenterOfMassPositions returns molecule atoms only; framework atoms are untouched.
    std::copy(newPositions.second.begin(), newPositions.second.end(), system.spanOfMoleculeAtoms().begin());

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    system.writeAcceptedEnergyLog("anisotropic-volume", std::numeric_limits<std::size_t>::max(), newTotalEnergy, {},
                                  energyDifference, acceptanceProbability);

    return newTotalEnergy - oldTotalEnergy;
  }

  system.forceField.cutOffFrameworkVDW = cutOffFrameworkVDW_stored;
  system.forceField.cutOffMoleculeVDW = cutOffMoleculeVDW_stored;
  system.forceField.cutOffCoulomb = cutOffCoulomb_stored;
  system.forceField.EwaldAlpha = ewald_alpha_stored;
  system.forceField.numberOfWaveVectors = ewald_k_stored;
  system.forceField.reciprocalIntegerCutOffSquared = reciprocalCutOffSquared_stored;

  return std::nullopt;
}
