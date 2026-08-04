module;

module system;

import std;

import archive;
import randomnumbers;
import stringutils;
import int3;
import uint3;
import double3;
import double3x3;
import double3x3x3;
import simd_quatd;
import cubic;
import atom;
import framework;
import component;
import cbmc_move_statistics;
import simulationbox;
import forcefield;
import units;
import property_loading;
import averages;
import skparser;
import skposcarparser;
import skstructure;
import skatom;
import skcell;
import sample_movies;
import property_enthalpy;
import property_pressure;
import energy_dudlambda;
import energy_status;
import energy_status_inter;
import energy_status_intra;
import property_simulationbox;
import average_energy_type;
import property_energy;
import property_partial_molar_properties;
import property_lambda_probability_histogram;
import property_widom;
import property_temperature;
import property_msd;
import running_energy;
import threadpool;
// import isotherm;
// import multi_site_isotherm;
// import pressure_range;
import bond_potential;
import move_statistics;
import mc_moves_probabilities;
import mc_moves_move_types;
import mc_moves_cputime;
import reaction;
import reactions;
import cbmc;
import cbmc_chain_data;
import interactions_framework_molecule;
import interactions_framework_molecule_grid;
import interactions_intermolecular;
import interactions_pair_kernel;
import interactions_ewald;
import interactions_internal;
import interactions_external_field;
import interactions_external_field_grid;
#ifdef BUILD_MBX
import interactions_mbx;
#endif
import interactions_polarization_derivatives;
import equation_of_states;
import thermostat;
import thermobarostat;
import json;
import integrators;
import integrators_compute;
import integrators_update;
import interpolation_energy_grid;
import property_number_of_molecules_evolution;
import property_volume_evolution;
import property_conserved_energy_evolution;
import minimization_cell_layout;
#if !(defined(__has_include) && __has_include(<mdspan>))
// import mdspan;
#endif

// construct System programmatically
/*! \brief Brief description.
 *         Brief description continued.
 *
 *  Detailed description starts here.
 */
System::System(ForceField forcefield, std::optional<SimulationBox> box, bool hasExternalField, double T,
               std::optional<double> P, double heliumVoidFraction, std::optional<Framework> f, std::vector<Component> c,
               std::vector<std::vector<double3>> initialpositions, std::vector<std::size_t> initialNumberOfMolecules,
               std::size_t numberOfBlocks, const MCMoveProbabilities& systemProbabilities, bool useMBXCalculator,
               std::optional<std::string> mbxFilePath, bool writeEnergyLog,
               bool computeZeroLoadingHeatOfAdsorption)
    : temperature(T),
      pressure(P.value_or(0.0) / Units::PressureConversionFactor),
      input_pressure(P.value_or(0.0)),
      beta(1.0 / (Units::KB * T)),
      heliumVoidFraction(heliumVoidFraction),
      framework(f),
      components(c),
      loadings(c.size()),
      swappableComponents(),
      initialNumberOfMolecules(initialNumberOfMolecules),
      numberOfMoleculesPerComponent(c.size()),
      numberOfIntegerMoleculesPerComponent(c.size()),
      numberOfFractionalMoleculesPerComponent(c.size()),
      numberOfGCFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfPairGCFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfPairSwapFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfGibbsFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfParallelReactionFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC(c.size()),
      numberOfReactionFractionalMoleculesPerComponent_CFCMC(),
      idealGasEnergiesPerComponent(c.size()),
      forceField(forcefield),
      hasExternalField(hasExternalField),
      useMBX(useMBXCalculator),
      mbxSettingsFilePath(mbxFilePath.value_or(std::string{})),
      writeEnergyLog(writeEnergyLog),
      computeZeroLoadingHeatOfAdsorption(computeZeroLoadingHeatOfAdsorption),
      numberOfPseudoAtoms(c.size(), std::vector<std::size_t>(forceField.pseudoAtoms.size())),
      totalNumberOfPseudoAtoms(forceField.pseudoAtoms.size()),
      atomData({}),
      moleculeData({}),
      runningEnergies(),
      currentEnergyStatus(1, f.has_value() ? 1 : 0, c.size()),
      netChargePerComponent(c.size()),
      mc_moves_probabilities(systemProbabilities),
      mc_moves_statistics(),
      reactions(),
      tmmc(),
      averageEnergies(numberOfBlocks, 1, f.has_value() ? 1 : 0, c.size()),
      averageLoadings(numberOfBlocks, c.size()),
      averageEnthalpiesOfAdsorption(numberOfBlocks, c.size()),
      averagePartialMolarProperties(numberOfBlocks, c.size()),
      averageTemperature(numberOfBlocks),
      averageTranslationalTemperature(numberOfBlocks),
      averageRotationalTemperature(numberOfBlocks),
      averagePressure(numberOfBlocks),
      averageSimulationBox(numberOfBlocks),
      interpolationGrids(forceField.pseudoAtoms.size() + 1, std::nullopt)
{
  if (useMBX)
  {
#ifndef BUILD_MBX
    throw std::runtime_error(
        "[System]: UseMBX was requested, but this RASPA binary was built without BUILD_MBX support");
#else
    if (!mbxFilePath.has_value() || mbxFilePath->empty())
    {
      throw std::runtime_error("[System]: UseMBX requires MBXSettingsFile");
    }

    const std::filesystem::path settingsPath(*mbxFilePath);
    std::error_code error;
    if (!std::filesystem::is_regular_file(settingsPath, error))
    {
      throw std::runtime_error(std::format("[System]: MBX settings file '{}' is not a readable regular file{}",
                                           settingsPath.string(),
                                           error ? std::format(" ({})", error.message()) : std::string{}));
    }
    std::ifstream settings(settingsPath, std::ios::binary);
    if (!settings || settings.peek() == std::ifstream::traits_type::eof())
    {
      throw std::runtime_error(
          std::format("[System]: MBX settings file '{}' is unreadable or empty", settingsPath.string()));
    }

    const std::filesystem::path canonicalPath = std::filesystem::canonical(settingsPath, error);
    mbxSettingsFilePath =
        error ? std::filesystem::absolute(settingsPath).lexically_normal().string() : canonicalPath.string();
#endif
  }
  currentEnergyStatus.useMBX = useMBX;

  input_pressureTensorDiagonal = double3(input_pressure, input_pressure, input_pressure);
  pressureTensorDiagonal = double3(pressure, pressure, pressure);

  // Temperature-dependent potentials (Feynman-Hibbs) require the external temperature;
  // recompute the derived constants, shifts, and tail-corrections with the system temperature.
  if (forceField.temperature != T)
  {
    forceField.temperature = T;
    forceField.preComputeDerivedParameters();
    forceField.preComputePotentialShift();
    forceField.preComputeTailCorrection();
  }

  if (box.has_value())
  {
    simulationBox = box.value();
  }

  removeRedundantMoves();
  determineSwappableComponents();
  determineFractionalComponents();
  assignDUdlambdaGroups();
  rescaleMoveProbabilities();
  rescaleMolarFractions();
  computeNumberOfPseudoAtoms();
  computeTailCorrectionCounts();

  createFrameworks();
  if (framework.has_value())
  {
    simulationBox = framework->simulationBox.scaled(framework->numberOfUnitCells);
  }

  if (useMBX) preComputeElecPermFrameworkMBX();

  forceField.initializeEwaldParameters(simulationBox);

  CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      eik_x, eik_y, eik_z, eik_xy, forceField, simulationBox, double3(0.0, 0.0, 0.0), 1.0);

  precomputeTotalRigidEnergy();

  translationalCenterOfMassConstraint = 0;
  translationalDegreesOfFreedom = 0;
  rotationalDegreesOfFreedom = 0;
  if (framework && framework->hasMobileAtoms())
  {
    if (framework->isMixed())
    {
      translationalDegreesOfFreedom += 3 * framework->flexibleAtomCount + 3 * framework->numberOfRigidGroups();
      for (const FrameworkGroup& group : framework->groups)
      {
        if (group.isRigidBody()) rotationalDegreesOfFreedom += group.rotationalDegreesOfFreedom;
      }
    }
    else
    {
      translationalDegreesOfFreedom += 3 * numberOfFrameworkAtoms;
    }
  }

  createInitialMolecules(initialpositions);
  initializeFixedLambdaFractionalMolecules();
  computeTailCorrectionCounts();

  if (computeZeroLoadingHeatOfAdsorption)
  {
    if (numberOfMolecules() != 0)
    {
      throw std::runtime_error(
          "[System]: ComputeZeroLoadingHeatOfAdsorption requires zero initial guest molecules");
    }
    if (framework.has_value() && framework->hasMobileAtoms())
    {
      throw std::runtime_error(
          "[System]: ComputeZeroLoadingHeatOfAdsorption currently requires a rigid framework");
    }
    if (mc_moves_probabilities.getProbability(Move::Types::VolumeChange) > 0.0 ||
        mc_moves_probabilities.getProbability(Move::Types::AnisotropicVolumeChange) > 0.0)
    {
      throw std::runtime_error(
          "[System]: ComputeZeroLoadingHeatOfAdsorption currently requires a fixed simulation volume");
    }
    bool hasWidomComponent = false;
    for (const Component& component : components)
    {
      if (component.mc_moves_probabilities.getProbability(Move::Types::Widom) <= 0.0) continue;
      hasWidomComponent = true;
      if (!component.rigid)
      {
        throw std::runtime_error(
            "[System]: ComputeZeroLoadingHeatOfAdsorption currently requires rigid Widom components");
      }
    }
    if (!hasWidomComponent)
    {
      throw std::runtime_error(
          "[System]: ComputeZeroLoadingHeatOfAdsorption requires a positive WidomProbability");
    }
  }

  // Build the per-component ideal-gas conformation reservoirs used to seed CBMC growth. Done after the
  // initial molecules are placed so their placement keeps using the cold-start seed (unchanged initial
  // geometry); the reservoir is only consulted by the production Monte-Carlo moves.
  buildConformationReservoirs();

  equationOfState = EquationOfState(EquationOfState::Type::PengRobinson, EquationOfState::MixingRules::VanDerWaals, T,
                                    P.value_or(0.0), simulationBox, heliumVoidFraction, components);

  averageEnthalpiesOfAdsorption.resize(swappableComponents.size());
  averagePartialMolarProperties.resize(swappableComponents.size());
}

void System::createFrameworks()
{
  netChargeFramework = 0.0;
  if (framework.has_value())
  {
    const std::vector<Atom>& atoms = framework->atoms;
    for (const Atom& atom : atoms)
    {
      atomData.push_back(atom);
      atomDynamics.push_back(AtomDynamics{});
      electricPotential.push_back(0.0);
      electricField.push_back(double3(0.0, 0.0, 0.0));
      electricFieldNew.push_back(double3(0.0, 0.0, 0.0));
    }
    numberOfFrameworkAtoms += atoms.size();
    // Lab-fixed prefix for Ewald: all atoms when fully rigid, Fixed groups when mixed, else none.
    numberOfRigidFrameworkAtoms += framework->numberOfFixedAtoms();
    netChargeFramework += framework->netCharge;
    netCharge += framework->netCharge;
  }
}

void System::rebuildForFramework(const Framework& newFramework, const SimulationBox& newSimulationBox)
{
  const std::size_t previousTranslationalFrameworkDof =
      (framework.has_value() && framework->hasMobileAtoms())
          ? (framework->isMixed() ? 3 * framework->flexibleAtomCount + 3 * framework->numberOfRigidGroups()
                                  : 3 * numberOfFrameworkAtoms)
          : 0;
  std::size_t previousRotationalFrameworkDof = 0;
  if (framework.has_value() && framework->isMixed())
  {
    for (const FrameworkGroup& group : framework->groups)
    {
      if (group.isRigidBody()) previousRotationalFrameworkDof += group.rotationalDegreesOfFreedom;
    }
  }
  const std::size_t previousFrameworkAtoms = numberOfFrameworkAtoms;

  // Detach any guest-molecule storage (the suffix after the framework prefix) so it can be re-appended once
  // the new framework atoms are in place. For a framework-only system these are empty.
  const auto moleculeOffset = static_cast<std::vector<Atom>::difference_type>(previousFrameworkAtoms);
  std::vector<Atom> moleculeAtoms(atomData.begin() + moleculeOffset, atomData.end());
  std::vector<AtomDynamics> moleculeDynamics(atomDynamics.begin() + moleculeOffset, atomDynamics.end());

  framework = newFramework;
  simulationBox = newSimulationBox;

  const std::vector<Atom>& frameworkAtoms = framework->atoms;
  numberOfFrameworkAtoms = frameworkAtoms.size();
  numberOfRigidFrameworkAtoms = framework->numberOfFixedAtoms();

  // Rebuild the per-atom storage as [new framework atoms] ++ [preserved molecule atoms]; the field/potential
  // buffers are reset to the new size (they are recomputed on the next energy/field evaluation).
  atomData.assign(frameworkAtoms.begin(), frameworkAtoms.end());
  atomData.insert(atomData.end(), moleculeAtoms.begin(), moleculeAtoms.end());
  atomDynamics.assign(frameworkAtoms.size(), AtomDynamics{});
  atomDynamics.insert(atomDynamics.end(), moleculeDynamics.begin(), moleculeDynamics.end());
  electricPotential.assign(atomData.size(), 0.0);
  electricField.assign(atomData.size(), double3(0.0, 0.0, 0.0));
  electricFieldNew.assign(atomData.size(), double3(0.0, 0.0, 0.0));

  // Swap the old framework net-charge contribution for the new one (adsorbate contribution is unchanged).
  netCharge -= netChargeFramework;
  netChargeFramework = framework->netCharge;
  netCharge += netChargeFramework;

  // Adjust framework degrees of freedom for the replacement host.
  translationalDegreesOfFreedom -= previousTranslationalFrameworkDof;
  rotationalDegreesOfFreedom -= previousRotationalFrameworkDof;
  if (framework->hasMobileAtoms())
  {
    if (framework->isMixed())
    {
      translationalDegreesOfFreedom += 3 * framework->flexibleAtomCount + 3 * framework->numberOfRigidGroups();
      for (const FrameworkGroup& group : framework->groups)
      {
        if (group.isRigidBody()) rotationalDegreesOfFreedom += group.rotationalDegreesOfFreedom;
      }
    }
    else
    {
      translationalDegreesOfFreedom += 3 * numberOfFrameworkAtoms;
    }
  }

  // Note: pseudo-atom counts intentionally track only component (guest) atoms, matching the constructor which
  // computes them before the framework atoms are appended; the preserved guest counts stay valid, so framework
  // atoms are not (re)counted here.

  forceField.initializeEwaldParameters(simulationBox);
  CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      eik_x, eik_y, eik_z, eik_xy, forceField, simulationBox, double3(0.0, 0.0, 0.0), 1.0);

  precomputeTotalRigidEnergy();
  if (useMBX) preComputeElecPermFrameworkMBX();
}

void System::determineSwappableComponents()
{
  // Satellite components of a group-swap move are swappable even when they carry no swap
  // probability themselves: the group move inserts/deletes their molecules.
  for (const Component& component : components)
  {
    if (component.mc_moves_probabilities.getProbability(Move::Types::GroupSwap) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GroupSwapCBMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GroupSwapCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GroupSwapCBCFCMC) > 0.0)
    {
      for (std::size_t satelliteComponentId : component.groupComponentIds)
      {
        if (satelliteComponentId < components.size())
        {
          components[satelliteComponentId].swappable = true;
        }
      }
    }
  }

  for (std::size_t componentId{0}; Component& component : components)
  {
    if (component.mc_moves_probabilities.getProbability(Move::Types::Swap) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::SwapCBMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::PairSwapCBMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::PairSwap) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::PairSwapCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::PairSwapCBCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GroupSwap) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GroupSwapCBMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GroupSwapCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GroupSwapCBCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::SwapCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::SwapCBCFCMC) > 0.0)
    {
      component.swappable = true;
    }

    if (component.mc_moves_probabilities.getProbability(Move::Types::GibbsSwapCBMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GibbsSwapCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GibbsSwapCBCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GibbsConventionalCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GibbsConventionalCBCFCMC) > 0.0)
    {
      component.swappable = true;
    }

    if (component.swappable)
    {
      swappableComponents.push_back(componentId);
    }

    ++componentId;
  }
}

// determine the required number of fractional molecules

void System::rescaleMoveProbabilities()
{
  for (Component& component : components)
  {
    component.mc_moves_probabilities.join(mc_moves_probabilities);
  }
}

void System::removeRedundantMoves()
{
  for (Component& component : components)
  {
    component.mc_moves_probabilities.removeRedundantMoves();
  }
}

void System::optimizeMCMoves()
{
  mc_moves_statistics.optimizeMCMoves();
  for (Component& component : components)
  {
    component.mc_moves_statistics.optimizeMCMoves();

    // Adapt the internal CBMC / ring-closure Monte-Carlo step sizes (per bead) towards their target
    // acceptance ratios.
    for (CBMCMoveStatistics& cbmcStatistics : component.cbmc_moves_statistics)
    {
      cbmcStatistics.optimize();
    }
  }
}

void System::rescaleMolarFractions()
{
  double totalMolfraction = 0.0;
  double numberOfSwappableComponents = 0.0;
  for (const Component& component : components)
  {
    if (component.swappable)
    {
      totalMolfraction += component.molFraction;
      numberOfSwappableComponents += 1.0;
    }
  }

  if (totalMolfraction > 0.0)
  {
    for (Component& component : components)
    {
      if (component.swappable)
      {
        component.molFraction /= totalMolfraction;
      }
    }
  }
  else
  {
    for (Component& component : components)
    {
      if (component.swappable)
      {
        component.molFraction /= numberOfSwappableComponents;
      }
    }
  }
}

void System::computeNumberOfPseudoAtoms()
{
  for (std::size_t i = 0; i != components.size(); ++i)
  {
    std::fill(numberOfPseudoAtoms[i].begin(), numberOfPseudoAtoms[i].end(), 0);
  }
  std::fill(totalNumberOfPseudoAtoms.begin(), totalNumberOfPseudoAtoms.end(), 0);

  for (const Atom& atom : atomData)
  {
    std::size_t componentId = static_cast<std::size_t>(atom.componentId);
    std::size_t type = static_cast<std::size_t>(atom.type);
    numberOfPseudoAtoms[componentId][type] += 1;
    totalNumberOfPseudoAtoms[type] += 1;
  }
}

void System::computeTailCorrectionCounts()
{
  std::size_t numberOfPseudoAtomTypes = forceField.pseudoAtoms.size();

  effectiveNumberOfPseudoAtomsVDW.assign(numberOfPseudoAtomTypes, 0.0);
  for (std::size_t group = 0; group < maximumNumberOfDUDlambdaGroups; ++group)
  {
    fractionalPseudoAtomCountsPerGroup[group].assign(numberOfPseudoAtomTypes, 0.0);
  }

  for (const Atom& atom : spanOfMoleculeAtoms())
  {
    std::size_t type = static_cast<std::size_t>(atom.type);
    effectiveNumberOfPseudoAtomsVDW[type] += atom.scalingVDW;
    if (atom.groupId != 0)
    {
      fractionalPseudoAtomCountsPerGroup[static_cast<std::size_t>(atom.groupId) - 1][type] += 1.0;
    }
  }
}

void System::addAtomToTailCorrectionCounts(const Atom& atom)
{
  std::size_t type = static_cast<std::size_t>(atom.type);
  effectiveNumberOfPseudoAtomsVDW[type] += atom.scalingVDW;
  if (atom.groupId != 0)
  {
    fractionalPseudoAtomCountsPerGroup[static_cast<std::size_t>(atom.groupId) - 1][type] += 1.0;
  }
}

void System::removeAtomFromTailCorrectionCounts(const Atom& atom)
{
  std::size_t type = static_cast<std::size_t>(atom.type);
  effectiveNumberOfPseudoAtomsVDW[type] -= atom.scalingVDW;
  if (atom.groupId != 0)
  {
    fractionalPseudoAtomCountsPerGroup[static_cast<std::size_t>(atom.groupId) - 1][type] -= 1.0;
  }
}

void System::sampleProperties(std::size_t systemId, std::size_t currentBlock, std::size_t currentCycle)
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  double w = weight();

  averageSimulationBox.addSample(currentBlock, simulationBox, w);

  double translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(
      moleculeData, spanOfMoleculeAtoms(), spanOfMoleculeDynamics(), components, framework, spanOfFrameworkAtoms(),
      spanOfFrameworkDynamics(), &forceField, spanOfGroupData());
  double translationalTemperature =
      2.0 * translationalKineticEnergy /
      (Units::KB * static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint));
  averageTranslationalTemperature.addSample(currentBlock, translationalTemperature, w);

  double rotationalKineticEnergy =
      Integrators::computeRotationalKineticEnergy(moleculeData, components, spanOfGroupData());
  double rotationalTemperature =
      rotationalDegreesOfFreedom > 0
          ? 2.0 * rotationalKineticEnergy / (Units::KB * static_cast<double>(rotationalDegreesOfFreedom))
          : 0.0;
  averageRotationalTemperature.addSample(currentBlock, rotationalTemperature, w);

  double overallTemperature =
      2.0 * (translationalKineticEnergy + rotationalKineticEnergy) /
      (Units::KB * static_cast<double>(translationalDegreesOfFreedom - translationalCenterOfMassConstraint +
                                       rotationalDegreesOfFreedom));
  averageTemperature.addSample(currentBlock, overallTemperature, w);

  loadings = LoadingData(components.size(), numberOfIntegerMoleculesPerComponent, simulationBox);
  averageLoadings.addSample(currentBlock, loadings, w);

  EnthalpyOfAdsorptionTerms enthalpyTerms = EnthalpyOfAdsorptionTerms(
      swappableComponents, numberOfIntegerMoleculesPerComponent, runningEnergies.potentialEnergy(), temperature);
  averageEnthalpiesOfAdsorption.addSample(currentBlock, enthalpyTerms, w);

  PartialMolarPropertiesTerms partialMolarTerms =
      PartialMolarPropertiesTerms(swappableComponents, numberOfIntegerMoleculesPerComponent,
                                  runningEnergies.potentialEnergy(), simulationBox.volume);
  averagePartialMolarProperties.addSample(currentBlock, partialMolarTerms, w);

  std::size_t numberOfMolecules =
      std::accumulate(numberOfIntegerMoleculesPerComponent.begin(), numberOfIntegerMoleculesPerComponent.end(), 0uz);
  double currentIdealPressure = static_cast<double>(numberOfMolecules) / (beta * simulationBox.volume);

  averagePressure.addSample(currentBlock, currentIdealPressure, currentExcessPressureTensor, w);

  for (std::size_t componentId{0}; Component& component : components)
  {
    double componentDensity =
        static_cast<double>(numberOfIntegerMoleculesPerComponent[componentId]) / simulationBox.volume;

    double lambda = component.lambdaGC.lambdaValue();
    double dudlambda = currentDUdlambda(lambda, component.lambdaGC.dUdlambdaGroupId);
    component.lambdaGC.sampleHistogram(currentBlock, componentDensity, dudlambda, containsTheFractionalMolecule, w);

    if (usesGibbsConventionalCFCMC())
    {
      const double gibbsLambda = component.lambdaGibbs.lambdaValue();
      const double gibbsDudlambda = currentDUdlambda(gibbsLambda, component.lambdaGibbs.dUdlambdaGroupId);
      component.lambdaGibbs.sampleHistogram(currentBlock, componentDensity, gibbsDudlambda, true, w);
    }

    if (componentDrivesPairSwapLambda(componentId, Move::Types::PairSwapCFCMC))
    {
      const double pairLambda = component.lambdaPairSwap.lambdaValue();
      component.lambdaPairSwap.sampleHistogram(currentBlock, componentDensity,
                                               currentDUdlambda(pairLambda, component.lambdaPairSwap.dUdlambdaGroupId),
                                               containsTheFractionalMolecule, w);
    }

    if (componentDrivesPairSwapLambda(componentId, Move::Types::PairSwapCBCFCMC))
    {
      const double pairLambda = component.lambdaPairSwapCB.lambdaValue();
      component.lambdaPairSwapCB.sampleHistogram(
          currentBlock, componentDensity, currentDUdlambda(pairLambda, component.lambdaPairSwapCB.dUdlambdaGroupId),
          containsTheFractionalMolecule, w);
    }

    if (componentDrivesGroupSwapLambda(componentId, Move::Types::GroupSwapCFCMC))
    {
      const double groupLambda = component.lambdaGroupSwap.lambdaValue();
      component.lambdaGroupSwap.sampleHistogram(
          currentBlock, componentDensity, currentDUdlambda(groupLambda, component.lambdaGroupSwap.dUdlambdaGroupId),
          containsTheFractionalMolecule, w);
    }

    if (componentDrivesGroupSwapLambda(componentId, Move::Types::GroupSwapCBCFCMC))
    {
      const double groupLambda = component.lambdaGroupSwapCB.lambdaValue();
      component.lambdaGroupSwapCB.sampleHistogram(
          currentBlock, componentDensity, currentDUdlambda(groupLambda, component.lambdaGroupSwapCB.dUdlambdaGroupId),
          containsTheFractionalMolecule, w);
    }

    ++componentId;
  }

  reactionLambdaSampleProductionHistograms(currentBlock, w);

  updateSamplePDBMovie(systemId, currentCycle);

  if (writeLammpsData.has_value())
  {
    writeLammpsData->update(currentCycle, components, atomData, atomDynamics, moleculeData, simulationBox, forceField,
                            numberOfIntegerMoleculesPerComponent, framework);
  }

  if (propertyConventionalRadialDistributionFunction.has_value())
  {
    propertyConventionalRadialDistributionFunction->sample(simulationBox, spanOfFrameworkAtoms(), spanOfMoleculeAtoms(),
                                                           currentCycle, currentBlock);
  }

  // Force-based RDF is not sampled here. Monte Carlo calls sampleForceBasedRDFWithFullGradients()
  // (full U, including intramolecular). Molecular dynamics reuses integrator forces via
  // sampleForceBasedRDFFromCurrentGradients().

  if (propertyMoleculeProperties.has_value())
  {
    propertyMoleculeProperties->sample(components, numberOfMoleculesPerComponent, spanOfMoleculeAtoms(), currentCycle,
                                       currentBlock);
  }

  if (averageEnergyHistogram.has_value())
  {
    averageEnergyHistogram->addSample(
        currentBlock, currentCycle,
        {runningEnergies.potentialEnergy(), runningEnergies.frameworkMoleculeVDW + runningEnergies.moleculeMoleculeVDW,
         runningEnergies.frameworkMoleculeCharge + runningEnergies.moleculeMoleculeCharge +
             runningEnergies.ewald_fourier + runningEnergies.ewald_self + runningEnergies.ewald_exclusion,
         runningEnergies.polarization},
        w);
  }

  if (averageNumberOfMoleculesHistogram.has_value())
  {
    averageNumberOfMoleculesHistogram->addSample(currentBlock, currentCycle, numberOfIntegerMoleculesPerComponent, w);
  }

  if (propertyMSD.has_value())
  {
    propertyMSD->addSample(currentCycle, moleculeData);
  }

  if (propertyVACF.has_value())
  {
    propertyVACF->addSample(currentCycle, moleculeData);
  }

  if (propertyDensityGrid.has_value())
  {
    propertyDensityGrid->sample(framework, simulationBox, spanOfMoleculeAtoms(), currentCycle);
  }

  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

  mc_moves_cputime.propertySampling += (t2 - t1);
}

void System::samplePropertiesEvolution(std::size_t absoluteCurrentCycle)
{
  if (propertyNumberOfMoleculesEvolution.has_value())
  {
    propertyNumberOfMoleculesEvolution->addSample(absoluteCurrentCycle, numberOfIntegerMoleculesPerComponent);
  }
  if (propertyVolumeEvolution.has_value())
  {
    propertyVolumeEvolution->addSample(absoluteCurrentCycle, simulationBox.volume);
  }
  if (propertyConservedEnergyEvolution.has_value())
  {
    propertyConservedEnergyEvolution->addSample(absoluteCurrentCycle, runningEnergies);
  }
}

void System::precomputeTotalRigidEnergy() noexcept
{
  Interactions::precomputeEwaldFourierRigid(eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, forceField,
                                            simulationBox, spanOfRigidFrameworkAtoms());
}

void System::preComputeElecPermFrameworkMBX()
{
#ifdef BUILD_MBX
  elecPermFrameworkMBX = useMBX && framework.has_value() ? Interactions::computeFrameworkElecPermMBXEnergy(
                                                               *this, simulationBox, framework, spanOfFrameworkAtoms())
                                                         : 0.0;
#else
  if (useMBX)
  {
    throw std::runtime_error(
        "[System]: MBX state cannot be initialized because this RASPA binary was built without BUILD_MBX support");
  }
  elecPermFrameworkMBX = 0.0;
#endif
}

void System::validateMBXMonteCarloConfiguration() const
{
  if (!useMBX) return;

  if (tmmc.doTMMC)
  {
    throw std::runtime_error(
        "[System MBX validation]: transition-matrix Monte Carlo has not been forward-ported to MBX; "
        "use the regular MonteCarlo driver without TMMC state");
  }
  if (forceField.computePolarization)
  {
    throw std::runtime_error(
        "[System MBX validation]: 'ComputePolarization' must be false with UseMBX; MBX already provides the "
        "retained induced-electrostatic energy");
  }
  if (!reactions.list.empty())
  {
    throw std::runtime_error("[System MBX validation]: reaction moves have not been forward-ported to MBX");
  }
  if (propertyRadialDistributionFunction.has_value())
  {
    throw std::runtime_error(
        "[System MBX validation]: force-based RDF sampling is not supported with UseMBX because its force path "
        "uses the classical gradient energy; conventional geometry-only RDF sampling remains supported");
  }

  for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
  {
    const Component& component = components[componentId];
    if (!component.rigid && component.mc_moves_probabilities.getProbability(Move::Types::Swap) > 0.0)
    {
      throw std::runtime_error(std::format(
          "[System MBX validation, component {} ('{}')]: conventional insertion/deletion is supported only for "
          "rigid components; use SwapCBMC for flexible or semi-flexible components",
          componentId, component.name));
    }
  }

  const auto moveIsSupported = [this](Move::Types move)
  {
    switch (move)
    {
      case Move::Types::Translation:
      case Move::Types::Rotation:
      case Move::Types::ReinsertionCBMC:
      case Move::Types::Swap:
      case Move::Types::SwapCBMC:
      case Move::Types::Widom:
        return true;
      case Move::Types::VolumeChange:
        return !framework.has_value() && !hasExternalField;
      default:
        return false;
    }
  };

  for (std::size_t moveIndex = 0; moveIndex < std::to_underlying(Move::Types::Count); ++moveIndex)
  {
    const Move::Types move = static_cast<Move::Types>(moveIndex);
    double probability = mc_moves_probabilities.getProbability(move);
    for (const Component& component : components)
    {
      probability = std::max(probability, component.mc_moves_probabilities.getProbability(move));
    }
    if (probability <= 0.0 || moveIsSupported(move)) continue;

    if (move == Move::Types::VolumeChange && (framework.has_value() || hasExternalField))
    {
      throw std::runtime_error(
          "[System MBX validation]: isotropic volume moves are supported only for framework-free boxes without "
          "an external-field interpolation grid");
    }
    throw std::runtime_error(std::format(
        "[System MBX validation]: move '{}' has nonzero probability but has not been forward-ported to MBX. "
        "Supported moves are translation, rotation, reinsertion CBMC, conventional insertion/deletion for rigid "
        "components, CBMC insertion/deletion, Widom, and isotropic volume changes in framework-free boxes without "
        "an external field",
        Move::moveNames[moveIndex]));
  }
}

void System::precomputeTotalGradients() noexcept
{
  runningEnergies = Integrators::updateGradients(
      moleculeData, spanOfMoleculeAtoms(), spanOfMoleculeDynamics(), spanOfFrameworkAtoms(), forceField, simulationBox,
      components, eik_x, eik_y, eik_z, eik_xy, trialEik, fixedFrameworkStoredEik, interpolationGrids,
      numberOfMoleculesPerComponent, framework, spanOfFrameworkDynamics());
}

RunningEnergy System::computeTotalEnergies(std::span<double> mbxEnergyTerms)
{
  std::span<const Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<Atom> moleculeAtomPositions = spanOfMoleculeAtoms();

  if (useMBX)
  {
#ifndef BUILD_MBX
    throw std::runtime_error(
        "[System]: MBX energy was requested, but this RASPA binary was built without BUILD_MBX support");
#else
    // MBX owns every guest intra/intermolecular and electrostatic contribution. RASPA contributes only
    // framework-guest VDW (including its tail correction) and the configured external field. In particular,
    // no classical guest intra/intermolecular, Ewald, or polarization energy is added here.
    const RunningEnergy frameworkMoleculeEnergy = Interactions::computeFrameworkMoleculeEnergy(
        forceField, simulationBox, interpolationGrids, framework, frameworkAtomPositions, moleculeAtomPositions);
    const RunningEnergy frameworkMoleculeTailEnergy = Interactions::computeFrameworkMoleculeTailEnergy(
        forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions);

    RunningEnergy externalFieldEnergy;
    Interactions::computeExternalFieldEnergy(hasExternalField, forceField, simulationBox, moleculeAtomPositions,
                                             externalFieldEnergy, externalFieldInterpolationGrid);

    RunningEnergy result = Interactions::computeMBXEnergySystem(*this, components, simulationBox, framework,
                                                                frameworkAtomPositions, moleculeAtomPositions,
                                                                mbxEnergyTerms);
    result.frameworkMoleculeVDW = frameworkMoleculeEnergy.frameworkMoleculeVDW;
    result.tail = frameworkMoleculeTailEnergy.tail;
    result.externalFieldVDW = externalFieldEnergy.externalFieldVDW;
    result.externalFieldCharge = externalFieldEnergy.externalFieldCharge;
    return result;
#endif
  }

  if (!mbxEnergyTerms.empty())
  {
    throw std::invalid_argument("MBX energy subterms were requested for a system that does not use MBX");
  }

  if (fixedFrameworkStoredEik.empty())
  {
    precomputeTotalRigidEnergy();
  }

  RunningEnergy runningIntraEnergy{};
  std::size_t index = 0;
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (numberOfMoleculesPerComponent[i] > 0)
    {
      std::span<const Molecule> span_molecules = {&moleculeData[index], numberOfMoleculesPerComponent[i]};
      runningIntraEnergy += Interactions::computeIntraMolecularEnergy(components[i].intraMolecularPotentials,
                                                                      span_molecules, spanOfMoleculeAtoms());
    }

    index += numberOfMoleculesPerComponent[i];
  }
  if (framework && framework->hasMobileAtoms())
  {
    runningIntraEnergy += Interactions::computeFrameworkIntraMolecularEnergy(forceField, *framework, simulationBox,
                                                                             frameworkAtomPositions);
  }

  if (forceField.computePolarization)
  {
    std::span<double3> moleculeElectricField = spanOfMoleculeElectricField();

    std::fill(moleculeElectricField.begin(), moleculeElectricField.end(), double3(0.0, 0.0, 0.0));

    RunningEnergy frameworkMoleculeEnergy = Interactions::computeFrameworkMoleculeElectricField(
        forceField, simulationBox, moleculeElectricField, frameworkAtomPositions, moleculeAtomPositions);

    // When molecule-molecule polarization is requested, the inter-molecular electric field has to be added
    // to the stored field so that the polarization energy is consistent with the incremental Monte-Carlo moves
    // (which maintain the same field). Otherwise only the framework field contributes and the plain
    // inter-molecular energy is sufficient.
    RunningEnergy intermolecularEnergy;
    if (forceField.omitInterPolarization)
    {
      intermolecularEnergy =
          Interactions::computeInterMolecularEnergy(forceField, simulationBox, moleculeAtomPositions);
    }
    else
    {
      intermolecularEnergy = Interactions::computeInterMolecularElectricField(
          forceField, simulationBox, moleculeElectricField, moleculeAtomPositions);
    }

    RunningEnergy frameworkMoleculeTailEnergy = Interactions::computeFrameworkMoleculeTailEnergy(
        forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions);
    RunningEnergy intermolecularTailEnergy =
        Interactions::computeInterMolecularTailEnergy(forceField, simulationBox, moleculeAtomPositions);

    RunningEnergy ewaldEnergy = Interactions::computeEwaldFourierElectricField(
        eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, forceField, simulationBox,
        moleculeElectricField, components, numberOfMoleculesPerComponent, moleculeAtomPositions);

    RunningEnergy polarizationEnergy = computePolarizationEnergy();

    RunningEnergy externalFieldEnergy;
    Interactions::computeExternalFieldEnergy(hasExternalField, forceField, simulationBox, moleculeAtomPositions,
                                             externalFieldEnergy, externalFieldInterpolationGrid);

    return frameworkMoleculeEnergy + intermolecularEnergy + frameworkMoleculeTailEnergy + intermolecularTailEnergy +
           ewaldEnergy + polarizationEnergy + runningIntraEnergy + externalFieldEnergy;
  }
  else
  {
    RunningEnergy frameworkMoleculeEnergy = Interactions::computeFrameworkMoleculeEnergy(
        forceField, simulationBox, interpolationGrids, framework, frameworkAtomPositions, moleculeAtomPositions);
    RunningEnergy intermolecularEnergy =
        Interactions::computeInterMolecularEnergy(forceField, simulationBox, moleculeAtomPositions);

    RunningEnergy frameworkMoleculeTailEnergy = Interactions::computeFrameworkMoleculeTailEnergy(
        forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions);
    RunningEnergy intermolecularTailEnergy =
        Interactions::computeInterMolecularTailEnergy(forceField, simulationBox, moleculeAtomPositions);

    RunningEnergy ewaldEnergy = Interactions::computeEwaldFourierEnergy(
        eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, forceField, simulationBox, components,
        numberOfMoleculesPerComponent, moleculeAtomPositions, netChargeFramework);

    RunningEnergy externalFieldEnergy;
    Interactions::computeExternalFieldEnergy(hasExternalField, forceField, simulationBox, moleculeAtomPositions,
                                             externalFieldEnergy, externalFieldInterpolationGrid);

    return frameworkMoleculeEnergy + intermolecularEnergy + frameworkMoleculeTailEnergy + intermolecularTailEnergy +
           ewaldEnergy + runningIntraEnergy + externalFieldEnergy;
  }
}

RunningEnergy System::computePolarizationEnergy() noexcept
{
  RunningEnergy energy{};

  std::span<const Atom> moleculeAtomPositions = spanOfMoleculeAtoms();
  std::span<double3> moleculeElectricField = spanOfMoleculeElectricField();

  for (std::size_t i = 0; i < moleculeAtomPositions.size(); ++i)
  {
    std::size_t type = moleculeAtomPositions[i].type;
    // The polarization coupling is scaled by the atom's Coulomb scaling so that a fractional (CFCMC)
    // molecule decouples from the field as lambda decreases (matching the incremental moves).
    double polarizability = moleculeAtomPositions[i].scalingCoulomb * forceField.pseudoAtoms[type].polarizability /
                            Units::CoulombicConversionFactor;
    energy.polarization -= 0.5 * polarizability * double3::dot(moleculeElectricField[i], moleculeElectricField[i]);
  }

  return energy;
}

void System::computeTotalElectrostaticPotential() noexcept
{
  if (fixedFrameworkStoredEik.empty())
  {
    precomputeTotalRigidEnergy();
  }

  std::span<Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<Atom> moleculeAtomPositions = spanOfMoleculeAtoms();
  std::span<double> moleculeElectrostaticPotential = spanOfMoleculeElectrostaticPotential();

  std::fill(moleculeElectrostaticPotential.begin(), moleculeElectrostaticPotential.end(), 0.0);

  Interactions::computeInterMolecularElectrostaticPotential(forceField, simulationBox, moleculeElectrostaticPotential,
                                                            moleculeAtomPositions);

  Interactions::computeFrameworkMoleculeElectrostaticPotential(
      forceField, simulationBox, moleculeElectrostaticPotential, frameworkAtomPositions, moleculeAtomPositions);

  Interactions::computeEwaldFourierElectrostaticPotential(
      eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, moleculeElectrostaticPotential, forceField,
      simulationBox, components, numberOfMoleculesPerComponent, moleculeAtomPositions);
}

void System::computeTotalElectricField() noexcept
{
  if (fixedFrameworkStoredEik.empty())
  {
    precomputeTotalRigidEnergy();
  }

  std::span<Atom> frameworkAtomPositions = spanOfFrameworkAtoms();
  std::span<Atom> moleculeAtomPositions = spanOfMoleculeAtoms();
  std::span<double3> moleculeElectricField = spanOfMoleculeElectricField();

  std::fill(moleculeElectricField.begin(), moleculeElectricField.end(), double3(0.0, 0.0, 0.0));

  Interactions::computeInterMolecularElectricField(forceField, simulationBox, moleculeElectricField,
                                                   moleculeAtomPositions);

  Interactions::computeFrameworkMoleculeElectricField(forceField, simulationBox, moleculeElectricField,
                                                      frameworkAtomPositions, moleculeAtomPositions);

  Interactions::computeEwaldFourierElectricField(eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik,
                                                 forceField, simulationBox, moleculeElectricField, components,
                                                 numberOfMoleculesPerComponent, moleculeAtomPositions);
}

std::pair<EnergyStatus, double3x3> System::computeMolecularPressure()
{
  // Scratch buffer so molecular-pressure sampling does not mutate live MD site gradients.
  // Strain-derivative routines write intermolecular/framework forces here for the atomic-to-molecular
  // virial correction. Intramolecular bonded forces are intentionally omitted (they cancel in the
  // molecular virial) and must not overwrite AtomDynamics used by Velocity-Verlet.
  std::vector<AtomDynamics> pressureMoleculeDynamics(spanOfMoleculeAtoms().size());
  std::span<AtomDynamics> pressureDynamics(pressureMoleculeDynamics);

  const std::span<const Atom> moleculeAtoms = spanOfMoleculeAtoms();

  // Polarization rides along with the real-space Coulomb strain loops: they gather the per-atom electric
  // field and its cell-strain response (molecular COM scaling) in one pass, so the pair walk is not
  // repeated. The Ewald reciprocal framework field and the final contraction into the polarization energy
  // and strain tensor happen afterwards in computePolarizationMolecularPressureStrain. The polarization
  // forces are intentionally kept out of 'pressureDynamics': the polarization strain below already uses
  // COM arms, so it must not participate in the atomic-to-molecular virial correction.
  const bool gatherPolarization = forceField.computePolarization && forceField.useCharge && !moleculeAtoms.empty();
  std::vector<double3> polarizationField;
  std::vector<std::array<double3, 9>> polarizationFieldStrain;
  std::vector<double3> polarizationComOffset;
  std::vector<double> polarizationPolarizability;
  Interactions::PolarizationFieldStrain polarizationGather{};
  const Interactions::PolarizationFieldStrain* frameworkGather = nullptr;
  const Interactions::PolarizationFieldStrain* interGather = nullptr;
  if (gatherPolarization)
  {
    polarizationField.assign(moleculeAtoms.size(), double3(0.0, 0.0, 0.0));
    polarizationFieldStrain.assign(moleculeAtoms.size(), {});
    polarizationComOffset.assign(moleculeAtoms.size(), double3(0.0, 0.0, 0.0));
    polarizationPolarizability.resize(moleculeAtoms.size());
    for (std::size_t i = 0; i < moleculeAtoms.size(); ++i)
    {
      // Scaled by the atom's Coulomb scaling: fractional (CFCMC) molecules decouple from the field.
      polarizationPolarizability[i] =
          moleculeAtoms[i].scalingCoulomb *
          forceField.pseudoAtoms[static_cast<std::size_t>(moleculeAtoms[i].type)].polarizability /
          Units::CoulombicConversionFactor;
    }

    // Mass-weighted COM offsets from the current atom positions for every molecule (rigid and flexible),
    // matching the COM-scaling volume move.
    for (const Molecule& molecule : moleculeData)
    {
      double totalMass = 0.0;
      double3 com(0.0, 0.0, 0.0);
      for (std::size_t k = 0; k < molecule.numberOfAtoms; ++k)
      {
        const Atom& atom = moleculeAtoms[molecule.atomIndex + k];
        const double mass = forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].mass;
        com += mass * atom.position;
        totalMass += mass;
      }
      com = com / totalMass;
      for (std::size_t k = 0; k < molecule.numberOfAtoms; ++k)
      {
        polarizationComOffset[molecule.atomIndex + k] = moleculeAtoms[molecule.atomIndex + k].position - com;
      }
    }

    polarizationGather = {
        std::span<double3>(polarizationField), std::span<std::array<double3, 9>>(polarizationFieldStrain),
        std::span<const double3>(polarizationComOffset), std::span<const double>(polarizationPolarizability)};
    frameworkGather = &polarizationGather;
    // The inter-molecular field obeys the omit flags of the polarization model; the routine's own
    // omitInterInteractions early-out covers the remaining case.
    if (!forceField.omitInterPolarization) interGather = &polarizationGather;
  }

  std::pair<EnergyStatus, double3x3> pressureInfo = Interactions::computeFrameworkMoleculeEnergyStrainDerivative(
      forceField, framework, interpolationGrids, components, simulationBox, spanOfFrameworkAtoms(),
      spanOfMoleculeAtoms(), pressureDynamics, frameworkGather);

  pressureInfo.first.translationalKineticEnergy = runningEnergies.translationalKineticEnergy;
  pressureInfo.first.rotationalKineticEnergy = runningEnergies.rotationalKineticEnergy;
  pressureInfo.first.noseHooverEnergy = runningEnergies.NoseHooverEnergy;

  pressureInfo = pairSum(
      pressureInfo, Interactions::computeInterMolecularEnergyStrainDerivative(
                        forceField, components, simulationBox, spanOfMoleculeAtoms(), pressureDynamics, interGather));

  pressureInfo =
      pairSum(pressureInfo, Interactions::computeEwaldFourierEnergyStrainDerivative(
                                eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, storedEik, forceField,
                                simulationBox, framework, components, numberOfMoleculesPerComponent,
                                spanOfMoleculeAtoms(), pressureDynamics, netChargeFramework, netChargePerComponent));

  std::size_t molecule_index = 0;
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (numberOfMoleculesPerComponent[i] > 0)
    {
      std::span<const Molecule> span_molecules = {&moleculeData[molecule_index], numberOfMoleculesPerComponent[i]};
      RunningEnergy runningIntraEnergy = Interactions::computeIntraMolecularEnergy(
          components[i].intraMolecularPotentials, span_molecules, spanOfMoleculeAtoms());

      pressureInfo.first.intraComponentEnergies[i].bond += runningIntraEnergy.bond;
      pressureInfo.first.intraComponentEnergies[i].ureyBradley += runningIntraEnergy.ureyBradley;
      pressureInfo.first.intraComponentEnergies[i].bend += runningIntraEnergy.bend;
      pressureInfo.first.intraComponentEnergies[i].inversionBend += runningIntraEnergy.inversionBend;
      pressureInfo.first.intraComponentEnergies[i].outOfPlaneBend += runningIntraEnergy.outOfPlaneBend;
      pressureInfo.first.intraComponentEnergies[i].torsion += runningIntraEnergy.torsion;
      pressureInfo.first.intraComponentEnergies[i].improperTorsion += runningIntraEnergy.improperTorsion;
      pressureInfo.first.intraComponentEnergies[i].bondBond += runningIntraEnergy.bondBond;
      pressureInfo.first.intraComponentEnergies[i].bondBend += runningIntraEnergy.bondBend;
      pressureInfo.first.intraComponentEnergies[i].bondTorsion += runningIntraEnergy.bondTorsion;
      pressureInfo.first.intraComponentEnergies[i].bendBend += runningIntraEnergy.bendBend;
      pressureInfo.first.intraComponentEnergies[i].bendTorsion += runningIntraEnergy.bendTorsion;
      pressureInfo.first.intraComponentEnergies[i].vanDerWaals += runningIntraEnergy.intraVDW;
      pressureInfo.first.intraComponentEnergies[i].coulomb += runningIntraEnergy.intraCoul;

      // Intramolecular potentials (bonds, bends, torsions, ...) must NOT contribute to the molecular
      // (center-of-mass based) pressure: internal forces sum to zero over each molecule and cancel in the
      // molecular virial. Their strain derivative is therefore intentionally not accumulated here, and their
      // gradients are deliberately kept out of 'pressureDynamics' so they cannot pollute the
      // atomic-to-molecular correction term computed below. Only the energy is recorded (above) for reporting.
    }

    molecule_index += numberOfMoleculesPerComponent[i];
  }

  if (gatherPolarization)
  {
    // Complete the gathered real-space field with the Ewald reciprocal framework contribution and contract
    // into the polarization energy and its (unsymmetrized) strain-derivative tensor; the symmetrization at
    // the end of this function averages the off-diagonal pairs, matching the symmetric strain generators.
    const auto [polarizationEnergy, polarizationStrain] = Interactions::computePolarizationMolecularPressureStrain(
        *this, polarizationField, polarizationFieldStrain, polarizationComOffset, polarizationPolarizability);

    pressureInfo.first.polarizationEnergy = EnergyDuDlambda(polarizationEnergy, 0.0);
    pressureInfo.second += polarizationStrain;
  }

  if (useMBX)
  {
#ifdef BUILD_MBX
    // The pressure tensor remains the classical molecular-pressure estimate. Expose the current MBX energy
    // alongside it for reporting without erasing the classical decomposition used to construct the tensor.
    const RunningEnergy mbx = Interactions::computeMBXEnergySystem(*this, components, simulationBox, framework,
                                                                   spanOfFrameworkAtoms(), spanOfMoleculeAtoms());
    pressureInfo.first.useMBX = true;
    pressureInfo.first.mbxEnergy = mbx.mbxEnergy;
#else
    throw std::runtime_error(
        "[System]: MBX pressure reporting was requested, but this binary was built without BUILD_MBX support");
#endif
  }

  pressureInfo.first.sumTotal();

  double pressureTailCorrection = 0.0;
  // Tail correction to the (excess) pressure virial. The per-pair integral 'tailCorrectionPressure' equals
  // Integrate[U'(r) r^3, {r, rc, Inf}]; the isotropic virial tail correction is -(2 pi / 3 V) Sum_ij n_i n_j <integral>
  // (see RASPA2 CalculateTailCorrection). The overall minus sign and the 1/3 are essential: the van der Waals tail is
  // attractive and must LOWER the pressure. The extra global negation of 'pressureInfo.second' below turns the
  // '-= pressureTailCorrection' into the correct negative diagonal contribution.
  double preFactor = -2.0 * std::numbers::pi / (3.0 * simulationBox.volume);
  for (std::vector<Atom>::iterator it1 = atomData.begin(); it1 != atomData.end(); ++it1)
  {
    std::size_t typeA = static_cast<std::size_t>(it1->type);
    double scalingVDWA = it1->scalingVDW;

    pressureTailCorrection += scalingVDWA * scalingVDWA * preFactor * forceField(typeA, typeA).tailCorrectionPressure;

    for (std::vector<Atom>::iterator it2 = it1 + 1; it2 != atomData.end(); ++it2)
    {
      std::size_t typeB = static_cast<std::size_t>(it2->type);
      double scalingVDWB = it2->scalingVDW;

      pressureTailCorrection +=
          scalingVDWA * scalingVDWB * 2.0 * preFactor * forceField(typeA, typeB).tailCorrectionPressure;
    }
  }

  pressureInfo.second.ax -= pressureTailCorrection;
  pressureInfo.second.by -= pressureTailCorrection;
  pressureInfo.second.cz -= pressureTailCorrection;

  // Correct rigid molecule contribution using the constraints forces.
  // Molecule::atomIndex indexes the molecule-atom span (framework atoms excluded).
  double3x3 correctionTerm{};
  for (Molecule& molecule : moleculeData)
  {
    const std::span<const Atom> span = moleculeAtoms.subspan(molecule.atomIndex, molecule.numberOfAtoms);
    const std::span<AtomDynamics> spanDynamics = {&pressureMoleculeDynamics[molecule.atomIndex],
                                                  molecule.numberOfAtoms};

    double totalMass = 0.0;
    double3 com(0.0, 0.0, 0.0);
    for (const Atom& atom : span)
    {
      double mass = forceField.pseudoAtoms[static_cast<std::size_t>(atom.type)].mass;
      com += mass * atom.position;
      totalMass += mass;
    }
    com = com / totalMass;

    for (std::size_t k = 0; k < span.size(); ++k)
    {
      const double3 position = span[k].position;
      const double3 gradient = spanDynamics[k].gradient;

      correctionTerm.ax += (position.x - com.x) * gradient.x;
      correctionTerm.ay += (position.x - com.x) * gradient.y;
      correctionTerm.az += (position.x - com.x) * gradient.z;

      correctionTerm.bx += (position.y - com.y) * gradient.x;
      correctionTerm.by += (position.y - com.y) * gradient.y;
      correctionTerm.bz += (position.y - com.y) * gradient.z;

      correctionTerm.cx += (position.z - com.z) * gradient.x;
      correctionTerm.cy += (position.z - com.z) * gradient.y;
      correctionTerm.cz += (position.z - com.z) * gradient.z;
    }
  }

  pressureInfo.second = -(pressureInfo.second - correctionTerm);

  double temp = 0.5 * (pressureInfo.second.ay + pressureInfo.second.bx);
  pressureInfo.second.ay = pressureInfo.second.bx = temp;
  temp = 0.5 * (pressureInfo.second.az + pressureInfo.second.cx);
  pressureInfo.second.az = pressureInfo.second.cx = temp;
  temp = 0.5 * (pressureInfo.second.bz + pressureInfo.second.cy);
  pressureInfo.second.bz = pressureInfo.second.cy = temp;

  return pressureInfo;
}

void System::checkCartesianPositions()
{
  std::span<Atom> moleculeAtomPositions = spanOfMoleculeAtoms();

  std::size_t index{};
  for (Molecule& molecule : moleculeData)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
    if (components[molecule.componentId].rigid)
    {
      simd_quatd q = molecule.orientation;
      double3x3 M = double3x3::buildRotationMatrixInverse(q);

      for (std::size_t i = 0; i != span.size(); i++)
      {
        double3 expandedPosition =
            molecule.centerOfMassPosition + M * components[molecule.componentId].atoms[i].position;
        if ((std::abs(span[i].position.x - expandedPosition.x) > 1e-5) ||
            (std::abs(span[i].position.y - expandedPosition.y) > 1e-5) ||
            (std::abs(span[i].position.z - expandedPosition.z) > 1e-5))
        {
          throw std::runtime_error(
              std::format("Difference detected between atom position ({} {} {}) and position generated from "
                          "quaternion ({} {} {})\n",
                          span[i].position.x, span[i].position.y, span[i].position.z, expandedPosition.x,
                          expandedPosition.y, expandedPosition.z));
        }
      }
    }
    index += molecule.numberOfAtoms;
  }
}

void System::setThermostat(const std::optional<Thermostat>& thermo)
{
  if (thermo.has_value())
  {
    thermostat = Thermostat(temperature, timeStep, translationalDegreesOfFreedom, rotationalDegreesOfFreedom,
                            thermo->thermostatChainLength, thermo->numberOfYoshidaSuzukiSteps,
                            thermo->timeScaleParameterThermostat);
  }
}

void System::setThermobarostat(const std::optional<Thermobarostat>& barostat)
{
  if (!barostat.has_value())
  {
    thermobarostat.reset();
    return;
  }
  molecularDynamicsEnsemble = barostat->ensemble;
  thermobarostat = Thermobarostat(barostat->ensemble, barostat->cellType, barostat->monoclinicAngle, temperature,
                                  pressure, timeStep, translationalDegreesOfFreedom, barostat->chainLength,
                                  barostat->numberOfYoshidaSuzukiSteps, barostat->timeScaleParameterBarostat);
  thermobarostat->numberOfRespaSteps = barostat->numberOfRespaSteps;
}

void System::setSamplePDBMovie(const std::optional<SampleMovie>& movie)
{
  if (movie.has_value())
  {
    samplePDBMovie = movie;
  }
}

void System::updateSamplePDBMovie(std::size_t systemId, std::size_t currentCycle)
{
  if (samplePDBMovie.has_value())
  {
    samplePDBMovie->update(forceField, systemId, simulationBox, spanOfMoleculeAtoms(), components,
                           numberOfMoleculesPerComponent, currentCycle, spanOfFrameworkAtoms());
  }
}

void System::setNumberOfMoleculesHistogram(const std::optional<PropertyNumberOfMoleculesHistogram>& hist)
{
  if (hist.has_value())
  {
    averageNumberOfMoleculesHistogram = PropertyNumberOfMoleculesHistogram(
        hist->numberOfBlocks, components.size(), hist->range, hist->sampleEvery, hist->writeEvery);
  }
}

void System::setAverageEnergyHistogram(const std::optional<PropertyEnergyHistogram>& hist)
{
  if (hist.has_value())
  {
    averageEnergyHistogram = hist;
  }
}

void System::setPropertyDensityGrid(const std::optional<PropertyDensityGrid>& grid)
{
  if (grid.has_value())
  {
    propertyDensityGrid = grid;
  }
}

void System::setPropertyNumberOfMoleculesEvolution(std::optional<PropertyNumberOfMoleculesEvolution> property)
{
  if (property.has_value())
  {
    propertyNumberOfMoleculesEvolution = property;
  }
}

void System::setPropertyVolumeEvolution(std::optional<PropertyVolumeEvolution> property)
{
  if (property.has_value())
  {
    propertyVolumeEvolution = property;
  }
}

void System::setPropertyConservedEnergyEvolution(std::optional<PropertyConservedEnergyEvolution> property)
{
  if (property.has_value())
  {
    propertyConservedEnergyEvolution = property;
  }
}

void System::setPropertyConventionalRDF(const std::optional<PropertyConventionalRadialDistributionFunction>& rdf)
{
  if (rdf.has_value())
  {
    propertyConventionalRadialDistributionFunction = PropertyConventionalRadialDistributionFunction(
        5, forceField.pseudoAtoms.size(), rdf->numberOfBins, 12.0, rdf->sampleEvery, rdf->writeEvery);
  }
}

void System::setPropertyRDF(const std::optional<PropertyRadialDistributionFunction>& rdf)
{
  if (rdf.has_value())
  {
    propertyRadialDistributionFunction = PropertyRadialDistributionFunction(
        5, forceField.pseudoAtoms.size(), rdf->numberOfBins, rdf->range, rdf->sampleEvery, rdf->writeEvery);
  }
}

bool System::forceBasedRDFSampleDue(std::size_t currentCycle) const
{
  return propertyRadialDistributionFunction.has_value() && propertyRadialDistributionFunction->sampleEvery > 0uz &&
         (currentCycle % propertyRadialDistributionFunction->sampleEvery == 0uz);
}

void System::sampleForceBasedRDFFromCurrentGradients(std::size_t currentCycle, std::size_t currentBlock)
{
  if (!propertyRadialDistributionFunction.has_value()) return;
  propertyRadialDistributionFunction->sample(simulationBox, spanOfFrameworkAtoms(), spanOfFrameworkDynamics(),
                                             moleculeData, spanOfMoleculeAtoms(), spanOfMoleculeDynamics(),
                                             currentCycle, currentBlock);
}

void System::sampleForceBasedRDFWithFullGradients(std::size_t currentCycle, std::size_t currentBlock)
{
  if (!forceBasedRDFSampleDue(currentCycle)) return;
  precomputeTotalGradients();
  sampleForceBasedRDFFromCurrentGradients(currentCycle, currentBlock);
}

void System::setPropertyMSD(const std::optional<PropertyMeanSquaredDisplacement>& msd)
{
  if (msd.has_value())
  {
    propertyMSD = PropertyMeanSquaredDisplacement(numberOfMoleculesPerComponent, moleculeData.size(), timeStep,
                                                  msd->numberOfBlockElementsMSD, msd->sampleEvery, msd->writeEvery);
  }
}

void System::setPropertyVACF(const std::optional<PropertyVelocityAutoCorrelationFunction>& vacf)
{
  if (vacf.has_value())
  {
    propertyVACF = PropertyVelocityAutoCorrelationFunction(numberOfMoleculesPerComponent, moleculeData.size(), timeStep,
                                                           vacf->numberOfBuffersVACF, vacf->bufferLengthVACF,
                                                           vacf->sampleEvery, vacf->writeEvery);
  }
}

std::string System::writeMBXStatus() const
{
  if (!useMBX) return {};

  const double frameworkPermanentInternal = elecPermFrameworkMBX / Units::EnergyToKCalPerMol;
  const std::string energyTermsDestination =
      !writeEnergyLog ? std::string{"disabled"}
                      : (energyTermsLogSink.filePath.empty() ? std::string{"enabled"} : energyTermsLogSink.filePath);
  const bool hasConventionalWidom = std::ranges::any_of(components, [](const Component& component) {
    return component.mc_moves_probabilities.getProbability(Move::Types::Widom) > 0.0;
  });
  const std::string widomEnergyTermsDestination =
      !writeEnergyLog || !hasConventionalWidom
          ? std::string{"disabled"}
          : (energyTermsLogSink.widomFilePath.empty() ? std::string{"enabled"}
                                                       : energyTermsLogSink.widomFilePath);
  return std::format(
      "MBX energy model\n"
      "===============================================================================\n"
      "Settings file: {}\n"
      "Accepted-move energy terms: {}\n"
      "Widom-trial energy terms: {}\n"
      "Zero-loading Widom heat: {}\n"
      "Current retained MBX guest energy: {: .6e} [{}]\n"
      "Bare-framework permanent electrostatics removed from MBX: {: .6e} [{}]\n\n",
      mbxSettingsFilePath, energyTermsDestination, widomEnergyTermsDestination,
      computeZeroLoadingHeatOfAdsorption ? "enabled" : "disabled",
      runningEnergies.mbxEnergy * Units::EnergyToKelvin,
      Units::displayedUnitOfEnergyString, frameworkPermanentInternal * Units::EnergyToKelvin,
      Units::displayedUnitOfEnergyString);
}

void System::configureEnergyTermsLog(std::string filePath, bool append)
{
  energyTermsLogSink = EnergyTermsLogSink{};
  if (!writeEnergyLog) return;

  const std::filesystem::path path = std::filesystem::path(std::move(filePath)).lexically_normal();
  if (path.empty())
  {
    throw std::invalid_argument("Energy-terms log path must not be empty");
  }

  std::error_code error;
  if (!path.parent_path().empty())
  {
    std::filesystem::create_directories(path.parent_path(), error);
    if (error)
    {
      throw std::runtime_error(std::format("Unable to create energy-terms output directory '{}': {}",
                                           path.parent_path().string(), error.message()));
    }
  }

  bool existingContent = false;
  if (append)
  {
    const bool exists = std::filesystem::exists(path, error);
    if (error)
    {
      throw std::runtime_error(
          std::format("Unable to inspect energy-terms log '{}': {}", path.string(), error.message()));
    }
    if (exists)
    {
      existingContent = std::filesystem::file_size(path, error) > 0;
      if (error)
      {
        throw std::runtime_error(
            std::format("Unable to inspect energy-terms log '{}': {}", path.string(), error.message()));
      }
    }
  }

  const std::ios::openmode mode = std::ios::out | (append ? std::ios::app : std::ios::trunc);
  auto stream = std::make_unique<std::ofstream>(path, mode);
  if (!*stream)
  {
    throw std::runtime_error(std::format("Unable to open energy-terms log '{}'", path.string()));
  }

  energyTermsLogSink.filePath = path.string();
  energyTermsLogSink.headerWritten = existingContent;
  energyTermsLogSink.stream = std::move(stream);

  if (!energyTermsLogSink.headerWritten)
  {
    if (useMBX)
    {
      // 1B is intentionally excluded from both the retained MBX total and the historical RASPA-MBX CSV schema.
      *energyTermsLogSink.stream
          << "type,component,N,total,hg_VDW,hg_tail,mbx_tot,E2b,E3b,E4b,Edisp,Eelec_perm,Eelec_ind,E_diff,Pacc\n";
    }
    else
    {
      *energyTermsLogSink.stream
          << "type,component,N,total,hg_VDW,gg_VDW,tail,hg_Charge,gg_Charge,E_ewald,E_diff,Pacc\n";
    }
    std::flush(*energyTermsLogSink.stream);
    energyTermsLogSink.headerWritten = true;
  }

  const bool hasConventionalWidom = std::ranges::any_of(components, [](const Component& component) {
    return component.mc_moves_probabilities.getProbability(Move::Types::Widom) > 0.0;
  });
  if (!hasConventionalWidom) return;

  const std::filesystem::path widomPath = path.parent_path() / std::format("widom_{}", path.filename().string());
  bool existingWidomContent = false;
  if (append)
  {
    const bool exists = std::filesystem::exists(widomPath, error);
    if (error)
    {
      throw std::runtime_error(
          std::format("Unable to inspect Widom energy-terms log '{}': {}", widomPath.string(), error.message()));
    }
    if (exists)
    {
      existingWidomContent = std::filesystem::file_size(widomPath, error) > 0;
      if (error)
      {
        throw std::runtime_error(
            std::format("Unable to inspect Widom energy-terms log '{}': {}", widomPath.string(), error.message()));
      }
    }
  }

  auto widomStream = std::make_unique<std::ofstream>(widomPath, mode);
  if (!*widomStream)
  {
    throw std::runtime_error(std::format("Unable to open Widom energy-terms log '{}'", widomPath.string()));
  }
  energyTermsLogSink.widomFilePath = widomPath.string();
  energyTermsLogSink.widomHeaderWritten = existingWidomContent;
  energyTermsLogSink.widomStream = std::move(widomStream);

  if (!energyTermsLogSink.widomHeaderWritten)
  {
    if (useMBX)
    {
      *energyTermsLogSink.widomStream
          << "component,N,trial_N,current_total,trial_total,hg_VDW,hg_tail,mbx_tot,E2b,E3b,E4b,Edisp,"
             "Eelec_perm,Eelec_ind,E_insert,weight,weighted_E_insert\n";
    }
    else
    {
      *energyTermsLogSink.widomStream
          << "component,N,trial_N,current_total,trial_total,hg_VDW,gg_VDW,tail,hg_Charge,gg_Charge,E_ewald,"
             "E_insert,weight,weighted_E_insert\n";
    }
    std::flush(*energyTermsLogSink.widomStream);
    energyTermsLogSink.widomHeaderWritten = true;
  }
}

void System::flushEnergyTermsLog() const
{
  if (energyTermsLogSink.stream) std::flush(*energyTermsLogSink.stream);
  if (energyTermsLogSink.widomStream) std::flush(*energyTermsLogSink.widomStream);
}

void System::writeAcceptedEnergyLog(std::string_view moveType, std::size_t componentId,
                                    const RunningEnergy& totalEnergy, std::span<const double> mbxTerms,
                                    double energyDifference, double acceptanceProbability) const
{
  if (!writeEnergyLog) return;
  if (useMBX && mbxTerms.size() != 7)
  {
    throw std::invalid_argument(
        std::format("Accepted MBX energy log requires 7 subterms, received {}", mbxTerms.size()));
  }

  std::ostream& target = energyTermsLogSink.stream ? static_cast<std::ostream&>(*energyTermsLogSink.stream) : std::cerr;
  std::osyncstream output(target);
  if (!energyTermsLogSink.headerWritten)
  {
    if (useMBX)
    {
      // 1B is intentionally excluded from both the retained MBX total and the historical RASPA-MBX CSV schema.
      output << "type,component,N,total,hg_VDW,hg_tail,mbx_tot,E2b,E3b,E4b,Edisp,Eelec_perm,Eelec_ind,"
                "E_diff,Pacc\n";
    }
    else
    {
      output << "type,component,N,total,hg_VDW,gg_VDW,tail,hg_Charge,gg_Charge,E_ewald,E_diff,Pacc\n";
    }
    energyTermsLogSink.headerWritten = true;
  }

  const bool allComponents = componentId == std::numeric_limits<std::size_t>::max();
  const std::string component = allComponents ? std::string{"all"} : std::to_string(componentId);
  const std::size_t moleculeCount =
      allComponents ? numberOfMolecules() : numberOfIntegerMoleculesPerComponent.at(componentId);
  output << std::setprecision(17) << moveType << ',' << component << ',' << moleculeCount << ','
         << totalEnergy.potentialEnergy() << ',';

  if (useMBX)
  {
    output << totalEnergy.frameworkMoleculeVDW << ',' << totalEnergy.tail << ',' << totalEnergy.mbxEnergy << ','
           << mbxTerms[1] / Units::EnergyToKCalPerMol << ',' << mbxTerms[2] / Units::EnergyToKCalPerMol << ','
           << mbxTerms[3] / Units::EnergyToKCalPerMol << ',' << mbxTerms[4] / Units::EnergyToKCalPerMol << ','
           << mbxTerms[5] / Units::EnergyToKCalPerMol << ',' << mbxTerms[6] / Units::EnergyToKCalPerMol << ',';
  }
  else
  {
    output << totalEnergy.frameworkMoleculeVDW << ',' << totalEnergy.moleculeMoleculeVDW << ',' << totalEnergy.tail
           << ',' << totalEnergy.frameworkMoleculeCharge << ',' << totalEnergy.moleculeMoleculeCharge << ','
           << totalEnergy.ewald_fourier + totalEnergy.ewald_self + totalEnergy.ewald_exclusion << ',';
  }
  output << energyDifference << ',' << acceptanceProbability << '\n' << std::flush_emit;
}

void System::writeWidomEnergyLog(std::size_t componentId, const RunningEnergy& trialTotalEnergy,
                                 std::span<const double> mbxTerms, double insertionEnergy,
                                 double widomWeight) const
{
  if (!writeEnergyLog) return;
  if (useMBX && mbxTerms.size() != 7)
  {
    throw std::invalid_argument(
        std::format("Widom MBX energy log requires 7 subterms, received {}", mbxTerms.size()));
  }

  std::ostream& target = energyTermsLogSink.widomStream
                             ? static_cast<std::ostream&>(*energyTermsLogSink.widomStream)
                             : std::cerr;
  std::osyncstream output(target);
  if (!energyTermsLogSink.widomHeaderWritten)
  {
    if (useMBX)
    {
      output << "component,N,trial_N,current_total,trial_total,hg_VDW,hg_tail,mbx_tot,E2b,E3b,E4b,Edisp,"
                "Eelec_perm,Eelec_ind,E_insert,weight,weighted_E_insert\n";
    }
    else
    {
      output << "component,N,trial_N,current_total,trial_total,hg_VDW,gg_VDW,tail,hg_Charge,gg_Charge,E_ewald,"
                "E_insert,weight,weighted_E_insert\n";
    }
    energyTermsLogSink.widomHeaderWritten = true;
  }

  const std::size_t moleculeCount = numberOfIntegerMoleculesPerComponent.at(componentId);
  if (moleculeCount == std::numeric_limits<std::size_t>::max())
  {
    throw std::overflow_error("Widom trial molecule count cannot be represented");
  }
  output << std::setprecision(17) << componentId << ',' << moleculeCount << ',' << moleculeCount + 1 << ','
         << runningEnergies.potentialEnergy() << ',' << trialTotalEnergy.potentialEnergy() << ',';

  if (useMBX)
  {
    output << trialTotalEnergy.frameworkMoleculeVDW << ',' << trialTotalEnergy.tail << ','
           << trialTotalEnergy.mbxEnergy << ',' << mbxTerms[1] / Units::EnergyToKCalPerMol << ','
           << mbxTerms[2] / Units::EnergyToKCalPerMol << ',' << mbxTerms[3] / Units::EnergyToKCalPerMol << ','
           << mbxTerms[4] / Units::EnergyToKCalPerMol << ',' << mbxTerms[5] / Units::EnergyToKCalPerMol << ','
           << mbxTerms[6] / Units::EnergyToKCalPerMol << ',';
  }
  else
  {
    output << trialTotalEnergy.frameworkMoleculeVDW << ',' << trialTotalEnergy.moleculeMoleculeVDW << ','
           << trialTotalEnergy.tail << ',' << trialTotalEnergy.frameworkMoleculeCharge << ','
           << trialTotalEnergy.moleculeMoleculeCharge << ','
           << trialTotalEnergy.ewald_fourier + trialTotalEnergy.ewald_self + trialTotalEnergy.ewald_exclusion << ',';
  }
  output << insertionEnergy << ',' << widomWeight << ',' << insertionEnergy * widomWeight << '\n'
         << std::flush_emit;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const System& s)
{
  archive << s.versionNumber;

  archive << s.temperature;
  archive << s.pressure;
  archive << s.input_pressure;
  archive << s.pressureTensorDiagonal;
  archive << s.input_pressureTensorDiagonal;
  archive << s.beta;
  archive << static_cast<std::uint8_t>(s.cellMinimizationType);
  archive << static_cast<std::uint8_t>(s.monoclinicAngleType);

  archive << s.heliumVoidFraction;

  archive << s.numberOfFrameworks;
  archive << s.numberOfFrameworkAtoms;
  archive << s.numberOfRigidFrameworkAtoms;

  archive << s.framework;
  archive << s.components;

  archive << s.equationOfState;
  archive << static_cast<std::uint8_t>(s.molecularDynamicsEnsemble);
  archive << s.thermostat;
  archive << s.thermobarostat;

  archive << s.loadings;

  archive << s.swappableComponents;
  archive << s.initialNumberOfMolecules;

  archive << s.numberOfMoleculesPerComponent;
  archive << s.numberOfIntegerMoleculesPerComponent;
  archive << s.numberOfFractionalMoleculesPerComponent;
  archive << s.numberOfGCFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfPairGCFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfPairSwapFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfGibbsFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfParallelReactionFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC;
  archive << s.numberOfReactionFractionalMoleculesPerComponent_CFCMC;

  archive << s.idealGasEnergiesPerComponent;

  archive << s.forceField;
  archive << s.hasExternalField;
  archive << s.useMBX;
  archive << s.mbxSettingsFilePath;
  archive << s.writeEnergyLog;
  archive << s.computeZeroLoadingHeatOfAdsorption;
  archive << s.elecPermFrameworkMBX;

  archive << s.numberOfPseudoAtoms;
  archive << s.totalNumberOfPseudoAtoms;

  archive << s.translationalCenterOfMassConstraint;
  archive << s.translationalDegreesOfFreedom;
  archive << s.rotationalDegreesOfFreedom;

  archive << s.timeStep;

  archive << s.simulationBox;
  archive << s.containsTheFractionalMolecule;

  archive << s.atomData;
  archive << s.atomDynamics;
  archive << s.moleculeData;
  archive << s.electricPotential;
  archive << s.electricField;
  archive << s.electricFieldNew;

  archive << s.conservedEnergy;
  archive << s.referenceEnergy;
  archive << s.accumulatedDrift;

  archive << s.rigidEnergies;
  archive << s.runningEnergies;

  archive << s.currentExcessPressureTensor;
  archive << s.currentEnergyStatus;

  archive << s.numberOfHybridMCSteps;

  archive << s.eik_xy;
  archive << s.eik_x;
  archive << s.eik_y;
  archive << s.eik_z;
  archive << s.storedEik;
  archive << s.fixedFrameworkStoredEik;
  archive << s.trialEik;
  archive << s.CoulombicFourierEnergySingleIon;
  archive << s.netCharge;
  archive << s.netChargeFramework;
  archive << s.netChargeAdsorbates;
  archive << s.netChargePerComponent;

  archive << s.mc_moves_probabilities;
  archive << s.mc_moves_statistics;
  archive << s.mc_moves_cputime;

  archive << s.reactions;
  archive << s.tmmc;

  archive << s.averageEnergies;
  archive << s.averageLoadings;
  archive << s.averageEnthalpiesOfAdsorption;
  archive << s.averagePartialMolarProperties;
  archive << s.averageTemperature;
  archive << s.averageTranslationalTemperature;
  archive << s.averageRotationalTemperature;
  archive << s.averagePressure;
  archive << s.averageSimulationBox;
  archive << s.propertyElasticConstantsFluctuation;
  archive << s.elasticConstantsSampleEvery;

  archive << s.samplePDBMovie;

  archive << s.propertyConventionalRadialDistributionFunction;
  archive << s.propertyRadialDistributionFunction;
  archive << s.propertyDensityGrid;
  archive << s.averageEnergyHistogram;
  archive << s.averageNumberOfMoleculesHistogram;
  archive << s.propertyMoleculeProperties;
  archive << s.propertyMSD;
  archive << s.propertyVACF;
  archive << s.writeLammpsData;

  archive << s.propertyNumberOfMoleculesEvolution;
  archive << s.propertyVolumeEvolution;
  archive << s.propertyConservedEnergyEvolution;

  // 'interpolationGrids' and 'externalFieldInterpolationGrid' are intentionally not serialized:
  // they are derived data that can dominate the checkpoint size and are rebuilt deterministically
  // on restart (MonteCarlo/MolecularDynamics::createInterpolationGrids).

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, System& s)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > s.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'System' at line {} in file {}\n", location.line(), location.file_name()));
  }

  archive >> s.temperature;
  archive >> s.pressure;
  archive >> s.input_pressure;
  archive >> s.pressureTensorDiagonal;
  archive >> s.input_pressureTensorDiagonal;
  archive >> s.beta;
  {
    std::uint8_t cellMinimizationType;
    std::uint8_t monoclinicAngleType;
    archive >> cellMinimizationType;
    archive >> monoclinicAngleType;
    s.cellMinimizationType = static_cast<CellMinimizationType>(cellMinimizationType);
    s.monoclinicAngleType = static_cast<MonoclinicAngleType>(monoclinicAngleType);
  }

  archive >> s.heliumVoidFraction;

  archive >> s.numberOfFrameworks;
  archive >> s.numberOfFrameworkAtoms;
  archive >> s.numberOfRigidFrameworkAtoms;

  archive >> s.framework;
  s.numberOfRigidFrameworkAtoms = s.framework && s.framework->rigid ? s.numberOfFrameworkAtoms : 0;
  archive >> s.components;

  archive >> s.equationOfState;
  {
    std::uint8_t molecularDynamicsEnsemble;
    archive >> molecularDynamicsEnsemble;
    s.molecularDynamicsEnsemble = static_cast<MolecularDynamicsEnsemble>(molecularDynamicsEnsemble);
  }
  archive >> s.thermostat;
  archive >> s.thermobarostat;

  archive >> s.loadings;

  archive >> s.swappableComponents;
  archive >> s.initialNumberOfMolecules;

  archive >> s.numberOfMoleculesPerComponent;
  archive >> s.numberOfIntegerMoleculesPerComponent;
  archive >> s.numberOfFractionalMoleculesPerComponent;
  archive >> s.numberOfGCFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfPairGCFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfPairSwapFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfGibbsFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfParallelReactionFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC;
  archive >> s.numberOfReactionFractionalMoleculesPerComponent_CFCMC;

  archive >> s.idealGasEnergiesPerComponent;

  archive >> s.forceField;
  archive >> s.hasExternalField;
  if (versionNumber >= 2)
  {
    archive >> s.useMBX;
    archive >> s.mbxSettingsFilePath;
    archive >> s.writeEnergyLog;
    if (versionNumber >= 3)
    {
      archive >> s.computeZeroLoadingHeatOfAdsorption;
    }
    else
    {
      s.computeZeroLoadingHeatOfAdsorption = false;
    }
    archive >> s.elecPermFrameworkMBX;
  }
  else
  {
    // Upstream version-1 System archives predate MBX state. Preserve stream alignment and default safely.
    s.useMBX = false;
    s.mbxSettingsFilePath.clear();
    s.writeEnergyLog = true;
    s.computeZeroLoadingHeatOfAdsorption = false;
    s.elecPermFrameworkMBX = 0.0;
  }
  s.energyTermsLogSink = System::EnergyTermsLogSink{};

  archive >> s.numberOfPseudoAtoms;
  archive >> s.totalNumberOfPseudoAtoms;

  archive >> s.translationalCenterOfMassConstraint;
  archive >> s.translationalDegreesOfFreedom;
  archive >> s.rotationalDegreesOfFreedom;

  archive >> s.timeStep;

  archive >> s.simulationBox;
  archive >> s.containsTheFractionalMolecule;

  archive >> s.atomData;
  archive >> s.atomDynamics;
  archive >> s.moleculeData;
  archive >> s.electricPotential;
  archive >> s.electricField;
  archive >> s.electricFieldNew;

  archive >> s.conservedEnergy;
  archive >> s.referenceEnergy;
  archive >> s.accumulatedDrift;

  archive >> s.rigidEnergies;
  archive >> s.runningEnergies;

  archive >> s.currentExcessPressureTensor;
  archive >> s.currentEnergyStatus;

  archive >> s.numberOfHybridMCSteps;

  archive >> s.eik_xy;
  archive >> s.eik_x;
  archive >> s.eik_y;
  archive >> s.eik_z;
  archive >> s.storedEik;
  archive >> s.fixedFrameworkStoredEik;
  archive >> s.trialEik;
  archive >> s.CoulombicFourierEnergySingleIon;
  archive >> s.netCharge;
  archive >> s.netChargeFramework;
  archive >> s.netChargeAdsorbates;
  archive >> s.netChargePerComponent;

  archive >> s.mc_moves_probabilities;
  archive >> s.mc_moves_statistics;
  archive >> s.mc_moves_cputime;

  archive >> s.reactions;
  archive >> s.tmmc;

  archive >> s.averageEnergies;
  archive >> s.averageLoadings;
  archive >> s.averageEnthalpiesOfAdsorption;
  archive >> s.averagePartialMolarProperties;
  archive >> s.averageTemperature;
  archive >> s.averageTranslationalTemperature;
  archive >> s.averageRotationalTemperature;
  archive >> s.averagePressure;
  archive >> s.averageSimulationBox;
  archive >> s.propertyElasticConstantsFluctuation;
  archive >> s.elasticConstantsSampleEvery;

  archive >> s.samplePDBMovie;

  archive >> s.propertyConventionalRadialDistributionFunction;
  archive >> s.propertyRadialDistributionFunction;
  archive >> s.propertyDensityGrid;
  archive >> s.averageEnergyHistogram;
  archive >> s.averageNumberOfMoleculesHistogram;
  archive >> s.propertyMoleculeProperties;
  archive >> s.propertyMSD;
  archive >> s.propertyVACF;
  archive >> s.writeLammpsData;

  archive >> s.propertyNumberOfMoleculesEvolution;
  archive >> s.propertyVolumeEvolution;
  archive >> s.propertyConservedEnergyEvolution;

  // archive >> s.columnNumberOfGridPoints;
  // archive >> s.columnTotalPressure;
  // archive >> s.columnPressureGradient;
  // archive >> s.columnVoidFraction;
  // archive >> s.columnParticleDensity;
  // archive >> s.columnEntranceVelocity;
  // archive >> s.columnLength;
  // archive >> s.columnTimeStep;
  // archive >> s.columnNumberOfTimeSteps;
  // archive >> s.columnAutoNumberOfTimeSteps;
  // archive >> s.mixturePredictionMethod;
  // archive >> s.pressure_range;
  // archive >> s.numberOfCarrierGases;
  // archive >> s.carrierGasComponent;
  // archive >> s.maxIsothermTerms;

  // The interpolation grids are not stored in the archive; size the empty containers here and let
  // the driver's restart path rebuild them (createInterpolationGrids).
  s.interpolationGrids =
      std::vector<std::optional<InterpolationEnergyGrid>>(s.forceField.pseudoAtoms.size() + 1, std::nullopt);
  s.externalFieldInterpolationGrid = std::nullopt;

  // Derived quantities: rebuild the aggregated tail-correction counts from the restored atoms.
  s.computeTailCorrectionCounts();

  s.currentEnergyStatus.useMBX = s.useMBX;
  if (s.useMBX)
  {
#ifndef BUILD_MBX
    throw std::runtime_error(
        "[System restart]: checkpoint requires MBX, but this RASPA binary was built without BUILD_MBX support");
#else
    std::error_code error;
    if (!std::filesystem::is_regular_file(s.mbxSettingsFilePath, error))
    {
      throw std::runtime_error(std::format("[System restart]: MBX settings file '{}' is unavailable{}",
                                           s.mbxSettingsFilePath,
                                           error ? std::format(" ({})", error.message()) : std::string{}));
    }
#endif
  }

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("System: Error in binary restart\n"));
  }
#endif

  return archive;
}

void System::writeRestartFile(std::size_t systemId)
{
  nlohmann::json json;

  if (!framework.has_value())
  {
    json["SimulationBox"] = simulationBox;
  }

  for (std::size_t component_id = 0; component_id < components.size(); ++component_id)
  {
    std::vector<double3> positions{};

    if (numberOfMoleculesPerComponent[component_id] > 0)
    {
      positions = spanOfIntegerAtomsOfComponent(component_id) | std::views::transform(&Atom::position) |
                  std::ranges::to<std::vector>();
    }
    json[components[component_id].name] = positions;
  }

  std::string fileNameString = std::format("output/restart_{}_{}.s{}.json", temperature, input_pressure, systemId);
  std::ofstream file(fileNameString);
  file << json.dump(2);
}

std::string System::repr() const { return std::string("system test"); }
