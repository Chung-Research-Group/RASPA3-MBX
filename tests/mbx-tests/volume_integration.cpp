#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import atom;
import molecule;
import double3;
import double3x3;
import int3;
import simd_quatd;
import pseudo_atom;
import vdwparameters;
import forcefield;
import component;
import system;
import simulationbox;
import running_energy;
import randomnumbers;
import json;
import interactions_mbx;
import interactions_intermolecular;
import interactions_ewald;
import mc_moves_volume;
import mc_moves_statistics;
import mc_moves_move_types;
import move_statistics;

namespace
{
// These are isolated branch fixtures, not an NPT trajectory or an ensemble test.
// Execute with OMP_NUM_THREADS=1. All MBX calls contain two valid rigid waters.
constexpr double energyTolerance = 1e-6;  // Native RASPA energy units.
constexpr double coordinateTolerance = 2e-12;
constexpr Move::Types volumeMoveType = Move::Types::VolumeChange;

constexpr std::string_view settings = R"({
  "MBX": {
    "box": [],
    "twobody_cutoff": 9.0,
    "threebody_cutoff": 5.0,
    "dipole_tolerance": 1e-9,
    "dipole_method": "cg",
    "alpha_ewald_elec": 0.6,
    "grid_density_elec": 2.5,
    "spline_order_elec": 6,
    "alpha_ewald_disp": 0.6,
    "grid_density_disp": 2.5,
    "spline_order_disp": 6,
    "ignore_2b_poly": [],
    "ignore_3b_poly": []
  }
})";

// Observe object bytes through volatile character loads. A by-value double
// classifier (even one using bit_cast) can itself fold to "finite" under the
// release finite-math assumptions. This oracle has no floating-value argument
// or floating arithmetic and cannot discard the required byte observations.
bool finiteBits(const double& value)
{
  static_assert(sizeof(double) == sizeof(std::uint64_t));
  const volatile unsigned char* bytes = reinterpret_cast<const volatile unsigned char*>(&value);
  std::array<unsigned char, sizeof(double)> observed{};
  for (std::size_t i = 0; i < observed.size(); ++i) observed[i] = bytes[i];
  return (std::bit_cast<std::uint64_t>(observed) & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000);
}

void expectVectorNear(const double3& actual, const double3& expected)
{
  EXPECT_NEAR(actual.x, expected.x, coordinateTolerance);
  EXPECT_NEAR(actual.y, expected.y, coordinateTolerance);
  EXPECT_NEAR(actual.z, expected.z, coordinateTolerance);
}

struct Cutoffs
{
  double frameworkVDW;
  double moleculeVDW;
  double coulomb;
  double alpha;
  int3 waveVectors;
  std::size_t reciprocalIntegerSquared;

  explicit Cutoffs(const ForceField& forceField)
      : frameworkVDW(forceField.cutOffFrameworkVDW),
        moleculeVDW(forceField.cutOffMoleculeVDW),
        coulomb(forceField.cutOffCoulomb),
        alpha(forceField.EwaldAlpha),
        waveVectors(forceField.numberOfWaveVectors),
        reciprocalIntegerSquared(forceField.reciprocalIntegerCutOffSquared)
  {
  }

  void expectEqual(const ForceField& actual) const
  {
    EXPECT_EQ(actual.cutOffFrameworkVDW, frameworkVDW);
    EXPECT_EQ(actual.cutOffMoleculeVDW, moleculeVDW);
    EXPECT_EQ(actual.cutOffCoulomb, coulomb);
    EXPECT_EQ(actual.EwaldAlpha, alpha);
    EXPECT_EQ(actual.numberOfWaveVectors.x, waveVectors.x);
    EXPECT_EQ(actual.numberOfWaveVectors.y, waveVectors.y);
    EXPECT_EQ(actual.numberOfWaveVectors.z, waveVectors.z);
    EXPECT_EQ(actual.reciprocalIntegerCutOffSquared, reciprocalIntegerSquared);
  }
};

class WaterBox
{
 public:
  TemporaryDirectory workspace;
  std::unique_ptr<System> system;

  explicit WaterBox(double pressurePa = 1e4)
  {
    workspace.write("mbx.json", settings);
    ForceField forceField({PseudoAtom("Ow", false, 15.9996, -0.84760, 0.0, 8, true),
                           PseudoAtom("Hw", false, 1.0008, 0.42380, 0.0, 1, true)},
                          {VDWParameters(78.19743111, 3.16555789), VDWParameters(0.0, 1.0)},
                          ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false, true);
    forceField.cutOffFrameworkVDWAutomatic = true;
    forceField.cutOffMoleculeVDWAutomatic = true;
    forceField.cutOffCoulombAutomatic = true;
    forceField.automaticEwald = true;
    // At L=18 and a +0.12 log-volume proposal this crosses k=4 -> 5.
    // Thus the added reciprocalIntegerCutOffSquared restoration is observable.
    forceField.EwaldPrecision = 5.357966575133417e-5;

    Component water(forceField, "h2o", 647.096, 22064000.0, 0.3442920843,
                    {Atom(double3(0.0, 0.0, 0.0), -0.84760, 1.0, 0, 0, 0, false, false),
                     Atom(double3(-0.7569503272636612, 0.0, -0.5858822766182950), 0.42380, 1.0, 0, 1, 0, false, false),
                     Atom(double3(0.7569503272636612, 0.0, -0.5858822766182950), 0.42380, 1.0, 0, 1, 0, false, false)},
                    {}, {}, 5, 21);

    system = std::make_unique<System>(forceField, SimulationBox(18.0, 18.0, 18.0), false, 300.0, pressurePa, 1.0,
                                      std::nullopt, std::vector<Component>{water}, std::vector<std::vector<double3>>{},
                                      std::vector<std::size_t>{2}, 5);
    system->useMBX = true;
    system->writeEnergyLog = false;
    system->mbxSettingsFilePath = (workspace.path() / "mbx.json").string();
    const std::array<double3, 2> centers{double3(6.5, 9.0, 9.0), double3(11.5, 9.0, 9.0)};
    for (std::size_t m = 0; m < 2; ++m)
    {
      auto atoms = system->spanOfMolecule(0, m);
      double3 center(0.0);
      double mass = 0.0;
      for (std::size_t a = 0; a < 3; ++a)
      {
        atoms[a].position = centers[m] + system->components[0].atoms[a].position;
        const double atomMass = system->forceField.pseudoAtoms[atoms[a].type].mass;
        center += atomMass * atoms[a].position;
        mass += atomMass;
      }
      system->moleculeData[m].centerOfMassPosition = center / mass;
      system->moleculeData[m].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
    }
    system->forceField.initializeAutomaticCutOff(system->simulationBox);
    system->precomputeTotalRigidEnergy();
    // Give rejection tests a genuine nonempty retained Ewald state to protect.
    Interactions::computeEwaldFourierEnergy(
        system->eik_x, system->eik_y, system->eik_z, system->eik_xy, system->fixedFrameworkStoredEik, system->trialEik,
        system->forceField, system->simulationBox, system->components, system->numberOfMoleculesPerComponent,
        system->spanOfMoleculeAtoms(), system->netChargeFramework);
    Interactions::acceptEwaldMove(system->forceField, system->storedEik, system->trialEik);
  }

  void initializeEnergy() { system->runningEnergies = system->computeTotalEnergies(); }
};

struct RetainedState
{
  SimulationBox box;
  std::vector<Atom> atoms;
  std::vector<Molecule> molecules;
  RunningEnergy energy;
  Cutoffs cutoffs;
  decltype(System::storedEik) eik;
  std::vector<std::size_t> moleculeCounts;
  std::vector<std::size_t> integerCounts;

  explicit RetainedState(const System& system)
      : box(system.simulationBox),
        atoms(system.spanOfMoleculeAtoms().begin(), system.spanOfMoleculeAtoms().end()),
        molecules(system.moleculeData),
        energy(system.runningEnergies),
        cutoffs(system.forceField),
        eik(system.storedEik),
        moleculeCounts(system.numberOfMoleculesPerComponent),
        integerCounts(system.numberOfIntegerMoleculesPerComponent)
  {
  }

  void expectUnchanged(const System& system) const
  {
    EXPECT_EQ(system.simulationBox, box);
    cutoffs.expectEqual(system.forceField);
    EXPECT_EQ(system.storedEik, eik);
    EXPECT_EQ(system.runningEnergies.jsonMC(), energy.jsonMC());
    EXPECT_EQ(system.numberOfMoleculesPerComponent, moleculeCounts);
    EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent, integerCounts);
    ASSERT_EQ(system.spanOfMoleculeAtoms().size(), atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
      const Atom& actual = system.spanOfMoleculeAtoms()[i];
      EXPECT_EQ(actual.position.x, atoms[i].position.x);
      EXPECT_EQ(actual.position.y, atoms[i].position.y);
      EXPECT_EQ(actual.position.z, atoms[i].position.z);
      EXPECT_EQ(actual.charge, atoms[i].charge);
      EXPECT_EQ(actual.scalingVDW, atoms[i].scalingVDW);
      EXPECT_EQ(actual.scalingCoulomb, atoms[i].scalingCoulomb);
      EXPECT_EQ(actual.type, atoms[i].type);
      EXPECT_EQ(actual.moleculeId, atoms[i].moleculeId);
      EXPECT_EQ(actual.componentId, atoms[i].componentId);
      EXPECT_EQ(static_cast<unsigned>(actual.groupId), static_cast<unsigned>(atoms[i].groupId));
      EXPECT_EQ(static_cast<unsigned>(actual.isFractional), static_cast<unsigned>(atoms[i].isFractional));
    }
    ASSERT_EQ(system.moleculeData.size(), molecules.size());
    for (std::size_t i = 0; i < molecules.size(); ++i)
    {
      EXPECT_EQ(nlohmann::json(system.moleculeData[i]), nlohmann::json(molecules[i]));
    }
  }
};

struct Proposal
{
  std::size_t seed;
  double secondUniform;
  double scale;
  SimulationBox box;
  std::vector<Atom> atoms;
  Cutoffs cutoffs;
  double classicalVDW;
};

Proposal prepareProposal(System& system, double logVolumeChange, bool smallAcceptanceUniform = true)
{
  // Search RNG values only: no physical evaluations, rerolled trajectories, or
  // acceptance-based seed selection. This bounded search is library-portable.
  std::size_t chosenSeed = 0;
  double first = 0.0;
  double second = 0.0;
  for (std::size_t seed = 1; seed <= 65536; ++seed)
  {
    RandomNumber probe(seed);
    first = probe.uniform();
    second = probe.uniform();
    const bool direction = logVolumeChange > 0.0 ? first > 0.7 && first < 0.9 : first > 0.1 && first < 0.3;
    const bool acceptanceDraw = smallAcceptanceUniform ? second > 1e-5 && second < 1e-3 : second > 0.2 && second < 0.8;
    if (direction && acceptanceDraw)
    {
      chosenSeed = seed;
      break;
    }
  }
  if (chosenSeed == 0) throw std::runtime_error("No deterministic volume fixture seed in bounded search");

  const double maxChange = logVolumeChange / (2.0 * first - 1.0);
  system.mc_moves_statistics.setMaxChange(volumeMoveType, maxChange);
  const double newVolume = std::exp(std::log(system.simulationBox.volume) + maxChange * (2.0 * first - 1.0));
  const double scale = std::pow(newVolume / system.simulationBox.volume, 1.0 / 3.0);
  SimulationBox newBox = system.simulationBox.scaled(scale);
  std::vector<Atom> newAtoms(system.spanOfMoleculeAtoms().begin(), system.spanOfMoleculeAtoms().end());
  // Independent coordinate construction for this cubic fixture. Do not call
  // the production scaledCenterOfMassPositions helper to form the expectation.
  for (const Molecule& molecule : system.moleculeData)
  {
    double3 center(0.0);
    double mass = 0.0;
    for (std::size_t a = molecule.atomIndex; a < molecule.atomIndex + molecule.numberOfAtoms; ++a)
    {
      const double atomMass = system.forceField.pseudoAtoms[newAtoms[a].type].mass;
      center += atomMass * newAtoms[a].position;
      mass += atomMass;
    }
    const double3 displacement = (scale - 1.0) * (center / mass);
    for (std::size_t a = molecule.atomIndex; a < molecule.atomIndex + molecule.numberOfAtoms; ++a)
      newAtoms[a].position += displacement;
  }
  ForceField trialForceField = system.forceField;
  trialForceField.initializeAutomaticCutOff(newBox);
  const double classicalVDW =
      Interactions::computeInterMolecularEnergy(trialForceField, newBox, newAtoms).moleculeMoleculeVDW;
  return {chosenSeed, second, scale, newBox, std::move(newAtoms), Cutoffs(trialForceField), classicalVDW};
}

void expectCounts(const System& system, double trials, double constructed, double accepted)
{
  const auto& counts = std::get<MoveStatistics<double>>(system.mc_moves_statistics[volumeMoveType]);
  EXPECT_EQ(counts.counts, trials);
  EXPECT_EQ(counts.totalCounts, trials);
  EXPECT_EQ(counts.constructed, constructed);
  EXPECT_EQ(counts.totalConstructed, constructed);
  EXPECT_EQ(counts.accepted, accepted);
  EXPECT_EQ(counts.totalAccepted, accepted);
}

RunningEnergy trialEnergy(const System& system, const Proposal& proposal)
{
  return Interactions::computeMBXEnergySystem(system, system.components, proposal.box, system.framework, {},
                                              proposal.atoms);
}

double logAcceptance(const System& system, const Proposal& proposal, const RunningEnergy& trial)
{
  const double n = static_cast<double>(system.numberOfMolecules());
  return (n + 1.0) * std::log(proposal.box.volume / system.simulationBox.volume) -
         system.beta * (system.pressure * (proposal.box.volume - system.simulationBox.volume) +
                        trial.potentialEnergy() - system.runningEnergies.potentialEnergy());
}

void expectRetainedEnergy(System& system)
{
  const RunningEnergy full = system.computeTotalEnergies();
  ASSERT_TRUE(finiteBits(full.potentialEnergy()));
  EXPECT_NEAR(system.runningEnergies.potentialEnergy(), full.potentialEnergy(), energyTolerance);
  EXPECT_NEAR(system.runningEnergies.mbxEnergy, full.mbxEnergy, energyTolerance);
  EXPECT_NEAR(system.runningEnergies.potentialEnergy(), system.runningEnergies.mbxEnergy, energyTolerance);
}

void exerciseAccepted(System& system, const Proposal& proposal)
{
  ASSERT_TRUE(finiteBits(proposal.classicalVDW));
  ASSERT_LE(proposal.classicalVDW, system.forceField.energyOverlapCriteria);
  const RetainedState before(system);
  const RunningEnergy expected = trialEnergy(system, proposal);
  ASSERT_TRUE(finiteBits(expected.potentialEnergy()));
  ASSERT_LT(std::log(proposal.secondUniform), logAcceptance(system, proposal, expected));
  RandomNumber random(proposal.seed);
  const auto delta = MC_Moves::volumeMove(random, system);
  ASSERT_TRUE(delta.has_value());
  EXPECT_EQ(random.count, 2u);
  // The move returns a difference; the native MC caller performs this addition.
  EXPECT_EQ(system.runningEnergies.jsonMC(), before.energy.jsonMC());
  EXPECT_NEAR(delta->potentialEnergy(), expected.potentialEnergy() - before.energy.potentialEnergy(), energyTolerance);
  system.runningEnergies += *delta;
  EXPECT_NEAR(system.simulationBox.volume, proposal.box.volume, coordinateTolerance);
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j)
    {
      EXPECT_NEAR(system.simulationBox.cell.mm[i][j], proposal.box.cell.mm[i][j], coordinateTolerance);
      EXPECT_NEAR(system.simulationBox.inverseCell.mm[i][j], proposal.box.inverseCell.mm[i][j], coordinateTolerance);
    }
  proposal.cutoffs.expectEqual(system.forceField);
  EXPECT_EQ(system.storedEik, system.trialEik);
  EXPECT_EQ(system.numberOfMoleculesPerComponent, before.moleculeCounts);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent, before.integerCounts);
  for (std::size_t i = 0; i < proposal.atoms.size(); ++i)
  {
    const Atom& actual = system.spanOfMoleculeAtoms()[i];
    expectVectorNear(actual.position, proposal.atoms[i].position);
    EXPECT_EQ(actual.type, before.atoms[i].type);
    EXPECT_EQ(actual.moleculeId, before.atoms[i].moleculeId);
    EXPECT_EQ(actual.componentId, before.atoms[i].componentId);
  }
  for (std::size_t m = 0; m < system.moleculeData.size(); ++m)
  {
    expectVectorNear(system.moleculeData[m].centerOfMassPosition,
                     proposal.scale * before.molecules[m].centerOfMassPosition);
    for (std::size_t a = 0; a < 3; ++a)
      for (std::size_t b = a + 1; b < 3; ++b)
        EXPECT_NEAR(
            (system.spanOfMoleculeAtoms()[3 * m + a].position - system.spanOfMoleculeAtoms()[3 * m + b].position)
                .length(),
            (before.atoms[3 * m + a].position - before.atoms[3 * m + b].position).length(), coordinateTolerance);
  }
  expectRetainedEnergy(system);
}
}  // namespace

TEST(mbx_volume_integration, AcceptedExpansionCommitsMBXEnergyCoordinatesAndAutomaticCutoffs)
{
  WaterBox fixture;
  System& system = *fixture.system;
  fixture.initializeEnergy();
  const Proposal proposal = prepareProposal(system, 0.12);
  ASSERT_NE(proposal.cutoffs.reciprocalIntegerSquared, system.forceField.reciprocalIntegerCutOffSquared);
  exerciseAccepted(system, proposal);
  expectCounts(system, 1, 1, 1);
}

TEST(mbx_volume_integration, AcceptedCompressionCommitsMBXEnergyCoordinatesAndAutomaticCutoffs)
{
  WaterBox fixture;
  System& system = *fixture.system;
  fixture.initializeEnergy();
  const Proposal proposal = prepareProposal(system, -0.12);
  exerciseAccepted(system, proposal);
  expectCounts(system, 1, 1, 1);
}

TEST(mbx_volume_integration, OrdinaryRejectionRestoresAllSixCutoffsAndAllowsNextMove)
{
  // High but finite positive pressure guarantees rejection of this expansion.
  // Pressure is a fixture control, unchanged for the following compression.
  WaterBox fixture(1e10);
  System& system = *fixture.system;
  fixture.initializeEnergy();
  const Proposal proposal = prepareProposal(system, 0.12, false);
  ASSERT_TRUE(finiteBits(proposal.classicalVDW));
  ASSERT_LE(proposal.classicalVDW, system.forceField.energyOverlapCriteria);
  ASSERT_NE(proposal.cutoffs.reciprocalIntegerSquared, system.forceField.reciprocalIntegerCutOffSquared);
  const RunningEnergy expected = trialEnergy(system, proposal);
  ASSERT_GT(std::log(proposal.secondUniform), logAcceptance(system, proposal, expected));
  const RetainedState before(system);
  ASSERT_FALSE(before.eik.empty());
  RandomNumber random(proposal.seed);
  EXPECT_FALSE(MC_Moves::volumeMove(random, system).has_value());
  EXPECT_EQ(random.count, 2u);
  expectCounts(system, 1, 1, 0);
  before.expectUnchanged(system);
  expectRetainedEnergy(system);
  exerciseAccepted(system, prepareProposal(system, -0.12));
  expectCounts(system, 2, 2, 1);
}

TEST(mbx_volume_integration, FiniteOverlapScreenRejectsBeforeConstructionAndAllowsNextMove)
{
  WaterBox fixture;
  System& system = *fixture.system;
  fixture.initializeEnergy();
  // Compress the O--O separation from 5 A to 1 A. Only the conventional
  // pre-screen sees this trial; no trial MBX energy is requested by this test.
  const Proposal proposal = prepareProposal(system, 3.0 * std::log(0.2));
  ASSERT_TRUE(finiteBits(proposal.classicalVDW));
  ASSERT_GT(proposal.classicalVDW, system.forceField.energyOverlapCriteria);
  ASSERT_NE(proposal.cutoffs.moleculeVDW, system.forceField.cutOffMoleculeVDW);
  const RetainedState before(system);
  RandomNumber random(proposal.seed);
  EXPECT_FALSE(MC_Moves::volumeMove(random, system).has_value());
  EXPECT_EQ(random.count, 1u);
  expectCounts(system, 1, 0, 0);
  before.expectUnchanged(system);
  expectRetainedEnergy(system);
  exerciseAccepted(system, prepareProposal(system, 0.12));
  expectCounts(system, 2, 1, 1);
}

TEST(mbx_volume_integration, NonfiniteConventionalScreenRejectsValidMBXGeometryAndAllowsNextMove)
{
  // Check the observer before using it as a fault-injection precondition.
  const double finiteControl = 1.0;
  const double infinityControl = std::bit_cast<double>(UINT64_C(0x7ff0000000000000));
  const double nanControl = std::bit_cast<double>(UINT64_C(0x7ff8000000000042));
  ASSERT_TRUE(finiteBits(finiteControl));
  ASSERT_FALSE(finiteBits(infinityControl));
  ASSERT_FALSE(finiteBits(nanControl));
  WaterBox fixture;
  System& system = *fixture.system;
  fixture.initializeEnergy();
  // Numerical-fault injection is confined to the conventional O--O screen.
  // The old/trial MBX geometry, MBX settings and MBX model remain valid.
  const double originalEpsilon = system.forceField.data[0].parameters.x;
  system.forceField.data[0].parameters.x = nanControl;
  ASSERT_FALSE(finiteBits(system.forceField.data[0].parameters.x));
  const Proposal proposal = prepareProposal(system, 0.12);
  ASSERT_FALSE(finiteBits(proposal.classicalVDW));
  ASSERT_NE(proposal.cutoffs.reciprocalIntegerSquared, system.forceField.reciprocalIntegerCutOffSquared);
  const RetainedState before(system);
  RandomNumber random(proposal.seed);
  const auto delta = MC_Moves::volumeMove(random, system);
  EXPECT_FALSE(delta.has_value()) << "Nonfinite pre-trial MBX screen was bypassed";
  EXPECT_EQ(random.count, 1u);
  expectCounts(system, 1, 0, 0);
  before.expectUnchanged(system);
  // A bypass is retained as a failing result; do not hide it by continuing
  // from an unexpectedly committed state or changing the production engine.
  if (delta.has_value()) return;
  expectRetainedEnergy(system);
  system.forceField.data[0].parameters.x = originalEpsilon;
  exerciseAccepted(system, prepareProposal(system, 0.12));
  expectCounts(system, 2, 1, 1);
}

TEST(mbx_volume_integration, DiagnosticSwitchPreservesProposalAcceptanceStateAndRNG)
{
  WaterBox quiet;
  WaterBox logged;
  quiet.initializeEnergy();
  // The two fixtures have identical fixed physical inputs. Reuse this known
  // old-state energy to avoid a redundant MBX evaluation.
  logged.system->runningEnergies = quiet.system->runningEnergies;
  const Proposal a = prepareProposal(*quiet.system, 0.12);
  const Proposal b = prepareProposal(*logged.system, 0.12);
  ASSERT_EQ(a.seed, b.seed);
  const auto quietPath = quiet.workspace.path() / "energy_terms.csv";
  const auto loggedPath = logged.workspace.path() / "energy_terms.csv";
  quiet.system->configureEnergyTermsLog(quietPath.string(), false);
  logged.system->writeEnergyLog = true;
  logged.system->configureEnergyTermsLog(loggedPath.string(), false);
  RandomNumber rngQuiet(a.seed);
  RandomNumber rngLogged(b.seed);
  const auto deltaQuiet = MC_Moves::volumeMove(rngQuiet, *quiet.system);
  const auto deltaLogged = MC_Moves::volumeMove(rngLogged, *logged.system);
  ASSERT_TRUE(deltaQuiet.has_value());
  ASSERT_TRUE(deltaLogged.has_value());
  EXPECT_EQ(rngQuiet, rngLogged);
  EXPECT_EQ(deltaQuiet->jsonMC(), deltaLogged->jsonMC());
  quiet.system->runningEnergies += *deltaQuiet;
  logged.system->runningEnergies += *deltaLogged;
  RetainedState(*quiet.system).expectUnchanged(*logged.system);
  EXPECT_FALSE(std::filesystem::exists(quietPath));
  ASSERT_TRUE(std::filesystem::exists(loggedPath));
  std::ifstream input(loggedPath);
  const std::string contents{std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
  EXPECT_NE(contents.find("\nvolume,all,2,"), std::string::npos);
  expectCounts(*quiet.system, 1, 1, 1);
  expectCounts(*logged.system, 1, 1, 1);
  expectRetainedEnergy(*quiet.system);
}
