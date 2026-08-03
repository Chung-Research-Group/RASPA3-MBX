#include <gtest/gtest.h>

#include "../raspa3-tests/input_reader_fixtures.hpp"
#include "../raspa3-tests/molecule_fixtures.hpp"
#include "../test_support.hpp"

import std;

import int3;
import double3;
import json;
import units;
import atom;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import running_energy;
import interactions_mbx;
import input_reader;
import monte_carlo;
import monte_carlo_transition_matrix;
import simulation_schedule;
import energy_evaluation;
import archive;
import randomnumbers;

namespace
{
constexpr double staticMbxEnergy = 444.5797675540443;

constexpr double oldRigidMoveMbxEnergy = 444.1757421121683;
constexpr double newRigidMoveMbxEnergy = 435.5616147453709;

constexpr double directRigidMoveDifference = -8.614127366797334;
constexpr double firstRigidMoveDifference = -3.20172544713563;
constexpr double secondRigidMoveDifference = -5.412401968369863;
constexpr double sequentialRigidMoveDifference = -8.614127415505493;

constexpr double mbxEnergyTolerance = 1e-6;

class ScopedCurrentPath
{
 public:
  explicit ScopedCurrentPath(const std::filesystem::path& path) : originalPath(std::filesystem::current_path())
  {
    std::filesystem::current_path(path);
  }

  ~ScopedCurrentPath() { std::filesystem::current_path(originalPath); }

 private:
  std::filesystem::path originalPath;
};

TemporaryDirectory makeInputReaderWorkspace()
{
  TemporaryDirectory directory;
  directory.write("force_field.json", input_reader_fixtures::kBoxForceFieldJson);
  directory.write("methane.json", molecule_fixtures::kMethaneJson);
  return directory;
}

nlohmann::json energyTermsInput()
{
  nlohmann::json input = nlohmann::json::parse(input_reader_fixtures::kNvtSimulationJson);
  input["SimulationType"] = "MonteCarlo";
  input["Components"][0]["CreateNumberOfMolecules"] = 1;
  return input;
}

const nlohmann::json mbxSettings = {{"Note", "This is an MBX v1.0 configuration file"},
                                    {"MBX",
                                     {{"box", nlohmann::json::array()},
                                      {"twobody_cutoff", 9.0},
                                      {"threebody_cutoff", 5.0},
                                      {"dipole_tolerance", 1E-9},
                                      {"dipole_method", "cg"},
                                      {"alpha_ewald_elec", 0.6},
                                      {"grid_density_elec", 2.5},
                                      {"spline_order_elec", 6},
                                      {"alpha_ewald_disp", 0.60},
                                      {"grid_density_disp", 2.5},
                                      {"spline_order_disp", 6},
                                      {"ignore_2b_poly", nlohmann::json::array()},
                                      {"ignore_3b_poly", nlohmann::json::array()}}}};

void setTwoCO2MFIPositions(std::span<Atom> atomData)
{
  atomData[0].position = double3(10.011, 4.97475 + 2.0, 0.0);
  atomData[1].position = double3(10.011, 4.97475 + 2.0, 1.149);
  atomData[2].position = double3(10.011, 4.97475 + 2.0, -1.149);

  atomData[3].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomData[4].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomData[5].position = double3(10.011, 4.97475 - 2.0, -1.149);
}

double3 centroid(std::span<const Atom> moleculeAtoms) noexcept
{
  double3 center(0.0);

  for (const Atom& atom : moleculeAtoms)
  {
    center += atom.position;
  }

  return center / static_cast<double>(moleculeAtoms.size());
}

void translateAndRotateMolecule(std::span<Atom> moleculeAtoms, const double3& displacement,
                                const double3& rotationVector) noexcept
{
  const double3 center = centroid(moleculeAtoms);
  const double theta = rotationVector.length();

  if (theta < 1e-14)
  {
    for (Atom& atom : moleculeAtoms)
    {
      atom.position += displacement;
    }
    return;
  }

  const double3 axis = rotationVector / theta;

  for (Atom& atom : moleculeAtoms)
  {
    const double3 relative = atom.position - center;
    atom.position = center + axis.rotateAroundAxis(relative, theta) + displacement;
  }
}

std::vector<Atom> moleculeSlice(const std::vector<Atom>& moleculeAtoms, std::size_t moleculeId, std::size_t numAtoms)
{
  const std::size_t moleculeOffset = moleculeId * numAtoms;
  const auto first = moleculeAtoms.begin() + static_cast<std::ptrdiff_t>(moleculeOffset);
  return std::vector<Atom>(first, first + static_cast<std::ptrdiff_t>(numAtoms));
}

void setupFrameworkPermanentEnergy(System& system, std::span<const Atom> frameworkAtoms)
{
  system.elecPermFrameworkMBX =
      Interactions::computeFrameworkElecPermMBXEnergy(system, system.simulationBox, system.framework, frameworkAtoms);
}

std::filesystem::path makeTemporaryMbxSettingsFilePath()
{
  std::random_device randomDevice;
  const auto timeStamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  const std::uint64_t randomSuffix =
      (static_cast<std::uint64_t>(randomDevice()) << 32U) ^ static_cast<std::uint64_t>(randomDevice());

  return std::filesystem::temp_directory_path() /
         ("_json_" + std::to_string(timeStamp) + "_" + std::to_string(randomSuffix) + ".mbx");
}

class TemporaryMbxSettingsFile
{
 public:
  TemporaryMbxSettingsFile() : filePath(makeTemporaryMbxSettingsFilePath())
  {
    std::ofstream file(filePath);
    if (!file)
    {
      throw std::runtime_error("Could not write MBX settings file: " + filePath.string());
    }

    file << mbxSettings.dump(3) << '\n';
  }

  ~TemporaryMbxSettingsFile()
  {
    std::error_code errorCode;
    std::filesystem::remove(filePath, errorCode);
  }

  const std::filesystem::path& path() const noexcept { return filePath; }

 private:
  std::filesystem::path filePath;
};

const TemporaryMbxSettingsFile& mbxSettingsFile()
{
  static const TemporaryMbxSettingsFile settingsFile;
  return settingsFile;
}

}  // namespace

TEST(mbx_static_energy, Test_2_CO2_in_MFI_2x2x2_system_matches_direct_MBX_energy)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component c = Component::makeCO2(forceField, 0, true);

  // mbx defines component names in lower case and takes COO type order
  c.name = "co2";
  std::swap(c.definedAtoms[0], c.definedAtoms[1]);
  std::vector<Component> components{c};

  Framework f = Framework::makeMFI(forceField, int3(2, 2, 2));

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, components, {}, {2}, 5);
  system.useMBX = true;
  system.mbxSettingsFilePath = mbxSettingsFile().path().string();

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  setTwoCO2MFIPositions(atomData);
  setupFrameworkPermanentEnergy(system, frameworkAtoms);

  RunningEnergy energySystem = Interactions::computeMBXEnergySystem(system, components, system.simulationBox,
                                                                    system.framework, frameworkAtoms, atomData);

  std::array<double, 7> reportedTerms{};
  const RunningEnergy evaluatedEnergy = system.computeTotalEnergies(reportedTerms);

  RunningEnergy energyDirect = Interactions::computeMBXEnergy(system, components, system.simulationBox,
                                                              system.framework, 0, frameworkAtoms, atomData, {}, false);

  EXPECT_NEAR(staticMbxEnergy, energySystem.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(staticMbxEnergy, energyDirect.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(energySystem.mbxEnergy, energyDirect.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(energySystem.mbxEnergy, evaluatedEnergy.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(evaluatedEnergy.mbxEnergy * Units::EnergyToKCalPerMol,
              std::accumulate(reportedTerms.begin() + 1, reportedTerms.end(), 0.0), mbxEnergyTolerance);
}

TEST(mbx_static_energy, Test_2_CO2_in_MFI_2x2x2_selected_molecule_reinsert_matches_system_energy)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component c = Component::makeCO2(forceField, 0, true);

  // mbx defines component names in lower case and takes COO type order
  c.name = "co2";
  std::swap(c.definedAtoms[0], c.definedAtoms[1]);
  std::vector<Component> components{c};

  Framework f = Framework::makeMFI(forceField, int3(2, 2, 2));

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, components, {}, {2}, 5);
  system.useMBX = true;
  system.mbxSettingsFilePath = mbxSettingsFile().path().string();

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  setTwoCO2MFIPositions(atomData);
  setupFrameworkPermanentEnergy(system, frameworkAtoms);

  const std::size_t numAtoms = components[0].atoms.size();
  std::span<const Atom> selectedMoleculeAtoms(atomData.data(), numAtoms);

  RunningEnergy energySystem = Interactions::computeMBXEnergySystem(system, components, system.simulationBox,
                                                                    system.framework, frameworkAtoms, atomData);

  RunningEnergy energySelectedIncluded =
      Interactions::computeMBXEnergy(system, components, system.simulationBox, system.framework, 0, frameworkAtoms,
                                     atomData, selectedMoleculeAtoms, true);

  EXPECT_NEAR(staticMbxEnergy, energySystem.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(staticMbxEnergy, energySelectedIncluded.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(energySystem.mbxEnergy, energySelectedIncluded.mbxEnergy, mbxEnergyTolerance);
}

TEST(mbx_static_energy, Test_nonfirst_component_global_molecule_id_replacement_and_deletion)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component firstComponent = Component::makeCO2(forceField, 0, true);
  Component secondComponent = Component::makeCO2(forceField, 1, true);

  // Use the same MBX monomer model in two distinct RASPA components. RASPA assigns the second molecule
  // global molecule ID 1, although it is component-local molecule 0 in the second component.
  for (Component* component : {&firstComponent, &secondComponent})
  {
    component->name = "co2";
    std::swap(component->definedAtoms[0], component->definedAtoms[1]);
  }
  std::vector<Component> components{firstComponent, secondComponent};

  Framework framework = Framework::makeMFI(forceField, int3(2, 2, 2));
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, components, {}, {1, 1}, 5);
  system.useMBX = true;
  system.mbxSettingsFilePath = mbxSettingsFile().path().string();

  std::span<Atom> moleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  setTwoCO2MFIPositions(moleculeAtoms);
  setupFrameworkPermanentEnergy(system, frameworkAtoms);

  const std::size_t atomsPerMolecule = components[0].atoms.size();
  ASSERT_EQ(2 * atomsPerMolecule, moleculeAtoms.size());
  ASSERT_EQ(0uz, static_cast<std::size_t>(moleculeAtoms.front().moleculeId));
  ASSERT_EQ(1uz, static_cast<std::size_t>(moleculeAtoms[atomsPerMolecule].moleculeId));
  ASSERT_EQ(1uz, static_cast<std::size_t>(moleculeAtoms[atomsPerMolecule].componentId));

  std::vector<Atom> oldSelectedMolecule(moleculeAtoms.begin() + static_cast<std::ptrdiff_t>(atomsPerMolecule),
                                        moleculeAtoms.end());
  std::vector<Atom> newSelectedMolecule = oldSelectedMolecule;
  translateAndRotateMolecule(newSelectedMolecule, double3(0.037, -0.019, 0.011), double3(0.009, -0.006, 0.004));

  const RunningEnergy oldFullEnergy = Interactions::computeMBXEnergySystem(
      system, components, system.simulationBox, system.framework, frameworkAtoms, moleculeAtoms);

  const RunningEnergy reinsertedEnergy =
      Interactions::computeMBXEnergy(system, components, system.simulationBox, system.framework, 1, frameworkAtoms,
                                     moleculeAtoms, oldSelectedMolecule, true);
  EXPECT_NEAR(oldFullEnergy.mbxEnergy, reinsertedEnergy.mbxEnergy, mbxEnergyTolerance);

  std::vector<Atom> replacedMoleculeAtoms(moleculeAtoms.begin(), moleculeAtoms.end());
  std::copy(newSelectedMolecule.begin(), newSelectedMolecule.end(),
            replacedMoleculeAtoms.begin() + static_cast<std::ptrdiff_t>(atomsPerMolecule));
  const RunningEnergy newFullEnergy = Interactions::computeMBXEnergySystem(
      system, components, system.simulationBox, system.framework, frameworkAtoms, replacedMoleculeAtoms);

  const RunningEnergy replacementEnergy =
      Interactions::computeMBXEnergy(system, components, system.simulationBox, system.framework, 1, frameworkAtoms,
                                     moleculeAtoms, newSelectedMolecule, true);
  const RunningEnergy replacementDifference =
      Interactions::computeMBXEnergyDifference(system, components, system.simulationBox, system.framework, 1,
                                               frameworkAtoms, moleculeAtoms, newSelectedMolecule, oldSelectedMolecule);

  EXPECT_NEAR(newFullEnergy.mbxEnergy, replacementEnergy.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(newFullEnergy.mbxEnergy - oldFullEnergy.mbxEnergy, replacementDifference.mbxEnergy, mbxEnergyTolerance);

  const RunningEnergy deletedEnergy =
      Interactions::computeMBXEnergy(system, components, system.simulationBox, system.framework, 1, frameworkAtoms,
                                     moleculeAtoms, oldSelectedMolecule, false);

  System deletionReference = system;
  std::vector<Atom> atomsToDelete(oldSelectedMolecule);
  deletionReference.deleteMolecule(1, 0, std::span<Atom>(atomsToDelete));
  const RunningEnergy directDeletedEnergy = Interactions::computeMBXEnergySystem(
      deletionReference, components, deletionReference.simulationBox, deletionReference.framework,
      deletionReference.spanOfFrameworkAtoms(), deletionReference.spanOfMoleculeAtoms());

  EXPECT_NEAR(directDeletedEnergy.mbxEnergy, deletedEnergy.mbxEnergy, mbxEnergyTolerance);
  EXPECT_GT(std::abs(oldFullEnergy.mbxEnergy - deletedEnergy.mbxEnergy), 1e-14);
}

TEST(mbx_energy_log, Test_print_energy_terms_input_and_legacy_alias)
{
  TemporaryDirectory workspace = makeInputReaderWorkspace();
  ScopedCurrentPath currentPath(workspace.path());

  nlohmann::json input = energyTermsInput();
  input["Systems"][0]["PrintEnergyTerms"] = false;
  workspace.write("simulation.json", input.dump(2));
  InputReader disabledReader("simulation.json");
  ASSERT_EQ(disabledReader.systems.size(), 1uz);
  EXPECT_FALSE(disabledReader.systems[0].writeEnergyLog);

  input["Systems"][0]["PrintEnergyTerms"] = true;
  workspace.write("simulation.json", input.dump(2));
  InputReader enabledReader("simulation.json");
  EXPECT_TRUE(enabledReader.systems[0].writeEnergyLog);

  input["Systems"][0].erase("PrintEnergyTerms");
  input["Systems"][0]["WriteEnergyLog"] = false;
  workspace.write("simulation.json", input.dump(2));
  InputReader legacyReader("simulation.json");
  EXPECT_FALSE(legacyReader.systems[0].writeEnergyLog);

  input["Systems"][0]["PrintEnergyTerms"] = false;
  workspace.write("simulation.json", input.dump(2));
  InputReader matchingAliasesReader("simulation.json");
  EXPECT_FALSE(matchingAliasesReader.systems[0].writeEnergyLog);

  input["Systems"][0]["PrintEnergyTerms"] = true;
  workspace.write("simulation.json", input.dump(2));
  EXPECT_THROW(InputReader conflictingAliasesReader("simulation.json"), std::runtime_error);

  input["Systems"][0].erase("WriteEnergyLog");
  input["Systems"][0]["PrintEnergyTerms"] = "yes";
  workspace.write("simulation.json", input.dump(2));
  EXPECT_THROW(InputReader invalidTypeReader("simulation.json"), std::runtime_error);

  input["Systems"][0].erase("PrintEnergyTerms");
  workspace.write("simulation.json", input.dump(2));
  InputReader defaultReader("simulation.json");
  EXPECT_TRUE(defaultReader.systems[0].writeEnergyLog);
}

TEST(run_control, Test_write_restart_every_parser_and_driver_propagation)
{
  TemporaryDirectory workspace = makeInputReaderWorkspace();
  ScopedCurrentPath currentPath(workspace.path());

  nlohmann::json input = energyTermsInput();
  input["WriteRestartEvery"] = 17;
  workspace.write("simulation.json", input.dump(2));
  InputReader configuredReader("simulation.json");
  EXPECT_EQ(17uz, configuredReader.writeRestartEvery);
  MonteCarlo configuredMonteCarlo(configuredReader);
  EXPECT_EQ(17uz, configuredMonteCarlo.writeRestartEvery);

  SimulationSchedule schedule;
  schedule.writeRestartEvery = 23;
  MonteCarlo scheduledMonteCarlo(schedule, {}, 42uz, 5, false);
  EXPECT_EQ(23uz, scheduledMonteCarlo.writeRestartEvery);

  TemporaryDirectory archiveWorkspace;
  const std::filesystem::path archivePath = archiveWorkspace.path() / "restart_data.bin";
  writeBinaryRestartFile(scheduledMonteCarlo, archivePath.string());
  MonteCarlo restoredMonteCarlo;
  readBinaryRestartFile(restoredMonteCarlo, archivePath.string());
  EXPECT_EQ(23uz, restoredMonteCarlo.writeRestartEvery);

  input["WriteRestartEvery"] = 0;
  workspace.write("simulation.json", input.dump(2));
  InputReader disabledReader("simulation.json");
  EXPECT_EQ(0uz, disabledReader.writeRestartEvery);

  input["WriteRestartEvery"] = -1;
  workspace.write("simulation.json", input.dump(2));
  EXPECT_THROW(InputReader negativeReader("simulation.json"), std::runtime_error);

  input["WriteRestartEvery"] = "often";
  workspace.write("simulation.json", input.dump(2));
  EXPECT_THROW(InputReader invalidTypeReader("simulation.json"), std::runtime_error);

  for (const nlohmann::json& invalidValue : {nlohmann::json(1.5), nlohmann::json(true), nlohmann::json(nullptr)})
  {
    input["WriteRestartEvery"] = invalidValue;
    workspace.write("simulation.json", input.dump(2));
    EXPECT_THROW(InputReader invalidReader("simulation.json"), std::runtime_error);
  }

  input.erase("WriteRestartEvery");
  workspace.write("simulation.json", input.dump(2));
  InputReader defaultReader("simulation.json");
  EXPECT_EQ(5000uz, defaultReader.writeRestartEvery);
}

TEST(energy_evaluation, Test_restart_json_is_evaluated_once_without_monte_carlo_moves)
{
  TemporaryDirectory workspace;
  ScopedCurrentPath currentPath(workspace.path());
  workspace.write("force_field.json", input_reader_fixtures::kBoxForceFieldJson);
  workspace.write("methane.json", molecule_fixtures::kMethaneJson);

  nlohmann::json input = energyTermsInput();
  input["SimulationType"] = "EnergyEvaluation";
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  input["Systems"][0]["BoxLengths"] = {25.0, 25.0, 25.0};
  input["Systems"][0]["RestartFileName"] = "snapshot";
  workspace.write("case/simulation.json", input.dump(2));

  const nlohmann::json snapshot{
      {"SimulationBox",
       {{"length-a", 30.0},
        {"length-b", 30.0},
        {"length-c", 30.0},
        {"angle-alpha", 90.0},
        {"angle-beta", 90.0},
        {"angle-gamma", 90.0}}},
      {"methane", {{5.0, 5.0, 5.0}, {9.0, 5.0, 5.0}}},
  };
  workspace.write("case/snapshot.json", snapshot.dump(2));

  InputReader reader((workspace.path() / "case/simulation.json").string());
  ASSERT_EQ(InputReader::SimulationType::EnergyEvaluation, reader.simulationType);
  ASSERT_EQ(1uz, reader.systems.size());
  ASSERT_EQ(2uz, reader.systems[0].numberOfIntegerMoleculesPerComponent.at(0));
  ASSERT_EQ(2uz, reader.systems[0].spanOfMoleculeAtoms().size());
  EXPECT_DOUBLE_EQ(30.0, reader.systems[0].simulationBox.lengthA);
  const std::vector<double3> positionsBefore =
      reader.systems[0].spanOfMoleculeAtoms() | std::views::transform(&Atom::position) | std::ranges::to<std::vector>();

  testing::internal::CaptureStdout();
  EnergyEvaluation evaluation(reader);
  evaluation.run();
  const std::string standardOutput = testing::internal::GetCapturedStdout();
  EXPECT_NE(std::string::npos, standardOutput.find("Total potential energy"));

  const std::vector<double3> positionsAfter =
      reader.systems[0].spanOfMoleculeAtoms() | std::views::transform(&Atom::position) | std::ranges::to<std::vector>();
  EXPECT_EQ(positionsBefore, positionsAfter);

  const std::filesystem::path outputPath = workspace.path() / "output/energy_evaluation.json";
  ASSERT_TRUE(std::filesystem::is_regular_file(outputPath));
  std::ifstream outputFile(outputPath);
  const nlohmann::json output = nlohmann::json::parse(outputFile);
  ASSERT_EQ(1uz, output["systems"].size());
  EXPECT_EQ("Classical", output["systems"][0]["model"]);
  EXPECT_DOUBLE_EQ(30.0, output["systems"][0]["simulationBox"]["boxLengths"][0].get<double>());
  EXPECT_NEAR(-144.239545749195855, output["systems"][0]["energy"]["Total potential energy [K]"].get<double>(), 1.0e-9);
  EXPECT_NEAR(-144.239545749195855, output["systems"][0]["energy"]["molecule-molecule VDW [K]"].get<double>(), 1.0e-9);
  EXPECT_FALSE(std::filesystem::exists(workspace.path() / "output/energy_terms.s0.csv"));

  // The evaluator owns its System copies, so constructing it from a temporary reader is safe.
  testing::internal::CaptureStdout();
  EnergyEvaluation temporaryReaderEvaluation{InputReader((workspace.path() / "case/simulation.json").string())};
  EXPECT_NO_THROW(temporaryReaderEvaluation.run());
  testing::internal::GetCapturedStdout();
}

TEST(energy_evaluation, Test_restart_snapshot_rejects_generated_extra_molecules_and_malformed_json)
{
  TemporaryDirectory workspace = makeInputReaderWorkspace();
  ScopedCurrentPath currentPath(workspace.path());

  nlohmann::json input = energyTermsInput();
  input["SimulationType"] = "EnergyEvaluation";
  input["Systems"][0]["RestartFileName"] = "snapshot.json";
  workspace.write("snapshot.json", R"({"methane": [[5.0, 5.0, 5.0]]})");
  workspace.write("simulation.json", input.dump(2));
  EXPECT_THROW(InputReader generatedMoleculeReader("simulation.json"), std::runtime_error);

  input["Components"][0]["CreateNumberOfMolecules"] = 0;

  workspace.write("snapshot.json", R"({"methnae": [[5.0, 5.0, 5.0]]})");
  workspace.write("simulation.json", input.dump(2));
  EXPECT_THROW(InputReader unknownComponentReader("simulation.json"), std::runtime_error);

  nlohmann::json duplicateComponentInput = input;
  duplicateComponentInput["Components"].push_back(duplicateComponentInput["Components"][0]);
  duplicateComponentInput["Components"][1]["Name"] = "METHANE";
  workspace.write("snapshot.json", R"({"methane": [[5.0, 5.0, 5.0]]})");
  workspace.write("simulation.json", duplicateComponentInput.dump(2));
  EXPECT_THROW(InputReader duplicateComponentReader("simulation.json"), std::runtime_error);

  workspace.write("snapshot.json", "{not valid json");
  workspace.write("simulation.json", input.dump(2));
  EXPECT_THROW(InputReader malformedReader("simulation.json"), std::runtime_error);

  workspace.write("snapshot.json", R"({"methane": [[5.0, 5.0, 5.0]]})");
  InputReader fractionalReader("simulation.json");
  fractionalReader.systems[0].components[0].hasFractionalMolecule = true;
  fractionalReader.systems[0].numberOfFractionalMoleculesPerComponent[0] = 1;
  EnergyEvaluation fractionalEvaluation(fractionalReader);
  EXPECT_THROW(fractionalEvaluation.run(), std::runtime_error);

  input["RestartFromBinaryFile"] = true;
  workspace.write("simulation.json", input.dump(2));
  InputReader binaryReader("simulation.json");
  EnergyEvaluation binaryEvaluation(binaryReader);
  EXPECT_THROW(binaryEvaluation.run(), std::runtime_error);

  input["RestartFromBinaryFile"] = "true";
  workspace.write("simulation.json", input.dump(2));
  EXPECT_THROW(InputReader invalidBinaryFlagReader("simulation.json"), std::runtime_error);
}

TEST(energy_evaluation, Test_fractional_configuration_is_rejected_before_system_construction)
{
  TemporaryDirectory workspace = makeInputReaderWorkspace();
  ScopedCurrentPath currentPath(workspace.path());

  nlohmann::json input = energyTermsInput();
  input["SimulationType"] = "EnergyEvaluation";
  input["Components"][0]["CreateNumberOfMolecules"] = 0;
  input["Systems"][0]["RestartFileName"] = "snapshot.json";
  workspace.write("snapshot.json", R"({"methane": [[5.0, 5.0, 5.0]]})");

  const std::array componentKeywords{"CFCMC_SwapProbability",
                                     "CFCMC_CBMC_SwapProbability",
                                     "CFCMC_WidomProbability",
                                     "CFCMC_CBMC_WidomProbability",
                                     "GibbsSwapCFCMCProbability",
                                     "GibbsSwapCBCFCMCProbability",
                                     "GibbsConventionalCFCMCProbability",
                                     "GibbsConventionalCBCFCMCProbability",
                                     "CFCMC_PairSwapProbability",
                                     "CFCMC_CBMC_PairSwapProbability",
                                     "CFCMC_GroupSwapProbability",
                                     "CFCMC_CBMC_GroupSwapProbability"};
  for (const std::string_view keyword : componentKeywords)
  {
    SCOPED_TRACE(keyword);
    nlohmann::json fractionalInput = input;
    fractionalInput["Components"][0][keyword] = 1.0;
    workspace.write("simulation.json", fractionalInput.dump(2));
    EXPECT_THROW(InputReader fractionalMoveReader("simulation.json"), std::runtime_error);
  }

  nlohmann::json fixedLambdaInput = input;
  fixedLambdaInput["Components"][0]["LambdaBinIndex"] = 0;
  workspace.write("simulation.json", fixedLambdaInput.dump(2));
  EXPECT_THROW(InputReader fixedLambdaReader("simulation.json"), std::runtime_error);

  const std::array reactionKeywords{"ReactionConventionalCFCMCProbability", "ReactionConventionalCBCFCMCProbability",
                                    "ReactionCFCMCProbability", "ReactionCBCFCMCProbability"};
  for (const std::string_view keyword : reactionKeywords)
  {
    SCOPED_TRACE(keyword);
    nlohmann::json fractionalInput = input;
    fractionalInput["Systems"][0][keyword] = 1.0;
    workspace.write("simulation.json", fractionalInput.dump(2));
    EXPECT_THROW(InputReader fractionalReactionReader("simulation.json"), std::runtime_error);
  }
}

TEST(mbx_energy_log, Test_component_rows_use_component_loading_and_volume_rows_use_total_loading)
{
  TemporaryDirectory temporaryDirectory;
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component firstComponent = Component::makeCO2(forceField, 0, true);
  Component secondComponent = Component::makeCO2(forceField, 1, true);
  std::vector<Component> components{firstComponent, secondComponent};
  Framework framework = Framework::makeMFI(forceField, int3(1, 1, 1));

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, components, {}, {1, 1}, 5);
  system.useMBX = true;
  system.writeEnergyLog = true;
  const std::filesystem::path logPath = temporaryDirectory.path() / "output" / "energy_terms.s0.csv";
  system.configureEnergyTermsLog(logPath.string(), false);

  const RunningEnergy totalEnergy;
  const std::array<double, 7> mbxTerms{};

  testing::internal::CaptureStderr();
  system.writeAcceptedEnergyLog("translation", 1, totalEnergy, mbxTerms, 0.0, 1.0);
  system.writeAcceptedEnergyLog("volume", std::numeric_limits<std::size_t>::max(), totalEnergy, mbxTerms, 0.0, 1.0);
  EXPECT_TRUE(testing::internal::GetCapturedStderr().empty());

  std::ifstream input(logPath);
  const std::string contents{std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
  const std::string header =
      "type,component,N,total,hg_VDW,hg_tail,mbx_tot,E2b,E3b,E4b,Edisp,Eelec_perm,Eelec_ind,E_diff,Pacc\n";
  EXPECT_TRUE(contents.starts_with(header));
  EXPECT_NE(std::string::npos, contents.find("\ntranslation,1,1,"));
  EXPECT_NE(std::string::npos, contents.find("\nvolume,all,2,"));

  // A binary-restart resume appends to the existing file and must not duplicate its CSV header.
  system.configureEnergyTermsLog(logPath.string(), true);
  system.writeAcceptedEnergyLog("rotation", 0, totalEnergy, mbxTerms, 0.0, 1.0);
  std::ifstream appendedInput(logPath);
  const std::string appended{std::istreambuf_iterator<char>(appendedInput), std::istreambuf_iterator<char>()};
  EXPECT_EQ(appended.find(header), appended.rfind(header));
  EXPECT_NE(std::string::npos, appended.find("\nrotation,0,1,"));
}

TEST(mbx_energy_log, Test_disabled_energy_terms_create_no_file_or_stderr_output)
{
  TemporaryDirectory temporaryDirectory;
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component component = Component::makeCO2(forceField, 0, true);
  Framework framework = Framework::makeMFI(forceField, int3(1, 1, 1));
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {component}, {}, {1}, 5);
  system.useMBX = true;
  system.writeEnergyLog = false;

  const std::filesystem::path logPath = temporaryDirectory.path() / "output" / "energy_terms.s0.csv";
  system.configureEnergyTermsLog(logPath.string(), false);

  testing::internal::CaptureStderr();
  system.writeAcceptedEnergyLog("translation", 0, RunningEnergy{}, {}, 0.0, 1.0);
  EXPECT_TRUE(testing::internal::GetCapturedStderr().empty());
  EXPECT_FALSE(std::filesystem::exists(logPath));
}

TEST(mbx_energy_log, Test_transition_matrix_fresh_truncates_and_resume_appends_energy_log)
{
  TemporaryDirectory temporaryDirectory;
  ScopedCurrentPath currentPath(temporaryDirectory.path());
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component component = Component::makeCO2(forceField, 0, true);
  Framework framework = Framework::makeMFI(forceField, int3(1, 1, 1));
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {framework}, {component}, {}, {1}, 5);
  system.writeEnergyLog = true;

  temporaryDirectory.write("output/energy_terms.s0.csv", "stale-row\n");
  SimulationSchedule schedule;
  {
    std::vector<System> systems{system};
    RandomNumber random(std::optional<std::size_t>{42uz});
    MonteCarloTransitionMatrix transitionMatrix(schedule, systems, random, 5);
    transitionMatrix.createOutputFiles();
  }

  const std::filesystem::path logPath = temporaryDirectory.path() / "output/energy_terms.s0.csv";
  {
    std::ifstream input(logPath);
    const std::string contents{std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
    EXPECT_EQ(std::string::npos, contents.find("stale-row"));
    EXPECT_TRUE(contents.starts_with("type,component,N,total,"));
  }

  {
    std::ofstream marker(logPath, std::ios::app);
    marker << "resume-tail\n";
  }
  {
    std::vector<System> systems{system};
    RandomNumber random(std::optional<std::size_t>{42uz});
    MonteCarloTransitionMatrix transitionMatrix(schedule, systems, random, 5);
    transitionMatrix.createOutputFiles(true);
  }

  std::ifstream resumedInput(logPath);
  const std::string resumedContents{std::istreambuf_iterator<char>(resumedInput), std::istreambuf_iterator<char>()};
  EXPECT_NE(std::string::npos, resumedContents.find("resume-tail"));
}

TEST(mbx_static_energy, Test_2_CO2_in_MFI_2x2x2_energy_difference_two_rigid_moves)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component c = Component::makeCO2(forceField, 0, true);

  // mbx defines component names in lower case and takes COO type order
  c.name = "co2";
  std::swap(c.definedAtoms[0], c.definedAtoms[1]);
  std::vector<Component> components{c};

  Framework f = Framework::makeMFI(forceField, int3(2, 2, 2));

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, components, {}, {2}, 5);
  system.useMBX = true;
  system.mbxSettingsFilePath = mbxSettingsFile().path().string();

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  setTwoCO2MFIPositions(atomData);
  setupFrameworkPermanentEnergy(system, frameworkAtoms);

  const std::size_t numAtoms = components[0].atoms.size();

  std::vector<Atom> oldMoleculeAtoms(atomData.begin(), atomData.end());
  std::vector<Atom> newMoleculeAtoms(atomData.begin(), atomData.end());

  translateAndRotateMolecule(std::span<Atom>(oldMoleculeAtoms.data(), numAtoms), double3(0.010, -0.015, 0.020),
                             double3(0.004, -0.003, 0.002));
  translateAndRotateMolecule(std::span<Atom>(oldMoleculeAtoms.data() + numAtoms, numAtoms),
                             double3(-0.012, 0.008, -0.010), double3(-0.002, 0.005, -0.004));

  translateAndRotateMolecule(std::span<Atom>(newMoleculeAtoms.data(), numAtoms), double3(0.035, -0.006, 0.016),
                             double3(0.011, -0.007, 0.006));
  translateAndRotateMolecule(std::span<Atom>(newMoleculeAtoms.data() + numAtoms, numAtoms),
                             double3(-0.025, 0.021, -0.018), double3(-0.006, 0.009, -0.005));

  std::vector<Atom> oldatoms0 = moleculeSlice(oldMoleculeAtoms, 0, numAtoms);
  std::vector<Atom> newatoms0 = moleculeSlice(newMoleculeAtoms, 0, numAtoms);
  std::vector<Atom> oldatoms1 = moleculeSlice(oldMoleculeAtoms, 1, numAtoms);
  std::vector<Atom> newatoms1 = moleculeSlice(newMoleculeAtoms, 1, numAtoms);

  RunningEnergy oldEnergy =
      Interactions::computeMBXEnergySystem(system, components, system.simulationBox, system.framework, frameworkAtoms,
                                           std::span<const Atom>(oldMoleculeAtoms.data(), oldMoleculeAtoms.size()));

  RunningEnergy newEnergy =
      Interactions::computeMBXEnergySystem(system, components, system.simulationBox, system.framework, frameworkAtoms,
                                           std::span<const Atom>(newMoleculeAtoms.data(), newMoleculeAtoms.size()));

  RunningEnergy directDifference = newEnergy - oldEnergy;

  RunningEnergy firstDifference = Interactions::computeMBXEnergyDifference(
      system, components, system.simulationBox, system.framework, 0, frameworkAtoms,
      std::span<const Atom>(oldMoleculeAtoms.data(), oldMoleculeAtoms.size()),
      std::span<const Atom>(newatoms0.data(), newatoms0.size()),
      std::span<const Atom>(oldatoms0.data(), oldatoms0.size()));

  std::vector<Atom> midMoleculeAtoms = oldMoleculeAtoms;
  std::copy(newatoms0.begin(), newatoms0.end(), midMoleculeAtoms.begin());

  RunningEnergy secondDifference = Interactions::computeMBXEnergyDifference(
      system, components, system.simulationBox, system.framework, 0, frameworkAtoms,
      std::span<const Atom>(midMoleculeAtoms.data(), midMoleculeAtoms.size()),
      std::span<const Atom>(newatoms1.data(), newatoms1.size()),
      std::span<const Atom>(oldatoms1.data(), oldatoms1.size()));

  const double sequentialDifference = firstDifference.mbxEnergy + secondDifference.mbxEnergy;

  EXPECT_NEAR(oldRigidMoveMbxEnergy, oldEnergy.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(newRigidMoveMbxEnergy, newEnergy.mbxEnergy, mbxEnergyTolerance);

  EXPECT_NEAR(directRigidMoveDifference, directDifference.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(firstRigidMoveDifference, firstDifference.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(secondRigidMoveDifference, secondDifference.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(sequentialRigidMoveDifference, sequentialDifference, mbxEnergyTolerance);

  EXPECT_GT(std::abs(directDifference.mbxEnergy), 1e-14);
  EXPECT_NEAR(directDifference.mbxEnergy, sequentialDifference, mbxEnergyTolerance);
}
