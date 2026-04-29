#include <gtest/gtest.h>

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
import energy_factor;
import running_energy;
import interactions_mbx;

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
  return std::vector<Atom>(moleculeAtoms.begin() + moleculeOffset, moleculeAtoms.begin() + moleculeOffset + numAtoms);
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

  RunningEnergy energyDirect = Interactions::computeMBXEnergy(system, components, system.simulationBox,
                                                              system.framework, 0, frameworkAtoms, atomData, {}, false);

  EXPECT_NEAR(staticMbxEnergy, energySystem.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(staticMbxEnergy, energyDirect.mbxEnergy, mbxEnergyTolerance);
  EXPECT_NEAR(energySystem.mbxEnergy, energyDirect.mbxEnergy, mbxEnergyTolerance);
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