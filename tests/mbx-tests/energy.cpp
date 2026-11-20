#include <gtest/gtest.h>

import std;
import atom;
import component;
import double3;
import forcefield;
import framework;
import int3;
import running_energy;
import system;

TEST(MBX, TwoComponentEnergyAndExchange)
{
  ForceField ff = ForceField::makeZeoliteForceField(12.0, true, true, true);
  std::vector<Component> components;
  for (std::size_t i = 0; i < 2; ++i)
  {
    Component co2 = Component::makeCO2(ff, i, true);
    co2.name = "co2";
    std::swap(co2.definedAtoms[0], co2.definedAtoms[1]);
    components.push_back(co2);
  }
  System system(ff, std::nullopt, false, 300.0, 1e4, 1.0, Framework::makeMFI(ff, int3(2, 2, 2)), components, {}, {1, 1},
                5);
  system.configureMBX(MBX_SETTINGS_FILE);
  for (std::size_t c = 0; c < 2; ++c)
  {
    auto atoms = system.spanOfMolecule(c, 0);
    const double y = c == 0 ? 6.97475 : 2.97475;
    atoms[0].position = double3(10.011, y, 0.0);
    atoms[1].position = double3(10.011, y, 1.149);
    atoms[2].position = double3(10.011, y, -1.149);
  }

  system.runningEnergies = system.computeTotalEnergies();
  // Reference from the original RASPA–MBX two-CO2 MFI fixture, in RASPA internal units.
  EXPECT_NEAR(444.5797675540443, system.runningEnergies.mbxEnergy, 1e-6);
  auto molecule = system.spanOfMolecule(1, 0);
  EXPECT_NEAR(0.0, system.mbxTrialDifference(1, 0, molecule, {}).mbxEnergy, 1e-6);
  std::vector<Atom> copy(molecule.begin(), molecule.end());
  const auto deletion = system.mbxTrialDifference(1, 0, {}, {});
  const RunningEnergy old = system.runningEnergies;
  system.deleteMolecule(1, 0, molecule);
  system.runningEnergies = system.computeTotalEnergies();
  EXPECT_NEAR(system.runningEnergies.mbxEnergy - old.mbxEnergy, deletion.mbxEnergy, 1e-6);
  EXPECT_NEAR(system.runningEnergies.tail - old.tail, deletion.tail, 1e-8);
  EXPECT_NE(0.0, deletion.tail);
  const auto insertion = system.mbxTrialDifference(1, {}, copy, {});
  EXPECT_NEAR(-deletion.mbxEnergy, insertion.mbxEnergy, 1e-6);
  EXPECT_NEAR(-deletion.tail, insertion.tail, 1e-8);
}
