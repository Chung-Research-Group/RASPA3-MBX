#include <gtest/gtest.h>

import std;
import atom;
import pseudo_atom;
import json;
import units;
import component;
import forcefield;
import system;
import simulationbox;

namespace
{
System makeComponentValidationSystem(bool correctOrder)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component component = Component::makeCO2(forceField, 0, true);
  component.name = "co2";
  if (correctOrder) std::swap(component.definedAtoms[0], component.definedAtoms[1]);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 298.0, 1e5, 1.0, {}, {component}, {}, {0}, 5);
  system.useMBX = true;
  return system;
}
}  // namespace

TEST(mbx_component_validation, CorrectCO2OrderIsAccepted)
{
  System system = makeComponentValidationSystem(true);
  EXPECT_NO_THROW(system.validateMBXComponentDefinitions());
}

TEST(mbx_component_validation, ClassicalOCOOrderIsRejectedBeforeEvaluation)
{
  System system = makeComponentValidationSystem(false);
  EXPECT_THROW(system.validateMBXComponentDefinitions(), std::runtime_error);
  EXPECT_THROW(system.computeTotalEnergies(), std::runtime_error);
}

TEST(mbx_component_validation, FlexibleGuestIsRejectedEvenWithoutSwap)
{
  System system = makeComponentValidationSystem(true);
  system.components[0].rigid = false;
  EXPECT_THROW(system.validateMBXMonteCarloConfiguration(), std::runtime_error);
  EXPECT_THROW(system.computeTotalEnergies(), std::runtime_error);
}

TEST(mbx_component_validation, ClassicalConfigurationIsUnaffected)
{
  System system = makeComponentValidationSystem(false);
  system.useMBX = false;
  system.components[0].rigid = false;
  EXPECT_NO_THROW(system.validateMBXComponentDefinitions());
  EXPECT_NO_THROW(system.validateMBXMonteCarloConfiguration());
}

TEST(mbx_component_validation, WaterRequiresOHHOrder)
{
  System system = makeComponentValidationSystem(true);
  system.components[0].name = "h2o";
  const std::array<std::size_t, 3> types{8, 1, 1};
  for (std::size_t i = 0; i < types.size(); ++i)
  {
    system.components[0].definedAtoms[i].first.type = static_cast<std::uint16_t>(i);
    system.forceField.pseudoAtoms[i].atomicNumber = types[i];
  }
  EXPECT_NO_THROW(system.validateMBXComponentDefinitions());
  std::swap(system.components[0].definedAtoms[0], system.components[0].definedAtoms[1]);
  EXPECT_THROW(system.validateMBXComponentDefinitions(), std::runtime_error);
}

TEST(mbx_component_validation, PressureMetadataSeparatesImposedInputAndDiagnostic)
{
  System system = makeComponentValidationSystem(true);
  const double imposedPressure = system.pressure;
  const nlohmann::json status = system.jsonSystemStatus();
  EXPECT_TRUE(status.at("pressureIsImposedInput").get<bool>());
  EXPECT_EQ(status.at("pressureEstimator"), "classical virial diagnostic; MBX virial not implemented");
  EXPECT_DOUBLE_EQ(status.at("pressure").get<double>(), imposedPressure * Units::PressureConversionFactor);
  EXPECT_NE(system.writeMBXStatus().find("MBX virial is not implemented"), std::string::npos);
  EXPECT_DOUBLE_EQ(system.pressure, imposedPressure);
  system.useMBX = false;
  EXPECT_FALSE(system.jsonSystemStatus().contains("pressureEstimator"));
}
