module;

module energy_evaluation;

import std;

import component;
import forcefield;
import input_reader;
import json;
import running_energy;
import simulationbox;
import system;

namespace
{
nlohmann::json moleculeCounts(const System& system)
{
  nlohmann::json counts;
  counts["total"] = system.moleculeData.size();
  counts["perComponent"] = nlohmann::json::array();

  for (std::size_t componentId = 0; componentId != system.components.size(); ++componentId)
  {
    counts["perComponent"].push_back({
        {"componentId", componentId},
        {"name", system.components[componentId].name},
        {"total", system.numberOfMoleculesPerComponent.at(componentId)},
        {"integer", system.numberOfIntegerMoleculesPerComponent.at(componentId)},
        {"fractional", system.numberOfFractionalMoleculesPerComponent.at(componentId)},
    });
  }

  return counts;
}

nlohmann::json rawMBXTerms(const std::array<double, 7>& terms)
{
  return {
      {"oneBody (1B)", terms[0]},          {"twoBody (2B)", terms[1]}, {"threeBody (3B)", terms[2]},
      {"fourBody (4B)", terms[3]},         {"dispersion", terms[4]},   {"permanentElectrostatics", terms[5]},
      {"inducedElectrostatics", terms[6]},
  };
}
}  // namespace

EnergyEvaluation::EnergyEvaluation(const InputReader& inputReader)
    : restartFromBinary(inputReader.restartFromBinary), systems(inputReader.systems)
{
}

void EnergyEvaluation::run()
{
  if (restartFromBinary)
  {
    throw std::runtime_error(
        "[EnergyEvaluation]: RestartFromBinary is not supported for a one-shot energy evaluation; "
        "use the configured JSON restart coordinates instead");
  }

  for (std::size_t systemId = 0; systemId != systems.size(); ++systemId)
  {
    const System& system = systems[systemId];
    const bool hasFractionalConfiguration =
        std::ranges::any_of(system.components, [](const Component& component)
                            { return component.hasFractionalMolecule || component.fixedLambdaBin.has_value(); }) ||
        std::ranges::any_of(system.numberOfFractionalMoleculesPerComponent,
                            [](std::size_t count) { return count != 0uz; });
    if (hasFractionalConfiguration)
    {
      throw std::runtime_error(
          std::format("[EnergyEvaluation]: system {} contains a CFCMC or fixed-lambda fractional molecule; "
                      "the coordinate JSON format does not store lambda/scaling state, so only integer-molecule "
                      "snapshots can be evaluated safely",
                      systemId));
    }
  }

  std::error_code directoryError;
  std::filesystem::create_directories("output", directoryError);
  if (directoryError)
  {
    throw std::runtime_error(
        std::format("[EnergyEvaluation]: could not create output directory: {}", directoryError.message()));
  }

  nlohmann::json output{
      {"simulationType", "EnergyEvaluation"},
      {"systems", nlohmann::json::array()},
  };

  for (std::size_t systemId = 0; systemId != systems.size(); ++systemId)
  {
    // All evaluator-specific initialization is intentionally confined to the owned copy.
    System& system = systems[systemId];
    system.forceField.initializeAutomaticCutOff(system.simulationBox);
    system.forceField.initializeEwaldParameters(system.simulationBox);
    system.precomputeTotalRigidEnergy();

    std::array<double, 7> mbxTerms{};
    const RunningEnergy energy =
        system.useMBX ? system.computeTotalEnergies(std::span<double>(mbxTerms)) : system.computeTotalEnergies();

    std::print(std::cout, "Energy evaluation for system {}\n{}", systemId, energy.printMC());

    nlohmann::json systemOutput{
        {"systemId", systemId},
        {"model", system.useMBX ? "MBX" : "Classical"},
        {"simulationBox", system.simulationBox.jsonStatus()},
        {"moleculeCounts", moleculeCounts(system)},
        {"energy", energy.jsonMC()},
    };

    if (system.useMBX)
    {
      systemOutput["mbxRawEnergyTerms [kcal/mol]"] = rawMBXTerms(mbxTerms);
      systemOutput["mbxEnergyConvention"] =
          "The raw MBX one-body (1B) term is reported here but excluded from the RASPA-MBX total energy.";
    }

    output["systems"].push_back(std::move(systemOutput));
  }

  const std::filesystem::path outputPath{"output/energy_evaluation.json"};
  std::ofstream stream(outputPath, std::ios::out | std::ios::trunc);
  if (!stream)
  {
    throw std::runtime_error(std::format("[EnergyEvaluation]: could not open '{}' for writing", outputPath.string()));
  }

  stream << output.dump(2) << '\n';
  stream.flush();
  if (!stream)
  {
    throw std::runtime_error(std::format("[EnergyEvaluation]: failed while writing '{}'", outputPath.string()));
  }
}
