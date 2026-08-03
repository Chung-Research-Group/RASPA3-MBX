module;

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <iterator>
#include <limits>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "bblock/system.h"

module interactions_mbx;

import atom;
import component;
import framework;
import running_energy;
import simulationbox;
import system;
import units;

namespace
{
constexpr std::size_t numberOfMBXEnergyTerms{7};

struct MBXEnergyTerms
{
  double oneBody{};
  double twoBody{};
  double threeBody{};
  double fourBody{};
  double dispersion{};
  double permanentElectrostatics{};
  double inducedElectrostatics{};

  [[nodiscard]] double includedEnergy() const noexcept
  {
    // RASPA-MBX intentionally delegates guest interactions to the many-body terms below. The MBX 1B,
    // Buckingham, Lennard-Jones, and bare-framework permanent electrostatics are not part of this energy.
    return twoBody + threeBody + fourBody + dispersion + permanentElectrostatics + inducedElectrostatics;
  }

  void copyTo(std::span<double> output) const
  {
    if (output.empty()) return;
    if (output.size() != numberOfMBXEnergyTerms)
    {
      throw std::invalid_argument(
          std::format("MBX energy log requires {} entries, received {}", numberOfMBXEnergyTerms, output.size()));
    }
    std::ranges::copy(
        std::array{oneBody, twoBody, threeBody, fourBody, dispersion, permanentElectrostatics, inducedElectrostatics},
        output.begin());
  }
};

[[nodiscard]] std::string readMBXSettings(const System& system)
{
  if (system.mbxSettingsFilePath.empty())
  {
    throw std::runtime_error("[MBX]: no MBX settings file was configured");
  }

  const std::filesystem::path settingsPath(system.mbxSettingsFilePath);
  std::error_code error;
  if (!std::filesystem::is_regular_file(settingsPath, error))
  {
    throw std::runtime_error(std::format("[MBX]: settings file '{}' is not a readable regular file{}",
                                         settingsPath.string(),
                                         error ? std::format(" ({})", error.message()) : std::string{}));
  }

  std::ifstream input(settingsPath, std::ios::binary);
  if (!input)
  {
    throw std::runtime_error(std::format("[MBX]: failed to open settings file '{}'", settingsPath.string()));
  }

  std::string settings{std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
  if (input.bad())
  {
    throw std::runtime_error(std::format("[MBX]: failed while reading settings file '{}'", settingsPath.string()));
  }
  if (settings.empty())
  {
    throw std::runtime_error(std::format("[MBX]: settings file '{}' is empty", settingsPath.string()));
  }
  return settings;
}

void setPeriodicCell(bblock::System& mbx, const SimulationBox& simulationBox)
{
  std::vector<double> box{simulationBox.cell.mm[0][0], simulationBox.cell.mm[0][1], simulationBox.cell.mm[0][2],
                          simulationBox.cell.mm[1][0], simulationBox.cell.mm[1][1], simulationBox.cell.mm[1][2],
                          simulationBox.cell.mm[2][0], simulationBox.cell.mm[2][1], simulationBox.cell.mm[2][2]};

  for (double& value : box)
  {
    if (std::abs(value) < 1.0e-10) value = 0.0;
  }
  mbx.SetPBC(box);
}

void advanceTag(int& nextTag, std::size_t numberOfAtoms)
{
  if (numberOfAtoms > static_cast<std::size_t>(std::numeric_limits<int>::max() - nextTag))
  {
    throw std::overflow_error("[MBX]: atom tags exceed the range supported by MBX");
  }
  nextTag += static_cast<int>(numberOfAtoms);
}

void addMonomer(bblock::System& mbx, const Component& component, std::span<const Atom> atoms, int& nextTag)
{
  if (atoms.size() != component.atoms.size())
  {
    throw std::runtime_error(std::format("[MBX]: component '{}' expects {} atoms per monomer, received {}",
                                         component.name, component.atoms.size(), atoms.size()));
  }

  std::vector<double> coordinates(3 * atoms.size());
  std::vector<std::string> atomNames(atoms.size());
  for (std::size_t atomIndex = 0; atomIndex < atoms.size(); ++atomIndex)
  {
    coordinates[3 * atomIndex] = atoms[atomIndex].position.x;
    coordinates[3 * atomIndex + 1] = atoms[atomIndex].position.y;
    coordinates[3 * atomIndex + 2] = atoms[atomIndex].position.z;
    // The RASPA-MBX settings/database convention uses the numeric RASPA pseudo-atom type as the MBX site name.
    atomNames[atomIndex] = std::to_string(static_cast<std::size_t>(atoms[atomIndex].type));
  }

  mbx.AddMonomer(std::move(coordinates), std::move(atomNames), component.name, 1, nextTag);
  advanceTag(nextTag, atoms.size());
}

void addSystemMonomers(bblock::System& mbx, const System& system, const std::vector<Component>& components,
                       std::span<const Atom> moleculeAtoms, int& nextTag,
                       std::optional<std::pair<std::size_t, std::size_t>> excludedMolecule = std::nullopt)
{
  if (components.size() != system.numberOfMoleculesPerComponent.size())
  {
    throw std::runtime_error("[MBX]: component and molecule-count arrays have inconsistent sizes");
  }

  std::size_t moleculeOffset{};
  for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
  {
    const std::size_t atomsPerMolecule = components[componentId].atoms.size();
    const std::size_t moleculeCount = system.numberOfMoleculesPerComponent[componentId];
    if (atomsPerMolecule == 0 && moleculeCount != 0)
    {
      throw std::runtime_error(
          std::format("[MBX]: component '{}' has molecules but no atoms", components[componentId].name));
    }
    if (moleculeOffset > moleculeAtoms.size() ||
        (atomsPerMolecule != 0 && moleculeCount > (moleculeAtoms.size() - moleculeOffset) / atomsPerMolecule))
    {
      throw std::runtime_error(
          std::format("[MBX]: molecule atom storage is too short for component '{}'", components[componentId].name));
    }

    for (std::size_t moleculeId = 0; moleculeId < moleculeCount; ++moleculeId)
    {
      if (excludedMolecule == std::pair{componentId, moleculeId}) continue;
      const std::size_t base = moleculeOffset + moleculeId * atomsPerMolecule;
      addMonomer(mbx, components[componentId], moleculeAtoms.subspan(base, atomsPerMolecule), nextTag);
    }
    moleculeOffset += atomsPerMolecule * moleculeCount;
  }

  if (moleculeOffset != moleculeAtoms.size())
  {
    throw std::runtime_error(
        std::format("[MBX]: molecule atom storage has {} unaccounted atoms", moleculeAtoms.size() - moleculeOffset));
  }
}

void addFrameworkCharges(bblock::System& mbx, const std::optional<Framework>& framework,
                         std::span<const Atom> frameworkAtoms, int& nextTag)
{
  if (!framework.has_value())
  {
    if (!frameworkAtoms.empty())
    {
      throw std::runtime_error("[MBX]: framework atoms were supplied without a framework");
    }
    return;
  }
  if (frameworkAtoms.empty()) return;

  std::vector<double> coordinates(3 * frameworkAtoms.size());
  std::vector<double> charges(frameworkAtoms.size());
  std::vector<std::size_t> isLocal(frameworkAtoms.size(), 1);
  std::vector<int> tags(frameworkAtoms.size());
  for (std::size_t atomIndex = 0; atomIndex < frameworkAtoms.size(); ++atomIndex)
  {
    coordinates[3 * atomIndex] = frameworkAtoms[atomIndex].position.x;
    coordinates[3 * atomIndex + 1] = frameworkAtoms[atomIndex].position.y;
    coordinates[3 * atomIndex + 2] = frameworkAtoms[atomIndex].position.z;
    charges[atomIndex] = frameworkAtoms[atomIndex].charge;
    tags[atomIndex] = nextTag;
    advanceTag(nextTag, 1);
  }
  mbx.SetExternalChargesAndPositions(std::move(charges), std::move(coordinates), std::move(isLocal), std::move(tags));
}

[[nodiscard]] MBXEnergyTerms evaluateMBX(bblock::System& mbx, double bareFrameworkPermanentEnergy)
{
  if (mbx.GetNumMon() == 0) return {};

  constexpr bool calculateGradients{false};
  MBXEnergyTerms terms;
  terms.oneBody = mbx.OneBodyEnergy(calculateGradients);
  terms.twoBody = mbx.TwoBodyEnergy(calculateGradients);
  terms.threeBody = mbx.ThreeBodyEnergy(calculateGradients);
  terms.fourBody = mbx.FourBodyEnergy(calculateGradients);
  terms.dispersion = mbx.Dispersion(calculateGradients);
  // Buckingham and Lennard-Jones are deliberately neither evaluated nor included: RASPA supplies only
  // the framework-guest VDW term, while MBX supplies the retained many-body guest terms.
  mbx.Electrostatics(calculateGradients);
  terms.permanentElectrostatics = mbx.GetPermanentElectrostaticEnergy() - bareFrameworkPermanentEnergy;
  terms.inducedElectrostatics = mbx.GetInducedElectrostaticEnergy();
  return terms;
}

[[nodiscard]] RunningEnergy asRunningEnergy(const MBXEnergyTerms& terms)
{
  RunningEnergy energy;
  energy.mbxEnergy = terms.includedEnergy() / Units::EnergyToKCalPerMol;
  return energy;
}

[[nodiscard]] std::optional<std::pair<std::size_t, std::size_t>> selectedExistingMolecule(
    const System& system, std::size_t selectedComponent, std::span<const Atom> selectedMoleculeAtoms)
{
  if (selectedMoleculeAtoms.empty() || selectedComponent >= system.numberOfMoleculesPerComponent.size())
  {
    return std::nullopt;
  }
  if (selectedComponent >= system.components.size() ||
      selectedMoleculeAtoms.size() != system.components[selectedComponent].atoms.size())
  {
    return std::nullopt;
  }

  // Atom::moleculeId is global across the component-major molecule storage, whereas addSystemMonomers()
  // enumerates molecules locally within each component. Convert the validated global ID to that local index.
  std::size_t firstGlobalMoleculeId{};
  for (std::size_t componentId = 0; componentId < selectedComponent; ++componentId)
  {
    firstGlobalMoleculeId += system.numberOfMoleculesPerComponent[componentId];
  }

  const std::size_t globalMoleculeId = static_cast<std::size_t>(selectedMoleculeAtoms.front().moleculeId);
  if (globalMoleculeId < firstGlobalMoleculeId) return std::nullopt;

  const std::size_t componentMoleculeId = globalMoleculeId - firstGlobalMoleculeId;
  if (componentMoleculeId >= system.numberOfMoleculesPerComponent[selectedComponent]) return std::nullopt;
  if (!std::ranges::all_of(selectedMoleculeAtoms,
                           [selectedComponent, globalMoleculeId](const Atom& atom)
                           {
                             return static_cast<std::size_t>(atom.componentId) == selectedComponent &&
                                    static_cast<std::size_t>(atom.moleculeId) == globalMoleculeId;
                           }))
  {
    return std::nullopt;
  }
  return std::pair{selectedComponent, componentMoleculeId};
}
}  // namespace

double Interactions::computeFrameworkElecPermMBXEnergy(const System& system, const SimulationBox& simulationBox,
                                                       const std::optional<Framework>& framework,
                                                       std::span<const Atom> frameworkAtoms)
{
  if (!framework.has_value() || frameworkAtoms.empty()) return 0.0;

  bblock::System mbx;
  int nextTag{1};
  mbx.SetUpFromJson(readMBXSettings(system));
  addFrameworkCharges(mbx, framework, frameworkAtoms, nextTag);
  setPeriodicCell(mbx, simulationBox);
  constexpr bool calculateGradients{false};
  mbx.Electrostatics(calculateGradients);
  return mbx.GetPermanentElectrostaticEnergy();
}

RunningEnergy Interactions::computeMBXEnergySystem(const System& system, const std::vector<Component>& components,
                                                   const SimulationBox& simulationBox,
                                                   const std::optional<Framework>& framework,
                                                   std::span<const Atom> frameworkAtoms,
                                                   std::span<const Atom> moleculeAtoms, std::span<double> mbxEnergyLog)
{
  if (!mbxEnergyLog.empty() && mbxEnergyLog.size() != numberOfMBXEnergyTerms)
  {
    throw std::invalid_argument(
        std::format("MBX energy log requires {} entries, received {}", numberOfMBXEnergyTerms, mbxEnergyLog.size()));
  }

  bblock::System mbx;
  int nextTag{1};
  addSystemMonomers(mbx, system, components, moleculeAtoms, nextTag);
  mbx.SetUpFromJson(readMBXSettings(system));
  addFrameworkCharges(mbx, framework, frameworkAtoms, nextTag);
  setPeriodicCell(mbx, simulationBox);

  const MBXEnergyTerms terms = evaluateMBX(mbx, system.elecPermFrameworkMBX);
  terms.copyTo(mbxEnergyLog);
  return asRunningEnergy(terms);
}

RunningEnergy Interactions::computeMBXEnergy(const System& system, const std::vector<Component>& components,
                                             const SimulationBox& simulationBox,
                                             const std::optional<Framework>& framework, std::size_t selectedComponent,
                                             std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
                                             std::span<const Atom> selectedMoleculeAtoms,
                                             bool includeSelectedMoleculeAtoms, std::vector<double>* mbxEnergyLog)
{
  if (selectedComponent >= components.size())
  {
    throw std::out_of_range(std::format("[MBX]: selected component {} is out of range", selectedComponent));
  }
  if (mbxEnergyLog != nullptr && mbxEnergyLog->size() != numberOfMBXEnergyTerms)
  {
    throw std::invalid_argument(
        std::format("MBX energy log requires {} entries, received {}", numberOfMBXEnergyTerms, mbxEnergyLog->size()));
  }
  if (includeSelectedMoleculeAtoms && selectedMoleculeAtoms.empty())
  {
    throw std::invalid_argument("[MBX]: cannot include an empty selected molecule");
  }
  if (!selectedMoleculeAtoms.empty() && selectedMoleculeAtoms.size() != components[selectedComponent].atoms.size())
  {
    throw std::invalid_argument(std::format("[MBX]: selected component '{}' requires {} atoms, received {}",
                                            components[selectedComponent].name,
                                            components[selectedComponent].atoms.size(), selectedMoleculeAtoms.size()));
  }

  bblock::System mbx;
  int nextTag{1};
  addSystemMonomers(mbx, system, components, moleculeAtoms, nextTag,
                    selectedExistingMolecule(system, selectedComponent, selectedMoleculeAtoms));
  if (includeSelectedMoleculeAtoms)
  {
    addMonomer(mbx, components[selectedComponent], selectedMoleculeAtoms, nextTag);
  }
  mbx.SetUpFromJson(readMBXSettings(system));
  addFrameworkCharges(mbx, framework, frameworkAtoms, nextTag);
  setPeriodicCell(mbx, simulationBox);

  const MBXEnergyTerms terms = evaluateMBX(mbx, system.elecPermFrameworkMBX);
  if (mbxEnergyLog != nullptr) terms.copyTo(std::span<double>(*mbxEnergyLog));
  return asRunningEnergy(terms);
}

RunningEnergy Interactions::computeMBXEnergyDifference(
    const System& system, const std::vector<Component>& components, const SimulationBox& simulationBox,
    const std::optional<Framework>& framework, std::size_t selectedComponent, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, std::span<const Atom> newMoleculeAtoms, std::span<const Atom> oldMoleculeAtoms)
{
  if (oldMoleculeAtoms.empty() || newMoleculeAtoms.empty())
  {
    throw std::invalid_argument(
        "[MBX]: replacement energy differences require nonempty old and new molecule configurations; "
        "insertion, deletion, and Widom states must be evaluated explicitly");
  }

  const auto oldMolecule = selectedExistingMolecule(system, selectedComponent, oldMoleculeAtoms);
  const auto newMolecule = selectedExistingMolecule(system, selectedComponent, newMoleculeAtoms);
  if (!oldMolecule.has_value() || !newMolecule.has_value() || oldMolecule != newMolecule)
  {
    throw std::invalid_argument(
        "[MBX]: replacement energy differences require old and new configurations of the same existing molecule");
  }

  const RunningEnergy oldEnergy = computeMBXEnergy(system, components, simulationBox, framework, selectedComponent,
                                                   frameworkAtoms, moleculeAtoms, oldMoleculeAtoms, true);
  const RunningEnergy newEnergy = computeMBXEnergy(system, components, simulationBox, framework, selectedComponent,
                                                   frameworkAtoms, moleculeAtoms, newMoleculeAtoms, true);
  return newEnergy - oldEnergy;
}
