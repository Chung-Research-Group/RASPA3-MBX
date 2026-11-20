module;
#include "bblock/system.h"

module system;

import std;
import atom;
import component;
import framework;
import simulationbox;
import units;

namespace
{
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
};

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

}  // namespace

bool System::mbxAvailable() { return true; }

double System::mbxEnergy(std::size_t component, std::optional<std::size_t> replacedMolecule,
                         std::span<const Atom> trial, std::array<double, 7>* terms) const
{
  if ((replacedMolecule || !trial.empty()) && component >= components.size())
    throw std::out_of_range("MBX trial component is outside the system");
  if (replacedMolecule && *replacedMolecule >= numberOfMoleculesPerComponent.at(component))
    throw std::out_of_range("MBX trial molecule is outside the component");
  bblock::System mbx;
  int nextTag = 1;
  std::optional<std::pair<std::size_t, std::size_t>> excluded;
  if (replacedMolecule) excluded.emplace(component, *replacedMolecule);
  addSystemMonomers(mbx, *this, components, spanOfMoleculeAtoms(), nextTag, excluded);
  if (!trial.empty()) addMonomer(mbx, components.at(component), trial, nextTag);
  mbx.SetUpFromJson(mbxSettings);
  addFrameworkCharges(mbx, framework, spanOfFrameworkAtoms(), nextTag);
  setPeriodicCell(mbx, simulationBox);
  const auto energy = evaluateMBX(mbx, mbxFrameworkPermanentEnergy);
  if (terms)
    *terms = {energy.oneBody,
              energy.twoBody,
              energy.threeBody,
              energy.fourBody,
              energy.dispersion,
              energy.permanentElectrostatics,
              energy.inducedElectrostatics};
  return energy.includedEnergy() / Units::EnergyToKCalPerMol;
}

void System::configureMBX(const std::string& path)
{
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error(std::format("Could not open MBXSettingsFile '{}'", path));
  std::ostringstream buffer;
  buffer << input.rdbuf();
  std::string settings = std::move(buffer).str();
  if (!input || !buffer || settings.empty()) throw std::runtime_error("Could not read MBX settings");
  if (framework && framework->hasMobileAtoms()) throw std::runtime_error("MBX requires a fixed framework");
  bblock::System bare;
  bare.SetUpFromJson(settings);
  double permanentEnergy = 0.0;
  if (framework)
  {
    int nextTag = 1;
    addFrameworkCharges(bare, framework, spanOfFrameworkAtoms(), nextTag);
    setPeriodicCell(bare, simulationBox);
    bare.Electrostatics(false);
    permanentEnergy = bare.GetPermanentElectrostaticEnergy();
  }
  mbxSettings = std::move(settings);
  mbxFrameworkPermanentEnergy = permanentEnergy;
}
