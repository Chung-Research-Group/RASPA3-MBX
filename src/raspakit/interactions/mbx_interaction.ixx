module;

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

export module interactions_mbx;

import atom;
import component;
import framework;
import running_energy;
import simulationbox;
import system;

export namespace Interactions
{
/** Bare-framework permanent electrostatic energy in MBX units (kcal/mol). */
double computeFrameworkElecPermMBXEnergy(const System& system, const SimulationBox& simulationBox,
                                         const std::optional<Framework>& framework,
                                         std::span<const Atom> frameworkAtoms);

/**
 * Absolute MBX guest energy for a complete System. The optional seven-element output span is populated
 * in the order [1B, 2B, 3B, 4B, dispersion, permanent electrostatics, induced electrostatics], in kcal/mol.
 */
RunningEnergy computeMBXEnergySystem(const System& system, const std::vector<Component>& components,
                                     const SimulationBox& simulationBox, const std::optional<Framework>& framework,
                                     std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
                                     std::span<double> mbxEnergyLog = {});

/**
 * Absolute MBX energy with a selected molecule replaced, inserted, or removed. When supplied, mbxEnergyLog
 * must point to a seven-element vector using the same order and units as computeMBXEnergySystem().
 */
RunningEnergy computeMBXEnergy(const System& system, const std::vector<Component>& components,
                               const SimulationBox& simulationBox, const std::optional<Framework>& framework,
                               std::size_t selectedComponent, std::span<const Atom> frameworkAtoms,
                               std::span<const Atom> moleculeAtoms, std::span<const Atom> selectedMoleculeAtoms,
                               bool includeSelectedMoleculeAtoms, std::vector<double>* mbxEnergyLog = nullptr);

/**
 * MBX energy difference for replacing one existing molecule with a trial configuration of that same molecule.
 * Both spans must be nonempty and identify the same existing component/molecule. Insertion, deletion, and Widom
 * moves must instead evaluate their explicit absolute old/new states with computeMBXEnergy().
 */
[[nodiscard]] RunningEnergy computeMBXEnergyDifference(
    const System& system, const std::vector<Component>& components, const SimulationBox& simulationBox,
    const std::optional<Framework>& framework, std::size_t selectedComponent, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, std::span<const Atom> newMoleculeAtoms,
    std::span<const Atom> oldMoleculeAtoms);

}  // namespace Interactions
