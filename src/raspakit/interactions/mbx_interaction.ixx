module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#include "bblock/system.h"
#endif

export module interactions_mbx;

#ifndef USE_LEGACY_HEADERS
import <span>;
import <optional>;
import <tuple>;
import <vector>;
#include "bblock/system.h"
#endif

import double3;
import double3x3;
import atom;
import running_energy;
import energy_status;
import simulationbox;
import energy_factor;
import gradient_factor;
import component;
import system;
import component;
import framework;

export namespace Interactions
{
// MBX Initialization function
// bblock::System* initialize_mbx();
/**
 * \brief Computes the inter-molecular energy between atoms.
 *
 * Calculates the van der Waals and Coulombic energy contributions between all pairs of atoms
 * in \p moleculeAtoms, excluding interactions within the same molecule.
 *
 * \param forceField The force field parameters used for the energy calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of atoms for which to compute inter-molecular energies.
 * \return The total inter-molecular energy contributions.
 */

double computeFrameworkElecPermMBXEnergy(
        const System &system,
        const SimulationBox &box,
        const std::optional<Framework> &framework,
        std::span<const Atom> frameworkAtoms
) noexcept;

RunningEnergy computeMBXEnergy(
        const System &system,
        const std::vector<Component> &components,
        const SimulationBox &box,
        const std::optional<Framework> &framework,
        size_t selectedComponent,
        std::span<const Atom> frameworkAtoms,
        std::span<const Atom> moleculeAtoms,
        std::span<const Atom> selectedMoleculeAtoms,
        bool includeSelectedMoleculeAtoms,
        std::vector<double>* mbxEnergyLog = nullptr
) noexcept;


RunningEnergy computeMBXEnergySystem(
        const System &system,
        const std::vector<Component> &components,
        const SimulationBox &box,
        const std::optional<Framework> &framework,
        std::span<const Atom> frameworkAtoms,
        std::span<const Atom> moleculeAtoms,
        std::vector<double>* mbxEnergyLog = nullptr
) noexcept;
/**
 * \brief Computes the difference in inter-molecular energy due to atom changes.
 *
 * Calculates the energy difference resulting from replacing \p oldatoms with \p newatoms
 * in the system, considering interactions with \p moleculeAtoms. Excludes interactions within
 * the same molecule. Returns std::nullopt if an overlap is detected based on the force field's overlap criteria.
 *
 * \param forceField The force field parameters used for the energy calculations.
 * \param simulationBox The simulation box containing the atoms.
 * \param moleculeAtoms A span of existing atoms in the system.
 * \param newatoms A span of new atoms to be added to the system.
 * \param oldatoms A span of atoms to be removed from the system.
 * \return The energy difference due to the atom changes, or std::nullopt if an overlap occurs.
 */
[[nodiscard]] RunningEnergy computeMBXEnergyDifference(
        const System &system,
        const std::vector<Component> &components,
        const SimulationBox &simulationBox,
        const std::optional<Framework> &framework,
        size_t selectedComponent,
        std::span<const Atom> frameworkAtoms,
        std::span<const Atom> moleculeAtoms,
        std::span<const Atom> newatoms,
        std::span<const Atom> oldatoms
) noexcept;

}