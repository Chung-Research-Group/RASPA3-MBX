module;

export module mc_moves_widom;

import std;

import double3;
import randomnumbers;
import running_energy;
import atom;
import system;

export namespace MC_Moves
{
/**
 * \brief Result of one non-mutating conventional Widom trial.
 *
 * A failed construction has weight zero and no insertion energy.  A completed
 * trial contains the Rosenbluth/Widom sample weight and the potential-energy
 * difference between the hypothetical N+1 state and the live N-particle state.
 */
struct WidomTrialResult
{
  double weight{};
  std::optional<double> insertionEnergy{};
};

/**
 * \brief Performs a Widom insertion move for the specified component.
 *
 * Attempts to insert a molecule of the selected component into the system using
 * Configurational Bias Monte Carlo (CBMC) method. Calculates the energy differences
 * and computes the insertion weight used for chemical potential estimation.
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 * \param selectedComponent Index of the component to perform the Widom move on.
 * \return The sample weight and, for a completed ghost insertion, its insertion energy.
 */
WidomTrialResult WidomMove(RandomNumber& random, System& system, std::size_t selectedComponent);
}  // namespace MC_Moves
