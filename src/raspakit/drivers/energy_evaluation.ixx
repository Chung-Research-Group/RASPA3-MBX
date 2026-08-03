module;

export module energy_evaluation;

import std;

import input_reader;
import system;

/**
 * \brief Evaluates the configured systems once without performing simulation moves.
 *
 * Each system is evaluated through a private copy so initialization of cutoffs,
 * Ewald parameters, and rigid-framework energy leaves the InputReader state intact.
 * Results are printed with RunningEnergy::printMC() and collected in the truncating
 * output/energy_evaluation.json file. MBX-only raw terms are reported in kcal/mol.
 */
export struct EnergyEvaluation
{
  explicit EnergyEvaluation(const InputReader& reader);

  /// Perform exactly one total-energy calculation per owned system copy.
  /// \throws std::runtime_error if a binary restart was requested or output cannot be written.
  void run();

 private:
  bool restartFromBinary{};
  std::vector<System> systems{};
};
