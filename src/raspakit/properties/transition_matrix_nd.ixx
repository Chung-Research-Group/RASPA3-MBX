module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#include <vector>
#endif

export module transition_matrix_nd;

import std;

import archive;
import double3;

/**
 * \brief Represents a histogram-based biasing scheme for multidimensional macrostate sampling.
 *
 * The TransitionMatrixMultidimensional struct stores a multidimensional macrostate histogram,
 * the current bias, and the corresponding logarithmic probability estimate \f$\ln \pi\f$.
 * Bias updates are computed directly from histogram counts.
 *
 * In addition, the struct tracks a collection matrix in the same layout as the older
 * TMMC-style implementation. That matrix is retained for diagnostics, restart data,
 * and text output, but it is not used in the bias update.
 */
export struct TransitionMatrixMultidimensional
{
  std::uint64_t versionNumber{0};  ///< Version number for serialization compatibility.

  TransitionMatrixMultidimensional() = default;

  /**
   * \brief Constructs a multidimensional histogram-bias object.
   *
   * The macrostate space is defined by inclusive bounds \p minMacrostate and
   * \p maxMacrostate for each tracked component. The \p componentIds array maps
   * local bias dimensions onto the full system component numbering.
   *
   * The collection-matrix storage is initialized as a per-state block with
   * \f$1 + 2 n_\mathrm{comp}\f$ entries:
   * - one "stay" entry,
   * - one decrement entry per component,
   * - one increment entry per component.
   *
   * \param minMacrostate Inclusive lower macrostate bounds.
   * \param maxMacrostate Inclusive upper macrostate bounds.
   * \param componentIds Mapping from local bias dimensions to system components.
   * \param sampleTMMCEvery Frequency of histogram-based bias updates.
   * \param writeTMMCEvery Frequency of writing statistics to file.
   * \param useBias Whether bias factors should be applied in acceptance rules.
   */
  TransitionMatrixMultidimensional(std::vector<std::size_t> minMacrostate, std::vector<std::size_t> maxMacrostate,
                                   std::vector<std::size_t> componentIds, std::size_t sampleTMMCEvery,
                                   std::size_t writeTMMCEvery, bool useBias)
      : numberOfComponents(minMacrostate.size()),
        numberOfDirections(2 * numberOfComponents + 1),
        minMacrostate(minMacrostate),
        maxMacrostate(maxMacrostate),
        componentIds(componentIds),
        numberOfStates(minMacrostate.size()),
        strides(minMacrostate.size(), 1),
        sampleTMMCEvery(sampleTMMCEvery),
        writeTMMCEvery(writeTMMCEvery),
        useBias(useBias)
  {
    numberOfStates.resize(numberOfComponents);
    strides.assign(numberOfComponents, 1);

    totalNumberOfStates = 1;
    for (std::size_t i = 0; i < numberOfComponents; ++i)
    {
      if (maxMacrostate[i] < minMacrostate[i])
      {
        throw std::runtime_error("Error: histogram bias invalid macrostate bounds.");
      }

      numberOfStates[i] = maxMacrostate[i] - minMacrostate[i] + 1;
      totalNumberOfStates *= numberOfStates[i];
    }

    for (std::size_t i = numberOfComponents; i > 1; --i)
    {
      strides[i - 2] = strides[i - 1] * numberOfStates[i - 1];
    }

    cmatrix.resize(totalNumberOfStates * numberOfDirections, 0.0);
    bias.resize(totalNumberOfStates, 0.0);
    lnpi.resize(totalNumberOfStates, 0.0);
    histogram.resize(totalNumberOfStates, 0);
    overallHistogram.resize(totalNumberOfStates, 0);
    sumPotentialEnergy.resize(totalNumberOfStates, 0.0);

    for (std::size_t s = 0; s < totalNumberOfStates; ++s)
    {
      auto local = unravelIndex(s);

      std::vector<std::size_t> absolute(numberOfComponents);
      for (std::size_t i = 0; i < numberOfComponents; ++i)
      {
        absolute[i] = minMacrostate[i] + local[i];
      }

      if (selectMatrixIndex(absolute) != s)
      {
        throw std::runtime_error("Error: histogram bias index/unravel mismatch.");
      }
    }
  };

  bool operator==(TransitionMatrixMultidimensional const &) const = default;

  std::size_t numberOfComponents{};         ///< Number of biased macrostate dimensions.
  std::size_t numberOfDirections{};         ///< Per-state collection-matrix entries: stay plus \f$\pm 1\f$ moves.
  std::size_t totalNumberOfStates{};        ///< Total number of discrete macrostates.
  std::vector<std::size_t> minMacrostate;   ///< Inclusive minimum macrostate value per component.
  std::vector<std::size_t> maxMacrostate;   ///< Inclusive maximum macrostate value per component.
  std::vector<std::size_t> componentIds;    ///< Mapping from local macrostate dimensions to system components.
  std::vector<std::size_t> numberOfStates;  ///< Number of discrete states per macrostate dimension.
  std::vector<std::size_t> strides;         ///< Row-major strides for flattening multidimensional indices.

  std::vector<double> cmatrix;                ///< Collection matrix tracked for output/restart purposes only.
  std::vector<double> bias;                   ///< Bias value associated with each macrostate.
  std::vector<double> lnpi;                   ///< Estimated \f$\ln \pi\f$ for each macrostate.
  std::vector<std::size_t> histogram;         ///< Histogram accumulated since the last bias update.
  std::vector<std::size_t> overallHistogram;  ///< Histogram accumulated over the current stage.
  std::vector<double> sumPotentialEnergy;     ///< Sum of potential energy samples per macrostate.

  std::size_t sampleTMMCEvery{};  ///< Update interval for histogram-based bias refinement.
  std::size_t writeTMMCEvery{};   ///< Output interval for writing statistics to file.
  bool useBias{true};             ///< Enables or disables the multiplicative bias factor.
  bool equilibration{true};       ///< Whether the object is still in equilibration mode.

  double maxBias{30.0};  ///< Upper cap applied to the bias during equilibration.
  double gamma{0.5};     ///< Relaxation parameter used in histogram-based bias updates.

  std::size_t numberOfUpdates = {0};   ///< Number of recorded matrix/histogram updates.
  std::size_t numberOfComputes = {0};  ///< Number of bias recomputations.

  bool rejectOutofBound = {true};            ///< Whether out-of-range moves should be rejected by the caller.
  bool rezeroAfterInitialization = {false};  ///< Whether initialization statistics should be reset before production.

  /**
   * \brief Converts absolute particle numbers into a flattened macrostate index.
   *
   * \param particleNumbers Absolute particle numbers for the tracked components.
   * \return Flattened row-major state index.
   */
  std::size_t selectMatrixIndex(const std::vector<std::size_t> &particleNumbers) const;

  /**
   * \brief Converts a flattened macrostate index into local per-dimension indices.
   *
   * The returned indices are local coordinates in the range
   * \f$[0,\,\texttt{numberOfStates}[i])\f$.
   *
   * \param state Flattened state index.
   * \return Local multidimensional index.
   */
  std::vector<std::size_t> unravelIndex(std::size_t state) const;

  /**
   * \brief Converts local per-dimension indices into absolute particle numbers.
   *
   * \param local Local multidimensional index.
   * \return Absolute particle numbers.
   */
  std::vector<std::size_t> indexToParticleNumbers(const std::vector<std::size_t> &local) const;

  /**
   * \brief Extracts the subset of particle numbers that belongs to the biased components.
   *
   * \param particleNumbers Full particle-number vector.
   * \return Particle numbers restricted to the tracked macrostate dimensions.
   */
  std::vector<std::size_t> pruneParticleNumbers(const std::vector<std::size_t> &particleNumbers) const;

  /**
   * \brief Maps a global component id onto its local macrostate index.
   *
   * \param selectedComponent Global component id.
   * \return Local component index within the tracked macrostate.
   */
  std::size_t getLocalIndex(const std::size_t selectedComponent) const;

  /**
   * \brief Updates the tracked collection matrix.
   *
   * The input \p Pacc is interpreted as deletion, stay, and insertion probabilities
   * in the order \p x, \p y, \p z. The selected component determines which pair of
   * decrement/increment directions is updated inside the per-state collection-matrix block.
   *
   * This matrix is tracked and written to file, but it is not used in \ref adjustBias.
   *
   * \param Pacc Acceptance probabilities for deletion (\p x), stay (\p y), and insertion (\p z).
   * \param particleNumbers Current particle numbers.
   * \param selectedComponent Global component id of the attempted move.
   */
  void updateMatrix(double3 Pacc, const std::vector<std::size_t> &particleNumbers, std::size_t selectedComponent);

  /**
   * \brief Updates the current and overall histograms of macrostate visits.
   *
   * \param particleNumbers Current particle numbers.
   * \param potentialEnergy Current potential energy sample.
   */
  void updateHistogram(const std::vector<std::size_t> &particleNumbers, double potentialEnergy);

  /**
   * \brief Computes the multiplicative bias factor for a neighboring macrostate move.
   *
   * For an insertion or deletion attempt in the selected component, this returns
   * \f$\exp(\Delta \mathrm{bias})\f$ between the destination and source macrostates.
   * If the attempted move would leave the allowed macrostate window, zero is returned.
   *
   * \param particleNumbers Current particle numbers.
   * \param selectedComponent Global component id of the attempted move.
   * \param insert True for insertion, false for deletion.
   * \return Multiplicative bias factor for the proposed move.
   */
  double biasFactor(const std::vector<std::size_t> &particleNumbers, std::size_t selectedComponent, bool insert);

  /**
   * \brief Refines the bias from the accumulated histogram.
   *
   * The current histogram is converted into an estimate of \f$\ln \pi\f$, normalized,
   * and mixed into the existing bias with relaxation parameter \ref gamma. After each
   * update, the current histogram is reset to zero. During equilibration, the bias
   * is capped by \ref maxBias.
   *
   * \param currentCycle Current simulation cycle.
   */
  void adjustBias(std::size_t currentCycle);

  /**
   * \brief Clears initialization statistics if enabled.
   *
   * Resets the collection matrix, histograms, \f$\ln \pi\f$, bias, and counters when
   * \ref rezeroAfterInitialization is true.
   */
  void clearCMatrix();

  /**
   * \brief Switches from equilibration to production sampling.
   *
   * Clears both histograms and disables equilibration. The tracked collection matrix
   * is intentionally left unchanged.
   */
  void startProduction();

  /**
   * \brief Writes the current statistics to a text file.
   *
   * The output includes macrostate bounds, collection-matrix entries, bias,
   * \f$\ln \pi\f$, the current histogram, the overall histogram, and the average
   * potential energy per macrostate.
   *
   * \param currentCycle Current simulation cycle.
   */
  void writeStatistics(std::size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const TransitionMatrixMultidimensional &m);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, TransitionMatrixMultidimensional &m);
};
