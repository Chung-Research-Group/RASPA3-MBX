module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <print>
#include <source_location>
#include <utility>
#include <vector>
#endif

module transition_matrix;

#ifdef USE_STD_IMPORT
import std;
#endif

import archive;
import double3;

// C(No -> Nn) += p(o -> n)
// C(No -> No) += 1 − p(o -> n)
//
// translation: double3(0.0, 1.0, 0.0)
// insertion: double3(Pacc, 1.0 - Pacc, 0.0)
// insertion overlap-detected: double3(0.0, 1.0, 0.0)
// deletion: double3(0.0, 1.0 - Pacc, Pacc)
// deletion overlap-detected: double3(0.0, 1.0, 0.0)
void TransitionMatrix::updateMatrix(double3 Pacc, std::size_t oldN)
{
  numberOfUpdates++;
  Pacc.clamp(0.0, 1.0);
  cmatrix[oldN - minMacrostate] += Pacc;
};

void TransitionMatrix::updateHistogram(std::size_t N)
{
  if ((N > maxMacrostate) || (N < minMacrostate)) return;
  histogram[N - minMacrostate]++;
}

// return the biasing Factor
double TransitionMatrix::biasFactor(std::size_t newN, std::size_t oldN)
{
  if (!useBias || newN > maxMacrostate) return 1.0;
  double TMMCBias = bias[newN - minMacrostate] - bias[oldN - minMacrostate];
  return std::exp(TMMCBias);
};

// From Vince Shen's pseudo code//
void TransitionMatrix::adjustBias(std::size_t currentCycle)
{
  if ((currentCycle % sampleTMMCEvery != 0) || currentCycle == 0) return;

  numberOfComputes++;

  // get the lowest and highest visited states in terms of loading
  std::size_t minVisitedN = static_cast<std::size_t>(std::distance(
      histogram.begin(), std::find_if(histogram.begin(), histogram.end(), [](const std::size_t &i) { return i; })));
  std::size_t maxVisitedN = static_cast<std::size_t>(
      std::distance(histogram.begin(),
                    std::find_if(histogram.rbegin(), histogram.rend(), [](const std::size_t &i) { return i; }).base()) -
      1);
  [[maybe_unused]] std::size_t nonzeroCount = maxVisitedN - minVisitedN + 1;

  lnpi[minVisitedN] = 0.0;
  double maxlnpi = lnpi[minVisitedN];
  // Update the lnpi for the sampled region//
  // x: -1; y: 0; z: +1//
  for (std::size_t i = minVisitedN; i < maxVisitedN; i++)
  {
    // Zhao's note: add protection to avoid numerical issues
    if (cmatrix[i].z != 0)
    {
      lnpi[i + 1] =
          lnpi[i] + std::log(cmatrix[i].z) - std::log(cmatrix[i].x + cmatrix[i].y + cmatrix[i].z);  // Forward//
    }
    forward_lnpi[i + 1] = lnpi[i + 1];
    if (cmatrix[i + 1].x != 0)
    {
      lnpi[i + 1] = lnpi[i + 1] - std::log(cmatrix[i + 1].x) +
                    std::log(cmatrix[i + 1].x + cmatrix[i + 1].y + cmatrix[i + 1].z);  // Reverse//
    }
    reverse_lnpi[i + 1] = lnpi[i + 1];
    if (lnpi[i + 1] > maxlnpi)
    {
      maxlnpi = lnpi[i + 1];
    }
  }

  // For the unsampled states, fill them with the minVisitedN/maxVisitedN stats
  for (std::size_t i = 0; i < minVisitedN; ++i)
  {
    lnpi[i] = lnpi[minVisitedN];
  }
  for (std::size_t i = maxVisitedN; i < maxMacrostate - minMacrostate + 1; ++i)
  {
    lnpi[i] = lnpi[maxVisitedN];
  }

  // Normalize
  for (std::size_t i = 0; i < maxMacrostate - minMacrostate + 1; ++i)
  {
    lnpi[i] -= maxlnpi;
  }
  double sumExps =
      std::accumulate(lnpi.begin(), lnpi.end(), 0.0, [](double sum, double item) { return sum + std::exp(item); });
  double normalFactor = -std::log(sumExps);

  for (std::size_t i = 0; i < maxMacrostate - minMacrostate + 1; ++i)
  {
    lnpi[i] += normalFactor;  // Zhao's note: mind the sign
    bias[i] = -lnpi[i];
  }
};

// Clear Collection matrix stats (used after initialization cycles)
void TransitionMatrix::clearCMatrix()
{
  if (!rezeroAfterInitialization) return;

  numberOfUpdates = 0;
  numberOfComputes = 0;
  double3 temp = {0.0, 0.0, 0.0};
  std::fill(cmatrix.begin(), cmatrix.end(), temp);
  std::fill(histogram.begin(), histogram.end(), 0);
  std::fill(lnpi.begin(), lnpi.end(), 0.0);
  std::fill(bias.begin(), bias.end(), 1.0);
};

void TransitionMatrix::writeStatistics(std::size_t currentCycle)
{
  if ((currentCycle % writeTMMCEvery != 0) || currentCycle == 0) return;

  std::ofstream textTMMCFile{};
  std::filesystem::path cwd = std::filesystem::current_path();

  std::string dirname = "tmmc/";
  std::string fname = dirname + "/" + "tmmc_statistics.txt";

  std::filesystem::path directoryName = cwd / dirname;
  std::filesystem::path fileName = cwd / fname;
  std::filesystem::create_directories(directoryName);
  textTMMCFile = std::ofstream(fileName, std::ios::out);

  std::print(textTMMCFile, "# performed: {} steps\n", numberOfUpdates);
  std::print(textTMMCFile, "# collection matrix updated: {} times\n", numberOfComputes);
  std::print(textTMMCFile, "# minimum microstate: {}\n", minMacrostate);
  std::print(textTMMCFile, "# maximum microstate: {}\n", maxMacrostate);
  std::print(textTMMCFile, "# column 1: N\n");
  std::print(textTMMCFile, "# column 2: CM[-1]\n");
  std::print(textTMMCFile, "# column 3: CM[ 0]0\n");
  std::print(textTMMCFile, "# column 4: CM[+1]\n");
  std::print(textTMMCFile, "# column 5: bias\n");
  std::print(textTMMCFile, "# column 6: lnpi\n");
  std::print(textTMMCFile, "# column 7: forward lnpi\n");
  std::print(textTMMCFile, "# column 8: reverse lnpi\n");
  std::print(textTMMCFile, "# column 9: histogram\n");
  std::print(textTMMCFile, "N CM[-1] CM[0] CM[1] bias lnpi Forward_lnpi Reverse_lnpi histogram\n");
  for (std::size_t j = minMacrostate; j < maxMacrostate + 1; j++)
  {
    std::size_t newj = j - minMacrostate;
    std::print(textTMMCFile, "{} {:8.5f} {:8.5f} {:8.5f} {} {} {} {} {}\n", j, cmatrix[newj].x, cmatrix[newj].y,
               cmatrix[newj].z, bias[newj], lnpi[newj], forward_lnpi[newj], reverse_lnpi[newj], histogram[newj]);
  }
};

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const TransitionMatrix &m)
{
  archive << m.versionNumber;

  archive << m.cmatrix;
  archive << m.bias;
  archive << m.lnpi;
  archive << m.forward_lnpi;
  archive << m.reverse_lnpi;
  archive << m.histogram;

  archive << m.minMacrostate;
  archive << m.maxMacrostate;
  archive << m.numberOfStates;

  archive << m.sampleTMMCEvery;
  archive << m.writeTMMCEvery;
  archive << m.subSampling;
  archive << m.useBias;

  archive << m.numberOfUpdates;
  archive << m.numberOfComputes;

  archive << m.rejectOutofBound;
  archive << m.rezeroAfterInitialization;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, TransitionMatrix &m)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > m.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'TransitionMatrix' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> m.cmatrix;
  archive >> m.bias;
  archive >> m.lnpi;
  archive >> m.forward_lnpi;
  archive >> m.reverse_lnpi;
  archive >> m.histogram;

  archive >> m.minMacrostate;
  archive >> m.maxMacrostate;
  archive >> m.numberOfStates;

  archive >> m.sampleTMMCEvery;
  archive >> m.writeTMMCEvery;
  archive >> m.subSampling;
  archive >> m.useBias;

  archive >> m.numberOfUpdates;
  archive >> m.numberOfComputes;

  archive >> m.rejectOutofBound;
  archive >> m.rezeroAfterInitialization;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("TransitionMatrix: Error in binary restart\n"));
  }
#endif

  return archive;
}
