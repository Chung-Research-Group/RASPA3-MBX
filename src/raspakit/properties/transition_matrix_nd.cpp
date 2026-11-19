module;

module transition_matrix_nd;

import std;

import archive;
import double3;

std::vector<std::size_t> TransitionMatrixMultidimensional::unravelIndex(std::size_t state) const
{
  if (state >= totalNumberOfStates)
  {
    throw std::runtime_error("Error: histogram bias state index exceeds bounds\n");
  }

  std::vector<std::size_t> indices(numberOfComponents, 0);
  for (std::size_t i = 0; i < numberOfComponents; ++i)
  {
    indices[i] = (state / strides[i]) % numberOfStates[i];
  }

  return indices;
}

std::size_t TransitionMatrixMultidimensional::collectionRow(const std::vector<std::size_t> &particleNumbers) const
{
  std::size_t row = 0;
  for (std::size_t i = 0; i < numberOfComponents; ++i)
  {
    const std::size_t n = particleNumbers.at(componentIds[i]);
    if (n < minMacrostate[i] || n > maxMacrostate[i])
      throw std::out_of_range("Population is outside the joint TMMC window");
    row += (n - minMacrostate[i]) * strides[i];
  }
  return row;
}

std::size_t TransitionMatrixMultidimensional::getLocalIndex(const std::size_t selectedComponent) const
{
  for (std::size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    if (componentIds[comp] == selectedComponent) return comp;
  }
  return numberOfComponents;  // An untracked component leaves the tracked particle numbers unchanged.
}

void TransitionMatrixMultidimensional::updateMatrix(double3 Pacc, std::size_t oldRow, std::size_t selectedComponent)
{
  numberOfUpdates++;
  Pacc.clamp(0.0, 1.0);

  std::size_t index = numberOfDirections * oldRow;
  std::size_t localComponentIndex = getLocalIndex(selectedComponent);

  if (localComponentIndex == numberOfComponents)
  {
    cmatrix[index] += 1.0;
    return;
  }
  cmatrix[index] += Pacc.y;
  cmatrix[index + 2 * localComponentIndex + 1] += Pacc.x;
  cmatrix[index + 2 * localComponentIndex + 2] += Pacc.z;
}

void TransitionMatrixMultidimensional::updateHistogram(const std::vector<std::size_t> &particleNumbers,
                                                       double potentialEnergy)
{
  std::size_t index = collectionRow(particleNumbers);

  histogram[index]++;
  overallHistogram[index]++;
  sumPotentialEnergy[index] += potentialEnergy;
}

double TransitionMatrixMultidimensional::biasFactor(const std::vector<std::size_t> &particleNumbers,
                                                    std::size_t selectedComponent, bool insert)
{
  const std::size_t axis = getLocalIndex(selectedComponent);
  if (axis == numberOfComponents) return 1.0;
  const std::size_t n = particleNumbers[selectedComponent];
  if (insert ? n >= maxMacrostate[axis] : n <= minMacrostate[axis]) return 0.0;
  if (!useBias) return 1.0;
  const std::size_t originalIndex = collectionRow(particleNumbers);
  const std::size_t newIndex = insert ? originalIndex + strides[axis] : originalIndex - strides[axis];
  double TMMCBias = bias[newIndex] - bias[originalIndex];
  return std::exp(TMMCBias);
}

void TransitionMatrixMultidimensional::adjustBias(std::size_t currentCycle)
{
  if (!useBias || (currentCycle % sampleTMMCEvery != 0) || currentCycle == 0) return;

  numberOfComputes++;

  double maxLnpi = -std::numeric_limits<double>::infinity();

  for (std::size_t state = 0; state < totalNumberOfStates; ++state)
  {
    double count = static_cast<double>(histogram[state]);
    if (count <= 0.0) count = 1.0;

    lnpi[state] = std::log(count) - bias[state];
    maxLnpi = std::max(maxLnpi, lnpi[state]);
  }

  double sumExps = 0.0;
  for (std::size_t state = 0; state < totalNumberOfStates; ++state)
  {
    sumExps += std::exp(lnpi[state] - maxLnpi);
  }

  const double normalFactor = -(maxLnpi + std::log(sumExps));
  for (std::size_t state = 0; state < totalNumberOfStates; ++state)
  {
    lnpi[state] += normalFactor;
  }

  for (std::size_t state = 0; state < totalNumberOfStates; ++state)
  {
    bias[state] = (1.0 - gamma) * bias[state] - gamma * lnpi[state];
  }

  if (equilibration)
  {
    for (std::size_t state = 0; state < totalNumberOfStates; ++state)
    {
      bias[state] = std::min(maxBias, bias[state]);
    }
  }

  std::fill(histogram.begin(), histogram.end(), 0);
}

void TransitionMatrixMultidimensional::startProduction()
{
  // Physical collection, visits and energy averages in the production file must share a baseline.
  std::fill(cmatrix.begin(), cmatrix.end(), 0.0);
  numberOfUpdates = 0;
  std::fill(histogram.begin(), histogram.end(), 0);
  std::fill(overallHistogram.begin(), overallHistogram.end(), 0);
  std::fill(sumPotentialEnergy.begin(), sumPotentialEnergy.end(), 0.0);
  equilibration = false;
}

void TransitionMatrixMultidimensional::writeStatistics(const std::string &filename) const
{
  const std::filesystem::path path(filename);
  if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
  std::ofstream textTMMCFile(path);
  if (!textTMMCFile) throw std::runtime_error("Cannot open TMMC statistics: " + filename);

  std::print(textTMMCFile, "# performed: {} steps\n", numberOfUpdates);
  std::print(textTMMCFile, "# bias updated: {} times\n", numberOfComputes);
  std::print(textTMMCFile, "# equilibration: {}\n", equilibration);
  std::print(textTMMCFile, "# max_bias: {}\n", maxBias);

  for (std::size_t comp = 0; comp < numberOfComponents; comp++)
  {
    std::print(textTMMCFile, "# axis {} component id: {}\n", comp, componentIds[comp]);
    std::print(textTMMCFile, "# minimum macrostate comp {}: {}\n", comp, minMacrostate[comp]);
    std::print(textTMMCFile, "# maximum macrostate comp {}: {}\n", comp, maxMacrostate[comp]);
  }

  std::size_t columnNumber = 1;
  std::print(textTMMCFile, "# column {}: State\n", columnNumber++);
  for (std::size_t comp = 0; comp < numberOfComponents; comp++)
  {
    std::print(textTMMCFile, "# column {}: N_{}\n", columnNumber++, comp);
  }
  std::print(textTMMCFile, "# column {}: CM\n", columnNumber++);
  for (std::size_t comp = 0; comp < numberOfComponents; comp++)
  {
    std::print(textTMMCFile, "# column {}: CM_{}^-\n", columnNumber++, comp);
    std::print(textTMMCFile, "# column {}: CM_{}^+\n", columnNumber++, comp);
  }
  std::print(textTMMCFile, "# column {}: bias\n", columnNumber++);
  std::print(textTMMCFile, "# column {}: lnpi\n", columnNumber++);
  std::print(textTMMCFile, "# column {}: histogram_overall\n", columnNumber++);
  std::print(textTMMCFile, "# column {}: average_potential_energy\n", columnNumber++);

  for (std::size_t state = 0; state < totalNumberOfStates; ++state)
  {
    std::vector<std::size_t> indices = unravelIndex(state);
    std::print(textTMMCFile, "{} ", state);
    for (std::size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      std::print(textTMMCFile, "{} ", minMacrostate[comp] + indices[comp]);
    }

    // The collection matrix is stored as one contiguous block per state:
    // [stay, comp0-, comp0+, comp1-, comp1+, ...].
    std::print(textTMMCFile, "{} ", cmatrix[state * numberOfDirections]);
    for (std::size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      std::print(textTMMCFile, "{} ", cmatrix[state * numberOfDirections + 2 * comp + 1]);
      std::print(textTMMCFile, "{} ", cmatrix[state * numberOfDirections + 2 * comp + 2]);
    }

    double averagePotentialEnergy = 0.0;
    if (overallHistogram[state] > 0)
    {
      averagePotentialEnergy = sumPotentialEnergy[state] / static_cast<double>(overallHistogram[state]);
    }

    std::print(textTMMCFile, "{} {} {} {}\n", bias[state], lnpi[state], overallHistogram[state],
               averagePotentialEnergy);
  }
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const TransitionMatrixMultidimensional &m)
{
  archive << m.versionNumber;

  archive << m.cmatrix;
  archive << m.bias;
  archive << m.lnpi;
  archive << m.histogram;
  archive << m.overallHistogram;
  archive << m.sumPotentialEnergy;

  archive << m.minMacrostate;
  archive << m.maxMacrostate;
  archive << m.componentIds;
  archive << m.numberOfStates;
  archive << m.strides;

  archive << m.numberOfComponents;
  archive << m.numberOfDirections;
  archive << m.totalNumberOfStates;

  archive << m.sampleTMMCEvery;
  archive << m.writeTMMCEvery;
  archive << m.useBias;
  archive << m.equilibration;
  archive << m.maxBias;
  archive << m.gamma;

  archive << m.numberOfUpdates;
  archive << m.numberOfComputes;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, TransitionMatrixMultidimensional &m)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > m.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'TransitionMatrixMultidimensional' at line {} in file {}\n",
                    location.line(), location.file_name()));
  }

  archive >> m.cmatrix;
  archive >> m.bias;
  archive >> m.lnpi;
  archive >> m.histogram;
  archive >> m.overallHistogram;
  archive >> m.sumPotentialEnergy;

  archive >> m.minMacrostate;
  archive >> m.maxMacrostate;
  archive >> m.componentIds;
  archive >> m.numberOfStates;
  archive >> m.strides;

  archive >> m.numberOfComponents;
  archive >> m.numberOfDirections;
  archive >> m.totalNumberOfStates;

  archive >> m.sampleTMMCEvery;
  archive >> m.writeTMMCEvery;
  archive >> m.useBias;
  archive >> m.equilibration;
  archive >> m.maxBias;
  archive >> m.gamma;

  archive >> m.numberOfUpdates;
  archive >> m.numberOfComputes;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("TransitionMatrixMultidimensional: Error in binary restart\n"));
  }
#endif
  return archive;
}
