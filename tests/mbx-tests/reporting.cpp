#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import archive;
import component;
import energy_status;
import framework;
import property_energy;
import units;

namespace
{
// Synthetic samples only: no System, molecular evaluator, random draws or simulation.
PropertyEnergy syntheticEnergies(const std::array<double, 3>& mbxValues)
{
  PropertyEnergy energies(3, 0, 0, 0);
  for (std::size_t block = 0; block < mbxValues.size(); ++block)
  {
    EnergyStatus sample(0, 0, 0, true);
    sample.mbxEnergy = mbxValues[block];
    sample.totalEnergy.energy = mbxValues[block] + 10.0;
    energies.addSample(block, sample, 2.0);
  }
  return energies;
}

std::string writeStatistics(const PropertyEnergy& energies, bool useMBX)
{
  std::optional<Framework> framework;
  std::vector<Component> components;
  return energies.writeAveragesStatistics(false, framework, components, useMBX);
}

std::string archiveBytes(const PropertyEnergy& energies, const std::filesystem::path& path)
{
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << energies;
  }
  std::ifstream stream(path, std::ios::binary);
  return std::string(std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>());
}

void expectMbxSection(const std::string& report, const std::array<double, 3>& blocks, double mean, double confidence)
{
  const std::size_t begin = report.find("MBX:\n");
  ASSERT_NE(begin, std::string::npos);
  const std::size_t end = report.find("Kinetic Energies:\n", begin);
  ASSERT_NE(end, std::string::npos);
  const std::string section = report.substr(begin, end - begin);
  for (std::size_t block = 0; block < blocks.size(); ++block)
  {
    EXPECT_NE(section.find(std::format("Block[ {:2d}] {: .6e}", block, Units::EnergyToKelvin * blocks[block])),
              std::string::npos);
  }
  EXPECT_NE(section.find(std::format("Average  {: .6e} +/- {: .6e} [{}]", Units::EnergyToKelvin * mean,
                                     Units::EnergyToKelvin * confidence, Units::displayedUnitOfEnergyString)),
            std::string::npos);
}

void expectTotalSection(const std::string& report, double mean, double confidence)
{
  const std::size_t begin = report.find("Total energy");
  ASSERT_NE(begin, std::string::npos);
  EXPECT_NE(
      report.substr(begin).find(std::format("Average  {: .6e} +/- {: .6e} [{}]", Units::EnergyToKelvin * mean,
                                            Units::EnergyToKelvin * confidence, Units::displayedUnitOfEnergyString)),
      std::string::npos);
}
}  // namespace

TEST(mbx_reporting, authoritative_model_shows_nonzero_statistics_without_mutating_accumulators)
{
  const std::array<double, 3> values{2.0, 4.0, 6.0};
  const PropertyEnergy energies = syntheticEnergies(values);
  ASSERT_FALSE(energies.zeroSample.useMBX);
  for (const auto& block : energies.bookKeeping) ASSERT_FALSE(block.first.useMBX);
  const auto [mean, confidence] = energies.average();
  ASSERT_FALSE(mean.useMBX);
  // Three equally weighted block means: sample standard deviation 2, Student-t(2, 95%) = 4.303.
  const double expectedConfidence = 4.303 * 2.0 / std::sqrt(3.0);
  EXPECT_DOUBLE_EQ(mean.mbxEnergy, 4.0);
  EXPECT_DOUBLE_EQ(mean.totalEnergy.energy, 14.0);
  EXPECT_NEAR(confidence.mbxEnergy, expectedConfidence, 1e-12);
  EXPECT_NEAR(confidence.totalEnergy.energy, expectedConfidence, 1e-12);

  TemporaryDirectory directory;
  const std::string before = archiveBytes(energies, directory.path() / "before.bin");
  const std::string report = writeStatistics(energies, true);
  expectMbxSection(report, values, 4.0, expectedConfidence);
  expectTotalSection(report, 14.0, expectedConfidence);
  EXPECT_EQ(archiveBytes(energies, directory.path() / "after.bin"), before);
}

TEST(mbx_reporting, authoritative_model_shows_exactly_zero_mbx_mean)
{
  const std::array<double, 3> values{-2.0, 0.0, 2.0};
  const PropertyEnergy energies = syntheticEnergies(values);
  ASSERT_DOUBLE_EQ(energies.average().first.mbxEnergy, 0.0);
  const double expectedConfidence = 4.303 * 2.0 / std::sqrt(3.0);
  const std::string report = writeStatistics(energies, true);
  expectMbxSection(report, values, 0.0, expectedConfidence);
  expectTotalSection(report, 10.0, expectedConfidence);
}

TEST(mbx_reporting, classical_model_suppresses_mbx_even_with_stale_true_metadata)
{
  PropertyEnergy energies = syntheticEnergies({2.0, 4.0, 6.0});
  energies.zeroSample.useMBX = true;
  for (auto& block : energies.bookKeeping) block.first.useMBX = true;
  ASSERT_TRUE(energies.average().first.useMBX);
  TemporaryDirectory directory;
  const std::string before = archiveBytes(energies, directory.path() / "before.bin");
  const std::string report = writeStatistics(energies, false);
  EXPECT_EQ(report.find("MBX:\n"), std::string::npos);
  expectTotalSection(report, 14.0, 4.303 * 2.0 / std::sqrt(3.0));
  EXPECT_EQ(archiveBytes(energies, directory.path() / "after.bin"), before);
}

TEST(mbx_reporting, restored_false_metadata_does_not_hide_mbx_statistics)
{
  const std::array<double, 3> values{2.0, 4.0, 6.0};
  const PropertyEnergy original = syntheticEnergies(values);
  TemporaryDirectory directory;
  const std::filesystem::path path = directory.path() / "legacy_false_metadata.bin";
  const std::string originalBytes = archiveBytes(original, path);
  PropertyEnergy restored;
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  // Exercise persisted false metadata produced by the pre-repair accumulator.
  // This is a synthetic archive in the current format, not a historical full restart.
  ASSERT_EQ(restored.numberOfBlocks, 3u);
  ASSERT_FALSE(restored.zeroSample.useMBX);
  for (const auto& block : restored.bookKeeping) ASSERT_FALSE(block.first.useMBX);
  ASSERT_FALSE(restored.average().first.useMBX);
  const std::string report = writeStatistics(restored, true);
  expectMbxSection(report, values, 4.0, 4.303 * 2.0 / std::sqrt(3.0));
  expectTotalSection(report, 14.0, 4.303 * 2.0 / std::sqrt(3.0));
  EXPECT_EQ(archiveBytes(restored, directory.path() / "after.bin"), originalBytes);
}
