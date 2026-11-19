#include <gtest/gtest.h>

import std;
import double3;
import transition_matrix;
import transition_matrix_nd;

TEST(JointTMMC, EachAttemptUpdatesThePreMoveRowAndASeparateRetainedState)
{
  TransitionMatrix matrix;
  matrix.doTMMC = true;
  matrix.joint.emplace(std::vector<std::size_t>{2, 4}, std::vector<std::size_t>{3, 6},
                       std::vector<std::size_t>{2, 0}, 1, 1, true);
  auto& joint = *matrix.joint;
  std::vector<std::size_t> counts{5, 0, 2};
  const auto oldRow = matrix.collectionRow(counts, 2);
  EXPECT_EQ(oldRow, 1U);
  joint.bias[oldRow + 3] = std::log(3.0);
  EXPECT_NEAR(matrix.exchangeBias(counts, 2, true), 3.0, 1e-12);
  ++counts[2];
  matrix.recordMove(double3(0.0, 0.7, 0.3), oldRow, 2);
  joint.updateHistogram(counts, -12.0);
  EXPECT_DOUBLE_EQ(joint.cmatrix[5 * oldRow], 0.7);
  EXPECT_DOUBLE_EQ(joint.cmatrix[5 * oldRow + 2], 0.3);
  EXPECT_EQ(joint.overallHistogram[oldRow], 0U);
  EXPECT_EQ(joint.overallHistogram[oldRow + 3], 1U);
  EXPECT_DOUBLE_EQ(joint.sumPotentialEnergy[oldRow + 3], -12.0);
  joint.startProduction();
  EXPECT_EQ(joint.numberOfUpdates, 0U);
  EXPECT_EQ(joint.overallHistogram[oldRow + 3], 0U);
  EXPECT_DOUBLE_EQ(joint.cmatrix[5 * oldRow + 2], 0.0);
  EXPECT_DOUBLE_EQ(joint.bias[oldRow + 3], std::log(3.0));
}

TEST(JointTMMC, FixedWindowRetainsOutwardWeightsWhileUntrackedSpeciesMoveFreely)
{
  TransitionMatrix matrix;
  matrix.doTMMC = true;
  matrix.joint.emplace(std::vector<std::size_t>{4, 1}, std::vector<std::size_t>{4, 1},
                       std::vector<std::size_t>{2, 0}, 1, 1, false);
  std::vector<std::size_t> counts{1, 7, 4};
  EXPECT_DOUBLE_EQ(matrix.exchangeBias(counts, 0, false), 0.0);
  EXPECT_DOUBLE_EQ(matrix.exchangeBias(counts, 2, true), 0.0);
  EXPECT_DOUBLE_EQ(matrix.exchangeBias(counts, 1, true), 1.0);
  matrix.recordMove(double3(0.0, -1.0, 2.0), 0, 2);
  matrix.recordMove(double3(0.25, 0.75, 0.0), 0, 0);
  matrix.recordMove(double3(0.0, 0.4, 0.6), 0, 1);
  EXPECT_EQ(matrix.joint->cmatrix, (std::vector<double>{1.75, 0.0, 1.0, 0.25, 0.0}));
}
