#include <gtest/gtest.h>

import std;
import finite_double;

TEST(FINITE_DOUBLE, classifies_ieee_encodings)
{
  // Inspect the supplied bit patterns, including subnormal numbers, signed zero, infinity and NaNs.
  for (std::uint64_t bits : {0x0000000000000000ULL, 0x8000000000000000ULL, 0x0000000000000001ULL,
                             0x3ff0000000000000ULL, 0x7fefffffffffffffULL, 0xffefffffffffffffULL})
  {
    EXPECT_TRUE(isFiniteDouble(std::bit_cast<double>(bits))) << std::hex << bits;
  }
  for (std::uint64_t bits : {0x7ff0000000000000ULL, 0xfff0000000000000ULL, 0x7ff8000000000001ULL,
                             0xfff8000000000001ULL, 0x7ff0000000000001ULL})
  {
    EXPECT_FALSE(isFiniteDouble(std::bit_cast<double>(bits))) << std::hex << bits;
  }
}
