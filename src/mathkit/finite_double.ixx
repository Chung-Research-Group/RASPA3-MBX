module;

export module finite_double;

import std;

// Check the IEEE exponent directly so NaN and infinity remain observable under -ffast-math.
// Volatile character reads preserve the object representation; a floating-value
// isfinite/bit_cast expression can otherwise be optimized to unconditional true.
export [[nodiscard]] inline bool isFiniteDouble(const double &value)
{
  static_assert(std::numeric_limits<double>::is_iec559 && sizeof(double) == sizeof(std::uint64_t));
  const volatile unsigned char *bytes = reinterpret_cast<const volatile unsigned char *>(&value);
  std::array<unsigned char, sizeof(double)> observed{};
  for (std::size_t i = 0; i < observed.size(); ++i) observed[i] = bytes[i];
  constexpr std::uint64_t exponentMask{0x7ff0000000000000ULL};
  return (std::bit_cast<std::uint64_t>(observed) & exponentMask) != exponentMask;
}
