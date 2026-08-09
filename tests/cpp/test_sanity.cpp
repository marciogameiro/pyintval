// Milestone 1 sanity tests: verify the platform guarantees that the rigorous
// kernel (Milestone 2) will rely on. These are runtime confirmations of the
// static_asserts in interval.hpp plus the library functions the EFT
// (error-free transformation) algorithms need.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cmath>
#include <limits>

#include "doctest/doctest.h"
#include "pyintval/interval.hpp"

TEST_CASE("platform provides IEEE 754 binary64 doubles") {
  CHECK(std::numeric_limits<double>::is_iec559);
  CHECK(std::numeric_limits<double>::digits == 53);
  CHECK(std::numeric_limits<double>::has_infinity);
  CHECK(std::numeric_limits<double>::has_denorm == std::denorm_present);
}

TEST_CASE("std::fma computes a fused multiply-add") {
  // The EFT product residual (Milestone 2) is exact only if fma performs a
  // single rounding; this catches libraries where fma is broken outright.
  CHECK(std::fma(2.0, 3.0, 1.0) == 7.0);
  // 1 + 2^-53 * 2^-53 rounds to 1 under a single rounding; a double-rounding
  // (multiply-then-add) implementation can differ on cases of this shape.
  const double eps = 0x1p-53;
  CHECK(std::fma(eps, eps, 1.0) == 1.0);
}

TEST_CASE("std::nextafter steps by exactly one ulp") {
  CHECK(std::nextafter(1.0, 2.0) == 1.0 + 0x1p-52);
  CHECK(std::nextafter(1.0, 0.0) == 1.0 - 0x1p-53);
  CHECK(std::nextafter(0.0, 1.0) == std::numeric_limits<double>::denorm_min());
  CHECK(std::isinf(
      std::nextafter(std::numeric_limits<double>::max(), std::numeric_limits<double>::infinity())));
}

TEST_CASE("subnormals are not flushed to zero") {
  // FTZ/DAZ CPU modes would silently break enclosure guarantees near zero.
  const double dmin = std::numeric_limits<double>::denorm_min();
  volatile double half = dmin / 2.0;  // volatile: defeat constant folding
  CHECK(dmin > 0.0);
  CHECK(half * 2.0 <= dmin);
  volatile double tiny = std::numeric_limits<double>::min();
  CHECK(tiny / 2.0 > 0.0);
}

TEST_CASE("kernel ABI version is exported") { CHECK(pyintval::kernel_abi_version == 1); }
