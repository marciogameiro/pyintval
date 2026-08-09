// Deterministic edge-table tests for the directed-rounding primitives.
// The exhaustive differential check against hardware rounding lives in
// test_rounding_hw.cpp; this file pins specific, human-checkable values and
// the structural invariants (RD <= RU, exactness when representable).

#include <cmath>
#include <limits>

#include "doctest/doctest.h"
#include "pyintval/detail/rounding.hpp"

using namespace pyintval::detail;

namespace {
constexpr double kInfv = std::numeric_limits<double>::infinity();
constexpr double kMax = std::numeric_limits<double>::max();
constexpr double kDMin = std::numeric_limits<double>::denorm_min();
constexpr double kMinNormal = std::numeric_limits<double>::min();

// Value-equality that treats +0 and -0 alike (the kernel documents zero signs
// as non-load-bearing) and NaN as never appearing.
bool same(double a, double b) { return a == b || (a == 0.0 && b == 0.0); }
}  // namespace

TEST_CASE("EFT building blocks are exact") {
  auto s = two_sum(1.0, 0x1p-53);
  CHECK(s.val == 1.0);
  CHECK(s.err == 0x1p-53);
  auto p = two_prod(0x1.5p0, 0x1.9p0);
  CHECK(p.val + p.err == 0x1.5p0 * 0x1.9p0);  // reconstructs exactly (small values)
  auto f = fast_two_sum(1e300, 1.0);
  CHECK(f.val == 1e300 + 1.0);
}

TEST_CASE("directed ops agree and are exact when the result is representable") {
  struct Row {
    double a, b;
  };
  const Row rows[] = {{1.0, 1.0}, {3.0, 4.0},  {6.0, 2.0},    {2.25, 0.0},
                      {0.5, 0.5}, {-2.0, 8.0}, {1024.0, 0.5}, {5.0, -3.0}};
  for (auto r : rows) {
    CHECK(add_rd(r.a, r.b) == add_ru(r.a, r.b));
    CHECK(add_rd(r.a, r.b) == r.a + r.b);
    CHECK(sub_rd(r.a, r.b) == sub_ru(r.a, r.b));
    CHECK(sub_rd(r.a, r.b) == r.a - r.b);
    CHECK(mul_rd(r.a, r.b) == mul_ru(r.a, r.b));
    CHECK(mul_rd(r.a, r.b) == r.a * r.b);
  }
  CHECK(div_rd(6.0, 2.0) == 3.0);
  CHECK(div_ru(6.0, 2.0) == 3.0);
  CHECK(sqrt_rd(4.0) == 2.0);
  CHECK(sqrt_ru(4.0) == 2.0);
  CHECK(sqrt_rd(2.25) == 1.5);
  CHECK(sqrt_ru(2.25) == 1.5);
  CHECK(sqrt_rd(0.0) == 0.0);
  CHECK(fma_rd(2.0, 3.0, 1.0) == 7.0);
  CHECK(fma_ru(2.0, 3.0, 1.0) == 7.0);
}

TEST_CASE("directed ops straddle the exact value on inexact results") {
  // Each pair: RD < RU with RD == pred(RU), and the RN result lies between.
  CHECK(add_ru(0.1, 0.2) == succ(add_rd(0.1, 0.2)));
  CHECK(add_rd(0.1, 0.2) <= 0.1 + 0.2);
  CHECK(0.1 + 0.2 <= add_ru(0.1, 0.2));

  CHECK(add_rd(1.0, 0x1p-54) == 1.0);
  CHECK(add_ru(1.0, 0x1p-54) == 1.0 + 0x1p-52);

  CHECK(sub_ru(1.0, 0x1p-54) == 1.0);  // 1 - 2^-54 rounds up to 1
  CHECK(sub_rd(1.0, 0x1p-54) == 1.0 - 0x1p-53);

  CHECK(div_rd(1.0, 3.0) == pred(div_ru(1.0, 3.0)));
  CHECK(sqrt_rd(2.0) == pred(sqrt_ru(2.0)));
  CHECK(mul_rd(0.1, 0.2) == pred(mul_ru(0.1, 0.2)));
}

TEST_CASE("overflow: RD saturates to DBL_MAX, RU to +inf (and mirror)") {
  CHECK(add_ru(kMax, kMax) == kInfv);
  CHECK(add_rd(kMax, kMax) == kMax);
  CHECK(add_rd(-kMax, -kMax) == -kInfv);
  CHECK(add_ru(-kMax, -kMax) == -kMax);
  CHECK(mul_ru(kMax, 2.0) == kInfv);
  CHECK(mul_rd(kMax, 2.0) == kMax);
  CHECK(mul_rd(-kMax, 2.0) == -kInfv);
  CHECK(mul_ru(-kMax, 2.0) == -kMax);
  CHECK(sub_ru(kMax, -kMax) == kInfv);
  CHECK(sub_rd(kMax, -kMax) == kMax);
  CHECK(fma_ru(kMax, 2.0, kMax) == kInfv);
  CHECK(fma_rd(kMax, 2.0, kMax) == kMax);
}

TEST_CASE("genuine infinite operands pass through exactly") {
  CHECK(add_ru(-kInfv, 5.0) == -kInfv);
  CHECK(add_rd(-kInfv, 5.0) == -kInfv);
  CHECK(mul_rd(kInfv, 2.0) == kInfv);
  CHECK(mul_ru(kInfv, 2.0) == kInfv);
  CHECK(same(div_ru(1.0, kInfv), 0.0));
  CHECK(same(div_rd(1.0, kInfv), 0.0));
  CHECK(sqrt_ru(kInfv) == kInfv);
  CHECK(sqrt_rd(kInfv) == kInfv);
}

TEST_CASE("underflow and subnormal rounding") {
  CHECK(mul_rd(kDMin, 0.5) == 0.0);
  CHECK(mul_ru(kDMin, 0.5) == kDMin);
  // 2^-537 * 2^-537 == 2^-1074 exactly (both directions agree).
  CHECK(mul_rd(0x1p-537, 0x1p-537) == kDMin);
  CHECK(mul_ru(0x1p-537, 0x1p-537) == kDMin);
  // A product landing just under the 2^-969 exact-EFT gate stays ordered.
  CHECK(mul_rd(0x1.3p-500, 0x1.7p-480) <= mul_ru(0x1.3p-500, 0x1.7p-480));

  // Division with tiny dividend and huge divisor hitting the sign shortcut:
  // quotient magnitude ~ 2^-1920, far below denorm_min.
  CHECK(same(div_rd(0x1p-1000, 0x1p920), 0.0));
  CHECK(div_ru(0x1p-1000, 0x1p920) == kDMin);
  CHECK(div_rd(-0x1p-1000, 0x1p920) == -kDMin);
  CHECK(same(div_ru(-0x1p-1000, 0x1p920), 0.0));

  // sqrt of subnormals stays a rigorous bracket.
  CHECK(sqrt_rd(kDMin) <= sqrt_ru(kDMin));
  CHECK(sqrt_rd(3.0 * kDMin) <= sqrt_ru(3.0 * kDMin));
  CHECK(mul_rd(sqrt_rd(3.0 * kDMin), sqrt_rd(3.0 * kDMin)) <= 3.0 * kDMin);
  CHECK(mul_ru(sqrt_ru(3.0 * kDMin), sqrt_ru(3.0 * kDMin)) >= 3.0 * kDMin);
}

TEST_CASE("succ and pred at boundaries") {
  CHECK(succ(1.0) == 1.0 + 0x1p-52);
  CHECK(pred(1.0) == 1.0 - 0x1p-53);
  CHECK(succ(0.0) == kDMin);
  CHECK(pred(0.0) == -kDMin);
  CHECK(succ(-kDMin) == 0.0);
  CHECK(pred(kDMin) == 0.0);
  CHECK(succ(kMax) == kInfv);
  CHECK(pred(-kMax) == -kInfv);
  CHECK(succ(kInfv) == kInfv);
  CHECK(pred(-kInfv) == -kInfv);
  CHECK(succ(pred(kMinNormal)) == kMinNormal);  // crosses the normal/subnormal line
}

TEST_CASE("ordering invariant RD <= RU on a structured sweep") {
  const double vals[] = {0.0,      kDMin, 3.0 * kDMin, 0x1p-1000, kMinNormal, 0x1p-537,
                         0x1p-100, 0.1,   1.0,         1.5,       3.0,        0x1p100,
                         0x1p900,  kMax,  -0.1,        -1.0,      -3.0,       -0x1p900};
  for (double a : vals) {
    for (double b : vals) {
      CHECK(add_rd(a, b) <= add_ru(a, b));
      CHECK(sub_rd(a, b) <= sub_ru(a, b));
      CHECK(mul_rd(a, b) <= mul_ru(a, b));
      if (b != 0.0) CHECK(div_rd(a, b) <= div_ru(a, b));
      CHECK(fma_rd(a, b, 1.0) <= fma_ru(a, b, 1.0));
    }
    if (a >= 0.0) CHECK(sqrt_rd(a) <= sqrt_ru(a));
  }
}
