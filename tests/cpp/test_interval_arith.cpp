// Milestone 2: semantic tests for the bare-interval arithmetic kernel
// (include/pyintval/interval.hpp) against IEEE 1788-2015 set-based reference
// semantics.
//
// Structure:
//   * empty propagation sweep over every operation;
//   * add/sub: exactness, outward rounding, overflow, unbounded endpoints;
//   * mul: every sign-class cell (finite and infinite representatives), the
//     0 * inf hazard cells, and a systematic 13x13 class sweep asserting that
//     no NaN endpoint ever escapes;
//   * div: zero-free divisors (exact and 1-ulp directed cases), the full
//     zero-containing-divisor table, recip;
//   * sqr, sqrt, fma (including a class sweep with the fma-tighter-than-
//     mul-then-add subset property), pown, abs/min/max, floor/ceil/trunc/
//     round/sign;
//   * randomized property tests with fixed seeds (CI-reproducible):
//     pointwise containment, inclusion isotonicity, and algebraic identities
//     (neg involution, sub == add-of-negation, bit-identical).
//
// Mathematical intent: every interval operation must return an enclosure of
// the exact image hull{ f(x) : x in X restricted to dom(f) } with outward
// rounding; the deterministic tables below pin the tightest representable
// enclosures for dyadic data (where no rounding may occur at all), and the
// property tests check the enclosure and monotonicity laws on random data.

#include <algorithm>
#include <bit>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "pyintval/interval.hpp"

namespace pv = pyintval;

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

pv::Interval I(double lo, double hi) { return pv::make(lo, hi); }
pv::Interval P(double v) { return pv::point(v); }

// Bit-identical comparison (distinguishes -0.0 from +0.0; the kernel
// canonicalizes zero endpoints to +0.0, so canonical results must round-trip
// exactly).
bool bit_equal(const pv::Interval& a, const pv::Interval& b) {
  return std::bit_cast<std::uint64_t>(a.lo) == std::bit_cast<std::uint64_t>(b.lo) &&
         std::bit_cast<std::uint64_t>(a.hi) == std::bit_cast<std::uint64_t>(b.hi);
}

// Extended-real membership: z in [lo, hi] with infinities allowed as limits
// (an RN pointwise result that overflowed to +-inf certifies the exact value
// exceeded the largest finite double, so the matching endpoint must be
// infinite for the enclosure to hold).
bool encloses(const pv::Interval& r, double z) {
  return !pv::is_empty(r) && r.lo <= z && z <= r.hi;
}

// Structural invariant every operation result must satisfy: no NaN endpoint,
// and either the canonical empty representation or ordered bounds excluding
// the point-at-infinity forms.
bool well_formed(const pv::Interval& r) {
  if (std::isnan(r.lo) || std::isnan(r.hi)) return false;
  if (pv::is_empty(r)) return r.lo == kInf && r.hi == -kInf;
  return r.lo <= r.hi && r.lo < kInf && r.hi > -kInf;
}

void check_interval_eq(const pv::Interval& got, const pv::Interval& want) {
  CHECK_MESSAGE(pv::equal(got, want), "got " << got << ", want " << want);
}

// One representative per sign-class shape, mixing bounded and unbounded
// endpoints: zero, positive/negative/mixed with 0 as interior point,
// endpoint, or absent, and entire. Every mul/fma pairing of these exercises
// every branch of the sign-class case analysis, including all cells where a
// naive endpoint product would form 0 * inf = NaN.
std::vector<pv::Interval> sign_class_reps() {
  return {P(0.0),        I(1.0, 2.0),   I(0.0, 2.0),    I(2.0, kInf),  I(0.0, kInf),
          I(-2.0, -1.0), I(-2.0, 0.0),  I(-kInf, -2.0), I(-kInf, 0.0), I(-2.0, 3.0),
          I(-kInf, 3.0), I(-2.0, kInf), pv::entire()};
}

// Deterministic member points of a nonempty interval (finite endpoints,
// midpoint, and small integers when contained) -- every one is a member, so
// any pointwise image must land inside the interval result.
std::vector<double> fixed_members(const pv::Interval& x) {
  std::vector<double> m;
  if (pv::is_empty(x)) return m;
  if (std::isfinite(x.lo)) m.push_back(x.lo);
  if (std::isfinite(x.hi)) m.push_back(x.hi);
  m.push_back(pv::mid(x));
  for (double c : {0.0, 1.0, -1.0, 100.0, -100.0}) {
    if (pv::is_member(c, x)) m.push_back(c);
  }
  return m;
}

// --- Random-testing infrastructure ------------------------------------------
//
// Generation is derived directly from the mt19937_64 output (whose sequence is
// standardized) rather than from std::uniform_*_distribution, whose mapping is
// implementation-defined -- so these tests generate the SAME cases under
// libstdc++, libc++, and MSVC's STL, and a failure reproduces on every
// platform.

// Portable uniform double in [0, 1) from the top 53 bits.
double u01(std::mt19937_64& rng) { return static_cast<double>(rng() >> 11) * 0x1p-53; }

// Portable uniform integer in [lo, hi] (a negligible modulo bias is irrelevant
// for randomized testing).
long urange(std::mt19937_64& rng, long lo, long hi) {
  return lo + static_cast<long>(rng() % static_cast<std::uint64_t>(hi - lo + 1));
}

// Endpoint generator mixing zeros, infinities, subnormals, near-overflow
// magnitudes, small integers, and general finite doubles across ~700 binades.
double random_endpoint(std::mt19937_64& rng) {
  switch (urange(rng, 0, 15)) {
    case 0:
      return 0.0;
    case 1:
      return -kInf;
    case 2:
      return kInf;
    case 3: {  // subnormal / deep-underflow magnitudes
      const double v = std::ldexp(static_cast<double>(urange(rng, 1, (1L << 52) - 1)), -1074);
      return (rng() & 1) ? v : -v;
    }
    case 4:
    case 5:  // near-overflow magnitudes
      return (2.0 * u01(rng) - 1.0) * std::numeric_limits<double>::max();
    case 6:  // small integers (exact-arithmetic corner cases)
      return static_cast<double>(urange(rng, -10, 10));
    default: {  // general finite: mantissa in [1,2) over exponents [-350,350]
      const double v = std::ldexp(1.0 + u01(rng), static_cast<int>(urange(rng, -350, 350)));
      return (rng() & 1) ? v : -v;
    }
  }
}

pv::Interval random_interval(std::mt19937_64& rng) {
  for (;;) {
    double a = random_endpoint(rng);
    double b = random_endpoint(rng);
    if (a > b) std::swap(a, b);
    if (pv::valid_bounds(a, b)) return pv::make(a, b);  // rejects [+inf,+inf] etc.
  }
}

// Member samples: finite endpoints, midpoint, and random interior points
// (unbounded sides clamped to the largest finite double, which is a member).
void sample_members(const pv::Interval& x, std::mt19937_64& rng, std::vector<double>& out) {
  out.clear();
  if (pv::is_empty(x)) return;
  if (std::isfinite(x.lo)) out.push_back(x.lo);
  if (std::isfinite(x.hi)) out.push_back(x.hi);
  out.push_back(pv::mid(x));
  const double flo = std::max(x.lo, -std::numeric_limits<double>::max());
  const double fhi = std::min(x.hi, std::numeric_limits<double>::max());
  for (int k = 0; k < 2; ++k) {
    const double t = u01(rng);
    double m = flo * (1.0 - t) + fhi * t;  // may overflow to +-inf; clamp fixes it
    m = std::clamp(m, flo, fhi);
    out.push_back(m);
  }
}

// Random subinterval of x (all sampled members lie in x, so the result is a
// subset by construction); occasionally keeps x itself or one unbounded side.
pv::Interval shrink(const pv::Interval& x, std::mt19937_64& rng) {
  if (pv::is_empty(x)) return x;
  std::vector<double> m;
  sample_members(x, rng, m);
  const auto pick = [&]() { return m[urange(rng, 0, static_cast<long>(m.size()) - 1)]; };
  switch (urange(rng, 0, 3)) {
    case 0:
      return x;
    case 1:
      return pv::make(x.lo, pick());  // keep (possibly unbounded) lower side
    case 2:
      return pv::make(pick(), x.hi);  // keep (possibly unbounded) upper side
    default: {
      double a = pick();
      double b = pick();
      if (a > b) std::swap(a, b);
      return pv::make(a, b);
    }
  }
}

}  // namespace

// --- Empty propagation -------------------------------------------------------

TEST_CASE("empty operands propagate: every operation returns empty") {
  const pv::Interval e = pv::empty();
  const std::vector<pv::Interval> others = {pv::empty(),  P(0.0),        I(1.0, 2.0),
                                            I(-2.0, 3.0), I(-kInf, 2.0), pv::entire()};

  struct UnaryCase {
    const char* name;
    pv::Interval (*fn)(const pv::Interval&);
  };
  const UnaryCase unary[] = {
      {"neg", pv::neg},
      {"recip", pv::recip},
      {"sqr", pv::sqr},
      {"sqrt", pv::sqrt},
      {"abs", pv::abs},
      {"floor", pv::floor},
      {"ceil", pv::ceil},
      {"trunc", pv::trunc},
      {"round_ties_to_even", pv::round_ties_to_even},
      {"round_ties_to_away", pv::round_ties_to_away},
      {"sign", pv::sign},
  };
  for (const auto& u : unary) {
    CAPTURE(u.name);
    CHECK(pv::is_empty(u.fn(e)));
    CHECK(well_formed(u.fn(e)));
  }

  struct BinaryCase {
    const char* name;
    pv::Interval (*fn)(const pv::Interval&, const pv::Interval&);
  };
  const BinaryCase binary[] = {{"add", pv::add}, {"sub", pv::sub}, {"mul", pv::mul},
                               {"div", pv::div}, {"min", pv::min}, {"max", pv::max}};
  for (const auto& b : binary) {
    CAPTURE(b.name);
    for (const auto& o : others) {
      CAPTURE(o);
      CHECK(pv::is_empty(b.fn(e, o)));
      CHECK(pv::is_empty(b.fn(o, e)));
    }
  }

  // fma: empty in each of the three slots.
  const std::vector<pv::Interval> small = {P(0.0), I(1.0, 2.0), pv::entire()};
  for (const auto& o1 : small) {
    for (const auto& o2 : small) {
      CAPTURE(o1);
      CAPTURE(o2);
      CHECK(pv::is_empty(pv::fma(e, o1, o2)));
      CHECK(pv::is_empty(pv::fma(o1, e, o2)));
      CHECK(pv::is_empty(pv::fma(o1, o2, e)));
    }
  }

  // pown: every exponent, including 0 (empty wins over the u^0 = 1 rule).
  for (int n : {-3, -2, -1, 0, 1, 2, 3, 8}) {
    CAPTURE(n);
    CHECK(pv::is_empty(pv::pown(e, n)));
  }
}

// --- Addition / subtraction --------------------------------------------------

TEST_CASE("add/sub: exact dyadic endpoints incur no rounding") {
  check_interval_eq(pv::add(I(1, 2), I(3, 4)), I(4, 6));
  check_interval_eq(pv::sub(I(1, 2), I(3, 4)), I(-3, -1));
  check_interval_eq(pv::add(I(-1.5, 2.5), P(0.5)), I(-1, 3));
  check_interval_eq(pv::sub(I(4, 6), I(3, 4)), I(0, 3));
  check_interval_eq(pv::add(I(-2, 3), P(0.0)), I(-2, 3));
}

TEST_CASE("add/sub: inexact endpoints are rounded outward by exactly one ulp") {
  // 0.1 + 0.2 is inexact in binary64: the enclosure must be the 1-ulp
  // interval bracketing the exact sum, and the RN sum is a member.
  const pv::Interval s = pv::add(P(0.1), P(0.2));
  CHECK(s.lo < s.hi);
  CHECK(std::nextafter(s.lo, kInf) == s.hi);
  CHECK(encloses(s, 0.1 + 0.2));

  // 1.0 - 0.1 is likewise inexact; the enclosure brackets the exact
  // difference and therefore contains RN(0.9).
  const pv::Interval d = pv::sub(P(1.0), P(0.1));
  CHECK(d.lo < d.hi);
  CHECK(std::nextafter(d.lo, kInf) == d.hi);
  CHECK(encloses(d, 1.0 - 0.1));
  CHECK(encloses(d, 0.9));

  // Counterpoint: sub(P(0.3), P(0.1)) is EXACT in binary64 -- the exact
  // difference of the two nearest doubles is 0x1.9999999999999p-3 (about
  // 0.19999999999999998), which is representable. A 1788-tight kernel must
  // return that singleton and must NOT widen it to cover the real 0.2, which
  // is not in the image of these point operands.
  const pv::Interval x = pv::sub(P(0.3), P(0.1));
  CHECK(pv::is_singleton(x));
  CHECK(x.lo == 0.3 - 0.1);
  CHECK(encloses(x, 0x1.9999999999999p-3));
}

TEST_CASE("add/sub: overflow saturates outward on the overflowing side only") {
  check_interval_eq(pv::add(P(DBL_MAX), P(DBL_MAX)), I(DBL_MAX, kInf));
  check_interval_eq(pv::sub(P(-DBL_MAX), P(DBL_MAX)), I(-kInf, -DBL_MAX));
  check_interval_eq(pv::add(I(-DBL_MAX, DBL_MAX), I(-DBL_MAX, DBL_MAX)), pv::entire());
}

TEST_CASE("add/sub: unbounded endpoints pass through exactly") {
  check_interval_eq(pv::add(I(-kInf, 2), P(1)), I(-kInf, 3));
  check_interval_eq(pv::add(I(-kInf, 2), I(3, kInf)), pv::entire());
  check_interval_eq(pv::add(I(2, kInf), I(3, kInf)), I(5, kInf));
  check_interval_eq(pv::sub(P(1), I(-kInf, 2)), I(-1, kInf));
  check_interval_eq(pv::sub(I(-kInf, 2), I(3, kInf)), I(-kInf, -1));
  check_interval_eq(pv::add(pv::entire(), P(1)), pv::entire());
}

// --- Multiplication ----------------------------------------------------------

TEST_CASE("mul: every sign-class cell, finite and infinite representatives") {
  struct MulRow {
    pv::Interval a, b, want;
  };
  const MulRow rows[] = {
      // P*P
      {I(1, 2), I(3, 4), I(3, 8)},
      {I(1, 2), I(3, kInf), I(3, kInf)},
      {I(2, kInf), I(3, 4), I(6, kInf)},
      {I(2, kInf), I(3, kInf), I(6, kInf)},
      {I(0, 2), I(3, 4), I(0, 8)},
      {I(0, kInf), I(3, 4), I(0, kInf)},
      // P*N
      {I(1, 2), I(-4, -3), I(-8, -3)},
      {I(1, 2), I(-kInf, -3), I(-kInf, -3)},
      {I(2, kInf), I(-4, -3), I(-kInf, -6)},
      {I(0, 2), I(-4, -3), I(-8, 0)},
      {I(0, kInf), I(-kInf, -3), I(-kInf, 0)},
      // P*M
      {I(1, 2), I(-2, 3), I(-4, 6)},
      {I(1, 2), I(-kInf, 3), I(-kInf, 6)},
      {I(1, 2), I(-2, kInf), I(-4, kInf)},
      {I(2, kInf), I(-2, 3), pv::entire()},
      {I(0, 2), I(-2, 3), I(-4, 6)},
      {I(0, kInf), I(-2, 3), pv::entire()},
      // N*P (mirrors of P*N; commutativity is also asserted below)
      {I(-2, -1), I(3, 4), I(-8, -3)},
      {I(-kInf, -1), I(3, 4), I(-kInf, -3)},
      {I(-2, -1), I(3, kInf), I(-kInf, -3)},
      {I(-2, 0), I(3, 4), I(-8, 0)},
      // N*N
      {I(-2, -1), I(-4, -3), I(3, 8)},
      {I(-kInf, -1), I(-4, -3), I(3, kInf)},
      {I(-2, -1), I(-kInf, -3), I(3, kInf)},
      {I(-kInf, -1), I(-kInf, -3), I(3, kInf)},
      {I(-2, 0), I(-4, -3), I(0, 8)},
      // N*M
      {I(-2, -1), I(-5, 7), I(-14, 10)},
      {I(-kInf, -1), I(-5, 7), pv::entire()},
      {I(-2, -1), I(-5, kInf), I(-kInf, 10)},
      {I(-2, 0), I(-5, 7), I(-14, 10)},
      // M*P
      {I(-2, 3), I(3, 4), I(-8, 12)},
      {I(-kInf, 3), I(3, 4), I(-kInf, 12)},
      {I(-2, kInf), I(3, 4), I(-8, kInf)},
      {I(-2, 3), I(3, kInf), pv::entire()},
      // M*N
      {I(-2, 3), I(-4, -3), I(-12, 8)},
      {I(-kInf, 3), I(-4, -3), I(-12, kInf)},
      {I(-2, kInf), I(-4, -3), I(-kInf, 8)},
      {I(-2, 3), I(-kInf, -3), pv::entire()},
      // M*M: min over the cross products, max over the same-sign products:
      // lo = min((-2)*7, 3*(-5)) = -15, hi = max((-2)*(-5), 3*7) = 21.
      {I(-2, 3), I(-5, 7), I(-15, 21)},
      {I(-2, 3), I(-kInf, 7), pv::entire()},
      {I(-kInf, 3), I(-5, 7), pv::entire()},
      {I(-kInf, 3), I(-2, kInf), pv::entire()},
      // Zero classes: [0,0] absorbs everything, including unbounded operands.
      {P(0.0), pv::entire(), P(0.0)},
      {P(0.0), I(3, kInf), P(0.0)},
      {P(0.0), I(-kInf, -3), P(0.0)},
      {P(0.0), I(-2, 3), P(0.0)},
      {P(0.0), P(0.0), P(0.0)},
      // 0 * inf hazard cells: the exact image treats 0 * unbounded as 0, so a
      // zero endpoint must never poison the result with NaN.
      {I(0, 1), I(2, kInf), I(0, kInf)},
      {I(0, 1), I(-kInf, -2), I(-kInf, 0)},
      {I(-1, 0), I(2, kInf), I(-kInf, 0)},
      {I(-1, 1), I(0, kInf), pv::entire()},
      {I(0, kInf), I(0, kInf), I(0, kInf)},
      {I(0, kInf), I(-kInf, 0), I(-kInf, 0)},
      {I(-kInf, 0), I(-kInf, 0), I(0, kInf)},
      {I(-1, 0), I(-kInf, 0), I(0, kInf)},
  };
  for (const auto& r : rows) {
    CAPTURE(r.a);
    CAPTURE(r.b);
    check_interval_eq(pv::mul(r.a, r.b), r.want);
    // Multiplication of real sets is commutative and the tightest enclosure
    // is unique, so the kernel must agree under operand swap.
    check_interval_eq(pv::mul(r.b, r.a), r.want);
  }
}

TEST_CASE("mul: 13x13 sign-class sweep never yields NaN and encloses member products") {
  const auto reps = sign_class_reps();
  for (const auto& a : reps) {
    for (const auto& b : reps) {
      CAPTURE(a);
      CAPTURE(b);
      const pv::Interval r = pv::mul(a, b);
      CHECK(well_formed(r));
      CHECK_FALSE(pv::is_empty(r));  // product of nonempty sets is nonempty
      CHECK(pv::equal(r, pv::mul(b, a)));
      for (double x : fixed_members(a)) {
        for (double y : fixed_members(b)) {
          const double z = x * y;  // finite * finite: never NaN
          CHECK_MESSAGE(encloses(r, z), "x=" << x << " y=" << y << " x*y=" << z << " r=" << r);
        }
      }
    }
  }
}

// --- Division ----------------------------------------------------------------

TEST_CASE("div: zero-free divisors, exact and directed") {
  // Dyadic data: exact, no rounding allowed.
  check_interval_eq(pv::div(I(1, 2), I(4, 8)), I(0.125, 0.5));
  check_interval_eq(pv::div(I(-2, -1), I(4, 8)), I(-0.5, -0.125));
  check_interval_eq(pv::div(I(-2, 1), I(4, 8)), I(-0.5, 0.25));
  check_interval_eq(pv::div(I(1, 2), I(-8, -4)), I(-0.5, -0.125));
  check_interval_eq(pv::div(I(-2, -1), I(-8, -4)), I(0.125, 0.5));
  check_interval_eq(pv::div(I(-2, 1), I(-8, -4)), I(-0.25, 0.5));
  // Unbounded operands.
  check_interval_eq(pv::div(I(1, 2), I(2, kInf)), I(0, 1));
  check_interval_eq(pv::div(I(1, 2), I(-kInf, -2)), I(-1, 0));
  check_interval_eq(pv::div(I(3, kInf), I(4, 8)), I(0.375, kInf));
  check_interval_eq(pv::div(I(-kInf, -3), I(4, 8)), I(-kInf, -0.375));

  // [1,1]/[3,3]: 1/3 is not representable, so the tightest enclosure is the
  // 1-ulp interval [RD(1/3), RU(1/3)] and the RN quotient is a member.
  const pv::Interval t = pv::div(P(1), P(3));
  CHECK(t.lo == pv::detail::div_rd(1.0, 3.0));
  CHECK(t.hi == pv::detail::div_ru(1.0, 3.0));
  CHECK(t.lo < t.hi);
  CHECK(std::nextafter(t.lo, kInf) == t.hi);
  CHECK(encloses(t, 1.0 / 3.0));
}

TEST_CASE("div: zero-containing divisors, full table") {
  // X / [0,0] = empty for every X: the domain of division excludes v = 0, so
  // restricting to a divisor with no nonzero element leaves nothing to map.
  for (const auto& x :
       {P(0.0), I(1, 2), I(-2, -1), I(-2, 3), pv::entire(), I(2, kInf), I(-kInf, 0.0)}) {
    CAPTURE(x);
    CHECK(pv::is_empty(pv::div(x, P(0.0))));
  }
  // [0,0] / Y = [0,0] whenever Y has a nonzero element (0/v = 0 for v != 0).
  for (const auto& y :
       {I(1, 2), I(-3, -2), I(-2, 3), I(0, 4), I(-4, 0), I(2, kInf), pv::entire()}) {
    CAPTURE(y);
    check_interval_eq(pv::div(P(0.0), y), P(0.0));
  }
  // 0 interior to the divisor: both half-lines appear, hull is entire.
  check_interval_eq(pv::div(I(1, 2), I(-1, 1)), pv::entire());
  // 0 as a divisor endpoint: exactly one direction blows up.
  check_interval_eq(pv::div(I(1, 2), I(0, 1)), I(1, kInf));
  check_interval_eq(pv::div(I(1, 2), I(-1, 0)), I(-kInf, -1));
  check_interval_eq(pv::div(I(-2, -1), I(0, 1)), I(-kInf, -1));
  check_interval_eq(pv::div(I(-2, -1), I(-1, 0)), I(1, kInf));
  check_interval_eq(pv::div(I(-2, 1), I(0, 1)), pv::entire());
  check_interval_eq(pv::div(I(-2, 1), I(-1, 0)), pv::entire());
  check_interval_eq(pv::div(I(1, 2), I(0, 4)), I(0.25, kInf));
}

TEST_CASE("recip: reciprocal composes the division semantics") {
  check_interval_eq(pv::recip(I(2, 4)), I(0.25, 0.5));
  check_interval_eq(pv::recip(I(-4, -2)), I(-0.5, -0.25));
  check_interval_eq(pv::recip(I(0, 1)), I(1, kInf));
  check_interval_eq(pv::recip(I(-1, 0)), I(-kInf, -1));
  check_interval_eq(pv::recip(I(-1, 1)), pv::entire());
  CHECK(pv::is_empty(pv::recip(P(0.0))));
  check_interval_eq(pv::recip(pv::entire()), pv::entire());
  check_interval_eq(pv::recip(I(1, kInf)), I(0, 1));
}

// --- Square, square root -----------------------------------------------------

TEST_CASE("sqr: image of u^2, exact on dyadic data") {
  check_interval_eq(pv::sqr(I(-3, 2)), I(0, 9));
  check_interval_eq(pv::sqr(I(2, 3)), I(4, 9));
  check_interval_eq(pv::sqr(I(-3, -2)), I(4, 9));
  check_interval_eq(pv::sqr(pv::entire()), I(0, kInf));
  check_interval_eq(pv::sqr(I(-kInf, -2)), I(4, kInf));
  check_interval_eq(pv::sqr(P(0.0)), P(0.0));

  // Enclosure of member squares (including inexact ones), and never wider
  // than the dependency-blind product enclosure.
  for (const auto& x : {I(-3, 2), I(0.1, 0.3), I(-0.7, 0.2), I(2, kInf), P(0.1)}) {
    CAPTURE(x);
    const pv::Interval r = pv::sqr(x);
    CHECK(pv::subset(r, pv::mul(x, x)));
    for (double v : fixed_members(x)) {
      CHECK_MESSAGE(encloses(r, v * v), "v=" << v << " v*v=" << v * v << " r=" << r);
    }
  }
}

TEST_CASE("sqrt: domain-clipped image with directed endpoints") {
  check_interval_eq(pv::sqrt(I(4, 9)), I(2, 3));
  check_interval_eq(pv::sqrt(I(-4, 9)), I(0, 3));  // negative part clipped
  CHECK(pv::is_empty(pv::sqrt(I(-4, -1))));        // empty domain intersection
  check_interval_eq(pv::sqrt(I(-4, 0.0)), P(0.0));
  check_interval_eq(pv::sqrt(pv::entire()), I(0, kInf));
  check_interval_eq(pv::sqrt(I(0, kInf)), I(0, kInf));
  check_interval_eq(pv::sqrt(I(4, kInf)), I(2, kInf));

  // sqrt([0,2]): upper endpoint is the round-up of the irrational sqrt(2).
  const pv::Interval s02 = pv::sqrt(I(0, 2));
  CHECK(s02.lo == 0.0);
  CHECK(s02.hi == pv::detail::sqrt_ru(2.0));
  CHECK(s02.hi >= std::sqrt(2.0));

  // sqrt([2,2]) is the 1-ulp bracket of sqrt(2), containing the RN value.
  const pv::Interval s2 = pv::sqrt(P(2.0));
  CHECK(s2.lo < s2.hi);
  CHECK(std::nextafter(s2.lo, kInf) == s2.hi);
  CHECK(encloses(s2, std::sqrt(2.0)));
}

// --- Fused multiply-add ------------------------------------------------------

TEST_CASE("fma: exact cases and zero-class shortcuts") {
  check_interval_eq(pv::fma(I(1, 2), I(3, 4), I(5, 6)), I(8, 14));
  check_interval_eq(pv::fma(P(0.0), pv::entire(), I(1, 2)), I(1, 2));
  check_interval_eq(pv::fma(pv::entire(), P(0.0), I(1, 2)), I(1, 2));
  check_interval_eq(pv::fma(I(-2, -1), I(3, 4), P(0.0)), I(-8, -3));
  check_interval_eq(pv::fma(I(1, 2), I(3, 4), I(-kInf, 0.0)), I(-kInf, 8));
}

TEST_CASE("fma: 13x13 class sweep -- no NaN, subset of mul-then-add") {
  const auto reps = sign_class_reps();
  const pv::Interval zs[] = {P(0.0), I(1, 2), I(-2, 3), I(-kInf, -2), I(2, kInf), pv::entire()};
  for (const auto& a : reps) {
    for (const auto& b : reps) {
      for (const auto& z : zs) {
        CAPTURE(a);
        CAPTURE(b);
        CAPTURE(z);
        const pv::Interval r = pv::fma(a, b, z);
        CHECK(well_formed(r));
        CHECK_FALSE(pv::is_empty(r));
        // The fused form rounds once per endpoint, so it must be at least as
        // tight as multiplying and then adding (two outward roundings).
        const pv::Interval two_step = pv::add(pv::mul(a, b), z);
        CHECK_MESSAGE(pv::subset(r, two_step),
                      "fma " << r << " not within mul-then-add " << two_step);
      }
    }
  }
}

// --- Integer powers ----------------------------------------------------------

TEST_CASE("pown: n = 0 gives [1,1] on every nonempty input, including 0 and entire") {
  for (const auto& x :
       {P(0.0), P(5.0), I(1, 2), I(-2, 3), pv::entire(), I(0, kInf), I(-kInf, -2.0)}) {
    CAPTURE(x);
    check_interval_eq(pv::pown(x, 0), P(1.0));
  }
}

TEST_CASE("pown: n = 1 is the identity, bit-for-bit") {
  for (const auto& x : {P(0.0), I(1, 2), I(-2, 3), pv::entire(), I(-kInf, 0.0)}) {
    CAPTURE(x);
    CHECK(bit_equal(pv::pown(x, 1), x));
  }
}

TEST_CASE("pown: positive powers, even via |u| and odd monotone") {
  check_interval_eq(pv::pown(I(-2, 3), 2), I(0, 9));
  check_interval_eq(pv::pown(I(-2, 3), 3), I(-8, 27));
  check_interval_eq(pv::pown(I(2, 3), 2), I(4, 9));
  check_interval_eq(pv::pown(I(-3, -2), 2), I(4, 9));
  check_interval_eq(pv::pown(I(-3, -2), 3), I(-27, -8));
  check_interval_eq(pv::pown(pv::entire(), 2), I(0, kInf));
  check_interval_eq(pv::pown(pv::entire(), 3), pv::entire());
  check_interval_eq(pv::pown(I(-kInf, -2.0), 2), I(4, kInf));
  check_interval_eq(pv::pown(I(-kInf, -2.0), 3), I(-kInf, -8));
  check_interval_eq(pv::pown(P(2.0), 10), P(1024.0));  // dyadic: exact
}

TEST_CASE("pown: negative powers compose with the division semantics") {
  check_interval_eq(pv::pown(I(2, 3), -1), I(pv::detail::div_rd(1.0, 3.0), 0.5));
  check_interval_eq(pv::pown(I(-1, 2), -1), pv::entire());
  check_interval_eq(pv::pown(I(0, 2), -1), I(0.5, kInf));
  check_interval_eq(pv::pown(I(0, 2), -2), I(0.25, kInf));
  check_interval_eq(pv::pown(I(-2, -1), -2), I(0.25, 1));
  check_interval_eq(pv::pown(I(-2, 0.0), -3), I(-kInf, -0.125));
  check_interval_eq(pv::pown(pv::entire(), -1), pv::entire());
}

TEST_CASE("pown: large exponents stay within a few ulps of tightest") {
  // 10^25 is not representable; the enclosure must contain it (via the
  // nearest double 1e25) and be at most a few ulps wide despite the chained
  // directed multiplications of binary exponentiation.
  const pv::Interval p = pv::pown(P(10.0), 25);
  CHECK(p.lo <= 1e25);
  CHECK(1e25 <= p.hi);
  const double ulp = std::nextafter(1e25, kInf) - 1e25;
  CHECK(pv::wid(p) <= 4.0 * ulp);
  // Powers of two are exact at every step: zero width required.
  check_interval_eq(pv::pown(P(2.0), 100), P(0x1p100));
}

// --- abs / min / max ---------------------------------------------------------

TEST_CASE("abs/min/max: exact endpoint selection") {
  check_interval_eq(pv::abs(I(-3, 2)), I(0, 3));
  check_interval_eq(pv::abs(I(-kInf, -2.0)), I(2, kInf));
  check_interval_eq(pv::abs(I(2, 3)), I(2, 3));
  check_interval_eq(pv::abs(I(-3, -2)), I(2, 3));
  check_interval_eq(pv::abs(pv::entire()), I(0, kInf));
  check_interval_eq(pv::abs(P(0.0)), P(0.0));

  check_interval_eq(pv::min(I(1, 4), I(2, 3)), I(1, 3));
  check_interval_eq(pv::max(I(1, 4), I(2, 3)), I(2, 4));
  check_interval_eq(pv::min(I(1, 4), I(-kInf, 2.0)), I(-kInf, 2.0));
  check_interval_eq(pv::min(I(1, 4), I(2, kInf)), I(1, 4));
  check_interval_eq(pv::max(I(1, 4), I(2, kInf)), I(2, kInf));
  check_interval_eq(pv::max(I(-3, -1), I(-2, 0.0)), I(-2, 0.0));
}

// --- floor / ceil / trunc / round / sign -------------------------------------

TEST_CASE("floor/ceil/trunc/round/sign: exact step functions on both endpoints") {
  check_interval_eq(pv::floor(I(-1.5, 2.5)), I(-2, 2));
  check_interval_eq(pv::ceil(I(-1.5, 2.5)), I(-1, 3));
  check_interval_eq(pv::trunc(I(-1.7, 2.7)), I(-1, 2));
  check_interval_eq(pv::floor(I(2, 3)), I(2, 3));  // integers are fixed points

  check_interval_eq(pv::round_ties_to_even(I(0.5, 2.5)), I(0, 2));
  check_interval_eq(pv::round_ties_to_away(I(0.5, 2.5)), I(1, 3));
  check_interval_eq(pv::round_ties_to_even(I(-2.5, -0.5)), I(-2, 0.0));
  check_interval_eq(pv::round_ties_to_away(I(-2.5, -0.5)), I(-3, -1));
  check_interval_eq(pv::round_ties_to_even(I(1.5, 1.5)), P(2.0));

  check_interval_eq(pv::sign(I(-3, 2)), I(-1, 1));
  check_interval_eq(pv::sign(P(0.0)), P(0.0));
  check_interval_eq(pv::sign(I(2, 5)), P(1.0));
  check_interval_eq(pv::sign(I(-5, -2)), P(-1.0));
  check_interval_eq(pv::sign(I(0.0, 3)), I(0.0, 1));
  check_interval_eq(pv::sign(I(-3, 0.0)), I(-1, 0.0));
  check_interval_eq(pv::sign(pv::entire()), I(-1, 1));
}

TEST_CASE("floor/ceil/trunc/round: infinite endpoints pass through") {
  check_interval_eq(pv::floor(I(-kInf, 1.5)), I(-kInf, 1));
  check_interval_eq(pv::ceil(I(1.5, kInf)), I(2, kInf));
  check_interval_eq(pv::trunc(I(-kInf, -2.5)), I(-kInf, -2));
  check_interval_eq(pv::trunc(pv::entire()), pv::entire());
  check_interval_eq(pv::round_ties_to_even(I(-kInf, 0.5)), I(-kInf, 0.0));
  check_interval_eq(pv::round_ties_to_away(I(0.5, kInf)), I(1, kInf));
}

// --- Property tests (randomized, fixed seeds) --------------------------------

// (1) Containment: for members x in X, y in Y the correctly rounded pointwise
// result RN(f(x, y)) lies inside the interval result. RN never rounds past a
// double bound of the exact value, so this follows from -- and thus verifies
// -- the outward-rounding enclosure at every endpoint pairing.
TEST_CASE("property: pointwise RN results are members of interval results") {
  std::mt19937_64 rng(0x17882015ULL);  // fixed seed: CI-reproducible
  constexpr int kIters = 20000;
  int failures = 0;
  auto expect = [&](bool ok, const char* op, const pv::Interval& r, double x, double y, double z) {
    if (!ok && ++failures <= 8) {
      FAIL_CHECK("containment violated for " << op << ": x=" << x << " y=" << y << " f=" << z
                                             << " result=" << r);
    }
  };
  std::vector<double> xs, ys, zs;
  for (int it = 0; it < kIters; ++it) {
    const pv::Interval x = random_interval(rng);
    const pv::Interval y = random_interval(rng);
    sample_members(x, rng, xs);
    sample_members(y, rng, ys);

    const pv::Interval sum = pv::add(x, y);
    const pv::Interval dif = pv::sub(x, y);
    const pv::Interval prd = pv::mul(x, y);
    const pv::Interval quo = pv::div(x, y);
    for (double u : xs) {
      for (double v : ys) {
        expect(encloses(sum, u + v), "add", sum, u, v, u + v);
        expect(encloses(dif, u - v), "sub", dif, u, v, u - v);
        expect(encloses(prd, u * v), "mul", prd, u, v, u * v);
        if (v != 0.0) expect(encloses(quo, u / v), "div", quo, u, v, u / v);
      }
    }

    // fma over a third interval (limited fan-out keeps the runtime modest).
    const pv::Interval zi = random_interval(rng);
    sample_members(zi, rng, zs);
    const pv::Interval fr = pv::fma(x, y, zi);
    const std::size_t nx = std::min<std::size_t>(xs.size(), 3);
    const std::size_t ny = std::min<std::size_t>(ys.size(), 3);
    const std::size_t nz = std::min<std::size_t>(zs.size(), 3);
    for (std::size_t i = 0; i < nx; ++i) {
      for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t k = 0; k < nz; ++k) {
          const double w = std::fma(xs[i], ys[j], zs[k]);
          expect(encloses(fr, w), "fma", fr, xs[i], ys[j], w);
        }
      }
    }

    // Unary: sqr and sqrt via correctly rounded pointwise results; pown via
    // directed pointwise brackets [RD(u^n), RU(u^n)] (std::pow is not
    // correctly rounded, so it cannot serve as a reference).
    const pv::Interval sq = pv::sqr(x);
    const pv::Interval rt = pv::sqrt(x);
    const pv::Interval p2 = pv::pown(x, 2);
    const pv::Interval p3 = pv::pown(x, 3);
    const pv::Interval pm1 = pv::pown(x, -1);
    for (double u : xs) {
      expect(encloses(sq, u * u), "sqr", sq, u, 0.0, u * u);
      if (u >= 0.0) expect(encloses(rt, std::sqrt(u)), "sqrt", rt, u, 0.0, std::sqrt(u));
      for (int n : {2, 3}) {
        const pv::Interval& pn = n == 2 ? p2 : p3;
        const double d = pv::detail::pown_scalar(u, n, /*round_up=*/false);
        const double h = pv::detail::pown_scalar(u, n, /*round_up=*/true);
        expect(encloses(pn, d) && encloses(pn, h), n == 2 ? "pown2" : "pown3", pn, u, 0.0, d);
      }
      if (u != 0.0) {
        const double d = pv::detail::div_rd(1.0, u);
        const double h = pv::detail::div_ru(1.0, u);
        expect(encloses(pm1, d) && encloses(pm1, h), "pown-1", pm1, u, 0.0, d);
      }
    }
  }
  CHECK(failures == 0);
}

// (2) Inclusion isotonicity: X' subset X and Y' subset Y imply
// op(X', Y') subset op(X, Y) -- the defining monotonicity of set-valued
// extensions, which any valid enclosure scheme must preserve.
TEST_CASE("property: inclusion isotonicity under random shrinking") {
  std::mt19937_64 rng(0xC0FFEE1788ULL);  // fixed seed: CI-reproducible
  constexpr int kIters = 20000;
  int failures = 0;
  auto expect = [&](bool ok, const char* op, const pv::Interval& small, const pv::Interval& big) {
    if (!ok && ++failures <= 8) {
      FAIL_CHECK("isotonicity violated for " << op << ": " << small << " not within " << big);
    }
  };
  for (int it = 0; it < kIters; ++it) {
    const pv::Interval x = random_interval(rng);
    const pv::Interval y = random_interval(rng);
    const pv::Interval z = random_interval(rng);
    const pv::Interval xs = shrink(x, rng);
    const pv::Interval ys = shrink(y, rng);
    const pv::Interval zs = shrink(z, rng);
    REQUIRE(pv::subset(xs, x));
    REQUIRE(pv::subset(ys, y));

    expect(pv::subset(pv::add(xs, ys), pv::add(x, y)), "add", pv::add(xs, ys), pv::add(x, y));
    expect(pv::subset(pv::sub(xs, ys), pv::sub(x, y)), "sub", pv::sub(xs, ys), pv::sub(x, y));
    expect(pv::subset(pv::mul(xs, ys), pv::mul(x, y)), "mul", pv::mul(xs, ys), pv::mul(x, y));
    expect(pv::subset(pv::div(xs, ys), pv::div(x, y)), "div", pv::div(xs, ys), pv::div(x, y));
    expect(pv::subset(pv::min(xs, ys), pv::min(x, y)), "min", pv::min(xs, ys), pv::min(x, y));
    expect(pv::subset(pv::max(xs, ys), pv::max(x, y)), "max", pv::max(xs, ys), pv::max(x, y));
    expect(pv::subset(pv::fma(xs, ys, zs), pv::fma(x, y, z)), "fma", pv::fma(xs, ys, zs),
           pv::fma(x, y, z));
    expect(pv::subset(pv::neg(xs), pv::neg(x)), "neg", pv::neg(xs), pv::neg(x));
    expect(pv::subset(pv::sqr(xs), pv::sqr(x)), "sqr", pv::sqr(xs), pv::sqr(x));
    expect(pv::subset(pv::sqrt(xs), pv::sqrt(x)), "sqrt", pv::sqrt(xs), pv::sqrt(x));
    expect(pv::subset(pv::abs(xs), pv::abs(x)), "abs", pv::abs(xs), pv::abs(x));
    expect(pv::subset(pv::recip(xs), pv::recip(x)), "recip", pv::recip(xs), pv::recip(x));
    for (int n : {2, 3, -2}) {
      expect(pv::subset(pv::pown(xs, n), pv::pown(x, n)), "pown", pv::pown(xs, n), pv::pown(x, n));
    }
  }
  CHECK(failures == 0);
}

// (3) Algebraic identities that must hold bit-for-bit: negation is an exact
// involution (endpoint negation incurs no rounding), and subtraction is
// defined as addition of the negation, so both spellings must agree exactly.
TEST_CASE("property: neg involution and sub == add-of-negation, bit-identical") {
  std::mt19937_64 rng(0x20260808ULL);  // fixed seed: CI-reproducible
  constexpr int kIters = 20000;
  int failures = 0;

  // Deterministic seed cases first: class representatives and empty.
  std::vector<pv::Interval> cases = sign_class_reps();
  cases.push_back(pv::empty());
  for (const auto& x : cases) {
    CAPTURE(x);
    CHECK(bit_equal(pv::neg(pv::neg(x)), x));
    for (const auto& y : cases) {
      // sub and add-of-negation must denote the same set (they are in fact
      // bit-identical by construction, but set-equality is the rigor invariant
      // and is robust to benign signed-zero representation differences).
      CHECK(pv::equal(pv::sub(x, y), pv::add(x, pv::neg(y))));
    }
  }

  for (int it = 0; it < kIters; ++it) {
    const pv::Interval x = random_interval(rng);
    const pv::Interval y = random_interval(rng);
    if (!bit_equal(pv::neg(pv::neg(x)), x) && ++failures <= 8) {
      FAIL_CHECK("neg(neg(X)) != X for X=" << x);
    }
    if (!pv::equal(pv::sub(x, y), pv::add(x, pv::neg(y))) && ++failures <= 8) {
      FAIL_CHECK("sub(X,Y) != add(X,neg(Y)) for X=" << x << " Y=" << y);
    }
  }
  CHECK(failures == 0);
}
