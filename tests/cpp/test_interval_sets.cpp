// Milestone 2: exhaustive tests for the non-arithmetic interval kernel.
//
// Covers: factories/queries (make, point, valid_bounds, is_empty/is_entire/
// is_common/is_singleton/is_member/contains_zero), set operations
// (intersection, convex_hull), numeric functions (mid, rad, wid, mag, mig),
// all eight comparison predicates over a fixed pair table, cancel_minus /
// cancel_plus, and the reverse operations (mul_rev_to_pair, mul_rev, sqr_rev,
// abs_rev), plus a smoke test of the diagnostic printer.
//
// Semantics under test are IEEE 1788-2015 set-based conventions: the empty
// set absorbs, entire is the "no such z" value of cancelMinus, infinities are
// directions (never members), inf(empty) = +inf / sup(empty) = -inf, and
// every returned enclosure is a superset of the exact set. Reverse-operation
// tables are cross-validated against the *defining* set by exact dyadic
// sampling, so the tests do not simply mirror the implementation. All random
// property tests use fixed mt19937_64 seeds (the generator is fully specified
// by the standard, so runs are reproducible across platforms and CI).

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <utility>

#include "doctest/doctest.h"
#include "pyintval/interval.hpp"

namespace {

namespace pd = pyintval::detail;
using pyintval::Interval;

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kMax = std::numeric_limits<double>::max();
constexpr double kMinNormal = std::numeric_limits<double>::min();
constexpr double kDen = std::numeric_limits<double>::denorm_min();

Interval I(double lo, double hi) { return pyintval::make(lo, hi); }

const Interval kEmp = pyintval::empty();
const Interval kEnt = pyintval::entire();

// Bitwise identity (distinguishes -0.0 from +0.0, unlike equal()).
bool bit_identical(const Interval& x, const Interval& y) {
  return std::bit_cast<std::uint64_t>(x.lo) == std::bit_cast<std::uint64_t>(y.lo) &&
         std::bit_cast<std::uint64_t>(x.hi) == std::bit_cast<std::uint64_t>(y.hi);
}

// Uniform over finite bit patterns: exercises the full exponent range,
// including subnormals (rarely) and both zeros.
double rand_finite_bits(std::mt19937_64& rng) {
  for (;;) {
    const double v = std::bit_cast<double>(rng());
    if (std::isfinite(v)) return v;
  }
}

// Mixed-regime generator: full-range doubles, forced subnormals/zeros, small
// dyadics, and small integers -- the regimes where midpoint/width/cancel
// rounding traps live.
double rand_mixed(std::mt19937_64& rng) {
  switch (rng() & 3u) {
    case 0:
      return rand_finite_bits(rng);
    case 1:
      // Clear the exponent field: a subnormal (or zero) of either sign.
      return std::bit_cast<double>(rng() & 0x800FFFFFFFFFFFFFULL);
    case 2: {
      // Dyadic values in [-4, 4): exercise the "ordinary" magnitude range.
      const std::int64_t n = static_cast<std::int64_t>(rng() >> 32) - (std::int64_t{1} << 31);
      return static_cast<double>(n) * 0x1p-29;
    }
    default:
      return static_cast<double>(static_cast<std::int64_t>(rng() % 2001u) - 1000);
  }
}

// Random nonempty interval; optionally replaces endpoints with infinite
// directions (probability ~1/32 per side) to cover half-lines and entire.
Interval rand_interval(std::mt19937_64& rng, bool allow_unbounded) {
  double x = rand_mixed(rng);
  double y = rand_mixed(rng);
  if (x > y) std::swap(x, y);
  if (allow_unbounded) {
    const std::uint64_t sel = rng();
    if ((sel & 0x1Fu) == 0) x = -kInf;
    if (((sel >> 5) & 0x1Fu) == 0) y = kInf;
  }
  return pyintval::make(x, y);
}

// --- mul_rev_to_pair cross-validation helpers -------------------------------
//
// For a fixed x, {v*x : v in b} is the closed interval spanned by the two
// endpoint products (the map v -> v*x is linear in v). On the dyadic sample
// grid below every product is exact in binary64, so intersecting that range
// with c decides "exists v in b with v*x in c" with no rounding error: an
// oracle independent of the kernel's division-based construction.
bool mul_rev_has_solution(double x, const Interval& b, const Interval& c) {
  const double p = b.lo * x;
  const double q = b.hi * x;
  const double lo = std::min(p, q);
  const double hi = std::max(p, q);
  return hi >= c.lo && c.hi >= lo;
}

// One-ulp inward shrink, for the tightness direction of the sampling check
// (the kernel may legitimately overshoot each boundary by one outward ulp).
Interval shrink_one_ulp(const Interval& p) {
  if (pyintval::is_empty(p)) return p;
  return pyintval::make(pd::succ(p.lo), pd::pred(p.hi));
}

struct MulRevRow {
  const char* name;
  Interval b, c, p1, p2;
};

// Expected values derived from the defining set {x : exists v in b, v*x in c}
// (all division boundaries below are exactly representable, so the kernel's
// outward rounding must reproduce them exactly).
const MulRevRow kMulRevRows[] = {
    {"b=[-1,1] c=[1,2]", I(-1, 1), I(1, 2), I(-kInf, -1), I(1, kInf)},
    {"b=[-1,1] c=[-2,-1]", I(-1, 1), I(-2, -1), I(-kInf, -1), I(1, kInf)},
    {"b=[-2,3] c=[6,6]", I(-2, 3), I(6, 6), I(-kInf, -3), I(2, kInf)},
    {"b=[0,0] c=[-1,1]", I(0, 0), I(-1, 1), kEnt, kEmp},
    {"b=[0,0] c=[0,0]", I(0, 0), I(0, 0), kEnt, kEmp},
    {"b=[0,0] c=[1,2]", I(0, 0), I(1, 2), kEmp, kEmp},
    {"b=[2,4] c=[6,8]", I(2, 4), I(6, 8), I(1.5, 4), kEmp},
    {"b=[2,4] c=[0,0]", I(2, 4), I(0, 0), I(0, 0), kEmp},
    {"b=[2,4] c=[-1,1]", I(2, 4), I(-1, 1), I(-0.5, 0.5), kEmp},
    {"b=[1,2] c=[-4,-2]", I(1, 2), I(-4, -2), I(-4, -1), kEmp},
    {"b=[-2,-1] c=[2,4]", I(-2, -1), I(2, 4), I(-4, -1), kEmp},
    {"b=[0,1] c=[1,2]", I(0, 1), I(1, 2), I(1, kInf), kEmp},
    {"b=[-1,0] c=[1,2]", I(-1, 0), I(1, 2), I(-kInf, -1), kEmp},
    {"b=[0,2] c=[-4,-1]", I(0, 2), I(-4, -1), I(-kInf, -0.5), kEmp},
    {"b=[-2,0] c=[1,2]", I(-2, 0), I(1, 2), I(-kInf, -0.5), kEmp},
    {"b=[-2,0] c=[-4,-1]", I(-2, 0), I(-4, -1), I(0.5, kInf), kEmp},
    {"b=[-1,1] c=[0,2]", I(-1, 1), I(0, 2), kEnt, kEmp},
};

}  // namespace

// --- Factories and validity --------------------------------------------------

TEST_CASE("make/point: invalid bounds yield empty, canonical empty per 1788 9.2") {
  CHECK(pyintval::is_empty(I(2, 1)));
  CHECK(pyintval::is_empty(I(kNaN, 1)));
  CHECK(pyintval::is_empty(I(1, kNaN)));
  CHECK(pyintval::is_empty(I(kNaN, kNaN)));
  CHECK(pyintval::is_empty(I(kInf, kInf)));  // infinity is never a member
  CHECK(pyintval::is_empty(I(-kInf, -kInf)));
  CHECK(pyintval::is_empty(I(kInf, -kInf)));
  CHECK(pyintval::is_entire(I(-kInf, kInf)));
  CHECK(!pyintval::is_empty(I(-kInf, 3)));  // half-lines are valid
  CHECK(!pyintval::is_empty(I(3, kInf)));

  const Interval p = pyintval::point(3.0);
  CHECK(p.lo == 3.0);
  CHECK(p.hi == 3.0);
  CHECK(pyintval::is_singleton(p));
  CHECK(pyintval::is_empty(pyintval::point(kInf)));
  CHECK(pyintval::is_empty(pyintval::point(-kInf)));
  CHECK(pyintval::is_empty(pyintval::point(kNaN)));

  // Canonical empty {+inf, -inf}: the 1788 inf/sup conventions fall out.
  CHECK(pyintval::inf(kEmp) == kInf);
  CHECK(pyintval::sup(kEmp) == -kInf);
  CHECK(pyintval::inf(I(2, 1)) == kInf);
  CHECK(pyintval::sup(I(kNaN, kNaN)) == -kInf);
}

TEST_CASE("valid_bounds mirrors the factory rules") {
  CHECK(pyintval::valid_bounds(1.0, 2.0));
  CHECK(pyintval::valid_bounds(1.0, 1.0));
  CHECK(pyintval::valid_bounds(-kInf, kInf));
  CHECK(pyintval::valid_bounds(-kInf, 5.0));
  CHECK(pyintval::valid_bounds(5.0, kInf));
  CHECK(pyintval::valid_bounds(-0.0, 0.0));
  CHECK(pyintval::valid_bounds(0.0, -0.0));  // 0 == -0, so lo <= hi holds
  CHECK(!pyintval::valid_bounds(2.0, 1.0));
  CHECK(!pyintval::valid_bounds(kNaN, 1.0));
  CHECK(!pyintval::valid_bounds(1.0, kNaN));
  CHECK(!pyintval::valid_bounds(kNaN, kNaN));
  CHECK(!pyintval::valid_bounds(kInf, kInf));
  CHECK(!pyintval::valid_bounds(-kInf, -kInf));
  CHECK(!pyintval::valid_bounds(kInf, -kInf));
}

TEST_CASE("zero endpoints canonicalize to +0.0") {
  const Interval a = I(-0.0, 0.0);
  CHECK(!std::signbit(a.lo));
  CHECK(!std::signbit(a.hi));
  const Interval b = I(-0.0, -0.0);
  CHECK(!std::signbit(b.lo));
  CHECK(!std::signbit(b.hi));
  const Interval c = I(0.0, -0.0);
  CHECK(!pyintval::is_empty(c));
  CHECK(!std::signbit(c.lo));
  CHECK(!std::signbit(c.hi));
  const Interval d = I(-1.0, -0.0);
  CHECK(d.lo == -1.0);
  CHECK(!std::signbit(d.hi));
  const Interval e = pyintval::point(-0.0);
  CHECK(!std::signbit(e.lo));
  CHECK(!std::signbit(e.hi));
  // [-0, 0], [0, -0], and [0, 0] denote the same set and share one bit pattern.
  CHECK(bit_identical(a, I(0.0, 0.0)));
  CHECK(bit_identical(c, I(0.0, 0.0)));
}

// --- Basic queries -----------------------------------------------------------

TEST_CASE("is_empty / is_entire / is_common / is_singleton") {
  CHECK(pyintval::is_empty(kEmp));
  CHECK(!pyintval::is_empty(I(1, 2)));
  CHECK(!pyintval::is_empty(kEnt));

  CHECK(pyintval::is_entire(kEnt));
  CHECK(!pyintval::is_entire(kEmp));
  CHECK(!pyintval::is_entire(I(1, kInf)));
  CHECK(!pyintval::is_entire(I(-kInf, 1)));
  CHECK(!pyintval::is_entire(I(-kMax, kMax)));

  // "common" = nonempty and bounded.
  CHECK(pyintval::is_common(I(1, 2)));
  CHECK(pyintval::is_common(I(-kMax, kMax)));
  CHECK(pyintval::is_common(pyintval::point(0.0)));
  CHECK(!pyintval::is_common(kEnt));
  CHECK(!pyintval::is_common(I(1, kInf)));
  CHECK(!pyintval::is_common(I(-kInf, 1)));
  CHECK(!pyintval::is_common(kEmp));

  CHECK(pyintval::is_singleton(pyintval::point(3.0)));
  CHECK(pyintval::is_singleton(pyintval::point(0.0)));
  CHECK(pyintval::is_singleton(I(kMax, kMax)));
  CHECK(!pyintval::is_singleton(kEmp));
  CHECK(!pyintval::is_singleton(I(1, 2)));
  CHECK(!pyintval::is_singleton(kEnt));
}

TEST_CASE("is_member and contains_zero") {
  CHECK(pyintval::is_member(3.0, I(2, 4)));
  CHECK(pyintval::is_member(2.0, I(2, 4)));  // endpoints are members
  CHECK(pyintval::is_member(4.0, I(2, 4)));
  CHECK(!pyintval::is_member(4.5, I(2, 4)));
  CHECK(!pyintval::is_member(1.5, I(2, 4)));
  // Infinity is a direction, never a member -- even of its own half-line.
  CHECK(!pyintval::is_member(kInf, I(2, kInf)));
  CHECK(!pyintval::is_member(-kInf, I(-kInf, 2)));
  CHECK(!pyintval::is_member(kInf, kEnt));
  CHECK(!pyintval::is_member(-kInf, kEnt));
  CHECK(!pyintval::is_member(kNaN, kEnt));
  CHECK(!pyintval::is_member(0.0, kEmp));
  CHECK(pyintval::is_member(kMax, I(2, kInf)));
  // Signed zeros denote the same real number.
  CHECK(pyintval::is_member(-0.0, I(0, 1)));
  CHECK(pyintval::is_member(0.0, I(-1, -0.0)));
  CHECK(pyintval::is_member(kDen, I(0, kMinNormal)));

  CHECK(pyintval::contains_zero(I(-1, 1)));
  CHECK(pyintval::contains_zero(I(0, 2)));
  CHECK(pyintval::contains_zero(I(-2, -0.0)));
  CHECK(pyintval::contains_zero(pyintval::point(0.0)));
  CHECK(pyintval::contains_zero(kEnt));
  CHECK(!pyintval::contains_zero(I(1, 2)));
  CHECK(!pyintval::contains_zero(I(-2, -1)));
  CHECK(!pyintval::contains_zero(I(kDen, 1)));
  CHECK(!pyintval::contains_zero(kEmp));
}

// --- Set operations ----------------------------------------------------------

TEST_CASE("intersection: standard, touching, disjoint, empty-absorbing") {
  CHECK(pyintval::equal(pyintval::intersection(I(1, 3), I(2, 4)), I(2, 3)));
  CHECK(pyintval::equal(pyintval::intersection(I(2, 4), I(1, 3)), I(2, 3)));
  // Touching intervals meet in a single point.
  CHECK(pyintval::equal(pyintval::intersection(I(1, 2), I(2, 3)), I(2, 2)));
  CHECK(pyintval::is_singleton(pyintval::intersection(I(1, 2), I(2, 3))));
  CHECK(pyintval::is_empty(pyintval::intersection(I(1, 2), I(3, 4))));
  CHECK(pyintval::is_empty(pyintval::intersection(kEmp, I(1, 2))));
  CHECK(pyintval::is_empty(pyintval::intersection(I(1, 2), kEmp)));
  CHECK(pyintval::is_empty(pyintval::intersection(kEmp, kEmp)));
  CHECK(pyintval::is_empty(pyintval::intersection(kEmp, kEnt)));
  CHECK(pyintval::equal(pyintval::intersection(kEnt, I(1, 2)), I(1, 2)));
  CHECK(pyintval::equal(pyintval::intersection(I(1, 2), kEnt), I(1, 2)));
  CHECK(pyintval::equal(pyintval::intersection(I(-kInf, 2), I(0, kInf)), I(0, 2)));
  CHECK(pyintval::equal(pyintval::intersection(I(1, 2), I(0, 3)), I(1, 2)));
  // Zero canonicalization survives the meet.
  const Interval z = pyintval::intersection(I(-1, -0.0), I(0, 5));
  CHECK(pyintval::equal(z, I(0, 0)));
  CHECK(!std::signbit(z.lo));
  CHECK(!std::signbit(z.hi));
}

TEST_CASE("convex_hull: identity on empty, spans gaps of disjoint operands") {
  CHECK(pyintval::is_empty(pyintval::convex_hull(kEmp, kEmp)));
  CHECK(pyintval::equal(pyintval::convex_hull(kEmp, I(1, 2)), I(1, 2)));
  CHECK(pyintval::equal(pyintval::convex_hull(I(1, 2), kEmp), I(1, 2)));
  CHECK(pyintval::equal(pyintval::convex_hull(I(1, 2), I(3, 4)), I(1, 4)));
  CHECK(pyintval::equal(pyintval::convex_hull(I(3, 4), I(1, 2)), I(1, 4)));
  CHECK(pyintval::equal(pyintval::convex_hull(I(1, 3), I(2, 4)), I(1, 4)));
  CHECK(pyintval::equal(pyintval::convex_hull(I(1, 2), I(1.5, 1.75)), I(1, 2)));
  CHECK(pyintval::is_entire(pyintval::convex_hull(I(-kInf, 0), I(1, kInf))));
  CHECK(pyintval::is_entire(pyintval::convex_hull(kEnt, I(1, 2))));
  CHECK(pyintval::equal(pyintval::convex_hull(kEmp, kEnt), kEnt));
  // Hull is an upper bound of both operands.
  CHECK(pyintval::subset(I(1, 2), pyintval::convex_hull(I(1, 2), I(3, 4))));
  CHECK(pyintval::subset(I(3, 4), pyintval::convex_hull(I(1, 2), I(3, 4))));
}

// --- Numeric functions -------------------------------------------------------

TEST_CASE("mid: special values and one-ulp intervals") {
  CHECK(std::isnan(pyintval::mid(kEmp)));
  CHECK(pyintval::mid(kEnt) == 0.0);
  CHECK(pyintval::mid(I(-kInf, 2)) == -kMax);
  CHECK(pyintval::mid(I(2, kInf)) == kMax);
  CHECK(pyintval::mid(I(1, 2)) == 1.5);
  CHECK(pyintval::mid(I(kMax, kMax)) == kMax);  // lo+hi overflows; halves are exact
  CHECK(pyintval::mid(I(-kMax, kMax)) == 0.0);
  CHECK(pyintval::mid(I(-kDen, kDen)) == 0.0);
  CHECK(pyintval::mid(pyintval::point(-7.25)) == -7.25);

  // mid of [x, succ(x)] must land on one of the two endpoints (there is no
  // double strictly between), hence be a member.
  const double specials[] = {1.0,        -1.0,  0.0,  kDen,  -2.0 * kDen,
                             kMinNormal, 1e300, -3.5, -kMax, kMax};
  for (const double x : specials) {
    const Interval v = pyintval::make(x, pd::succ(x));
    const double m = pyintval::mid(v);
    CAPTURE(v);
    CAPTURE(m);
    CHECK((m == v.lo || m == v.hi));
    CHECK(pyintval::is_member(m, v));
  }
}

TEST_CASE("mid is a finite member of every nonempty interval (property, 20000 cases)") {
  std::mt19937_64 rng(0x17882015DEADBEEFULL);
  int failures = 0;
  for (int i = 0; i < 20000; ++i) {
    const Interval v = rand_interval(rng, /*allow_unbounded=*/true);
    const double m = pyintval::mid(v);
    const bool ok = std::isfinite(m) && pyintval::is_member(m, v);
    if (!ok) {
      ++failures;
      if (failures <= 3) {
        CAPTURE(i);
        CAPTURE(v);
        CAPTURE(m);
        CHECK(ok);
      }
    }
  }
  CHECK(failures == 0);
}

TEST_CASE("rad: special values") {
  CHECK(std::isnan(pyintval::rad(kEmp)));
  CHECK(pyintval::rad(I(1, 2)) == 0.5);
  CHECK(pyintval::rad(pyintval::point(3.0)) == 0.0);
  CHECK(pyintval::rad(pyintval::point(0.0)) == 0.0);
  CHECK(pyintval::rad(I(2, kInf)) == kInf);
  CHECK(pyintval::rad(I(-kInf, 2)) == kInf);
  CHECK(pyintval::rad(kEnt) == kInf);
  CHECK(pyintval::rad(I(1.0, pd::succ(1.0))) == 0x1p-52);
  CHECK(pyintval::rad(I(-kMax, kMax)) == kMax);
}

TEST_CASE("rad: [mid-rad, mid+rad] encloses the interval (property, 10000 cases)") {
  // For doubles m, r, lo, hi: the real inequality m - r <= lo is equivalent
  // to RU(m - r) <= lo (lo is itself a double), and hi <= m + r to
  // RD(m + r) >= hi. Testing with the directed primitives is therefore exact
  // and cannot mask a violation behind test-side rounding.
  std::mt19937_64 rng(0x52AD1CA1F00DULL);
  int failures = 0;
  for (int i = 0; i < 10000; ++i) {
    const Interval v = rand_interval(rng, /*allow_unbounded=*/true);
    const double m = pyintval::mid(v);
    const double r = pyintval::rad(v);
    const bool ok = r >= 0.0 && pd::sub_ru(m, r) <= v.lo && pd::add_rd(m, r) >= v.hi;
    if (!ok) {
      ++failures;
      if (failures <= 3) {
        CAPTURE(i);
        CAPTURE(v);
        CAPTURE(m);
        CAPTURE(r);
        CHECK(ok);
      }
    }
  }
  CHECK(failures == 0);
}

TEST_CASE("wid: special values, exactness at one ulp and subnormals") {
  CHECK(std::isnan(pyintval::wid(kEmp)));
  CHECK(pyintval::wid(I(1, 2)) == 1.0);
  CHECK(pyintval::wid(pyintval::point(7.0)) == 0.0);
  CHECK(pyintval::wid(I(1.0, pd::succ(1.0))) == 0x1p-52);
  CHECK(pyintval::wid(I(-pd::succ(1.0), -1.0)) == 0x1p-52);
  CHECK(pyintval::wid(kEnt) == kInf);
  CHECK(pyintval::wid(I(-kInf, 0)) == kInf);
  CHECK(pyintval::wid(I(0, kInf)) == kInf);
  CHECK(pyintval::wid(I(-kMax, kMax)) == kInf);  // overflows; RU gives +inf
  // Near subnormals the width must stay exact / strictly positive.
  CHECK(pyintval::wid(I(0, kDen)) == kDen);
  CHECK(pyintval::wid(I(-kDen, 0)) == kDen);
  CHECK(pyintval::wid(I(kDen, 2 * kDen)) == kDen);
  CHECK(pyintval::wid(I(-kDen, kDen)) == 2 * kDen);
}

TEST_CASE("wid of [x, succ(x)] is strictly positive (property, 5000 cases)") {
  std::mt19937_64 rng(0x51DEC0DE1234ULL);
  int failures = 0;
  for (int i = 0; i < 5000; ++i) {
    const double x = rand_mixed(rng);
    const Interval v = pyintval::make(x, pd::succ(x));
    const double w = pyintval::wid(v);
    const bool ok = w > 0.0;
    if (!ok) {
      ++failures;
      if (failures <= 3) {
        CAPTURE(i);
        CAPTURE(v);
        CAPTURE(w);
        CHECK(ok);
      }
    }
  }
  CHECK(failures == 0);
}

TEST_CASE("mag and mig") {
  CHECK(std::isnan(pyintval::mag(kEmp)));
  CHECK(pyintval::mag(I(-3, 2)) == 3.0);
  CHECK(pyintval::mag(I(-2, 3)) == 3.0);
  CHECK(pyintval::mag(kEnt) == kInf);
  CHECK(pyintval::mag(I(-1, kInf)) == kInf);
  CHECK(pyintval::mag(I(-kInf, 1)) == kInf);
  CHECK(pyintval::mag(pyintval::point(0.0)) == 0.0);
  CHECK(pyintval::mag(I(-0.0, 0.0)) == 0.0);
  CHECK(pyintval::mag(I(2, 5)) == 5.0);

  CHECK(std::isnan(pyintval::mig(kEmp)));
  CHECK(pyintval::mig(I(-3, 2)) == 0.0);  // straddles zero
  CHECK(pyintval::mig(I(-3, -2)) == 2.0);
  CHECK(pyintval::mig(I(2, 3)) == 2.0);
  CHECK(pyintval::mig(I(2, kInf)) == 2.0);
  CHECK(pyintval::mig(I(-kInf, -2)) == 2.0);
  CHECK(pyintval::mig(I(0, 5)) == 0.0);
  CHECK(pyintval::mig(I(-5, -0.0)) == 0.0);
  CHECK(pyintval::mig(kEnt) == 0.0);
}

// --- Comparison predicates ---------------------------------------------------

TEST_CASE("comparison predicates: full truth table over the reference pair set") {
  // Columns: equal, subset, interior, disjoint, less, strict_less, precedes,
  // strict_precedes. Conventions: empty is subset/interior of everything;
  // less/strict_less are true on (E,E) only among empty pairs;
  // precedes/strict_precedes are true whenever either side is empty; matching
  // infinities count as interior / strictly-less on their unbounded side.
  struct Row {
    const char* name;
    Interval a, b;
    bool eq, sub, itr, dis, ls, sls, prec, sprec;
  };
  const Row rows[] = {
      // name                     a           b        eq     sub    itr    dis    ls     sls prec
      // sprec
      {"(E,E)", kEmp, kEmp, true, true, true, true, true, true, true, true},
      {"(E,[1,2])", kEmp, I(1, 2), false, true, true, true, false, false, true, true},
      {"([1,2],E)", I(1, 2), kEmp, false, false, false, true, false, false, true, true},
      {"(E,R)", kEmp, kEnt, false, true, true, true, false, false, true, true},
      {"([1,2],[1,2])", I(1, 2), I(1, 2), true, true, false, false, true, false, false, false},
      {"([1,2],[0,3])", I(1, 2), I(0, 3), false, true, true, false, false, false, false, false},
      {"([0,3],[1,2])", I(0, 3), I(1, 2), false, false, false, false, false, false, false, false},
      {"([1,2],[2,3])", I(1, 2), I(2, 3), false, false, false, false, true, true, true, false},
      {"([1,2],[3,4])", I(1, 2), I(3, 4), false, false, false, true, true, true, true, true},
      {"([1,2],R)", I(1, 2), kEnt, false, true, true, false, false, false, false, false},
      {"(R,R)", kEnt, kEnt, true, true, true, false, true, true, false, false},
      {"([-inf,2],[-inf,3])", I(-kInf, 2), I(-kInf, 3), false, true, true, false, true, true, false,
       false},
      {"([1,inf],[0,inf])", I(1, kInf), I(0, kInf), false, true, true, false, false, false, false,
       false},
      {"([1,2],[1,3])", I(1, 2), I(1, 3), false, true, false, false, true, false, false, false},
      {"([-inf,2],R)", I(-kInf, 2), kEnt, false, true, true, false, true, true, false, false},
  };
  for (const Row& r : rows) {
    CAPTURE(r.name);
    CHECK(pyintval::equal(r.a, r.b) == r.eq);
    CHECK(pyintval::subset(r.a, r.b) == r.sub);
    CHECK(pyintval::interior(r.a, r.b) == r.itr);
    CHECK(pyintval::disjoint(r.a, r.b) == r.dis);
    CHECK(pyintval::less(r.a, r.b) == r.ls);
    CHECK(pyintval::strict_less(r.a, r.b) == r.sls);
    CHECK(pyintval::precedes(r.a, r.b) == r.prec);
    CHECK(pyintval::strict_precedes(r.a, r.b) == r.sprec);
    // Structural implications that must hold for every pair.
    CHECK(pyintval::disjoint(r.a, r.b) == pyintval::disjoint(r.b, r.a));
    CHECK(pyintval::equal(r.a, r.b) == pyintval::equal(r.b, r.a));
    if (pyintval::interior(r.a, r.b)) CHECK(pyintval::subset(r.a, r.b));
    if (pyintval::strict_less(r.a, r.b)) CHECK(pyintval::less(r.a, r.b));
    if (pyintval::strict_precedes(r.a, r.b)) CHECK(pyintval::precedes(r.a, r.b));
  }
}

// --- cancel_minus / cancel_plus ----------------------------------------------

TEST_CASE("cancel_minus: empty/unbounded conventions and exact cases") {
  // Empty a: z = empty works iff b is empty or bounded; otherwise no z exists.
  CHECK(pyintval::is_empty(pyintval::cancel_minus(kEmp, kEmp)));
  CHECK(pyintval::is_empty(pyintval::cancel_minus(kEmp, I(1, 2))));
  CHECK(pyintval::is_entire(pyintval::cancel_minus(kEmp, kEnt)));
  CHECK(pyintval::is_entire(pyintval::cancel_minus(kEmp, I(1, kInf))));
  CHECK(pyintval::is_entire(pyintval::cancel_minus(kEmp, I(-kInf, 1))));
  // Nonempty a with empty or unbounded operand: undefined -> entire.
  CHECK(pyintval::is_entire(pyintval::cancel_minus(I(1, 2), kEmp)));
  CHECK(pyintval::is_entire(pyintval::cancel_minus(kEnt, I(1, 2))));
  CHECK(pyintval::is_entire(pyintval::cancel_minus(kEnt, kEnt)));
  CHECK(pyintval::is_entire(pyintval::cancel_minus(kEnt, kEmp)));
  CHECK(pyintval::is_entire(pyintval::cancel_minus(I(1, kInf), I(1, 2))));
  CHECK(pyintval::is_entire(pyintval::cancel_minus(I(1, 2), I(1, kInf))));

  // Exact case: b + z == a exactly.
  const Interval z = pyintval::cancel_minus(I(2, 4), I(1, 2));
  CHECK(pyintval::equal(z, I(1, 2)));
  CHECK(pyintval::equal(pyintval::add(I(1, 2), z), I(2, 4)));

  // Width violation: wid(a) < wid(b) admits no z.
  CHECK(pyintval::is_entire(pyintval::cancel_minus(I(1, 2), I(0, 3))));

  // Equal widths (down to the last ulp): z collapses to a point.
  const Interval a_eq = I(1.0 + 0x1p-52, 2.0);
  const Interval b_eq = I(1.0, 2.0 - 0x1p-52);
  CHECK(pyintval::wid(a_eq) == pyintval::wid(b_eq));
  CHECK(pyintval::equal(pyintval::cancel_minus(a_eq, b_eq), I(0x1p-52, 0x1p-52)));
  // One ulp narrower on a: the exact width test must reject it.
  CHECK(pyintval::is_entire(pyintval::cancel_minus(I(1.0 + 0x1p-52, 2.0 - 0x1p-52), b_eq)));
}

TEST_CASE("cancel_minus: directed rounding on inexact endpoint differences") {
  // a.hi - b.hi = (1e300 + 1e284) - 1 is not representable: the upper bound
  // of z must round up so that b + z still covers a.
  const Interval a = I(1e300, 1e300 + 1e284);
  CHECK(a.hi > a.lo);  // 1e284 exceeds ulp(1e300)/2, so the sum rounded up
  const Interval b = I(0, 1);
  const Interval z = pyintval::cancel_minus(a, b);
  CHECK(!pyintval::is_entire(z));
  CHECK(pyintval::subset(a, pyintval::add(b, z)));
  // Tightness pins: lo is exact, hi rounds up to the next double.
  CHECK(z.lo == 1e300);
  CHECK(z.hi == a.hi);
}

TEST_CASE("cancel_minus: b + z encloses a whenever z is not entire (property, 10000 cases)") {
  std::mt19937_64 rng(0xCA9CE1B5A5A5ULL);
  int failures = 0;
  int informative = 0;
  for (int i = 0; i < 10000; ++i) {
    const Interval a = rand_interval(rng, /*allow_unbounded=*/false);
    const Interval b = rand_interval(rng, /*allow_unbounded=*/false);
    const Interval z = pyintval::cancel_minus(a, b);
    if (pyintval::is_entire(z)) continue;
    ++informative;
    const bool ok = pyintval::subset(a, pyintval::add(b, z));
    if (!ok) {
      ++failures;
      if (failures <= 3) {
        CAPTURE(i);
        CAPTURE(a);
        CAPTURE(b);
        CAPTURE(z);
        CHECK(ok);
      }
    }
  }
  CHECK(failures == 0);
  CHECK(informative >= 1000);  // the property must not be vacuous
}

TEST_CASE("cancel_minus: exact integer-scaled widths decide entire (property, 10000 cases)") {
  // Endpoints n * 2^k with |n| <= 2^40 and k in [-500, 500] are exactly
  // representable normal doubles, so interval widths are exactly (n2 - n1)
  // in units of 2^k and comparable in integer arithmetic -- an oracle for
  // the kernel's exact width test that involves no floating-point rounding.
  std::mt19937_64 rng(0x00E1D7EE2015ULL);
  int failures = 0;
  int entire_cases = 0;
  int proper_cases = 0;
  const std::int64_t span = std::int64_t{1} << 40;
  for (int i = 0; i < 10000; ++i) {
    const int k = static_cast<int>(rng() % 1001u) - 500;
    const auto ri = [&rng] {
      return static_cast<std::int64_t>(rng() % static_cast<std::uint64_t>(2 * span + 1)) - span;
    };
    std::int64_t n1 = ri(), n2 = ri(), m1 = ri(), m2 = ri();
    if (n1 > n2) std::swap(n1, n2);
    if (m1 > m2) std::swap(m1, m2);
    const Interval a = pyintval::make(std::ldexp(static_cast<double>(n1), k),
                                      std::ldexp(static_cast<double>(n2), k));
    const Interval b = pyintval::make(std::ldexp(static_cast<double>(m1), k),
                                      std::ldexp(static_cast<double>(m2), k));
    const __int128 wa = static_cast<__int128>(n2) - n1;
    const __int128 wb = static_cast<__int128>(m2) - m1;
    const Interval z = pyintval::cancel_minus(a, b);
    bool ok;
    if (wa < wb) {
      ok = pyintval::is_entire(z);
      ++entire_cases;
    } else {
      ok = !pyintval::is_entire(z) && pyintval::subset(a, pyintval::add(b, z));
      ++proper_cases;
    }
    if (!ok) {
      ++failures;
      if (failures <= 3) {
        CAPTURE(i);
        CAPTURE(a);
        CAPTURE(b);
        CAPTURE(z);
        CHECK(ok);
      }
    }
  }
  CHECK(failures == 0);
  CHECK(entire_cases >= 1000);
  CHECK(proper_cases >= 1000);
}

TEST_CASE("cancel_plus(a, b) == cancel_minus(a, -b) bit-identically (property, 5000 cases)") {
  CHECK(pyintval::equal(pyintval::cancel_plus(I(1, 2), I(-2, -1)), I(0, 0)));
  CHECK(pyintval::is_empty(pyintval::cancel_plus(kEmp, I(1, 2))));
  CHECK(pyintval::is_entire(pyintval::cancel_plus(I(1, 2), I(0, 3))));
  std::mt19937_64 rng(0xBEEFCAFE0042ULL);
  int failures = 0;
  for (int i = 0; i < 5000; ++i) {
    const Interval a = rand_interval(rng, /*allow_unbounded=*/false);
    const Interval b = rand_interval(rng, /*allow_unbounded=*/false);
    const Interval z1 = pyintval::cancel_plus(a, b);
    const Interval z2 = pyintval::cancel_minus(a, pyintval::neg(b));
    const bool ok = bit_identical(z1, z2);
    if (!ok) {
      ++failures;
      if (failures <= 3) {
        CAPTURE(i);
        CAPTURE(a);
        CAPTURE(b);
        CAPTURE(z1);
        CAPTURE(z2);
        CHECK(ok);
      }
    }
  }
  CHECK(failures == 0);
}

// --- Reverse operations ------------------------------------------------------

TEST_CASE("mul_rev_to_pair: table cases, slot discipline, and empty operands") {
  const auto [e1, e2] = pyintval::mul_rev_to_pair(kEmp, I(1, 2));
  CHECK(pyintval::is_empty(e1));
  CHECK(pyintval::is_empty(e2));
  const auto [f1, f2] = pyintval::mul_rev_to_pair(I(1, 2), kEmp);
  CHECK(pyintval::is_empty(f1));
  CHECK(pyintval::is_empty(f2));

  for (const MulRevRow& row : kMulRevRows) {
    CAPTURE(row.name);
    const auto [p1, p2] = pyintval::mul_rev_to_pair(row.b, row.c);
    CHECK(pyintval::equal(p1, row.p1));
    CHECK(pyintval::equal(p2, row.p2));
    // Slot discipline: the second slot is used only when the first is, and
    // when both are used the first holds the smaller (disjoint) piece.
    if (pyintval::is_empty(p1)) CHECK(pyintval::is_empty(p2));
    if (!pyintval::is_empty(p2)) {
      CHECK(!pyintval::is_empty(p1));
      CHECK(p1.hi <= p2.lo);
      CHECK(pyintval::disjoint(p1, p2));
    }
  }
}

TEST_CASE("mul_rev_to_pair: sampled cross-check against the defining set") {
  // Independent verification: for each table case, walk an exact dyadic grid
  // of x values. Soundness: every x with a witness v in b (v*x in c) must lie
  // in the returned union (outward rounding can only widen it). Tightness:
  // every grid x strictly inside a returned piece -- one ulp in from each
  // boundary -- must possess a witness (grid spacing 1/8 dwarfs one ulp, so
  // only exact boundary points are excused).
  for (const MulRevRow& row : kMulRevRows) {
    CAPTURE(row.name);
    const auto [p1, p2] = pyintval::mul_rev_to_pair(row.b, row.c);
    const Interval s1 = shrink_one_ulp(p1);
    const Interval s2 = shrink_one_ulp(p2);
    for (int gi = -64; gi <= 64; ++gi) {
      const double x = static_cast<double>(gi) / 8.0;
      const bool solution = mul_rev_has_solution(x, row.b, row.c);
      const bool in_union = pyintval::is_member(x, p1) || pyintval::is_member(x, p2);
      const bool in_shrunk = pyintval::is_member(x, s1) || pyintval::is_member(x, s2);
      CAPTURE(x);
      CHECK((!solution || in_union));   // soundness
      CHECK((!in_shrunk || solution));  // tightness up to one ulp
    }
  }
}

TEST_CASE("mul_rev: hull of the pair; the x overload intersects before hulling") {
  CHECK(pyintval::equal(pyintval::mul_rev(I(2, 4), I(6, 8)), I(1.5, 4)));
  // Two half-line pieces: the bare hull spans the gap and becomes entire.
  CHECK(pyintval::is_entire(pyintval::mul_rev(I(-1, 1), I(1, 2))));
  // Intersecting with x first keeps only the reachable piece: tighter than
  // hulling first ([1,3] vs [0.5,3]).
  CHECK(pyintval::equal(pyintval::mul_rev(I(-1, 1), I(1, 2), I(0.5, 3)), I(1, 3)));
  CHECK(pyintval::equal(pyintval::mul_rev(I(-1, 1), I(1, 2), I(-3, -0.5)), I(-3, -1)));
  CHECK(pyintval::is_empty(pyintval::mul_rev(I(2, 4), I(6, 8), I(0, 1))));
  CHECK(pyintval::is_empty(pyintval::mul_rev(I(0, 0), I(1, 2))));
  CHECK(pyintval::is_empty(pyintval::mul_rev(kEmp, I(1, 2))));
  CHECK(pyintval::is_entire(pyintval::mul_rev(I(0, 0), I(-1, 1))));
  // Consistency: the bare overload equals the hull of the pair on all rows.
  for (const MulRevRow& row : kMulRevRows) {
    CAPTURE(row.name);
    const auto [p1, p2] = pyintval::mul_rev_to_pair(row.b, row.c);
    CHECK(pyintval::equal(pyintval::mul_rev(row.b, row.c), pyintval::convex_hull(p1, p2)));
  }
}

TEST_CASE("sqr_rev") {
  CHECK(pyintval::equal(pyintval::sqr_rev(I(4, 9)), I(-3, 3)));
  CHECK(pyintval::equal(pyintval::sqr_rev(I(4, 9), I(0, 5)), I(2, 3)));
  CHECK(pyintval::equal(pyintval::sqr_rev(I(4, 9), I(-5, 0)), I(-3, -2)));
  CHECK(pyintval::equal(pyintval::sqr_rev(I(4, 9), I(2.5, 10)), I(2.5, 3)));
  CHECK(pyintval::equal(pyintval::sqr_rev(I(4, 9), I(-2.5, 2.5)), I(-2.5, 2.5)));
  CHECK(pyintval::is_empty(pyintval::sqr_rev(I(-2, -1))));  // u^2 is never negative
  CHECK(pyintval::equal(pyintval::sqr_rev(I(-1, 4)), I(-2, 2)));
  CHECK(pyintval::equal(pyintval::sqr_rev(I(0, 4)), I(-2, 2)));
  CHECK(pyintval::is_empty(pyintval::sqr_rev(kEmp)));
  CHECK(pyintval::is_empty(pyintval::sqr_rev(I(4, 9), kEmp)));
  CHECK(pyintval::is_entire(pyintval::sqr_rev(kEnt)));

  // sqrt(2) is irrational: the enclosure must widen one ulp outward, giving
  // the symmetric hull [-sqrt_ru(2), sqrt_ru(2)].
  const Interval s2 = pyintval::sqr_rev(pyintval::point(2.0));
  const double r2u = pd::sqrt_ru(2.0);
  CHECK(pd::sqrt_rd(2.0) < r2u);  // genuinely inexact, so the pair straddles
  CHECK(s2.hi == r2u);
  CHECK(s2.lo == -r2u);
  CHECK(pyintval::is_member(std::sqrt(2.0), s2));
  CHECK(pyintval::subset(pyintval::point(2.0), pyintval::sqr(s2)));
}

TEST_CASE("abs_rev") {
  CHECK(pyintval::equal(pyintval::abs_rev(I(1, 2)), I(-2, 2)));
  CHECK(pyintval::equal(pyintval::abs_rev(I(1, 2), I(0, kInf)), I(1, 2)));
  CHECK(pyintval::equal(pyintval::abs_rev(I(1, 2), I(-kInf, 0)), I(-2, -1)));
  CHECK(pyintval::is_empty(pyintval::abs_rev(I(-3, -1))));  // |u| is never negative
  CHECK(pyintval::is_empty(pyintval::abs_rev(I(-3, -1), kEnt)));
  CHECK(pyintval::is_empty(pyintval::abs_rev(I(-3, -1), I(-10, 10))));
  CHECK(pyintval::equal(pyintval::abs_rev(I(-1, 2)), I(-2, 2)));
  CHECK(pyintval::equal(pyintval::abs_rev(I(-1, 2), I(-1, 0)), I(-1, 0)));
  CHECK(pyintval::equal(pyintval::abs_rev(I(0, 3), I(-1, 5)), I(-1, 3)));
  CHECK(pyintval::is_empty(pyintval::abs_rev(kEmp)));
  CHECK(pyintval::is_empty(pyintval::abs_rev(I(1, 2), kEmp)));
  CHECK(pyintval::is_entire(pyintval::abs_rev(kEnt)));
}

// --- Diagnostics -------------------------------------------------------------

TEST_CASE("operator<< prints [empty] for the empty set, brackets otherwise") {
  std::ostringstream oss_empty;
  oss_empty << kEmp;
  CHECK(oss_empty.str() == "[empty]");
  std::ostringstream oss;
  oss << I(1, 2);
  const std::string s = oss.str();
  CHECK(s == "[1, 2]");
  std::ostringstream oss_ent;
  oss_ent << kEnt;
  const std::string t = oss_ent.str();
  CHECK(t.front() == '[');
  CHECK(t.back() == ']');
  CHECK(t.find(',') != std::string::npos);
}
