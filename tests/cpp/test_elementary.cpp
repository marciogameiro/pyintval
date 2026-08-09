// Structural and known-value tests for the interval elementary functions.
// The authoritative rigorous enclosure check (random inputs vs a high-precision
// mpmath oracle) lives in tests/python/test_elementary.py; this file pins the
// mathematical structure that must hold on every platform without an oracle:
// known values, domain handling, periodic extrema, unbounded limits, empty
// propagation, and tightness.

#include <cmath>
#include <limits>

#include "doctest/doctest.h"
#include "pyintval/elementary.hpp"

using namespace pyintval;

namespace {
constexpr double kInfv = std::numeric_limits<double>::infinity();
bool holds(double v, const Interval& x) { return is_member(v, x); }
}  // namespace

TEST_CASE("known values are enclosed") {
  CHECK(holds(1.0, exp(point(0.0))));
  CHECK(exp(point(1.0)).lo <= 2.718281828459045);
  CHECK(2.718281828459045 <= exp(point(1.0)).hi);
  CHECK(holds(0.0, log(point(1.0))));
  CHECK(holds(1.0, exp2(point(0.0))));
  CHECK(holds(8.0, exp2(point(3.0))));
  CHECK(holds(3.0, log2(point(8.0))));
  CHECK(holds(2.0, log10(point(100.0))));
  CHECK(holds(0.0, expm1(point(0.0))));
  CHECK(holds(0.0, log1p(point(0.0))));
  CHECK(holds(2.0, cbrt(point(8.0))));
  CHECK(holds(-3.0, cbrt(point(-27.0))));
  CHECK(holds(0.0, sin(point(0.0))));
  CHECK(holds(1.0, cos(point(0.0))));
  CHECK(holds(0.0, tan(point(0.0))));
  CHECK(holds(0.0, sinh(point(0.0))));
  CHECK(holds(1.0, cosh(point(0.0))));
  CHECK(holds(0.0, tanh(point(0.0))));
  CHECK(holds(0.0, asin(point(0.0))));
  CHECK(holds(0.0, atan(point(0.0))));
  CHECK(holds(0.0, asinh(point(0.0))));
  CHECK(holds(0.0, acosh(point(1.0))));
  CHECK(holds(0.0, atanh(point(0.0))));
  CHECK(holds(0.0, erf(point(0.0))));
  CHECK(holds(1.0, erfc(point(0.0))));
  CHECK(holds(5.0, hypot(point(3.0), point(4.0))));
}

TEST_CASE("bounded ranges are respected") {
  CHECK(subset(tanh(entire()), make(-1.0, 1.0)));
  CHECK(subset(erf(entire()), make(-1.0, 1.0)));
  CHECK(subset(atan(entire()), make(-2.0, 2.0)));  // within (-pi/2, pi/2)
  CHECK(subset(sin(make(-100.0, 100.0)), make(-1.0, 1.0)));
  CHECK(subset(cos(make(-100.0, 100.0)), make(-1.0, 1.0)));
  CHECK(cosh(make(-2.0, 3.0)).lo >= 1.0);
  CHECK(erfc(entire()).lo >= 0.0);
  CHECK(erfc(entire()).hi <= 2.0);
}

TEST_CASE("periodic extrema are captured") {
  // Interval spanning a full period -> full range.
  CHECK(equal(sin(make(0.0, 7.0)), make(-1.0, 1.0)));
  CHECK(equal(cos(make(0.0, 7.0)), make(-1.0, 1.0)));
  // Interval straddling pi/2 reaches the maximum 1.
  CHECK(sin(make(1.5, 1.7)).hi == 1.0);
  // Interval straddling pi reaches the minimum -1 for cos.
  CHECK(cos(make(3.0, 3.3)).lo == -1.0);
  // A narrow interval away from any extremum does NOT pin to +-1.
  CHECK(sin(make(0.1, 0.2)).hi < 1.0);
  CHECK(sin(make(0.1, 0.2)).lo > -1.0);
}

TEST_CASE("tan asymptotes give the entire line") {
  CHECK(is_entire(tan(make(1.5, 1.7))));  // straddles pi/2
  CHECK(is_entire(tan(make(0.0, 4.0))));  // spans more than a period
  CHECK(is_entire(tan(entire())));
  // A branch-interior interval stays bounded and increasing.
  Interval t = tan(make(0.1, 0.2));
  CHECK(!is_entire(t));
  CHECK(t.lo <= t.hi);
}

TEST_CASE("domains restrict toward empty") {
  CHECK(is_empty(log(make(-2.0, -1.0))));
  CHECK(is_empty(log2(make(-2.0, -1.0))));
  CHECK(is_empty(log1p(make(-3.0, -2.0))));
  CHECK(is_empty(asin(make(2.0, 3.0))));
  CHECK(is_empty(acos(make(-3.0, -2.0))));
  CHECK(is_empty(acosh(make(0.0, 0.5))));
  CHECK(is_empty(atanh(make(-3.0, -1.5))));
  // Partial-domain input is clipped, not rejected.
  CHECK(log(make(-1.0, 4.0)).lo == -kInfv);  // approaches 0+ -> -inf
  CHECK(equal(asin(make(-2.0, 0.0)), asin(make(-1.0, 0.0))));
  CHECK(acosh(make(0.5, 2.0)).lo >= 0.0);
}

TEST_CASE("unbounded inputs and limits") {
  CHECK(exp(make(-kInfv, 0.0)).lo == 0.0);
  CHECK(exp(make(-kInfv, 0.0)).hi <= 1.0 + 1e-15);
  CHECK(exp(make(0.0, kInfv)).hi == kInfv);
  CHECK(expm1(make(-kInfv, 0.0)).lo == -1.0);
  CHECK(log(make(0.0, kInfv)).lo == -kInfv);
  CHECK(atanh(make(0.5, 1.0)).hi == kInfv);
  CHECK(sinh(make(-kInfv, kInfv)).lo == -kInfv);
  CHECK(cosh(make(-kInfv, kInfv)).hi == kInfv);
}

TEST_CASE("empty propagates through every function") {
  const Interval e = empty();
  CHECK(is_empty(exp(e)));
  CHECK(is_empty(log(e)));
  CHECK(is_empty(sin(e)));
  CHECK(is_empty(cos(e)));
  CHECK(is_empty(tan(e)));
  CHECK(is_empty(sqrt(e)));
  CHECK(is_empty(atan2(e, point(1.0))));
  CHECK(is_empty(hypot(e, point(1.0))));
  CHECK(is_empty(pow(e, point(2.0))));
}

TEST_CASE("results are tight (about two ulps) at points") {
  auto narrow = [](const Interval& x) {
    return wid(x) <= 4.0 * std::ldexp(std::max(mag(x), 1.0), -52);
  };
  CHECK(narrow(exp(point(1.0))));
  CHECK(narrow(log(point(3.0))));
  CHECK(narrow(sin(point(1.0))));
  CHECK(narrow(cos(point(1.0))));
  CHECK(narrow(atan(point(0.5))));
  CHECK(narrow(cbrt(point(5.0))));
}

TEST_CASE("pow: integer exponents match pown; real powers via nonnegative base") {
  CHECK(equal(pow(make(-2.0, 3.0), point(2.0)), pown(make(-2.0, 3.0), 2)));
  CHECK(equal(pow(make(-2.0, 3.0), point(3.0)), pown(make(-2.0, 3.0), 3)));
  CHECK(holds(8.0, pow(point(2.0), point(3.0))));
  CHECK(holds(std::sqrt(2.0), pow(point(2.0), point(0.5))));
  // Negative base with non-integer exponent is undefined -> empty.
  CHECK(is_empty(pow(make(-3.0, -1.0), point(0.5))));
}

TEST_CASE("atan2 quadrants and branch cut") {
  CHECK(holds(std::atan2(1.0, 1.0), atan2(point(1.0), point(1.0))));    // ~ pi/4
  CHECK(holds(std::atan2(1.0, -1.0), atan2(point(1.0), point(-1.0))));  // ~ 3pi/4
  CHECK(holds(std::atan2(-1.0, -1.0), atan2(point(-1.0), point(-1.0))));
  // Box straddling the negative x-axis: conservative full-angle enclosure.
  Interval wrap = atan2(make(-1.0, 1.0), make(-2.0, -1.0));
  CHECK(wrap.lo <= -3.14);
  CHECK(wrap.hi >= 3.14);
}
