// Regression tests for enclosure/soundness bugs found by the adversarial rigor
// audit. Each pins a concrete counterexample so the bug cannot silently return.

#include <cmath>
#include <limits>

#include "doctest/doctest.h"
#include "pyintval/elementary.hpp"
#include "pyintval/interval.hpp"
#include "pyintval/text.hpp"

using namespace pyintval;

namespace {
constexpr double kMax = std::numeric_limits<double>::max();
constexpr double kInfv = std::numeric_limits<double>::infinity();
}  // namespace

// Bug: Knuth TwoSum's `bb = s - a` overflowed to +-inf near maxDouble, yielding
// a NaN error that disabled the directed-rounding adjustment, so add/sub
// returned a degenerate interval NOT containing the true sum.
TEST_CASE("regression: two_sum does not overflow internally near maxDouble") {
  const auto p = detail::two_sum(-0x1.909ef1ec45b4cp+1020, kMax);
  CHECK(p.err == p.err);          // not NaN
  CHECK(std::isfinite(p.err));    // exact finite error
  CHECK(p.val + p.err == p.val);  // val already absorbs the tiny error here
}

TEST_CASE("regression: add near overflow is a rigorous outward bound") {
  const double a = -0x1.d3c6ee130d306p+1021;
  const Interval r = add(make(a, a), make(kMax, kMax));
  CHECK(r.lo < r.hi);                 // nondegenerate: the true sum is enclosed
  CHECK(detail::succ(r.lo) == r.hi);  // tight: exactly one ulp wide
  // The exact sum a + kMax lies strictly inside (verified with rationals in the
  // Python audit); here we pin the exact endpoints the fix must produce.
  CHECK(r.lo == 0x1.8b0e447b3cb3dp+1023);
  CHECK(r.hi == 0x1.8b0e447b3cb3ep+1023);
  // sub() shares the code path.
  const Interval s = sub(make(-a, -a), make(kMax, kMax));
  CHECK(s.lo <= s.hi);
}

// Bug: cancel_minus's exact width test used TwoSum results that overflowed,
// corrupting exact_sum_sign and returning a bounded interval instead of Entire.
TEST_CASE("regression: cancel_minus returns Entire when a sum overflows") {
  const Interval r =
      cancel_minus(make(-kMax, -kMax), make(-0x1.15346c2bf6e36p+970, -0x1.4e5e32df92bc2p-33));
  CHECK(is_entire(r));
}

// Bug: sin/cos/tan of a large-magnitude POINT collapsed to [-1,1] / entire
// because the periodic-extremum test fired on the loose argument reduction,
// discarding the exact CORE-MATH value.
TEST_CASE("regression: sin/cos of a large-magnitude point stay tight") {
  // A point interval yields [pred(cr), succ(cr)] -- two ulps wide, NOT [-1, 1].
  auto narrow = [](const Interval& x) { return detail::succ(detail::succ(x.lo)) == x.hi; };
  const Interval s = sin(point(0x1p53));
  CHECK(s.lo < s.hi);
  CHECK(narrow(s));
  CHECK(s.lo >= -1.0);
  CHECK(s.hi <= 1.0);
  CHECK(is_member(cr_sin(0x1p53), s));  // encloses the correctly rounded value

  const Interval c = cos(point(0x1p60));
  CHECK(narrow(c));
  CHECK(is_member(cr_cos(0x1p60), c));

  const Interval t = tan(point(0x1p53));  // not near an asymptote -> finite, tight
  CHECK(!is_entire(t));
  CHECK(narrow(t));
  CHECK(is_member(cr_tan(0x1p53), t));
}

// Bug: a literal whose exact value fell in (DBL_MAX, DBL_MAX + 1/2 ulp) sent
// text_to_interval into an infinite loop via an out-of-range double->uint64 cast
// on +inf. (If this regressed, the test would hang rather than fail.)
TEST_CASE("regression: near-overflow literal terminates and brackets correctly") {
  Interval x;
  REQUIRE(text_to_interval("1.79769313486231577e308", x));
  CHECK(x.lo == kMax);
  CHECK(x.hi == kInfv);
  REQUIRE(text_to_interval("-1.79769313486231577e308", x));
  CHECK(x.lo == -kInfv);
  CHECK(x.hi == -kMax);
  REQUIRE(text_to_interval("0x1.fffffffffffff8p1023", x));  // hex just above DBL_MAX
  CHECK(x.lo == kMax);
  CHECK(x.hi == kInfv);
}
