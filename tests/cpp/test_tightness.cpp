// Tightness (overestimation-control) tests for *composed* interval expressions,
// exercised directly on the header-only C++ kernel.
//
// Enclosure tests prove a result CONTAINS the truth; these prove a result is
// not needlessly WIDE. For each fixed evaluation form the tightest correctly-
// rounded interval is unique -- and for the polynomial cases below it has exact
// integer endpoints, so a literal endpoint check IS a tightness check (the exact
// same minimal intervals the Python suite pins bit-for-bit against a rational
// oracle in tests/python/test_tightness.py). The cases also document the two
// irreducible overestimation sources -- the dependency problem and the wrapping
// effect -- and confirm the kernel still rigorously encloses the truth.
//
// These are written as deeply nested free-function calls on purpose: fused
// interval expressions are exactly the pattern that must survive aggressive
// optimization on the header-only path (see PYINTVAL_NOINLINE_GCC), so composed
// expressions belong in the C++ suite and not only behind the Python bindings.

#include "doctest/doctest.h"
#include "pyintval/interval.hpp"

namespace pv = pyintval;

namespace {
pv::Interval I(double lo, double hi) { return pv::make(lo, hi); }
pv::Interval P(double v) { return pv::point(v); }
bool enc(const pv::Interval& r, double z) { return !pv::is_empty(r) && r.lo <= z && z <= r.hi; }
bool eq(const pv::Interval& r, double lo, double hi) { return r.lo == lo && r.hi == hi; }
}  // namespace

TEST_CASE("dependency: x - x is not zero, but is exactly as tight as the form allows") {
  using namespace pv;
  const Interval x = I(0.0, 1.0);
  CHECK(eq(sub(x, x), -1.0, 1.0));           // [lo-hi, hi-lo], the dependency width
  CHECK(enc(sub(x, x), 0.0));                // the true value {0} is still enclosed
  CHECK(eq(sub(P(2.5), P(2.5)), 0.0, 0.0));  // point input: no dependency, collapses
}

TEST_CASE("dependency: logistic x*(1-x) vs x - x*x on [0,1], true range [0, 1/4]") {
  using namespace pv;
  const Interval x = I(0.0, 1.0);
  const Interval f1 = mul(x, sub(P(1.0), x));  // 4x too wide up top
  const Interval f2 = sub(x, mul(x, x));       // 8x too wide, wrong sign below
  CHECK(eq(f1, 0.0, 1.0));
  CHECK(eq(f2, -1.0, 1.0));
  CHECK(subset(f1, f2));  // the x*(1-x) form is strictly tighter
  CHECK(enc(f1, 0.25));   // both enclose the true max 1/4 ...
  CHECK(enc(f2, 0.25));
  CHECK((f1.lo <= 0.0 && f2.lo <= 0.0));  // ... and the true min 0
}

TEST_CASE("evaluation order: quartic x^4-4x^3+4x^2 on [0,3], exact range [0, 9]") {
  using namespace pv;
  const Interval x = I(0.0, 3.0);
  const Interval naive = add(sub(pown(x, 4), mul(P(4.0), pown(x, 3))), mul(P(4.0), pown(x, 2)));
  const Interval horner = mul(mul(add(mul(sub(x, P(4.0)), x), P(4.0)), x), x);
  const Interval factored = mul(sqr(x), sqr(sub(x, P(2.0))));  // x^2 (x-2)^2
  CHECK(eq(naive, -108.0, 117.0));
  CHECK(eq(horner, -72.0, 36.0));
  CHECK(eq(factored, 0.0, 36.0));
  CHECK(subset(factored, horner));  // tighter rearrangements nest in looser ones
  CHECK(subset(horner, naive));
  for (const Interval& f : {naive, horner, factored}) {
    CHECK(f.lo <= 0.0);  // every form encloses the exact range [0, 9]
    CHECK(9.0 <= f.hi);
  }
}

TEST_CASE("evaluation order: cubic x^3-x-1 on [1,2], Horner attains exact range [-1, 5]") {
  using namespace pv;
  const Interval x = I(1.0, 2.0);
  const Interval naive = sub(sub(pown(x, 3), x), P(1.0));
  const Interval horner = sub(mul(sub(mul(x, x), P(1.0)), x), P(1.0));
  CHECK(eq(naive, -2.0, 6.0));
  CHECK(eq(horner, -1.0, 5.0));  // == the exact range
  CHECK(subset(horner, naive));
  CHECK(enc(naive, 0.0));  // both enclose the real root (plastic number)
  CHECK(enc(horner, 0.0));
}

TEST_CASE("Rump's example: enclosure survives catastrophic cancellation") {
  using namespace pv;
  const Interval A = P(77617.0), B = P(33096.0);
  const Interval b2 = pown(B, 2), b4 = pown(B, 4), b6 = pown(B, 6), b8 = pown(B, 8);
  const Interval a2 = pown(A, 2);
  const Interval f =
      add(add(add(mul(P(333.75), b6),
                  mul(a2, sub(sub(sub(mul(mul(P(11.0), a2), b2), b6), mul(P(121.0), b4)), P(2.0)))),
              mul(P(5.5), b8)),
          div(A, mul(P(2.0), B)));
  const double truth = -0.8273960599468214;  // -54767/66192
  CHECK(enc(f, truth));                      // rigorous enclosure of the true value ...
  CHECK(wid(f) > 1.0);                       // ... but honestly wide: cannot even certify the sign
}

TEST_CASE("wrapping effect: eight 45-degree rotations of [-1,1]^2 grow ~16x per axis") {
  using namespace pv;
  const Interval s = recip(sqrt(P(2.0)));  // 1/sqrt(2)
  Interval x = I(-1.0, 1.0), y = I(-1.0, 1.0);
  for (int k = 0; k < 8; ++k) {  // 8 * 45 = 360 degrees: the true image is the input
    const Interval nx = mul(sub(x, y), s);
    const Interval ny = mul(add(x, y), s);
    x = nx;
    y = ny;
  }
  CHECK(subset(I(-1.0, 1.0), x));  // truth stays enclosed after a full turn
  CHECK(subset(I(-1.0, 1.0), y));
  CHECK(wid(x) > 8.0);   // wrapping: far wider than the true width 2 ...
  CHECK(wid(x) < 33.0);  // ... but bounded by (sqrt2)^8 = 16x growth
  CHECK(wid(y) < 33.0);
}
