// Tests for the IEEE 1788 decoration system: that each operation reports the
// correct local decoration, and that composition propagates the weakest one.

#include "doctest/doctest.h"
#include "pyintval/decoration.hpp"

using namespace pyintval;

namespace {
Decoration dec_of(const DecoratedInterval& d) { return d.dec; }
DecoratedInterval D(double lo, double hi) { return decorate(make(lo, hi)); }
}  // namespace

TEST_CASE("construction assigns the tightest decoration") {
  CHECK(D(1.0, 2.0).dec == Decoration::com);                       // bounded, nonempty
  CHECK(decorate(make(0.0, detail::kInf)).dec == Decoration::dac);  // unbounded
  CHECK(decorate(entire()).dec == Decoration::dac);
  CHECK(decorate(empty()).dec == Decoration::trv);
  CHECK(nai().dec == Decoration::ill);
}

TEST_CASE("continuous operations preserve com/dac") {
  CHECK(add(D(1, 2), D(3, 4)).dec == Decoration::com);
  CHECK(mul(D(1, 2), D(3, 4)).dec == Decoration::com);
  CHECK(exp(D(0, 1)).dec == Decoration::com);
  CHECK(sin(D(0, 1)).dec == Decoration::com);
  CHECK(sqr(D(-2, 3)).dec == Decoration::com);
  // An unbounded (but still defined & continuous) input drops com to dac.
  CHECK(add(D(1, 2), decorate(make(0.0, detail::kInf))).dec == Decoration::dac);
  CHECK(exp(decorate(make(-detail::kInf, 0.0))).dec == Decoration::dac);
}

TEST_CASE("square root domain") {
  CHECK(sqrt(D(1, 4)).dec == Decoration::com);   // strictly inside domain
  CHECK(sqrt(D(0, 4)).dec == Decoration::com);   // 0 is in the domain, continuous
  CHECK(sqrt(D(-1, 4)).dec == Decoration::trv);  // domain violated
  CHECK(sqrt(D(-4, -1)).dec == Decoration::trv);  // fully outside -> empty, trv
}

TEST_CASE("logarithms require strictly positive input") {
  CHECK(log(D(1, 2)).dec == Decoration::com);
  CHECK(log(D(0, 2)).dec == Decoration::trv);   // log undefined at 0
  CHECK(log(D(-1, 2)).dec == Decoration::trv);
  CHECK(log1p(D(-0.5, 2.0)).dec == Decoration::com);
  CHECK(log1p(D(-1.0, 2.0)).dec == Decoration::trv);
}

TEST_CASE("division and reciprocal by a zero-containing divisor are trv") {
  CHECK(div(D(1, 2), D(3, 4)).dec == Decoration::com);
  CHECK(div(D(1, 2), D(-1, 1)).dec == Decoration::trv);
  CHECK(div(D(1, 2), D(0, 1)).dec == Decoration::trv);
  CHECK(recip(D(2, 4)).dec == Decoration::com);
  CHECK(recip(D(-1, 1)).dec == Decoration::trv);
}

TEST_CASE("inverse trig and hyperbolic domains") {
  CHECK(asin(D(-0.5, 0.5)).dec == Decoration::com);
  CHECK(asin(D(-1.0, 1.0)).dec == Decoration::com);   // closed domain endpoints ok
  CHECK(asin(D(-2.0, 0.5)).dec == Decoration::trv);
  CHECK(acos(D(0.0, 1.0)).dec == Decoration::com);
  CHECK(acos(D(-1.5, 0.0)).dec == Decoration::trv);
  CHECK(acosh(D(1.0, 3.0)).dec == Decoration::com);
  CHECK(acosh(D(0.5, 3.0)).dec == Decoration::trv);
  CHECK(atanh(D(-0.5, 0.5)).dec == Decoration::com);
  CHECK(atanh(D(-1.0, 0.5)).dec == Decoration::trv);   // open domain endpoint
}

TEST_CASE("tan asymptotes drop to trv") {
  CHECK(tan(D(0.1, 0.2)).dec == Decoration::com);
  CHECK(tan(D(1.5, 1.7)).dec == Decoration::trv);   // straddles pi/2
  CHECK(tan(D(0.0, 4.0)).dec == Decoration::trv);   // spans a full period
}

TEST_CASE("step functions are def unless locally constant") {
  CHECK(floor(D(0.2, 0.8)).dec == Decoration::com);   // constant on the cell
  CHECK(floor(D(0.5, 1.5)).dec == Decoration::def);   // crosses an integer
  CHECK(ceil(D(0.2, 0.8)).dec == Decoration::com);
  CHECK(trunc(D(-0.5, 0.5)).dec == Decoration::com);  // trunc == 0 on (-1,1): continuous
  CHECK(trunc(D(0.5, 1.5)).dec == Decoration::def);   // crosses the jump at 1
  CHECK(sign(D(1.0, 2.0)).dec == Decoration::com);
  CHECK(sign(D(-1.0, 1.0)).dec == Decoration::def);
}

TEST_CASE("atan2 origin, branch cut, and clean regions") {
  CHECK(atan2(D(1, 2), D(1, 2)).dec == Decoration::com);       // first quadrant
  CHECK(atan2(D(-1, 1), D(-2, -1)).dec == Decoration::def);    // crosses branch cut
  CHECK(atan2(D(-1, 1), D(-1, 1)).dec == Decoration::trv);     // contains origin
  CHECK(hypot(D(3, 4), D(-1, 2)).dec == Decoration::com);      // continuous everywhere
}

TEST_CASE("pow domain") {
  CHECK(pow(D(2, 3), decorate(point(0.5))).dec == Decoration::com);   // positive base
  CHECK(pow(D(-1, 3), decorate(point(0.5))).dec == Decoration::trv);  // base can be < 0
  CHECK(pow(D(-2, 3), decorate(point(2.0))).dec == Decoration::com);  // integer exponent
  CHECK(pow(D(-2, 3), decorate(point(-1.0))).dec == Decoration::trv);  // n<0 with 0 in base
}

TEST_CASE("composition propagates the weakest decoration") {
  // A single domain violation anywhere poisons the whole certificate.
  CHECK(exp(sqrt(D(-1, 4))).dec == Decoration::trv);
  CHECK(add(log(D(0, 2)), D(1, 2)).dec == Decoration::trv);
  // A clean composition certifies defined-and-continuous.
  DecoratedInterval f = div(sin(D(0, 1)), add(sqr(D(0, 1)), decorate(point(1.0))));
  CHECK(static_cast<int>(f.dec) >= static_cast<int>(Decoration::dac));
  // ill poisons regardless of the other operand.
  CHECK(add(nai(), D(1, 2)).dec == Decoration::ill);
  CHECK(mul(D(1, 2), nai()).dec == Decoration::ill);
}
