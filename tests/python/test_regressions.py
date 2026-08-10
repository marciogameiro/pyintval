"""Regression tests for bugs found by the adversarial rigor audit.

Each test pins a concrete counterexample so the defect cannot silently return.
Kernel-level regressions (TwoSum overflow, text hang, sin/cos tightness) live in
the C++ suite; these cover the Python binding layer.
"""

import math
from fractions import Fraction

import pytest

import pyintval as iv

I = iv.Interval
DI = iv.DecoratedInterval


# Bug: Python ints >= 2**53 were promoted to a round-to-nearest double point,
# so operations returned a degenerate interval NOT containing the true value.
class TestBigIntPromotion:
    @pytest.mark.parametrize("n", [2**53 + 1, 2**60 + 1, 10**20, -(3**42), 2**70 - 3])
    def test_scalar_promotion_encloses_exactly(self, n):
        r = I(1) * n  # promote n on the right
        assert Fraction(r.lo) <= n <= Fraction(r.hi)
        r = n + I(0)  # reflected
        assert Fraction(r.lo) <= n <= Fraction(r.hi)

    @pytest.mark.parametrize("n", [2**53 + 1, 10**30, -(2**64) + 5])
    def test_point_construction_encloses(self, n):
        p = I(n)
        assert Fraction(p.lo) <= n <= Fraction(p.hi)
        assert n in p

    def test_explicit_bounds_round_outward(self):
        x = I(2**53 + 1, 2**53 + 3)
        assert Fraction(x.lo) <= 2**53 + 1
        assert 2**53 + 3 <= Fraction(x.hi)

    def test_small_ints_stay_exact_points(self):
        for n in range(-100, 101):
            assert I(n).lo == I(n).hi == float(n)

    def test_huge_int_overflows_to_unbounded(self):
        x = I(10**400)
        assert x.hi == math.inf
        assert Fraction(x.lo) <= 10**400
        x = I(-(10**400))
        assert x.lo == -math.inf

    def test_decorated_big_int_promotion(self):
        r = DI(1) * (2**60 + 1)
        assert Fraction(r.lo) <= 2**60 + 1 <= Fraction(r.hi)
        assert r.decoration == "com"


# Bug: the 'com' decoration was assigned on unbounded results (should be 'dac').
class TestDecorationComOverclaim:
    def test_overflow_results_are_dac_not_com(self):
        d = iv.exp(DI(709.8, 710.0))
        assert d.hi == math.inf
        assert d.decoration == "dac"

    def test_cosh_reaching_infinity_is_dac(self):
        d = iv.cosh(DI(710.0, 711.0))
        assert d.hi == math.inf
        assert d.decoration == "dac"

    def test_bounded_continuous_stays_com(self):
        assert iv.exp(DI(0, 1)).decoration == "com"
        assert (DI(1, 2) * DI(3, 4)).decoration == "com"

    def test_recip_near_zero_unbounded_is_dac_or_trv(self):
        # 1/[0.5, inf) style results are unbounded -> at most dac.
        d = iv.recip(DI(1e-320, 1.0))  # tiny lower bound -> huge upper
        assert d.decoration in ("dac", "trv")


# Bug: DecoratedInterval.from_parts minted false certificates.
class TestFromPartsValidation:
    def test_com_requires_common_interval(self):
        with pytest.raises(ValueError):
            DI.from_parts(iv.entire(), "com")
        with pytest.raises(ValueError):
            DI.from_parts(I(0, math.inf), "com")

    def test_defined_requires_nonempty(self):
        with pytest.raises(ValueError):
            DI.from_parts(I("[empty]"), "dac")
        with pytest.raises(ValueError):
            DI.from_parts(I("[empty]"), "def")

    def test_valid_combinations_accepted(self):
        assert DI.from_parts(I(1, 2), "com").decoration == "com"
        assert DI.from_parts(iv.entire(), "dac").decoration == "dac"
        assert DI.from_parts(I("[empty]"), "trv").decoration == "trv"
        assert DI.from_parts(I(1, 2), "trv").decoration == "trv"  # weakening is allowed


# Bug: DecoratedInterval repr used shortest decimals (widen on reparse); it must
# be exact so eval(repr(x)) reconstructs bit-for-bit.
def test_decorated_repr_is_exact():
    for d in [DI("0.1"), DI(1, 2), DI(0, math.inf), iv.sqrt(DI(-1, 4)), DI.nai()]:
        ns = {"DecoratedInterval": DI, "Interval": iv.Interval}
        assert eval(repr(d), ns) == d


# Cross-check: the near-overflow parse (that used to hang) is reachable and safe
# from the Python API too.
def test_near_overflow_literal_does_not_hang():
    x = I("1.79769313486231577e308")
    assert x.lo == 1.7976931348623157e308
    assert x.hi == math.inf


# Bug (F2, found by the ITF1788 conformance suite): pow(x, y) with a base
# touching 0 returned [empty] via exp(y * log(0)) instead of the IEEE 1788
# continuity value. pow is defined for base x >= 0; the x = 0 boundary
# contributes limit values, except that for a base of exactly {0} the point 0^0
# is undefined and 0^{y<=0} is excluded.
class TestPowBaseZero:
    def test_pure_zero_base(self):
        assert iv.pow(I(0, 0), I(0, 1)) == I(0, 0)  # 0^{y>0} = 0
        assert iv.pow(I(0, 0), I(0, 0)).is_empty  # 0^0 undefined
        assert iv.pow(I(0, 0), I(-1, -1)).is_empty  # 0^-1 undefined
        assert iv.pow(I(0, 0), I(0, 2.5)) == I(0, 0)

    def test_base_spanning_zero_encloses(self):
        # Was [empty] before the fix; must now enclose the true range.
        r = iv.pow(I(0, 0.5), I(2.5, 2.5))
        assert 0.0 in r and 0.5**2.5 in r
        assert iv.pow(I(0, 0.5), I(0, 0)) == I(1, 1)  # x^0 = 1
        assert iv.pow(I(0, 0.5), I(-1, -1)) == I(2, math.inf)  # 0^-1 -> +inf

    def test_base_below_one_negative_infinite_exponent(self):
        # x < 1 with y -> -inf gives +inf; the corner must not vanish.
        assert iv.pow(I(0.1, 0.5), I(-math.inf, 0.1)).hi == math.inf
        assert iv.pow(I(0.1, 0.5), I.entire()) == I(0, math.inf)
