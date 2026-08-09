"""Arithmetic operators, scalar promotion, set ops, and predicates."""

import math

import pytest

import pyintval as iv

I = iv.Interval


def test_addition_and_subtraction():
    assert (I(1, 2) + I(3, 4)).endpoints == (4.0, 6.0)
    assert (I(1, 2) - I(3, 4)).endpoints == (-3.0, -1.0)


def test_multiplication_sign_classes():
    assert (I(2, 3) * I(4, 5)).endpoints == (8.0, 15.0)
    assert (I(-2, 3) * I(-5, 7)).endpoints == (-15.0, 21.0)
    assert (I(0, 1) * I(2, math.inf)).endpoints == (0.0, math.inf)


def test_division_set_based():
    assert (I(1, 2) / I(4, 8)).endpoints == (0.125, 0.5)
    assert (I(1, 2) / I(-1, 1)).is_entire
    assert (I(1, 2) / I(0, 1)).endpoints == (1.0, math.inf)
    assert (I(1, 2) / I(0, 0)).is_empty  # empty domain


def test_scalar_promotion_both_sides():
    assert (2 * I(3, 4)).endpoints == (6.0, 8.0)
    assert (I(3, 4) * 2).endpoints == (6.0, 8.0)
    assert (10 - I(1, 2)).endpoints == (8.0, 9.0)
    assert (I(1, 2) + 0.5).endpoints == (1.5, 2.5)
    assert (1 / I(2, 4)).endpoints == (0.25, 0.5)


def test_scalar_promotion_equivalent_to_point():
    # A scalar operand is promoted to a point at its double value, so operating
    # with it equals operating with the explicit point interval. (The sum still
    # rounds outward, so point + point may itself be a nondegenerate interval.)
    assert (I(1) + 0.1) == (I(1) + I(0.1))
    assert (2 * I(3, 4)) == (I(2) * I(3, 4))
    assert (I(5) - 3) == I(2)


def test_non_finite_scalar_rejected():
    with pytest.raises(ValueError):
        _ = I(1, 2) + math.inf


def test_unary_and_abs():
    assert (-I(1, 2)).endpoints == (-2.0, -1.0)
    assert abs(I(-3, 2)).endpoints == (0.0, 3.0)
    assert (+I(1, 2)) == I(1, 2)


def test_integer_powers():
    assert (I(-2, 3) ** 2).endpoints == (0.0, 9.0)
    assert (I(-2, 3) ** 3).endpoints == (-8.0, 27.0)
    assert (I(2, 2) ** 10).endpoints == (1024.0, 1024.0)
    assert (I(2, 4) ** -1).endpoints == (0.25, 0.5)
    assert (I(2, 3) ** 0) == I(1, 1)


def test_float_and_real_powers():
    assert (I(2, 3) ** 2.0).endpoints == (4.0, 9.0)  # integral float -> pown
    assert 2**0.5 in (I(2) ** 0.5)  # real power via exp/log
    assert 2.0 in (I(4) ** 0.5)
    assert (I(2, 3) ** I(2)).endpoints == (4.0, 9.0)  # interval exponent


def test_set_operations():
    assert (I(1, 2) & I(1.5, 3)).endpoints == (1.5, 2.0)
    assert (I(1, 2) & I(5, 6)).is_empty
    assert (I(1, 2) | I(5, 6)).endpoints == (1.0, 6.0)
    assert I(1, 2).intersection(I(1.5, 3)).endpoints == (1.5, 2.0)
    assert I(1, 2).hull(I(5, 6)).endpoints == (1.0, 6.0)


def test_membership_and_containment():
    x = I(0.0, 1.0)
    assert 0.5 in x
    assert 2.0 not in x
    assert math.inf not in I(0.0, math.inf)  # infinity is never a member
    assert I(0.2, 0.8) in x  # subset via `in`
    assert x.contains(0.5)
    assert x.contains(I(0.2, 0.8))


def test_predicates():
    assert I(1, 2).subset(I(0, 3))
    assert I(0, 3).superset(I(1, 2))
    assert I(1, 2).is_interior_to(I(0, 3))
    assert not I(1, 2).is_interior_to(I(1, 3))
    assert I(1, 2).is_disjoint(I(3, 4))
    assert I(1, 2).overlaps(I(1.5, 3))
    assert I(1, 2).precedes(I(2, 3))
    assert not I(1, 2).strict_precedes(I(2, 3))
    assert I(1, 2).strict_precedes(I(3, 4))
    assert iv.Interval.empty().subset(I(1, 2))


def test_equality_with_non_interval_is_false():
    assert (I(1, 1) == 1.0) is False
    assert (I(1, 1) != 1.0) is True


def test_module_math_functions():
    assert iv.sqrt(I(4, 9)).endpoints == (2.0, 3.0)
    assert iv.sqrt(I(-4, 9)).endpoints == (0.0, 3.0)  # domain restricted
    assert iv.sqrt(I(-4, -1)).is_empty
    assert iv.sqr(I(-3, 2)).endpoints == (0.0, 9.0)
    assert iv.recip(I(2, 4)).endpoints == (0.25, 0.5)
    assert iv.fma(I(1, 2), I(3, 4), I(5, 6)).endpoints == (8.0, 14.0)
    assert iv.pown(I(2, 3), 2).endpoints == (4.0, 9.0)
    assert iv.floor(I(-1.5, 2.5)).endpoints == (-2.0, 2.0)
    assert iv.ceil(I(-1.5, 2.5)).endpoints == (-1.0, 3.0)
    assert iv.trunc(I(-1.7, 2.7)).endpoints == (-1.0, 2.0)
    assert iv.sign(I(-3, 2)).endpoints == (-1.0, 1.0)
    assert iv.min(I(1, 4), I(2, 3)).endpoints == (1.0, 3.0)
    assert iv.max(I(1, 4), I(2, 3)).endpoints == (2.0, 4.0)


def test_reverse_operations():
    lo, hi = iv.mul_rev(I(-1, 1), I(1, 2)).endpoints
    assert lo == -math.inf and hi == math.inf  # hull of two half-lines
    assert iv.sqr_rev(I(4, 9)).endpoints == (-3.0, 3.0)
    assert iv.sqr_rev(I(4, 9), I(0, 5)).endpoints == (2.0, 3.0)
    assert iv.abs_rev(I(1, 2)).endpoints == (-2.0, 2.0)


def test_cancellative_subtraction():
    z = iv.cancel_minus(I(2, 4), I(1, 2))
    assert z.endpoints == (1.0, 2.0)
    assert (I(1, 2) + z).subset(I(2, 4)) or I(2, 4).subset(I(1, 2) + z)
    assert iv.cancel_minus(I(1, 2), I(0, 3)).is_entire  # width violation


def test_empty_propagates_through_operations():
    e = iv.Interval.empty()
    assert (e + I(1, 2)).is_empty
    assert (I(1, 2) * e).is_empty
    assert iv.sqrt(e).is_empty
    assert (e**2).is_empty
