"""Tightness (overestimation-control) tests for *composed* interval expressions.

The enclosure tests elsewhere prove a result CONTAINS the truth. These tests
prove a result is not needlessly WIDE. For a fixed evaluation form -- a fixed
sequence of interval operations -- a correct kernel must return the *tightest*
correctly-rounded interval that sequence can produce. That interval is unique,
and we compute it independently with exact rational arithmetic
(``fractions.Fraction``) plus directed rounding, then assert pyintval
reproduces it bit-for-bit (``Interval.endpoints``).

The same cases document the two irreducible sources of interval overestimation
that *no* library can remove:

* the **dependency problem** -- a variable occurring more than once is treated
  as independent copies (``x - x`` is not ``0``; ``x*(1-x)`` beats ``x - x*x``),
* the **wrapping effect** -- an axis-aligned box cannot represent a rotated box,
  so error accumulates under iteration.

The point is not that pyintval makes these tight (it cannot -- that is inherent
to interval arithmetic on boxes), but that (a) it adds no slop beyond the
minimum the chosen form allows, and (b) it still rigorously encloses the truth.
"""

import math
from fractions import Fraction

import pyintval as iv

I = iv.Interval


# --------------------------------------------------------------------------- #
# Exact "tightest correctly-rounded" reference interval evaluator.
#
# Every endpoint here is a double and every operation is +, -, *, / or an
# integer power -- all exact over the rationals. Rounding each exact endpoint
# OUTWARD to the nearest double yields the unique minimal interval a correct
# kernel must return for that operation sequence: the tightness oracle. Because
# the reference replays the *same* sequence, a bit-for-bit match certifies the
# kernel introduces no excess width; any mismatch is a real tightness defect.
# --------------------------------------------------------------------------- #


def _rd(q: Fraction) -> float:
    """Largest double <= q (round toward -inf)."""
    d = float(q)
    return d if Fraction(d) <= q else math.nextafter(d, -math.inf)


def _ru(q: Fraction) -> float:
    """Smallest double >= q (round toward +inf)."""
    d = float(q)
    return d if Fraction(d) >= q else math.nextafter(d, math.inf)


class Ref:
    """Minimal correctly-rounded interval for a fixed operation sequence."""

    __slots__ = ("hi", "lo")

    def __init__(self, lo, hi):
        self.lo = float(lo)
        self.hi = float(hi)

    @classmethod
    def _c(cls, o):
        return o if isinstance(o, Ref) else cls(o, o)

    def _q(self):
        return Fraction(self.lo), Fraction(self.hi)

    def __add__(self, o):
        (a, b), (c, d) = self._q(), Ref._c(o)._q()
        return Ref(_rd(a + c), _ru(b + d))

    __radd__ = __add__

    def __sub__(self, o):
        (a, b), (c, d) = self._q(), Ref._c(o)._q()
        return Ref(_rd(a - d), _ru(b - c))

    def __rsub__(self, o):
        return Ref._c(o).__sub__(self)

    def __mul__(self, o):
        (a, b), (c, d) = self._q(), Ref._c(o)._q()
        corners = (a * c, a * d, b * c, b * d)
        return Ref(_rd(min(corners)), _ru(max(corners)))

    __rmul__ = __mul__

    def __truediv__(self, o):
        (a, b), (c, d) = self._q(), Ref._c(o)._q()
        assert c > 0 or d < 0, "reference divisor must not straddle zero"
        corners = (a / c, a / d, b / c, b / d)
        return Ref(_rd(min(corners)), _ru(max(corners)))

    def __neg__(self):
        return Ref(-self.hi, -self.lo)

    def sqr(self):
        a, b = self._q()
        lo = Fraction(0) if a < 0 < b else min(a * a, b * b)
        return Ref(_rd(lo), _ru(max(a * a, b * b)))

    def pown(self, n):
        a, b = self._q()
        if n % 2 == 1:  # monotone increasing
            return Ref(_rd(a**n), _ru(b**n))
        lo = Fraction(0) if a < 0 < b else min(a**n, b**n)
        return Ref(_rd(lo), _ru(max(a**n, b**n)))

    @property
    def endpoints(self):
        return (self.lo, self.hi)


def _match(result: iv.Interval, ref: Ref, label: str):
    """The kernel result must equal the tightest correctly-rounded interval."""
    assert result.endpoints == ref.endpoints, (
        f"{label}: pyintval {result.endpoints} is not the tightest "
        f"correctly-rounded result {ref.endpoints}"
    )


# --------------------------------------------------------------------------- #
# Dependency problem
# --------------------------------------------------------------------------- #


def test_x_minus_x_is_dependency_not_zero():
    # The two occurrences of x are independent copies, so the tightest result of
    # THIS form is [lo-hi, hi-lo], never the true range {0}. Assert the kernel is
    # exactly that tight (no wider) and still encloses 0.
    for lo, hi in [(0.0, 1.0), (-2.0, 3.0), (0.1, 0.2), (1.5, 1.5)]:
        x = I(lo, hi)
        _match(x - x, Ref(lo, hi) - Ref(lo, hi), f"x-x on [{lo},{hi}]")
        assert 0.0 in (x - x)
    # A point input has no dependency: the result collapses to exactly {0}.
    assert (I(2.5) - I(2.5)).endpoints == (0.0, 0.0)
    # A width-2 input yields width 2 -- the full dependency overestimation.
    assert (I(0.0, 1.0) - I(0.0, 1.0)).endpoints == (-1.0, 1.0)


def test_logistic_two_forms():
    # x*(1-x) and x - x*x are algebraically identical; on [0,1] the true range is
    # [0, 1/4], yet the two forms give different, both-loose enclosures.
    x, xr = I(0.0, 1.0), Ref(0.0, 1.0)

    form1 = x * (I(1.0) - x)  # (1-x) in [0,1] -> product [0,1]: 4x too wide up top
    _match(form1, xr * (Ref(1.0, 1.0) - xr), "x*(1-x)")
    assert form1.endpoints == (0.0, 1.0)

    form2 = x - x * x  # [0,1] - [0,1] = [-1,1]: 8x too wide, wrong sign below
    _match(form2, xr - xr * xr, "x - x*x")
    assert form2.endpoints == (-1.0, 1.0)

    # Both rigorously enclose the exact range [0, 1/4]; form1 is strictly tighter.
    for f in (form1, form2):
        assert f.lo <= 0.0 and 0.25 <= f.hi
    assert form1.subset(form2)


# --------------------------------------------------------------------------- #
# Evaluation-order sensitivity (same polynomial, different rearrangements)
# --------------------------------------------------------------------------- #


def test_quartic_evaluation_order():
    # p(x) = x^4 - 4x^3 + 4x^2 = x^2 (x-2)^2 on [0,3]; exact range [0, 9].
    x, xr = I(0.0, 3.0), Ref(0.0, 3.0)

    naive = iv.pown(x, 4) - 4 * iv.pown(x, 3) + 4 * iv.pown(x, 2)
    _match(naive, xr.pown(4) - 4 * xr.pown(3) + 4 * xr.pown(2), "quartic naive")
    assert naive.endpoints == (-108.0, 117.0)

    horner = (((x - 4) * x + 4) * x) * x
    _match(horner, (((xr - 4) * xr + 4) * xr) * xr, "quartic Horner")
    assert horner.endpoints == (-72.0, 36.0)

    factored = iv.sqr(x) * iv.sqr(x - 2)
    _match(factored, xr.sqr() * (xr - 2).sqr(), "quartic factored")
    assert factored.endpoints == (0.0, 36.0)

    # Tighter rearrangements nest inside looser ones; all enclose [0, 9].
    assert factored.subset(horner)
    assert horner.subset(naive)
    for f in (naive, horner, factored):
        assert f.lo <= 0.0 and 9.0 <= f.hi


def test_cubic_horner_reaches_exact_range():
    # x^3 - x - 1 on [1,2]; f increases there, so the exact range is [-1, 5].
    # Horner attains it exactly; the naive monomial form does not.
    x, xr = I(1.0, 2.0), Ref(1.0, 2.0)

    naive = iv.pown(x, 3) - x - 1
    _match(naive, xr.pown(3) - xr - 1, "cubic naive")
    assert naive.endpoints == (-2.0, 6.0)

    horner = (x * x - 1) * x - 1
    _match(horner, (xr * xr - 1) * xr - 1, "cubic Horner")
    assert horner.endpoints == (-1.0, 5.0)  # == the exact range

    assert horner.subset(naive)
    assert 0.0 in naive and 0.0 in horner  # both enclose the root (plastic number)


# --------------------------------------------------------------------------- #
# Catastrophic cancellation: enclosure must survive where floats do not
# --------------------------------------------------------------------------- #


def test_rump_encloses_truth_despite_cancellation():
    # Rump's example. Naive float64 returns a large positive number (wrong sign);
    # the true value is -54767/66192 ~ -0.827396. A rigorous kernel must still
    # ENCLOSE the truth -- and honestly reports a huge width here, because naive
    # double-precision interval evaluation cannot even certify the sign.
    a, b = 77617.0, 33096.0
    A, B = I(a), I(b)
    b2, b4, b6, b8, a2 = (iv.pown(B, 2), iv.pown(B, 4), iv.pown(B, 6), iv.pown(B, 8), iv.pown(A, 2))
    f = (
        I(333.75) * b6
        + a2 * (I(11.0) * a2 * b2 - b6 - I(121.0) * b4 - I(2.0))
        + I(5.5) * b8
        + A / (I(2.0) * B)
    )

    fa, fb = Fraction(a), Fraction(b)
    true = (
        Fraction(1335, 4) * fb**6
        + fa**2 * (11 * fa**2 * fb**2 - fb**6 - 121 * fb**4 - 2)
        + Fraction(11, 2) * fb**8
        + fa / (2 * fb)
    )
    assert true == Fraction(-54767, 66192)  # the classic exact value

    # Rigorous enclosure of the exact value ...
    assert Fraction(f.lo) <= true <= Fraction(f.hi)
    # ... but the enclosure is uselessly wide, not a false point estimate.
    assert f.wid > 1.0

    # Contrast: naive float64 is nowhere near the truth (catastrophic loss).
    naive_float = (
        333.75 * b**6 + a**2 * (11 * a**2 * b**2 - b**6 - 121 * b**4 - 2) + 5.5 * b**8 + a / (2 * b)
    )
    assert abs(naive_float - float(true)) > 1.0


# --------------------------------------------------------------------------- #
# Wrapping effect
# --------------------------------------------------------------------------- #


def test_wrapping_effect_rotation():
    # Rotate the box [-1,1]^2 by 45 degrees eight times: eight 45-degree turns is
    # a full revolution, so the true image is again exactly [-1,1]^2. But each
    # interval rotation must hull the rotated diamond, so a naive box grows by
    # ~sqrt(2) per axis per step -> ~[-16,16]^2 after eight steps. That growth is
    # the wrapping effect; it is irreducible for plain interval boxes.
    s = iv.recip(iv.sqrt(I(2.0)))  # 1/sqrt(2), a ~1-ulp interval
    sr = Ref(s.lo, s.hi)  # reuse the SAME constant so we test the rotation, not sqrt

    x = y = I(-1.0, 1.0)
    xr = yr = Ref(-1.0, 1.0)
    for _ in range(8):
        x, y = (x - y) * s, (x + y) * s
        xr, yr = (xr - yr) * sr, (xr + yr) * sr
        _match(x, xr, "rotation x")
        _match(y, yr, "rotation y")

    # The truth (full turn = original box) stays enclosed throughout.
    assert I(-1.0, 1.0).subset(x)
    assert I(-1.0, 1.0).subset(y)
    # Wrapping actually occurred: far wider than the true width 2 ...
    assert x.wid > 8.0
    # ... yet bounded by the theoretical (sqrt2)^8 = 16x growth per axis.
    assert x.wid < 33.0 and y.wid < 33.0
