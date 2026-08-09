"""Property-based tests (hypothesis): the fundamental interval guarantees,
exercised through the Python binding layer.

The core contract of a rigorous interval library is the *containment property*:
for an operation f and inputs X, Y, the interval result f(X, Y) contains f(x, y)
for every x in X, y in Y. These tests sample members and assert membership,
plus algebraic identities and inclusion monotonicity.
"""

import math
import operator

from hypothesis import assume, given
from hypothesis import strategies as st

import pyintval as iv

I = iv.Interval

# Finite doubles spanning a wide dynamic range, including subnormals and zero.
finite = st.floats(allow_nan=False, allow_infinity=False,
                   min_value=-1e100, max_value=1e100)


@st.composite
def intervals(draw, allow_unbounded=True):
    a = draw(finite)
    b = draw(finite)
    lo, hi = (a, b) if a <= b else (b, a)
    if allow_unbounded:
        kind = draw(st.integers(min_value=0, max_value=3))
        if kind == 1:
            lo = -math.inf
        elif kind == 2:
            hi = math.inf
        elif kind == 3 and draw(st.booleans()):
            return iv.Interval.entire()
    return I(lo, hi)


def members(x, rng_vals):
    """A few representative members of a (nonempty, bounded-enough) interval."""
    out = []
    if math.isfinite(x.lo):
        out.append(x.lo)
    if math.isfinite(x.hi):
        out.append(x.hi)
    m = x.mid
    if math.isfinite(m):
        out.append(m)
    for t in rng_vals:
        if math.isfinite(x.lo) and math.isfinite(x.hi):
            out.append(x.lo + t * (x.hi - x.lo))
    return [v for v in out if v in x]


@given(intervals(), intervals(), st.lists(st.floats(0, 1), min_size=1, max_size=3))
def test_containment_add_sub_mul(x, y, ts):
    for op, f in [(operator.add, iv.Interval.__add__),
                  (operator.sub, iv.Interval.__sub__),
                  (operator.mul, iv.Interval.__mul__)]:
        r = f(x, y)
        for a in members(x, ts):
            for b in members(y, ts):
                val = op(a, b)
                if math.isfinite(val):
                    assert val in r


@given(intervals(), intervals(), st.lists(st.floats(0, 1), min_size=1, max_size=3))
def test_containment_division(x, y, ts):
    r = x / y
    if r.is_empty:
        return
    for a in members(x, ts):
        for b in members(y, ts):
            if b != 0.0:
                val = a / b
                if math.isfinite(val):
                    assert val in r


@given(intervals(allow_unbounded=False), st.integers(min_value=0, max_value=8),
       st.lists(st.floats(0, 1), min_size=1, max_size=3))
def test_containment_pown(x, n, ts):
    r = x ** n
    for a in members(x, ts):
        try:
            val = a ** n
        except OverflowError:
            continue  # pointwise value exceeds DBL_MAX; r is unbounded there
        if math.isfinite(val):
            assert val in r


@given(intervals())
def test_sqrt_contains_pointwise(x):
    r = iv.sqrt(x)
    for a in (x.lo, x.hi, x.mid):
        if math.isfinite(a) and a >= 0.0:
            assert math.sqrt(a) in r


@given(intervals())
def test_neg_involution_and_sub_identity(x):
    assert -(-x) == x
    y = I(-1.0, 2.0)
    assert (x - y) == (x + (-y))


@given(intervals(), intervals())
def test_inclusion_isotonicity_add(x, y):
    # Shrinking either operand cannot enlarge the result.
    if x.is_empty or y.is_empty:
        return
    xm, ym = x.mid, y.mid
    if not (math.isfinite(xm) and math.isfinite(ym)):
        return
    xs, ys = I(xm, xm), I(ym, ym)
    assert (xs + ys).subset(x + y)


@given(intervals(), intervals())
def test_intersection_hull_relations(x, y):
    inter = x & y
    hull = x | y
    assert inter.subset(x)
    assert inter.subset(y)
    if not x.is_empty:
        assert x.subset(hull)
    if not y.is_empty:
        assert y.subset(hull)


@given(finite)
def test_string_roundtrip_encloses(v):
    p = I(v)
    reparsed = I(str(p))
    assert p.subset(reparsed)


@given(st.text(min_size=0, max_size=6))
def test_arbitrary_text_never_crashes(s):
    try:
        I(s)
    except ValueError:
        pass  # malformed input must raise cleanly, never crash


@given(intervals())
def test_mid_is_a_member_when_nonempty(x):
    if x.is_empty:
        return
    # mid returns a finite double guaranteed to lie inside the interval, even
    # for unbounded intervals (where it is the largest-magnitude finite double).
    assert x.mid in x


@given(intervals(allow_unbounded=False))
def test_bounded_mid_rad_cover(x):
    if x.is_empty:
        return
    m, r = x.mid, x.rad
    assert m - r <= x.lo + 1e-300 or math.isclose(m - r, x.lo, rel_tol=1e-12, abs_tol=1e-300)
    assert x.hi - 1e-300 <= m + r or math.isclose(m + r, x.hi, rel_tol=1e-12, abs_tol=1e-300)
