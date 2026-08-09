"""Rigorous enclosure tests for the elementary functions, against a
high-precision mpmath oracle.

For each function and many deterministically-seeded input intervals, the true
mathematical function is sampled across the interval at ~200-digit precision;
every sample must lie within the computed interval result. This is the
authoritative check that the interval transcendentals are rigorous enclosures.
Structural/known-value checks live in the C++ suite.
"""

import math

import pytest

import pyintval as iv

mp = pytest.importorskip("mpmath")
mp.mp.prec = 220

I = iv.Interval
INF = math.inf


def _encloses(out: iv.Interval, value) -> bool:
    v = mp.mpf(value)
    lo = -mp.inf if out.lo == -INF else mp.mpf(out.lo)
    hi = mp.inf if out.hi == INF else mp.mpf(out.hi)
    return lo <= v <= hi


def _sample(a: float, b: float, n: int):
    # mpf arithmetic: k=0 -> exactly a, k=n-1 -> exactly b, all strictly inside.
    lo, hi = mp.mpf(a), mp.mpf(b)
    return [lo + (hi - lo) * mp.mpf(k) / (n - 1) for k in range(n)]


def _real_cbrt(t):
    return mp.sign(t) * mp.cbrt(abs(t))  # mp.cbrt gives the complex principal root


# (name, interval-fn, mpmath-fn, safe-domain-lo, safe-domain-hi)
UNARY = [
    ("exp", iv.exp, mp.exp, -700.0, 700.0),
    ("expm1", iv.expm1, mp.expm1, -700.0, 700.0),
    ("exp2", iv.exp2, lambda t: mp.power(2, t), -1000.0, 1000.0),
    ("exp10", iv.exp10, lambda t: mp.power(10, t), -300.0, 300.0),
    ("cbrt", iv.cbrt, _real_cbrt, -1e6, 1e6),
    ("sinh", iv.sinh, mp.sinh, -700.0, 700.0),
    ("cosh", iv.cosh, mp.cosh, -700.0, 700.0),
    ("tanh", iv.tanh, mp.tanh, -1e6, 1e6),
    ("asinh", iv.asinh, mp.asinh, -1e6, 1e6),
    ("atan", iv.atan, mp.atan, -1e6, 1e6),
    ("erf", iv.erf, mp.erf, -1e3, 1e3),
    ("erfc", iv.erfc, mp.erfc, -100.0, 100.0),
    ("sin", iv.sin, mp.sin, -1e3, 1e3),
    ("cos", iv.cos, mp.cos, -1e3, 1e3),
    ("tan", iv.tan, mp.tan, -1e3, 1e3),
    ("log", iv.log, mp.log, 1e-300, 1e300),
    ("log2", iv.log2, lambda t: mp.log(t, 2), 1e-300, 1e300),
    ("log10", iv.log10, mp.log10, 1e-300, 1e300),
    ("log1p", iv.log1p, mp.log1p, -0.9999, 1e300),
    ("asin", iv.asin, mp.asin, -1.0, 1.0),
    ("acos", iv.acos, mp.acos, -1.0, 1.0),
    ("atanh", iv.atanh, mp.atanh, -0.9999, 0.9999),
    ("acosh", iv.acosh, mp.acosh, 1.0, 1e6),
]


def _gen_interval(rng, dlo, dhi):
    """A random [lo, hi] strictly inside the safe domain [dlo, dhi]."""
    c = rng.uniform(dlo, dhi)
    span = dhi - dlo
    w = span * 10 ** rng.uniform(-12, -0.5)
    lo = max(c - w * rng.random(), dlo)
    hi = min(c + w * rng.random(), dhi)
    if lo > hi:
        lo, hi = hi, lo
    return lo, hi


@pytest.mark.parametrize("name,ivf,mpf,dlo,dhi", UNARY, ids=[u[0] for u in UNARY])
def test_unary_enclosure(name, ivf, mpf, dlo, dhi):
    import random

    rng = random.Random(1000 + UNARY.index((name, ivf, mpf, dlo, dhi)))
    for _ in range(500):
        lo, hi = _gen_interval(rng, dlo, dhi)
        out = ivf(I(lo, hi))
        assert not out.is_empty
        # tan asymptotes legitimately yield the entire line; nothing to enclose.
        if out.is_entire:
            continue
        for p in _sample(lo, hi, 15):
            assert _encloses(out, mpf(p)), f"{name}({float(p)!r}) not in {out}"


@pytest.mark.parametrize("name,ivf,mpf", [
    ("sin", iv.sin, mp.sin), ("cos", iv.cos, mp.cos),
], ids=["sin", "cos"])
def test_periodic_dense(name, ivf, mpf):
    import random

    rng = random.Random(0xC05 + (name == "cos"))
    for _ in range(2000):
        a = rng.uniform(-30, 30)
        b = a + 10 ** rng.uniform(-4, 1.2)
        out = ivf(I(a, b))
        for p in _sample(a, b, 60):  # dense: must catch interior peaks/troughs
            assert _encloses(out, mpf(p)), f"{name} peak missed on [{a},{b}] -> {out}"


def test_two_argument_enclosure():
    import random

    rng = random.Random(0xA7A2)
    for _ in range(3000):
        x1, x2 = sorted(rng.uniform(-30, 30) for _ in range(2))
        y1, y2 = sorted(rng.uniform(-30, 30) for _ in range(2))
        hy = iv.hypot(I(x1, x2), I(y1, y2))
        at = iv.atan2(I(y1, y2), I(x1, x2))
        for xx in _sample(x1, x2, 6):
            for yy in _sample(y1, y2, 6):
                assert _encloses(hy, mp.hypot(xx, yy))
                assert _encloses(at, mp.atan2(yy, xx))  # mpmath order is atan2(y, x)
        # pow with strictly positive base
        bx1, bx2 = sorted(rng.uniform(0.01, 20) for _ in range(2))
        pw = iv.pow(I(bx1, bx2), I(y1, y2))
        if pw.is_entire:
            continue
        for xx in _sample(bx1, bx2, 5):
            for yy in _sample(y1, y2, 5):
                assert _encloses(pw, mp.power(xx, yy))


def test_tightness_at_points():
    # A point interval must yield an enclosure only a couple of ulps wide.
    for f, val in [(iv.exp, 1.3), (iv.log, 3.7), (iv.sin, 0.9), (iv.cos, 2.1),
                   (iv.atan, 5.0), (iv.tanh, 0.6), (iv.cbrt, 9.0), (iv.erf, 0.4)]:
        out = f(I(val))
        assert out.wid <= 4 * math.ulp(max(abs(out.hi), abs(out.lo), 1e-300))


def test_known_identities():
    for x in [0.5, 1.0, 2.0, 10.0, 100.0]:
        assert x in iv.exp(iv.log(I(x)))
    for t in [0.3, 1.0, 2.5, 4.0]:
        s = iv.sin(I(t))
        c = iv.cos(I(t))
        assert 1.0 in (s * s + c * c)
    assert math.e in iv.e()
    assert math.pi in iv.pi()


def test_constants_are_tight():
    assert iv.pi().wid <= 4 * math.ulp(3.14159)
    assert iv.e().wid <= 4 * math.ulp(2.71828)
