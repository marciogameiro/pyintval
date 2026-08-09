"""Rigorous interval arithmetic for Python.

pyintval implements IEEE 1788-2015 set-based interval arithmetic with
correctly rounded double-precision endpoints: every operation returns an
interval mathematically guaranteed to contain the true result.

Example
-------
>>> import pyintval as iv
>>> x = iv.Interval("0.1")          # rigorously encloses the decimal 0.1
>>> y = iv.sqrt(x) + iv.Interval(2) * x
>>> y.lo <= y.hi
True
>>> iv.Interval(1) / iv.Interval(-1, 1)   # division straddling zero
Interval('[entire]')
"""

from __future__ import annotations

from pyintval._core import (
    KERNEL_ABI_VERSION,
    Interval,
    __version__,
    abs,
    abs_rev,
    build_info,
    cancel_minus,
    cancel_plus,
    ceil,
    empty,
    entire,
    floor,
    fma,
    hull,
    intersection,
    max,
    min,
    mul_rev,
    pi,
    pown,
    recip,
    round,
    round_ties_to_away,
    sign,
    sqr,
    sqr_rev,
    sqrt,
    trunc,
)

__all__ = [
    "Interval",
    "KERNEL_ABI_VERSION",
    "__version__",
    "abs",
    "abs_rev",
    "build_info",
    "cancel_minus",
    "cancel_plus",
    "ceil",
    "empty",
    "entire",
    "floor",
    "fma",
    "hull",
    "intersection",
    "max",
    "min",
    "mul_rev",
    "pi",
    "pown",
    "recip",
    "round",
    "round_ties_to_away",
    "sign",
    "sqr",
    "sqr_rev",
    "sqrt",
    "trunc",
]
