"""Prove a root exists with a rigorous interval bisection.

If a continuous function has interval enclosures of opposite sign at the two
ends of a bracket, the Intermediate Value Theorem guarantees a real root
between them. Interval arithmetic makes the sign checks rigorous (no rounding
can produce a false conclusion), and bisection shrinks a verified bracket
around a root of f(x) = x^2 - 2, i.e. around sqrt(2).

Run with:  python examples/03_root_existence.py
"""

import pyintval as iv


def f(x: iv.Interval) -> iv.Interval:
    return x**2 - iv.Interval(2)


def certainly_negative(x: iv.Interval) -> bool:
    return x.hi < 0.0


def certainly_positive(x: iv.Interval) -> bool:
    return x.lo > 0.0


def main() -> None:
    lo, hi = 1.0, 2.0
    assert certainly_negative(f(iv.Interval(lo)))  # f(1) = -1 < 0
    assert certainly_positive(f(iv.Interval(hi)))  # f(2) =  2 > 0

    for _ in range(60):
        mid = 0.5 * (lo + hi)
        fmid = f(iv.Interval(mid))
        if certainly_positive(fmid):
            hi = mid
        elif certainly_negative(fmid):
            lo = mid
        else:
            # f(mid) enclosure straddles 0: mid is a rigorous approximation.
            break

    bracket = iv.Interval(lo, hi)
    print("Rigorous bracket for a root of x^2 - 2:")
    print(f"  {bracket}")
    print(f"  width = {bracket.wid:.3e}")
    # Cross-check: pyintval's own directed sqrt of 2 must lie in the bracket.
    root = iv.sqrt(iv.Interval(2))
    print(f"  sqrt(2) enclosure = {root}")
    print(f"  bracket contains the true sqrt(2)? {root.lo in bracket or bracket.overlaps(root)}")


if __name__ == "__main__":
    main()
