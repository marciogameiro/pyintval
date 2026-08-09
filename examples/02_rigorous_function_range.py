"""Rigorously enclose the range of a function over a box.

Evaluating a function in interval arithmetic yields a guaranteed superset of
its true range: every real output for every real input in the box is contained
in the interval result. This is the workhorse of computer-assisted proofs and
rigorous dynamics (e.g. showing a map cannot vanish, or bounding its image).

Run with:  python examples/02_rigorous_function_range.py
"""

import pyintval as iv


def f(x: iv.Interval) -> iv.Interval:
    # f(x) = x^2 - 2x + 1 = (x - 1)^2 >= 0 for all real x.
    return x**2 - 2 * x + 1


def main() -> None:
    box = iv.Interval(-1, 3)
    enclosure = f(box)
    print(f"f(x) = x^2 - 2x + 1 over x in {box}")
    print(f"  guaranteed range enclosure: {enclosure}")
    print(f"  proven nonnegative on the box? {enclosure.lo >= 0.0}")

    # A tighter enclosure by subdividing the box (interval arithmetic is
    # subdistributive, so smaller boxes give sharper bounds).
    n = 100
    lo = box.lo
    width = (box.hi - box.lo) / n
    hull = iv.Interval.empty()
    for i in range(n):
        piece = iv.Interval(lo + i * width, lo + (i + 1) * width)
        hull = hull.hull(f(piece))
    print(f"  refined enclosure ({n} subdivisions): {hull}")
    print(f"  refined lower bound closer to the true min 0: {hull.lo:.6f}")


if __name__ == "__main__":
    main()
