"""Basic usage of pyintval: constructing intervals and doing arithmetic.

Run with:  python examples/01_basic_usage.py
"""

import pyintval as iv


def main() -> None:
    # Construct from bounds, from a number, or from a string. A string literal
    # is parsed with correct OUTWARD rounding, so it provably encloses the exact
    # decimal value -- note that 0.1 is not representable in binary64:
    a = iv.Interval(1, 2)
    b = iv.Interval("0.1")
    print("a        =", a)
    print("b        =", b, "  (a nondegenerate enclosure of the real 1/10)")
    print("b.lo < b.hi:", b.lo < b.hi)

    # Arithmetic returns intervals guaranteed to contain the true result.
    print("a + b    =", a + b)
    print("a * b    =", a * b)
    print("a ** 3   =", a**3)
    print("2 * a - 1 =", 2 * a - 1)  # scalars are promoted to point intervals

    # Set-based division never raises: dividing by an interval that straddles
    # zero yields an unbounded result instead of a ZeroDivisionError.
    print("a / [-1, 1] =", a / iv.Interval(-1, 1))
    print("a / [0, 2]  =", a / iv.Interval(0, 2))

    # Query the interval.
    print("mid, rad, wid of a:", a.mid, a.rad, a.wid)
    print("float 0.1 in Interval('0.1')?", 0.1 in iv.Interval("0.1"))
    print("a subset of [0, 3]?", a.subset(iv.Interval(0, 3)))


if __name__ == "__main__":
    main()
