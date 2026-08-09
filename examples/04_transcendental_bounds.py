"""Rigorous bounds involving transcendental functions.

pyintval's elementary functions return guaranteed enclosures, so you can bound
and prove inequalities about expressions mixing exp, log, and trig. Each
single evaluation is at most a couple of ulps wider than the exact range -- but
when a variable appears more than once, naive evaluation is conservative (the
"dependency problem"), and subdividing the input sharpens the bound. This
example shows both effects.

Run with:  python examples/04_transcendental_bounds.py
"""

import pyintval as iv


def prove_positive(f, a, b, pieces=200):
    """Rigorously certify f(x) > 0 for all x in [a, b] by subdividing."""
    width = (b - a) / pieces
    worst = iv.entire()
    for i in range(pieces):
        piece = iv.Interval(a + i * width, a + (i + 1) * width)
        val = f(piece)
        worst = val if i == 0 else worst.hull(val)
        if val.lo <= 0.0:
            return False, worst
    return True, worst


def main() -> None:
    # Tight, certified enclosures of famous constants.
    print("pi in", iv.pi())
    print("e  in", iv.e())

    # sin(x)^2 + cos(x)^2 == 1. The enclosure contains 1, confirming the
    # identity -- but it is loose because x appears twice and interval
    # arithmetic treats the two occurrences as independent.
    x = iv.Interval(0.3, 0.9)
    s, c = iv.sin(x), iv.cos(x)
    print("\nsin(x)^2 + cos(x)^2 over x in", x)
    print("  enclosure:", s * s + c * c, " contains 1?", 1.0 in (s * s + c * c))

    # Prove exp(x) > x + 1 for all x in [0.1, 2]. A single interval evaluation
    # of expm1(x) - x is far too loose to conclude a sign because x appears
    # twice; refining the subdivision drives the enclosure to the true, strictly
    # positive range (its minimum near x = 0.1 is only ~0.005).
    def g(xi):
        return iv.expm1(xi) - xi  # = exp(x) - 1 - x

    single = g(iv.Interval(0.1, 2.0))
    print("\nexp(x) - (x + 1) over x in [0.1, 2]:")
    print("  single-shot enclosure:", single, " -> can conclude > 0?", single.lo > 0.0)
    for pieces in (50, 500, 5000):
        proven, refined = prove_positive(g, 0.1, 2.0, pieces)
        print(f"  {pieces:5d} pieces -> lower bound {refined.lo:+.6f}, proven > 0? {proven}")

    # Rigorously bound a transcendental expression over a box.
    box = iv.Interval(0.5, 1.5)
    f = iv.exp(-(box**2)) * iv.cos(box)
    print("\nf(x) = exp(-x^2) * cos(x) over x in", box)
    print("  guaranteed range enclosure:", f)


if __name__ == "__main__":
    main()
