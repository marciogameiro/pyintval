"""Certifying that a map is defined and continuous on a box.

Many computer-assisted proofs in dynamics (isolation, Conley index, Morse graph
validity) require the map under study to be defined and continuous on each box
of a grid. A DecoratedInterval tracks this automatically: after evaluating an
arbitrary composed formula on a box, a decoration of 'dac' or 'com' is a
machine-checked certificate of "defined and continuous here" -- no per-function
hand analysis, and it holds through any composition.

Run with:  python examples/05_continuity_certificate.py
"""

import pyintval as iv

DI = iv.DecoratedInterval


def f(x: DI, y: DI) -> DI:
    # A map mixing a square root and a division -- both have domain conditions.
    #   f(x, y) = sqrt(x) / (y - 1) + log(x + y)
    return iv.sqrt(x) / (y - 1) + iv.log(x + y)


def certified(box_x, box_y) -> bool:
    x = DI(*box_x)
    y = DI(*box_y)
    result = f(x, y)
    ok = result.is_defined_and_continuous
    print(
        f"  f on x={box_x}, y={box_y}: enclosure {result.interval}, "
        f"decoration '{result.decoration}' -> continuous? {ok}"
    )
    return ok


def main() -> None:
    print("f(x, y) = sqrt(x) / (y - 1) + log(x + y)\n")

    print("Boxes where f is provably defined and continuous:")
    certified((1.0, 2.0), (2.0, 3.0))  # sqrt ok, y-1 > 0, x+y > 0
    certified((0.5, 4.0), (1.5, 5.0))

    print("\nBoxes where a domain condition may fail (decoration drops to 'trv'):")
    certified((-0.5, 2.0), (2.0, 3.0))  # sqrt of a possibly-negative x
    certified((1.0, 2.0), (0.5, 1.5))  # y - 1 straddles 0: division by zero
    certified((1.0, 2.0), (-3.0, -1.0))  # x + y can be <= 0: log undefined

    print(
        "\nThe certificate composes: one bad sub-expression poisons the whole "
        "result, so you never need to check the pieces by hand."
    )


if __name__ == "__main__":
    main()
