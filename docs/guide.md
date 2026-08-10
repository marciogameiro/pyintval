# User guide

## The enclosure guarantee

An interval $[a, b]$ represents the set of all real numbers between $a$ and $b$
inclusive. For a function $f$, pyintval computes an interval that is guaranteed
to **contain** the true image of $f$ over the input:

$$\{\, f(x) : x \in [a, b] \,\} \subseteq \texttt{f}([a, b]).$$

This is the defining property of rigorous interval arithmetic. Because the
result is always a superset of the truth, any conclusion you can draw from the
interval (a sign, a bound, non-containment of a value) holds for the real
computation despite floating-point rounding.

## Correct rounding and directed endpoints

The endpoints are IEEE 754 double-precision numbers. To keep the enclosure
valid, the lower endpoint is always rounded **toward $-\infty$** and the upper
endpoint **toward $+\infty$** ("outward rounding"). pyintval achieves this
without ever changing the CPU rounding mode: it stays in round-to-nearest and
recovers the exact directed result from error-free transformations (TwoSum,
the FMA-based TwoProduct) plus a one-ulp `nextafter` adjustment. This makes the
primitives thread-safe and immune to compiler reordering, and they are
cross-checked bit-for-bit against the hardware's own directed rounding.

For the elementary functions, the correctly rounded CORE-MATH value is within
half an ulp of the true result, so widening it one ulp on each side yields a
rigorous enclosure at most about two ulps wide for a point input.

## Reading decimals rigorously

A Python `float` such as `0.1` is already the nearest double to the decimal
$1/10$; using it as an interval endpoint means "the double `0.1`", not the
exact decimal. To enclose the exact decimal value, pass a **string**:

```python
iv.Interval(0.1)      # the point {0.1-as-a-double}
iv.Interval("0.1")    # a 1-ulp interval that provably contains the real 1/10
```

The string parser converts each literal with correct directed rounding, using
exact big-integer comparison, so `"0.5"` yields the exact point `[0.5, 0.5]`
while `"0.1"` yields the tightest straddling interval.

## Text literal forms

`Interval(...)` accepts the IEEE 1788 text grammar, always rounding outward so
the parsed interval provably encloses the intended value:

| Form | Example | Result |
|---|---|---|
| decimal / hex-float point | `"0.1"`, `"0x1.8p1"` | tightest enclosing point |
| inf-sup | `"[1, 2]"`, `"[1.2e-3, 4]"` | `[RD(lo), RU(hi)]` |
| single number | `"[1.5]"` | `[1.5, 1.5]` |
| half-bounded | `"[1,]"`, `"[,2]"`, `"[,]"` | `[1, +inf]`, `[-inf, 2]`, entire |
| empty / entire | `"[empty]"`, `"[]"`, `"[entire]"` | the empty set / entire |
| uncertain | `"3.56?1"`, `"-10?u"` | midpoint ± radius in ulps |

The **uncertain form** `m?rad` means the midpoint `m` plus or minus `rad` units
in `m`'s last decimal place (half a unit if `rad` is omitted); a trailing `u`/`d`
makes it one-sided, and a trailing exponent scales the whole thing:

```python
iv.Interval("3.56?1")      # [3.55, 3.57]
iv.Interval("-10?")        # [-10.5, -9.5]
iv.Interval("3.56?1e2")    # [355, 357]
```

A `DecoratedInterval` literal may carry an explicit `_dec` suffix (case
insensitive); the decoration is validated — it may not over-claim — and `"[nai]"`
or an invalid decoration yields NaI. With no suffix the tightest decoration is
assigned:

```python
iv.DecoratedInterval("[1, 2]_com")   # [1, 2]_com
iv.DecoratedInterval("[entire]")     # [entire]_dac
iv.DecoratedInterval("[empty]_com")  # NaI — empty cannot be 'com'
```

## Set-based semantics

pyintval implements the IEEE 1788-2015 *set-based* flavor. Operations never
raise mid-computation; instead:

- The **empty set** models an impossible/undefined result and propagates through
  every subsequent operation (`iv.sqrt(iv.Interval(-4, -1))` is empty).
- **Unbounded** intervals (half-lines, or the whole line `[entire]`) model
  results that grow without bound; e.g. dividing by an interval that contains
  zero yields an unbounded result rather than a `ZeroDivisionError`.
- A domain that is only partly violated is **clipped**: `iv.sqrt(iv.Interval(-4, 9))`
  returns `[0, 3]`.

This composability is what lets you sweep a map over a large grid of boxes
without a single bad box aborting the run.

## The dependency problem

When a variable appears more than once in an expression, interval arithmetic
treats the occurrences as independent, so the result can be wider than the true
range. For example, over $x \in [-1, 3]$:

```python
x = iv.Interval(-1, 3)
print(x**2 - 2*x + 1)     # [-5, 12], much wider than the true [0, 16]
```

The enclosure is still rigorous — it contains the true range — just not tight.
Subdividing the input and taking the hull of the pieces sharpens the bound, and
is the standard way to drive interval methods to a conclusion:

```python
hull = iv.Interval.empty()
lo, n = -1.0, 200
w = 4.0 / n
for i in range(n):
    piece = iv.Interval(lo + i*w, lo + (i+1)*w)
    hull = hull.hull(piece**2 - 2*piece + 1)
print(hull)               # much closer to [0, 16]
```

## Decorations: certifying continuity

A `DecoratedInterval` pairs an interval with an IEEE 1788 *decoration* recording
the strongest property known of the function evaluated so far:

| Decoration | Meaning |
|---|---|
| `com` | defined, continuous, and bounded on a common (bounded, nonempty) input |
| `dac` | defined and continuous on the input |
| `def` | defined, but continuity not guaranteed (e.g. `floor`, `sign` over a step) |
| `trv` | only the enclosure is guaranteed; a domain violation may have occurred |
| `ill` | not an interval (poisons every computation) |

Decorations propagate as the minimum of each operation's local decoration and
its inputs', so after evaluating an arbitrary composed formula on a box, a
result of `dac` or `com` is a machine-checked certificate that the whole
expression is **defined and continuous on that box** — with no per-function
hand analysis:

```python
D = iv.DecoratedInterval
def f(x, y):
    return iv.sqrt(x) / (y - 1) + iv.log(x + y)

print(f(D(1, 2), D(2, 3)).is_defined_and_continuous)   # True
print(f(D(1, 2), D(0.5, 1.5)).decoration)              # 'trv': y-1 straddles 0
```

This is exactly the hypothesis needed by many computer-assisted proofs
(isolation, Conley index, Morse graph validity).

## Standards conformance

pyintval implements the set-based flavor of **IEEE Std 1788-2015**. It is
validated against [ITF1788](https://github.com/oheim/ITF1788), the reference
conformance suite for the standard (~7,200 cases): every result is checked to
*enclose* the standard's tightest interval, and each release gates on that rigor
property in continuous integration. The transcendental functions are correctly
rounded and then widened one ulp per endpoint, so they always enclose but are
intentionally not bit-for-bit *tightest*.
