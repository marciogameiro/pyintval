# pyintval

**Rigorous interval arithmetic for Python.** Every operation returns an interval
that is mathematically guaranteed to contain the true result, using correctly
rounded double-precision endpoints and IEEE 1788-2015 set-based semantics.

```{toctree}
:hidden:
:maxdepth: 2

guide
api
```

## Installation

```sh
pip install pyintval
```

Wheels are provided for Linux (x86-64, aarch64), macOS (arm64, x86-64), and
Windows (AMD64) on CPython 3.10–3.14. Building from source requires a C++20
compiler; on Windows that must be **clang-cl** (the vendored CORE-MATH kernels
use features MSVC lacks).

## Quickstart

```python
import pyintval as iv

# A string literal is parsed with correct OUTWARD rounding, so it provably
# encloses the exact decimal value (0.1 is not representable in binary64):
x = iv.Interval("0.1")
print(x.lo < x.hi)              # True: a nondegenerate 1-ulp enclosure of 1/10

# Arithmetic returns guaranteed enclosures.
y = iv.sqrt(x) + iv.Interval(2) * x
print(y.lo, y.hi)

# Set-based division never raises; it widens toward the unbounded result.
print(iv.Interval(1) / iv.Interval(-1, 1))     # Interval('[entire]')

# Elementary functions are correctly rounded, widened one ulp per side.
print(iv.exp(iv.log(iv.Interval(5))))          # encloses 5
print(3.141592653589793 in iv.pi())            # True

# A DecoratedInterval certifies "defined and continuous on this box".
d = iv.DecoratedInterval(1.0, 4.0)
print(iv.sqrt(d).decoration)                   # 'com'
print(iv.sqrt(iv.DecoratedInterval(-1.0, 4.0)).decoration)  # 'trv'
```

## What makes it rigorous

- **Correctly rounded arithmetic.** `+`, `-`, `*`, `/`, `sqrt`, and `fma` are
  computed with error-free transformations and directed one ulp outward — no
  rounding-mode switching, thread-safe, and cross-validated bit-for-bit against
  hardware directed rounding.
- **Correctly rounded elementary functions.** `exp`, `log`, trigonometric and
  hyperbolic functions and their inverses, `pow`, `atan2`, `hypot`, `erf`, and
  more are built on the [CORE-MATH](https://core-math.gitlabpages.inria.fr/)
  kernels and widened one ulp per endpoint, giving enclosures at most a couple
  of ulps wider than optimal. They are verified against a high-precision mpmath
  oracle over millions of inputs.
- **Set-based semantics.** Empty and unbounded intervals propagate instead of
  raising, so a single out-of-domain box never aborts a bulk computation.
- **Decorations.** An opt-in `DecoratedInterval` tracks an IEEE 1788 decoration
  that certifies, through any composition, whether the evaluated function is
  defined and continuous on its input — the hypothesis many computer-assisted
  proofs require.

See the {doc}`guide` for the guarantees in detail, or the {doc}`api` for the
full reference.
