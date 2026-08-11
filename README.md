# pyintval

[![CI](https://github.com/marciogameiro/pyintval/actions/workflows/ci.yml/badge.svg)](https://github.com/marciogameiro/pyintval/actions/workflows/ci.yml)
[![Docs](https://readthedocs.org/projects/pyintval/badge/?version=latest)](https://pyintval.readthedocs.io/)
[![PyPI](https://img.shields.io/pypi/v/pyintval)](https://pypi.org/project/pyintval/)

**Rigorous interval arithmetic for Python** — every operation returns an
interval that is mathematically guaranteed to contain the true result.

> **Released on PyPI**, with prebuilt wheels for Linux, macOS, and Windows on
> CPython 3.10–3.14. Validated against the [ITF1788](https://github.com/oheim/ITF1788)
> reference conformance suite for IEEE 1788-2015.

## What it provides

- IEEE 1788-2015 **set-based interval arithmetic** over double-precision
  endpoints: empty and unbounded intervals propagate instead of raising, so
  bulk computations never abort mid-sweep.
- **Correctly rounded endpoints** for `+`, `-`, `*`, `/`, `sqrt`, `fma` via
  error-free transformations (no rounding-mode switching, thread-safe), and
  **elementary functions** (`exp`, `log`, `sin`, `cos`, `tan`, their inverses,
  hyperbolics, `pow`, `atan2`, `hypot`, `cbrt`, `expm1`, `log1p`, `erf`,
  `erfc`, ...) built on the correctly rounded
  [CORE-MATH](https://core-math.gitlabpages.inria.fr/) kernels, widened one ulp
  per side — enclosures at most a couple of ulps wider than optimal.
- Correct handling of the hard cases: periodic extrema of `sin`/`cos`/`tan`
  across wrapping intervals, `atan2`'s branch cut, and domain restriction
  (e.g. `sqrt` of a partly-negative interval, `log` toward zero).
- Verified against a high-precision **mpmath oracle** across millions of
  inputs, and the rounding primitives cross-checked bit-for-bit against
  hardware directed rounding.

Also provides **decorated intervals** — a machine-checked certificate that a
composed function is defined and continuous on its input box. NumPy-style
interval arrays are planned.

## At a glance

```python
import pyintval as iv

x = iv.Interval("0.1")               # rigorously encloses the decimal 0.1
y = iv.sqrt(x) + iv.sin(x) * 2       # guaranteed enclosure of the true image
print(y.lo, y.hi)

iv.Interval(1) / iv.Interval(-1, 1)  # -> Interval('[entire]'), never raises
iv.exp(iv.log(iv.Interval(5)))       # encloses 5
math_pi_in = 3.141592653589793 in iv.pi()
```

## Installation

Once released:

```sh
pip install pyintval
```

From source (requires a C++20 compiler):

```sh
pip install .
```

On Windows, building from source requires **clang-cl** (the vendored CORE-MATH
kernels use features MSVC lacks); prebuilt wheels have no such requirement.

## Documentation

📖 **[pyintval.readthedocs.io](https://pyintval.readthedocs.io/)** — the user
guide and full API reference.

Built with Sphinx from `docs/`; build it locally with:

```sh
pip install ".[docs]"
sphinx-build -b html docs docs/_build/html
```

## Development

```sh
# Python tests
pip install -e ".[test]"
pytest

# C++ unit tests
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure

# Lint / format
pip install pre-commit
pre-commit run --all-files
```

See [RELEASING.md](RELEASING.md) for the wheel-building and PyPI publishing
process.

## Examples

The [`examples/`](examples/) directory has runnable scripts: basic usage,
rigorous range enclosure, an intermediate-value-theorem root proof, certified
transcendental bounds, and a continuity certificate via decorated intervals.

## License

MIT. Vendored components: [doctest](https://github.com/doctest/doctest) (MIT)
and the [CORE-MATH](https://core-math.gitlabpages.inria.fr/) correctly-rounded
kernels (MIT); see `third_party/`.
