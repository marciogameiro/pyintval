# pyintval

[![CI](https://github.com/marciogameiro/pyintval/actions/workflows/ci.yml/badge.svg)](https://github.com/marciogameiro/pyintval/actions/workflows/ci.yml)

**Rigorous interval arithmetic for Python** — every operation returns an
interval that is mathematically guaranteed to contain the true result.

> **Status: under active development.** The packaging skeleton is in place;
> the arithmetic kernel is being built milestone by milestone.

## What it will provide

- IEEE 1788-2015 **set-based interval arithmetic** over double-precision
  endpoints: empty and unbounded intervals propagate instead of raising, so
  bulk computations never abort mid-sweep.
- **Correctly rounded endpoints** for `+`, `-`, `*`, `/`, `sqrt` via
  error-free transformations (no rounding-mode switching, thread-safe), and
  transcendental functions (`exp`, `log`, `sin`, `cos`, `pow`, ...) built on
  the correctly rounded [CORE-MATH](https://core-math.gitlabpages.inria.fr/)
  routines, widened by 1 ulp per side.
- **Decorated intervals**: a machine-checked certificate that a composed
  function is defined and continuous on its input box — the hypothesis needed
  by many computer-assisted proofs.
- Full IEEE 1788 operation coverage plus extras (`cbrt`, `hypot`, `expm1`,
  `log1p`, `erf`, `erfc`), verified against the ITF1788 standard test suite
  and cross-checked against MPFR.

## Planned API (illustrative)

```python
import pyintval as iv

x = iv.Interval("0.1")          # rigorously encloses decimal 0.1
y = iv.sqrt(x) + iv.sin(x) * 2  # guaranteed enclosure of the true image
y.lo, y.hi, y in iv.Interval(0, 1)

d = iv.DecoratedInterval(-1, 1)
iv.sqrt(d).decoration            # 'trv': domain was violated — not continuous
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

## License

MIT. Vendored components: [doctest](https://github.com/doctest/doctest) (MIT);
CORE-MATH routines (MIT) will be vendored in a later milestone.
