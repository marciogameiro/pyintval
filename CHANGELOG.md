# Changelog

All notable changes to pyintval are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and the project follows
[Semantic Versioning](https://semver.org/).

## [0.1.0] — 2026-08-10

First public release.

### Added

- **Rigorous interval arithmetic** over IEEE 754 double-precision endpoints, in
  the IEEE 1788-2015 set-based flavor: empty and unbounded intervals propagate
  instead of raising, so bulk computations never abort mid-sweep.
- **Correctly rounded basic operations** — `+`, `-`, `*`, `/`, `sqrt`, `fma` —
  via error-free transformations with directed one-ulp rounding. No
  rounding-mode switching (thread-safe), and cross-checked bit-for-bit against
  hardware directed rounding over millions of inputs.
- **Correctly rounded elementary functions** — `exp`/`exp2`/`exp10`/`expm1`,
  `log`/`log2`/`log10`/`log1p`, `sin`/`cos`/`tan` and their inverses,
  hyperbolics and their inverses, `pow`, `atan2`, `hypot`, `cbrt`, `erf`/`erfc`
  — built on the MIT-licensed [CORE-MATH](https://core-math.gitlabpages.inria.fr/)
  kernels, widened one ulp per endpoint, with correct handling of periodic
  extrema, branch cuts, and domain restriction. Verified against a
  high-precision mpmath oracle.
- **Decorated intervals** (`DecoratedInterval`) implementing the IEEE 1788
  decoration system (`com`/`dac`/`def`/`trv`/`ill`) — a machine-checked
  certificate that a composed function is defined and continuous on its input.
- **Correctly rounded text I/O**: `Interval("0.1")` provably encloses the exact
  decimal 1/10, and formatting rounds outward.
- **Pythonic API**: operator overloading with scalar promotion, set operations,
  the full comparison-predicate set, rigorous `pi()`/`e()`, reverse
  (constraint-propagation) operations, `pickle`/`repr` round-tripping, and PEP
  561 type stubs. No runtime dependencies.
- Prebuilt wheels for Linux (x86-64, aarch64), macOS (arm64, x86-64), and
  Windows (AMD64) on CPython 3.10–3.14.

[0.1.0]: https://github.com/marciogameiro/pyintval/releases/tag/v0.1.0
