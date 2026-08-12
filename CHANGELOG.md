# Changelog

All notable changes to pyintval are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and the project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- **More IEEE 1788 text literals.** The uncertain form now accepts an infinite
  radius (`Interval("0.0??")`, `"2.5??u"`), and an overflowing radius or exponent
  saturates to an unbounded interval instead of failing. A decorated literal whose
  finite value overflows to unbounded downgrades an over-claimed `com` to `dac`
  (e.g. `DecoratedInterval("[1.0E+400]_com")`), while an explicitly unbounded
  literal stated as `com` remains NaI.
- **Rational text endpoints.** `Interval("2/3")`, `"[-4/2, 10/5]"` and the like
  parse `p/q` literals with correct outward rounding, decided by exact
  big-integer comparison, so the result rigorously encloses the exact rational.
  With this, pyintval **encloses the IEEE tightest result on all 7,236 ITF1788
  conformance tests** (transcendentals remain intentionally ~1 ulp wider than
  tightest).
- **ITF1788 conformance for the C++ kernel** (`tests/itf1788/run_conformance_cpp.py`):
  the reference corpus is now also applied to the header-only C++ API directly —
  no Python or pybind11 in the result path — via a doctest test-framework plugin
  (no Boost dependency) and a C++ arithmetic plugin. The kernel encloses the IEEE
  tightest interval on all 11,069 assertions (5,370 cases) of the 12 core corpus
  files. Runs as a CI job. See `tests/itf1788/README.md`.

### Fixed

- **Decorations of the step functions** `floor`/`ceil`/`trunc`/`round`/`sign`
  are now capped at `dac` and never `com`, matching IEEE 1788: even on a common
  input with a locally-constant (singleton) result, a step function is
  discontinuous in every neighborhood, so the previous `com` over-claimed
  continuity. Found by the ITF1788 conformance suite.
- **Decoration of `pow(x, y)` at a base touching zero.** When the base contains
  0 and the exponent reaches 0 (`0^0`) or below, the operation hits an undefined
  point on the box, so the result decoration is now `trv` instead of an
  over-claimed `com`/`dac`.

## [0.2.0] — 2026-08-10

### Added

- **Tightness regression tests** for composed interval expressions
  (`tests/{python,cpp}/test_tightness.py|cpp`): a fractions-based reference
  computes the tightest correctly-rounded result of a fixed evaluation form and
  asserts pyintval matches it bit-for-bit, over the dependency problem
  (`x*(1-x)` vs `x - x*x`), evaluation-order (Horner vs naive), Rump's example,
  and the wrapping effect.
- **IEEE 1788-2015 conformance testing via the ITF1788 reference suite**
  (`tests/itf1788/`): a plugin plus runner that generate and run the ~7,200-case
  corpus against pyintval and gate on the enclosure (rigor) standard. pyintval
  encloses the tightest result on 99.5% of tests exercising an implemented
  operation. See `tests/itf1788/README.md` for the documented remaining gaps.
- **Decorated overloads for the non-arithmetic operations**: `intersection`,
  `hull`, `cancel_minus`, `cancel_plus`, and the reverse operations `mul_rev`,
  `sqr_rev`, `abs_rev` now accept `DecoratedInterval` arguments (result carries
  the trivial `trv` decoration, and NaI propagates), as IEEE 1788 requires.
- **Expanded text parsing** toward the IEEE 1788 literal grammar: the uncertain
  form (`Interval("3.56?1")`, `"-10?u"`, `"3.56?1e2"`), half-bounded and
  empty-bound inf-sup literals (`"[1,]"`, `"[,2]"`, `"[,]"`), and native
  decorated literals (`DecoratedInterval("[1,2]_com")`, case-insensitive, with
  `"[nai]"` and over-claiming or unknown decorations resolving to NaI). Rational
  endpoints (`"[2/3, 1]"`) are not yet supported.

### Fixed

- **`pow(x, y)` at a base touching zero.** The two-argument power previously
  returned `[empty]` for a base interval containing 0 (it evaluated
  `exp(y·log x)` with `log 0 = -inf`). It now implements the IEEE 1788 corner
  semantics — `0^{y>0}=0`, `0^0` undefined only for a base of exactly `{0}`,
  limit values at the `x=0` boundary otherwise — so such results are correct
  rigorous enclosures. Found by the ITF1788 conformance suite.
- **Soundness: `cancel_minus`/`cancel_plus` near ±DBL_MAX on Linux.** GCC's `-O2`
  optimizer miscompiled the cancellative operations' overflow path in the
  header-only kernel, so the gcc-built (manylinux) wheels returned the empty set
  for `cancel_minus([max, max], [-max, -max])` — whose true result is
  `[max, +inf]` — instead of a valid enclosure. Pinned the two functions to `-O1`
  on GCC (clang, and hence the macOS/Windows wheels, were never affected). Also
  found by the ITF1788 conformance suite.

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

[0.2.0]: https://github.com/marciogameiro/pyintval/releases/tag/v0.2.0
[0.1.0]: https://github.com/marciogameiro/pyintval/releases/tag/v0.1.0
