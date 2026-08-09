# Vendored CORE-MATH routines

The subdirectories here contain correctly-rounded binary64 elementary-function
implementations from the **CORE-MATH** project, used by pyintval as the
round-to-nearest kernels that its interval wrappers widen outward by one ulp to
obtain rigorous enclosures.

- Upstream: https://core-math.gitlabpages.inria.fr/
- Repository: https://gitlab.inria.fr/core-math/core-math
- Pinned commit: `07cf01e12a42b82cc478341982936cad7f3f9bdc` (branch `master`)
- Retrieved: 2026-08-09
- License: MIT (see `LICENSE` in this directory)

Each function directory (`exp/`, `log/`, `sin/`, ...) holds the upstream source
file `<fn>.c` plus the private headers it `#include`s (e.g. `dint.h`, `qint.h`,
`pow.h`, `tint.h`). The `#include "..."` directives resolve relative to each
source file, so the directories are self-contained and compile with no extra
`-I` flags. Note that the several `dint.h` copies are intentionally *not*
identical across functions — do not deduplicate them onto a shared include path.

Only the `cr_<fn>` production entry points are vendored; the upstream test
harnesses (`check_worst.c`, Makefiles) are not included. The sources are
compiled unmodified. They require a compiler with `__builtin_*`, `__int128`,
and (on x86-64) `<x86intrin.h>` support — i.e. Clang or GCC on every platform,
and clang-cl (not MSVC) on Windows; the upstream `#if defined(__x86_64__)`
guards supply portable fallbacks on AArch64.

CORE-MATH's correct-rounding guarantee assumes the FE_TONEAREST rounding mode,
which pyintval never leaves (it emulates directed rounding without touching the
FPU mode), so the assumption always holds.
