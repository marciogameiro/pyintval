#pragma once

// pyintval kernel — header-only C++ core.
//
// Milestone 1 placeholder: establishes the namespace, the ABI version, and the
// platform requirements that the entire rigor argument of this library rests
// on. The rigorous interval type and its directed-rounding arithmetic land in
// Milestone 2.

#include <limits>

namespace pyintval {

// Structural requirement: every enclosure proof in this library assumes that
// `double` is an IEEE 754 (IEC 559) binary64 type — 53-bit significand, with
// infinities, signed zeros, and subnormals. Refuse to compile anywhere this
// does not hold rather than silently produce non-rigorous results.
static_assert(std::numeric_limits<double>::is_iec559,
              "pyintval requires IEEE 754 (IEC 559) double-precision floats");
static_assert(std::numeric_limits<double>::digits == 53,
              "pyintval requires binary64 doubles (53-bit significand)");
static_assert(std::numeric_limits<double>::has_infinity, "pyintval requires IEEE 754 infinities");

// Bumped whenever the in-memory layout of the kernel types changes.
inline constexpr int kernel_abi_version = 0;

}  // namespace pyintval
