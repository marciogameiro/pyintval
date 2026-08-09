#pragma once

// Interval elementary functions, built on the vendored CORE-MATH correctly
// rounded binary64 kernels (detail/coremath.hpp).
//
// Rigor principle. CORE-MATH's cr_f(x) is the double nearest to the true
// f(x), so |cr_f(x) - f(x)| <= 1/2 ulp(cr_f(x)). Therefore
//     pred(cr_f(x)) <= f(x) <= succ(cr_f(x))
// holds for every x (the half-ulp error can never reach the neighbor, and at a
// power-of-two boundary meets it exactly). We call pred/succ below "lb"/"ub":
// each transcendental endpoint is the correctly rounded value nudged one ulp
// outward, giving a rigorous enclosure at most ~1 ulp wider than optimal per
// endpoint.
//
// For each function we then apply its shape:
//   - monotone on its domain          -> map the (clipped) endpoints;
//   - even with an interior minimum    -> use mignitude/magnitude (cosh);
//   - periodic (sin, cos, tan)         -> map the endpoints and additionally
//     force the value to the extremum (+-1, or +-inf for tan) whenever the
//     interval might contain a critical point, decided by a rigorous, and
//     deliberately conservative, "contains an angle mod period" test;
//   - two-argument (atan2, hypot, pow) -> region/corner analysis.
//
// Domain violations shrink toward empty and unbounded limits widen toward the
// entire line, exactly as the algebraic kernel does.

#include <algorithm>
#include <cmath>

#include "pyintval/detail/coremath.hpp"
#include "pyintval/interval.hpp"

namespace pyintval {

namespace detail {

// One-ulp-outward directed bounds of a correctly rounded value.
inline double lb(double cr) noexcept { return pred(cr); }  // <= true value
inline double ub(double cr) noexcept { return succ(cr); }  // >= true value

// Rigorous enclosures of pi-related constants: the nearest double is below the
// true value, its successor above. These are constexpr aggregate literals --
// NOT computed via make()/add() -- so they carry no runtime static initializer.
// (A computed static initializer here ran fused kernel code at module load,
// which GCC 10 miscompiled into a crash when the extension was dlopen'd; a
// compile-time constant has nothing to miscompile.) kThreeHalfPi is the exact
// value add(kPi, kHalfPi) produces; kZero is [0, 0].
inline constexpr Interval kPi{0x1.921fb54442d18p1, 0x1.921fb54442d19p1};
inline constexpr Interval kHalfPi{0x1.921fb54442d18p0, 0x1.921fb54442d19p0};
inline constexpr Interval kTwoPi{0x1.921fb54442d18p2, 0x1.921fb54442d19p2};
inline constexpr Interval kThreeHalfPi{0x1.2d97c7f3321d2p2, 0x1.2d97c7f3321d3p2};
inline constexpr Interval kZero{0.0, 0.0};

// Conservative test: might the interval x contain a real number congruent to
// `target` modulo `period`? Computes n = (x - target) / period as an interval
// and asks whether its enclosure contains an integer. Because the enclosure is
// a superset of the true set of n-values, a "false" is certain (no such point
// exists) while a "true" may be a harmless over-report -- which only ever
// widens the resulting range, never breaks the enclosure.
inline bool maybe_contains_mod(const Interval& x, const Interval& target,
                               const Interval& period) noexcept {
  const Interval n = div(sub(x, target), period);
  if (is_empty(n)) return false;
  if (n.lo == -kInf || n.hi == kInf) return true;
  return std::ceil(n.lo) <= std::floor(n.hi);
}

}  // namespace detail

// --- Exponentials (monotone increasing, range depends on function) ----------

inline Interval exp(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  double lo = (x.lo == -detail::kInf) ? 0.0 : detail::lb(cr_exp(x.lo));
  const double hi = (x.hi == detail::kInf) ? detail::kInf : detail::ub(cr_exp(x.hi));
  if (lo < 0.0) lo = 0.0;  // exp > 0
  return make(lo, hi);
}

inline Interval exp2(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  double lo = (x.lo == -detail::kInf) ? 0.0 : detail::lb(cr_exp2(x.lo));
  const double hi = (x.hi == detail::kInf) ? detail::kInf : detail::ub(cr_exp2(x.hi));
  if (lo < 0.0) lo = 0.0;
  return make(lo, hi);
}

inline Interval exp10(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  double lo = (x.lo == -detail::kInf) ? 0.0 : detail::lb(cr_exp10(x.lo));
  const double hi = (x.hi == detail::kInf) ? detail::kInf : detail::ub(cr_exp10(x.hi));
  if (lo < 0.0) lo = 0.0;
  return make(lo, hi);
}

inline Interval expm1(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  double lo = (x.lo == -detail::kInf) ? -1.0 : detail::lb(cr_expm1(x.lo));
  const double hi = (x.hi == detail::kInf) ? detail::kInf : detail::ub(cr_expm1(x.hi));
  if (lo < -1.0) lo = -1.0;  // expm1 > -1
  return make(lo, hi);
}

// --- Logarithms (monotone increasing on (0, inf); (-1, inf) for log1p) ------

inline Interval log(const Interval& x) noexcept {
  if (is_empty(x) || x.hi <= 0.0) return empty();
  const double lo = (x.lo <= 0.0) ? -detail::kInf : detail::lb(cr_log(x.lo));
  const double hi = (x.hi == detail::kInf) ? detail::kInf : detail::ub(cr_log(x.hi));
  return make(lo, hi);
}

inline Interval log2(const Interval& x) noexcept {
  if (is_empty(x) || x.hi <= 0.0) return empty();
  const double lo = (x.lo <= 0.0) ? -detail::kInf : detail::lb(cr_log2(x.lo));
  const double hi = (x.hi == detail::kInf) ? detail::kInf : detail::ub(cr_log2(x.hi));
  return make(lo, hi);
}

inline Interval log10(const Interval& x) noexcept {
  if (is_empty(x) || x.hi <= 0.0) return empty();
  const double lo = (x.lo <= 0.0) ? -detail::kInf : detail::lb(cr_log10(x.lo));
  const double hi = (x.hi == detail::kInf) ? detail::kInf : detail::ub(cr_log10(x.hi));
  return make(lo, hi);
}

inline Interval log1p(const Interval& x) noexcept {
  if (is_empty(x) || x.hi <= -1.0) return empty();
  const double lo = (x.lo <= -1.0) ? -detail::kInf : detail::lb(cr_log1p(x.lo));
  const double hi = (x.hi == detail::kInf) ? detail::kInf : detail::ub(cr_log1p(x.hi));
  return make(lo, hi);
}

// --- Roots and monotone odd functions on all of R ---------------------------

inline Interval cbrt(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  const double lo = (x.lo == -detail::kInf) ? -detail::kInf : detail::lb(cr_cbrt(x.lo));
  const double hi = (x.hi == detail::kInf) ? detail::kInf : detail::ub(cr_cbrt(x.hi));
  return make(lo, hi);
}

inline Interval sinh(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  const double lo = (x.lo == -detail::kInf) ? -detail::kInf : detail::lb(cr_sinh(x.lo));
  const double hi = (x.hi == detail::kInf) ? detail::kInf : detail::ub(cr_sinh(x.hi));
  return make(lo, hi);
}

inline Interval asinh(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  const double lo = (x.lo == -detail::kInf) ? -detail::kInf : detail::lb(cr_asinh(x.lo));
  const double hi = (x.hi == detail::kInf) ? detail::kInf : detail::ub(cr_asinh(x.hi));
  return make(lo, hi);
}

inline Interval tanh(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  double lo = (x.lo == -detail::kInf) ? -1.0 : detail::lb(cr_tanh(x.lo));
  double hi = (x.hi == detail::kInf) ? 1.0 : detail::ub(cr_tanh(x.hi));
  if (lo < -1.0) lo = -1.0;
  if (hi > 1.0) hi = 1.0;  // tanh in (-1, 1)
  return make(lo, hi);
}

inline Interval erf(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  double lo = (x.lo == -detail::kInf) ? -1.0 : detail::lb(cr_erf(x.lo));
  double hi = (x.hi == detail::kInf) ? 1.0 : detail::ub(cr_erf(x.hi));
  if (lo < -1.0) lo = -1.0;
  if (hi > 1.0) hi = 1.0;  // erf in (-1, 1)
  return make(lo, hi);
}

// erfc is monotone DECREASING: erfc(-inf) = 2, erfc(+inf) = 0.
inline Interval erfc(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  double lo = (x.hi == detail::kInf) ? 0.0 : detail::lb(cr_erfc(x.hi));
  double hi = (x.lo == -detail::kInf) ? 2.0 : detail::ub(cr_erfc(x.lo));
  if (lo < 0.0) lo = 0.0;
  if (hi > 2.0) hi = 2.0;  // erfc in (0, 2)
  return make(lo, hi);
}

// --- Bounded monotone inverses (domain-clipped) -----------------------------

inline Interval atan(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  double lo = (x.lo == -detail::kInf) ? -detail::kHalfPi.hi : detail::lb(cr_atan(x.lo));
  double hi = (x.hi == detail::kInf) ? detail::kHalfPi.hi : detail::ub(cr_atan(x.hi));
  return make(lo, hi);
}

inline Interval asin(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  const double a = std::max(x.lo, -1.0), b = std::min(x.hi, 1.0);
  if (a > b) return empty();
  return make(detail::lb(cr_asin(a)), detail::ub(cr_asin(b)));
}

// acos is monotone DECREASING on [-1, 1].
inline Interval acos(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  const double a = std::max(x.lo, -1.0), b = std::min(x.hi, 1.0);
  if (a > b) return empty();
  return make(detail::lb(cr_acos(b)), detail::ub(cr_acos(a)));
}

inline Interval acosh(const Interval& x) noexcept {
  if (is_empty(x) || x.hi < 1.0) return empty();
  const double a = std::max(x.lo, 1.0);
  double lo = detail::lb(cr_acosh(a));
  const double hi = (x.hi == detail::kInf) ? detail::kInf : detail::ub(cr_acosh(x.hi));
  if (lo < 0.0) lo = 0.0;  // acosh >= 0
  return make(lo, hi);
}

// atanh on (-1, 1); +-inf as the argument approaches +-1.
inline Interval atanh(const Interval& x) noexcept {
  if (is_empty(x) || x.hi <= -1.0 || x.lo >= 1.0) return empty();
  const double lo = (x.lo <= -1.0) ? -detail::kInf : detail::lb(cr_atanh(x.lo));
  const double hi = (x.hi >= 1.0) ? detail::kInf : detail::ub(cr_atanh(x.hi));
  return make(lo, hi);
}

// --- cosh: even, minimum 1 at 0, increasing in |x| --------------------------

inline Interval cosh(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  const double near0 = mig(x);  // closest point to 0 (0 if the interval straddles 0)
  const double far = mag(x);    // farthest point from 0
  double lo = detail::lb(cr_cosh(near0));
  if (lo < 1.0) lo = 1.0;  // cosh >= 1
  const double hi = (far == detail::kInf) ? detail::kInf : detail::ub(cr_cosh(far));
  return make(lo, hi);
}

// --- Periodic: sin, cos, tan ------------------------------------------------

inline Interval sin(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  if (!is_common(x)) return make(-1.0, 1.0);  // unbounded: full oscillation
  double lo = std::min(detail::lb(cr_sin(x.lo)), detail::lb(cr_sin(x.hi)));
  double hi = std::max(detail::ub(cr_sin(x.lo)), detail::ub(cr_sin(x.hi)));
  // A point interval cannot contain an interior extremum, so skip the extremum
  // test (which, for large arguments, would loosely collapse the result to
  // [-1, 1] even though cr_sin evaluates the exact value).
  if (x.lo != x.hi) {
    if (detail::maybe_contains_mod(x, detail::kHalfPi, detail::kTwoPi)) hi = 1.0;
    if (detail::maybe_contains_mod(x, detail::kThreeHalfPi, detail::kTwoPi)) lo = -1.0;
  }
  if (lo < -1.0) lo = -1.0;
  if (hi > 1.0) hi = 1.0;
  return make(lo, hi);
}

inline Interval cos(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  if (!is_common(x)) return make(-1.0, 1.0);
  double lo = std::min(detail::lb(cr_cos(x.lo)), detail::lb(cr_cos(x.hi)));
  double hi = std::max(detail::ub(cr_cos(x.lo)), detail::ub(cr_cos(x.hi)));
  if (x.lo != x.hi) {
    if (detail::maybe_contains_mod(x, detail::kZero, detail::kTwoPi))
      hi = 1.0;  // cos=1 at 0 mod 2pi
    if (detail::maybe_contains_mod(x, detail::kPi, detail::kTwoPi))
      lo = -1.0;  // cos=-1 at pi mod 2pi
  }
  if (lo < -1.0) lo = -1.0;
  if (hi > 1.0) hi = 1.0;
  return make(lo, hi);
}

inline Interval tan(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  // Unbounded, or spanning at least a full period, or possibly crossing an
  // asymptote (pi/2 mod pi): the set-based image is the whole real line. A point
  // interval never straddles an asymptote (no double equals pi/2 + k*pi), so
  // skip the (large-argument-loose) asymptote test for singletons.
  if (!is_common(x)) return entire();
  if (x.lo != x.hi && detail::maybe_contains_mod(x, detail::kHalfPi, detail::kPi)) return entire();
  // No asymptote inside: tan is increasing on the single branch.
  return make(detail::lb(cr_tan(x.lo)), detail::ub(cr_tan(x.hi)));
}

// --- Two-argument functions -------------------------------------------------

// hypot(x, y) = sqrt(x^2 + y^2): increasing in |x| and |y|.
inline Interval hypot(const Interval& x, const Interval& y) noexcept {
  if (is_empty(x) || is_empty(y)) return empty();
  const double xn = mig(x), yn = mig(y);  // closest to the axes
  const double xf = mag(x), yf = mag(y);  // farthest
  double lo = detail::lb(cr_hypot(xn, yn));
  if (lo < 0.0) lo = 0.0;
  const double hi =
      (xf == detail::kInf || yf == detail::kInf) ? detail::kInf : detail::ub(cr_hypot(xf, yf));
  return make(lo, hi);
}

// atan2(y, x): the angle of (x, y) in (-pi, pi]. The branch cut lies along the
// negative x-axis (y = 0, x < 0). When the box avoids both the origin and that
// cut, atan2 is smooth and monotone in each variable separately, so its
// extrema over the box occur at the four corners. Otherwise the image can wrap
// through +-pi, and we return the conservative enclosure [-pi, pi].
inline Interval atan2(const Interval& y, const Interval& x) noexcept {
  using detail::kInf;
  using detail::kPi;
  if (is_empty(y) || is_empty(x)) return empty();
  const bool avoids_cut = (x.lo > 0.0) || (y.lo > 0.0) || (y.hi < 0.0);
  if (!avoids_cut) return make(-kPi.hi, kPi.hi);
  const double xs[2] = {x.lo, x.hi};
  const double ys[2] = {y.lo, y.hi};
  double lo = kInf, hi = -kInf;
  for (double yy : ys) {
    for (double xx : xs) {
      const double v = cr_atan2(yy, xx);
      lo = std::min(lo, v);
      hi = std::max(hi, v);
    }
  }
  return make(detail::lb(lo), detail::ub(hi));
}

// pow(x, y) = x^y. A point integer exponent reduces to pown (which handles
// zero and negative bases exactly). Otherwise a real power requires a
// nonnegative base, so the domain is restricted to x >= 0 and the value is the
// rigorous composition exp(y * log(x)).
inline Interval pow(const Interval& x, const Interval& y) noexcept {
  if (is_empty(x) || is_empty(y)) return empty();
  if (y.lo == y.hi && std::floor(y.lo) == y.lo && std::fabs(y.lo) <= 1.0e9) {
    return pown(x, static_cast<int>(y.lo));
  }
  const Interval xp = intersection(x, make(0.0, detail::kInf));
  if (is_empty(xp)) return empty();  // negative base with non-integer exponent: undefined
  return exp(mul(y, log(xp)));
}

}  // namespace pyintval
