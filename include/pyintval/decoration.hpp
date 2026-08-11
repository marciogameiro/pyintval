#pragma once

// Decorated intervals (IEEE 1788-2015 §8): a bare interval paired with a
// decoration that certifies, purely from the computation's history, the
// strongest property known to hold of the function evaluated so far.
//
//   com  ("common")   : defined, continuous, and bounded on the input, and the
//                       input was a common (bounded, nonempty) interval.
//   dac  ("dac")      : defined and continuous on the input.
//   def  ("defined")  : defined on the input, but continuity is not guaranteed
//                       (e.g. floor, sign crossing a step).
//   trv  ("trivial")  : only the enclosure is guaranteed; a domain violation
//                       may have occurred somewhere on the input.
//   ill  ("ill")      : the decorated interval is not well-formed (NaI); this
//                       poisons every computation it touches.
//
// Propagation rule (§8.8.4): for an operation phi with bare inputs X_i carrying
// decorations d_i, the output decoration is
//        min( d_local(phi, X_1..X_n),  min_i d_i )
// where d_local is the strongest property phi itself guarantees on the bare
// input box. Because the decorations are totally ordered (ill < trv < def <
// dac < com) and we always take the minimum, a decoration is a rigorous LOWER
// bound on the truth: it never over-claims. Where a test must be conservative
// (e.g. "does this interval contain a tan asymptote?"), reporting the weaker
// decoration stays sound.
//
// After evaluating an arbitrary composed formula F on a box X, a result
// decoration of dac or com is a machine-checked certificate that F is defined
// and continuous on X -- exactly the hypothesis many computer-assisted proofs
// require -- with no per-function hand analysis.

#include <algorithm>

#include "pyintval/elementary.hpp"
#include "pyintval/interval.hpp"

namespace pyintval {

enum class Decoration : int { ill = 0, trv = 1, def = 2, dac = 3, com = 4 };

inline Decoration dmin(Decoration a, Decoration b) noexcept {
  return static_cast<int>(a) < static_cast<int>(b) ? a : b;
}

inline const char* decoration_name(Decoration d) noexcept {
  switch (d) {
    case Decoration::com:
      return "com";
    case Decoration::dac:
      return "dac";
    case Decoration::def:
      return "def";
    case Decoration::trv:
      return "trv";
    default:
      return "ill";
  }
}

struct DecoratedInterval {
  Interval x;
  Decoration dec;
};

namespace detail {

// The strongest decoration a bare interval can carry when freshly constructed.
inline Decoration newdec(const Interval& x) noexcept {
  if (is_empty(x)) return Decoration::trv;
  if (is_common(x)) return Decoration::com;
  return Decoration::dac;  // nonempty but unbounded
}

// Local decoration of a defined-and-continuous operation: com requires every
// input to be common AND the result to be common (bounded, nonempty). An
// unbounded result -- e.g. exp near overflow, or cosh([709.8, 710]) reaching
// +inf -- is only dac even from common inputs. Empty result (only from an empty
// input) degrades to trv.
inline Decoration cont_local(std::initializer_list<Interval> inputs, const Interval& y) noexcept {
  if (is_empty(y)) return Decoration::trv;
  for (const Interval& in : inputs) {
    if (!is_common(in)) return Decoration::dac;
  }
  return is_common(y) ? Decoration::com : Decoration::dac;
}

}  // namespace detail

// --- Construction -----------------------------------------------------------

inline DecoratedInterval decorate(const Interval& x) noexcept { return {x, detail::newdec(x)}; }
inline DecoratedInterval decorate(const Interval& x, Decoration d) noexcept { return {x, d}; }

// A Not-an-Interval marker: ill-formed, poisons every computation.
inline DecoratedInterval nai() noexcept { return {empty(), Decoration::ill}; }

// --- Helpers to assemble a decorated result ---------------------------------

namespace detail {
inline DecoratedInterval finish(const Interval& y, Decoration local, Decoration d0) noexcept {
  return {y, dmin(local, d0)};
}
inline DecoratedInterval finish(const Interval& y, Decoration local, Decoration d0,
                                Decoration d1) noexcept {
  return {y, dmin(local, dmin(d0, d1))};
}
inline DecoratedInterval finish(const Interval& y, Decoration local, Decoration d0, Decoration d1,
                                Decoration d2) noexcept {
  return {y, dmin(local, dmin(d0, dmin(d1, d2)))};
}
}  // namespace detail

// --- Arithmetic (defined and continuous everywhere, except division) --------

inline DecoratedInterval neg(const DecoratedInterval& a) {
  const Interval y = neg(a.x);
  return detail::finish(y, detail::cont_local({a.x}, y), a.dec);
}
inline DecoratedInterval add(const DecoratedInterval& a, const DecoratedInterval& b) {
  const Interval y = add(a.x, b.x);
  return detail::finish(y, detail::cont_local({a.x, b.x}, y), a.dec, b.dec);
}
inline DecoratedInterval sub(const DecoratedInterval& a, const DecoratedInterval& b) {
  const Interval y = sub(a.x, b.x);
  return detail::finish(y, detail::cont_local({a.x, b.x}, y), a.dec, b.dec);
}
inline DecoratedInterval mul(const DecoratedInterval& a, const DecoratedInterval& b) {
  const Interval y = mul(a.x, b.x);
  return detail::finish(y, detail::cont_local({a.x, b.x}, y), a.dec, b.dec);
}
inline DecoratedInterval fma(const DecoratedInterval& a, const DecoratedInterval& b,
                             const DecoratedInterval& c) {
  const Interval y = fma(a.x, b.x, c.x);
  return detail::finish(y, detail::cont_local({a.x, b.x, c.x}, y), a.dec, b.dec, c.dec);
}
inline DecoratedInterval sqr(const DecoratedInterval& a) {
  const Interval y = sqr(a.x);
  return detail::finish(y, detail::cont_local({a.x}, y), a.dec);
}
inline DecoratedInterval abs(const DecoratedInterval& a) {
  const Interval y = abs(a.x);
  return detail::finish(y, detail::cont_local({a.x}, y), a.dec);
}
inline DecoratedInterval min(const DecoratedInterval& a, const DecoratedInterval& b) {
  const Interval y = min(a.x, b.x);
  return detail::finish(y, detail::cont_local({a.x, b.x}, y), a.dec, b.dec);
}
inline DecoratedInterval max(const DecoratedInterval& a, const DecoratedInterval& b) {
  const Interval y = max(a.x, b.x);
  return detail::finish(y, detail::cont_local({a.x, b.x}, y), a.dec, b.dec);
}

// Division is undefined where the divisor is zero.
inline DecoratedInterval div(const DecoratedInterval& a, const DecoratedInterval& b) {
  const Interval y = div(a.x, b.x);
  Decoration local = contains_zero(b.x) ? Decoration::trv : detail::cont_local({a.x, b.x}, y);
  return detail::finish(y, local, a.dec, b.dec);
}
inline DecoratedInterval recip(const DecoratedInterval& a) {
  const Interval y = recip(a.x);
  Decoration local = contains_zero(a.x) ? Decoration::trv : detail::cont_local({a.x}, y);
  return detail::finish(y, local, a.dec);
}

// pown: defined and continuous everywhere for n >= 0; for n < 0 undefined at 0.
inline DecoratedInterval pown(const DecoratedInterval& a, int n) {
  const Interval y = pown(a.x, n);
  Decoration local = (n < 0 && contains_zero(a.x)) ? Decoration::trv : detail::cont_local({a.x}, y);
  return detail::finish(y, local, a.dec);
}

// --- Discontinuous but everywhere-defined step functions --------------------

namespace detail {
// A step function (floor/ceil/trunc/round/sign) is locally constant -- hence
// continuous -- exactly when its bare result collapses to a single value, giving
// dac; otherwise it jumped across a discontinuity, giving def. It is never com:
// even on a common input a step function is discontinuous in every neighborhood,
// so IEEE 1788 caps its local decoration at dac.
inline Decoration step_local(const Interval& y) noexcept {
  if (is_empty(y)) return Decoration::trv;
  if (is_singleton(y)) return Decoration::dac;
  return Decoration::def;
}
}  // namespace detail

inline DecoratedInterval floor(const DecoratedInterval& a) {
  const Interval y = floor(a.x);
  return detail::finish(y, detail::step_local(y), a.dec);
}
inline DecoratedInterval ceil(const DecoratedInterval& a) {
  const Interval y = ceil(a.x);
  return detail::finish(y, detail::step_local(y), a.dec);
}
inline DecoratedInterval trunc(const DecoratedInterval& a) {
  const Interval y = trunc(a.x);
  return detail::finish(y, detail::step_local(y), a.dec);
}
inline DecoratedInterval round_ties_to_even(const DecoratedInterval& a) {
  const Interval y = round_ties_to_even(a.x);
  return detail::finish(y, detail::step_local(y), a.dec);
}
inline DecoratedInterval round_ties_to_away(const DecoratedInterval& a) {
  const Interval y = round_ties_to_away(a.x);
  return detail::finish(y, detail::step_local(y), a.dec);
}
inline DecoratedInterval sign(const DecoratedInterval& a) {
  const Interval y = sign(a.x);
  return detail::finish(y, detail::step_local(y), a.dec);
}

// --- Domain-restricted elementary functions ---------------------------------

namespace detail {
// Local decoration for a function defined & continuous on a domain: trv if the
// input escapes the domain, else com/dac by commonness. `violated` is true iff
// some point of the input lies outside the (continuity) domain.
inline Decoration domain_local(bool violated, const Interval& x, const Interval& y) noexcept {
  if (violated) return Decoration::trv;
  return cont_local({x}, y);
}
}  // namespace detail

inline DecoratedInterval sqrt(const DecoratedInterval& a) {
  const Interval y = sqrt(a.x);
  const bool violated = is_empty(a.x) ? false : (a.x.lo < 0.0);
  return detail::finish(y, detail::domain_local(violated, a.x, y), a.dec);
}

// Continuous everywhere on R.
#define PYINTVAL_DEC_CONT_UNARY(NAME)                              \
  inline DecoratedInterval NAME(const DecoratedInterval& a) {      \
    const Interval y = NAME(a.x);                                  \
    return detail::finish(y, detail::cont_local({a.x}, y), a.dec); \
  }
PYINTVAL_DEC_CONT_UNARY(exp)
PYINTVAL_DEC_CONT_UNARY(exp2)
PYINTVAL_DEC_CONT_UNARY(exp10)
PYINTVAL_DEC_CONT_UNARY(expm1)
PYINTVAL_DEC_CONT_UNARY(cbrt)
PYINTVAL_DEC_CONT_UNARY(sin)
PYINTVAL_DEC_CONT_UNARY(cos)
PYINTVAL_DEC_CONT_UNARY(atan)
PYINTVAL_DEC_CONT_UNARY(sinh)
PYINTVAL_DEC_CONT_UNARY(cosh)
PYINTVAL_DEC_CONT_UNARY(tanh)
PYINTVAL_DEC_CONT_UNARY(asinh)
PYINTVAL_DEC_CONT_UNARY(erf)
PYINTVAL_DEC_CONT_UNARY(erfc)
#undef PYINTVAL_DEC_CONT_UNARY

// Defined & continuous on (0, inf).
#define PYINTVAL_DEC_LOG(NAME)                                               \
  inline DecoratedInterval NAME(const DecoratedInterval& a) {                \
    const Interval y = NAME(a.x);                                            \
    const bool violated = is_empty(a.x) ? false : (a.x.lo <= 0.0);           \
    return detail::finish(y, detail::domain_local(violated, a.x, y), a.dec); \
  }
PYINTVAL_DEC_LOG(log)
PYINTVAL_DEC_LOG(log2)
PYINTVAL_DEC_LOG(log10)
#undef PYINTVAL_DEC_LOG

inline DecoratedInterval log1p(const DecoratedInterval& a) {
  const Interval y = log1p(a.x);
  const bool violated = is_empty(a.x) ? false : (a.x.lo <= -1.0);
  return detail::finish(y, detail::domain_local(violated, a.x, y), a.dec);
}

// Defined & continuous on [-1, 1] (closed; the endpoints are in the domain).
#define PYINTVAL_DEC_ASINCOS(NAME)                                                 \
  inline DecoratedInterval NAME(const DecoratedInterval& a) {                      \
    const Interval y = NAME(a.x);                                                  \
    const bool violated = is_empty(a.x) ? false : (a.x.lo < -1.0 || a.x.hi > 1.0); \
    return detail::finish(y, detail::domain_local(violated, a.x, y), a.dec);       \
  }
PYINTVAL_DEC_ASINCOS(asin)
PYINTVAL_DEC_ASINCOS(acos)
#undef PYINTVAL_DEC_ASINCOS

inline DecoratedInterval acosh(const DecoratedInterval& a) {
  const Interval y = acosh(a.x);
  const bool violated = is_empty(a.x) ? false : (a.x.lo < 1.0);
  return detail::finish(y, detail::domain_local(violated, a.x, y), a.dec);
}

inline DecoratedInterval atanh(const DecoratedInterval& a) {
  const Interval y = atanh(a.x);
  const bool violated = is_empty(a.x) ? false : (a.x.lo <= -1.0 || a.x.hi >= 1.0);
  return detail::finish(y, detail::domain_local(violated, a.x, y), a.dec);
}

// tan: undefined at the asymptotes pi/2 + k*pi. Use the same conservative
// angle-mod-period test as the bare version (over-reporting an asymptote only
// weakens the decoration to trv, which is sound).
inline DecoratedInterval tan(const DecoratedInterval& a) {
  const Interval y = tan(a.x);
  bool violated = true;
  if (!is_empty(a.x) && is_common(a.x)) {
    violated = detail::maybe_contains_mod(a.x, detail::kHalfPi, detail::kPi);
  }
  return detail::finish(y, detail::domain_local(violated, a.x, y), a.dec);
}

// hypot: defined & continuous on all of R^2.
inline DecoratedInterval hypot(const DecoratedInterval& a, const DecoratedInterval& b) {
  const Interval y = hypot(a.x, b.x);
  return detail::finish(y, detail::cont_local({a.x, b.x}, y), a.dec, b.dec);
}

// atan2(y, x): undefined at the origin; defined but discontinuous across the
// negative x-axis branch cut; continuous elsewhere.
inline DecoratedInterval atan2(const DecoratedInterval& y, const DecoratedInterval& x) {
  const Interval r = atan2(y.x, x.x);
  Decoration local;
  if (is_empty(y.x) || is_empty(x.x)) {
    local = Decoration::trv;
  } else {
    const bool origin = x.x.lo <= 0.0 && 0.0 <= x.x.hi && y.x.lo <= 0.0 && 0.0 <= y.x.hi;
    const bool avoids_cut = (x.x.lo > 0.0) || (y.x.lo > 0.0) || (y.x.hi < 0.0);
    if (origin) {
      local = Decoration::trv;  // undefined at (0,0)
    } else if (avoids_cut) {
      local = detail::cont_local({y.x, x.x}, r);
    } else {
      local = Decoration::def;  // crosses the branch cut: defined but not continuous
    }
  }
  return detail::finish(r, local, y.dec, x.dec);
}

// pow(x, y): real power requires x > 0 to be defined & continuous. A point
// integer exponent behaves like pown (defined for all x, minus 0 for n < 0).
inline DecoratedInterval pow(const DecoratedInterval& a, const DecoratedInterval& b) {
  const Interval y = pow(a.x, b.x);
  Decoration local;
  if (is_empty(a.x) || is_empty(b.x)) {
    local = Decoration::trv;
  } else if (b.x.lo == b.x.hi && std::floor(b.x.lo) == b.x.lo && std::fabs(b.x.lo) <= 1.0e9) {
    const bool neg_exp_at_zero = b.x.lo < 0.0 && contains_zero(a.x);
    local = neg_exp_at_zero ? Decoration::trv : detail::cont_local({a.x}, y);
  } else {
    local = (a.x.lo > 0.0) ? detail::cont_local({a.x, b.x}, y) : Decoration::trv;
  }
  return detail::finish(y, local, a.dec, b.dec);
}

// --- Set operations, cancellative subtraction, reverse operations -----------
//
// These are not point functions, so IEEE 1788 (§11.11) gives their result the
// trivial decoration trv regardless of the inputs' local behaviour. Passing trv
// as the local decoration through finish() also propagates NaI: an ill input
// (decoration 0) is the minimum, so it poisons the result exactly as required.

inline DecoratedInterval intersection(const DecoratedInterval& a, const DecoratedInterval& b) {
  return detail::finish(intersection(a.x, b.x), Decoration::trv, a.dec, b.dec);
}
inline DecoratedInterval convex_hull(const DecoratedInterval& a, const DecoratedInterval& b) {
  return detail::finish(convex_hull(a.x, b.x), Decoration::trv, a.dec, b.dec);
}
inline DecoratedInterval cancel_minus(const DecoratedInterval& a, const DecoratedInterval& b) {
  return detail::finish(cancel_minus(a.x, b.x), Decoration::trv, a.dec, b.dec);
}
inline DecoratedInterval cancel_plus(const DecoratedInterval& a, const DecoratedInterval& b) {
  return detail::finish(cancel_plus(a.x, b.x), Decoration::trv, a.dec, b.dec);
}
inline DecoratedInterval mul_rev(const DecoratedInterval& b, const DecoratedInterval& c) {
  return detail::finish(mul_rev(b.x, c.x), Decoration::trv, b.dec, c.dec);
}
inline DecoratedInterval mul_rev(const DecoratedInterval& b, const DecoratedInterval& c,
                                 const DecoratedInterval& x) {
  return detail::finish(mul_rev(b.x, c.x, x.x), Decoration::trv, b.dec, c.dec, x.dec);
}
inline DecoratedInterval sqr_rev(const DecoratedInterval& c) {
  return detail::finish(sqr_rev(c.x), Decoration::trv, c.dec);
}
inline DecoratedInterval sqr_rev(const DecoratedInterval& c, const DecoratedInterval& x) {
  return detail::finish(sqr_rev(c.x, x.x), Decoration::trv, c.dec, x.dec);
}
inline DecoratedInterval abs_rev(const DecoratedInterval& c) {
  return detail::finish(abs_rev(c.x), Decoration::trv, c.dec);
}
inline DecoratedInterval abs_rev(const DecoratedInterval& c, const DecoratedInterval& x) {
  return detail::finish(abs_rev(c.x, x.x), Decoration::trv, c.dec, x.dec);
}

}  // namespace pyintval
