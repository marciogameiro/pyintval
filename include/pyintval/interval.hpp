#pragma once

// pyintval kernel: IEEE 1788-2015 set-based bare intervals over binary64.
//
// An interval is a closed, connected subset of the reals: bounded [a,b],
// half-lines [a,+inf) / (-inf,b], the whole line ("entire"), or the empty
// set. Every operation f returns an enclosure of the exact image
//   hull{ f(x) : x in X restricted to dom(f) }
// with endpoints rounded outward (lower endpoints toward -inf, upper toward
// +inf) using the directed-rounding primitives in detail/rounding.hpp.
// Nothing here raises or returns NaN intervals mid-computation: domain
// violations shrink toward the empty set, unboundedness widens toward
// entire, and both propagate silently -- decorations (Milestone 4) recover
// the "was the domain respected" information when the caller needs it.
//
// Representation invariants (enforced by the factories):
//   nonempty:  lo <= hi, lo < +inf, hi > -inf; zero endpoints stored as +0.0
//   empty:     canonically {+inf, -inf}  (so inf(empty) = +inf and
//              sup(empty) = -inf fall out for free, as 1788 9.2 specifies)
// The struct is a plain aggregate and every operation is a free function on
// values: the same kernels will later back an IntervalArray without change.

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <ostream>
#include <utility>

#include "pyintval/detail/rounding.hpp"

namespace pyintval {

// The entire rigor argument assumes IEEE 754 binary64 semantics; refuse to
// compile anywhere this does not hold.
static_assert(std::numeric_limits<double>::is_iec559,
              "pyintval requires IEEE 754 (IEC 559) double-precision floats");
static_assert(std::numeric_limits<double>::digits == 53,
              "pyintval requires binary64 doubles (53-bit significand)");

// Bumped whenever the in-memory layout of the kernel types changes.
inline constexpr int kernel_abi_version = 1;

struct Interval {
  double lo;
  double hi;
};

namespace detail {
// Zero endpoints are stored as +0.0 so that bitwise-identical canonical forms
// exist ([-0.0, 0.0] and [0.0, 0.0] denote the same set); comparisons below
// use values, but tests and hashing want one representation.
inline double canon_zero(double v) noexcept { return v == 0.0 ? 0.0 : v; }
}  // namespace detail

// --- Factories --------------------------------------------------------------

inline Interval empty() noexcept { return {detail::kInf, -detail::kInf}; }
inline Interval entire() noexcept { return {-detail::kInf, detail::kInf}; }

// Valid bounds: lo <= hi, excluding the point-at-infinity cases [+inf,+inf]
// and [-inf,-inf] (infinities are allowed only as *open* directions, never as
// members). NaN fails every comparison and is therefore invalid.
inline bool valid_bounds(double lo, double hi) noexcept {
  return lo <= hi && lo < detail::kInf && hi > -detail::kInf;
}

// Construct from bounds; invalid bounds yield the empty set (the Python layer
// raises instead, after calling valid_bounds itself).
inline Interval make(double lo, double hi) noexcept {
  if (!valid_bounds(lo, hi)) return empty();
  return {detail::canon_zero(lo), detail::canon_zero(hi)};
}

inline Interval point(double v) noexcept { return make(v, v); }

// --- Basic queries ----------------------------------------------------------

inline bool is_empty(const Interval& x) noexcept { return !(x.lo <= x.hi); }
inline double inf(const Interval& x) noexcept { return x.lo; }
inline double sup(const Interval& x) noexcept { return x.hi; }
inline bool is_entire(const Interval& x) noexcept {
  return x.lo == -detail::kInf && x.hi == detail::kInf;
}
// 1788 "common" interval: nonempty and bounded.
inline bool is_common(const Interval& x) noexcept {
  return !is_empty(x) && x.lo > -detail::kInf && x.hi < detail::kInf;
}
inline bool is_singleton(const Interval& x) noexcept { return !is_empty(x) && x.lo == x.hi; }
inline bool is_member(double v, const Interval& x) noexcept {
  return std::isfinite(v) && !is_empty(x) && x.lo <= v && v <= x.hi;
}
inline bool contains_zero(const Interval& x) noexcept {
  return !is_empty(x) && x.lo <= 0.0 && 0.0 <= x.hi;
}

// --- Set operations ---------------------------------------------------------

inline Interval intersection(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a) || is_empty(b)) return empty();
  return make(std::max(a.lo, b.lo), std::min(a.hi, b.hi));  // invalid -> empty
}

// Interval hull of the union (the union itself need not be an interval).
inline Interval convex_hull(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a)) return b;
  if (is_empty(b)) return a;
  return make(std::min(a.lo, b.lo), std::max(a.hi, b.hi));
}

// --- Numeric functions (1788 10.5.9) ----------------------------------------

// Midpoint: a finite double guaranteed to lie inside the interval.
// Conventions: mid(empty) = NaN, mid(entire) = 0, half-lines map to the
// finite double of largest magnitude on their unbounded side.
inline double mid(const Interval& x) noexcept {
  if (is_empty(x)) return std::numeric_limits<double>::quiet_NaN();
  if (is_entire(x)) return 0.0;
  if (x.lo == -detail::kInf) return -detail::kMaxDouble;
  if (x.hi == detail::kInf) return detail::kMaxDouble;
  double m = 0.5 * (x.lo + x.hi);
  if (!std::isfinite(m)) m = 0.5 * x.lo + 0.5 * x.hi;  // lo+hi overflowed; halves are exact
  // Clamp: RN can push the naive midpoint of subnormal endpoints just outside.
  return std::min(std::max(m, x.lo), x.hi);
}

// Radius: rad(x) is the smallest double with [mid - rad, mid + rad] >= x.
inline double rad(const Interval& x) noexcept {
  if (is_empty(x)) return std::numeric_limits<double>::quiet_NaN();
  const double m = mid(x);
  return std::max(detail::sub_ru(m, x.lo), detail::sub_ru(x.hi, m));
}

inline double wid(const Interval& x) noexcept {
  if (is_empty(x)) return std::numeric_limits<double>::quiet_NaN();
  return detail::sub_ru(x.hi, x.lo);
}

// Magnitude sup|x| and mignitude min|x| over the interval.
inline double mag(const Interval& x) noexcept {
  if (is_empty(x)) return std::numeric_limits<double>::quiet_NaN();
  return std::max(std::fabs(x.lo), std::fabs(x.hi));
}

inline double mig(const Interval& x) noexcept {
  if (is_empty(x)) return std::numeric_limits<double>::quiet_NaN();
  if (contains_zero(x)) return 0.0;
  return std::min(std::fabs(x.lo), std::fabs(x.hi));
}

// --- Comparison predicates (1788 10.5.10, with the empty-set conventions) ---

inline bool equal(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a) || is_empty(b)) return is_empty(a) && is_empty(b);
  return a.lo == b.lo && a.hi == b.hi;
}

inline bool subset(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a)) return true;
  if (is_empty(b)) return false;
  return b.lo <= a.lo && a.hi <= b.hi;
}

// a is in the topological interior of b (relative to the extended reals:
// matching infinities are interior on their unbounded side).
inline bool interior(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a)) return true;
  if (is_empty(b)) return false;
  const bool lo_ok = b.lo < a.lo || (b.lo == -detail::kInf && a.lo == -detail::kInf);
  const bool hi_ok = a.hi < b.hi || (a.hi == detail::kInf && b.hi == detail::kInf);
  return lo_ok && hi_ok;
}

inline bool disjoint(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a) || is_empty(b)) return true;
  return a.hi < b.lo || b.hi < a.lo;
}

// Weak order "a is (weakly) to the left of b": both bounds no larger.
inline bool less(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a) || is_empty(b)) return is_empty(a) && is_empty(b);
  return a.lo <= b.lo && a.hi <= b.hi;
}

inline bool strict_less(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a) || is_empty(b)) return is_empty(a) && is_empty(b);
  const bool lo_ok = a.lo < b.lo || (a.lo == -detail::kInf && b.lo == -detail::kInf);
  const bool hi_ok = a.hi < b.hi || (a.hi == detail::kInf && b.hi == detail::kInf);
  return lo_ok && hi_ok;
}

inline bool precedes(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a) || is_empty(b)) return true;
  return a.hi <= b.lo;
}

inline bool strict_precedes(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a) || is_empty(b)) return true;
  return a.hi < b.lo;
}

// --- Arithmetic -------------------------------------------------------------

inline Interval neg(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  return make(-x.hi, -x.lo);
}

// Image of {u + v}: endpoint-wise monotone, so lower+lower and upper+upper.
// -inf + finite and +inf + finite are exact in the primitives; -inf and +inf
// can never meet (lo < +inf and hi > -inf on nonempty operands).
inline Interval add(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a) || is_empty(b)) return empty();
  return make(detail::add_rd(a.lo, b.lo), detail::add_ru(a.hi, b.hi));
}

inline Interval sub(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a) || is_empty(b)) return empty();
  return make(detail::sub_rd(a.lo, b.hi), detail::sub_ru(a.hi, b.lo));
}

namespace detail {
// Sign class of a nonempty interval, used to pick the extremal endpoint
// pairing for products/quotients without ever forming 0 * inf (which is NaN
// in FP but contributes 0 to the exact image over interval members).
enum class SignClass { kZero, kPos, kNeg, kMix };  // [0,0] | >=0 | <=0 | straddles

inline SignClass classify(const Interval& x) noexcept {
  if (x.lo == 0.0 && x.hi == 0.0) return SignClass::kZero;
  if (x.lo >= 0.0) return SignClass::kPos;
  if (x.hi <= 0.0) return SignClass::kNeg;
  return SignClass::kMix;
}
}  // namespace detail

// Image of {u * v}. Case analysis over sign classes; in every case the chosen
// endpoint pairs are the extremal ones and never pair a zero endpoint with an
// infinite one (kPos pairs lo-lo/hi-hi, kNeg pairs are both-negative, etc.).
inline Interval mul(const Interval& a, const Interval& b) noexcept {
  using detail::mul_rd;
  using detail::mul_ru;
  using detail::SignClass;
  if (is_empty(a) || is_empty(b)) return empty();
  const SignClass ca = detail::classify(a);
  const SignClass cb = detail::classify(b);
  if (ca == SignClass::kZero || cb == SignClass::kZero) return point(0.0);

  switch (ca) {
    case SignClass::kPos:
      switch (cb) {
        case SignClass::kPos:
          return make(mul_rd(a.lo, b.lo), mul_ru(a.hi, b.hi));
        case SignClass::kNeg:
          return make(mul_rd(a.hi, b.lo), mul_ru(a.lo, b.hi));
        default:  // kMix
          return make(mul_rd(a.hi, b.lo), mul_ru(a.hi, b.hi));
      }
    case SignClass::kNeg:
      switch (cb) {
        case SignClass::kPos:
          return make(mul_rd(a.lo, b.hi), mul_ru(a.hi, b.lo));
        case SignClass::kNeg:
          return make(mul_rd(a.hi, b.hi), mul_ru(a.lo, b.lo));
        default:  // kMix
          return make(mul_rd(a.lo, b.hi), mul_ru(a.lo, b.lo));
      }
    default:  // kMix
      switch (cb) {
        case SignClass::kPos:
          return make(mul_rd(a.lo, b.hi), mul_ru(a.hi, b.hi));
        case SignClass::kNeg:
          return make(mul_rd(a.hi, b.lo), mul_ru(a.lo, b.lo));
        default:  // kMix: both cross-pairs compete for each bound
          return make(std::min(mul_rd(a.lo, b.hi), mul_rd(a.hi, b.lo)),
                      std::max(mul_ru(a.lo, b.lo), mul_ru(a.hi, b.hi)));
      }
  }
}

// Image of {u / v} over v != 0 (1788 set-based division). When 0 is interior
// to b the image has two unbounded components and the hull is entire; when 0
// is an endpoint of b only one direction blows up; division by [0,0] has
// empty domain.
inline Interval div(const Interval& a, const Interval& b) noexcept {
  using detail::div_rd;
  using detail::div_ru;
  using detail::kInf;
  using detail::SignClass;
  if (is_empty(a) || is_empty(b)) return empty();
  const SignClass cb = detail::classify(b);
  if (cb == SignClass::kZero) return empty();
  const SignClass ca = detail::classify(a);
  if (ca == SignClass::kZero) return point(0.0);

  if (cb == SignClass::kPos && b.lo > 0.0) {  // b strictly positive
    switch (ca) {
      case SignClass::kPos:
        return make(div_rd(a.lo, b.hi), div_ru(a.hi, b.lo));
      case SignClass::kNeg:
        return make(div_rd(a.lo, b.lo), div_ru(a.hi, b.hi));
      default:
        return make(div_rd(a.lo, b.lo), div_ru(a.hi, b.lo));
    }
  }
  if (cb == SignClass::kNeg && b.hi < 0.0) {  // b strictly negative
    switch (ca) {
      case SignClass::kPos:
        return make(div_rd(a.hi, b.hi), div_ru(a.lo, b.lo));
      case SignClass::kNeg:
        return make(div_rd(a.hi, b.lo), div_ru(a.lo, b.hi));
      default:
        return make(div_rd(a.hi, b.hi), div_ru(a.lo, b.hi));
    }
  }
  // 0 is a member of b (endpoint or interior); a has nonzero elements.
  if (cb == SignClass::kMix) return entire();
  if (cb == SignClass::kPos) {  // b = [0, b.hi], b.hi > 0
    switch (ca) {
      case SignClass::kPos:
        return make(div_rd(a.lo, b.hi), kInf);
      case SignClass::kNeg:
        return make(-kInf, div_ru(a.hi, b.hi));
      default:
        return entire();
    }
  }
  // b = [b.lo, 0], b.lo < 0
  switch (ca) {
    case SignClass::kPos:
      return make(-kInf, div_ru(a.lo, b.lo));
    case SignClass::kNeg:
      return make(div_rd(a.hi, b.lo), kInf);
    default:
      return entire();
  }
}

inline Interval recip(const Interval& x) noexcept { return div(point(1.0), x); }

// Image of {u^2}: tighter than mul(x, x) because u ranges over one variable.
inline Interval sqr(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  const double lo_abs = mig(x);
  const double hi_abs = mag(x);
  const double lo = lo_abs == 0.0 ? 0.0 : detail::mul_rd(lo_abs, lo_abs);
  const double hi = hi_abs == detail::kInf ? detail::kInf : detail::mul_ru(hi_abs, hi_abs);
  return make(lo, hi);
}

// Image of {sqrt(u)} over u >= 0: domain restriction clips the negative part.
inline Interval sqrt(const Interval& x) noexcept {
  if (is_empty(x) || x.hi < 0.0) return empty();
  const double lo = x.lo <= 0.0 ? 0.0 : detail::sqrt_rd(x.lo);
  return make(lo, detail::sqrt_ru(x.hi));
}

// Image of {u*v + w} over the box: since w is independent, each bound is the
// corresponding product-bound fused with the matching bound of z in a single
// directed rounding -- the same sign-class pairing as mul.
inline Interval fma(const Interval& a, const Interval& b, const Interval& z) noexcept {
  using detail::fma_rd;
  using detail::fma_ru;
  using detail::SignClass;
  if (is_empty(a) || is_empty(b) || is_empty(z)) return empty();
  const SignClass ca = detail::classify(a);
  const SignClass cb = detail::classify(b);
  if (ca == SignClass::kZero || cb == SignClass::kZero) return z;

  switch (ca) {
    case SignClass::kPos:
      switch (cb) {
        case SignClass::kPos:
          return make(fma_rd(a.lo, b.lo, z.lo), fma_ru(a.hi, b.hi, z.hi));
        case SignClass::kNeg:
          return make(fma_rd(a.hi, b.lo, z.lo), fma_ru(a.lo, b.hi, z.hi));
        default:
          return make(fma_rd(a.hi, b.lo, z.lo), fma_ru(a.hi, b.hi, z.hi));
      }
    case SignClass::kNeg:
      switch (cb) {
        case SignClass::kPos:
          return make(fma_rd(a.lo, b.hi, z.lo), fma_ru(a.hi, b.lo, z.hi));
        case SignClass::kNeg:
          return make(fma_rd(a.hi, b.hi, z.lo), fma_ru(a.lo, b.lo, z.hi));
        default:
          return make(fma_rd(a.lo, b.hi, z.lo), fma_ru(a.lo, b.lo, z.hi));
      }
    default:
      switch (cb) {
        case SignClass::kPos:
          return make(fma_rd(a.lo, b.hi, z.lo), fma_ru(a.hi, b.hi, z.hi));
        case SignClass::kNeg:
          return make(fma_rd(a.hi, b.lo, z.lo), fma_ru(a.lo, b.lo, z.hi));
        default:
          return make(std::min(fma_rd(a.lo, b.hi, z.lo), fma_rd(a.hi, b.lo, z.lo)),
                      std::max(fma_ru(a.lo, b.lo, z.hi), fma_ru(a.hi, b.hi, z.hi)));
      }
  }
}

// --- Pointwise monotone / piecewise-exact functions -------------------------

inline Interval abs(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  return make(mig(x), mag(x));  // |.| is exact on doubles
}

inline Interval min(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a) || is_empty(b)) return empty();
  return make(std::min(a.lo, b.lo), std::min(a.hi, b.hi));
}

inline Interval max(const Interval& a, const Interval& b) noexcept {
  if (is_empty(a) || is_empty(b)) return empty();
  return make(std::max(a.lo, b.lo), std::max(a.hi, b.hi));
}

// floor/ceil/trunc/round are exact, monotone step functions: apply to both
// endpoints. (Infinite endpoints pass through unchanged.)
inline Interval floor(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  return make(std::floor(x.lo), std::floor(x.hi));
}

inline Interval ceil(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  return make(std::ceil(x.lo), std::ceil(x.hi));
}

inline Interval trunc(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  return make(std::trunc(x.lo), std::trunc(x.hi));
}

inline Interval round_ties_to_even(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  // nearbyint honors the current rounding mode, which pyintval never changes
  // from round-to-nearest-even.
  return make(std::nearbyint(x.lo), std::nearbyint(x.hi));
}

inline Interval round_ties_to_away(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  return make(std::round(x.lo), std::round(x.hi));
}

inline Interval sign(const Interval& x) noexcept {
  if (is_empty(x)) return empty();
  const auto sgn = [](double v) noexcept { return v > 0.0 ? 1.0 : (v < 0.0 ? -1.0 : 0.0); };
  return make(sgn(x.lo), sgn(x.hi));
}

// --- Integer powers ---------------------------------------------------------

namespace detail {
// v^n (n >= 1) with directed rounding, by binary exponentiation over the
// magnitude. For v >= 0 every factor is nonnegative, so chaining mul_rd
// (resp. mul_ru) at each step keeps the running product a lower (resp.
// upper) bound of the exact power: valid, within a few ulps of tightest.
inline double pown_scalar(double v, int n, bool round_up) noexcept {
  const bool negative = v < 0.0;
  double base = std::fabs(v);
  // A negative base with odd n flips the result sign, which flips the
  // rounding direction of the magnitude computation.
  const bool odd = (n % 2) != 0;
  const bool sign_flip = negative && odd;
  const bool mag_up = sign_flip ? !round_up : round_up;
  double acc = 1.0;
  bool acc_used = false;
  while (n > 0) {
    if (n & 1) {
      acc = acc_used ? (mag_up ? mul_ru(acc, base) : mul_rd(acc, base)) : base;
      acc_used = true;
    }
    n >>= 1;
    if (n > 0) base = mag_up ? mul_ru(base, base) : mul_rd(base, base);
  }
  return sign_flip ? -acc : acc;
}
}  // namespace detail

// Image of {u^n} for integer n (1788 pown). n = 0 gives [1,1] on any
// nonempty input (the standard defines pown(u, 0) = 1 including at u = 0);
// odd n is monotone; even n factors through |u|; negative n composes with
// the division semantics (so 0 in x widens toward unbounded results).
inline Interval pown(const Interval& x, int n) noexcept {
  if (is_empty(x)) return empty();
  if (n == 0) return point(1.0);
  if (n < 0) {
    // n == INT_MIN would overflow -n; such exponents are rejected upstream.
    return div(point(1.0), pown(x, -n));
  }
  if (n == 1) return x;
  if (n % 2 != 0) {
    return make(detail::pown_scalar(x.lo, n, /*round_up=*/false),
                detail::pown_scalar(x.hi, n, /*round_up=*/true));
  }
  const double lo_abs = mig(x);
  const double hi_abs = mag(x);
  const double lo = lo_abs == 0.0 ? 0.0 : detail::pown_scalar(lo_abs, n, false);
  const double hi = hi_abs == detail::kInf ? detail::kInf : detail::pown_scalar(hi_abs, n, true);
  return make(lo, hi);
}

// --- Cancellative addition/subtraction (1788 9.2) ---------------------------

// cancelMinus(a, b): the tightest z with b + z == a, when a is a translate of
// a superset of b (wid(a) >= wid(b)); otherwise no such z exists and the
// result is entire. Unbounded operands never qualify.
inline Interval cancel_minus(const Interval& a, const Interval& b) noexcept {
  // Empty a: z = empty satisfies b + z == a whenever b is empty or bounded;
  // an unbounded b admits no z at all. Nonempty a with empty or unbounded
  // operand likewise admits no z: the undefined cases yield entire.
  if (is_empty(a)) return (is_empty(b) || is_common(b)) ? empty() : entire();
  if (is_empty(b) || !is_common(a) || !is_common(b)) return entire();
  // Exact test wid(a) < wid(b), i.e. a.hi + b.lo < a.lo + b.hi, decided
  // exactly via TwoSum expansions (directed arithmetic could hide a
  // borderline violation inside its rounding).
  const detail::EftPair lhs = detail::two_sum(a.hi, b.lo);
  const detail::EftPair rhs = detail::two_sum(a.lo, b.hi);
  // If either sum overflows (operands near +-max), the exact comparison is
  // undecidable here; entire is the sound fallback (it is a superset of any
  // valid cancelMinus result).
  if (!std::isfinite(lhs.val) || !std::isfinite(rhs.val)) return entire();
  if (detail::exact_sum_sign({lhs.val, lhs.err, -rhs.val, -rhs.err}) < 0) return entire();
  const double lo = detail::sub_rd(a.lo, b.lo);
  const double hi = detail::sub_ru(a.hi, b.hi);
  if (lo > hi) return entire();  // defensive; the exact test should preclude this
  return make(lo, hi);
}

inline Interval cancel_plus(const Interval& a, const Interval& b) noexcept {
  return cancel_minus(a, neg(b));
}

// --- Reverse operations (constraint propagation) ----------------------------

// mulRevToPair(b, c): the set {x : exists v in b with v*x in c}, which is a
// union of at most two intervals (two half-lines when 0 is interior to b but
// not in c). Returned in increasing order; unused slots are empty.
inline std::pair<Interval, Interval> mul_rev_to_pair(const Interval& b,
                                                     const Interval& c) noexcept {
  using detail::div_rd;
  using detail::div_ru;
  using detail::kInf;
  if (is_empty(b) || is_empty(c)) return {empty(), empty()};
  const bool zero_in_c = contains_zero(c);
  if (b.lo == 0.0 && b.hi == 0.0) {
    // v = 0 forces v*x = 0: every x works iff 0 in c.
    return {zero_in_c ? entire() : empty(), empty()};
  }
  if (!contains_zero(b)) return {div(c, b), empty()};
  if (zero_in_c) return {entire(), empty()};
  // 0 in b (as endpoint or interior), 0 not in c: solutions require |x| large
  // enough that some nonzero v in b reaches c.
  if (b.lo < 0.0 && b.hi > 0.0) {
    if (c.lo > 0.0) {
      return {make(-kInf, div_ru(c.lo, b.lo)), make(div_rd(c.lo, b.hi), kInf)};
    }
    return {make(-kInf, div_ru(c.hi, b.hi)), make(div_rd(c.hi, b.lo), kInf)};
  }
  if (b.lo == 0.0) {  // b = [0, b.hi], b.hi > 0
    if (c.lo > 0.0) return {make(div_rd(c.lo, b.hi), kInf), empty()};
    return {make(-kInf, div_ru(c.hi, b.hi)), empty()};
  }
  // b = [b.lo, 0], b.lo < 0
  if (c.lo > 0.0) return {make(-kInf, div_ru(c.lo, b.lo)), empty()};
  return {make(div_rd(c.hi, b.lo), kInf), empty()};
}

// mulRev(b, c, x): interval hull of x intersected with {x : b*x meets c}.
inline Interval mul_rev(const Interval& b, const Interval& c, const Interval& x) noexcept {
  const auto [p1, p2] = mul_rev_to_pair(b, c);
  return convex_hull(intersection(p1, x), intersection(p2, x));
}
inline Interval mul_rev(const Interval& b, const Interval& c) noexcept {
  return mul_rev(b, c, entire());
}

// sqrRev(c, x): hull of x intersected with {u : u^2 in c}. The solution set
// is symmetric: |u| in sqrt(c intersect [0, inf)), widened outward so the
// enclosure is a superset of the exact solution set.
inline Interval sqr_rev(const Interval& c, const Interval& x) noexcept {
  const Interval cc = intersection(c, make(0.0, detail::kInf));
  if (is_empty(cc)) return empty();
  const double r_lo = cc.lo <= 0.0 ? 0.0 : detail::sqrt_rd(cc.lo);
  const double r_hi = detail::sqrt_ru(cc.hi);
  const Interval pos = make(r_lo, r_hi);
  const Interval negp = make(-r_hi, -r_lo);
  return convex_hull(intersection(negp, x), intersection(pos, x));
}
inline Interval sqr_rev(const Interval& c) noexcept { return sqr_rev(c, entire()); }

// absRev(c, x): hull of x intersected with {u : |u| in c}.
inline Interval abs_rev(const Interval& c, const Interval& x) noexcept {
  const Interval cc = intersection(c, make(0.0, detail::kInf));
  if (is_empty(cc)) return empty();
  const Interval pos = cc;
  const Interval negp = make(-cc.hi, -cc.lo);
  return convex_hull(intersection(negp, x), intersection(pos, x));
}
inline Interval abs_rev(const Interval& c) noexcept { return abs_rev(c, entire()); }

// --- Diagnostics ------------------------------------------------------------

// Debug/diagnostic output for tests; user-facing formatting lives in text.hpp.
inline std::ostream& operator<<(std::ostream& os, const Interval& x) {
  if (is_empty(x)) return os << "[empty]";
  const auto old_flags = os.flags();
  os.precision(17);
  os << "[" << x.lo << ", " << x.hi << "]";
  os.flags(old_flags);
  return os;
}

}  // namespace pyintval
