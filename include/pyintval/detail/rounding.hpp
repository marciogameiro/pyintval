#pragma once

// Directed-rounding primitives via error-free transformations (EFT).
//
// Strategy: stay in the default round-to-nearest (RN) mode and recover the
// exact rounded-down (RD) / rounded-up (RU) result of each basic operation
// from the RN result plus the *sign* of the rounding error, which an EFT
// computes exactly. No fesetround, no global floating-point state: the
// primitives are thread-safe and immune to compiler reordering around mode
// switches.
//
// Provenance of algorithms and thresholds:
//  - TwoSum: D.E. Knuth (TAOCP vol. 2, 4.2.2); the error of RN addition is
//    itself a double and TwoSum computes it exactly, with no precondition on
//    the magnitudes of the operands. Boldo/Graillat/Muller showed no
//    intermediate step overflows when the rounded sum is finite.
//  - TwoProdFMA: Ogita/Rump/Oishi, "Accurate sum and dot product",
//    SIAM J. Sci. Comput. 2005; exact when the product does not underflow.
//  - Underflow-safe scaling branches and thresholds: S.M. Rump's rounding
//    emulation scheme as implemented in M. Kashiwagi's kv library
//    (rdouble-nohwround.hpp, MIT) and matsueushi/RoundingEmulator.jl (MIT).
//    Each threshold is explained where used.
//
// The unit tests cross-validate every primitive bit-for-bit against hardware
// directed rounding (fesetround) over deterministic edge-case tables and
// large random samples, including subnormals, overflow, and signed zeros.
//
// Convention used throughout: RD(x op y) is the largest double <= the exact
// result, RU(x op y) the smallest double >= it. Signs of zero results may
// differ from a true directed-rounding unit; interval code canonicalizes
// zero endpoints, so only the *value* is load-bearing.

#include <array>
#include <cfloat>
#include <cmath>
#include <limits>

#if defined(__FAST_MATH__)
#error "pyintval: -ffast-math / -Ofast breaks the enclosure guarantees; compile without it."
#endif
#if defined(FLT_EVAL_METHOD) && FLT_EVAL_METHOD != 0 && FLT_EVAL_METHOD != 1
// FLT_EVAL_METHOD == 2 (x87 extended-precision evaluation) double-rounds
// every operation, which invalidates the EFT proofs.
#error \
    "pyintval: double arithmetic must be evaluated in double (FLT_EVAL_METHOD 0 or 1). On 32-bit x86 use -msse2 -mfpmath=sse."
#endif

namespace pyintval::detail {

inline constexpr double kInf = std::numeric_limits<double>::infinity();
inline constexpr double kMaxDouble = std::numeric_limits<double>::max();
inline constexpr double kDenormMin = std::numeric_limits<double>::denorm_min();

// Exact successor/predecessor on the double grid. std::nextafter is used
// deliberately instead of Rump's inline succ/pred: Rump's formula contains a
// contractible a*b + c shape that -ffp-contract could fuse into an fma,
// silently changing results; nextafter has no such hazard. (succ(+inf) stays
// +inf and pred(-inf) stays -inf, exactly what unbounded endpoints need.)
inline double succ(double x) noexcept { return std::nextafter(x, kInf); }
inline double pred(double x) noexcept { return std::nextafter(x, -kInf); }

struct EftPair {
  double val;  // RN result of the operation
  double err;  // exact error: true result == val + err (when the EFT is exact)
};

// TwoSum (Knuth): a + b == val + err exactly, for all finite a, b.
inline EftPair two_sum(double a, double b) noexcept {
  const double s = a + b;
  const double bb = s - a;
  const double err = (a - (s - bb)) + (b - bb);
  return {s, err};
}

// FastTwoSum (Dekker): requires |a| >= |b| (or a == 0); 3 flops.
inline EftPair fast_two_sum(double a, double b) noexcept {
  const double s = a + b;
  const double err = b - (s - a);
  return {s, err};
}

// TwoProdFMA (Ogita/Rump/Oishi): a * b == val + err exactly, provided the
// exact product does not fall below 2^-1022 - 2^-1074 in magnitude (the err
// term would underflow). Callers gate on |val| >= 2^-969 which is safely
// inside the exact regime, or rescale first.
inline EftPair two_prod(double a, double b) noexcept {
  const double p = a * b;
  const double err = std::fma(a, b, -p);
  return {p, err};
}

// ---------------------------------------------------------------------------
// Addition / subtraction.
//
// Structure (all four): the RN sum plus the sign of the exact TwoSum error
// determines the directed result. Overflow: if RN(a+b) == +inf but a, b are
// finite, the exact sum lies in (maxDouble, +inf), so RD is maxDouble and RU
// is +inf; mirrored for -inf. Infinite *operands* pass their (exact) result
// through.
// ---------------------------------------------------------------------------

inline double add_ru(double a, double b) noexcept {
  const EftPair s = two_sum(a, b);
  if (s.val == kInf) return s.val;
  if (s.val == -kInf) {
    // RN overflowed to -inf: exact sum is finite iff both operands are.
    return (a == -kInf || b == -kInf) ? -kInf : -kMaxDouble;
  }
  return s.err > 0.0 ? succ(s.val) : s.val;
}

inline double add_rd(double a, double b) noexcept {
  const EftPair s = two_sum(a, b);
  if (s.val == -kInf) return s.val;
  if (s.val == kInf) {
    return (a == kInf || b == kInf) ? kInf : kMaxDouble;
  }
  return s.err < 0.0 ? pred(s.val) : s.val;
}

inline double sub_ru(double a, double b) noexcept { return add_ru(a, -b); }
inline double sub_rd(double a, double b) noexcept { return add_rd(a, -b); }

// ---------------------------------------------------------------------------
// Multiplication.
//
// |p| >= 2^-969: TwoProdFMA error is exact (2^-969 = 2^(-1074 + 106 - 1)
// keeps the error term, bounded by ulp(p)/2, above the subnormal floor), so
// its sign decides the adjustment.
//
// |p| < 2^-969 (result near/below the subnormal range, where the fma residual
// itself can be rounded): rescale both operands by 2^537 (squares to 2^1074,
// lifting the product fully out of the subnormal range without overflowing,
// since |a|,|b| < 2^-416 here). The scaled product s + s2 is then exact, and
// comparing it against the scaled-back candidate t = (p * 2^537) * 2^537
// (exact: p is a small power-of-two multiple) decides the direction.
// Following kv's mul_up/mul_down.
// ---------------------------------------------------------------------------

inline double mul_ru(double a, double b) noexcept {
  const EftPair p = two_prod(a, b);
  if (p.val == kInf) return p.val;
  if (p.val == -kInf) {
    return (std::fabs(a) == kInf || std::fabs(b) == kInf) ? -kInf : -kMaxDouble;
  }
  constexpr double kErrorFreeGate = 0x1p-969;
  constexpr double kRescale = 0x1p537;
  if (std::fabs(p.val) >= kErrorFreeGate) {
    return p.err > 0.0 ? succ(p.val) : p.val;
  }
  const EftPair s = two_prod(a * kRescale, b * kRescale);
  const double t = (p.val * kRescale) * kRescale;
  if (t < s.val || (t == s.val && s.err > 0.0)) return succ(p.val);
  return p.val;
}

inline double mul_rd(double a, double b) noexcept {
  const EftPair p = two_prod(a, b);
  if (p.val == -kInf) return p.val;
  if (p.val == kInf) {
    return (std::fabs(a) == kInf || std::fabs(b) == kInf) ? kInf : kMaxDouble;
  }
  constexpr double kErrorFreeGate = 0x1p-969;
  constexpr double kRescale = 0x1p537;
  if (std::fabs(p.val) >= kErrorFreeGate) {
    return p.err < 0.0 ? pred(p.val) : p.val;
  }
  const EftPair s = two_prod(a * kRescale, b * kRescale);
  const double t = (p.val * kRescale) * kRescale;
  if (t > s.val || (t == s.val && s.err < 0.0)) return pred(p.val);
  return p.val;
}

// ---------------------------------------------------------------------------
// Division.
//
// Exact zero / infinite / NaN operands pass through x / y (the FP quotient is
// exact or correctly signed-infinite in those cases). Otherwise normalize the
// divisor positive, and if the dividend is small enough that the correctness
// check below could underflow (|x| < 2^-969), rescale both by 2^105 -- unless
// the divisor is huge (|y| >= 2^918 = 2^(1023 - 105)), where rescaling could
// overflow it; there the quotient magnitude is below 2^-969-2^-918 ~ 2^-1887,
// deep under the subnormal floor, so the directed results are +-denorm_min
// and +-0 by sign alone.
//
// Correctness check: d = RN(x/y); the sign of x - d*y (computed exactly via
// TwoProdFMA, valid because x was rescaled out of the danger zone) tells
// whether d is below or above the exact quotient. Following kv's
// div_up/div_down.
// ---------------------------------------------------------------------------

inline double div_ru(double x, double y) noexcept {
  if (x == 0.0 || y == 0.0 || std::fabs(x) == kInf || std::fabs(y) == kInf || x != x || y != y) {
    return x / y;
  }
  double xn = x, yn = y;
  if (y < 0.0) {
    xn = -x;
    yn = -y;
  }
  constexpr double kSmallDividend = 0x1p-969;
  constexpr double kHugeDivisor = 0x1p918;
  constexpr double kRescale = 0x1p105;
  if (std::fabs(xn) < kSmallDividend) {
    if (std::fabs(yn) < kHugeDivisor) {
      xn *= kRescale;
      yn *= kRescale;
    } else {
      // |quotient| < 2^-1887: rounds to the subnormal floor by sign.
      return xn < 0.0 ? -0.0 : kDenormMin;
    }
  }
  const double d = xn / yn;
  if (d == kInf) return d;
  if (d == -kInf) return -kMaxDouble;
  const EftPair r = two_prod(d, yn);
  if (r.val < xn || (r.val == xn && r.err < 0.0)) return succ(d);
  return d;
}

inline double div_rd(double x, double y) noexcept {
  if (x == 0.0 || y == 0.0 || std::fabs(x) == kInf || std::fabs(y) == kInf || x != x || y != y) {
    return x / y;
  }
  double xn = x, yn = y;
  if (y < 0.0) {
    xn = -x;
    yn = -y;
  }
  constexpr double kSmallDividend = 0x1p-969;
  constexpr double kHugeDivisor = 0x1p918;
  constexpr double kRescale = 0x1p105;
  if (std::fabs(xn) < kSmallDividend) {
    if (std::fabs(yn) < kHugeDivisor) {
      xn *= kRescale;
      yn *= kRescale;
    } else {
      return xn < 0.0 ? -kDenormMin : 0.0;
    }
  }
  const double d = xn / yn;
  if (d == -kInf) return d;
  if (d == kInf) return kMaxDouble;
  const EftPair r = two_prod(d, yn);
  if (r.val > xn || (r.val == xn && r.err > 0.0)) return pred(d);
  return d;
}

// ---------------------------------------------------------------------------
// Square root (x >= 0; negative inputs are excluded at the interval layer).
//
// d = RN(sqrt(x)); the sign of d*d - x (exact via TwoProdFMA) tells whether d
// overshoots or undershoots. For x < 2^-969 the square d*d would fall into
// the inexact-EFT zone, so compare (d * 2^53)^2 against x * 2^106 instead --
// both scalings are exact (power-of-two, no overflow at these magnitudes).
// Following kv's sqrt_up/sqrt_down.
// ---------------------------------------------------------------------------

inline double sqrt_ru(double x) noexcept {
  const double d = std::sqrt(x);
  if (d == kInf) return d;  // sqrt(+inf) or sqrt overflow cannot occur otherwise
  constexpr double kSmallGate = 0x1p-969;
  if (x < kSmallGate) {
    const double x2 = x * 0x1p106;
    const double d2 = d * 0x1p53;
    const EftPair r = two_prod(d2, d2);
    if (r.val < x2 || (r.val == x2 && r.err < 0.0)) return succ(d);
    return d;
  }
  const EftPair r = two_prod(d, d);
  if (r.val < x || (r.val == x && r.err < 0.0)) return succ(d);
  return d;
}

inline double sqrt_rd(double x) noexcept {
  const double d = std::sqrt(x);
  // sqrt of a finite double never overflows, so d == +inf only for x == +inf,
  // where the result is exact.
  if (d == kInf) return d;
  constexpr double kSmallGate = 0x1p-969;
  if (x < kSmallGate) {
    const double x2 = x * 0x1p106;
    const double d2 = d * 0x1p53;
    const EftPair r = two_prod(d2, d2);
    if (r.val > x2 || (r.val == x2 && r.err > 0.0)) return pred(d);
    return d;
  }
  const EftPair r = two_prod(d, d);
  if (r.val > x || (r.val == x && r.err > 0.0)) return pred(d);
  return d;
}

// ---------------------------------------------------------------------------
// Exact sign of a four-term sum.
//
// Distills {v0,v1,v2,v3} with TwoSum sweeps (Ogita-Rump-Oishi VecSum) to a
// fixed point: a non-overlapping expansion in increasing magnitude order,
// whose leading nonzero term then dominates the tail and carries the sign of
// the exact sum. Callers must ensure no intermediate partial sum overflows;
// every use below feeds terms that cancel massively by construction.
// ---------------------------------------------------------------------------

inline int exact_sum_sign(std::array<double, 4> v) noexcept {
  // For two TwoSum pairs convergence takes 2-3 sweeps; the cap is defensive.
  for (int sweep = 0; sweep < 6; ++sweep) {
    bool changed = false;
    for (int i = 1; i < 4; ++i) {
      const EftPair p = two_sum(v[i - 1], v[i]);
      changed = changed || p.val != v[i] || p.err != v[i - 1];
      v[i - 1] = p.err;
      v[i] = p.val;
    }
    if (!changed) break;
  }
  for (int i = 3; i >= 0; --i) {
    if (v[i] > 0.0) return 1;
    if (v[i] < 0.0) return -1;
  }
  return 0;
}

// ---------------------------------------------------------------------------
// Fused multiply-add a*b + c, directed and TIGHT in every regime.
//
// r = RN(fma(a,b,c)) is correctly rounded, so the directed results are r
// itself or its neighbor, decided by the sign of delta = (a*b + c) - r.
// fma_delta_sign computes that sign EXACTLY by reducing to a four-term exact
// sum in a regime-appropriate scaling (every rescaling below is by a power of
// two acting on operands proven large enough that no bits are lost):
//
//  - huge product (|a*b| >= 2^1019, possibly overflowing the double range,
//    though never past 2^1025 since r is finite): scale the larger factor by
//    2^-56 so TwoProdFMA is exact; the huge partner among {c, r} (at least
//    one is huge, since they jointly cancel the product to within ulp(r))
//    joins at the same scale; the cancelled residue is scaled back up
//    exactly and compared against the remaining small term.
//  - ordinary regime (2^-969 <= |a*b| < 2^1019): TwoProdFMA is exact and
//    r - c = a*b - delta is finite, so {u1, u2, -(r-c) pair} distills
//    directly (u1 cancels against r-c, keeping partial sums bounded).
//  - tiny product (|a*b| < 2^-969, where the TwoProdFMA residual would
//    underflow): rescale the product by 2^537 per factor (exact, lifting it
//    out of the subnormal range). If |c| >= 2^-950 then r = RN(c + eps) with
//    |eps| < 2^-968, so c - r is exact by Sterbenz and either dominates the
//    product outright or joins the scaled comparison. Otherwise everything
//    lives below 2^-949 and the whole identity rescales by 2^1074 exactly.
//    Products below 2^-2043 are dwarfed by any nonzero c - r (a multiple of
//    2^-1074) and only tie-break its sign when it is exactly zero.
//
// Tightness matters beyond bragging rights: IEEE 1788 requires tightest
// accuracy for fma, and only a tight (hence monotone) rounding makes the
// interval fma inclusion-isotone.
// ---------------------------------------------------------------------------

inline int sign_of(double x) noexcept { return x > 0.0 ? 1 : (x < 0.0 ? -1 : 0); }

// Sign of the (nonzero) exact product a*b; requires a, b nonzero.
inline int product_sign(double a, double b) noexcept {
  return std::signbit(a) == std::signbit(b) ? 1 : -1;
}

// Exact sign of (a*b + c) - r. Preconditions: a, b nonzero and finite,
// c finite and nonzero, r = RN(a*b + c) finite.
inline int fma_delta_sign(double a, double b, double c, double r) noexcept {
  constexpr double kGate = 0x1p-969;  // exact-TwoProdFMA gate
  constexpr double kS = 0x1p-56;
  const double p = a * b;
  if (!(std::fabs(p) < 0x1p1019)) {                     // huge or overflowing product
    const bool scale_a = std::fabs(a) >= std::fabs(b);  // larger factor is >= 2^509
    const EftPair u = two_prod(scale_a ? a * kS : a, scale_a ? b : b * kS);
    const bool c_huge = std::fabs(c) >= 0x1p-966;  // c * kS exact
    const bool r_huge = std::fabs(r) >= 0x1p-966;  // r * kS exact
    if (c_huge && r_huge) return exact_sum_sign({u.val, u.err, c * kS, -(r * kS)});
    if (c_huge) {  // r below 2^-966: c cancels the product, residue rescales up
      const EftPair t = two_sum(u.val, c * kS);
      return exact_sum_sign(
          {std::ldexp(t.val, 56), std::ldexp(t.err, 56), std::ldexp(u.err, 56), -r});
    }
    // c below 2^-966 forces r to be the huge cancelling partner.
    const EftPair t = two_sum(u.val, -(r * kS));
    return exact_sum_sign({std::ldexp(t.val, 56), std::ldexp(t.err, 56), std::ldexp(u.err, 56), c});
  }
  if (std::fabs(p) >= kGate) {
    const EftPair u = two_prod(a, b);  // exact
    const EftPair d = two_sum(r, -c);  // finite: |r - c| = |a*b - delta| < 2^1020
    return exact_sum_sign({u.val, u.err, -d.val, -d.err});
  }
  // |a*b| below the exact-TwoProdFMA range: rescale the product by 2^1074.
  const EftPair us = two_prod(a * 0x1p537, b * 0x1p537);
  if (std::fabs(c) >= 0x1p-950) {
    const double dd = c - r;  // exact by Sterbenz: r = RN(c + eps), |eps| < 2^-968
    if (std::fabs(dd) > 0x1p-966) return sign_of(dd);  // product (< 2^-968) cannot compete
    if (std::fabs(us.val) < kGate) {
      // |a*b| < 2^-2043 is dwarfed by any nonzero dd (a multiple of 2^-1074).
      return dd != 0.0 ? sign_of(dd) : product_sign(a, b);
    }
    return exact_sum_sign({us.val, us.err, std::ldexp(dd, 1074), 0.0});
  }
  // Everything tiny (|c| < 2^-950 implies |r| < 2^-949): rescale by 2^1074.
  if (std::fabs(us.val) < kGate) {
    const EftPair g = two_sum(c, -r);  // exact; value is a multiple of 2^-1074
    if (g.val != 0.0) return sign_of(g.val);
    if (g.err != 0.0) return sign_of(g.err);
    return product_sign(a, b);
  }
  return exact_sum_sign({us.val, us.err, std::ldexp(c, 1074), -std::ldexp(r, 1074)});
}

inline double fma_ru(double a, double b, double c) noexcept {
  const double r = std::fma(a, b, c);
  if (r == kInf) return r;
  if (r == -kInf) {
    return (std::fabs(a) == kInf || std::fabs(b) == kInf || c == -kInf) ? -kInf : -kMaxDouble;
  }
  // A zero factor makes the fma exact (a*b + c == c up to zero sign); a zero
  // addend reduces to the product, whose directed rounding is exact in every
  // range via mul's rescaling.
  if (a == 0.0 || b == 0.0) return r;
  if (c == 0.0) return mul_ru(a, b);
  return fma_delta_sign(a, b, c, r) > 0 ? succ(r) : r;
}

inline double fma_rd(double a, double b, double c) noexcept {
  const double r = std::fma(a, b, c);
  if (r == -kInf) return r;
  if (r == kInf) {
    return (std::fabs(a) == kInf || std::fabs(b) == kInf || c == kInf) ? kInf : kMaxDouble;
  }
  if (a == 0.0 || b == 0.0) return r;
  if (c == 0.0) return mul_rd(a, b);
  return fma_delta_sign(a, b, c, r) < 0 ? pred(r) : r;
}

}  // namespace pyintval::detail
