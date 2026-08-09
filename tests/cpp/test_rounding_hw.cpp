// THE ORACLE: cross-validate every directed-rounding primitive bit-for-bit
// against the hardware's own directed rounding (fesetround). The EFT
// primitives stay in round-to-nearest and must reproduce exactly what the CPU
// produces in FE_UPWARD / FE_DOWNWARD. Any mismatch is an enclosure-critical
// kernel bug.
//
// This translation unit is compiled with strict FP semantics
// (-frounding-math / /fp:strict; see tests/cpp/CMakeLists.txt) so the compiler
// neither constant-folds nor reorders arithmetic across the mode switches.

#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

#include "doctest/doctest.h"
#include "pyintval/detail/rounding.hpp"

#if defined(__clang__)
#pragma STDC FENV_ACCESS ON
#elif defined(_MSC_VER)
#pragma fenv_access(on)
#endif
// GCC: no pragma; the -frounding-math compile flag governs this TU.

using namespace pyintval::detail;

namespace {

constexpr double kInfv = std::numeric_limits<double>::infinity();
constexpr double kMax = std::numeric_limits<double>::max();
constexpr double kDMin = std::numeric_limits<double>::denorm_min();
constexpr double kMinNormal = std::numeric_limits<double>::min();

bool same(double a, double b) { return (a == b) || (a == 0.0 && b == 0.0); }

// Hardware directed operations via volatile operands so the compiler emits the
// actual instruction under the active rounding mode.
double hw_add(double a, double b, int mode) {
  volatile double va = a, vb = b;
  std::fesetround(mode);
  volatile double r = va + vb;
  std::fesetround(FE_TONEAREST);
  return r;
}
double hw_sub(double a, double b, int mode) {
  volatile double va = a, vb = b;
  std::fesetround(mode);
  volatile double r = va - vb;
  std::fesetround(FE_TONEAREST);
  return r;
}
double hw_mul(double a, double b, int mode) {
  volatile double va = a, vb = b;
  std::fesetround(mode);
  volatile double r = va * vb;
  std::fesetround(FE_TONEAREST);
  return r;
}
double hw_div(double a, double b, int mode) {
  volatile double va = a, vb = b;
  std::fesetround(mode);
  volatile double r = va / vb;
  std::fesetround(FE_TONEAREST);
  return r;
}
double hw_sqrt(double a, int mode) {
  volatile double va = a;
  std::fesetround(mode);
  volatile double r = std::sqrt(va);
  std::fesetround(FE_TONEAREST);
  return r;
}
double hw_fma(double a, double b, double c, int mode) {
  volatile double va = a, vb = b, vc = c;
  std::fesetround(mode);
  volatile double r = std::fma(va, vb, vc);
  std::fesetround(FE_TONEAREST);
  return r;
}

// Deterministic xorshift64 for reproducible CI.
struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed ? seed : 0x1234567890abcdefULL) {}
  uint64_t next() {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
  }
  double dbl() {
    uint64_t b = next();
    double d;
    std::memcpy(&d, &b, sizeof(d));
    return d;
  }
  // Random finite double with unbiased exponent in [lo, hi].
  double scaled(int lo, int hi) {
    double m = 1.0 + static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0);
    int e = lo + static_cast<int>(next() % static_cast<uint64_t>(hi - lo + 1));
    double v = std::ldexp(m, e);
    return (next() & 1) ? v : -v;
  }
};

const double kSpecials[] = {
    0.0,        -0.0,       kDMin,       -kDMin,   3.0 * kDMin, -3.0 * kDMin,   0x1p-1073,
    -0x1p-1073, kMinNormal, -kMinNormal, 0x1p-969, -0x1p-969,   pred(0x1p-969), succ(0x1p-969),
    0x1p-537,   -0x1p-537,  1.0,         -1.0,     succ(1.0),   pred(1.0),      1.5,
    -1.5,       0x1p53,     -0x1p53,     0x1p918,  -0x1p918,    0x1p1023,       -0x1p1023,
    kMax,       -kMax,      kInfv,       -kInfv,   0.1,         -0.1,           3.0,
    -3.0};

}  // namespace

TEST_CASE("succ/pred match std::nextafter") {
  for (double v : kSpecials) {
    if (v != v) continue;
    CHECK(same(succ(v), std::nextafter(v, kInfv)));
    CHECK(same(pred(v), std::nextafter(v, -kInfv)));
  }
  Rng rng(0xABCDEF01);
  for (int i = 0; i < 200000; ++i) {
    double v = rng.dbl();
    if (v != v) continue;
    CHECK(same(succ(v), std::nextafter(v, kInfv)));
    CHECK(same(pred(v), std::nextafter(v, -kInfv)));
  }
}

TEST_CASE("add/sub/mul/div match hardware directed rounding on specials") {
  int mism = 0;
  for (double a : kSpecials) {
    for (double b : kSpecials) {
      double hu, hd;

      hu = hw_add(a, b, FE_UPWARD);
      hd = hw_add(a, b, FE_DOWNWARD);
      if (hu == hu && !same(add_ru(a, b), hu)) ++mism;
      if (hd == hd && !same(add_rd(a, b), hd)) ++mism;

      hu = hw_sub(a, b, FE_UPWARD);
      hd = hw_sub(a, b, FE_DOWNWARD);
      if (hu == hu && !same(sub_ru(a, b), hu)) ++mism;
      if (hd == hd && !same(sub_rd(a, b), hd)) ++mism;

      hu = hw_mul(a, b, FE_UPWARD);
      hd = hw_mul(a, b, FE_DOWNWARD);
      if (hu == hu && !same(mul_ru(a, b), hu)) ++mism;
      if (hd == hd && !same(mul_rd(a, b), hd)) ++mism;

      hu = hw_div(a, b, FE_UPWARD);
      hd = hw_div(a, b, FE_DOWNWARD);
      if (hu == hu && !same(div_ru(a, b), hu)) ++mism;
      if (hd == hd && !same(div_rd(a, b), hd)) ++mism;
    }
  }
  CHECK(mism == 0);
}

TEST_CASE("sqrt matches hardware directed rounding on specials") {
  int mism = 0;
  for (double a : kSpecials) {
    if (!(a >= 0.0)) continue;
    double hu = hw_sqrt(a, FE_UPWARD), hd = hw_sqrt(a, FE_DOWNWARD);
    if (hu == hu && !same(sqrt_ru(a), hu)) ++mism;
    if (hd == hd && !same(sqrt_rd(a), hd)) ++mism;
  }
  CHECK(mism == 0);
}

TEST_CASE("binary ops match hardware over random bit patterns") {
  Rng rng(0x51501788);
  int mism = 0;
  for (int i = 0; i < 200000; ++i) {
    double a = rng.dbl(), b = rng.dbl();
    if (a != a || b != b) continue;
    if (!same(add_ru(a, b), hw_add(a, b, FE_UPWARD))) ++mism;
    if (!same(add_rd(a, b), hw_add(a, b, FE_DOWNWARD))) ++mism;
    if (!same(sub_ru(a, b), hw_sub(a, b, FE_UPWARD))) ++mism;
    if (!same(sub_rd(a, b), hw_sub(a, b, FE_DOWNWARD))) ++mism;
    if (!same(mul_ru(a, b), hw_mul(a, b, FE_UPWARD))) ++mism;
    if (!same(mul_rd(a, b), hw_mul(a, b, FE_DOWNWARD))) ++mism;
    double du = hw_div(a, b, FE_UPWARD), dd = hw_div(a, b, FE_DOWNWARD);
    if (du == du && !same(div_ru(a, b), du)) ++mism;
    if (dd == dd && !same(div_rd(a, b), dd)) ++mism;
  }
  CHECK(mism == 0);
}

TEST_CASE("mul/div match hardware in the subnormal-result regime") {
  Rng rng(0x5AB0);
  int mism = 0;
  for (int i = 0; i < 200000; ++i) {
    double a = rng.scaled(-600, -400), b = rng.scaled(-600, -450);
    if (!same(mul_ru(a, b), hw_mul(a, b, FE_UPWARD))) ++mism;
    if (!same(mul_rd(a, b), hw_mul(a, b, FE_DOWNWARD))) ++mism;
    double num = rng.scaled(-1074, -1000), den = rng.scaled(-40, 40);
    double du = hw_div(num, den, FE_UPWARD), dd = hw_div(num, den, FE_DOWNWARD);
    if (du == du && !same(div_ru(num, den), du)) ++mism;
    if (dd == dd && !same(div_rd(num, den), dd)) ++mism;
  }
  CHECK(mism == 0);
}

TEST_CASE("sqrt matches hardware over random nonnegative doubles") {
  Rng rng(0x59477);
  int mism = 0;
  for (int i = 0; i < 200000; ++i) {
    double a = std::fabs(rng.dbl());
    if (a != a) continue;
    if (!same(sqrt_ru(a), hw_sqrt(a, FE_UPWARD))) ++mism;
    if (!same(sqrt_rd(a), hw_sqrt(a, FE_DOWNWARD))) ++mism;
  }
  CHECK(mism == 0);
}

TEST_CASE("fma matches hardware, incl. cancellation and subnormal regimes") {
  Rng rng(0xF3A0);
  int mism = 0;
  // Uniform random triples.
  for (int i = 0; i < 200000; ++i) {
    double a = rng.dbl(), b = rng.dbl(), c = rng.dbl();
    if (a != a || b != b || c != c) continue;
    double hu = hw_fma(a, b, c, FE_UPWARD), hd = hw_fma(a, b, c, FE_DOWNWARD);
    if (hu == hu && !same(fma_ru(a, b, c), hu)) ++mism;
    if (hd == hd && !same(fma_rd(a, b, c), hd)) ++mism;
  }
  // Massive cancellation: c near -a*b, plus huge and subnormal product regimes.
  for (int i = 0; i < 200000; ++i) {
    double a = rng.scaled(480, 520), b = rng.scaled(480, 520);
    double p = a * b;
    double c = -p;
    for (int k = 0; k < 2; ++k) {
      double hu = hw_fma(a, b, c, FE_UPWARD), hd = hw_fma(a, b, c, FE_DOWNWARD);
      if (hu == hu && !same(fma_ru(a, b, c), hu)) ++mism;
      if (hd == hd && !same(fma_rd(a, b, c), hd)) ++mism;
      c = std::nextafter(c, (rng.next() & 1) ? kInfv : -kInfv);
    }
    double a2 = rng.scaled(-600, -420), b2 = rng.scaled(-600, -450);
    double c2 = rng.scaled(-1074, -930);
    double hu2 = hw_fma(a2, b2, c2, FE_UPWARD), hd2 = hw_fma(a2, b2, c2, FE_DOWNWARD);
    if (hu2 == hu2 && !same(fma_ru(a2, b2, c2), hu2)) ++mism;
    if (hd2 == hd2 && !same(fma_rd(a2, b2, c2), hd2)) ++mism;
  }
  CHECK(mism == 0);
}
