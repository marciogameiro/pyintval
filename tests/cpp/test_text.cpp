// Tests for correctly-rounded text <-> interval conversion (text.hpp).
// The central claim under test: literal endpoints are converted with correct
// DIRECTED rounding, so the interval provably contains the exact real value.

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <string>

#include "doctest/doctest.h"
#include "pyintval/text.hpp"

using namespace pyintval;

namespace {
constexpr double kInfv = std::numeric_limits<double>::infinity();
constexpr double kMax = std::numeric_limits<double>::max();
constexpr double kDMin = std::numeric_limits<double>::denorm_min();

Interval parse(const char* s) {
  Interval x;
  bool ok = text_to_interval(s, x);
  REQUIRE(ok);
  return x;
}
bool rejects(const char* s) {
  Interval x;
  return !text_to_interval(s, x);
}
bool bit_equal(const Interval& a, const Interval& b) {
  uint64_t la, lb, ha, hb;
  std::memcpy(&la, &a.lo, 8);
  std::memcpy(&lb, &b.lo, 8);
  std::memcpy(&ha, &a.hi, 8);
  std::memcpy(&hb, &b.hi, 8);
  // Treat the two zero bit-patterns as equal (kernel canonicalizes to +0).
  auto z = [](uint64_t& u, double d) {
    if (d == 0.0) u = 0;
  };
  z(la, a.lo);
  z(lb, b.lo);
  z(ha, a.hi);
  z(hb, b.hi);
  return la == lb && ha == hb;
}
}  // namespace

TEST_CASE("exact literals produce exact (zero-width) endpoints") {
  Interval x = parse("0.5");
  CHECK(x.lo == 0.5);
  CHECK(x.hi == 0.5);
  x = parse("[0.25, 2]");
  CHECK(x.lo == 0.25);
  CHECK(x.hi == 2.0);
  x = parse("[-3, 3]");
  CHECK(x.lo == -3.0);
  CHECK(x.hi == 3.0);
}

TEST_CASE("the canonical 0.1 straddles 1/10 by exactly one ulp") {
  Interval x = parse("0.1");
  CHECK(x.lo == 0x1.9999999999999p-4);
  CHECK(x.hi == 0x1.999999999999ap-4);
  CHECK(x.lo < x.hi);
  CHECK(detail::succ(x.lo) == x.hi);
  // "[0.1, 0.1]" denotes the same rigorous 1-ulp enclosure.
  CHECK(bit_equal(parse("[0.1, 0.1]"), x));
  CHECK(bit_equal(parse("[0.1]"), x));
}

TEST_CASE("overflow literals round outward to [DBL_MAX, inf]") {
  Interval x = parse("1e400");
  CHECK(x.lo == kMax);
  CHECK(x.hi == kInfv);
  x = parse("-1e400");
  CHECK(x.lo == -kInfv);
  CHECK(x.hi == -kMax);
  x = parse("1e309");
  CHECK(x.lo == kMax);
  CHECK(x.hi == kInfv);
  // Just below the overflow boundary: still finite, both endpoints finite.
  x = parse("9.9e307");
  CHECK(std::isfinite(x.lo));
  CHECK(std::isfinite(x.hi));
  CHECK(x.lo <= x.hi);
}

TEST_CASE("underflow literals round outward across denorm_min") {
  Interval x = parse("1e-400");
  CHECK(x.lo == 0.0);
  CHECK(x.hi == kDMin);
  x = parse("-1e-400");
  CHECK(x.lo == -kDMin);
  CHECK(x.hi == 0.0);
  // 4.9e-324 is just below denorm_min (~4.9407e-324): bracket [0, denorm_min].
  x = parse("4.9e-324");
  CHECK(x.lo == 0.0);
  CHECK(x.hi == kDMin);
}

TEST_CASE("subnormal literals convert with correct direction") {
  // 2.5e-323 is not representable, so the bracket is exactly one ulp wide with
  // both endpoints positive subnormals. (Comparing against the C++ literal
  // 2.5e-323 would be meaningless -- the compiler rounds it too.)
  Interval x = parse("2.5e-323");
  CHECK(x.lo > 0.0);
  CHECK(x.hi > x.lo);
  CHECK(detail::succ(x.lo) == x.hi);
  CHECK(x.hi < 1e-322);
  x = parse("1e-323");  // 2 * denorm_min exactly representable? 1e-323 is not
  CHECK(x.lo <= x.hi);
  CHECK(x.lo >= 0.0);
}

TEST_CASE("hex-float literals: exact when representable, bracketed otherwise") {
  Interval x = parse("0x1.8p1");  // 3.0 exactly
  CHECK(x.lo == 3.0);
  CHECK(x.hi == 3.0);
  x = parse("0x1p-1074");  // denorm_min exactly
  CHECK(x.lo == kDMin);
  CHECK(x.hi == kDMin);
  x = parse("-0x1.fffffffffffffp+1023");  // -DBL_MAX exactly
  CHECK(x.lo == -kMax);
  CHECK(x.hi == -kMax);
  // More bits than fit in a double: must straddle by one ulp.
  x = parse("0x1.00000000000008p0");  // 1 + 2^-53, the midpoint above 1
  CHECK(x.lo == 1.0);
  CHECK(x.hi == detail::succ(1.0));
  CHECK(x.lo < x.hi);
}

TEST_CASE("long digit strings force the exact bignum path") {
  // 100 digits of 1/3: value is below 1/3, so the nearest double may lie either
  // side; the enclosure must contain the true value and be at most one ulp.
  std::string third = "0.";
  third.append(100, '3');
  Interval x = parse(third.c_str());
  CHECK(x.lo <= x.hi);
  CHECK(x.hi <= detail::succ(x.lo));
  CHECK(x.lo < 1.0 / 3.0 + 1e-16);
  // 1e-300 written with a very long exact tail still brackets correctly.
  std::string big = "1";
  big.append(300, '0');
  Interval y = parse(big.c_str());  // 10^300, representable range
  CHECK(y.lo <= 1e300);
  CHECK(1e300 <= y.hi);
}

TEST_CASE("exponent and sign spelling variants") {
  CHECK(bit_equal(parse("1e5"), parse("1E5")));
  CHECK(bit_equal(parse("1e5"), parse("+1e5")));
  CHECK(bit_equal(parse("100000"), parse("1e5")));
  CHECK(parse("1e5").lo == 100000.0);
  CHECK(parse(".5").lo == 0.5);
  CHECK(parse("1.").lo == 1.0);
  CHECK(parse("0").lo == 0.0);
  CHECK(parse("0").hi == 0.0);
  CHECK(parse("-0").lo == 0.0);  // canonicalized to +0
  CHECK(parse("0e999").lo == 0.0);
}

TEST_CASE("keyword forms") {
  CHECK(is_entire(parse("[entire]")));
  CHECK(is_empty(parse("[empty]")));
  CHECK(is_empty(parse("[]")));
  CHECK(is_empty(parse("[ ]")));
  CHECK(is_entire(parse("[-inf, inf]")));
  Interval x = parse("[3, inf]");
  CHECK(x.lo == 3.0);
  CHECK(x.hi == kInfv);
  x = parse("[  -inf ,  -2 ]");  // whitespace tolerance
  CHECK(x.lo == -kInfv);
  CHECK(x.hi == -2.0);
}

TEST_CASE("malformed inputs are rejected") {
  const char* bad[] = {"",        "inf",  "nan",   "[inf,inf]", "[-inf,-inf]", "[2,1]",
                       "1.2.3",   "1e",   "1e+",   "0x",        ".",           "e5",
                       "[1,2",    "1,2]", "--1",   "1 2",       "0x1.p",       "abc",
                       "[1,2,3]", "1/3",  "3.5?x", "[nan]",     "[1, nan]"};
  for (const char* s : bad) {
    CAPTURE(s);
    CHECK(rejects(s));
  }
}

TEST_CASE("uncertain form d?ruE (F4)") {
  // Exactly representable endpoints are pinned bit-for-bit.
  CHECK(bit_equal(parse("-10?"), make(-10.5, -9.5)));   // half-ulp default radius
  CHECK(bit_equal(parse("-10?u"), make(-10.0, -9.5)));  // one-sided up
  CHECK(bit_equal(parse("-10?12"), make(-22.0, 2.0)));  // radius of 12 ulps
  CHECK(bit_equal(parse("3.56?1e2"), make(355.0, 357.0)));
  CHECK(parse("0.0?u").lo == 0.0);  // up starts exactly at the midpoint
  CHECK(parse("0.0?d").hi == 0.0);
  // Inexact endpoints: the exact midpoint +- radius is rigorously enclosed.
  const Interval x = parse("3.56?1");  // [3.55, 3.57]
  CHECK((x.lo <= 3.55 && 3.57 <= x.hi));
  CHECK(is_member(3.56, x));
  // Infinite radius (double '?'): the whole line, or a half-line at the midpoint.
  CHECK(is_entire(parse("0.0??")));
  CHECK(bit_equal(parse("2.5??u"), make(2.5, detail::kInf)));
  CHECK(bit_equal(parse("-10??d"), make(-detail::kInf, -10.0)));
  // A radius with too many digits to represent saturates to an infinite radius.
  const std::string huge = "10?" + std::string(320, '9');
  CHECK(is_entire(parse(huge.c_str())));
}

TEST_CASE("half-bounded and empty inf-sup literals (F4)") {
  CHECK(is_entire(parse("[,]")));
  CHECK(bit_equal(parse("[1,]"), make(1.0, detail::kInf)));
  CHECK(bit_equal(parse("[,2]"), make(-detail::kInf, 2.0)));
  CHECK(bit_equal(parse("[-1,]"), make(-1.0, detail::kInf)));
}

// Enclosure is the rigorous I/O contract: re-parsing formatted text yields a
// superset of the original (outward rounding applies to serialization too),
// widening each finite endpoint by at most one ulp. Exactly-representable
// endpoints round-trip with no widening at all.
namespace {
bool encloses_tightly(const Interval& orig, const Interval& re) {
  if (!subset(orig, re)) return false;
  if (!is_empty(orig)) {
    if (re.lo != -kInfv && re.lo < detail::pred(orig.lo)) return false;
    if (re.hi != kInfv && re.hi > detail::succ(orig.hi)) return false;
  }
  return true;
}
}  // namespace

TEST_CASE("round-trip through text encloses the original within one ulp") {
  const Interval table[] = {make(0.0, 0.0),
                            make(1.0, 2.0),
                            make(-3.5, 3.5),
                            make(0.1, 0.2),
                            empty(),
                            entire(),
                            make(-kInfv, 2.0),
                            make(2.0, kInfv),
                            make(kDMin, 3.0 * kDMin),
                            make(-kMax, kMax),
                            make(kMax, kInfv),
                            point(0.1),
                            point(kDMin),
                            point(-kDMin),
                            make(0.0, kDMin),
                            point(std::nextafter(1.0, 2.0))};
  for (const Interval& x : table) {
    std::string s = interval_to_text(x);
    Interval y;
    CAPTURE(s);
    REQUIRE(text_to_interval(s, y));
    CHECK(encloses_tightly(x, y));
  }
  // Endpoints whose shortest decimal is the exact value round-trip with no
  // widening (small integers and short dyadic fractions; NOT huge values like
  // DBL_MAX whose decimal form is inexact).
  for (const Interval& x : {make(1.0, 2.0), make(-3.5, 3.5), make(0.5, 0.25 + 1.75), point(0.0),
                            make(2.0, kInfv), entire(), empty()}) {
    Interval y;
    REQUIRE(text_to_interval(interval_to_text(x).c_str(), y));
    CHECK(bit_equal(x, y));
  }
}

TEST_CASE("randomized enclosure round-trip on point intervals (fixed seed)") {
  std::mt19937_64 rng(0x7E5701788ULL);
  int checked = 0;
  for (int i = 0; i < 10000; ++i) {
    uint64_t bits = rng();
    double d;
    std::memcpy(&d, &bits, 8);
    if (!std::isfinite(d)) continue;
    ++checked;
    Interval p = point(d);
    std::string s = interval_to_text(p);
    Interval y;
    CAPTURE(s);
    REQUIRE(text_to_interval(s, y));
    CHECK(encloses_tightly(p, y));
  }
  CHECK(checked > 9000);
}

TEST_CASE("hex-float serialization round-trips bit-exactly") {
  // The exact (lossless) serialization path: hex floats are dyadic, so
  // directed parsing recovers them with zero widening.
  std::mt19937_64 rng(0x1788C0FFEEULL);
  for (int i = 0; i < 10000; ++i) {
    uint64_t bits = rng();
    double d;
    std::memcpy(&d, &bits, 8);
    if (!std::isfinite(d)) continue;
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%a", d);  // exact hex float
    Interval y;
    CAPTURE(buf);
    REQUIRE(text_to_interval(buf, y));
    CHECK(y.lo == d);
    CHECK(y.hi == d);
  }
}

TEST_CASE("enclosure property: parsed interval contains the printed midpoint") {
  // For a spread of decimal strings, the bracket must contain the value that
  // strtod (correctly rounded) produces, and be at most one ulp wide.
  const char* lits[] = {"0.1",
                        "0.2",
                        "0.3",
                        "3.14159265358979",
                        "2.718281828459045",
                        "1.1",
                        "123.456",
                        "9.999999999999999e-2",
                        "6.022e23",
                        "1.602176634e-19",
                        "0.0000001",
                        "987654321.123456789"};
  for (const char* s : lits) {
    Interval x = parse(s);
    double mid = std::strtod(s, nullptr);
    CAPTURE(s);
    CHECK(x.lo <= mid);
    CHECK(mid <= x.hi);
    CHECK(x.hi <= detail::succ(x.lo));  // at most one ulp wide
  }
}
