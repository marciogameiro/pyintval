#pragma once

// Correctly-rounded text <-> interval conversion.
//
// The rigor requirement: a decimal (or hex-float) literal endpoint must be
// converted with CORRECT DIRECTED ROUNDING -- a lower endpoint toward -inf, an
// upper endpoint toward +inf -- so that the resulting interval provably
// contains the exact real value the literal denotes. "0.5" yields the exact
// point [0.5, 0.5]; "0.1" yields the tightest 1-ulp interval that truly
// straddles the real number 1/10.
//
// Method: parse the literal's structure exactly into an integer significand
// times 2^p2 * 5^p5 (decimal: p2 = p5 = decimal exponent; hex: p5 = 0). Obtain
// a round-to-nearest guess r from strtod (which also signals overflow ->
// +-inf and underflow -> 0). Then determine the EXACT sign of (value - r) by
// comparing two big integers -- no floating point in the decision -- and step
// r with nextafter to the correctly directed neighbors. strtod need only land
// within a couple of ulps: the exact comparison and the (defensive) stepping
// loop correct any rounding weakness in the platform's strtod.

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

#include "pyintval/interval.hpp"

#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
#include <charconv>
#define PYINTVAL_HAVE_TO_CHARS 1
#else
#include <cstdio>
#define PYINTVAL_HAVE_TO_CHARS 0
#endif

// Locale-independent strtod: a comma-locale must not change how "1.5" parses.
#if defined(_WIN32)
#include <locale.h>
#elif defined(__has_include)
#if __has_include(<xlocale.h>)
#include <xlocale.h>
#endif
#include <locale.h>
#else
#include <locale.h>
#endif

namespace pyintval {
namespace textdetail {

// --- Minimal big unsigned integer (base 2^32 limbs, little-endian) ----------
// Only the operations the exact comparison needs: build from digits, multiply
// by a power of five, shift left (multiply by a power of two), and compare.

struct BigU {
  std::vector<uint32_t> limb;  // little-endian; no trailing zero limbs (0 == empty)

  bool is_zero() const noexcept { return limb.empty(); }
  void trim() noexcept {
    while (!limb.empty() && limb.back() == 0) limb.pop_back();
  }
  static BigU from_u64(uint64_t v) {
    BigU b;
    while (v) {
      b.limb.push_back(static_cast<uint32_t>(v));
      v >>= 32;
    }
    return b;
  }
  // this = this * mul + add, with mul, add < 2^32.
  void mul_add(uint32_t mul, uint32_t add) noexcept {
    uint64_t carry = add;
    for (size_t i = 0; i < limb.size(); ++i) {
      const uint64_t cur = static_cast<uint64_t>(limb[i]) * mul + carry;
      limb[i] = static_cast<uint32_t>(cur);
      carry = cur >> 32;
    }
    while (carry) {
      limb.push_back(static_cast<uint32_t>(carry));
      carry >>= 32;
    }
    trim();
  }
  // this <<= bits  (multiply by 2^bits).
  void shl_bits(unsigned bits) {
    if (is_zero() || bits == 0) return;
    const unsigned whole = bits / 32, frac = bits % 32;
    if (frac) {
      uint32_t carry = 0;
      for (size_t i = 0; i < limb.size(); ++i) {
        const uint64_t cur = (static_cast<uint64_t>(limb[i]) << frac) | carry;
        limb[i] = static_cast<uint32_t>(cur);
        carry = static_cast<uint32_t>(cur >> 32);
      }
      if (carry) limb.push_back(carry);
    }
    if (whole) limb.insert(limb.begin(), whole, 0u);
  }
  static int cmp(const BigU& a, const BigU& b) noexcept {
    if (a.limb.size() != b.limb.size()) return a.limb.size() < b.limb.size() ? -1 : 1;
    for (size_t i = a.limb.size(); i-- > 0;) {
      if (a.limb[i] != b.limb[i]) return a.limb[i] < b.limb[i] ? -1 : 1;
    }
    return 0;
  }
};

// Multiply x by 5^k in chunks (5^13 = 1220703125 < 2^31 fits a uint32).
inline void mul_pow5(BigU& x, int k) {
  static constexpr uint32_t kPow5[14] = {1u,       5u,        25u,        125u,       625u,
                                         3125u,    15625u,    78125u,     390625u,    1953125u,
                                         9765625u, 48828125u, 244140625u, 1220703125u};
  while (k >= 13) {
    x.mul_add(kPow5[13], 0);
    k -= 13;
  }
  if (k > 0) x.mul_add(kPow5[k], 0);
}

// --- Exact representation of a parsed real literal --------------------------
// value = sign * significand * 2^p2 * 5^p5  (significand a nonnegative BigU).

struct ParsedReal {
  bool ok = false;
  int sign = 1;
  bool is_inf = false;
  bool is_zero = false;  // exact zero significand
  BigU sig;
  int p2 = 0;
  int p5 = 0;
  std::string mag_token;  // unsigned magnitude, a valid C literal for strtod
};

// Locale-independent strtod on a NUL-terminated unsigned magnitude token.
inline double c_strtod(const char* s, char** end) {
#if defined(_WIN32)
  static _locale_t loc = _create_locale(LC_ALL, "C");
  return _strtod_l(s, end, loc);
#elif defined(LC_ALL_MASK)
  static locale_t loc = ::newlocale(LC_ALL_MASK, "C", static_cast<locale_t>(0));
  return ::strtod_l(s, end, loc);
#else
  return std::strtod(s, end);
#endif
}

inline int hex_val(char c) noexcept {
  if (c >= '0' && c <= '9') return c - '0';
  if (c >= 'a' && c <= 'f') return c - 'a' + 10;
  if (c >= 'A' && c <= 'F') return c - 'A' + 10;
  return -1;
}

inline std::string_view trim(std::string_view s) noexcept {
  size_t a = 0, b = s.size();
  while (a < b && (s[a] == ' ' || s[a] == '\t')) ++a;
  while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t')) --b;
  return s.substr(a, b - a);
}

inline bool iequal(std::string_view a, const char* lit) noexcept {
  size_t i = 0;
  for (; i < a.size() && lit[i]; ++i) {
    char c = a[i];
    if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
    if (c != lit[i]) return false;
  }
  return i == a.size() && lit[i] == '\0';
}

// Parse a single real literal (decimal or hex-float) exactly. Whitespace must
// already be trimmed by the caller.
inline ParsedReal parse_real(std::string_view s) {
  ParsedReal pr;
  if (s.empty()) return pr;
  size_t i = 0;
  if (s[i] == '+' || s[i] == '-') {
    pr.sign = (s[i] == '-') ? -1 : 1;
    ++i;
  }
  std::string_view body = s.substr(i);
  if (body.empty()) return pr;
  if (iequal(body, "inf") || iequal(body, "infinity")) {
    pr.ok = true;
    pr.is_inf = true;
    return pr;
  }
  if (iequal(body, "nan")) return pr;  // NaN is never a valid endpoint

  pr.mag_token.assign(body);

  // Hex float: 0x<hex>[.<hex>][p<decexp>].
  if (body.size() > 2 && body[0] == '0' && (body[1] == 'x' || body[1] == 'X')) {
    size_t j = 2;
    bool any = false;
    int frac_hex = 0;
    bool seen_dot = false;
    while (j < body.size()) {
      const char c = body[j];
      if (c == '.') {
        if (seen_dot) return pr;
        seen_dot = true;
        ++j;
        continue;
      }
      const int hv = hex_val(c);
      if (hv < 0) break;
      pr.sig.mul_add(16, static_cast<uint32_t>(hv));
      if (seen_dot) ++frac_hex;
      any = true;
      ++j;
    }
    if (!any) return pr;
    int bexp = 0;
    if (j < body.size() && (body[j] == 'p' || body[j] == 'P')) {
      ++j;
      int esign = 1;
      if (j < body.size() && (body[j] == '+' || body[j] == '-')) {
        esign = (body[j] == '-') ? -1 : 1;
        ++j;
      }
      if (j >= body.size() || body[j] < '0' || body[j] > '9') return pr;
      long e = 0;
      while (j < body.size() && body[j] >= '0' && body[j] <= '9') {
        e = e * 10 + (body[j] - '0');
        if (e > 1000000) e = 1000000;  // clamp; strtod handles the actual overflow
        ++j;
      }
      bexp = static_cast<int>(esign * e);
    }
    if (j != body.size()) return pr;  // trailing garbage
    pr.sig.trim();
    pr.is_zero = pr.sig.is_zero();
    pr.p5 = 0;
    pr.p2 = bexp - 4 * frac_hex;
    pr.ok = true;
    return pr;
  }

  // Decimal: [digits][.digits][(e|E)[sign]digits].
  bool any_digit = false;
  int frac_digits = 0;
  bool seen_dot = false;
  size_t j = 0;
  while (j < body.size()) {
    const char c = body[j];
    if (c == '.') {
      if (seen_dot) return pr;
      seen_dot = true;
      ++j;
      continue;
    }
    if (c < '0' || c > '9') break;
    pr.sig.mul_add(10, static_cast<uint32_t>(c - '0'));
    if (seen_dot) ++frac_digits;
    any_digit = true;
    ++j;
  }
  if (!any_digit) return pr;
  int dexp = 0;
  if (j < body.size() && (body[j] == 'e' || body[j] == 'E')) {
    ++j;
    int esign = 1;
    if (j < body.size() && (body[j] == '+' || body[j] == '-')) {
      esign = (body[j] == '-') ? -1 : 1;
      ++j;
    }
    if (j >= body.size() || body[j] < '0' || body[j] > '9') return pr;
    long e = 0;
    while (j < body.size() && body[j] >= '0' && body[j] <= '9') {
      e = e * 10 + (body[j] - '0');
      if (e > 1000000) e = 1000000;
      ++j;
    }
    dexp = static_cast<int>(esign * e);
  }
  if (j != body.size()) return pr;  // trailing garbage
  pr.sig.trim();
  pr.is_zero = pr.sig.is_zero();
  const int e = dexp - frac_digits;  // value = sig * 10^e
  pr.p2 = e;
  pr.p5 = e;
  pr.ok = true;
  return pr;
}

// Exact sign of (magnitude - v) for v > 0 finite, magnitude = sig*2^p2*5^p5.
// Returns -1, 0, or +1.
inline int compare_mag(const BigU& sig, int p2, int p5, double v) {
  if (sig.is_zero()) return v > 0.0 ? -1 : 0;
  // A finite decimal magnitude is always below +inf. Guard this BEFORE frexp /
  // the uint64_t cast (both UB on infinity), and so the defensive stepping loop
  // in bracket() terminates when it steps hi up to +inf instead of spinning.
  if (v == detail::kInf) return -1;
  int ev = 0;
  const double f = std::frexp(v, &ev);                           // v = f * 2^ev, f in [0.5, 1)
  const uint64_t mv = static_cast<uint64_t>(std::ldexp(f, 53));  // exact integer
  ev -= 53;                                                      // v = mv * 2^ev
  BigU L, R;
  if (p5 >= 0) {
    L = sig;
    mul_pow5(L, p5);
    R = BigU::from_u64(mv);
  } else {
    L = sig;
    R = BigU::from_u64(mv);
    mul_pow5(R, -p5);
  }
  // Compare L * 2^p2 vs R * 2^ev by lifting the smaller-exponent side.
  if (p2 >= ev) {
    L.shl_bits(static_cast<unsigned>(p2 - ev));
  } else {
    R.shl_bits(static_cast<unsigned>(ev - p2));
  }
  return BigU::cmp(L, R);
}

// Directed bracket [rd, ru] of the exact value denoted by pr, with rd rounded
// toward -inf and ru toward +inf. Returns false if pr is not ok.
inline bool bracket(const ParsedReal& pr, double& rd, double& ru) {
  using detail::kDenormMin;
  using detail::kInf;
  using detail::kMaxDouble;
  using detail::pred;
  using detail::succ;
  if (!pr.ok) return false;
  if (pr.is_inf) {
    rd = ru = pr.sign < 0 ? -kInf : kInf;
    return true;
  }
  if (pr.is_zero) {
    rd = ru = 0.0;
    return true;
  }
  char* end = nullptr;
  const double guess = c_strtod(pr.mag_token.c_str(), &end);  // magnitude, >= 0
  double lo, hi;                                              // bracket of the positive magnitude
  if (guess == kInf) {                                        // overflow: value magnitude > DBL_MAX
    lo = kMaxDouble;
    hi = kInf;
  } else if (guess == 0.0) {  // underflow: 0 < magnitude <= denorm_min/2
    lo = 0.0;
    hi = kDenormMin;
  } else {
    const int c = compare_mag(pr.sig, pr.p2, pr.p5, guess);
    if (c == 0) {
      lo = hi = guess;
    } else if (c > 0) {  // magnitude > guess
      lo = guess;
      hi = succ(guess);
      while (compare_mag(pr.sig, pr.p2, pr.p5, hi) > 0) {  // defensive
        lo = hi;
        hi = succ(hi);
      }
    } else {  // magnitude < guess
      hi = guess;
      lo = pred(guess);
      while (lo > 0.0 && compare_mag(pr.sig, pr.p2, pr.p5, lo) < 0) {  // defensive
        hi = lo;
        lo = pred(lo);
      }
    }
  }
  if (pr.sign < 0) {
    rd = -hi;
    ru = -lo;
  } else {
    rd = lo;
    ru = hi;
  }
  return true;
}

// IEEE 1788 uncertain form: "m?rad<dir><exp>" -- m is a decimal midpoint, the
// radius is `rad` ulps of m's last decimal place (or half an ulp if `rad` is
// absent), `dir` (u/d) makes it one-sided, and a trailing exponent scales the
// whole thing by 10^exp. Examples: 3.56?1 -> [3.55,3.57], -10? -> [-10.5,-9.5],
// 0.0?2u -> [0.0,0.02], 3.56?1e2 -> [355,357]. Both endpoints are exact decimals,
// so we build "<int>e<exp>" strings and reuse parse_real/bracket for correct
// outward rounding (lower endpoint down, upper endpoint up).
inline bool parse_uncertain(std::string_view s, double& lo, double& hi) {
  const size_t q = s.find('?');
  if (q == std::string_view::npos) return false;
  const std::string_view mant = s.substr(0, q);
  const std::string_view rest = s.substr(q + 1);
  constexpr long long kGuard = 900000000000000000LL;  // keep int64 arithmetic safe

  // Midpoint mantissa: [sign] digits [. digits], with no exponent of its own.
  int sign = 1;
  size_t i = 0;
  if (!mant.empty() && (mant[i] == '+' || mant[i] == '-')) {
    sign = (mant[i] == '-') ? -1 : 1;
    ++i;
  }
  long long m = 0;
  int frac = 0;
  bool seen_dot = false, any = false;
  for (; i < mant.size(); ++i) {
    const char c = mant[i];
    if (c == '.') {
      if (seen_dot) return false;
      seen_dot = true;
      continue;
    }
    if (c < '0' || c > '9' || m > kGuard) return false;
    m = m * 10 + (c - '0');
    if (seen_dot) ++frac;
    any = true;
  }
  if (!any) return false;

  size_t j = 0;
  long long rad = 0;
  bool have_rad = false;
  while (j < rest.size() && rest[j] >= '0' && rest[j] <= '9') {
    if (rad > kGuard) return false;
    rad = rad * 10 + (rest[j] - '0');
    have_rad = true;
    ++j;
  }
  char dir = 0;
  if (j < rest.size() && (rest[j] == 'u' || rest[j] == 'd' || rest[j] == 'U' || rest[j] == 'D')) {
    dir = static_cast<char>(rest[j] | 0x20);  // to lower
    ++j;
  }
  int exp = 0;
  if (j < rest.size() && (rest[j] == 'e' || rest[j] == 'E')) {
    ++j;
    int esign = 1;
    if (j < rest.size() && (rest[j] == '+' || rest[j] == '-')) {
      esign = (rest[j] == '-') ? -1 : 1;
      ++j;
    }
    if (j >= rest.size() || rest[j] < '0' || rest[j] > '9') return false;
    long e = 0;
    while (j < rest.size() && rest[j] >= '0' && rest[j] <= '9') {
      e = e * 10 + (rest[j] - '0');
      if (e > 1000000) e = 1000000;
      ++j;
    }
    exp = static_cast<int>(esign * e);
  }
  if (j != rest.size()) return false;  // trailing garbage

  // A half-ulp default radius is handled by scaling the whole thing by 10 so the
  // radius (5) and midpoint stay integers; the unit's decimal exponent absorbs it.
  const long long scale = have_rad ? 1 : 10;
  const long long r = have_rad ? rad : 5;
  const long long mid = sign * m * scale;
  const int unit_exp = -frac - (have_rad ? 0 : 1) + exp;
  const long long lo_val = (dir == 'u') ? mid : mid - r;
  const long long hi_val = (dir == 'd') ? mid : mid + r;

  const ParsedReal plo = parse_real(std::to_string(lo_val) + "e" + std::to_string(unit_exp));
  const ParsedReal phi = parse_real(std::to_string(hi_val) + "e" + std::to_string(unit_exp));
  double lrd, lru, urd, uru;
  if (!bracket(plo, lrd, lru) || !bracket(phi, urd, uru)) return false;
  lo = lrd;  // lower endpoint rounds down (outward)
  hi = uru;  // upper endpoint rounds up (outward)
  return true;
}

// --- Formatting -------------------------------------------------------------

inline std::string format_double(double x) {
  if (x == detail::kInf) return "inf";
  if (x == -detail::kInf) return "-inf";
  if (x == 0.0) return "0";
#if PYINTVAL_HAVE_TO_CHARS
  char buf[64];
  auto res = std::to_chars(buf, buf + sizeof(buf), x);  // shortest round-trip
  return std::string(buf, res.ptr);
#else
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.17g", x);
  return std::string(buf);
#endif
}

}  // namespace textdetail

// Parse text into an interval. Accepts (case-insensitive, whitespace-tolerant):
//   "[empty]" / "[]"          -> empty
//   "[entire]"                -> entire
//   "[l, u]"                  -> [RD(l), RU(u)]
//   "[x]" or bare "x"         -> [RD(x), RU(x)]  (a rigorously rounded point)
// Returns false (and sets out to empty) on any malformed input, on NaN, or on
// bounds that cannot form a valid interval (e.g. "[inf,inf]", "[2,1]").
inline bool text_to_interval(std::string_view s, Interval& out) {
  using namespace textdetail;
  out = empty();
  std::string_view t = trim(s);
  if (t.empty()) return false;

  if (t.front() == '[') {
    if (t.back() != ']') return false;
    std::string_view inner = trim(t.substr(1, t.size() - 2));
    if (inner.empty() || iequal(inner, "empty")) {
      out = empty();
      return true;
    }
    if (iequal(inner, "entire")) {
      out = entire();
      return true;
    }
    const size_t comma = inner.find(',');
    if (comma == std::string_view::npos) {
      const ParsedReal p = parse_real(trim(inner));
      double rd, ru;
      if (!bracket(p, rd, ru)) return false;
      if (!valid_bounds(rd, ru)) return false;
      out = make(rd, ru);
      return true;
    }
    // Inf-sup form; an empty endpoint is unbounded, so "[,]" is entire and
    // "[1,]" is [1, +inf]. Each present endpoint rounds outward.
    const std::string_view lo_tok = trim(inner.substr(0, comma));
    const std::string_view hi_tok = trim(inner.substr(comma + 1));
    double lo = -detail::kInf, hi = detail::kInf;
    if (!lo_tok.empty()) {
      double rd, ru;
      if (!bracket(parse_real(lo_tok), rd, ru)) return false;
      lo = rd;  // lower endpoint rounds down (outward)
    }
    if (!hi_tok.empty()) {
      double rd, ru;
      if (!bracket(parse_real(hi_tok), rd, ru)) return false;
      hi = ru;  // upper endpoint rounds up (outward)
    }
    if (!valid_bounds(lo, hi)) return false;
    out = make(lo, hi);
    return true;
  }

  // Bare literal: an uncertain form "m?..." or a rigorously rounded point.
  if (t.find('?') != std::string_view::npos) {
    double lo, hi;
    if (!parse_uncertain(t, lo, hi) || !valid_bounds(lo, hi)) return false;
    out = make(lo, hi);
    return true;
  }
  const ParsedReal p = parse_real(t);
  double rd, ru;
  if (!bracket(p, rd, ru)) return false;
  if (!valid_bounds(rd, ru)) return false;
  out = make(rd, ru);
  return true;
}

// Render an interval as text that text_to_interval reproduces exactly. Each
// finite endpoint uses the shortest round-tripping decimal.
inline std::string interval_to_text(const Interval& x) {
  if (is_empty(x)) return "[empty]";
  if (is_entire(x)) return "[entire]";
  return "[" + textdetail::format_double(x.lo) + ", " + textdetail::format_double(x.hi) + "]";
}

}  // namespace pyintval
