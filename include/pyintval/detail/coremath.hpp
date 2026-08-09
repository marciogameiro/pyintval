#pragma once

// Declarations of the vendored CORE-MATH correctly-rounded binary64 entry
// points (see third_party/core-math/). Each cr_* function returns the
// round-to-nearest correctly-rounded value of the corresponding elementary
// function -- i.e. the true result is within half an ulp of the return value.
// The interval layer (elementary.hpp) turns these into rigorous enclosures by
// widening one ulp per endpoint.
//
// These are C functions; declare them with C linkage so the C++ interval code
// can call them. The definitions are compiled from the vendored .c sources and
// linked into the extension.

extern "C" {
double cr_exp(double);
double cr_exp2(double);
double cr_exp10(double);
double cr_expm1(double);
double cr_log(double);
double cr_log2(double);
double cr_log10(double);
double cr_log1p(double);
double cr_cbrt(double);
double cr_sin(double);
double cr_cos(double);
double cr_tan(double);
double cr_asin(double);
double cr_acos(double);
double cr_atan(double);
double cr_atan2(double, double);  // (y, x)
double cr_sinh(double);
double cr_cosh(double);
double cr_tanh(double);
double cr_asinh(double);
double cr_acosh(double);
double cr_atanh(double);
double cr_hypot(double, double);
double cr_pow(double, double);
double cr_erf(double);
double cr_erfc(double);
}
