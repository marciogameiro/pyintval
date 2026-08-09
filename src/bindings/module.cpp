// pybind11 bindings: expose the C++20 interval kernel to Python as
// pyintval._core. The Python-facing API (operators, properties, module-level
// math functions) is assembled here on top of the header-only kernel; all
// rigor lives in the headers, this file only marshals types.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <limits>
#include <string>

#include "pyintval/decoration.hpp"
#include "pyintval/elementary.hpp"
#include "pyintval/interval.hpp"
#include "pyintval/text.hpp"

namespace py = pybind11;
using pyintval::DecoratedInterval;
using pyintval::Decoration;
using pyintval::Interval;

#ifndef PYINTVAL_VERSION
#define PYINTVAL_VERSION "0.0.0"
#endif

namespace {

py::object not_implemented() {
  return py::reinterpret_borrow<py::object>(py::handle(Py_NotImplemented));
}

// Promote a Python operand to an interval for arithmetic: intervals pass
// through, finite real scalars become points. Anything else -> false, so the
// caller can return NotImplemented and let Python try the reflected operation.
bool promote(py::handle h, Interval& out) {
  if (py::isinstance<Interval>(h)) {
    out = py::cast<Interval>(h);
    return true;
  }
  if (py::isinstance<py::float_>(h) || py::isinstance<py::int_>(h)) {
    const double v = py::cast<double>(h);
    if (!std::isfinite(v)) {
      throw py::value_error("cannot use a non-finite scalar as an interval operand");
    }
    out = pyintval::point(v);
    return true;
  }
  return false;
}

// Exact, bit-lossless textual form of one endpoint (hex float), used by repr.
std::string hex_endpoint(double v) {
  if (v == pyintval::detail::kInf) return "inf";
  if (v == -pyintval::detail::kInf) return "-inf";
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%a", v);
  return std::string(buf);
}

std::string repr(const Interval& x) {
  if (pyintval::is_empty(x)) return "Interval('[empty]')";
  if (pyintval::is_entire(x)) return "Interval('[entire]')";
  return "Interval('[" + hex_endpoint(x.lo) + ", " + hex_endpoint(x.hi) + "]')";
}

Interval make_checked(double lo, double hi) {
  if (!pyintval::valid_bounds(lo, hi)) {
    throw py::value_error("invalid interval bounds: require lo <= hi with lo < +inf and hi > -inf");
  }
  return pyintval::make(lo, hi);
}

Interval from_object(py::handle a, py::handle b) {
  if (b.is_none()) {
    if (py::isinstance<py::str>(a)) {
      Interval out;
      if (!pyintval::text_to_interval(py::cast<std::string>(a), out)) {
        throw py::value_error("malformed interval literal: '" + py::cast<std::string>(a) + "'");
      }
      return out;
    }
    if (py::isinstance<Interval>(a)) return py::cast<Interval>(a);
    const double v = py::cast<double>(a);  // TypeError for non-numbers
    return make_checked(v, v);
  }
  return make_checked(py::cast<double>(a), py::cast<double>(b));
}

Interval pow_int(const Interval& a, long n) {
  if (n < -1000000000L || n > 1000000000L) {
    throw py::value_error("integer exponent out of supported range");
  }
  return pyintval::pown(a, static_cast<int>(n));
}

// Promote a Python operand to a decorated interval: decorated pass through,
// bare intervals and finite scalars are freshly decorated.
bool promote_dec(py::handle h, DecoratedInterval& out) {
  if (py::isinstance<DecoratedInterval>(h)) {
    out = py::cast<DecoratedInterval>(h);
    return true;
  }
  if (py::isinstance<Interval>(h)) {
    out = pyintval::decorate(py::cast<Interval>(h));
    return true;
  }
  if (py::isinstance<py::float_>(h) || py::isinstance<py::int_>(h)) {
    const double v = py::cast<double>(h);
    if (!std::isfinite(v)) {
      throw py::value_error("cannot use a non-finite scalar as an interval operand");
    }
    out = pyintval::decorate(pyintval::point(v));
    return true;
  }
  return false;
}

Decoration decoration_from_str(const std::string& s) {
  if (s == "com") return Decoration::com;
  if (s == "dac") return Decoration::dac;
  if (s == "def") return Decoration::def;
  if (s == "trv") return Decoration::trv;
  if (s == "ill") return Decoration::ill;
  throw py::value_error("unknown decoration '" + s + "' (expected com/dac/def/trv/ill)");
}

DecoratedInterval dec_pow_int(const DecoratedInterval& a, long n) {
  if (n < -1000000000L || n > 1000000000L) {
    throw py::value_error("integer exponent out of supported range");
  }
  return pyintval::pown(a, static_cast<int>(n));
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() =
      "Compiled core of pyintval: rigorous IEEE 1788-2015 set-based interval "
      "arithmetic with correctly rounded double-precision endpoints.";
  m.attr("__version__") = PYINTVAL_VERSION;
  m.attr("KERNEL_ABI_VERSION") = pyintval::kernel_abi_version;

  m.def(
      "build_info",
      []() {
        py::dict info;
        info["ieee754_doubles"] = std::numeric_limits<double>::is_iec559;
#if defined(_MSVC_LANG)
        info["cxx_standard"] = static_cast<long>(_MSVC_LANG);
#else
        info["cxx_standard"] = static_cast<long>(__cplusplus);
#endif
#if defined(__VERSION__)
        info["compiler"] = std::string(__VERSION__);
#elif defined(_MSC_VER)
        info["compiler"] = "MSVC " + std::to_string(_MSC_VER);
#else
        info["compiler"] = "unknown";
#endif
        return info;
      },
      "Return information about how the extension module was compiled.");

  py::class_<Interval> iv(m, "Interval", R"doc(
A closed, connected set of real numbers with double-precision endpoints.

Construct from bounds ``Interval(lo, hi)``, from a single real number
``Interval(x)`` (a degenerate point interval), or from a string
``Interval("[0.1, 0.2]")`` / ``Interval("0.1")`` which is parsed with correct
outward rounding so the result provably encloses the exact decimal value.

Every operation returns an interval guaranteed to contain the true image of
the operation over its inputs. Division by an interval containing zero yields
an unbounded interval (or the empty set) rather than raising.
)doc");

  iv.def(py::init([](py::object a, py::object b) { return from_object(a, b); }), py::arg("lo"),
         py::arg("hi") = py::none());

  // --- Class/static constructors -------------------------------------------
  iv.def_static("empty", []() { return pyintval::empty(); }, "The empty set.");
  iv.def_static("entire", []() { return pyintval::entire(); }, "The whole real line (-inf, +inf).");

  // --- Properties -----------------------------------------------------------
  iv.def_property_readonly(
      "lo", [](const Interval& x) { return x.lo; }, "Lower endpoint (+inf for the empty set).");
  iv.def_property_readonly(
      "hi", [](const Interval& x) { return x.hi; }, "Upper endpoint (-inf for the empty set).");
  iv.def_property_readonly("mid", &pyintval::mid, "A finite midpoint inside the interval.");
  iv.def_property_readonly("rad", &pyintval::rad, "Radius (half width, rounded up).");
  iv.def_property_readonly("wid", &pyintval::wid, "Width, rounded up.");
  iv.def_property_readonly("mag", &pyintval::mag, "Magnitude sup{|x| : x in self}.");
  iv.def_property_readonly("mig", &pyintval::mig, "Mignitude inf{|x| : x in self}.");
  iv.def_property_readonly("is_empty", &pyintval::is_empty);
  iv.def_property_readonly("is_entire", &pyintval::is_entire);
  iv.def_property_readonly("is_singleton", &pyintval::is_singleton,
                           "True if the interval is a single point.");
  iv.def_property_readonly("is_common", &pyintval::is_common,
                           "True if the interval is nonempty and bounded.");
  iv.def_property_readonly(
      "endpoints", [](const Interval& x) { return py::make_tuple(x.lo, x.hi); },
      "The pair (lo, hi).");

  // --- Arithmetic operators -------------------------------------------------
  using IBin = Interval (*)(const Interval&, const Interval&);
  auto binop = [](IBin f, bool reflected) {
    return [f, reflected](const Interval& self, py::handle other) -> py::object {
      Interval o;
      if (!promote(other, o)) return not_implemented();
      return py::cast(reflected ? f(o, self) : f(self, o));
    };
  };
  iv.def("__add__", binop(static_cast<IBin>(pyintval::add), false));
  iv.def("__radd__", binop(static_cast<IBin>(pyintval::add), true));
  iv.def("__sub__", binop(static_cast<IBin>(pyintval::sub), false));
  iv.def("__rsub__", binop(static_cast<IBin>(pyintval::sub), true));
  iv.def("__mul__", binop(static_cast<IBin>(pyintval::mul), false));
  iv.def("__rmul__", binop(static_cast<IBin>(pyintval::mul), true));
  iv.def("__truediv__", binop(static_cast<IBin>(pyintval::div), false));
  iv.def("__rtruediv__", binop(static_cast<IBin>(pyintval::div), true));
  iv.def("__neg__", [](const Interval& x) { return pyintval::neg(x); });
  iv.def("__pos__", [](const Interval& x) { return x; });
  iv.def("__abs__", [](const Interval& x) { return pyintval::abs(x); });
  iv.def(
      "__pow__",
      [](const Interval& a, py::object e, py::object mod) -> py::object {
        if (!mod.is_none()) return not_implemented();
        if (py::isinstance<py::int_>(e)) return py::cast(pow_int(a, py::cast<long>(e)));
        if (py::isinstance<Interval>(e)) return py::cast(pyintval::pow(a, py::cast<Interval>(e)));
        if (py::isinstance<py::float_>(e)) {
          const double d = py::cast<double>(e);
          if (std::isfinite(d) && d == std::floor(d))
            return py::cast(pow_int(a, static_cast<long>(d)));
          return py::cast(pyintval::pow(a, pyintval::point(d)));  // real power via exp/log
        }
        return not_implemented();
      },
      py::arg("exponent"), py::arg("modulo") = py::none());

  // --- Set operations -------------------------------------------------------
  iv.def("__and__",
         [](const Interval& a, const Interval& b) { return pyintval::intersection(a, b); });
  iv.def("__or__",
         [](const Interval& a, const Interval& b) { return pyintval::convex_hull(a, b); });
  iv.def("intersection", &pyintval::intersection, py::arg("other"));
  iv.def("hull", &pyintval::convex_hull, py::arg("other"),
         "Interval hull of the union of self and other.");

  // --- Comparisons and predicates ------------------------------------------
  iv.def("__eq__", [](const Interval& a, py::handle o) -> py::object {
    if (!py::isinstance<Interval>(o)) return not_implemented();
    return py::cast(pyintval::equal(a, py::cast<Interval>(o)));
  });
  iv.def("__ne__", [](const Interval& a, py::handle o) -> py::object {
    if (!py::isinstance<Interval>(o)) return not_implemented();
    return py::cast(!pyintval::equal(a, py::cast<Interval>(o)));
  });
  iv.def("__contains__", [](const Interval& a, py::handle o) -> bool {
    if (py::isinstance<Interval>(o)) return pyintval::subset(py::cast<Interval>(o), a);
    return pyintval::is_member(py::cast<double>(o), a);
  });
  iv.def("__hash__", [](const Interval& x) -> Py_hash_t {
    std::uint64_t a, b;
    std::memcpy(&a, &x.lo, 8);
    std::memcpy(&b, &x.hi, 8);
    const std::uint64_t h = a * 1000003ULL ^ (b + 0x9e3779b97f4a7c15ULL);
    Py_hash_t r = static_cast<Py_hash_t>(h);
    return r == -1 ? -2 : r;
  });
  iv.def(
      "subset", [](const Interval& a, const Interval& b) { return pyintval::subset(a, b); },
      py::arg("other"), "True if self is a subset of other.");
  iv.def(
      "superset", [](const Interval& a, const Interval& b) { return pyintval::subset(b, a); },
      py::arg("other"), "True if self is a superset of other.");
  iv.def(
      "is_interior_to",
      [](const Interval& a, const Interval& b) { return pyintval::interior(a, b); },
      py::arg("other"));
  iv.def("is_disjoint", &pyintval::disjoint, py::arg("other"));
  iv.def(
      "overlaps", [](const Interval& a, const Interval& b) { return !pyintval::disjoint(a, b); },
      py::arg("other"));
  iv.def(
      "contains",
      [](const Interval& a, py::handle o) -> bool {
        if (py::isinstance<Interval>(o)) return pyintval::subset(py::cast<Interval>(o), a);
        return pyintval::is_member(py::cast<double>(o), a);
      },
      py::arg("value"), "True if value (a number or interval) lies in self.");
  iv.def("precedes", &pyintval::precedes, py::arg("other"));
  iv.def("strict_precedes", &pyintval::strict_precedes, py::arg("other"));
  iv.def("less", &pyintval::less, py::arg("other"));
  iv.def("strict_less", &pyintval::strict_less, py::arg("other"));

  // --- Text ----------------------------------------------------------------
  iv.def("__repr__", &repr);
  iv.def("__str__", [](const Interval& x) { return pyintval::interval_to_text(x); });

  // --- Pickling / copy ------------------------------------------------------
  iv.def(py::pickle([](const Interval& x) { return py::make_tuple(x.lo, x.hi); },
                    [](py::tuple t) {
                      if (t.size() != 2) throw std::runtime_error("invalid interval pickle");
                      return Interval{t[0].cast<double>(), t[1].cast<double>()};
                    }));
  iv.def("__copy__", [](const Interval& x) { return x; });
  iv.def("__deepcopy__", [](const Interval& x, py::handle) { return x; }, py::arg("memo"));

  // --- Module-level functions ----------------------------------------------
  // Math functions are polymorphic: an Interval argument yields an Interval, a
  // DecoratedInterval yields a DecoratedInterval (carrying the propagated
  // decoration). These macros bind both overloads under one name.
#define PYINTVAL_FN1(PYNAME, FN)                                                  \
  m.def(PYNAME, [](const Interval& x) { return pyintval::FN(x); }, py::arg("x")); \
  m.def(PYNAME, [](const DecoratedInterval& x) { return pyintval::FN(x); }, py::arg("x"))
#define PYINTVAL_FN2(PYNAME, FN, A, B)                                                             \
  m.def(                                                                                           \
      PYNAME, [](const Interval& a, const Interval& b) { return pyintval::FN(a, b); }, py::arg(A), \
      py::arg(B));                                                                                 \
  m.def(                                                                                           \
      PYNAME,                                                                                      \
      [](const DecoratedInterval& a, const DecoratedInterval& b) { return pyintval::FN(a, b); },   \
      py::arg(A), py::arg(B))

  m.def("empty", []() { return pyintval::empty(); });
  m.def("entire", []() { return pyintval::entire(); });
  PYINTVAL_FN1("sqrt", sqrt);
  PYINTVAL_FN1("sqr", sqr);
  PYINTVAL_FN1("abs", abs);
  PYINTVAL_FN1("recip", recip);
  m.def(
      "fma",
      [](const Interval& x, const Interval& y, const Interval& z) {
        return pyintval::fma(x, y, z);
      },
      py::arg("x"), py::arg("y"), py::arg("z"),
      "Fused multiply-add: a tight enclosure of x*y + z.");
  m.def(
      "fma",
      [](const DecoratedInterval& x, const DecoratedInterval& y, const DecoratedInterval& z) {
        return pyintval::fma(x, y, z);
      },
      py::arg("x"), py::arg("y"), py::arg("z"));
  m.def(
      "pown", [](const Interval& x, long n) { return pow_int(x, n); }, py::arg("x"), py::arg("n"));
  m.def(
      "pown", [](const DecoratedInterval& x, long n) { return dec_pow_int(x, n); }, py::arg("x"),
      py::arg("n"));
  PYINTVAL_FN1("floor", floor);
  PYINTVAL_FN1("ceil", ceil);
  PYINTVAL_FN1("trunc", trunc);
  PYINTVAL_FN1("round", round_ties_to_even);
  PYINTVAL_FN1("round_ties_to_away", round_ties_to_away);
  PYINTVAL_FN1("sign", sign);
  PYINTVAL_FN2("min", min, "x", "y");
  PYINTVAL_FN2("max", max, "x", "y");
  m.def("hull", &pyintval::convex_hull, py::arg("x"), py::arg("y"));
  m.def("intersection", &pyintval::intersection, py::arg("x"), py::arg("y"));
  m.def("cancel_minus", &pyintval::cancel_minus, py::arg("a"), py::arg("b"));
  m.def("cancel_plus", &pyintval::cancel_plus, py::arg("a"), py::arg("b"));
  m.def(
      "mul_rev", [](const Interval& b, const Interval& c) { return pyintval::mul_rev(b, c); },
      py::arg("b"), py::arg("c"), "{x : b*x meets c} as an interval hull.");
  m.def(
      "mul_rev",
      [](const Interval& b, const Interval& c, const Interval& x) {
        return pyintval::mul_rev(b, c, x);
      },
      py::arg("b"), py::arg("c"), py::arg("x"));
  m.def("sqr_rev", [](const Interval& c) { return pyintval::sqr_rev(c); }, py::arg("c"));
  m.def(
      "sqr_rev", [](const Interval& c, const Interval& x) { return pyintval::sqr_rev(c, x); },
      py::arg("c"), py::arg("x"));
  m.def("abs_rev", [](const Interval& c) { return pyintval::abs_rev(c); }, py::arg("c"));
  m.def(
      "abs_rev", [](const Interval& c, const Interval& x) { return pyintval::abs_rev(c, x); },
      py::arg("c"), py::arg("x"));

  // --- Elementary (transcendental) functions -------------------------------
  // Each returns a rigorous enclosure of the true image over its input, built
  // on correctly-rounded CORE-MATH kernels widened one ulp per endpoint.
  PYINTVAL_FN1("exp", exp);
  PYINTVAL_FN1("exp2", exp2);
  PYINTVAL_FN1("exp10", exp10);
  PYINTVAL_FN1("expm1", expm1);
  PYINTVAL_FN1("log", log);
  PYINTVAL_FN1("log2", log2);
  PYINTVAL_FN1("log10", log10);
  PYINTVAL_FN1("log1p", log1p);
  PYINTVAL_FN1("cbrt", cbrt);
  PYINTVAL_FN1("sin", sin);
  PYINTVAL_FN1("cos", cos);
  PYINTVAL_FN1("tan", tan);
  PYINTVAL_FN1("asin", asin);
  PYINTVAL_FN1("acos", acos);
  PYINTVAL_FN1("atan", atan);
  PYINTVAL_FN2("atan2", atan2, "y", "x");
  PYINTVAL_FN1("sinh", sinh);
  PYINTVAL_FN1("cosh", cosh);
  PYINTVAL_FN1("tanh", tanh);
  PYINTVAL_FN1("asinh", asinh);
  PYINTVAL_FN1("acosh", acosh);
  PYINTVAL_FN1("atanh", atanh);
  PYINTVAL_FN2("hypot", hypot, "x", "y");
  PYINTVAL_FN2("pow", pow, "x", "y");
  PYINTVAL_FN1("erf", erf);
  PYINTVAL_FN1("erfc", erfc);
#undef PYINTVAL_FN1
#undef PYINTVAL_FN2

  // Rigorous enclosures of common constants.
  m.def(
      "pi", []() { return pyintval::detail::kPi; },
      "A tight interval enclosing the mathematical constant pi.");
  m.def(
      "e", []() { return pyintval::exp(pyintval::point(1.0)); },
      "A tight interval enclosing Euler's number e.");

  // --- DecoratedInterval ----------------------------------------------------
  py::class_<DecoratedInterval> di(m, "DecoratedInterval", R"doc(
An interval paired with an IEEE 1788 decoration certifying, from the
computation's history, the strongest property known of the function evaluated
so far: 'com' (defined, continuous, bounded on a common input), 'dac' (defined
and continuous), 'def' (defined), 'trv' (only the enclosure is guaranteed), or
'ill' (not an interval). A result decoration of 'dac' or 'com' is a
machine-checked certificate that the whole composed expression is defined and
continuous on its input box.
)doc");

  di.def(py::init([](py::object a, py::object b) { return pyintval::decorate(from_object(a, b)); }),
         py::arg("lo"), py::arg("hi") = py::none());
  di.def_static(
      "from_parts",
      [](const Interval& x, const std::string& d) {
        return pyintval::decorate(x, decoration_from_str(d));
      },
      py::arg("interval"), py::arg("decoration"),
      "Build a decorated interval from a bare interval and an explicit decoration.");
  di.def_static("nai", []() { return pyintval::nai(); }, "The Not-an-Interval (ill) value.");

  di.def_property_readonly(
      "interval", [](const DecoratedInterval& d) { return d.x; }, "The bare interval.");
  di.def_property_readonly(
      "decoration", [](const DecoratedInterval& d) { return std::string(decoration_name(d.dec)); },
      "The decoration as one of 'com'/'dac'/'def'/'trv'/'ill'.");
  di.def_property_readonly("lo", [](const DecoratedInterval& d) { return d.x.lo; });
  di.def_property_readonly("hi", [](const DecoratedInterval& d) { return d.x.hi; });
  di.def_property_readonly("is_common",
                           [](const DecoratedInterval& d) { return d.dec == Decoration::com; });
  di.def_property_readonly("is_defined_and_continuous", [](const DecoratedInterval& d) {
    return static_cast<int>(d.dec) >= static_cast<int>(Decoration::dac);
  });
  di.def_property_readonly("is_defined", [](const DecoratedInterval& d) {
    return static_cast<int>(d.dec) >= static_cast<int>(Decoration::def);
  });
  di.def_property_readonly("is_nai",
                           [](const DecoratedInterval& d) { return d.dec == Decoration::ill; });

  using DBin = DecoratedInterval (*)(const DecoratedInterval&, const DecoratedInterval&);
  auto dbinop = [](DBin f, bool reflected) {
    return [f, reflected](const DecoratedInterval& self, py::handle other) -> py::object {
      DecoratedInterval o;
      if (!promote_dec(other, o)) return not_implemented();
      return py::cast(reflected ? f(o, self) : f(self, o));
    };
  };
  di.def("__add__", dbinop(static_cast<DBin>(pyintval::add), false));
  di.def("__radd__", dbinop(static_cast<DBin>(pyintval::add), true));
  di.def("__sub__", dbinop(static_cast<DBin>(pyintval::sub), false));
  di.def("__rsub__", dbinop(static_cast<DBin>(pyintval::sub), true));
  di.def("__mul__", dbinop(static_cast<DBin>(pyintval::mul), false));
  di.def("__rmul__", dbinop(static_cast<DBin>(pyintval::mul), true));
  di.def("__truediv__", dbinop(static_cast<DBin>(pyintval::div), false));
  di.def("__rtruediv__", dbinop(static_cast<DBin>(pyintval::div), true));
  di.def("__neg__", [](const DecoratedInterval& x) { return pyintval::neg(x); });
  di.def("__pos__", [](const DecoratedInterval& x) { return x; });
  di.def("__abs__", [](const DecoratedInterval& x) { return pyintval::abs(x); });
  di.def(
      "__pow__",
      [](const DecoratedInterval& a, py::object e, py::object mod) -> py::object {
        if (!mod.is_none()) return not_implemented();
        if (py::isinstance<py::int_>(e)) return py::cast(dec_pow_int(a, py::cast<long>(e)));
        if (py::isinstance<DecoratedInterval>(e))
          return py::cast(pyintval::pow(a, py::cast<DecoratedInterval>(e)));
        if (py::isinstance<py::float_>(e)) {
          const double d = py::cast<double>(e);
          if (std::isfinite(d) && d == std::floor(d))
            return py::cast(dec_pow_int(a, static_cast<long>(d)));
          return py::cast(pyintval::pow(a, pyintval::decorate(pyintval::point(d))));
        }
        return not_implemented();
      },
      py::arg("exponent"), py::arg("modulo") = py::none());

  di.def("__repr__", [](const DecoratedInterval& d) {
    return "DecoratedInterval.from_parts(Interval('" +
           (pyintval::is_empty(d.x) ? std::string("[empty]") : pyintval::interval_to_text(d.x)) +
           "'), '" + decoration_name(d.dec) + "')";
  });
  di.def("__str__", [](const DecoratedInterval& d) {
    return pyintval::interval_to_text(d.x) + "_" + decoration_name(d.dec);
  });
  di.def("__eq__", [](const DecoratedInterval& a, py::handle o) -> py::object {
    if (!py::isinstance<DecoratedInterval>(o)) return not_implemented();
    const auto& b = py::cast<DecoratedInterval>(o);
    return py::cast(pyintval::equal(a.x, b.x) && a.dec == b.dec);
  });
  di.def(py::pickle(
      [](const DecoratedInterval& d) {
        return py::make_tuple(d.x.lo, d.x.hi, static_cast<int>(d.dec));
      },
      [](py::tuple t) {
        if (t.size() != 3) throw std::runtime_error("invalid decorated-interval pickle");
        return DecoratedInterval{Interval{t[0].cast<double>(), t[1].cast<double>()},
                                 static_cast<Decoration>(t[2].cast<int>())};
      }));
}
