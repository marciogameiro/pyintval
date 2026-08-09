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

#include "pyintval/interval.hpp"
#include "pyintval/text.hpp"

namespace py = pybind11;
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
  auto binop = [](const std::function<Interval(const Interval&, const Interval&)>& f,
                  bool reflected) {
    return [f, reflected](const Interval& self, py::handle other) -> py::object {
      Interval o;
      if (!promote(other, o)) return not_implemented();
      return py::cast(reflected ? f(o, self) : f(self, o));
    };
  };
  iv.def("__add__", binop(pyintval::add, false));
  iv.def("__radd__", binop(pyintval::add, true));
  iv.def("__sub__", binop(pyintval::sub, false));
  iv.def("__rsub__", binop(pyintval::sub, true));
  iv.def("__mul__", binop(pyintval::mul, false));
  iv.def("__rmul__", binop(pyintval::mul, true));
  iv.def("__truediv__", binop(pyintval::div, false));
  iv.def("__rtruediv__", binop(pyintval::div, true));
  iv.def("__neg__", [](const Interval& x) { return pyintval::neg(x); });
  iv.def("__pos__", [](const Interval& x) { return x; });
  iv.def("__abs__", [](const Interval& x) { return pyintval::abs(x); });
  iv.def(
      "__pow__",
      [](const Interval& a, py::object e, py::object mod) -> py::object {
        if (!mod.is_none()) return not_implemented();
        if (py::isinstance<py::int_>(e)) return py::cast(pow_int(a, py::cast<long>(e)));
        if (py::isinstance<py::float_>(e)) {
          const double d = py::cast<double>(e);
          if (std::isfinite(d) && d == std::floor(d))
            return py::cast(pow_int(a, static_cast<long>(d)));
          return not_implemented();  // real powers arrive with the transcendentals (M4)
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
  m.def("empty", []() { return pyintval::empty(); });
  m.def("entire", []() { return pyintval::entire(); });
  m.def("sqrt", &pyintval::sqrt, py::arg("x"));
  m.def("sqr", &pyintval::sqr, py::arg("x"));
  m.def("abs", &pyintval::abs, py::arg("x"));
  m.def("recip", &pyintval::recip, py::arg("x"));
  m.def("fma", &pyintval::fma, py::arg("x"), py::arg("y"), py::arg("z"),
        "Fused multiply-add: a tight enclosure of x*y + z.");
  m.def(
      "pown", [](const Interval& x, long n) { return pow_int(x, n); }, py::arg("x"), py::arg("n"));
  m.def("floor", &pyintval::floor, py::arg("x"));
  m.def("ceil", &pyintval::ceil, py::arg("x"));
  m.def("trunc", &pyintval::trunc, py::arg("x"));
  m.def("round", &pyintval::round_ties_to_even, py::arg("x"), "Round half to even.");
  m.def("round_ties_to_away", &pyintval::round_ties_to_away, py::arg("x"));
  m.def("sign", &pyintval::sign, py::arg("x"));
  m.def("min", &pyintval::min, py::arg("x"), py::arg("y"));
  m.def("max", &pyintval::max, py::arg("x"), py::arg("y"));
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

  // Rigorous enclosure of pi (nearest double is below pi; widen up by one ulp).
  m.def(
      "pi", []() { return pyintval::make(0x1.921fb54442d18p1, 0x1.921fb54442d19p1); },
      "A tight interval enclosing the mathematical constant pi.");
}
