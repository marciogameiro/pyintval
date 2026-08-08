#include <pybind11/pybind11.h>

#include <limits>
#include <string>

#include "pyintval/interval.hpp"

namespace py = pybind11;

#ifndef PYINTVAL_VERSION
#define PYINTVAL_VERSION "0.0.0"
#endif

PYBIND11_MODULE(_core, m) {
  m.doc() =
      "Compiled core of pyintval — rigorous interval arithmetic with correctly "
      "rounded double-precision endpoints (Milestone 1 skeleton).";

  m.attr("__version__") = PYINTVAL_VERSION;
  m.attr("KERNEL_ABI_VERSION") = pyintval::kernel_abi_version;

  m.def(
      "build_info",
      []() {
        py::dict info;
        info["ieee754_doubles"] = std::numeric_limits<double>::is_iec559;
    // MSVC keeps __cplusplus at 199711L unless /Zc:__cplusplus is passed;
    // _MSVC_LANG always carries the actual language standard.
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
}
