cmake_minimum_required(VERSION 3.26)
project(pyintval LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# The vendored CORE-MATH transcendental kernels use __builtin_*, __int128 and
# (on x86-64) <x86intrin.h>. Those are provided by Clang and GCC on every
# platform, but NOT by MSVC. On Windows the extension must therefore be built
# with clang-cl. Fail early with an actionable message rather than deep in a
# CORE-MATH source file.
if(MSVC AND NOT CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "GNU"
   AND NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  message(FATAL_ERROR
    "pyintval's correctly-rounded transcendental kernels (CORE-MATH) require "
    "Clang/GCC features unavailable in MSVC. On Windows, build with clang-cl "
    "(e.g. set CMAKE_ARGS=\"-T ClangCL\" or CC/CXX to clang-cl).")
endif()

# Two build modes share this file:
#  - via scikit-build-core (pip install): build the Python extension, skip C++ tests
#  - standalone CMake (development):      build the C++ tests, skip the extension
if(DEFINED SKBUILD)
  set(_pyintval_bindings_default ON)
  set(_pyintval_tests_default OFF)
else()
  set(_pyintval_bindings_default OFF)
  set(_pyintval_tests_default ON)
endif()
option(PYINTVAL_BUILD_BINDINGS "Build the pybind11 extension module" ${_pyintval_bindings_default})
option(PYINTVAL_BUILD_TESTS "Build the C++ unit tests" ${_pyintval_tests_default})

# Vendored CORE-MATH kernels. Each <fn>.c is self-contained: its #include "..."
# directives resolve relative to its own directory, so no include paths are
# needed (and must not be added -- the several dint.h copies are not identical).
file(GLOB PYINTVAL_COREMATH_SOURCES CONFIGURE_DEPENDS
     ${PROJECT_SOURCE_DIR}/third_party/core-math/*/*.c)

# Correct rounding assumes standard IEEE semantics: never fast-math.
if(NOT MSVC OR CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "GNU")
  set_source_files_properties(${PYINTVAL_COREMATH_SOURCES} PROPERTIES
      COMPILE_OPTIONS "-fno-fast-math")
endif()

if(PYINTVAL_BUILD_BINDINGS)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG REQUIRED)

  pybind11_add_module(_core src/bindings/module.cpp ${PYINTVAL_COREMATH_SOURCES})
  target_include_directories(_core PRIVATE ${PROJECT_SOURCE_DIR}/include)
  if(DEFINED SKBUILD_PROJECT_VERSION_FULL)
    target_compile_definitions(_core PRIVATE PYINTVAL_VERSION="${SKBUILD_PROJECT_VERSION_FULL}")
  elseif(DEFINED SKBUILD_PROJECT_VERSION)
    target_compile_definitions(_core PRIVATE PYINTVAL_VERSION="${SKBUILD_PROJECT_VERSION}")
  endif()
  install(TARGETS _core LIBRARY DESTINATION pyintval)
endif()

if(PYINTVAL_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests/cpp)
endif()
