cmake_minimum_required(VERSION 3.26)
project(pyintval LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

if(PYINTVAL_BUILD_BINDINGS)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG REQUIRED)

  pybind11_add_module(_core src/bindings/module.cpp)
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
