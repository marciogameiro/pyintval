# The elementary-function tests call the interval transcendentals, which are
# backed by the vendored CORE-MATH C kernels; link them into the test binary
# (${PYINTVAL_COREMATH_SOURCES} is defined in the top-level CMakeLists).
add_executable(pyintval_cpp_tests
    test_sanity.cpp
    test_rounding.cpp
    test_rounding_hw.cpp
    test_interval_arith.cpp
    test_interval_sets.cpp
    test_text.cpp
    test_elementary.cpp
    test_decoration.cpp
    ${PYINTVAL_COREMATH_SOURCES})
target_include_directories(pyintval_cpp_tests PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/third_party)
target_compile_features(pyintval_cpp_tests PRIVATE cxx_std_20)

# test_rounding_hw.cpp flips the hardware rounding mode (fesetround) to
# cross-validate the EFT primitives bit-for-bit; that translation unit alone
# must be compiled with strict floating-point semantics so the compiler does
# not fold or reorder FP operations across the mode switches.
if(MSVC)
  set_source_files_properties(test_rounding_hw.cpp PROPERTIES COMPILE_OPTIONS "/fp:strict")
else()
  set_source_files_properties(test_rounding_hw.cpp PROPERTIES COMPILE_OPTIONS "-frounding-math")
endif()

add_test(NAME cpp_unit COMMAND pyintval_cpp_tests)
