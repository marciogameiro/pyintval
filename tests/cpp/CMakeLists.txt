add_executable(pyintval_cpp_tests test_sanity.cpp)
target_include_directories(pyintval_cpp_tests PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/third_party)
target_compile_features(pyintval_cpp_tests PRIVATE cxx_std_20)

add_test(NAME cpp_unit COMMAND pyintval_cpp_tests)
