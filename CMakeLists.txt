cmake_minimum_required(VERSION 3.16)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(modularFFT
    VERSION 1.0.0
    DESCRIPTION "Header-only modular fast Fourier transform (NTT) over Z/p"
    LANGUAGES CXX)

# Header-only interface library; GCC or Clang on 64-bit is required
# (__int128). Consumers can add this directory with add_subdirectory().
add_library(modularfft INTERFACE)
target_compile_features(modularfft INTERFACE cxx_std_17)
target_include_directories(modularfft INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

option(MODULARFFT_BUILD_EXAMPLES "Build the example" ON)
option(MODULARFFT_BUILD_TESTS "Build the tests" ON)
option(MODULARFFT_BUILD_BENCHMARKS "Build the benchmark" ON)

set(MODULARFFT_WARNING_FLAGS -Wall -Wextra)

if(MODULARFFT_BUILD_EXAMPLES)
    add_executable(multiply examples/multiply.cpp)
    target_link_libraries(multiply PRIVATE modularfft)
    target_compile_options(multiply PRIVATE ${MODULARFFT_WARNING_FLAGS})
endif()

if(MODULARFFT_BUILD_TESTS)
    enable_testing()
    add_executable(test_modularfft tests/test_modularfft.cpp)
    target_link_libraries(test_modularfft PRIVATE modularfft)
    target_compile_options(test_modularfft PRIVATE ${MODULARFFT_WARNING_FLAGS})
    add_test(NAME unit_tests COMMAND test_modularfft)
endif()

if(MODULARFFT_BUILD_BENCHMARKS)
    add_executable(bench_multiply benchmarks/bench_multiply.cpp)
    target_link_libraries(bench_multiply PRIVATE modularfft)
    target_compile_options(bench_multiply PRIVATE ${MODULARFFT_WARNING_FLAGS})
endif()
