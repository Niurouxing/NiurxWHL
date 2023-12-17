cmake_minimum_required(VERSION 3.4...3.18)
project(mimo)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_RELEASE_FLAGS "-O3")

add_subdirectory(pybind11)
include_directories(${PROJECT_SOURCE_DIR}/include)
pybind11_add_module(mimo src/bindings.cpp)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(mimo PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
