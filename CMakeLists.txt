cmake_minimum_required(VERSION 3.4...3.18)

SET(MODNAME "mimo")

PROJECT(${MODNAME})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_RELEASE_FLAGS "-O3")

add_subdirectory(pybind11)
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB CPP_FILES "src/*.cpp")

pybind11_add_module(${MODNAME}
    MODULE
        ${CPP_FILES}
        ${PROJECT_SOURCE_DIR}/binding.cpp
)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(${MODNAME} PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})


# -----------------------------------------------------------
# cmake_minimum_required(VERSION 3.10)
# project(MimoTest)

# set(CMAKE_CXX_STANDARD 17)

# if(APPLE)
#     link_directories("/opt/homebrew/lib")
#     include_directories("/opt/homebrew/include")
# endif()

# # include the headers files in ./include/
# include_directories(${PROJECT_SOURCE_DIR}/include)

# # include all the cpp files in ./src except the binding.cpp
# file(GLOB CPP_FILES "src/*.cpp")

# add_executable(MimoTest test/test.cpp ${CPP_FILES})
# target_link_libraries(MimoTest -larmadillo)