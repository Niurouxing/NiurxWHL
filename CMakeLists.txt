cmake_minimum_required(VERSION 3.4...3.18)

SET(MODNAME "mimo")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")


PROJECT(${MODNAME})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


 

add_subdirectory(pybind11)
include_directories(${PROJECT_SOURCE_DIR}/include)


file(GLOB_RECURSE CPP_FILES "src/*.cpp")

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
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")

# project(MimoTest)


# set(CMAKE_CXX_STANDARD 17)

# if(APPLE)
#     link_directories("/opt/homebrew/lib")
#     include_directories("/opt/homebrew/include")
# endif()

# # include the headers files in ./include/
# include_directories(${PROJECT_SOURCE_DIR}/include)

# file(GLOB_RECURSE CPP_FILES "src/*.cpp")

# add_executable(MimoTest test/test.cpp ${CPP_FILES})
# target_link_libraries(MimoTest -larmadillo)