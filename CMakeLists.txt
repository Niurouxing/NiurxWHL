CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

SET(MODNAME "mimo")

PROJECT(${MODNAME})


set(BUILD_SHARED_LIBS OFF CACHE BOOL "build shared library" FORCE)

ADD_SUBDIRECTORY(external/pybind11)
ADD_SUBDIRECTORY(external/armadillo-code)

# 开启o2优化
set(PROJECT_RELEASE_FLAGS "-O2")
 
if(APPLE)
    link_directories("/opt/homebrew/lib")
    include_directories("/opt/homebrew/include")
endif()

if($ENV{CIBUILDWHEEL})
    target_compile_definitions(armadillo PRIVATE -DARMA_USE_LAPACK)
    target_compile_definitions(armadillo PRIVATE -DARMA_USE_BLAS)
endif()

# Add -fPIC for Armadillo (and OpenBLAS if compiled)
if(NOT MSVC)
    # clang on Windows does not support -fPIC
    if(NOT WIN32)
        target_compile_options(armadillo PRIVATE -fPIC)
    endif()
endif()

ADD_SUBDIRECTORY(external/carma)

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB CPP_FILES "src/*.cpp")

pybind11_add_module(${MODNAME}
    MODULE
        ${CPP_FILES}
)

TARGET_LINK_LIBRARIES(${MODNAME} PUBLIC carma::carma)

TARGET_INCLUDE_DIRECTORIES(${MODNAME}
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

TARGET_COMPILE_OPTIONS(${MODNAME}
    PUBLIC
        "$<$<CONFIG:RELEASE>:${PROJECT_RELEASE_FLAGS}>"
)

TARGET_COMPILE_DEFINITIONS(${MODNAME}
    PUBLIC
        "$<$<CONFIG:RELEASE>:${PROJECT_RELEASE_DEFINITIONS}>"
)

INSTALL(TARGETS ${MODNAME} DESTINATION .)
 













# #############################################################################################################################

# cmake_minimum_required(VERSION 3.10)
# project(MimoTest)

# set(CMAKE_CXX_STANDARD 17)

# if(APPLE)
#     link_directories("/opt/homebrew/lib")
#     include_directories("/opt/homebrew/include")
# endif()

# # include the headers files in ./include/
# include_directories(${PROJECT_SOURCE_DIR}/include)

# # include the utils.cpp in ./src/
# file(GLOB CPP_FILES "src/utils.cpp")


# add_executable(MimoTest test/test.cpp ${CPP_FILES})
# target_link_libraries(MimoTest -larmadillo)

