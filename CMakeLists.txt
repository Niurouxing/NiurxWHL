CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

SET(MODNAME "mimo")

PROJECT(${MODNAME})

set(BUILD_SHARED_LIBS OFF CACHE BOOL "build shared library" FORCE)

ADD_SUBDIRECTORY(external/pybind11)
ADD_SUBDIRECTORY(external/armadillo-code)


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

# ##############################################################################
#                                  EXECUTABLE                                  #
# ##############################################################################
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

# ##############################################################################
#                                   INSTALL                                    #
# ##############################################################################
INSTALL(TARGETS ${MODNAME} DESTINATION .)
#FILE(GLOB PY_EXAMPLE_FILES "${PROJECT_SOURCE_DIR}/examples/*.py")
#INSTALL(FILES ${PY_EXAMPLE_FILES} DESTINATION examples)

# ##############################################################################
#                                  EXAMPLES                                    #
# ##############################################################################
#ADD_TEST(NAME example
#        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#        COMMAND ${PYTHON_EXECUTABLE} carma_examples.py)
#SET_PROPERTY(TEST example
#        PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")
