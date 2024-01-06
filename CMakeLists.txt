cmake_minimum_required(VERSION 3.18)

SET(MODNAME "mimo")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")



PROJECT(${MODNAME})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE)
    find_package(OpenBLAS PATHS /opt/homebrew/opt/openblas/lib/cmake/OpenBLAS/) 
    include_directories(${OpenBLAS_INCLUDE_DIRS})
endif()

if(UNIX AND NOT APPLE)
    set(OpenBLAS_DIR /home/niurx/桌面/OpenBLAS-0.3.25/lib/cmake/openblas/) # 去手动编译openblas并指定路径
    find_package(OpenBLAS REQUIRED)
    include_directories(${OpenBLAS_INCLUDE_DIRS})
endif()

if(WIN32)
    set(OpenBLAS_DIR ${PROJECT_SOURCE_DIR}/external/openblas/lib/cmake/openblas/)
    find_package(OpenBLAS REQUIRED)
    include_directories( ${PROJECT_SOURCE_DIR}/external/openblas/include)
endif()

add_subdirectory(external/pybind11)


# 将当前目录及子目录下所有的 .h 文件添加到 HEADER_FILES 变量中
file(GLOB_RECURSE HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")

# 创建一个空的变量来存储包含目录
set(INCLUDE_DIRS "")

# 对于每个头文件
foreach(HEADER ${HEADER_FILES})
    # 获取该头文件的目录
    get_filename_component(DIR ${HEADER} DIRECTORY)
    # 将目录添加到包含目录列表中
    set(INCLUDE_DIRS ${INCLUDE_DIRS} ${DIR})
endforeach()

# 删除重复的目录
list(REMOVE_DUPLICATES INCLUDE_DIRS)

# 将这些目录添加到包含目录中
include_directories(${INCLUDE_DIRS})


file(GLOB_RECURSE CPP_FILES "src/*.cpp")

pybind11_add_module(${MODNAME}
    MODULE
        ${CPP_FILES}
        ${PROJECT_SOURCE_DIR}/binding.cpp
)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(${MODNAME} PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})


add_executable(MimoTest test/test.cpp ${CPP_FILES})

if(APPLE)
    target_link_libraries(${MODNAME} PRIVATE ${OpenBLAS_LIBRARIES})
    target_link_libraries(MimoTest PRIVATE ${OpenBLAS_LIBRARIES})
endif()

if(UNIX AND NOT APPLE)
    target_link_libraries(${MODNAME} PRIVATE ${OpenBLAS_LIBRARIES})
    target_link_libraries(MimoTest PRIVATE ${OpenBLAS_LIBRARIES})
endif()

if(WIN32)
    target_link_libraries(${MODNAME} PRIVATE ${PROJECT_SOURCE_DIR}/external/openblas/bin/libopenblas.dll)
    target_link_libraries(MimoTest PRIVATE ${PROJECT_SOURCE_DIR}/external/openblas/bin/libopenblas.dll)
endif()





