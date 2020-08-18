cmake_minimum_required(VERSION 3.12)
project(StaticBTESolver VERSION 1.0 LANGUAGES CXX C Fortran)

# for cpu profiling
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake-modules" "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

option(USE_GPU "compile as GPU version" off)
option(USE_GMSH "use gmsh to generate mesh" off)
option(BUILD_TESTING "build tests folder" off)
option(USE_TIME "time each process" on)

if(DEFINED CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
else()
    SET (CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
    include(CTest)
endif()

add_subdirectory(src)

if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING) AND BUILD_TESTING)
    add_subdirectory(tests)
endif()