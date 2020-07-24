cmake_minimum_required(VERSION 3.12)
project(StaticBTESolver VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake-modules" "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

option(USE_GPU "build a GPU version" off)
option(BULD_TESTING "build tests folder" on)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
    include(CTest)
endif()

add_subdirectory(src)

if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING) AND BUILD_TESTING)
    add_subdirectory(tests)
endif()