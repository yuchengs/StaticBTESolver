
# Try to find gmsh
#
# Once done this will define
#
#  Gmsh_FOUND
#  Gmsh_INCLUDES
#  Gmsh_LIBRARIES
#
#  Usage:
#  find_package(gmsh)
#

cmake_policy(VERSION 3.10)

find_path(Gmsh_INCLUDE_DIR
        NAMES gmsh.h gmsh/gmsh.h Gmsh/gmsh.h
        HINTS /usr/local/include /usr/include /opt/include)

find_library(Gmsh_LIBRARY
        NAMES gmsh libgmsh libgmsh.so libgmsh.a
        HINTS /usr/local/lib /usr/lib /opt/lib)

mark_as_advanced(Gmsh_INCLUDE_DIR Gmsh_LIBRARY)

set(Gmsh_LIBRARIES ${Gmsh_LIBRARY})
set(Gmsh_INCLUDE_DIRS ${Gmsh_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Gmsh
        REQUIRED_VARS Gmsh_LIBRARIES Gmsh_INCLUDE_DIRS
        FAIL_MESSAGE "gmsh could not be found. Try to resolve this by adding hints to find_path and find_library")