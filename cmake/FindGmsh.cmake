
# Try to find gmsh
#
# Once done this will define
#
#  GMSH_FOUND
#  GMSH_INCLUDES
#  GMSH_LIBRARIES
#
#  Usage:
#  find_package(gmsh)
#

cmake_policy(VERSION 3.10)

find_path(GMSH_INCLUDE_DIR
        NAMES gmsh.h gmsh/gmsh.h Gmsh/gmsh.h
        HINTS /usr/local/include /usr/include /opt/include)

find_library(GMSH_LIBRARY
        NAMES gmsh libgmsh libgmsh.so libgmsh.a
        HINTS /usr/local/lib /usr/lib /opt/lib)

mark_as_advanced(GMSH_INCLUDE_DIR GMSH_LIBRARY)

set(GMSH_LIBRARIES ${GMSH_LIBRARY})
set(GMSH_INCLUDE_DIRS ${GMSH_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(gmsh
        REQUIRED_VARS GMSH_LIBRARIES GMSH_INCLUDE_DIRS
        FAIL_MESSAGE "gmsh could not be found. Try to resolve this by adding hints to find_path and find_library")