file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${BTESolver_SOURCE_DIR}/includes/BTESolver/*.h")
add_library(BTESolver
            BTEMesh.cpp
            BTEBand.cpp
            BTEBoundaryCondition.cpp
        StaticBTESolver.cpp
            BTEGeometry.cpp
            utility.cpp
            ${HEADER_LIST})

target_include_directories(BTESolver PUBLIC "../includes" "${PETSC_INCLUDES}" "${GMSH_INCLUDE_DIRS}")
target_link_libraries(BTESolver "${PETSC_LIBRARIES}" "${GMSH_LIBRARIES}")
target_compile_features(BTESolver PUBLIC cxx_std_14)

add_executable(BTEcmd
               BTEMain.cpp
               BTEMesh.cpp
               BTEBand.cpp
               BTEBoundaryCondition.cpp
        StaticBTESolver.cpp
               BTEGeometry.cpp
               utility.cpp
               ${HEADER_LIST})
target_include_directories(BTEcmd PUBLIC "../includes" "${PETSC_INCLUDES}" "${GMSH_INCLUDE_DIRS}")
target_link_libraries(BTEcmd "${PETSC_LIBRARIES}" "${GMSH_LIBRARIES}")


