find_package(Gmsh REQUIRED)

if (USE_GPU)
    include(CheckLanguage)
    check_language(CUDA)
    if (CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
    else()
        message(FATAL_ERROR "No CUDA support")
    endif()
else()
    find_package(PETSc REQUIRED)
endif()


add_library(BTEutility utility.cpp "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESolver/utility.h")
target_include_directories(BTEutility PUBLIC "${StaticBTESolver_SOURCE_DIR}/includes")
target_compile_features(BTEutility PUBLIC cxx_std_14)

add_library(BTEBand BTEBand.cpp "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESolver/BTEBand.h")
target_include_directories(BTEBand PUBLIC "${StaticBTESolver_SOURCE_DIR}/includes")
target_link_libraries(BTEBand PUBLIC BTEutility)
target_compile_features(BTEBand PUBLIC cxx_std_14)

add_library(BTEBoundaryCondition BTEBoundaryCondition.cpp "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESolver/BTEBoundaryCondition.h")
target_include_directories(BTEBoundaryCondition PUBLIC "${StaticBTESolver_SOURCE_DIR}/includes")
target_link_libraries(BTEBoundaryCondition PUBLIC BTEutility)
target_compile_features(BTEBoundaryCondition PUBLIC cxx_std_14)

add_library(BTEMesh BTEMesh.cpp "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESolver/BTEMesh.h")
target_include_directories(BTEMesh PUBLIC "${StaticBTESolver_SOURCE_DIR}/includes")
target_link_libraries(BTEMesh PUBLIC BTEutility)
target_compile_features(BTEMesh PUBLIC cxx_std_14)

add_library(BTEGeometry BTEGeometry.cpp "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESolver/BTEGeometry.h")
target_include_directories(BTEGeometry PUBLIC "${StaticBTESolver_SOURCE_DIR}/includes"
                                   PRIVATE "${GMSH_INCLUDE_DIRS}")
target_link_libraries(BTEGeometry PUBLIC BTEMesh
                              PRIVATE "${GMSH_LIBRARIES}")
target_compile_features(BTEGeometry PUBLIC cxx_std_14)

add_library(StaticBTESolver StaticBTESolver.cpp "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESolver/StaticBTESolver.h")
if (USE_GPU)
    target_include_directories(StaticBTESolver PUBLIC "${StaticBTESolver_SOURCE_DIR}/includes")
    target_link_libraries(StaticBTESolver PUBLIC BTEBand BTEBoundaryCondition BTEMesh BTEGeometry
            PRIVATE CUDA::toolkit CUDA::cusolver CUDA::cusparse CUDA::cublas)
else()
    target_include_directories(StaticBTESolver PUBLIC "${StaticBTESolver_SOURCE_DIR}/includes"
                                               PRIVATE "${PETSC_INCLUDES}")
    target_link_libraries(StaticBTESolver PUBLIC BTEBand BTEBoundaryCondition BTEMesh BTEGeometry
                                          PRIVATE "${PETSC_LIBRARIES}")
endif()

add_executable(BTEcmd
        BTEMain.cpp)
target_link_libraries(BTEcmd PUBLIC StaticBTESolver)

