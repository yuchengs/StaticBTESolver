find_package(Gmsh REQUIRED)

if (USE_GPU)
    find_package(MPI REQUIRED)
    find_package(CUDAToolkit REQUIRED)
    enable_language(CUDA)
    set_source_files_properties(StaticBTESolver.cpp PROPERTIES LANGUAGE CUDA)
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
                                   PRIVATE "${Gmsh_INCLUDE_DIRS}")
target_link_libraries(BTEGeometry PUBLIC BTEMesh
                              PRIVATE "${Gmsh_LIBRARIES}")
target_compile_features(BTEGeometry PUBLIC cxx_std_14)

add_library(StaticBTESolver StaticBTESolver.cpp "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESolver/StaticBTESolver.h")
if (USE_GPU)
    target_include_directories(StaticBTESolver PUBLIC "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESolver"
                                               PRIVATE "${StaticBTESolver_SOURCE_DIR}/includes/viennacl")
    set_target_properties(StaticBTESolver
            PROPERTIES
            CUDA_SEPARABLE_COMPILATION on)
    set_target_properties(StaticBTESolver
            PROPERTIES
            COMPILE_DEFINITIONS "USE_GPU")
    target_link_libraries(StaticBTESolver PUBLIC BTEBand BTEBoundaryCondition BTEMesh BTEGeometry
            PRIVATE MPI::MPI_CXX CUDA::cudart)
else()
    target_include_directories(StaticBTESolver PUBLIC "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESolver"
                                               PRIVATE "${PETSC_INCLUDES}")
    target_link_libraries(StaticBTESolver PUBLIC BTEBand BTEBoundaryCondition BTEMesh BTEGeometry
                                          PRIVATE "${PETSC_LIBRARIES}")
endif()

add_executable(BTEcmd
        BTEMain.cpp)
if (USE_GPU)
    target_link_libraries(BTEcmd PUBLIC StaticBTESolver MPI::MPI_CXX)
    set_target_properties(BTEcmd PROPERTIES COMPILE_DEFINITIONS "USE_GPU")
else()
    target_link_libraries(BTEcmd PUBLIC StaticBTESolver)
endif()

