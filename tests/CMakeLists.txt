file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mesh2D.mphtxt
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mesh3D.mphtxt
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_executable(testBTEMesh testBTEMesh.cpp)
target_compile_features(testBTEMesh PRIVATE cxx_std_14)
target_include_directories(testBTEMesh PUBLIC ../includes)
target_link_libraries(testBTEMesh PRIVATE BTESolver)

add_test(NAME testBTEMesh COMMAND testBTEMesh)


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Input-dispersion-relation-fp.dat
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_executable(testBTEBand testBTEBand.cpp)
target_compile_features(testBTEBand PRIVATE cxx_std_14)
target_include_directories(testBTEBand PUBLIC ../includes)
target_link_libraries(testBTEBand PRIVATE BTESolver)

add_test(NAME testBTEBand COMMAND testBTEBand)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/inputbc2D.dat
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/inputbc3D.dat
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_executable(testBTEBoundaryCondition testBTEBoundaryCondition.cpp)
target_compile_features(testBTEBoundaryCondition PRIVATE cxx_std_14)
target_include_directories(testBTEBoundaryCondition PUBLIC ../includes)
target_link_libraries(testBTEBoundaryCondition PRIVATE BTESolver)

add_test(NAME testBTEBoundaryCondition COMMAND testBTEBoundaryCondition)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/inputdata2D.dat
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(testBTESolver2D testBTESolver2D.cpp)
target_compile_features(testBTESolver2D PRIVATE cxx_std_14)
target_include_directories(testBTESolver2D PUBLIC ../includes)
target_link_libraries(testBTESolver2D PRIVATE BTESolver)

add_test(NAME testBTESolver2D COMMAND testBTESolver2D)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/inputdata3D.dat
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(testBTESolver3D testBTESolver3D.cpp)
target_compile_features(testBTESolver3D PRIVATE cxx_std_14)
target_include_directories(testBTESolver3D PUBLIC ../includes)
target_link_libraries(testBTESolver3D PRIVATE BTESolver)

add_test(NAME testBTESolver3D COMMAND testBTESolver3D)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/cube.geo
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(testBTEGeometry3D testBTEGeometry3D.cpp)
target_compile_features(testBTEGeometry3D PRIVATE cxx_std_14)
target_include_directories(testBTEGeometry3D PUBLIC ../includes)
target_link_libraries(testBTEGeometry3D PRIVATE BTESolver)

add_test(NAME testBTEGeometry3D COMMAND testBTEGeometry3D)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/2DM1DG
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/3DM1DG
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/2DM2DG
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/3DM2DG
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/3DM3DG
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})