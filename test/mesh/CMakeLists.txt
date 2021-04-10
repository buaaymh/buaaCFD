add_executable(test_mesh_dim2 dim2.cpp)
target_include_directories(test_mesh_dim2 PRIVATE ${EIGEN_INC})
set_target_properties(test_mesh_dim2 PROPERTIES OUTPUT_NAME dim2)
add_test(NAME TestMeshDim2 COMMAND dim2)