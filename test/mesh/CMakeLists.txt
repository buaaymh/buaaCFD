add_executable(test_mesh_point point.cpp)
target_include_directories(test_mesh_point PRIVATE ${EIGEN_INC})
set_target_properties(test_mesh_point PROPERTIES OUTPUT_NAME point)
add_test(NAME TestMeshPoint COMMAND point)

# add_executable(test_mesh_edge edge.cpp)
# target_include_directories(test_mesh_edge PRIVATE ${EIGEN_INC})
# set_target_properties(test_mesh_edge PROPERTIES OUTPUT_NAME edge)
# add_test(NAME TestMeshEdge COMMAND edge)

# add_executable(test_mesh_triangle triangle.cpp)
# target_include_directories(test_mesh_triangle PRIVATE ${EIGEN_INC})
# set_target_properties(test_mesh_triangle PROPERTIES OUTPUT_NAME triangle)
# add_test(NAME TestMeshTriangle COMMAND triangle)

# add_executable(test_mesh_dim2 dim2.cpp)
# target_include_directories(test_mesh_dim2 PRIVATE ${EIGEN_INC})
# set_target_properties(test_mesh_dim2 PROPERTIES OUTPUT_NAME dim2)
# add_test(NAME TestMeshDim2 COMMAND dim2)
