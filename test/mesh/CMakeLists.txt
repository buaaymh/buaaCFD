add_executable(test_mesh_basis basis.cpp)
set_target_properties(test_mesh_basis PROPERTIES OUTPUT_NAME basis)
add_test(NAME TestMeshBasis COMMAND basis)

# add_executable(test_mesh_dim2 dim2.cpp)
# set_target_properties(test_mesh_dim2 PROPERTIES OUTPUT_NAME dim2)
# add_test(NAME TestMeshDim2 COMMAND dim2)

# add_executable(test_mesh_edge edge.cpp)
# set_target_properties(test_mesh_edge PROPERTIES OUTPUT_NAME edge)
# add_test(NAME TestMeshEdge COMMAND edge)

# if (${PROJECT_NAME}_ENABLE_VTK)
#   link_libraries(${VTK_LIBRARIES})
#   add_executable(test_mesh_vtk vtk.cpp)
#   set_target_properties(test_mesh_vtk PROPERTIES OUTPUT_NAME vtk)
#   add_test(NAME TestMeshVtk COMMAND vtk)
# endif (${PROJECT_NAME}_ENABLE_VTK)