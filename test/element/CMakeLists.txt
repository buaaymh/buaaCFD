add_executable(test_element_node node.cpp)
target_include_directories(test_element_node PRIVATE ${EIGEN_INC})
set_target_properties(test_element_node PROPERTIES OUTPUT_NAME node)
add_test(NAME TestElementNode COMMAND node)

add_executable(test_element_edge edge.cpp)
target_include_directories(test_element_edge PRIVATE ${EIGEN_INC})
set_target_properties(test_element_edge PROPERTIES OUTPUT_NAME edge)
add_test(NAME TestElementEdge COMMAND edge)

add_executable(test_element_triangle triangle.cpp)
target_include_directories(test_element_triangle PRIVATE ${EIGEN_INC})
set_target_properties(test_element_triangle PROPERTIES OUTPUT_NAME triangle)
add_test(NAME TestElementTriangle COMMAND triangle)

