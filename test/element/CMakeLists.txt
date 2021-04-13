add_executable(test_element_point point.cpp)
set_target_properties(test_element_point PROPERTIES OUTPUT_NAME point)
add_test(NAME TestElementPoint COMMAND point)

add_executable(test_element_node node.cpp)
set_target_properties(test_element_node PROPERTIES OUTPUT_NAME node)
add_test(NAME TestElementNode COMMAND node)

add_executable(test_element_edge edge.cpp)
set_target_properties(test_element_edge PROPERTIES OUTPUT_NAME edge)
add_test(NAME TestElementEdge COMMAND edge)

add_executable(test_element_triangle triangle.cpp)
set_target_properties(test_element_triangle PROPERTIES OUTPUT_NAME triangle)
add_test(NAME TestElementTriangle COMMAND triangle)

