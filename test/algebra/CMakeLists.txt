add_executable(test_algebra_column column.cpp)
target_include_directories(test_algebra_column PRIVATE ${EIGEN_INC})
set_target_properties(test_algebra_column PROPERTIES OUTPUT_NAME column)
add_test(NAME TestAlgebraColumn COMMAND column)
