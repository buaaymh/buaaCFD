add_executable(test_algebra_column column.cpp)
set_target_properties(test_algebra_column PROPERTIES OUTPUT_NAME column)
add_test(NAME TestAlgebraColumn COMMAND column)
