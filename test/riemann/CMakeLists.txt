add_executable(test_riemann_types types.cpp)
set_target_properties(test_riemann_types PROPERTIES OUTPUT_NAME types)
add_test(NAME TestRiemannTypes COMMAND types)

add_executable(test_riemann_linear linear.cpp)
set_target_properties(test_riemann_linear PROPERTIES OUTPUT_NAME linear)
add_test(NAME TestRiemannLinear COMMAND linear)

add_executable(test_riemann_ausm ausm.cpp)
set_target_properties(test_riemann_ausm PROPERTIES OUTPUT_NAME ausm)
add_test(NAME TestRiemannAusm COMMAND ausm)

add_executable(test_riemann_rotated rotated.cpp)
set_target_properties(test_riemann_rotated PROPERTIES OUTPUT_NAME rotated)
add_test(NAME TestRiemannRotated COMMAND rotated)
