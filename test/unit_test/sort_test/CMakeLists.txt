set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_sort_int_exe
  test_sort_int.f90
)
target_link_libraries(test_sort_int_exe dcaspt2_module)
add_executable(test_sort_int_reverse_exe
  test_sort_int_reverse.f90
)
target_link_libraries(test_sort_int_reverse_exe dcaspt2_module)
add_executable(test_sort_real_exe
  test_sort_real.f90
)
target_link_libraries(test_sort_real_exe dcaspt2_module)
add_executable(test_sort_real_reverse_exe
  test_sort_real_reverse.f90
)
target_link_libraries(test_sort_real_reverse_exe dcaspt2_module)

set_target_properties(test_sort_int_exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set_target_properties(test_sort_int_reverse_exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set_target_properties(test_sort_real_exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set_target_properties(test_sort_real_reverse_exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
