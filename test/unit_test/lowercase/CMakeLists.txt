set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_lowercase_exe
  test_lowercase.f90
)
target_link_libraries(test_lowercase_exe dcaspt2_module)
set_target_properties(test_lowercase_exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
