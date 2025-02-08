set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_uppercase_exe
  test_uppercase.f90
)
target_link_libraries(test_uppercase_exe dcaspt2_module)
set_target_properties(test_uppercase_exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
