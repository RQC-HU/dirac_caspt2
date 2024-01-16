set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_replace_exe
  test_replace.f90
)
target_link_libraries(test_replace_exe dcaspt2_module)
