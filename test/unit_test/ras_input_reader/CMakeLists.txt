set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_ras_input_reader_exe
  ${CMAKE_SOURCE_DIR}/src/four_caspt2_module.f90
  ${CMAKE_SOURCE_DIR}/src/mem.f90
  ${CMAKE_SOURCE_DIR}/src/module_file_manager.f90
  ${CMAKE_SOURCE_DIR}/src/module_sort_swap.f90
  ${CMAKE_SOURCE_DIR}/src/read_input_module.f90
  test_ras_input_reader.f90
)
