set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ras3_bitcheck_exe
  ${CMAKE_SOURCE_DIR}/src/four_caspt2_module.f90
  ${CMAKE_SOURCE_DIR}/src/mem.f90
  ${CMAKE_SOURCE_DIR}/src/module_file_manager.f90
  ${CMAKE_SOURCE_DIR}/src/module_sort_swap.f90
  ${CMAKE_SOURCE_DIR}/src/read_input_module.f90
  ${CMAKE_SOURCE_DIR}/src/ras_det_check.f90
  test_ras3_bitcheck.f90
)
