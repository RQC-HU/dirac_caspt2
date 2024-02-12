set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_time_calc_exe
  test_time_calc.f90
)
target_link_libraries(test_time_calc_exe dcaspt2_module)
