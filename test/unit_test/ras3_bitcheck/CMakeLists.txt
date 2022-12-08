set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ras3_bitcheck_exe
  test_ras3_bitcheck.f90
)
target_link_libraries(ras3_bitcheck_exe dcaspt2_module)
