cmake_minimum_required(VERSION 3.7)

message(STATUS "CMAKE_SOURCE_DIR        =${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_BINARY_DIR        =${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ras_input_reader_exe
    ${CMAKE_SOURCE_DIR}/src/four_caspt2_module.f90
    ${CMAKE_SOURCE_DIR}/src/mem.f90
    ${CMAKE_SOURCE_DIR}/src/input_reader.f90
    test_input_reader.f90
)
