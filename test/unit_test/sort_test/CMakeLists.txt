cmake_minimum_required(VERSION 3.7)

message(STATUS "CMAKE_SOURCE_DIR        =${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_BINARY_DIR        =${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(sort_test_int
    ${CMAKE_SOURCE_DIR}/src/module_sort.f90
    sort_test_int.f90
)
add_executable(sort_test_real
    ${CMAKE_SOURCE_DIR}/src/module_sort.f90
    sort_test_real.f90
)
