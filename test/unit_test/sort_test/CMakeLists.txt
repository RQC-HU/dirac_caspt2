cmake_minimum_required(VERSION 3.7)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_sort_int_exe
    ${CMAKE_SOURCE_DIR}/src/module_sort_swap.f90
    test_sort_int.f90
)
add_executable(test_sort_int_reverse_exe
    ${CMAKE_SOURCE_DIR}/src/module_sort_swap.f90
    test_sort_int_reverse.f90
)
add_executable(test_sort_real_exe
    ${CMAKE_SOURCE_DIR}/src/module_sort_swap.f90
    test_sort_real.f90
)
add_executable(test_sort_real_reverse_exe
    ${CMAKE_SOURCE_DIR}/src/module_sort_swap.f90
    test_sort_real_reverse.f90
)
