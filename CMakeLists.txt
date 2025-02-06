cmake_minimum_required(VERSION 3.14)
project(dirac_caspt2 Fortran)
message(STATUS "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

# Build options
option(MPI "Enable MPI?" OFF)
option(OPENMP "Enable OpenMP?" OFF)
option(MKL "Enable MKL?" ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules) # Set module output path

# Deny In-source build ( Ref : https://github.com/eigenteam/eigen-git-mirror/blob/36b95962756c1fce8e29b1f8bc45967f30773c00/CMakeLists.txt#L7-L9)
if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# Set default Install path if -DCMAKE_INSTALL_PREFIX is not specified
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/usr/local/bin" CACHE PATH "..." FORCE)
endif()

if(MPI)
  find_package(MPI REQUIRED) # Find MPI
  add_compile_options(-DHAVE_MPI) # Add MPI preprocessor flag
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel OR CMAKE_Fortran_COMPILER_ID STREQUAL IntelLLVM)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    set(CMAKE_Fortran_FLAGS "-mp1")
  endif()

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DINTEL -g -traceback -cpp -i8 -I$ENV{MKLROOT}/include -pad -integer-size 64 -unroll -warn nounused -nogen-interface")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -debug extended -debug-parameters -warn -DDEBUG")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  link_libraries(-i8)

  if(MKL)
    link_libraries(-L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl)
  endif()

  if(OPENMP)
    find_package(OpenMP REQUIRED)
    add_compile_options(-qopenmp)
    link_libraries(-qopenmp)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DGNU -g -fbacktrace -cpp -fdefault-integer-8 -m64 -I$ENV{MKLROOT}/include -Wall -Wno-unused-variable")

  # if gfortran version >= 10.0, Use -fallow-argument-mismatch option to treat the mismatch between the actual and dummy arguments as a warning instead of an error.
  # (Ref: https://gcc.gnu.org/gcc-10/changes.html#:~:text=Use%20the%20new%20option%20%2Dfallow%2Dargument%2Dmismatch%20to%20turn%20these%20errors%20into%20warnings)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()

  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -fcheck=all -DDEBUG")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  link_libraries(-fdefault-integer-8 -m64)

  if(MKL)
    link_libraries(-L$ENV{MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_core -lmkl_gnu_thread -lgomp -lpthread -lm -ldl)
  endif()

  if(OPENMP)
    find_package(OpenMP REQUIRED)
    add_compile_options(-fopenmp)
    link_libraries(-fopenmp)
  endif()
else()
  message(FATAL_ERROR "Unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

add_subdirectory(src)
add_subdirectory(test)
configure_file(${PROJECT_SOURCE_DIR}/tools/dcaspt2_input ${PROJECT_BINARY_DIR}/dcaspt2 COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/tools/gen_dcaspt2_restart ${PROJECT_BINARY_DIR}/gen_dcaspt2_restart COPYONLY)
install(PROGRAMS ${PROJECT_BINARY_DIR}/dcaspt2 DESTINATION ${CMAKE_INSTALL_PREFIX})
install(PROGRAMS ${PROJECT_BINARY_DIR}/gen_dcaspt2_restart DESTINATION ${CMAKE_INSTALL_PREFIX})

# Create a custom target named commit_hash
add_custom_target(commit_hash ALL)

# Add a custom command to the target (commit_hash)
add_custom_command(
  TARGET commit_hash
  PRE_BUILD
  COMMAND ${CMAKE_COMMAND}
  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
  -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
  -P ${PROJECT_SOURCE_DIR}/cmake/commit_hash.cmake
  ALWAYS
)
install(FILES ${PROJECT_BINARY_DIR}/.commit_hash DESTINATION ${CMAKE_INSTALL_PREFIX})
