add_library(dcaspt2_module_global STATIC
  module_dict.F90
  module_error.F90
  module_global_variables.F90
  module_mpi.F90
  module_time.F90
)

add_library(dcaspt2_module STATIC
  module_2integrals.F90
  module_blas.F90
  module_essential_input.F90
  module_file_manager.F90
  module_index_utils.F90
  module_intra.F90
  module_realonly.F90
  module_restart_file.F90
  module_sort_swap.F90
  ras_det_check.F90
  read_input_module.F90
  takekr.F90

  # CASPT2 modules
  ulambda_s_half.F90

  # IVO modules
  module_ivo_consistency_check.F90
)
target_link_libraries(dcaspt2_module dcaspt2_module_global)

add_executable(r4dcaspt2exe

  # General
  create_newmdcint.F90
  cutoff.F90
  dcaspt2_deallocate.F90
  dcaspt2_main.F90
  dcaspt2_run_subprograms.F90
  density_matrix.F90
  diag.F90
  e0test.F90
  get_filename.F90
  mem.F90
  one_e_exct.F90
  print_head.F90
  read_cidata.F90
  read_mrconee.F90
  tramo.F90

  # IVO
  fockivo.F90
  print_head_ivo.F90
  r4divo.F90

  # CASCI
  calce0.F90
  casci.F90
  casmat.F90
  checkdgc.F90
  print_head_casci.F90
  r4dcasci.F90
  search_cas_configuration.F90

  # CASPT2
  divide_2_elec_integral_into_subspaces.F90
  e0after_tra.F90
  fock_matrix_of_hf.F90
  fockcasci.F90
  fockdiag.F90
  print_head_caspt2.F90
  r4dcaspt2_tra.F90
  sdet.F90
  solve_A_subspace.F90
  solve_B_subspace.F90
  solve_C_subspace.F90
  solve_D_subspace.F90
  solve_E_subspace.F90
  solve_F_subspace.F90
  solve_G_subspace.F90
  solve_H_subspace.F90
  trac.F90
  prtoutfock.F90
)
target_link_libraries(r4dcaspt2exe dcaspt2_module)

# r4dcaspt2exe is placed in the PROJECT_BINARY_DIR on build time
set_target_properties(r4dcaspt2exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

install(TARGETS r4dcaspt2exe RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
