add_library(dcaspt2_module_global STATIC
  module_dict.f90
  module_error.f90
  module_global_variables.f90
  module_mpi.f90
  module_time.f90
)

add_library(dcaspt2_module STATIC
  module_2integrals.f90
  module_blas.f90
  module_essential_input.f90
  module_file_manager.f90
  module_index_utils.f90
  module_intra.f90
  module_realonly.f90
  module_restart_file.f90
  module_sort_swap.f90
  ras_det_check.f90
  read_input_module.f90
  takekr.f90

  # CASPT2 modules
  ulambda_s_half.f90

  # IVO modules
  module_ivo_consistency_check.f90
)
target_link_libraries(dcaspt2_module dcaspt2_module_global)

add_executable(r4dcaspt2exe

  # General
  create_newmdcint.f90
  cutoff.f90
  dcaspt2_deallocate.f90
  dcaspt2_main.f90
  dcaspt2_run_subprograms.f90
  density_matrix.f90
  diag.f90
  e0test.f90
  get_filename.f90
  mem.f90
  one_e_exct.f90
  print_head.f90
  read_cidata.f90
  read_mrconee.f90
  tramo.f90

  # IVO
  fockivo.f90
  print_head_ivo.f90
  r4divo.f90

  # CASCI
  calce0.f90
  casci.f90
  casmat.f90
  checkdgc.f90
  print_head_casci.f90
  r4dcasci.f90
  search_cas_configuration.f90

  # CASPT2
  divide_2_elec_integral_into_subspaces.f90
  e0after_tra.f90
  fock_matrix_of_hf.f90
  fockcasci.f90
  fockdiag.f90
  print_head_caspt2.f90
  r4dcaspt2_tra.f90
  sdet.f90
  solve_A_subspace.f90
  solve_B_subspace.f90
  solve_C_subspace.f90
  solve_D_subspace.f90
  solve_E_subspace.f90
  solve_F_subspace.f90
  solve_G_subspace.f90
  solve_H_subspace.f90
  trac.f90
  prtoutfock.f90
)
target_link_libraries(r4dcaspt2exe dcaspt2_module)

install(TARGETS r4dcaspt2exe RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
