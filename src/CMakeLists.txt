add_library(dcaspt2_module STATIC
  four_caspt2_module.f90
  module_error.f90
  module_file_manager.f90
  module_mpi.f90
  module_sort_swap.f90
  ras_det_check.f90
  read_input_module.f90
)
add_executable(r4dcascicoexe
  calce0.f90
  casci.f90
  casdet.f90
  casmat.f90
  checkdgc.f90
  create_binmdcint.f90
  cutoff.f90
  diag.f90
  dim.f90
  e0after_tra.f90
  e0test.f90
  fockcasci.f90
  fockdiag.f90
  fockhf1.f90
  get_filename.f90
  mem.f90
  one_e_exct.f90
  prtoutfock.f90
  r4dcasci.f90
  read1mo.f90
  readorb_enesym.f90
  readint2_casci.f90
  sdet.f90
  takekr.f90
  timing.f90
  trac.f90
  tramo.f90
  uramda_s_half.f90
)
target_link_libraries(r4dcascicoexe dcaspt2_module)

add_executable(r4dcaspt2ocoexe
  calce0.f90
  checkdgc.f90
  cutoff.f90
  diag.f90
  dim.f90
  e0test.f90
  get_filename.f90
  intra.f90
  mem.f90
  one_e_exct.f90
  r4dcaspt2_tra.f90
  ras_det_check.f90
  readorb_enesym.f90
  read1mo.f90
  readint2_ord.f90
  solve_A_subspace.f90
  solve_B_subspace.f90
  solve_C_subspace.f90
  solve_D_subspace.f90
  solve_E_subspace.f90
  solve_F_subspace.f90
  solve_G_subspace.f90
  solve_H_subspace.f90
  takekr.f90
  timing.f90
  tramo.f90
  uramda_s_half.f90
)
target_link_libraries(r4dcaspt2ocoexe dcaspt2_module)

add_executable(r4divocoexe
  create_binmdcint.f90
	diag.f90
	dim.f90
	e0test.f90
	fockivo.f90
	get_filename.f90
	mem.f90
	one_e_exct.f90
	r4divo.f90
	read1mo.f90
	readint2_casci.f90
	readorb_enesym.f90
	takekr.f90
	timing.f90
	uramda_s_half.f90
)
target_link_libraries(r4divocoexe dcaspt2_module)

install(TARGETS r4divocoexe RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS r4dcascicoexe RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS r4dcaspt2ocoexe RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
