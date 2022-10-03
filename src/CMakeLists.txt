cmake_minimum_required(VERSION 3.7)

add_executable(r4dcascicoexe
	calce0.f90
	casci_ty.f90
	casdet_ty.f90
	casmat.f90
	checkdgc.f90
	create_binmdcint.f90
	cutoff.f90
	diag.f90
	dim.f90
	e0after_tra_ty.f90
	e0test_v2.f90
	fockcasci_ty.f90
	fockdiag_ty.f90
	fockhf1_ty.f90
	four_caspt2_module.f90
	get_filename.f90
	mem.f90
	module_file_manager.f90
	module_sort_swap.f90
	one_e_exct.f90
	prtoutfock.f90
	r4dcasci_co.f90
	ras_det_check.f90
	read_input_module.f90
	read1mo_co.f90
	readorb_enesym_co.f90
	readint2_casci_co.f90
	sdet.f90
	takekr.f90
	timing.f90
	trac.f90
	tramo_ty.f90
	uramda_s_half.f90
)

add_executable(r4dcaspt2ocoexe
	calce0.f90
	checkdgc.f90
	cutoff.f90
	diag.f90
	dim.f90
	e0test_v2.f90
	four_caspt2_module.f90
	get_filename.f90
	intra.f90
	mem.f90
	module_file_manager.f90
	module_sort_swap.f90
	one_e_exct.f90
	r4dcaspt2_tra_co.f90
	ras_det_check.f90
	read_input_module.f90
	readorb_enesym_co.f90
	read1mo_co.f90
	readint2_ord_co.f90
	solvall_A_ord_ty.f90
	solvall_B_ord_ty.f90
	solvall_C_ord_ty.f90
	solvall_D_ord_ty.f90
	solvall_E_ord_ty.f90
	solvall_F_ord_ty.f90
	solvall_G_ord_ty.f90
	solvall_H_ord_ty.f90
	takekr.f90
	timing.f90
	tramo_ty.f90
	uramda_s_half.f90
)

install(TARGETS r4dcascicoexe RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS r4dcaspt2ocoexe RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})

