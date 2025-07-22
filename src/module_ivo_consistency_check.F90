module module_ivo_consistency_check
    implicit none
    ! This module contains subroutines to check the consistency of the input and DFPCMO data
    ! It is called from the main program
    ! Author : Kohei Noda

contains
    subroutine ivo_consistency_check
        use module_cmo_variables
        use module_cmo_handler, only: ivo_cmo_read
        use module_global_variables, only: irpamo, dirac_version, integrated_caspt2, ninact, nact, nsec, nsymrpa, &
                                           occ_mo_num, vcut_mo_num, rank
        use module_file_manager
        use module_error
        implicit none
        integer :: unit_dfpcmo, iostat
        logical :: end_of_file
        character(150) :: line
        character(:), allocatable :: filename
        integer :: isym, nv_input, nv_dfpcmo, start_isym, end_isym, isym_for_syminfo
        integer :: start_idx_input, end_idx_input, start_idx_dfpcmo, end_idx_dfpcmo
        integer :: i

        if (rank == 0) print *, "Start checking the consistency of your input and DFPCMO data"

        call ivo_cmo_read
        if (rank == 0) print *, "syminfo", syminfo
        do i = 1, A
            ! Define the indices of the virtual MOs in the input and DFPCMO data
            start_idx_input = ninact + nact + 1
            end_idx_input = ninact + nact + nsec
            if (i == 1) then
                start_idx_dfpcmo = positronic_mo(i) + occ_mo_num(i) + 1
                end_idx_dfpcmo = mo(i) - vcut_mo_num(i)
                start_isym = 1
                end_isym = nsymrpa/2
            else
                start_idx_dfpcmo = mo(1) + positronic_mo(i) + occ_mo_num(i) + 1
                end_idx_dfpcmo = mo(1) + mo(i) - vcut_mo_num(i)
                start_isym = nsymrpa/2 + 1
                end_isym = nsymrpa
            end if
            if (rank == 0) then
                ! Print the virtual MOs syminfo
                print *, "Virtual syminfo", syminfo(start_idx_dfpcmo:end_idx_dfpcmo)
                ! Check the consistency of the input and DFPCMO data
                ! At any isym(irreducible representation)
                ! the number of virtual MOs in the DFPCMO file must be equal to the number of virtual MOs in the input file
                print *, "irpamo", irpamo(start_idx_input:end_idx_input)
            end if
            do isym = start_isym, end_isym, 2
                if (i == 1) then
                    isym_for_syminfo = isym
                else
                    isym_for_syminfo = isym - nsymrpa/2
                end if
                if (all(syminfo(:) == 0)) then
                    nv_dfpcmo = end_idx_dfpcmo - start_idx_dfpcmo + 1  ! All virtual MOs are in the same irreducible representation
                else
                    nv_dfpcmo = count(abs(syminfo(start_idx_dfpcmo:end_idx_dfpcmo)) == isym_for_syminfo)  ! Number of virtual MOs corresponding to isym in the DFPCMO file
                end if
                nv_input = count(irpamo(start_idx_input:end_idx_input) == isym)  ! Number of virtual MOs corresponding to isym in the input file
                if (rank == 0) print *, "isym", isym, "isym_f_s", isym_for_syminfo, "nv_dfpcmo", nv_dfpcmo, "nv_input", nv_input
                if (nv_input /= nv_dfpcmo) then
                    if (rank == 0) then
                        print *, "isym =", isym, "syminfo =", syminfo(start_idx_dfpcmo:end_idx_dfpcmo)
                        print *, "The number of virtual MOs in the DFPCMO file is not equal", &
                            "to the number of virtual MOs in the input file.", &
                            "isym = ", isym, "nv_input = ", nv_input, "nv_dfpcmo = ", nv_dfpcmo
                        print *, "Please check your input file."
                        print *, "Maybe you forgot to set the nvcut(g,u) parameter in the input file?"
                    end if
                    call stop_with_errorcode(1)
                end if
            end do
        end do
        call reset_cmo_variables

    contains
        subroutine check_end_of_file
            implicit none
            call check_iostat(iostat=iostat, file=filename, end_of_file_reached=end_of_file)
            if (end_of_file) then
                if (rank == 0) print *, "Error: The DFPCMO file contains less data than expected."
                call stop_with_errorcode(1)
            end if
        end subroutine check_end_of_file

    end subroutine ivo_consistency_check
end module module_ivo_consistency_check
