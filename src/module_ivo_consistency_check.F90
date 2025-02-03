module module_ivo_consistency_check
    implicit none
    ! This module contains subroutines to check the consistency of the input and DFPCMO data
    ! It is called from the main program
    ! Author : Kohei Noda

contains
    subroutine ivo_consistency_check
        use module_global_variables, only: irpamo, dirac_version, ninact, nact, nsec, nsymrpa, occ_mo_num, vcut_mo_num, rank
        use module_file_manager
        use module_error
        implicit none
        integer :: unit_dfpcmo, iostat
        logical :: end_of_file
        character(150) :: line
        character(:), allocatable :: filename
        integer :: A, B
        real(8), allocatable :: evals(:), BUF(:)
        integer, allocatable :: supersym(:), kappa(:)
        integer :: isym, nv_input, nv_dfpcmo, start_isym, end_isym, isym_for_supersym
        integer :: start_idx_input, end_idx_input, start_idx_dfpcmo, end_idx_dfpcmo
        integer :: i, total_mo, total_ao
        integer :: positronic_mo(2), electronic_mo(2), basis_ao(2), basis_all(2), mo(2)

        if (rank == 0) print *, "Start checking the consistency of your input and DFPCMO data"
        positronic_mo(:) = 0; electronic_mo(:) = 0; basis_ao(:) = 0; basis_all(:) = 0; mo(:) = 0
        filename = "DFPCMO"
        call open_formatted_file(unit=unit_dfpcmo, file=filename, status="old")
        ! If DIRAC version is >= 21, additional information is printed in the DFPCMO file
        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A)', iostat=iostat) line ! INFO
            call check_end_of_file
        end if

        read (unit_dfpcmo, '(A)', iostat=iostat) line ! (e.g.) UO2                                                Sun Apr 16 12:55:22 2023
        call check_end_of_file

        if (dirac_version >= 21) then
            read (unit_dfpcmo, *, iostat=iostat) A, B, (positronic_mo(i), electronic_mo(i), basis_ao(i), i=1, A)
        else ! DIRAC version < 21
            read (unit_dfpcmo, *, iostat=iostat) A, (positronic_mo(i), electronic_mo(i), basis_ao(i), i=1, A)
        end if
        basis_all = (positronic_mo + electronic_mo)*basis_ao
        mo = positronic_mo + electronic_mo
        total_mo = sum(positronic_mo) + sum(electronic_mo)
        total_ao = sum(basis_all)
        call check_end_of_file

        read (unit_dfpcmo, '(A)', iostat=iostat) line ! (e.g.) -0.2818820677301684E+05
        call check_end_of_file

        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A)', iostat=iostat) line ! COEFF
            call check_end_of_file
        end if

        ! Read MO coefficients
        allocate (BUF(total_ao))
        read (unit_dfpcmo, *, iostat=iostat) BUF

        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A)', iostat=iostat) line ! EVALS
            call check_end_of_file
        end if

        ! Read MO energies
        allocate (evals(total_mo))
        read (unit_dfpcmo, *, iostat=iostat) evals
        call check_end_of_file

        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A)', iostat=iostat) line ! SUPERSYM
            call check_end_of_file
        end if

        ! Read MO supersymmetry
        allocate (supersym(total_mo))
        read (unit_dfpcmo, *, iostat=iostat) supersym
        if (rank == 0) print *, "supersym", supersym

        if (.not. is_eof(unit=unit_dfpcmo, file=filename, is_formatted=.true.)) then
            allocate (kappa(total_mo))
            ! Read KAPPA info
            if (dirac_version >= 21) then
                read (unit_dfpcmo, '(A)', iostat=iostat) line ! KAPPA
                print *, line
                call check_end_of_file
            end if
            read (unit_dfpcmo, *, iostat=iostat) kappa
            call check_iostat(iostat, filename, end_of_file)
        end if
        close (unit_dfpcmo) ! Close the DFPCMO file

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
                ! Print the virtual MOs supersymmetry
                print *, "Virtual supersym", supersym(start_idx_dfpcmo:end_idx_dfpcmo)
                ! Check the consistency of the input and DFPCMO data
                ! At any isym(irreducible representation)
                ! the number of virtual MOs in the DFPCMO file must be equal to the number of virtual MOs in the input file
                print *, "irpamo", irpamo(start_idx_input:end_idx_input)
            end if
            do isym = start_isym, end_isym, 2
                if (i == 1) then
                    isym_for_supersym = isym
                else
                    isym_for_supersym = isym - nsymrpa/2
                end if
                if (all(supersym(:) == 0)) then
                    nv_dfpcmo = end_idx_dfpcmo - start_idx_dfpcmo + 1  ! All virtual MOs are in the same irreducible representation
                else
                    nv_dfpcmo = count(abs(supersym(start_idx_dfpcmo:end_idx_dfpcmo)) == isym_for_supersym)  ! Number of virtual MOs corresponding to isym in the DFPCMO file
                end if
                nv_input = count(irpamo(start_idx_input:end_idx_input) == isym)  ! Number of virtual MOs corresponding to isym in the input file
                if (rank == 0) print *, "isym", isym, "isym_f_s", isym_for_supersym, "nv_dfpcmo", nv_dfpcmo, "nv_input", nv_input
                if (nv_input /= nv_dfpcmo) then
                    if (rank == 0) then
                        print *, "isym =", isym, "supersym =", supersym(start_idx_dfpcmo:end_idx_dfpcmo)
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
