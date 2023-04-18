module module_ivo_consistency_check
    implicit none
    ! This module contains subroutines to check the consistency of the input and DFPCMO data
    ! It is called from the main program
    ! Author : Kohei Noda

contains
    subroutine ivo_consistency_check
        use module_global_variables, only: irpamo, dirac_version, ninact, nact, nsec, noccg, noccu, nvcutg, nvcutu, nsymrpa
        use module_file_manager
        use module_error
        implicit none
        integer :: unit_dfpcmo, iostat
        logical :: end_of_file
        character(150) :: line
        character(:), allocatable :: filename
        integer :: A, B, npg, neg, nbasg, npu, neu, nbasu
        integer :: ao_all, mo_all
        real(8), allocatable :: evals(:), BUF(:)
        integer, allocatable :: supersym(:)
        integer :: isym, nv_input, nv_dfpcmo
        integer :: start_idx_input, end_idx_input, start_idx_gerade, end_idx_gerade, start_idx_ungerade, end_idx_ungerade

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
            read (unit_dfpcmo, *, iostat=iostat) A, B, npg, neg, nbasg, npu, neu, nbasu
        else ! DIRAC version < 21
            read (unit_dfpcmo, *, iostat=iostat) A, npg, neg, nbasg, npu, neu, nbasu
        end if
        call check_end_of_file

        read (unit_dfpcmo, '(A)', iostat=iostat) line ! (e.g.) -0.2818820677301684E+05
        call check_end_of_file

        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A)', iostat=iostat) line ! COEFF
            call check_end_of_file
        end if

        ! Read MO coefficients
        ao_all = (npg + neg)*nbasg + (npu + neu)*nbasu
        allocate (BUF(ao_all))
        read (unit_dfpcmo, *, iostat=iostat) BUF

        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A)', iostat=iostat) line ! EVALS
            call check_end_of_file
        end if

        ! Read MO energies
        mo_all = (npg + neg) + (npu + neu)
        allocate (evals(mo_all))
        read (unit_dfpcmo, *, iostat=iostat) evals
        call check_end_of_file

        if (dirac_version >= 21) then
            read (unit_dfpcmo, '(A)', iostat=iostat) line ! SUPERSYM
            call check_end_of_file
        end if

        ! Read MO supersymmetry
        allocate (supersym(mo_all))
        read (unit_dfpcmo, *, iostat=iostat) supersym
        print *, "supersym", supersym

        ! End of file check
        read (unit_dfpcmo, '(A)', iostat=iostat) line ! iostat should be -1 (end of file)
        call check_iostat(iostat=iostat, file=filename, end_of_file_reached=end_of_file)
        if (.not. end_of_file) then
            ! Additional information is printed in the DFPCMO file
            print *, "Error: The DFPCMO file is not correctly formatted"
            call stop_with_errorcode(1)
        end if

        close (unit_dfpcmo) ! Close the DFPCMO file

        ! Define the indices of the virtual MOs in the input and DFPCMO data
        start_idx_input = ninact + nact + 1
        end_idx_input = ninact + nact + nsec
        start_idx_gerade = npg + noccg + 1
        end_idx_gerade = npg + neg - nvcutg
        start_idx_ungerade = npg + neg + npu + noccu + 1
        end_idx_ungerade = npg + neg + npu + neu - nvcutu

        ! Print the virtual MOs supersymmetry
        print *, "Gerade virtual supersym", supersym(start_idx_gerade:end_idx_gerade)
        print *, "Ungerade virtual supersym", supersym(start_idx_ungerade:end_idx_ungerade)
        ! Check the consistency of the input and DFPCMO data
        ! At any isym(irreducible representation)
        ! the number of virtual MOs in the DFPCMO file must be equal to the number of virtual MOs in the input file
        print *, "irpamo", irpamo(start_idx_input:end_idx_input)
        do isym = 1, nsymrpa, 2
            nv_input = count(irpamo(start_idx_input:end_idx_input) == isym)  ! Number of virtual MOs corresponding to isym in the input file

            if (isym <= nsymrpa/2) then  ! Gerade
                nv_dfpcmo = count(abs(supersym(start_idx_gerade:end_idx_gerade)) == isym)  ! Number of virtual MOs corresponding to isym in the DFPCMO file
                if (nv_input /= nv_dfpcmo) then
                    print *, "isym =", isym, "supersym =", supersym(start_idx_gerade:end_idx_gerade)
                    print *, "The number of virtual MOs in the DFPCMO file is not equal", &
                        "to the number of virtual MOs in the input file.", &
                        "isym = ", isym, "nv_input = ", nv_input, "nv_dfpcmo = ", nv_dfpcmo
                    call stop_with_errorcode(1)
                end if
            else ! Ungerade
                nv_dfpcmo = count(abs(supersym(start_idx_ungerade:end_idx_ungerade)) == isym - nsymrpa/2)  ! Number of virtual MOs corresponding to isym in the DFPCMO file
                if (nv_input /= nv_dfpcmo) then
                    print *, "isym =", isym, "supersym =", supersym(start_idx_ungerade:end_idx_ungerade)
                    print *, "The number of virtual MOs in the DFPCMO file is not equal", &
                        "to the number of virtual MOs in the input file.", &
                        "isym = ", isym, "nv_input = ", nv_input, "nv_dfpcmo = ", nv_dfpcmo
                    call stop_with_errorcode(1)
                end if
            end if
        end do

    contains
        subroutine check_end_of_file
            implicit none
            call check_iostat(iostat=iostat, file=filename, end_of_file_reached=end_of_file)
            if (end_of_file) then
                print *, "Error: The DFPCMO file contains less data than expected."
                call stop_with_errorcode(1)
            end if
        end subroutine check_end_of_file

    end subroutine ivo_consistency_check
end module module_ivo_consistency_check
