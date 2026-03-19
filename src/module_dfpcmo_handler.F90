module module_cmo_handler
    use module_cmo_variables
    use module_file_manager
    use module_global_variables
    implicit none

    private

    integer :: i, idx_irrep, unit_dfpcmo, iostat
    integer :: total_ao, total_mo
    character*150 :: line0, line1, line2, line3, line4, line5, line6, format_str

    public :: ivo_cmo_read, ivo_cmo_write
contains
    subroutine ivo_cmo_read
        implicit none
        ! read CMO data from DFPCMO file
        call open_formatted_file(unit=unit_dfpcmo, file='DFPCMO', status='old', optional_action='read')
        rewind (unit_dfpcmo)

        ! From DIRAC dirgp.F WRIPCMO (Write DHF-coefficients and eigenvalues )

        if (dirac_version >= 21 .or. integrated_caspt2) then
            read (unit_dfpcmo, '(A150)') line0
        end if
        read (unit_dfpcmo, '(A150)') line1
        if (dirac_version >= 21 .or. integrated_caspt2) then
            ! cmo_nfsym is nfsym2 in DIRAC (https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/dirac/dirgp.F#L77-78)
            ! cmo_nfsym is 1 (no inversion symmetry) or 2 (inversion symmetry)
            ! cmo_nz is 1 (real), 2 (complex) or 4 (quaterion)
            read (unit_dfpcmo, *) cmo_nfsym, cmo_nz, &
                (positronic_mo(idx_irrep), electronic_mo(idx_irrep), basis_ao(idx_irrep), idx_irrep=1, cmo_nfsym)
        else
            read (unit_dfpcmo, *) cmo_nfsym, &
                (positronic_mo(idx_irrep), electronic_mo(idx_irrep), basis_ao(idx_irrep), idx_irrep=1, cmo_nfsym)
            cmo_nz = mrconee_nz ! cannot read nz value from DFPCMO for old DIRAC (version < 21), assume that the MRCONEE nz is same as DFPCMO's
        end if
        read (unit_dfpcmo, '(A150)') line2

        if (debug .and. rank == 0) print *, 'end reading information, symmetry information and energy'

        basis_all = (positronic_mo + electronic_mo)*basis_ao
        total_ao = sum(basis_all)
        mo = positronic_mo + electronic_mo
        total_mo = sum(mo)

        allocate (eval(total_mo))
        allocate (syminfo(total_mo))
        Allocate (BUF(cmo_nz*total_ao))

        BUF = 0.0d+00
        if (dirac_version >= 21 .or. integrated_caspt2) then
            read (unit_dfpcmo, '(A150)') line3
        end if
        ! Read MO coefficient of DFPCMO
        read (unit_dfpcmo, *, iostat=iostat) BUF

        if (debug .and. rank == 0) print *, 'end reading MO coefficient'

        if (dirac_version >= 21 .or. integrated_caspt2) then
            read (unit_dfpcmo, '(A150)') line4
        end if

        Read (unit_dfpcmo, *) eval
        if (rank == 0) then
            if (debug) then
                Do i = 1, total_mo
                    print *, "eval(i)", eval(i)
                End do
            end if
            print *, 'end reading eigenvalue'
        end if
        if (dirac_version >= 21 .or. integrated_caspt2) then
            read (unit_dfpcmo, '(A150)') line5
        end if

        Read (unit_dfpcmo, *) syminfo
        if (rank == 0) then
            if (debug) then
                Do i = 1, total_mo
                    print *, "syminfo(i)", syminfo(i)
                End do
            end if
            print *, 'end reading symmetry information'
        end if

        if (.not. is_eof(unit=unit_dfpcmo, file='DFPCMO', is_formatted=.true.)) then
            ! KAPPA infomation
            allocate (kappa(total_mo))
            if (dirac_version >= 21 .or. integrated_caspt2) then
                read (unit_dfpcmo, '(A150)') line6
            end if
            read (unit_dfpcmo, *) kappa
        end if
        close (unit_dfpcmo)
    end subroutine ivo_cmo_read

    subroutine ivo_cmo_write
        ! Create new DFPCMO : DFPCMONEW
        if (rank == 0) then
            call open_formatted_file(unit=unit_dfpcmo, file='DFPCMONEW', status='replace', optional_action="write")
            if (dirac_version >= 21 .or. integrated_caspt2) then
                write (unit_dfpcmo, '(A150)') line0
            end if
            write (unit_dfpcmo, '(A150)') line1
            if (dirac_version >= 21 .or. integrated_caspt2) then
                if (cmo_nfsym == 1) then
                    format_str = '(5(X,I0))'
                else
                    format_str = '(8(X,I0))'
                end if
                write (unit_dfpcmo, format_str) cmo_nfsym, cmo_nz, (positronic_mo(i), electronic_mo(i), basis_ao(i), i=1, cmo_nfsym)
            else
                if (cmo_nfsym == 1) then
                    format_str = '(4(X,I0))'
                else
                    format_str = '(7(X,I0))'
                end if
                write (unit_dfpcmo, format_str) cmo_nfsym, (positronic_mo(i), electronic_mo(i), basis_ao(i), i=1, cmo_nfsym)
            end if
            write (unit_dfpcmo, '(A150)') line2
            if (dirac_version >= 21 .or. integrated_caspt2) then
                write (unit_dfpcmo, '(A150)') line3
            end if

            if (dirac_version >= 26 .or. integrated_caspt2) then
                write (unit_dfpcmo, '(6G25.17)') BUF
            else
                write (unit_dfpcmo, '(6F22.16)') BUF
            end if
            if (dirac_version >= 21 .or. integrated_caspt2) then
                write (unit_dfpcmo, '(A150)') line4
            end if

            if (dirac_version >= 26 .or. integrated_caspt2) then
                write (unit_dfpcmo, '(6G25.17)') eval
            else
                write (unit_dfpcmo, '(6E22.12)') eval
            end if
            if (dirac_version >= 21 .or. integrated_caspt2) then
                write (unit_dfpcmo, '(A150)') line5
            end if
            write (unit_dfpcmo, '(66(X,I0))') (syminfo(i), i=1, total_mo)

            if (allocated(kappa)) then
                if (dirac_version >= 21 .or. integrated_caspt2) then
                    write (unit_dfpcmo, '(A150)') line6
                end if
                write (unit_dfpcmo, '(66(X,I0))') kappa
            end if

            close (unit_dfpcmo)
        end if
    end subroutine ivo_cmo_write
end module module_cmo_handler
