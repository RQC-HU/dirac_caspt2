! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

subroutine r4divo_co   ! DO IVO CALC ONLY FOR SMALL BASIS SETS

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_file_manager
    use module_2integrals
    use module_realonly, only: realonly
    use module_time
    use read_input_module
    use module_ivo_consistency_check

    Implicit NONE
    integer                     :: input_unit, nuniq
    character(:), allocatable   :: filename
    character(*), parameter     :: int_input_form = '(1x,a,1x,i0)'
    character(len=30)           :: real_str

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    if (rank == 0) then
        call print_head_ivo
        print *, ''
        print *, 'START RELATIVISTIC IVO PROGRAM'
    end if

    tmem = 0.0d+00
    call write_allocated_memory_size
    call get_current_time(init_time); call print_time(init_time); start_time = init_time

    if (rank == 0) then
        print int_input_form, 'ninact     =', ninact
        print int_input_form, 'nact       =', nact
        print int_input_form, 'nsec       =', nsec
        print int_input_form, 'nelec      =', nelec
        print int_input_form, 'nroot      =', nroot
        print int_input_form, 'selectroot =', selectroot
        print int_input_form, 'totsym     =', totsym
        write (real_str, '(E20.10)') eshift
        print '(1x,a,1x,a)', 'eshift     =', trim(adjustl(real_str))          ! NO USE IN IVO BUT FOR CASCI AND CASPT2 IT IS USED
        print int_input_form, 'nhomo      =', nhomo
        if (inversion) then
            print int_input_form, "noccg      =", occ_mo_num(1)
            print int_input_form, "noccu      =", occ_mo_num(2)
            print int_input_form, "nvcutg     =", vcut_mo_num(1)
            print int_input_form, "nvcutu     =", vcut_mo_num(2)
        else
            print int_input_form, "nocc       =", occ_mo_num(1)
            print int_input_form, "nvcut      =", vcut_mo_num(1)
        end if
        print int_input_form, 'diracver   =', dirac_version
        print int_input_form, 'scheme     =', mdcint_scheme
        print *, 'debugprint =', debug
        if (enable_restart) print *, "restart    =", enable_restart
        print *, ''
    end if

    ! Check consistency of IVO input and DFPCMO file.
    call ivo_consistency_check

    ! Create UTChem type MDCINT file from Dirac MDCINT file
    if (rank == 0) print *, "Create_newmdcint"
    call create_newmdcint

    call get_mdcint_filename(0)
    ! Read UTChem type MDCINT files and expands the 2-electron integral in memory
    call readint2_casci(mdcintnew, nuniq)

    if (rank == 0) then
        call write_allocated_memory_size
        print *, ' '
        print *, '*******************************'
        print *, ' '
        ! print *, 'IREP IS ', repna(totsym)
        print *, ' '
        print *, '*******************************'
    end if
    iroot = selectroot

    if (rank == 0) then
        Allocate (fock_cmplx(nsec, nsec)); Call memplus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)

        fock_cmplx(:, :) = 0.0d+00

!! NOW MAKE FOCK MATRIX FOR IVO (only virtual spinors)
!! fij = hij + SIGUMA_a(ij|aa)-(ia|aj)}

        Call fockivo

        Call memminus(KIND(fock_cmplx), SIZE(fock_cmplx), 2); deallocate (fock_cmplx)
    end if

    ! Deallocate memory
    if (allocated(inttwi)) then
        Call memminus(KIND(inttwi), SIZE(inttwi), 1); deallocate (inttwi)
    end if
    if (allocated(inttwr)) then
        Call memminus(KIND(inttwr), SIZE(inttwr), 1); deallocate (inttwr)
    end if
    if (allocated(int2r_f1)) then
        Call memminus(KIND(int2r_f1), SIZE(int2r_f1), 1); deallocate (int2r_f1)
    end if
    if (allocated(int2r_f2)) then
        Call memminus(KIND(int2r_f2), SIZE(int2r_f2), 1); deallocate (int2r_f2)
    end if
    if (allocated(int2i_f1)) then
        Call memminus(KIND(int2i_f1), SIZE(int2i_f1), 1); deallocate (int2i_f1)
    end if
    if (allocated(int2i_f2)) then
        Call memminus(KIND(int2i_f2), SIZE(int2i_f2), 1); deallocate (int2i_f2)
    end if

    call write_allocated_memory_size
    call get_current_time_and_print_diff(init_time, end_time) ! print the total time
    if (rank == 0) print *, 'END OF RELATIVISTIC IVO PROGRAM'

END subroutine r4divo_co
