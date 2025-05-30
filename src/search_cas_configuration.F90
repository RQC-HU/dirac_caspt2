! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE search_cas_configuration

! Find the CAS configurations and store the index of the CAS configurations in cas_idx
! Also, store the reverse index of the CAS configurations in dict_cas_idx_reverse.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use, intrinsic :: iso_fortran_env, only: int64
    use module_global_variables
    use module_validation, only: validate_ndet
    use module_error, only: stop_with_errorcode
    use module_dict, only: add, destruct_dict
    Implicit NONE

    integer(kind=int64) :: idx, t, allow_det_num, current_det
    integer             :: cur_sym
    integer(kind=int64), allocatable :: det_cnt(:)

    if (debug .and. rank == 0) print *, 'Enter search_cas_configuration'
    allocate (det_cnt(2*nsymrpa))
    det_cnt(:) = 0
    allow_det_num = 0
    ndet = 0
    call destruct_dict(dict_cas_idx)
    call destruct_dict(dict_cas_idx_reverse)

    ! ========================================================================================
    ! Find the CASCI configuration and store the index of the CASCI configuration in dict_cas_idx.
    ! Also, store the reverse index of the CASCI configuration in dict_cas_idx_reverse.
    ! Loop over only the nelec electron configurations
    ! ========================================================================================
    current_det = 2**nelec - 1 ! First configuration that is the number of electrons are nact. (e.g. 1111 for 4 electrons)
    do while (current_det < 2**nact)

        if (.not. satisfy_ras_conditions()) then
            ! Do not satisfy the RAS condition
            current_det = find_next_configuration()
            cycle
        end if

        ! RAS conditions are satisfied
        allow_det_num = allow_det_num + 1

        ! Calculate the total symmetry of the configuration and count-up the number of configurations for calculated symmetry.
        cur_sym = cas_configuration_totsym()
        det_cnt(cur_sym) = det_cnt(cur_sym) + 1

        ! Check if the configuration is CASCI configuration
        if (is_cas_configuration(cur_sym)) then
            ndet = ndet + 1
            call add(dict_cas_idx, int(ndet, kind=int64), current_det)
            call add(dict_cas_idx_reverse, current_det, int(ndet, kind=int64))
        end if

        current_det = find_next_configuration()
    End do

    if (docountndet) then
        if (rank == 0) then
            print *, 'Numbers of CASCI configurations for each total symmetry'
            do idx = 1, 2*nsymrpa
                print '(2(a,i0))', "ndet( ", idx, " ) = ", det_cnt(idx)
            end do
        end if
        return ! skip the rest of validation
    end if

    if (rank == 0) then
        print *, 'Number of candidates for configuration of CASCI = ', allow_det_num
        print *, 'total symmetry of CASCI configuration = ', totsym
        print *, 'Number of CASCI configuration = ', ndet
    end if

    call check_det_cnt
    call validate_ndet
contains

    function find_next_configuration() result(next_configuration)
        implicit none
        integer(kind=int64) :: next_configuration
        ! Find the next configuration
        ! http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        ! is the reference implementation of the following code (public domain)
        ! Thanks to Dario Sneidermanis of Argentina, who provided this on November 28, 2009.
        t = or(current_det, current_det - 1) + 1
        next_configuration = or(t, ishft(and(t, -t)/and(current_det, -current_det), -1) - 1) ! Find the next configuration
    end function find_next_configuration

    logical function satisfy_ras_conditions()
        use ras_det_check
        implicit none
        if (ras1_size /= 0) then
            if (.not. satisfy_ras1_condition(current_det, int(ras1_max_hole, kind=int64))) then
                ! Do not satisfy the RAS1 condition
                satisfy_ras_conditions = .false.
                return
            end if
        end if
        if (ras3_size /= 0) then
            if (.not. satisfy_ras3_condition(current_det, int(ras3_max_elec, kind=int64))) then
                ! Do not satisfy the RAS3 condition
                satisfy_ras_conditions = .false.
                return
            end if
        end if
        ! Satisfy the RAS condition
        satisfy_ras_conditions = .true.
    end function satisfy_ras_conditions

    function cas_configuration_totsym() result(ret_isym)
        use module_index_utils, only: convert_active_to_global_idx
        implicit none
        integer :: i, j, jsym, ielec
        integer :: ret_isym  ! Return value

        ret_isym = 1 ! Default value
        if (nsymrpa == 1) then
            return  ! If nsymrpa == 1, all the configurations are irrep 1.
        else
            ielec = 0
            Do i = 1, nact
                if (btest(current_det, i - 1)) then
                    ielec = ielec + 1
                    j = convert_active_to_global_idx(i)
                    jsym = irpamo(j)
                    if (mod(ielec, 2) == 1) then
                        ret_isym = MULTB_DS(jsym, ret_isym) ! ret_isym will be double irrep: odd number of electron
                    else
                        if (mod(jsym, 2) == 1) then
                            ret_isym = MULTB_D(jsym + 1, ret_isym) ! ret_isym will be single irrep: even number of electron !MULTB_D is (fai*|fai)
                        else
                            ret_isym = MULTB_D(jsym - 1, ret_isym) ! ret_isym will be single irrep: even number of electron
                        end if
                    end if
                    if (ret_isym > nsymrpa .and. rank == 0) then
                        print *, 'ielec, current_det, ret_isym, jsym-1', &
                            ielec, current_det, ret_isym, jsym - 1
                    end if
                End if
            End do
            If (mod(ielec, 2) == 0) ret_isym = ret_isym + nsymrpa ! even number electronic system
#ifdef DEBUG
            if (rank == 0) print '(a,i20,1x,a,b50,1x,a,i5)', &
                "current_det:", current_det, " bit(current_det):", current_det, " isym:", ret_isym
#endif
        end if
    end function cas_configuration_totsym

    logical function is_cas_configuration(isym)
        implicit none
        integer, intent(in) :: isym
        ! Check if the configuration is allowed
        if (nsymrpa == 1) then
            is_cas_configuration = .true.  ! If nsymrpa == 1, all configurations are CASCI configuration
        else
            is_cas_configuration = isym == totsym  ! isym == totsym means that the configuration is CASCI configuration
        end if
    end function is_cas_configuration

    subroutine check_det_cnt()
        implicit none

        ! Stop the program if ndet == 0 because the number of CASCI configuration is 0.
        ! It means that subsequent calculations cannot be performed successfully.
        if (ndet == 0) then
            if (rank == 0) then
                print *, "[ERROR]: The number of CASCI configuration is 0. Therefore, subsequent calculations", &
                    " cannot be performed successfully and the program is terminated."
                print '(a,i0,a)', " Please check your totsym parameter in your input file. your totsym parameter: ", totsym, "."
                print *, "Following are the number of configurations for each totsym."
                do idx = 1, 2*nsymrpa
                    print '(2(a,i0))', "ndet( ", idx, " ) = ", det_cnt(idx)
                end do
            end if
            call stop_with_errorcode(1)
        end if
    end subroutine check_det_cnt

end subroutine search_cas_configuration
