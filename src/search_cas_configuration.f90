! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE search_cas_configuration

! Find the CAS configurations and store the index of the CAS configurations in cas_idx
! Also, store the reverse index of the CAS configurations in dict_cas_idx_reverse.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use module_global_variables
    use module_error, only: stop_with_errorcode
    use module_dict, only: add
    Implicit NONE

    integer :: t, allow_det_num, current_det

    if (rank == 0) print *, 'Enter search_cas_configuration'
    allow_det_num = 0
    ndet = 0

    ! ========================================================================================
    ! Find the CASCI determinant and store the index of the CASCI determinant in dict_cas_idx.
    ! Also, store the reverse index of the CASCI determinant in dict_cas_idx_reverse.
    ! Loop over only the nelec electron determinants
    ! ========================================================================================
    current_det = 2**nelec - 1 ! First determinant that is the number of electrons are nact. (e.g. 1111 for 4 electrons)
    do while (current_det < 2**nact)

        if (.not. satisfy_ras_conditions()) then
            ! Do not satisfy the RAS condition
            current_det = find_next_determinant()
            cycle
        end if

        ! RAS conditions are satisfied
        allow_det_num = allow_det_num + 1

        ! Check if the determinant is CASCI determinant
        if (is_cas_determinant()) then
            ndet = ndet + 1
            call add(dict_cas_idx, ndet, current_det)
            call add(dict_cas_idx_reverse, current_det, ndet)
        end if

        current_det = find_next_determinant()
    End do

    ! Stop the program if ndet == 0 because ndet == 0 means the number of CASCI determinant.
    if (ndet == 0) then
        if (rank == 0) then
            print *, "[ERROR]: The number of CASCI determinant is 0. Therefore, subsequent calculations", &
                " cannot be performed successfully and the program is terminated."
        end if
        call stop_with_errorcode(1)
    end if

    if (rank == 0) then
        print *, 'allow  = ', allow_det_num
        print *, 'totsym = ', totsym
        print *, 'ndet   = ', ndet
    end if

    call validate_ndet_and_input_parameters()
contains

    function find_next_determinant() result(next_determinant)
        implicit none
        integer :: next_determinant
        ! Find the next determinant
        ! http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        ! is the reference implementation of the following code (public domain)
        ! Thanks to Dario Sneidermanis of Argentina, who provided this on November 28, 2009.
        t = or(current_det, current_det - 1) + 1
        next_determinant = or(t, ishft(and(t, -t)/and(current_det, -current_det), -1) - 1) ! Find the next determinant
    end function find_next_determinant

    logical function satisfy_ras_conditions()
        use ras_det_check
        implicit none
        if (ras1_size /= 0) then
            if (.not. satisfy_ras1_condition(current_det, ras1_max_hole)) then
                ! Do not satisfy the RAS1 condition
                satisfy_ras_conditions = .false.
                return
            end if
        end if
        if (ras3_size /= 0) then
            if (.not. satisfy_ras3_condition(current_det, ras3_max_elec)) then
                ! Do not satisfy the RAS3 condition
                satisfy_ras_conditions = .false.
                return
            end if
        end if
        ! Satisfy the RAS condition
        satisfy_ras_conditions = .true.
    end function satisfy_ras_conditions

    logical function is_cas_determinant()
        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
        implicit none
        integer :: i, j, jsym, ielec, isym, isym1

        is_cas_determinant = .false. ! Default value
        if (nsymrpa == 1) then
            is_cas_determinant = .true. ! If nsymrpa == 1, all the determinants are in the same irrep.
        else
            isym = 1
            ielec = 0
            Do i = 1, nact
                if (btest(current_det, i - 1) .eqv. .true.) then
                    ielec = ielec + 1
                    j = convert_active_to_global_idx(i)
                    jsym = irpamo(j)
                    if (mod(ielec, 2) == 1) then
                        isym1 = MULTB_DS(jsym, isym) ! isym will be double irrep: odd number of electron
                        isym = isym1
                    else
                        if (mod(jsym, 2) == 1) then
                            isym1 = MULTB_D(jsym + 1, isym) ! isym will be single irrep: even number of electron !MULTB_D is (fai*|fai)
                            isym = isym1
                        else
                            isym1 = MULTB_D(jsym - 1, isym) ! isym will be single irrep: even number of electron
                            isym = isym1
                        end if
                    end if
                    if (isym1 > nsymrpa .and. rank == 0) print *, 'ielec, current_det, isym, jsym-1, isym1', &
                        ielec, current_det, isym, jsym - 1, isym1
                End if
            End do
            If (mod(ielec, 2) == 0) isym = isym + nsymrpa ! even number electronic system
#ifdef DEBUG
            if (rank == 0) print '(a,i20,1x,a,b50,1x,a,i5)', &
                "current_det:", current_det, "bit(current_det):", current_det, "isym:", isym
#endif
            ! Check if the determinant is allowed
            if (isym == totsym) then
                is_cas_determinant = .true.
            else
                is_cas_determinant = .false.
            end if
        end if
    end function is_cas_determinant

    subroutine validate_ndet_and_input_parameters()
        implicit none
        if (ndet < selectroot) then
            if (rank == 0) then
                print *, "ERROR: ndet < selectroot"
                print '(A,I0,A)', "Cannot calculate ", selectroot, "th RASCI/CASCI energy"
                print '(A,I0,A,I0)', "because the number of CASCI determinant is ", ndet, " and it is less than ", selectroot
                print *, "Please increase the number of active orbitals or the number of electrons"
                print *, "or decrease the number of selected root."
                print *, "Exit the program."
            end if
            call stop_with_errorcode(1)
        end if

        if (ndet < nroot) then
            if (rank == 0) then
                print *, "WARNING: ndet < nroot"
                print '(A,I0,A)', "Cannot print ", nroot, "th RASCI/CASCI energy"
                print '(A,I0,A,I0)', "because the number of CASCI determinant is ", ndet, " and it is less than ", nroot
                print *, "Therefore, explicitly replace nroot with the number of CASCI determinant."
                print *, "new nroot = ", ndet
            end if
            nroot = ndet
        end if
    end subroutine validate_ndet_and_input_parameters

end subroutine search_cas_configuration
