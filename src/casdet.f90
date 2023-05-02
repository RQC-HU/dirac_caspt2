! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE casdet

! Find the CASCI determinant and store the index of the CASCI determinant in cas_idx
! Also, store the reverse index of the CASCI determinant in cas_idx_reverse.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module
    use module_error, only: stop_with_errorcode
    use ras_det_check
    Implicit NONE

    integer :: t, isym, allow_det_num, current_det
    integer, allocatable  :: cas_idx0(:)
    logical :: is_det_allow

    allow_det_num = 0
    if (rank == 0) print *, 'Enter casdet'
    Allocate (cas_idx0(ndet))
    Allocate (cas_idx_reverse(2**nact - 1)); call memplus(kind(cas_idx_reverse), size(cas_idx_reverse), 1)
    cas_idx0 = 0
    cas_idx_reverse = 0
    ndet = 0

    ! ===============================================
    ! Find the CASCI determinant and store the index of the CASCI determinant in cas_idx0.
    ! Also, store the reverse index of the CASCI determinant in cas_idx_reverse.
    ! Loop over only the nelec electron determinants
    ! ===============================================
    current_det = 2**nelec - 1 ! First determinant that is the number of electrons are nact. (e.g. 1111 for 4 electrons)
    do while (current_det < 2**nact)
        if (ras1_size /= 0) then
            is_det_allow = ras1_det_check(current_det, ras1_max_hole)
            if (.not. is_det_allow) cycle
        end if

        if (ras3_size /= 0) then
            is_det_allow = ras3_det_check(current_det, ras3_max_elec)
            if (.not. is_det_allow) cycle
        end if

        allow_det_num = allow_det_num + 1
        if (nsymrpa == 1) then
            ndet = ndet + 1
            cas_idx0(ndet) = current_det
            cas_idx_reverse(current_det) = ndet
        else
            Call detsym(current_det, isym)
            if (rank == 0) print '(a,i20,a,b50,a,i5)', "i:", current_det, "bit(i)", current_det, "isym:", isym
            if (isym == totsym) then
                ndet = ndet + 1
                cas_idx0(ndet) = current_det
                cas_idx_reverse(current_det) = ndet
            end if
        End if

        ! http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        ! is the reference implementation of the following code (public domain)
        ! Thanks to Dario Sneidermanis of Argentina, who provided this on November 28, 2009.
        t = or(current_det, current_det - 1) + 1
        current_det = or(t, ishft(and(t,-t)/and(current_det, -current_det), -1) - 1)  ! current_det changed to the next permutation
    End do

    ! Stop the program if ndet == 0 because ndet == 0 means the number of CASCI determinant.
    if (ndet == 0) then
        if (rank == 0) then
            print *, "[ERROR]: The number of CASCI determinant is 0. Therefore, subsequent calculations", &
                " cannot be performed successfully and the program is terminated."
        end if
        call stop_with_errorcode(1)
    end if

    Allocate (cas_idx(ndet))
    cas_idx(1:ndet) = cas_idx0(1:ndet)
    if (rank == 0) then
        print *, 'allow  = ', allow_det_num
        print *, 'totsym = ', totsym
        print *, 'ndet   = ', ndet
    end if
    Deallocate (cas_idx0)

end subroutine casdet

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE detsym(ii, isym)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE

    integer, intent(in)  :: ii
    integer, intent(out) :: isym

    integer :: i, j, jsym, ielec, isym1

    isym = 1
    ielec = 0
    Do i = 1, nact
        if (btest(ii, i - 1) .eqv. .true.) then
            ielec = ielec + 1
            j = i + ninact
            jsym = irpamo(j)
            if (mod(ielec, 2) == 1) then
                isym1 = MULTB_DS(jsym, isym) ! isym will be double irrep: odd number of electron
                isym = isym1
                if (isym1 > nsymrpa .and. rank == 0) print *, 'ielec, ii, isym, jsym, isym1', ielec, ii, isym, jsym + 1, isym1
            else
                if (mod(jsym, 2) == 1) then
                    isym1 = MULTB_D(jsym + 1, isym) ! isym will be single irrep: even number of electron !MULTB_D is (fai*|fai)
                    isym = isym1
                    if (isym1 > nsymrpa .and. rank == 0) print *, 'ielec, ii, isym, jsym+1, isym1', ielec, ii, isym, jsym + 1, isym1
                else
                    isym1 = MULTB_D(jsym - 1, isym) ! isym will be single irrep: even number of electron
                    isym = isym1
                    if (isym1 > nsymrpa .and. rank == 0) print *, 'ielec, ii, isym, jsym-1, isym1', ielec, ii, isym, jsym - 1, isym1
                end if
            end if

        End if
    End do
    If (mod(ielec, 2) == 0) isym = isym + nsymrpa ! even number electronic system

end subroutine detsym
