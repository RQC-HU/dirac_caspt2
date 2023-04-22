! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE casdet_ty

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module
    use module_error, only: stop_with_errorcode
    use ras_det_check
    Implicit NONE

    integer :: i, isym
    integer, allocatable  :: cas_idx0(:)
    integer :: upper_allowed_hole, allow_det_num
    logical :: is_det_allow
    upper_allowed_hole = 1 ! RAS1の許容されるホール数

    ! e.g. RAS1が4spinorの場合RAS1のビット表現であるras1_bitはどのように表せるか?
    ! -> 00001111のようなものが表せていれば良い
    ! 00001111は10進数で15(=2^0+2^1+2^+2^3)であるので
    ! ras1_bit = 15

    ! <一般のspinor数nのとき>
    ! 4spinorのときの例からras1_bitは等比数列の和となるので
    ! 等比数列の和の公式 a(1-r^n)/(1-r) (ここでa:初項,n:項数,r:公比)より
    ! 一般のspinor数n(ただしnは自然数)についてras1_bitは
    ! ras1_bit = 1*(1-2^n)/(1-2) = 2^n - 1 となる
    ! 実際に4spinorのときを確かめると
    ! ras1_bit = 2^4 - 1 = 16 - 1 = 15 となり上の4spinorの例と一致する

    ! ras1_bit = 2**2 - 1 ! RAS1のビット表現
    allow_det_num = 0
    if (rank == 0) print *, 'Enter casdet_ty'
    Allocate (cas_idx0(ndet))
    Allocate (cas_idx_reverse(2**nact - 1)); call memplus(kind(cas_idx_reverse), size(cas_idx_reverse), 1)
    cas_idx0 = 0
    cas_idx_reverse = 0
    ndet = 0
    !    67108864* 8 / (1024^2) = 500MB, 26 spinor
    Do i = 1, 2**nact - 1
        if (POPCNT(i) == nelec) then
            is_det_allow = .true.
            if (is_ras1_configured) then
                is_det_allow = ras1_det_check(i, ras1_max_hole)
                if (.not. is_det_allow) cycle
            end if

            if (is_ras3_configured) then
                is_det_allow = ras3_det_check(i, ras3_max_elec)
                if (.not. is_det_allow) cycle
            end if

            allow_det_num = allow_det_num + 1
            if (nsymrpa == 1) then
                ndet = ndet + 1
                cas_idx0(ndet) = i
                cas_idx_reverse(i) = ndet
            else
                Call detsym_ty(i, isym)
                if (rank == 0) print '(a,i20,a,b50,a,i5)', "i:", i, "bit(i)", i, "isym:", isym
                if (isym == totsym) then
                    !if (rank == 0) print '(a,L,a,i20,a,b50)', 'is_det_allow', is_det_allow, ",i:", i, "bit(i)", i
                    ndet = ndet + 1
                    cas_idx0(ndet) = i
                    cas_idx_reverse(i) = ndet
                end if
            End if
        End if
    End do
    ! Stop the program if ndet == 0 because ndet == 0 means the number of CASCI determinant.
    if (ndet == 0) then
        if (rank == 0) then

            print *, "[ERROR]: The number of CASCI determinant is 0. Therefore, subsequent calculations  &
    &            cannot be performed successfully and the program is terminated."
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
!        print *,cas_idx(1:ndet)
    Deallocate (cas_idx0)

end subroutine casdet_ty

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE detsym_ty(ii, isym)

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

end subroutine detsym_ty
