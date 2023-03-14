! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE casdet_ty

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module
    use ras_det_check
    Implicit NONE

    integer :: i, isym
    integer, allocatable  :: idet0(:)
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
    Allocate (idet0(ndet))
    Allocate (idetr(2**nact - 1)); call memplus(kind(idetr), size(idetr), 1)
    idet0 = 0
    idetr = 0
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
            if (trim(ptgrp) == 'C1') then
                ndet = ndet + 1
                idet0(ndet) = i
                idetr(i) = ndet
            else
                Call detsym_ty(i, isym)
!                if (rank == 0) print '(a,L,a,i20,a,b50)', 'is_det_allow', is_det_allow, ",i:", i, "bit(i)", i
                if (isym == totsym) then
                    ndet = ndet + 1
                    idet0(ndet) = i
                    idetr(i) = ndet
                end if
            End if
        End if
    End do

    Allocate (idet(ndet))
    idet(1:ndet) = idet0(1:ndet)
    if (rank == 0) then
        print *, 'allow  = ', allow_det_num
        print *, 'totsym = ', totsym
        print *, 'ndet   = ', ndet
    end if
!        print *,idet(1:ndet)
    Deallocate (idet0)

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
                if (isym1 > nsymrp .and. rank == 0) print *, 'ielec, ii, isym, jsym, isym1', ielec, ii, isym, jsym + 1, isym1
            else
                if (mod(jsym, 2) == 1) then
                    isym1 = MULTB_D(jsym + 1, isym) ! isym will be single irrep: even number of electron !MULTB_D is (fai*|fai)
                    isym = isym1
                    if (isym1 > nsymrp .and. rank == 0) print *, 'ielec, ii, isym, jsym+1, isym1', ielec, ii, isym, jsym + 1, isym1
                else
                    isym1 = MULTB_D(jsym - 1, isym) ! isym will be single irrep: even number of electron
                    isym = isym1
                    if (isym1 > nsymrp .and. rank == 0) print *, 'ielec, ii, isym, jsym-1, isym1', ielec, ii, isym, jsym - 1, isym1
                end if
            end if

        End if
    End do
    If (mod(ielec, 2) == 0) isym = isym + nsymrp ! even number electronic system

end subroutine detsym_ty
