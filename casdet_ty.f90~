! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE casdet_ty(totsym)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        use four_caspt2_module

        Implicit NONE

        integer, intent(in)   :: totsym

        integer :: nbitsa
        integer :: i, isym
        integer, allocatable  :: idet0(:)

        write(*,*)'Enter casdet_ty'
        Allocate(idet0(ndet))
        idet0 = 0
        ndet  = 0

        Do i = 1, 2**nact-1
           if(nbitsa(i) == nelec) then
              if(trim(ptgrp)=='C1') then
                 ndet = ndet +1
                 idet0(ndet) = i
              else
                 Call detsym_ty(i, isym)
                 if(isym == totsym) then
                 ndet = ndet +1
                 idet0(ndet) = i
              endif
              Endif
           Endif
        End do

        Allocate(idet(ndet))
        idet(1:ndet) = idet0(1:ndet)
        write(*,*)'totsym = ',totsym
        write(*,*)'ndet   = ',ndet
!        write(*,*)idet(1:ndet)
        Deallocate(idet0)

 1000   end subroutine casdet_ty

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE detsym_ty(ii,isym)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        use four_caspt2_module

        Implicit NONE

        integer, intent(in)  :: ii
        integer, intent(out) :: isym

        integer :: i, j, jsym, ielec, isym1

      isym = 1
        ielec = 0
        Do i = 1, nact
           if(btest(ii, i-1) .eqv..true.) then
              ielec = ielec + 1
              j = i+ninact
              jsym = irpamo(j)
              if(mod(ielec,2)==1) then
                 isym1 = MULTB_DS(jsym, isym) ! isym will be double irrep: odd number of electron
                 if(isym1 > nsymrp) write(*,*)'ielec, ii, isym, jsym, isym1',ielec, ii, isym, jsym+1, isym1
                 isym = isym1
              else
                 if(mod(jsym, 2)==1) then
                    isym1 = MULTB_D (jsym+1, isym) ! isym will be single irrep: even number of electron !MULTB_D is (fai*|fai)
                    if(isym1 > nsymrp) write(*,*)'ielec, ii, isym, jsym+1, isym1',ielec, ii, isym, jsym+1, isym1
                    isym = isym1
                 else
                    isym1 = MULTB_D (jsym-1, isym) ! isym will be single irrep: even number of electron
                    if(isym1 > nsymrp) write(*,*)'ielec, ii, isym, jsym-1, isym1',ielec, ii, isym, jsym-1, isym1
                    isym = isym1
                 endif
              endif

           Endif
        End do
        If(mod(ielec, 2) == 0) isym = isym + nsymrp ! even number electronic system

 1000   end subroutine detsym_ty


