! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM r4dcaspt2_tra_co   ! DO CASPT2 CALC WITH MO TRANSFORMATION

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
    integer                 :: ii, jj, kk, ll, typetype, i0
    integer                 ::  j, i, k, l, nuniq
    integer                 :: k0, l0, nint, n, dimn, n0, n1, nspace(3, 3)
    integer                 ::  totsym, inisym, endsym, ieshift

    logical                 :: test, cutoff

    real*8                  :: i2r, i2i, dr, di, nsign, e0, e2, e2all, weight0
    complex*16              ::  cmplxint, dens, trace1, trace2
    complex*16, allocatable  :: ci(:)
    real*8, allocatable      :: ecas(:)

    character*50            :: filename

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!        debug = .TRUE.
    debug = .FALSE.
    thres = 1.0d-15
!        thres = 0.0d+00

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

    write (*, *) ''
    write (*, *) ' ENTER R4DCASPT2_TRA_TY PROGRAM written by M. Abe 2007.7.23'
    write (*, *) ''
    tmem = 0.0d+00

    val = 0
    Call DATE_AND_TIME(VALUES=val)
    Write (*, *) 'Year = ', val(1), 'Mon = ', val(2), 'Date = ', val(3)
    Write (*, *) 'Hour = ', val(5), 'Min = ', val(6), 'Sec = ', val(7), '.', val(8)

    totalsec = val(8)*(1.0d-03) + val(7) + val(6)*(6.0d+01) + val(5)*(6.0d+01)**2
    initdate = val(3)
    inittime = totalsec

    write (*, *) inittime

    Call timing(val(3), totalsec, date0, tsec)

    eshift = 0.0d+00
    ieshift = 0

    open (5, file='active.inp', form='formatted', status='old')
    read (5, '(I4)') ninact
    read (5, '(I4)') nact
    read (5, '(I4)') nsec
    read (5, '(I4)') nelec
    read (5, '(I4)') nroot
    read (5, '(I4)') selectroot
    read (5, '(I4)') totsym
    read (5, '(I4)') ncore
    read (5, '(I4)') nbas
    read (5, '(E8.2)') eshift
    read (5, '(A6)') ptgrp
    close (5)

    write (*, *) 'ninact     =', ninact
    write (*, *) 'nact       =', nact
    write (*, *) 'nsec       =', nsec
    write (*, *) 'nelec      =', nelec
    write (*, *) 'nroot      =', nroot
    write (*, *) 'selectroot =', selectroot
    write (*, *) 'totsym     =', totsym
    write (*, *) 'ncore      =', ncore
    write (*, *) 'nbas       =', nbas
    write (*, *) 'eshift     =', eshift
    write (*, *) 'ptgrp      =', ptgrp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    If (mod(nelec, 2) == 0) then
        evenelec = .true.
    else
        evenelec = .false.
    End if

!        write(*,*)' ENTER READ r4dmoint1'
    write (*, *) ' ENTER READ MRCONEE'

    filename = 'MRCONEE'

    call readorb_enesym_co(filename)
    call read1mo_co(filename)

!       write(*,*)' EXIT READ r4dmoint1'
    write (*, *) ' EXIT READ MRCONEE'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write (*, *) ' ENTER READ MDCINT'

    filename = 'MDCINTNEW'

    Call readint2_ord_co(filename)

    write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nmo = ninact + nact + nsec
    write (*, *) 'nmo        =', nmo

    open (10, file='CIMAT', form='unformatted', status='old')

    read (10) ndet
    Allocate (idet(1:ndet)); Call memplus(KIND(idet), SIZE(idet), 1)
    Allocate (ecas(1:ndet)); Call memplus(KIND(ecas), SIZE(ecas), 1)

    read (10) idet(1:ndet)
    read (10) ecas(1:ndet)

    close (10)

    Allocate (eigen(1:nroot)); Call memplus(KIND(eigen), SIZE(eigen), 1)
    eigen = 0.0d+00
    eigen(1:nroot) = ecas(1:nroot) + ecore

    Deallocate (ecas)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write (*, *) ' ENTER READ NEWCICOEFF', ndet

    Allocate (ci(1:ndet))
    ci = 0.0d+00

    open (10, file='NEWCICOEFF', form='unformatted', status='old')

    read (10) ci(1:ndet)

    close (10)

    Allocate (cir(1:ndet, selectroot:selectroot))
    Allocate (cii(1:ndet, selectroot:selectroot))

    cir(1:ndet, selectroot) = DBLE(ci(1:ndet))
    cii(1:ndet, selectroot) = DIMAG(ci(1:ndet))

!        Do i0 = 1, ndet
!           write(*,'(2E20.10)')cir(i0,selectroot),cii(i0,selectroot)
!        End do

!        Do i0 = 1, ndet
!           write(*,'(2E20.10)')ci(i0)
!        End do

    deallocate (ci)

!        write(*,*)cir(1:ndet,selectroot)
!        write(*,*)cii(1:ndet,selectroot)

    write (*, *) ' EXIT READ NEWCICOEFF'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    open (10, file='EPS', form='unformatted', status='old')

    read (10) nmo
    Allocate (eps(1:nmo)); Call memplus(KIND(eps), SIZE(eps), 1)
    eps = 0.0d+00
    read (10) eps(1:nmo)

    close (10)
!        Do i = 1, nmo
!           write(*,*)'eps(',i,')= ',eps(i)
!        Enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    open (10, file='TRANSFOCK', form='unformatted', status='old')

    read (10) nmo
    Allocate (f(nmo, nmo)); Call memplus(KIND(f), SIZE(f), 2)
    read (10) f

    close (10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write (*, *) ' '
    write (*, *) '*******************************'
    write (*, *) ' '
    write (*, *) 'IREP IS ', repna(totsym)
    write (*, *) ' '
    write (*, *) '*******************************'
    write (*, *) ' '

    realcvec = .TRUE.

!    This is test for bug fix about realc part

    write (*, *) realc, 'realc'
    write (*, *) realcvec, 'realcvec'

    realc = .FALSE.      !!!      realc =.TRUE.
    realcvec = .FALSE.   !!!      realcvec =.TRUE.

    write (*, *) 'FOR TEST WE DO (F,F)'
    write (*, *) realc, 'realc'
    write (*, *) realcvec, 'realcvec'

!!=============================================!
!                                              !
    iroot = selectroot                       !
!                                              !
!!=============================================!

!      write(*,*)'RECALCULATION OF CASCI ENERGY'
!      Call e0test_v2

    e2 = 0.0d+00

    Call calce0(e0)

    e2all = 0.0d+00

    date1 = initdate
    tsec1 = totalsec

    Call timing(date1, tsec1, date0, tsec0)
    date1 = date0
    tsec1 = tsec0

    write (*, *) 'Enter intra3 A1int'
    Call intra_3(2, 1, 2, 2, 'A1int')
    write (*, *) 'Enter intra3 A2int'
    Call intra_3(2, 1, 1, 1, 'A2int')
    sumc2local = 0.0d+00
    write (*, *) 'Enter solvA'
    Call solvA_ord_ty(e0, e2)
    e2all = e2all + e2
    write (*, *) e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_2(2, 1, 2, 1, 'Bint ')

    sumc2local = 0.0d+00
    Call solvB_ord_ty(e0, e2)
    e2all = e2all + e2
    write (*, *) e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_3(3, 2, 2, 2, 'C1int')
    Call intra_3(3, 2, 1, 1, 'C2int')
    Call intra_1(3, 1, 1, 2, 'C3int')
    sumc2local = 0.0d+00
    Call solvC_ord_ty(e0, e2)
    e2all = e2all + e2
!      write(*,*) 'e2c is not added at while'

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_3(3, 1, 2, 2, 'D1int')
    Call intra_1(3, 2, 2, 1, 'D2int')
    Call intra_3(3, 1, 1, 1, 'D3int')
    sumc2local = 0.0d+00
    Call solvD_ord_ty(e0, e2)
    e2all = e2all + e2
    write (*, *) e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_1(3, 1, 2, 1, 'Eint')

    sumc2local = 0.0d+00
    Call solvE_ord_ty(e0, e2)
    e2all = e2all + e2
    write (*, *) e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_2(3, 2, 3, 2, 'Fint ')

    sumc2local = 0.0d+00
    Call solvF_ord_ty(e0, e2)
    e2all = e2all + e2
    write (*, *) e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_1(3, 1, 3, 2, 'Gint ')

    sumc2local = 0.0d+00
    Call solvG_ord_ty(e0, e2)
    e2all = e2all + e2
    write (*, *) e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    Call intra_2(3, 1, 3, 1, 'Hint ')

    sumc2local = 0.0d+00
    Call solvH_ord_ty(e0, e2)
    e2all = e2all + e2
    write (*, *) e2all

    date1 = date0
    tsec1 = tsec0
    Call timing(date1, tsec1, date0, tsec0)

    write (*, '("c^2 ",F30.15)') sumc2
    weight0 = 1.0d+00/(1.0d+00 + sumc2)
    write (*, '("weight of 0th wave function is",F30.15)') weight0

    write (*, '("Total second order energy is ",F30.15," a.u.")') e2all - eshift*sumc2
    write (*, '(" ")')
    write (*, '("Total energy is ",F30.15," a.u.")') e2all + eigen(iroot) - eshift*sumc2

    deallocate (cir); Call memminus(KIND(cir), SIZE(cir), 1)
    deallocate (cii); Call memminus(KIND(cii), SIZE(cii), 1)
    deallocate (eigen); Call memminus(KIND(eigen), SIZE(eigen), 1)
    deallocate (eps); Call memminus(KIND(eps), SIZE(eps), 1)
    deallocate (idet); Call memminus(KIND(idet), SIZE(idet), 1)

!    End do  ! totsym
    deallocate (sp); Call memminus(KIND(sp), SIZE(sp), 1)
    deallocate (orb); Call memminus(KIND(orb), SIZE(orb), 1)
    deallocate (irpmo); Call memminus(KIND(irpmo), SIZE(irpmo), 1)
    deallocate (irpamo); Call memminus(KIND(irpamo), SIZE(irpamo), 1)
    deallocate (indmo); Call memminus(KIND(indmo), SIZE(indmo), 1)
    deallocate (indmor); Call memminus(KIND(indmor), SIZE(indmor), 1)
    deallocate (oner); Call memminus(KIND(oner), SIZE(oner), 1)
    deallocate (onei); Call memminus(KIND(onei), SIZE(onei), 1)
    deallocate (MULTB_S); Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1)
    deallocate (MULTB_D); Call memminus(KIND(MULTB_D), SIZE(MULTB_D), 1)
    deallocate (MULTB_DS); Call memminus(KIND(MULTB_DS), SIZE(MULTB_DS), 1)
    deallocate (MULTB_DF); Call memminus(KIND(MULTB_DF), SIZE(MULTB_DF), 1)
    deallocate (MULTB_DB); Call memminus(KIND(MULTB_DB), SIZE(MULTB_DB), 1)
    deallocate (MULTB_SB); Call memminus(KIND(MULTB_SB), SIZE(MULTB_SB), 1)

    Call timing(val(3), totalsec, date0, tsec0)
    write (*, *) 'End r4dcaspt2_tra_ty'

1000 continue
END program r4dcaspt2_tra_co

