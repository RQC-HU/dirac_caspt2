! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE readint2_ty(filename, nuniq) ! 2 electorn integrals created by typart in utchem

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE

    character*50, intent(in) :: filename

    character  :: datex*10, timex*8

    integer :: mdcint, nkr, idum, nuniq, nmom
    integer :: nz, type
    integer :: j0, i0, i1
    integer :: k0, l0, ii, jj, kk, ll, signind
    integer :: i, j, k, l, ikr, jkr, lkr, kkr
    integer :: SignIJ, SignKL, itr, jtr, ltr, ktr, inz, totalint

    integer, allocatable :: indk(:), indl(:), kr(:)

    real*8, allocatable :: rklr(:), rkli(:), int2rs(:), int2is(:)

    logical :: breit

    Allocate (int2rs(0:nmo**4)); Call memplus(KIND(int2rs), SIZE(int2rs), 1)
    Allocate (int2is(0:nmo**4)); Call memplus(KIND(int2is), SIZE(int2is), 1)

    Allocate (kr(-nmo/2:nmo/2)); Call memplus(KIND(kr), SIZE(kr), 1)
    Allocate (indtwr(nmo, nmo, nmo, nmo)); Call memplus(KIND(indtwr), SIZE(indtwr), 1)
    Allocate (indtwi(nmo, nmo, nmo, nmo)); Call memplus(KIND(indtwi), SIZE(indtwi), 1)

    kr = 0

    Allocate (indk((nmo/2)**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    Allocate (indl((nmo/2)**2)); Call memplus(KIND(indl), SIZE(indl), 1)
    Allocate (rklr((nmo/2)**2)); Call memplus(KIND(rklr), SIZE(rklr), 1)
    Allocate (rkli((nmo/2)**2)); Call memplus(KIND(rkli), SIZE(rkli), 1)

    write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024

    nuniq = 0
    indk(:) = 0
    indl(:) = 0
    rklr(:) = 0.0d+00
    rkli(:) = 0.0d+00
    int2r(:) = 0.0d+00
    int2i(:) = 0.0d+00
    indtwr = 0
    indtwi = 0

!###########################################################
!  THIS PART IS TAKEN FROM GOSCIP MOLFDIR PROGRAM PACKAGE
!###########################################################

    totalint = 0
    mdcint = 11
    open (mdcint, file=trim(filename), form='unformatted', status='unknown', err=10)

60  read (mdcint, ERR=40, END=50) i, j, nz, &
        (indk(inz), indl(inz), inz=1, nz), &
        (rklr(inz), rkli(inz), inz=1, nz)

    if (i == 0) goto 50

    totalint = totalint + nz

    itr = i + (-1)**(mod(i, 2) + 1)
    jtr = j + (-1)**(mod(j, 2) + 1)

    nmom = ninact + nact + nsec

    SignIJ = (-1)**(mod(i, 2) + mod(j, 2))
!                  SignIJ = SIGN(1,ikr) * SIGN(1,jkr)

    Do inz = 1, nz

        k = indk(inz)
        ktr = k + (-1)**(mod(k, 2) + 1)
        l = indl(inz)
        ltr = l + (-1)**(mod(l, 2) + 1)

        If (i > ninact + nact .and. j > ninact + nact .and. &
        &  k > ninact + nact .and. l > ninact + nact) goto 70

        SignKL = (-1)**(mod(k, 2) + mod(l, 2))
!                     SignKL = SIGN(1,kkr) * SIGN(1,lkr)
        nuniq = nuniq + 1

!=-> Original integral plus time-reversed partners
        INDTWR(I, J, K, L) = NUNIQ
        INDTWR(JTR, ITR, K, L) = NUNIQ*SignIJ
        INDTWR(I, J, LTR, KTR) = NUNIQ*SignKL
        INDTWR(JTR, ITR, LTR, KTR) = NUNIQ*SignIJ*SignKL
        INDTWI(I, J, K, L) = NUNIQ
        INDTWI(JTR, ITR, K, L) = NUNIQ*SignIJ
        INDTWI(I, J, LTR, KTR) = NUNIQ*SignKL
        INDTWI(JTR, ITR, LTR, KTR) = NUNIQ*SignIJ*SignKL
!=-> Complex conjugate plus time-reversed partners
        INDTWR(J, I, L, K) = NUNIQ
        INDTWR(ITR, JTR, L, K) = NUNIQ*SignIJ
        INDTWR(J, I, KTR, LTR) = NUNIQ*SignKL
        INDTWR(ITR, JTR, KTR, LTR) = NUNIQ*SignIJ*SignKL
        INDTWI(J, I, L, K) = -NUNIQ
        INDTWI(ITR, JTR, L, K) = -NUNIQ*SignIJ
        INDTWI(J, I, KTR, LTR) = -NUNIQ*SignKL
        INDTWI(ITR, JTR, KTR, LTR) = -NUNIQ*SignIJ*SignKL
!=-> Particle interchanged plus time-reversed partners
        INDTWR(K, L, I, J) = NUNIQ
        INDTWR(LTR, KTR, I, J) = NUNIQ*SignKL
        INDTWR(K, L, JTR, ITR) = NUNIQ*SignIJ
        INDTWR(LTR, KTR, JTR, ITR) = NUNIQ*SignIJ*SignKL
        INDTWI(K, L, I, J) = NUNIQ
        INDTWI(LTR, KTR, I, J) = NUNIQ*SignKL
        INDTWI(K, L, JTR, ITR) = NUNIQ*SignIJ
        INDTWI(LTR, KTR, JTR, ITR) = NUNIQ*SignIJ*SignKL
!=-> Particle interchanged and complex conjugated plus time-reversed partners
        INDTWR(L, K, J, I) = NUNIQ
        INDTWR(KTR, LTR, J, I) = NUNIQ*SignKL
        INDTWR(L, K, ITR, JTR) = NUNIQ*SignIJ
        INDTWR(KTR, LTR, ITR, JTR) = NUNIQ*SignIJ*SignKL
        INDTWI(L, K, J, I) = -NUNIQ
        INDTWI(KTR, LTR, J, I) = -NUNIQ*SignKL
        INDTWI(L, K, ITR, JTR) = -NUNIQ*SignIJ
        INDTWI(KTR, LTR, ITR, JTR) = -NUNIQ*SignIJ*SignKL

        int2rs(nuniq) = rklr(inz)
        int2is(nuniq) = rkli(inz)

!                     If(abs(rklr(inz))>1.0d-1) write(*,*)rklr(inz),rkli(inz), &
!                     & i, j, k, l

        if (abs(rkli(inz)) > thres) realc = .false.

5       FORMAT(4(4I3, 2I6))

70  End do

    indk(:) = 0
    indl(:) = 0
    rklr = 0.0d+00
    rkli = 0.0d+00

    Goto 60

10  write (*, *) 'error for opening mdcint 10'
    go to 100
20  write (*, *) 'error for reading mdcint 20'
    go to 100
30  write (*, *) 'end mdcint 30'
    go to 100
40  write (*, *) 'error for reading mdcint 40'
    go to 100
50  write (*, *) 'end mdcint 50 normal'
    go to 100

100 continue

    close (mdcint)
    write (*, *) nuniq, totalint

    Allocate (int2r(0:nuniq)); Call memplus(KIND(int2r), SIZE(int2r), 1)

    int2r(0:nuniq) = int2rs(0:nuniq)

    Deallocate (int2rs); Call memminus(KIND(int2rs), SIZE(int2rs), 1)

    Allocate (int2i(0:nuniq)); Call memplus(KIND(int2i), SIZE(int2i), 1)

    int2i(0:nuniq) = int2is(0:nuniq)

    Deallocate (int2is); Call memminus(KIND(int2is), SIZE(int2is), 1)

    deallocate (indk); Call memminus(KIND(indk), SIZE(indk), 1)
    deallocate (indl); Call memminus(KIND(indl), SIZE(indl), 1)
    deallocate (rklr); Call memminus(KIND(rklr), SIZE(rklr), 1)
    deallocate (rkli); Call memminus(KIND(rkli), SIZE(rkli), 1)
    deallocate (kr); Call memminus(KIND(kr), SIZE(kr), 1)

end subroutine readint2_ty
