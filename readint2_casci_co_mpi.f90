! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! SUBROUTINE readint2_casci_co(nuniq)  ! 2 electorn integrals created by typart in utchem
program readint2_casci_co  ! 2 electorn integrals created by typart in utchem

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module, &
        only: realonly, int2r_f1, int2i_f1, int2r_f2, int2i_f2, &
              indtwr, indtwi, int2r, int2i, tmem !, nmo, ninact, nact, nsec

    Implicit NONE

    character*50 :: filename

    character  :: datex*10, timex*8

    integer    :: mdcint, nkr, idum, nuniq, nmom, nmoc
!        integer ::  nkr, idum
!        integer    :: mdcint, nuniq, nmom, nmoc
    integer    :: nz, type
!        integer :: nz
!        integer    :: type
    integer    :: j0, i0, i1
    integer    :: k0, l0, ii, jj, kk, ll, signind
!        integer    :: i, j, k, l, jtr0, itr0
    integer    :: i, j, k, l
    integer    :: ikr, jkr, inz, kkr, lkr
    integer    :: jtr0, itr0

    integer    :: SignIJ, SignKL, itr, jtr, ltr, ktr, totalint, save, count

    complex*16 :: cint2

!        integer, allocatable :: indk(:), indl(:)
    integer, allocatable :: indk(:), indl(:), kr(:)
    real*8, allocatable  :: rklr(:), rkli(:), int2rs(:), int2is(:)
!        real*8, allocatable  :: int2rs(:), int2is(:)
!        double precision, allocatable  :: rklr(:), rkli(:)
    logical :: breit

    !@ comment : What is nuniq ? n unique
    ! For standalone mode. If you run whole casci/caspt2 code, comment
    ! out next block.
    real(8) :: thres = 1.0d-15
    logical :: realc = .false.
    integer :: ninact = 8, nact = 12, nsec = 172, nmo = 192
    integer, allocatable :: sp(:)

    integer :: ierr, MPI_COMM_WORLD, nprocs, rank ! variables for MPI

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_rank(MPI_COMM_WORLD, rank, ierr)

    ! For standalone mode. If you run whole casci/caspt2 code, comment
    ! out next block.
    Allocate (sp(1:nmo)); Call memplus(KIND(sp), SIZE(sp), 1)
    sp(1:ninact) = 1
    sp(ninact + 1:ninact + nact) = 2
    sp(ninact + nact + 1:ninact + nact + nsec) = 3
    sp(ninact + nact + nsec + 1:nmo) = 4
! Iwamuro modify
    realonly = .false.

    nmoc = ninact + nact
    nmom = ninact + nact + nsec

    Allocate (int2rs(0:nmoc**4)); Call memplus(KIND(int2rs), SIZE(int2rs), 1)
    Allocate (int2is(0:nmoc**4)); Call memplus(KIND(int2is), SIZE(int2is), 1)

    Allocate (int2r_f1(ninact + nact + 1:ninact + nact + nsec, ninact + nact + 1:ninact + nact + nsec, nmoc, nmoc))
    Allocate (int2i_f1(ninact + nact + 1:ninact + nact + nsec, ninact + nact + 1:ninact + nact + nsec, nmoc, nmoc))
    Allocate (int2r_f2(ninact + nact + 1:ninact + nact + nsec, nmoc, nmoc, ninact + nact + 1:ninact + nact + nsec))
    Allocate (int2i_f2(ninact + nact + 1:ninact + nact + nsec, nmoc, nmoc, ninact + nact + 1:ninact + nact + nsec))
    Call memplus(KIND(int2r_f1), SIZE(int2r_f1), 1)
    Call memplus(KIND(int2i_f1), SIZE(int2i_f1), 1)
    Call memplus(KIND(int2r_f2), SIZE(int2r_f2), 1)
    Call memplus(KIND(int2i_f2), SIZE(int2i_f2), 1)

    Allocate (indtwr(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(indtwr), SIZE(indtwr), 1)
    Allocate (indtwi(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(indtwi), SIZE(indtwi), 1)

!        Allocate(indk((nmo/2)**2)); Call memplus(KIND(indk),SIZE(indk),1)
!        Allocate(indl((nmo/2)**2)); Call memplus(KIND(indl),SIZE(indl),1)
!        Allocate(rklr((nmo/2)**2)); Call memplus(KIND(rklr),SIZE(rklr),1)
!        Allocate(rkli((nmo/2)**2)); Call memplus(KIND(rkli),SIZE(rkli),1)

    Allocate (indk(nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    Allocate (indl(nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    Allocate (rklr(nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    Allocate (rkli(nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)

!Iwamuro modify
    Allocate (kr(-nmo/2:nmo/2)); Call memplus(KIND(kr), SIZE(kr), 1)

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
    int2r_f1 = 0.0d+00
    int2i_f1 = 0.0d+00
    int2r_f2 = 0.0d+00
    int2i_f2 = 0.0d+00

    totalint = 0
    mdcint = 11
    open (mdcint, file=trim(filename), form='unformatted', status='old', err=10)

    read (mdcint, err=20, end=30) datex, timex, nkr, &
        (kr(i0), kr(-1*i0), i0=1, nkr)

    write (*, *) datex, timex
    write (*, *) 'nkr', nkr, 'kr(+),kr(-)', (kr(i0), kr(-1*i0), i0=1, nkr)

60  read (mdcint, END=50) i, j, nz, &
        (indk(inz), indl(inz), rklr(inz), rkli(inz), inz=1, nz)

!         do inz = 1,nz
!            write(*,'(I4)') i!,j,indk(inz),indl(inz),rklr(inz),rkli(inz)
!            write(*,'(5I4,2E20.10)') i,j,nz,indk(inz),indl(inz),rklr(inz),rkli(inz)
!         enddo

    if (i == 0) goto 50

    totalint = totalint + nz
    ! tr=trans クラマースペアの番号
    itr = i + (-1)**(mod(i, 2) + 1)
    jtr = j + (-1)**(mod(j, 2) + 1)

    i0 = i
    itr0 = itr
    j0 = j
    jtr0 = jtr
!Iwamuro debug
!        write(*,*) 'debug1'
    Do inz = 1, nz
        ! クラマースペアと、軌道エネルギー順
        ! MDCINTは正負の値でインデックスがついている
        i = i0
        itr = itr0
        j = j0
        jtr = jtr0

        k = indk(inz)
        ktr = k + (-1)**(mod(k, 2) + 1)
        l = indl(inz)
        ltr = l + (-1)**(mod(l, 2) + 1)
!Iwamuro debug
!        write(*,*) 'debug2'
!            write(*,*)sp(i),sp(j),sp(k),sp(l)
!            if(sp(l) == 0) write(*,*)i,j,k,l,'0'

        If (i > nmoc .and. j > nmoc .and. k > nmoc .and. l > nmoc) goto 70 ! (33|33) is ignored
!Iwamuro debug
!        write(*,*) 'debug3'
        If (i == j .and. k > l) goto 70
!Iwamuro debug
!        write(*,*) 'debug4', i, j, k, l
!            write(*,*)sp(i)!,sp(j),sp(k),sp(l)
        If (i <= nmoc .and. j <= nmoc .and. k <= nmoc .and. l <= nmoc) then
!               write(*,'("type 1",4I4,2E20.10)')i,j,k,l,rklr(inz),rkli(inz)
!Iwamuro debug
!        write(*,*) 'debug5'
            SignIJ = (-1)**(mod(i, 2) + mod(j, 2))
            SignKL = (-1)**(mod(k, 2) + mod(l, 2))
            ! nuniqをどうやってMPIでuniqueな番号にするのか
            nuniq = nuniq + 1
!Iwamuro debug
!        write(*,*) 'debug6'
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

!Iwamuro debug
!        write(*,*) 'debug7'
            if (abs(rkli(inz)) > thres) realc = .false.

        elseif (sp(i) == 3 .and. sp(j) == 3 .and. sp(k) < 3 .and. sp(l) == sp(k)) then !(33|11) or (33|22) type
!              write(*,'("type 2",4I4,2E20.10)')i,j,k,l,rklr(inz),rkli(inz)

            count = 0

11          if (mod(i, 2) == 0) then
                itr = i - 1
            else
                itr = i + 1
            end if

            if (mod(j, 2) == 0) then
                jtr = j - 1
            else
                jtr = j + 1
            end if

            if (mod(k, 2) == 0) then
                ktr = k - 1
            else
                ktr = k + 1
            end if

            if (mod(l, 2) == 0) then
                ltr = l - 1
            else
                ltr = l + 1
            end if

            SignIJ = (-1.0d+00)**mod(i + j, 2)
            SignKL = (-1.0d+00)**mod(k + l, 2)
!               write(*,*)'sign',signIJ,signKL

            int2r_f1(i, j, k, l) = rklr(inz)
            int2i_f1(i, j, k, l) = rkli(inz)

            int2r_f1(jtr, itr, k, l) = SignIJ*rklr(inz)
            int2i_f1(jtr, itr, k, l) = SignIJ*rkli(inz)

            int2r_f1(i, j, ltr, ktr) = SignKL*rklr(inz)
            int2i_f1(i, j, ltr, ktr) = SignKL*rkli(inz)

            int2r_f1(jtr, itr, ltr, ktr) = SignIJ*SignKL*rklr(inz)
            int2i_f1(jtr, itr, ltr, ktr) = SignIJ*SignKL*rkli(inz)

            count = count + 1
            cint2 = DCMPLX(rklr(inz), rkli(inz))

            if (count == 1) then
                Call takekr(i, j, k, l, cint2)              ! Consider Kramers pair
                rklr(inz) = DBLE(cint2)
                rkli(inz) = DIMAG(cint2)
                goto 11
            else
                goto 70
            end if

        elseif (sp(k) == 3 .and. sp(l) == 3 .and. sp(i) < 3 .and. sp(i) == sp(j)) then !(11|33) or (22|33) type
!               write(*,'("type 3",4I4,2E20.10)')i,j,k,l,rklr(inz),rkli(inz)

            count = 0

21          if (mod(i, 2) == 0) then
                itr = i - 1
            else
                itr = i + 1
            end if

            if (mod(j, 2) == 0) then
                jtr = j - 1
            else
                jtr = j + 1
            end if

            if (mod(k, 2) == 0) then
                ktr = k - 1
            else
                ktr = k + 1
            end if

            if (mod(l, 2) == 0) then
                ltr = l - 1
            else
                ltr = l + 1
            end if

            SignIJ = (-1.0d+00)**mod(i + j, 2)
            SignKL = (-1.0d+00)**mod(k + l, 2)
!               write(*,*)'sign',signIJ,signKL

            int2r_f1(k, l, i, j) = rklr(inz)
            int2i_f1(k, l, i, j) = rkli(inz)

            int2r_f1(k, l, jtr, itr) = SignIJ*rklr(inz)
            int2i_f1(k, l, jtr, itr) = SignIJ*rkli(inz)

            int2r_f1(ltr, ktr, i, j) = SignKL*rklr(inz)
            int2i_f1(ltr, ktr, i, j) = SignKL*rkli(inz)

            int2r_f1(ltr, ktr, jtr, itr) = SignIJ*SignKL*rklr(inz)
            int2i_f1(ltr, ktr, jtr, itr) = SignIJ*SignKL*rkli(inz)

            count = count + 1
            cint2 = DCMPLX(rklr(inz), rkli(inz))
            if (count == 1) then
                Call takekr(i, j, k, l, cint2)              ! Consider Kramers pair
                rklr(inz) = DBLE(cint2)
                rkli(inz) = DIMAG(cint2)
                goto 21
            else
                goto 70
            end if

        elseif (max(sp(i), sp(j)) == 3 .and. max(sp(k), sp(l)) == 3 .and. &
              &  min(sp(i), sp(j)) == min(sp(k), sp(l))) then                !(31|31) or (32|32) series

            count = 0

12          if (mod(i, 2) == 0) then
                itr = i - 1
            else
                itr = i + 1
            end if

            if (mod(j, 2) == 0) then
                jtr = j - 1
            else
                jtr = j + 1
            end if

            if (mod(k, 2) == 0) then
                ktr = k - 1
            else
                ktr = k + 1
            end if

            if (mod(l, 2) == 0) then
                ltr = l - 1
            else
                ltr = l + 1
            end if

            SignIJ = (-1.0d+00)**mod(i + j, 2)
            SignKL = (-1.0d+00)**mod(k + l, 2)

            if (i > j .and. k > l) then ! (31|31) or (32|32) ==> (31|13) or (32|23)

                int2r_f2(i, j, ltr, ktr) = signKL*rklr(inz)
                int2i_f2(i, j, ltr, ktr) = signKL*rkli(inz)

!                  write(*,*)i,j,ltr,ktr,int2r_f2(i,j,ltr,ktr),int2i_f2(i,j,ltr,ktr)

            elseif (i > j .and. k < l) then ! (31|13) or (32|23) ==> (31|13) or (32|23)

                int2r_f2(i, j, k, l) = rklr(inz)
                int2i_f2(i, j, k, l) = rkli(inz)

!                  write(*,*)i,j,k,l,int2r_f2(i,j,k,l),int2i_f2(i,j,k,l)

            elseif (i < j .and. k < l) then ! (13|13) or (23|23) ==> (31|13) or (32|23)

                int2r_f2(jtr, itr, k, l) = signIJ*rklr(inz)
                int2i_f2(jtr, itr, k, l) = signIJ*rkli(inz)

!                  write(*,*)jtr,itr,k,l,int2r_f2(jtr,itr,k,l),int2i_f2(jtr,itr,k,l)

            elseif (i < j .and. k > l) then ! (13|31) or (23|32) ==> (31|13) or (32|23)

                int2r_f2(jtr, itr, ltr, ktr) = signIJ*signKL*rklr(inz)
                int2i_f2(jtr, itr, ltr, ktr) = signIJ*signKL*rkli(inz)

!                  write(*,*)jtr,itr,ltr,ktr,int2r_f2(jtr,itr,ltr,ktr),int2i_f2(jtr,itr,ltr,ktr)

            end if

            count = count + 1
            cint2 = DCMPLX(rklr(inz), rkli(inz))
            if (count == 1 .or. count == 3) then
                Call takekr(i, j, k, l, cint2)              ! Consider Kramers pair
                rklr(inz) = DBLE(cint2)
                rkli(inz) = DIMAG(cint2)
                goto 12
            elseif (count == 2) then           ! variables exchange (AA|BB) => (BB|AA)
                save = i
                i = k
                k = save
                save = j
                j = l
                l = save
                goto 12
            else
                goto 70
            end if
        else
        end if

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
41  write (*, *) 'error for reading mdcint 41'
    go to 100
! 56      write(*,*)'error for reading mdcint 56'
!         go to 100

100 continue

    close (mdcint)

    write (*, *) nuniq, totalint

!         write(*,*) "debug1"

    Allocate (int2r(0:nuniq)); Call memplus(KIND(int2r), SIZE(int2r), 1)

    int2r(0:nuniq) = int2rs(0:nuniq)

!         write(*,*) "debug2"

    Deallocate (int2rs); Call memminus(KIND(int2rs), SIZE(int2rs), 1)

!         write(*,*) "debug3"

    Allocate (int2i(0:nuniq)); Call memplus(KIND(int2i), SIZE(int2i), 1)

    int2i(0:nuniq) = int2is(0:nuniq)

!         write(*,*) "debug4"

    Deallocate (int2is); Call memminus(KIND(int2is), SIZE(int2is), 1)

!         write(*,*) "debug5"

    deallocate (indk); Call memminus(KIND(indk), SIZE(indk), 1)
    deallocate (indl); Call memminus(KIND(indl), SIZE(indl), 1)
    deallocate (rklr); Call memminus(KIND(rklr), SIZE(rklr), 1)
    deallocate (rkli); Call memminus(KIND(rkli), SIZE(rkli), 1)
    deallocate (kr); Call memminus(KIND(kr), SIZE(kr), 1)
    ! For standalone mode. If you run whole casci/caspt2 code, comment
    ! out next block.
    Deallocate (sp)

    call MPI_Finalize(ierr)
!         write(*,*) "debug6"

end program readint2_casci_co
! end subroutine readint2_casci_co