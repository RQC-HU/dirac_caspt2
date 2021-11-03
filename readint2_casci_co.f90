! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE readint2_casci_co(filename, nuniq)  ! 2 electorn integrals created by typart in utchem

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
    include 'mpif.h'
    character*50, intent(in) :: filename

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
! Iwamuro modify
    realonly = .false.

    nmoc = ninact + nact
    nmom = ninact + nact + nsec
    if (rank == 0) then ! Process limits for output
        write (normaloutput, '(A,I8)') "Enter readint2_casci_co", rank
    end if
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
    Allocate (inttwr(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(indtwr), SIZE(indtwr), 1)
    Allocate (inttwi(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(indtwi), SIZE(indtwi), 1)

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

    if (rank == 0) then ! Process limits for output
        write (normaloutput, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024
    end if

    nuniq = 0
    indk(:) = 0
    indl(:) = 0
    rklr(:) = 0.0d+00
    rkli(:) = 0.0d+00
    int2rs(:) = 0.0d+00
    int2is(:) = 0.0d+00
    indtwr = 0
    indtwi = 0
    inttwr = 0.0d+00
    inttwi = 0.0d+00
    int2r_f1 = 0.0d+00
    int2i_f1 = 0.0d+00
    int2r_f2 = 0.0d+00
    int2i_f2 = 0.0d+00

    totalint = 0
    mdcint = 11
    ! open (mdcint, file=trim(filename), form='unformatted', status='old', err=10)
    open (mdcint, file=trim(filename), form='unformatted', status='old')

    read (mdcint, err=20, end=30) datex, timex, nkr, &
        (kr(i0), kr(-1*i0), i0=1, nkr)

    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) datex, timex
        write (normaloutput, *) 'readint2', 'nkr', nkr, 'kr(+),kr(-)', (kr(i0), kr(-1*i0), i0=1, nkr)
    end if
60  read (mdcint, END=50) i, j, nz, &
        (indk(inz), indl(inz), rklr(inz), rkli(inz), inz=1, nz)

!         do inz = 1,nz
!            write(*,'(I4)') i!,j,indk(inz),indl(inz),rklr(inz),rkli(inz)
!            write(*,'(5I4,2E20.10)') i,j,nz,indk(inz),indl(inz),rklr(inz),rkli(inz)
!         enddo

    if (i == 0) goto 50

    totalint = totalint + nz

    itr = i + (-1)**(mod(i, 2) + 1)
    jtr = j + (-1)**(mod(j, 2) + 1)

    i0 = i
    itr0 = itr
    j0 = j
    jtr0 = jtr
!Iwamuro debug
!        write(*,*) 'debug1'
    Do inz = 1, nz

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
            nuniq = nuniq + 1
!Iwamuro debug
!        write(*,*) 'debug6'
            !=-> Original integral plus time-reversed partners
            INTTWR(I, J, K, L) = rklr(inz)
            INTTWR(JTR, ITR, K, L) = rklr(inz)*SignIJ
            INTTWR(I, J, LTR, KTR) = rklr(inz)*SignKL
            INTTWR(JTR, ITR, LTR, KTR) = rklr(inz)*SignIJ*SignKL
            INTTWI(I, J, K, L) = rkli(inz)
            INTTWI(JTR, ITR, K, L) = rkli(inz)*SignIJ
            INTTWI(I, J, LTR, KTR) = rkli(inz)*SignKL
            INTTWI(JTR, ITR, LTR, KTR) = rkli(inz)*SignIJ*SignKL
            !=-> Complex conjugate plus time-reversed partners
            INTTWR(J, I, L, K) = rklr(inz)
            INTTWR(ITR, JTR, L, K) = rklr(inz)*SignIJ
            INTTWR(J, I, KTR, LTR) = rklr(inz)*SignKL
            INTTWR(ITR, JTR, KTR, LTR) = rklr(inz)*SignIJ*SignKL
            INTTWI(J, I, L, K) = -rkli(inz)
            INTTWI(ITR, JTR, L, K) = -rkli(inz)*SignIJ
            INTTWI(J, I, KTR, LTR) = -rkli(inz)*SignKL
            INTTWI(ITR, JTR, KTR, LTR) = -rkli(inz)*SignIJ*SignKL
            !=-> Particle interchanged plus time-reversed partners
            INTTWR(K, L, I, J) = rklr(inz)
            INTTWR(LTR, KTR, I, J) = rklr(inz)*SignKL
            INTTWR(K, L, JTR, ITR) = rklr(inz)*SignIJ
            INTTWR(LTR, KTR, JTR, ITR) = rklr(inz)*SignIJ*SignKL
            INTTWI(K, L, I, J) = rkli(inz)
            INTTWI(LTR, KTR, I, J) = rkli(inz)*SignKL
            INTTWI(K, L, JTR, ITR) = rkli(inz)*SignIJ
            INTTWI(LTR, KTR, JTR, ITR) = rkli(inz)*SignIJ*SignKL
            !=-> Particle interchanged and complex conjugated plus time-reversed partners
            INTTWR(L, K, J, I) = rklr(inz)
            INTTWR(KTR, LTR, J, I) = rklr(inz)*SignKL
            INTTWR(L, K, ITR, JTR) = rklr(inz)*SignIJ
            INTTWR(KTR, LTR, ITR, JTR) = rklr(inz)*SignIJ*SignKL
            INTTWI(L, K, J, I) = -rkli(inz)
            INTTWI(KTR, LTR, J, I) = -rkli(inz)*SignKL
            INTTWI(L, K, ITR, JTR) = -rkli(inz)*SignIJ
            INTTWI(KTR, LTR, ITR, JTR) = -rkli(inz)*SignIJ*SignKL

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

    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'error for opening mdcint 10'
    end if
    go to 100
20  if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'error for reading mdcint 20'
    end if
    go to 100
30  if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'end mdcint 30'
    end if
    go to 100
40  if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'error for reading mdcint 40'
    end if
    go to 100
50  if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'end mdcint 50 normal'
    end if
    go to 100
41  if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'error for reading mdcint 41'
    end if
    go to 100
! 56      write(*,*)'error for reading mdcint 56'
!         go to 100

100 continue

    close (mdcint)
    if (rank == 0) then
        write (normaloutput, *) nuniq, totalint
    end if
!         write(*,*) "debug1"

    Allocate (int2r(0:nuniq)); Call memplus(KIND(int2r), SIZE(int2r), 1)

    ! int2r(0:nuniq) = int2rs(0:nuniq)

!         write(*,*) "debug2"

    Deallocate (int2rs); Call memminus(KIND(int2rs), SIZE(int2rs), 1)

!         write(*,*) "debug3"

    Allocate (int2i(0:nuniq)); Call memplus(KIND(int2i), SIZE(int2i), 1)

    ! int2i(0:nuniq) = int2is(0:nuniq)

!         write(*,*) "debug4"

    Deallocate (int2is); Call memminus(KIND(int2is), SIZE(int2is), 1)

!         write(*,*) "debug5"

    deallocate (indk); Call memminus(KIND(indk), SIZE(indk), 1)
    deallocate (indl); Call memminus(KIND(indl), SIZE(indl), 1)
    deallocate (rklr); Call memminus(KIND(rklr), SIZE(rklr), 1)
    deallocate (rkli); Call memminus(KIND(rkli), SIZE(rkli), 1)
    deallocate (kr); Call memminus(KIND(kr), SIZE(kr), 1)

    call MPI_Allreduce(MPI_IN_PLACE, inttwr(1, 1, 1, 1), &
                       nmoc**4, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, inttwi(1, 1, 1, 1), &
                       nmoc**4, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, int2r_f1(ninact + nact + 1, ninact + nact + 1, 1, 1), &
                       nsec*nsec*nmoc*nmoc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, int2i_f1(ninact + nact + 1, ninact + nact + 1, 1, 1), &
                       nsec*nsec*nmoc*nmoc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, int2r_f2(ninact + nact + 1, 1, 1, ninact + nact + 1), &
                       nsec*nmoc*nmoc*nsec, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, int2i_f2(ninact + nact + 1, 1, 1, ninact + nact + 1), &
                       nsec*nmoc*nmoc*nsec, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
!         write(*,*) "debug6"

end subroutine readint2_casci_co
