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
    ! integer    :: nz, type
    !        integer :: nz
    integer    :: type
    integer    :: j0, i0, i1
    integer    :: k0, l0, ii, jj, kk, ll, signind
    !        integer    :: i, j, k, l, jtr0, itr0
    ! integer    :: i, j, k, l
    integer    :: k, l
    integer, allocatable    :: i(:), j(:), nz(:)
    integer    :: ikr, jkr, inz, kkr, lkr
    integer    :: jtr0, itr0

    integer    :: SignIJ, SignKL, itr, jtr, ltr, ktr, totalint, save, count

    complex*16 :: cint2

    !        integer, allocatable :: indk(:), indl(:)
    ! integer, allocatable :: indk(:), indl(:), kr(:)
    integer, allocatable :: indk(:, :), indl(:, :), kr(:)
    ! real*8, allocatable  :: rklr(:), rkli(:)
    real*8, allocatable  :: rklr(:, :), rkli(:, :)
    ! real(8), allocatable int2rs(:), int2is(:)
    !        real*8, allocatable  :: int2rs(:), int2is(:)
    !        double precision, allocatable  :: rklr(:), rkli(:)
    logical :: breit, continue_read = .true.
    integer :: idx, readmax = 1000
    ! Iwamuro modify
    realonly = .false.

    nmoc = ninact + nact
    nmom = ninact + nact + nsec
    if (rank == 0) then ! Process limits for output
        write (normaloutput, '(A,I8)') "Enter readint2_casci_co", rank
    end if
    ! Allocate (int2rs(0:nmoc**4)); Call memplus(KIND(int2rs), SIZE(int2rs), 1)
    ! Allocate (int2is(0:nmoc**4)); Call memplus(KIND(int2is), SIZE(int2is), 1)
    allocate (i(1000)); call memplus(kind(i), size(i), 1)
    allocate (j(1000)); call memplus(kind(j), size(j), 1)
    allocate (nz(1000)); call memplus(kind(nz), size(nz), 1)
    Allocate (int2r_f1(ninact + nact + 1:ninact + nact + nsec, ninact + nact + 1:ninact + nact + nsec, nmoc, nmoc))
    Allocate (int2i_f1(ninact + nact + 1:ninact + nact + nsec, ninact + nact + 1:ninact + nact + nsec, nmoc, nmoc))
    Allocate (int2r_f2(ninact + nact + 1:ninact + nact + nsec, nmoc, nmoc, ninact + nact + 1:ninact + nact + nsec))
    Allocate (int2i_f2(ninact + nact + 1:ninact + nact + nsec, nmoc, nmoc, ninact + nact + 1:ninact + nact + nsec))
    Call memplus(KIND(int2r_f1), SIZE(int2r_f1), 1)
    Call memplus(KIND(int2i_f1), SIZE(int2i_f1), 1)
    Call memplus(KIND(int2r_f2), SIZE(int2r_f2), 1)
    Call memplus(KIND(int2i_f2), SIZE(int2i_f2), 1)

    ! Allocate (indtwr(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(indtwr), SIZE(indtwr), 1)
    ! Allocate (indtwi(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(indtwi), SIZE(indtwi), 1)
    Allocate (inttwr(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(inttwr), SIZE(inttwr), 1)
    Allocate (inttwi(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(inttwi), SIZE(inttwi), 1)

    !        Allocate(indk((nmo/2)**2)); Call memplus(KIND(indk),SIZE(indk),1)
    !        Allocate(indl((nmo/2)**2)); Call memplus(KIND(indl),SIZE(indl),1)
    !        Allocate(rklr((nmo/2)**2)); Call memplus(KIND(rklr),SIZE(rklr),1)
    !        Allocate(rkli((nmo/2)**2)); Call memplus(KIND(rkli),SIZE(rkli),1)

    ! Allocate (indk(nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    ! Allocate (indl(nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    ! Allocate (rklr(nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    ! Allocate (rkli(nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    Allocate (indk(1000, nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    Allocate (indl(1000, nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    Allocate (rklr(1000, nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    Allocate (rkli(1000, nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    !Iwamuro modify
    Allocate (kr(-nmo/2:nmo/2)); Call memplus(KIND(kr), SIZE(kr), 1)

    if (rank == 0) then ! Process limits for output
        write (normaloutput, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024
    end if
    nuniq = 0
    i(:) = 0
    j(:) = 0
    indk(:, :) = 0
    indl(:, :) = 0
    rklr(:, :) = 0.0d+00
    rkli(:, :) = 0.0d+00
    ! int2rs(:) = 0.0d+00
    ! int2is(:) = 0.0d+00
    ! indtwr = 0
    ! indtwi = 0
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

    ! if (rank == 0) write (normaloutput, *) 'read'
60  do idx = 1, readmax
        ! if (rank == 0) write(normaloutput, *) 'idx',idx
        read (mdcint, END=81) i(idx), j(idx), nz(idx), &
            (indk(idx, inz), indl(idx, inz), rklr(idx, inz), rkli(idx, inz), inz=1, nz(idx))
    end do
    goto 82
81  continue_read = .false.
    readmax = idx
    !         do inz = 1,nz(idx)
    !            write(*,'(I4)') i!,j,indk(idx,inz),indl(idx,inz),rklr(idx,inz),rkli(idx,inz)
    !            write(*,'(5I4,2E20.10)') i,j,nz(idx),indk(idx,inz),indl(idx,inz),rklr(idx,inz),rkli(idx,inz)
    !         enddo
82  continue
    do idx = 1, readmax
        if (i(idx) == 0) goto 50

        totalint = totalint + nz(idx)

        itr = i(idx) + (-1)**(mod(i(idx), 2) + 1)
        jtr = j(idx) + (-1)**(mod(j(idx), 2) + 1)

        i0 = i(idx)
        itr0 = itr
        j0 = j(idx)
        jtr0 = jtr
        !Iwamuro debug
        !        write(*,*) 'debug1'
        Do inz = 1, nz(idx)

            i(idx) = i0
            itr = itr0
            j(idx) = j0
            jtr = jtr0

            k = indk(idx, inz)
            ktr = k + (-1)**(mod(k, 2) + 1)
            l = indl(idx, inz)
            ltr = l + (-1)**(mod(l, 2) + 1)
            !Iwamuro debug
            !        write(*,*) 'debug2'
            !            write(*,*)sp(i(idx)),sp(j(idx)),sp(k),sp(l)
            !            if(sp(l) == 0) write(*,*)i(idx),j(idx),k,l,'0'

            If (i(idx) > nmoc .and. j(idx) > nmoc .and. k > nmoc .and. l > nmoc) goto 70 ! (33|33) is ignored
            !Iwamuro debug
            !        write(*,*) 'debug3'
            If (i(idx) == j(idx) .and. k > l) goto 70
            !Iwamuro debug
            !        write(*,*) 'debug4', i(idx), j(idx), k, l
            !            write(*,*)sp(i(idx))!,sp(j(idx)),sp(k),sp(l)
            If (i(idx) <= nmoc .and. j(idx) <= nmoc .and. k <= nmoc .and. l <= nmoc) then
                !               write(*,'("type 1",4I4,2E20.10)')i(idx),j(idx),k,l,rklr(idx,inz),rkli(idx,inz)
                !Iwamuro debug
                !        write(*,*) 'debug5'
                SignIJ = (-1)**(mod(i(idx), 2) + mod(j(idx), 2))
                SignKL = (-1)**(mod(k, 2) + mod(l, 2))
                nuniq = nuniq + 1
                !Iwamuro debug
                !        write(*,*) 'debug6'
                !=-> Original integral plus time-reversed partners
                INTTWR(I(idx), J(idx), K, L) = rklr(idx, inz)
                INTTWR(JTR, ITR, K, L) = rklr(idx, inz)*SignIJ
                INTTWR(I(idx), J(idx), LTR, KTR) = rklr(idx, inz)*SignKL
                INTTWR(JTR, ITR, LTR, KTR) = rklr(idx, inz)*SignIJ*SignKL
                INTTWI(I(idx), J(idx), K, L) = rkli(idx, inz)
                INTTWI(JTR, ITR, K, L) = rkli(idx, inz)*SignIJ
                INTTWI(I(idx), J(idx), LTR, KTR) = rkli(idx, inz)*SignKL
                INTTWI(JTR, ITR, LTR, KTR) = rkli(idx, inz)*SignIJ*SignKL
                !=-> Complex conjugate plus time-reversed partners
                INTTWR(J(idx), I(idx), L, K) = rklr(idx, inz)
                INTTWR(ITR, JTR, L, K) = rklr(idx, inz)*SignIJ
                INTTWR(J(idx), I(idx), KTR, LTR) = rklr(idx, inz)*SignKL
                INTTWR(ITR, JTR, KTR, LTR) = rklr(idx, inz)*SignIJ*SignKL
                INTTWI(J(idx), I(idx), L, K) = -rkli(idx, inz)
                INTTWI(ITR, JTR, L, K) = -rkli(idx, inz)*SignIJ
                INTTWI(J(idx), I(idx), KTR, LTR) = -rkli(idx, inz)*SignKL
                INTTWI(ITR, JTR, KTR, LTR) = -rkli(idx, inz)*SignIJ*SignKL
                !=-> Particle interchanged plus time-reversed partners
                INTTWR(K, L, I(idx), J(idx)) = rklr(idx, inz)
                INTTWR(LTR, KTR, I(idx), J(idx)) = rklr(idx, inz)*SignKL
                INTTWR(K, L, JTR, ITR) = rklr(idx, inz)*SignIJ
                INTTWR(LTR, KTR, JTR, ITR) = rklr(idx, inz)*SignIJ*SignKL
                INTTWI(K, L, I(idx), J(idx)) = rkli(idx, inz)
                INTTWI(LTR, KTR, I(idx), J(idx)) = rkli(idx, inz)*SignKL
                INTTWI(K, L, JTR, ITR) = rkli(idx, inz)*SignIJ
                INTTWI(LTR, KTR, JTR, ITR) = rkli(idx, inz)*SignIJ*SignKL
                !=-> Particle interchanged and complex conjugated plus time-reversed partners
                INTTWR(L, K, J(idx), I(idx)) = rklr(idx, inz)
                INTTWR(KTR, LTR, J(idx), I(idx)) = rklr(idx, inz)*SignKL
                INTTWR(L, K, ITR, JTR) = rklr(idx, inz)*SignIJ
                INTTWR(KTR, LTR, ITR, JTR) = rklr(idx, inz)*SignIJ*SignKL
                INTTWI(L, K, J(idx), I(idx)) = -rkli(idx, inz)
                INTTWI(KTR, LTR, J(idx), I(idx)) = -rkli(idx, inz)*SignKL
                INTTWI(L, K, ITR, JTR) = -rkli(idx, inz)*SignIJ
                INTTWI(KTR, LTR, ITR, JTR) = -rkli(idx, inz)*SignIJ*SignKL

                ! int2rs(nuniq) = rklr(idx,inz)
                ! int2is(nuniq) = rkli(idx,inz)

                !Iwamuro debug
                !        write(*,*) 'debug7'
                if (abs(rkli(idx, inz)) > thres) realc = .false.

            elseif (sp(i(idx)) == 3 .and. sp(j(idx)) == 3 .and. sp(k) < 3 .and. sp(l) == sp(k)) then !(33|11) or (33|22) type
                !              write(*,'("type 2",4I4,2E20.10)')i(idx),j(idx),k,l,rklr(idx,inz),rkli(idx,inz)

                count = 0

11              if (mod(i(idx), 2) == 0) then
                    itr = i(idx) - 1
                else
                    itr = i(idx) + 1
                end if

                if (mod(j(idx), 2) == 0) then
                    jtr = j(idx) - 1
                else
                    jtr = j(idx) + 1
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

                SignIJ = (-1.0d+00)**mod(i(idx) + j(idx), 2)
                SignKL = (-1.0d+00)**mod(k + l, 2)
                !               write(*,*)'sign',signIJ,signKL

                int2r_f1(i(idx), j(idx), k, l) = rklr(idx, inz)
                int2i_f1(i(idx), j(idx), k, l) = rkli(idx, inz)

                int2r_f1(jtr, itr, k, l) = SignIJ*rklr(idx, inz)
                int2i_f1(jtr, itr, k, l) = SignIJ*rkli(idx, inz)

                int2r_f1(i(idx), j(idx), ltr, ktr) = SignKL*rklr(idx, inz)
                int2i_f1(i(idx), j(idx), ltr, ktr) = SignKL*rkli(idx, inz)

                int2r_f1(jtr, itr, ltr, ktr) = SignIJ*SignKL*rklr(idx, inz)
                int2i_f1(jtr, itr, ltr, ktr) = SignIJ*SignKL*rkli(idx, inz)

                count = count + 1
                cint2 = DCMPLX(rklr(idx, inz), rkli(idx, inz))

                if (count == 1) then
                    Call takekr(i(idx), j(idx), k, l, cint2)              ! Consider Kramers pair
                    rklr(idx, inz) = DBLE(cint2)
                    rkli(idx, inz) = DIMAG(cint2)
                    goto 11
                else
                    goto 70
                end if

            elseif (sp(k) == 3 .and. sp(l) == 3 .and. sp(i(idx)) < 3 .and. sp(i(idx)) == sp(j(idx))) then !(11|33) or (22|33) type
                !               write(*,'("type 3",4I4,2E20.10)')i(idx),j(idx),k,l,rklr(idx,inz),rkli(idx,inz)

                count = 0

21              if (mod(i(idx), 2) == 0) then
                    itr = i(idx) - 1
                else
                    itr = i(idx) + 1
                end if

                if (mod(j(idx), 2) == 0) then
                    jtr = j(idx) - 1
                else
                    jtr = j(idx) + 1
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

                SignIJ = (-1.0d+00)**mod(i(idx) + j(idx), 2)
                SignKL = (-1.0d+00)**mod(k + l, 2)
                !               write(*,*)'sign',signIJ,signKL

                int2r_f1(k, l, i(idx), j(idx)) = rklr(idx, inz)
                int2i_f1(k, l, i(idx), j(idx)) = rkli(idx, inz)

                int2r_f1(k, l, jtr, itr) = SignIJ*rklr(idx, inz)
                int2i_f1(k, l, jtr, itr) = SignIJ*rkli(idx, inz)

                int2r_f1(ltr, ktr, i(idx), j(idx)) = SignKL*rklr(idx, inz)
                int2i_f1(ltr, ktr, i(idx), j(idx)) = SignKL*rkli(idx, inz)

                int2r_f1(ltr, ktr, jtr, itr) = SignIJ*SignKL*rklr(idx, inz)
                int2i_f1(ltr, ktr, jtr, itr) = SignIJ*SignKL*rkli(idx, inz)

                count = count + 1
                cint2 = DCMPLX(rklr(idx, inz), rkli(idx, inz))
                if (count == 1) then
                    Call takekr(i(idx), j(idx), k, l, cint2)              ! Consider Kramers pair
                    rklr(idx, inz) = DBLE(cint2)
                    rkli(idx, inz) = DIMAG(cint2)
                    goto 21
                else
                    goto 70
                end if

            elseif (max(sp(i(idx)), sp(j(idx))) == 3 .and. max(sp(k), sp(l)) == 3 .and. &
                  &  min(sp(i(idx)), sp(j(idx))) == min(sp(k), sp(l))) then                !(31|31) or (32|32) series

                count = 0

12              if (mod(i(idx), 2) == 0) then
                    itr = i(idx) - 1
                else
                    itr = i(idx) + 1
                end if

                if (mod(j(idx), 2) == 0) then
                    jtr = j(idx) - 1
                else
                    jtr = j(idx) + 1
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

                SignIJ = (-1.0d+00)**mod(i(idx) + j(idx), 2)
                SignKL = (-1.0d+00)**mod(k + l, 2)

                if (i(idx) > j(idx) .and. k > l) then ! (31|31) or (32|32) ==> (31|13) or (32|23)

                    int2r_f2(i(idx), j(idx), ltr, ktr) = signKL*rklr(idx, inz)
                    int2i_f2(i(idx), j(idx), ltr, ktr) = signKL*rkli(idx, inz)

                    !                  write(*,*)i(idx),j(idx),ltr,ktr,int2r_f2(i(idx),j(idx),ltr,ktr),int2i_f2(i(idx),j(idx),ltr,ktr)

                elseif (i(idx) > j(idx) .and. k < l) then ! (31|13) or (32|23) ==> (31|13) or (32|23)

                    int2r_f2(i(idx), j(idx), k, l) = rklr(idx, inz)
                    int2i_f2(i(idx), j(idx), k, l) = rkli(idx, inz)

                    !                  write(*,*)i(idx),j(idx),k,l,int2r_f2(i(idx),j(idx),k,l),int2i_f2(i(idx),j(idx),k,l)

                elseif (i(idx) < j(idx) .and. k < l) then ! (13|13) or (23|23) ==> (31|13) or (32|23)

                    int2r_f2(jtr, itr, k, l) = signIJ*rklr(idx, inz)
                    int2i_f2(jtr, itr, k, l) = signIJ*rkli(idx, inz)

                    !                  write(*,*)jtr,itr,k,l,int2r_f2(jtr,itr,k,l),int2i_f2(jtr,itr,k,l)

                elseif (i(idx) < j(idx) .and. k > l) then ! (13|31) or (23|32) ==> (31|13) or (32|23)

                    int2r_f2(jtr, itr, ltr, ktr) = signIJ*signKL*rklr(idx, inz)
                    int2i_f2(jtr, itr, ltr, ktr) = signIJ*signKL*rkli(idx, inz)

                    !                  write(*,*)jtr,itr,ltr,ktr,int2r_f2(jtr,itr,ltr,ktr),int2i_f2(jtr,itr,ltr,ktr)

                end if

                count = count + 1
                cint2 = DCMPLX(rklr(idx, inz), rkli(idx, inz))
                if (count == 1 .or. count == 3) then
                    Call takekr(i(idx), j(idx), k, l, cint2)              ! Consider Kramers pair
                    rklr(idx, inz) = DBLE(cint2)
                    rkli(idx, inz) = DIMAG(cint2)
                    goto 12
                elseif (count == 2) then           ! variables exchange (AA|BB) => (BB|AA)
                    save = i(idx)
                    i(idx) = k
                    k = save
                    save = j(idx)
                    j(idx) = l
                    l = save
                    goto 12
                else
                    goto 70
                end if
            else
            end if

70      End do
    end do

    i(:) = 0
    ! j(:) = 0
    ! nz(:) = 0
    ! indk(:, :) = 0
    ! indl(:, :) = 0
    ! rklr = 0.0d+00
    ! rkli = 0.0d+00
    if (continue_read) then
        Goto 60
    else
        goto 50
    end if

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

    ! Allocate (int2r(0:nuniq)); Call memplus(KIND(int2r), SIZE(int2r), 1)

    ! int2r(0:nuniq) = int2rs(0:nuniq)

    !         write(*,*) "debug2"

    ! Deallocate (int2rs); Call memminus(KIND(int2rs), SIZE(int2rs), 1)

    !         write(*,*) "debug3"

    ! Allocate (int2i(0:nuniq)); Call memplus(KIND(int2i), SIZE(int2i), 1)

    ! int2i(0:nuniq) = int2is(0:nuniq)

    !         write(*,*) "debug4"

    ! Deallocate (int2is); Call memminus(KIND(int2is), SIZE(int2is), 1)

    !         write(*,*) "debug5"

    deallocate (indk); Call memminus(KIND(indk), SIZE(indk), 1)
    deallocate (indl); Call memminus(KIND(indl), SIZE(indl), 1)
    deallocate (rklr); Call memminus(KIND(rklr), SIZE(rklr), 1)
    deallocate (rkli); Call memminus(KIND(rkli), SIZE(rkli), 1)
    deallocate (kr); Call memminus(KIND(kr), SIZE(kr), 1)
    deallocate (i); Call memminus(KIND(i), SIZE(i), 1)
    deallocate (j); Call memminus(KIND(j), SIZE(j), 1)
    deallocate (nz); call memminus(kind(nz), size(nz), 1)

    call MPI_Allreduce(MPI_IN_PLACE, inttwr(1, 1, 1, 1), nmoc**4, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, inttwi(1, 1, 1, 1), nmoc**4, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
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
