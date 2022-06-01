! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE intra_1(spi, spj, spk, spl, fname)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer, intent(in)        :: spi, spj, spk, spl
    character(50), intent(in)    :: fname

    integer, allocatable :: indsym(:, :, :), nsym(:, :)
    complex*16, allocatable :: traint2(:, :, :, :)

    real*8                  :: thresd
    complex*16              :: cint2

    integer :: i, j, k, l, i1, j1, k1, l1, inew, jnew, knew, lnew
    integer :: ii, ji, ki, li, ie, je, ke, le
    integer :: inz, nmx, ini(3), end(3), isp, isym, imo, ired

    integer :: n_cnt
    logical :: is_opened
    thresd = 1.0d-15

    ini(1) = 1
    end(1) = ninact
    ini(2) = ninact + 1
    end(2) = ninact + nact
    ini(3) = ninact + nact + 1
    end(3) = ninact + nact + nsec

    nmx = max(ninact, nact, nsec)
    Allocate (indsym(3, nsymrpa, nmx)); Call memplus(KIND(indsym), SIZE(indsym), 1)
    Allocate (nsym(3, nsymrpa)); Call memplus(KIND(nsym), SIZE(nsym), 1)
    indsym = 0
    nsym = 0

    Do isym = 1, nsymrpa
        Do isp = 1, 3                           ! i th space( inact or act or sec)
            Do imo = ini(isp), end(isp)
                if (irpamo(imo) == isym) then
                    nsym(isp, isym) = nsym(isp, isym) + 1
                    indsym(isp, isym, nsym(isp, isym)) = imo
                end if
            End do
        End do
    End do

    ii = ini(spi)
    ie = end(spi)
    ji = ini(spj)
    je = end(spj)
    ki = ini(spk)
    ke = end(spk)
    li = ini(spl)
    le = end(spl)

    Allocate (traint2(ii:ie, ji:je, ki:ke, li:le)); Call memplus(KIND(traint2), SIZE(traint2), 2)

!debug
!        write(*,'("C1int",8I4)')ii,ie,ji,je,ki,ke,li,leE

    traint2 = 0.0d+00
    if (rank == 0) then ! Process limits for output
        write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024
    end if

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and first index transformation !
!                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='old', form='unformatted')  ! no symmetry about spi,spj,spk,spl

30  read (1, err=10, end=20) i, j, k, l, cint2

    isym = irpmo(l)

    Do lnew = 1, nsym(spl, isym)
        l1 = indsym(spl, isym, lnew)
        traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*f(l, l1)
    End do

    Call takekr(i, j, k, l, cint2)
    isym = irpmo(l)

    Do lnew = 1, nsym(spl, isym)
        l1 = indsym(spl, isym, lnew)
        traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*f(l, l1)
    End do

    goto 30 ! Continue to read 2-integrals

20  close (1)

#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, traint2(ii, ji, ki, li), (ie - ii + 1)*(je - ji + 1)*(ke - ki + 1)*(le - li + 1), &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='replace', form='unformatted')
    call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2(ii:ie, ji:je, ki:ke, li:le), thresd)
    close (1)

    traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and second index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='old', form='unformatted')
31  read (1, err=10, end=21) i, j, k, l, cint2

    isym = irpmo(k)

    Do knew = 1, nsym(spk, isym)
        k1 = indsym(spk, isym, knew)
        traint2(i, j, k1, l) = traint2(i, j, k1, l) + cint2*DCONJG(f(k, k1))
    End do

    goto 31 ! Continue to read integrals

21  close (1)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, traint2(ii, ji, ki, li), (ie - ii + 1)*(je - ji + 1)*(ke - ki + 1)*(le - li + 1), &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='replace', form='unformatted')
    call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2(ii:ie, ji:je, ki:ke, li:le), thresd)
    close (1)

    traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and third index transformation  !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='old', form='unformatted')

32  read (1, err=10, end=22) i, j, k, l, cint2

    isym = irpmo(j)

    Do jnew = 1, nsym(spj, isym)
        j1 = indsym(spj, isym, jnew)
        traint2(i, j1, k, l) = traint2(i, j1, k, l) + cint2*f(j, j1)
    End do

    goto 32 ! Continue to read 2-integrals

22  close (1)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, traint2(ii, ji, ki, li), (ie - ii + 1)*(je - ji + 1)*(ke - ki + 1)*(le - li + 1), &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='replace', form='unformatted')
    call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2(ii:ie, ji:je, ki:ke, li:le), thresd)
    close (1)

    traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and fourth index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='old', form='unformatted')

33  read (1, err=10, end=23) i, j, k, l, cint2

    isym = irpmo(i)

    Do inew = 1, nsym(spi, isym)
        i1 = indsym(spi, isym, inew)
        traint2(i1, j, k, l) = traint2(i1, j, k, l) + cint2*DCONJG(f(i, i1))
    End do

    goto 33 ! Continue to read 2-integrals

23  close (1)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, traint2(ii, ji, ki, li), (ie - ii + 1)*(je - ji + 1)*(ke - ki + 1)*(le - li + 1), &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='replace', form='unformatted')
    call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2(ii:ie, ji:je, ki:ke, li:le), thresd)
    close (1)

    goto 100
10  write (*, *) 'error opening file'
    inquire (1, opened=is_opened)
    if (is_opened .eqv. .true.) then
        close (1)
    end if

100 continue
    deallocate (traint2); Call memminus(KIND(traint2), SIZE(traint2), 2)

    deallocate (indsym); Call memminus(KIND(indsym), SIZE(indsym), 1)
    deallocate (nsym); Call memminus(KIND(nsym), SIZE(nsym), 1)

end subroutine intra_1

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE intra_2(spi, spj, spk, spl, fname)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer, intent(in)        :: spi, spj, spk, spl
    character(50), intent(in)    :: fname

    integer, allocatable :: indsym(:, :, :), nsym(:, :)
    complex*16, allocatable :: traint2(:, :, :, :)

    real*8                  :: thresd
    complex*16              :: cint2

    integer :: i, j, k, l
    integer ::i1, j1, k1, l1, inew, jnew, knew, lnew
    integer :: ii, ji, ki, li, ie, je, ke, le, lkr0, kkr0
    integer :: inz, nmx, ini(3), end(3), isp, isym, imo, ired, save
    integer :: n_cnt

    logical :: is_opened
    thresd = 1.0d-15

    if (.not. (spi == spk .and. spj == spl)) then
        write (*, *) 'error intra_2', spi, spj, spk, spl
        stop
    end if

    ini(1) = 1
    end(1) = ninact
    ini(2) = ninact + 1
    end(2) = ninact + nact
    ini(3) = ninact + nact + 1
    end(3) = ninact + nact + nsec

    nmx = max(ninact, nact, nsec)
    Allocate (indsym(3, nsymrpa, nmx)); Call memplus(KIND(indsym), SIZE(indsym), 1)
    Allocate (nsym(3, nsymrpa)); Call memplus(KIND(nsym), SIZE(nsym), 1)
    indsym = 0
    nsym = 0

    Do isym = 1, nsymrpa
        Do isp = 1, 3                           ! i th space( inact or act or sec)
            Do imo = ini(isp), end(isp)
                if (irpamo(imo) == isym) then
                    nsym(isp, isym) = nsym(isp, isym) + 1
                    indsym(isp, isym, nsym(isp, isym)) = imo
                end if
            End do
        End do
    End do

    ii = ini(spi)
    ie = end(spi)
    ji = ini(spj)
    je = end(spj)
    ki = ini(spk)
    ke = end(spk)
    li = ini(spl)
    le = end(spl)

    Allocate (traint2(ii:ie, ji:je, ki:ke, li:le)); Call memplus(KIND(traint2), SIZE(traint2), 2)

!        write(*,*)ii,ie,ji,je,ki,ke,li,le

    traint2 = 0.0d+00
    if (rank == 0) then ! Process limits for output
        write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024
    end if

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and first index transformation !
!                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='old', form='unformatted')
30  read (1, err=10, end=20) i, j, k, l, cint2

    isym = irpmo(l)

    Do lnew = 1, nsym(spl, isym)
        l1 = indsym(spl, isym, lnew)
        traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*f(l, l1)
    End do

    if (i == k .and. j == l) goto 50
    if (ABS(i - k) == 1 .and. ABS(j - l) == 1 .and. &
    & ABS(i/2 - k/2) == 1 .and. ABS(j/2 - l/2) == 1) goto 50

    ! swap indices i = k, j = l, k = i, l = j
    save = i
    i = k
    k = save
    save = j
    j = l
    l = save
    isym = irpmo(l)

    Do lnew = 1, nsym(spl, isym)
        l1 = indsym(spl, isym, lnew)
        traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*f(l, l1)
    End do

50  Call takekr(i, j, k, l, cint2)
    isym = irpmo(l)

    Do lnew = 1, nsym(spl, isym)
        l1 = indsym(spl, isym, lnew)
        traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*f(l, l1)
    End do

    if (i == k .and. j == l) goto 30 ! Continue to read 2-integrals
    if (ABS(i - k) == 1 .and. ABS(j - l) == 1 .and. &
    & ABS(i/2 - k/2) == 1 .and. ABS(j/2 - l/2) == 1) goto 30 ! Continue to read 2-integrals

    ! swap indecis i = k, j = l, k = i, l = j
    save = i
    i = k
    k = save
    save = j
    j = l
    l = save

    isym = irpmo(l)

    Do lnew = 1, nsym(spl, isym)
        l1 = indsym(spl, isym, lnew)
        traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*f(l, l1)
    End do

    goto 30 ! Continue to read 2-integrals

20  close (1)

#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, traint2(ii, ji, ki, li), (ie - ii + 1)*(je - ji + 1)*(ke - ki + 1)*(le - li + 1), &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='replace', form='unformatted')
    call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2(ii:ie, ji:je, ki:ke, li:le), thresd)
    close (1)

    traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and second index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='old', form='unformatted')
31  read (1, err=10, end=21) i, j, k, l, cint2

    isym = irpmo(k)

    Do knew = 1, nsym(spk, isym)
        k1 = indsym(spk, isym, knew)
        traint2(i, j, k1, l) = traint2(i, j, k1, l) + cint2*DCONJG(f(k, k1))
    End do

    goto 31 ! Continue to read 2-integrals

21  close (1)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, traint2(ii, ji, ki, li), (ie - ii + 1)*(je - ji + 1)*(ke - ki + 1)*(le - li + 1), &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='replace', form='unformatted')
    call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2(ii:ie, ji:je, ki:ke, li:le), thresd)
    close (1)

    traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and third index transformation  !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='old', form='unformatted')
32  read (1, err=10, end=22) i, j, k, l, cint2

    isym = irpmo(j)

    Do jnew = 1, nsym(spj, isym)
        j1 = indsym(spj, isym, jnew)
        traint2(i, j1, k, l) = traint2(i, j1, k, l) + cint2*f(j, j1)
    End do

    goto 32

22  close (1)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, traint2(ii, ji, ki, li), (ie - ii + 1)*(je - ji + 1)*(ke - ki + 1)*(le - li + 1), &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='replace', form='unformatted')
    call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2(ii:ie, ji:je, ki:ke, li:le), thresd)
    close (1)

    traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and fourth index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='old', form='unformatted')
33  read (1, err=10, end=23) i, j, k, l, cint2

    isym = irpmo(i)

    Do inew = 1, nsym(spi, isym)
        i1 = indsym(spi, isym, inew)
        traint2(i1, j, k, l) = traint2(i1, j, k, l) + cint2*DCONJG(f(i, i1))
    End do

    goto 33 ! Continue to read 2-integrals

23  close (1)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, traint2(ii, ji, ki, li), (ie - ii + 1)*(je - ji + 1)*(ke - ki + 1)*(le - li + 1), &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='replace', form='unformatted')
    call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2(ii:ie, ji:je, ki:ke, li:le), thresd)
    close (1)

    goto 100
10  write (*, *) 'error opening file'
    inquire (1, opened=is_opened)
    if (is_opened .eqv. .true.) then
        close (1)
    end if

100 continue
    deallocate (traint2); Call memminus(KIND(traint2), SIZE(traint2), 2)

    deallocate (indsym); Call memminus(KIND(indsym), SIZE(indsym), 1)
    deallocate (nsym); Call memminus(KIND(nsym), SIZE(nsym), 1)

end subroutine intra_2

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE intra_3(spi, spj, spk, spl, fname)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer, intent(in)        :: spi, spj, spk, spl
    character(50), intent(in)    :: fname

    integer, allocatable    :: indsym(:, :, :), nsym(:, :)
    complex*16, allocatable :: traint2(:, :, :, :)

    real*8                  :: thresd
    complex*16              :: cint2

    integer :: i, j, k, l
    integer :: initial_i, initial_j, initial_k, initial_l
    integer :: i1, j1, k1, l1, inew, jnew, knew, lnew
    integer :: ii, ji, ki, li, ie, je, ke, le, lkr0, kkr0
    integer :: inz, nmx, ini(3), end(3), isp, isym, imo, ired, save
    integer :: n_cnt

    logical :: is_opened
    thresd = 1.0d-15

    if (.not. (spk == spl)) then
        write (*, *) 'error intra_3', spi, spj, spk, spl
        stop
    end if
    ini(1) = 1
    end(1) = ninact
    ini(2) = ninact + 1
    end(2) = ninact + nact
    ini(3) = ninact + nact + 1
    end(3) = ninact + nact + nsec

    nmx = max(ninact, nact, nsec)
    Allocate (indsym(3, nsymrpa, nmx)); Call memplus(KIND(indsym), SIZE(indsym), 1)
    Allocate (nsym(3, nsymrpa)); Call memplus(KIND(nsym), SIZE(nsym), 1)
    indsym = 0
    nsym = 0

    Do isym = 1, nsymrpa
        Do isp = 1, 3                           ! i th space( inact or act or sec)
            Do imo = ini(isp), end(isp)
                if (irpamo(imo) == isym) then
                    nsym(isp, isym) = nsym(isp, isym) + 1
                    indsym(isp, isym, nsym(isp, isym)) = imo
                end if
            End do
        End do
    End do

    ii = ini(spi)
    ie = end(spi)
    ji = ini(spj)
    je = end(spj)
    ki = ini(spk)
    ke = end(spk)
    li = ini(spl)
    le = end(spl)

    Allocate (traint2(ii:ie, ji:je, ki:ke, li:le)); Call memplus(KIND(traint2), SIZE(traint2), 2)

!modify iwamuro
!        write(*,'("intra_3",8I4)')ii,ie,ji,je,ki,ke,li,le

    traint2 = 0.0d+00
    if (rank == 0) then ! Process limits for output
        write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024
    end if

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and first index transformation !
!                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='old', form='unformatted')
30  read (1, err=10, end=20) i, j, k, l, cint2
    ! save initial indices i,j,k,l to initial_i,initial_j,initial_k,initial_l, respectively.
    initial_i = i
    initial_j = j
    initial_k = k
    initial_l = l

    isym = irpamo(l)

    Do lnew = 1, nsym(spl, isym)
        l1 = indsym(spl, isym, lnew)
        traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*f(l, l1)
    End do

    Call takekr(i, j, k, l, cint2)
    isym = irpamo(l)
!Iwamuro modify
!           write(*,'("takekr",4I4,2E20.10)')i,j,k,l,cint2

    Do lnew = 1, nsym(spl, isym)
        l1 = indsym(spl, isym, lnew)
        traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*f(l, l1)
    End do

    if (ABS(k - l) == 1 .and. ABS(k/2 - l/2) == 1) goto 30 ! Continue to read 2-integrals

    i = initial_i
    j = initial_j
    if (mod(initial_k, 2) == 0) l = initial_k - 1
    if (mod(initial_k, 2) == 1) l = initial_k + 1
    if (mod(initial_l, 2) == 0) k = initial_l - 1
    if (mod(initial_l, 2) == 1) k = initial_l + 1
    cint2 = (-1.0d+00)**mod(initial_k + initial_l, 2)*cint2
    isym = irpamo(l)

    Do lnew = 1, nsym(spl, isym)
        l1 = indsym(spl, isym, lnew)
        traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*f(l, l1)
    End do

    Call takekr(i, j, k, l, cint2)
    isym = irpamo(l)

    Do lnew = 1, nsym(spl, isym)
        l1 = indsym(spl, isym, lnew)
        traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*f(l, l1)
    End do

    goto 30 ! Continue to read 2-integrals

20  close (1)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, traint2(ii, ji, ki, li), (ie - ii + 1)*(je - ji + 1)*(ke - ki + 1)*(le - li + 1), &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='replace', form='unformatted')
    call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2(ii:ie, ji:je, ki:ke, li:le), thresd)
    close (1)

    traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and second index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (rank == 0) then ! Process limits for output
        write (*, *) 'Read intergals  and second index transformation'
    end if
    open (1, file=trim(fname), status='old', form='unformatted')
31  read (1, err=10, end=21) i, j, k, l, cint2

    isym = irpamo(k)

    Do knew = 1, nsym(spk, isym)
        k1 = indsym(spk, isym, knew)
        traint2(i, j, k1, l) = traint2(i, j, k1, l) + cint2*DCONJG(f(k, k1))
    End do

    goto 31

21  close (1)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, traint2(ii, ji, ki, li), (ie - ii + 1)*(je - ji + 1)*(ke - ki + 1)*(le - li + 1), &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    open (1, file=trim(fname), status='replace', form='unformatted')
    call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2(ii:ie, ji:je, ki:ke, li:le), thresd)
    close (1)

    traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and third index transformation  !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='old', form='unformatted')
32  read (1, err=10, end=22) i, j, k, l, cint2

    isym = irpamo(j)

    Do jnew = 1, nsym(spj, isym)
        j1 = indsym(spj, isym, jnew)
        traint2(i, j1, k, l) = traint2(i, j1, k, l) + cint2*f(j, j1)
    End do

    goto 32

22  close (1)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, traint2(ii, ji, ki, li), (ie - ii + 1)*(je - ji + 1)*(ke - ki + 1)*(le - li + 1), &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='replace', form='unformatted')
    n_cnt = 0
    Do l = li, le
        Do k = ki, ke
            Do j = ji, je
                Do i = ii, ie
                    if (ABS(traint2(i, j, k, l)) > thresd) then
                        if (mod(n_cnt, nprocs) == rank) then
                            write (1) i, j, k, l, traint2(i, j, k, l)
                        end if
                        n_cnt = n_cnt + 1
                        ! write (1, '(4I4, 2e2.1)') i, j, k, l, traint2(i, j, k, l)
                    end if
                End do
            End do
        End do
    End do

    close (1)

    traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and fourth index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='old', form='unformatted')
33  read (1, err=10, end=23) i, j, k, l, cint2

    isym = irpamo(i)

    Do inew = 1, nsym(spi, isym)
        i1 = indsym(spi, isym, inew)
        traint2(i1, j, k, l) = traint2(i1, j, k, l) + cint2*DCONJG(f(i, i1))
    End do

    goto 33

23  close (1)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, traint2(ii, ji, ki, li), (ie - ii + 1)*(je - ji + 1)*(ke - ki + 1)*(le - li + 1), &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    open (1, file=trim(fname), status='replace', form='unformatted')
    call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2(ii:ie, ji:je, ki:ke, li:le), thresd)
    goto 100 ! No error in intra3

10  write (*, *) 'error opening file', rank
    inquire (1, opened=is_opened)
    if (is_opened .eqv. .true.) then
        close (1)
    end if
    goto 101

100 continue
    if (rank == 0) write (*, *) 'read and write file properly. filename : ', trim(fname)
101 continue
    deallocate (traint2); Call memminus(KIND(traint2), SIZE(traint2), 2)

    deallocate (indsym); Call memminus(KIND(indsym), SIZE(indsym), 1)
    deallocate (nsym); Call memminus(KIND(nsym), SIZE(nsym), 1)

end subroutine intra_3

subroutine write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, thresd)
    use four_caspt2_module, only: nprocs, rank
    implicit none
    integer                 :: n_cnt, i, j, k, l
    integer, intent(in)     :: ii, ie, ji, je, ki, ke, li, le
    real(8)                 :: thresd
    complex*16, intent(in)  :: traint2(ii:ie, ji:je, ki:ke, li:le)
    ! 4重ループを1重ループに変換する方法
    ! 元の記述 : do l=li,le;do k=ki,ke;do j=ji,je;do i=ii,ie;
    !-> newle=le-li+1;newke=ke-ki+1;newje=je-ji+1;newie=ie-ii+1
    !-> ここまでの処理でdo li,leがdo 1,newleに変換済み
    !-> i=mod(idx                      ,newie)+ii
    !-> j=mod(idx/newie                ,newje)+ji
    !-> k=mod(idx/(newie*newje)        ,newke)+ki
    !-> l=mod(idx/(newie*newje*newke)  ,newle)+li
    !-> これでwrite(1) i,j,k,l,traint2(i,j,k,l)の記述のまま同じ処理ができる
    !-> ほとんど速くならなかったのでrevertした

    n_cnt = 0
    Do l = li, le
        Do k = ki, ke
            Do j = ji, je
                Do i = ii, ie
                    if (ABS(traint2(i, j, k, l)) > thresd) then
                        if (mod(n_cnt, nprocs) == rank) then ! Averaging the size of the subspace 2-integral file per a MPI process
                            write (1) i, j, k, l, traint2(i, j, k, l)
                        end if
                        ! if traint2(i,j,k,l)>thresd, all MPI process need to count up n_cnt!!!
                        n_cnt = n_cnt + 1
                    end if
                End do
            End do
        End do
    End do
end subroutine write_traint2_to_disk
