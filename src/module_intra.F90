module module_intra
! 2 electron integrals transformation module
    use module_takekr, only: takekr
    implicit none
    private
    public intra_1, intra_2, intra_3
    interface write_traint2_to_disk_fourth
        module procedure write_traint2_to_disk_fourth_complex
        module procedure write_traint2_to_disk_fourth_real
    end interface write_traint2_to_disk_fourth

    interface write_traint2_to_disk
        module procedure write_traint2_to_disk_complex
        module procedure write_traint2_to_disk_real
    end interface write_traint2_to_disk
contains

    subroutine intra_1(spi, spj, spk, spl, fname)
        use module_realonly, only: realonly
        implicit none
        integer, intent(in)        :: spi, spj, spk, spl
        character(len=*), intent(in)    :: fname
        if (realonly%is_realonly()) then
            call intra_1_real(spi, spj, spk, spl, fname)
        else
            call intra_1_complex(spi, spj, spk, spl, fname)
        end if
    end subroutine intra_1
    subroutine intra_2(spi, spj, spk, spl, fname)
        use module_realonly, only: realonly
        implicit none
        integer, intent(in)        :: spi, spj, spk, spl
        character(len=*), intent(in)    :: fname
        if (realonly%is_realonly()) then
            call intra_2_real(spi, spj, spk, spl, fname)
        else
            call intra_2_complex(spi, spj, spk, spl, fname)
        end if
    end subroutine intra_2
    subroutine intra_3(spi, spj, spk, spl, fname)
        use module_realonly, only: realonly
        implicit none
        integer, intent(in)        :: spi, spj, spk, spl
        character(len=*), intent(in)    :: fname
        if (realonly%is_realonly()) then
            call intra_3_real(spi, spj, spk, spl, fname)
        else
            call intra_3_complex(spi, spj, spk, spl, fname)
        end if
    end subroutine intra_3

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE intra_1_complex(spi, spj, spk, spl, fname)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_file_manager
        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)        :: spi, spj, spk, spl
        character(len=*), intent(in)    :: fname
        logical :: is_end_of_file

        integer                 :: unit_int2_subspace
        integer, allocatable    :: indsym(:, :, :), nsym(:, :)
        complex*16, allocatable :: traint2(:, :, :, :)
        complex*16 :: readint2

        real(8)                  :: cutoff_threshold
        complex*16              :: cint2

        integer :: i, j, k, l, i1, j1, k1, l1, inew, jnew, knew, lnew
        integer :: ii, ji, ki, li, ie, je, ke, le, iostat
        integer :: nmx, ini(3), end(3), isp, isym, imo

        cutoff_threshold = 1.0d-15
        readint2 = 0.0d+00
        ini(1) = 1
        end(1) = ninact
        ini(2) = global_act_start
        end(2) = global_act_end
        ini(3) = global_sec_start
        end(3) = global_sec_end

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

        traint2 = 0.0d+00
        call write_allocated_memory_size

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and first index transformation !
!                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = readint2
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_cmplx(l, l1)
            End do

            Call takekr(i, j, k, l, cint2)
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_cmplx(l, l1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)

        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and second index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = readint2
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            isym = irpamo(k)

            Do knew = 1, nsym(spk, isym)
                k1 = indsym(spk, isym, knew)
                traint2(i, j, k1, l) = traint2(i, j, k1, l) + cint2*DCONJG(fock_cmplx(k, k1))
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and third index transformation  !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = readint2
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            isym = irpamo(j)

            Do jnew = 1, nsym(spj, isym)
                j1 = indsym(spj, isym, jnew)
                traint2(i, j1, k, l) = traint2(i, j1, k, l) + cint2*fock_cmplx(j, j1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and fourth index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = readint2
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            isym = irpamo(i)

            Do inew = 1, nsym(spi, isym)
                i1 = indsym(spi, isym, inew)
                traint2(i1, j, k, l) = traint2(i1, j, k, l) + cint2*DCONJG(fock_cmplx(i, i1))
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk_fourth(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        Call memminus(KIND(traint2), SIZE(traint2), 2); deallocate (traint2)

        Call memminus(KIND(indsym), SIZE(indsym), 1); deallocate (indsym)
        Call memminus(KIND(nsym), SIZE(nsym), 1); deallocate (nsym)

    end subroutine intra_1_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE intra_1_real(spi, spj, spk, spl, fname)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_file_manager
        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)        :: spi, spj, spk, spl
        character(len=*), intent(in)    :: fname
        logical :: is_end_of_file

        integer                 :: unit_int2_subspace
        integer, allocatable    :: indsym(:, :, :), nsym(:, :)
        real(8), allocatable :: traint2(:, :, :, :)
        complex*16 :: readint2

        real(8)                  :: cutoff_threshold
        real(8)              :: cint2

        integer :: i, j, k, l, i1, j1, k1, l1, inew, jnew, knew, lnew
        integer :: ii, ji, ki, li, ie, je, ke, le, iostat
        integer :: nmx, ini(3), end(3), isp, isym, imo

        cutoff_threshold = 1.0d-15
        readint2 = 0.0d+00
        ini(1) = 1
        end(1) = ninact
        ini(2) = global_act_start
        end(2) = global_act_end
        ini(3) = global_sec_start
        end(3) = global_sec_end

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

        Allocate (traint2(ii:ie, ji:je, ki:ke, li:le)); Call memplus(KIND(traint2), SIZE(traint2), 1)

        traint2 = 0.0d+00
        call write_allocated_memory_size

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and first index transformation !
!                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = real(readint2, kind=KIND(cint2))
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_real(l, l1)
            End do

            Call takekr(i, j, k, l, cint2)
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_real(l, l1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)

        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and second index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = real(readint2, kind=KIND(cint2))
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            isym = irpamo(k)

            Do knew = 1, nsym(spk, isym)
                k1 = indsym(spk, isym, knew)
                traint2(i, j, k1, l) = traint2(i, j, k1, l) + cint2*fock_real(k, k1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and third index transformation  !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = real(readint2, kind=KIND(cint2))
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            isym = irpamo(j)

            Do jnew = 1, nsym(spj, isym)
                j1 = indsym(spj, isym, jnew)
                traint2(i, j1, k, l) = traint2(i, j1, k, l) + cint2*fock_real(j, j1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and fourth index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = real(readint2, kind=KIND(cint2))
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            isym = irpamo(i)

            Do inew = 1, nsym(spi, isym)
                i1 = indsym(spi, isym, inew)
                traint2(i1, j, k, l) = traint2(i1, j, k, l) + cint2*fock_real(i, i1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk_fourth(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        Call memminus(KIND(traint2), SIZE(traint2), 1); deallocate (traint2)

        Call memminus(KIND(indsym), SIZE(indsym), 1); deallocate (indsym)
        Call memminus(KIND(nsym), SIZE(nsym), 1); deallocate (nsym)

    end subroutine intra_1_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE intra_2_complex(spi, spj, spk, spl, fname)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_error, only: stop_with_errorcode
        use module_file_manager
        use module_global_variables
        use module_sort_swap, only: swap
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)        :: spi, spj, spk, spl
        character(len=*), intent(in)    :: fname
        logical :: is_end_of_file

        integer                 :: unit_int2_subspace
        integer, allocatable    :: indsym(:, :, :), nsym(:, :)
        complex*16, allocatable :: traint2(:, :, :, :)

        real(8)                  :: cutoff_threshold
        complex*16              :: cint2
        complex*16              :: readint2

        integer :: i, j, k, l
        integer :: i1, j1, k1, l1, inew, jnew, knew, lnew
        integer :: ii, ji, ki, li, ie, je, ke, le
        integer :: nmx, ini(3), end(3), isp, isym, imo, iostat

        cutoff_threshold = 1.0d-15
        readint2 = 0.0d+00

        if (.not. (spi == spk .and. spj == spl)) then
            print *, 'error intra_2', spi, spj, spk, spl
            call stop_with_errorcode(1)
        end if

        ini(1) = 1
        end(1) = ninact
        ini(2) = global_act_start
        end(2) = global_act_end
        ini(3) = global_sec_start
        end(3) = global_sec_end

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

        traint2 = 0.0d+00
        call write_allocated_memory_size

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and first index transformation !
!                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = readint2
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_cmplx(l, l1)
            End do

            if (i == k .and. j == l) then
                continue
            else if (ABS(i - k) == 1 .and. ABS(j - l) == 1 .and. ABS(i/2 - k/2) == 1 .and. ABS(j/2 - l/2) == 1) then
                continue
            else

                call swap(i, k)
                call swap(j, l)
                isym = irpamo(l)

                Do lnew = 1, nsym(spl, isym)
                    l1 = indsym(spl, isym, lnew)
                    traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_cmplx(l, l1)
                End do
            end if
            Call takekr(i, j, k, l, cint2)
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_cmplx(l, l1)
            End do

            if (i == k .and. j == l) cycle ! Continue to read 2-integrals
            if (ABS(i - k) == 1 .and. ABS(j - l) == 1 .and. &
            & ABS(i/2 - k/2) == 1 .and. ABS(j/2 - l/2) == 1) cycle ! Continue to read 2-integrals

            call swap(i, k)
            call swap(j, l)
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_cmplx(l, l1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and second index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do

            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = readint2
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(k)

            Do knew = 1, nsym(spk, isym)
                k1 = indsym(spk, isym, knew)
                traint2(i, j, k1, l) = traint2(i, j, k1, l) + cint2*DCONJG(fock_cmplx(k, k1))
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and third index transformation  !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = readint2
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(j)

            Do jnew = 1, nsym(spj, isym)
                j1 = indsym(spj, isym, jnew)
                traint2(i, j1, k, l) = traint2(i, j1, k, l) + cint2*fock_cmplx(j, j1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and fourth index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = readint2
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(i)

            Do inew = 1, nsym(spi, isym)
                i1 = indsym(spi, isym, inew)
                traint2(i1, j, k, l) = traint2(i1, j, k, l) + cint2*DCONJG(fock_cmplx(i, i1))
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk_fourth(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        Call memminus(KIND(traint2), SIZE(traint2), 2); deallocate (traint2)
        Call memminus(KIND(indsym), SIZE(indsym), 1); deallocate (indsym)
        Call memminus(KIND(nsym), SIZE(nsym), 1); deallocate (nsym)

    end subroutine intra_2_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE intra_2_real(spi, spj, spk, spl, fname)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_error, only: stop_with_errorcode
        use module_file_manager
        use module_global_variables
        use module_sort_swap, only: swap
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)        :: spi, spj, spk, spl
        character(len=*), intent(in)    :: fname
        logical :: is_end_of_file

        integer                 :: unit_int2_subspace
        integer, allocatable    :: indsym(:, :, :), nsym(:, :)
        real(8), allocatable :: traint2(:, :, :, :)

        real(8)                  :: cutoff_threshold
        real(8)              :: cint2
        complex*16              :: readint2

        integer :: i, j, k, l
        integer :: i1, j1, k1, l1, inew, jnew, knew, lnew
        integer :: ii, ji, ki, li, ie, je, ke, le
        integer :: nmx, ini(3), end(3), isp, isym, imo, iostat

        cutoff_threshold = 1.0d-15
        readint2 = 0.0d+00

        if (.not. (spi == spk .and. spj == spl)) then
            print *, 'error intra_2', spi, spj, spk, spl
            call stop_with_errorcode(1)
        end if

        ini(1) = 1
        end(1) = ninact
        ini(2) = global_act_start
        end(2) = global_act_end
        ini(3) = global_sec_start
        end(3) = global_sec_end

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

        Allocate (traint2(ii:ie, ji:je, ki:ke, li:le)); Call memplus(KIND(traint2), SIZE(traint2), 1)

        traint2 = 0.0d+00
        call write_allocated_memory_size

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and first index transformation !
!                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = real(readint2, kind=KIND(cint2))
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_real(l, l1)
            End do

            if (i == k .and. j == l) then
                continue
            else if (ABS(i - k) == 1 .and. ABS(j - l) == 1 .and. ABS(i/2 - k/2) == 1 .and. ABS(j/2 - l/2) == 1) then
                continue
            else

                call swap(i, k)
                call swap(j, l)
                isym = irpamo(l)

                Do lnew = 1, nsym(spl, isym)
                    l1 = indsym(spl, isym, lnew)
                    traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_real(l, l1)
                End do
            end if
            Call takekr(i, j, k, l, cint2)
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_real(l, l1)
            End do

            if (i == k .and. j == l) cycle ! Continue to read 2-integrals
            if (ABS(i - k) == 1 .and. ABS(j - l) == 1 .and. &
            & ABS(i/2 - k/2) == 1 .and. ABS(j/2 - l/2) == 1) cycle ! Continue to read 2-integrals

            call swap(i, k)
            call swap(j, l)
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_real(l, l1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and second index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do

            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = real(readint2, kind=KIND(cint2))
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(k)

            Do knew = 1, nsym(spk, isym)
                k1 = indsym(spk, isym, knew)
                traint2(i, j, k1, l) = traint2(i, j, k1, l) + cint2*fock_real(k, k1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and third index transformation  !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = real(readint2, kind=KIND(cint2))
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(j)

            Do jnew = 1, nsym(spj, isym)
                j1 = indsym(spj, isym, jnew)
                traint2(i, j1, k, l) = traint2(i, j1, k, l) + cint2*fock_real(j, j1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and fourth index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = real(readint2, kind=KIND(cint2))
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(i)

            Do inew = 1, nsym(spi, isym)
                i1 = indsym(spi, isym, inew)
                traint2(i1, j, k, l) = traint2(i1, j, k, l) + cint2*fock_real(i, i1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk_fourth(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        Call memminus(KIND(traint2), SIZE(traint2), 1); deallocate (traint2)
        Call memminus(KIND(indsym), SIZE(indsym), 1); deallocate (indsym)
        Call memminus(KIND(nsym), SIZE(nsym), 1); deallocate (nsym)

    end subroutine intra_2_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE intra_3_complex(spi, spj, spk, spl, fname)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_error, only: stop_with_errorcode
        use module_file_manager
        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)        :: spi, spj, spk, spl
        character(len=*), intent(in)    :: fname
        logical :: is_end_of_file

        integer                 :: unit_int2_subspace
        integer, allocatable    :: indsym(:, :, :), nsym(:, :)
        complex*16, allocatable :: traint2(:, :, :, :)
        complex*16              :: readint2

        real(8)                  :: cutoff_threshold
        complex*16              :: cint2, initial_cint2

        integer :: i, j, k, l
        integer :: initial_i, initial_j, initial_k, initial_l
        integer :: i1, j1, k1, l1, inew, jnew, knew, lnew
        integer :: ii, ji, ki, li, ie, je, ke, le
        integer :: nmx, ini(3), end(3), isp, isym, imo, iostat

        cutoff_threshold = 1.0d-15
        readint2 = 0.0d+00

        if (.not. (spk == spl)) then
            print *, 'error intra_3', spi, spj, spk, spl
            call stop_with_errorcode(1)
        end if
        ini(1) = 1
        end(1) = ninact
        ini(2) = global_act_start
        end(2) = global_act_end
        ini(3) = global_sec_start
        end(3) = global_sec_end

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
        traint2 = 0.0d+00
        call write_allocated_memory_size

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and first index transformation !
!                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = readint2
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

! save initial indices i,j,k,l,cint2 to initial_i,initial_j,initial_k,initial_l,initial_cint2 respectively.
            initial_i = i
            initial_j = j
            initial_k = k
            initial_l = l
            initial_cint2 = cint2

            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_cmplx(l, l1)
            End do

            Call takekr(i, j, k, l, cint2)
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_cmplx(l, l1)
            End do

            if (ABS(k - l) == 1 .and. ABS(k/2 - l/2) == 1) cycle ! Continue to read 2-integrals

            i = initial_i
            j = initial_j
            k = initial_l + merge(1, -1, mod(initial_l, 2)==1) ! k = initial_l + 1 if initial_l is odd, otherwise k = initial_l - 1
            l = initial_k + merge(1, -1, mod(initial_k, 2)==1) ! l = initial_k + 1 if initial_k is odd, otherwise l = initial_k - 1
            cint2 = initial_cint2*merge(1, -1, mod(initial_k + initial_l, 2)==0) ! cint2 = initial_cint2 if initial_k + initial_l is even, otherwise cint2 = -initial_cint2
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_cmplx(l, l1)
            End do

            Call takekr(i, j, k, l, cint2)
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_cmplx(l, l1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and second index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = readint2
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(k)

            Do knew = 1, nsym(spk, isym)
                k1 = indsym(spk, isym, knew)
                traint2(i, j, k1, l) = traint2(i, j, k1, l) + cint2*DCONJG(fock_cmplx(k, k1))
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and third index transformation  !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = readint2
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(j)

            Do jnew = 1, nsym(spj, isym)
                j1 = indsym(spj, isym, jnew)
                traint2(i, j1, k, l) = traint2(i, j1, k, l) + cint2*fock_cmplx(j, j1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and fourth index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = readint2
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(i)

            Do inew = 1, nsym(spi, isym)
                i1 = indsym(spi, isym, inew)
                traint2(i1, j, k, l) = traint2(i1, j, k, l) + cint2*DCONJG(fock_cmplx(i, i1))
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk_fourth(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)

        if (debug .and. rank == 0) print *, 'read and write file properly. filename : ', trim(fname)
        Call memminus(KIND(traint2), SIZE(traint2), 2); deallocate (traint2)
        Call memminus(KIND(indsym), SIZE(indsym), 1); deallocate (indsym)
        Call memminus(KIND(nsym), SIZE(nsym), 1); deallocate (nsym)

    end subroutine intra_3_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE intra_3_real(spi, spj, spk, spl, fname)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_error, only: stop_with_errorcode
        use module_file_manager
        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)        :: spi, spj, spk, spl
        character(len=*), intent(in)    :: fname
        logical :: is_end_of_file

        integer                 :: unit_int2_subspace
        integer, allocatable    :: indsym(:, :, :), nsym(:, :)
        real(8), allocatable :: traint2(:, :, :, :)
        complex*16              :: readint2

        real(8)                  :: cutoff_threshold
        real(8)              :: cint2, initial_cint2

        integer :: i, j, k, l
        integer :: initial_i, initial_j, initial_k, initial_l
        integer :: i1, j1, k1, l1, inew, jnew, knew, lnew
        integer :: ii, ji, ki, li, ie, je, ke, le
        integer :: nmx, ini(3), end(3), isp, isym, imo, iostat

        cutoff_threshold = 1.0d-15
        readint2 = 0.0d+00

        if (.not. (spk == spl)) then
            print *, 'error intra_3', spi, spj, spk, spl
            call stop_with_errorcode(1)
        end if
        ini(1) = 1
        end(1) = ninact
        ini(2) = global_act_start
        end(2) = global_act_end
        ini(3) = global_sec_start
        end(3) = global_sec_end

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

        Allocate (traint2(ii:ie, ji:je, ki:ke, li:le)); Call memplus(KIND(traint2), SIZE(traint2), 1)
        traint2 = 0.0d+00
        call write_allocated_memory_size

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and first index transformation !
!                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = real(readint2, kind=KIND(cint2))
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

! save initial indices i,j,k,l,cint2 to initial_i,initial_j,initial_k,initial_l,initial_cint2 respectively.
            initial_i = i
            initial_j = j
            initial_k = k
            initial_l = l
            initial_cint2 = cint2

            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_real(l, l1)
            End do

            Call takekr(i, j, k, l, cint2)
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_real(l, l1)
            End do

            if (ABS(k - l) == 1 .and. ABS(k/2 - l/2) == 1) cycle ! Continue to read 2-integrals

            i = initial_i
            j = initial_j
            k = initial_l + merge(1, -1, mod(initial_l, 2)==1) ! k = initial_l + 1 if initial_l is odd, otherwise k = initial_l - 1
            l = initial_k + merge(1, -1, mod(initial_k, 2)==1) ! l = initial_k + 1 if initial_k is odd, otherwise l = initial_k - 1
            cint2 = initial_cint2*merge(1, -1, mod(initial_k + initial_l, 2)==0) ! cint2 = initial_cint2 if initial_k + initial_l is even, otherwise cint2 = -initial_cint2
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_real(l, l1)
            End do

            Call takekr(i, j, k, l, cint2)
            isym = irpamo(l)

            Do lnew = 1, nsym(spl, isym)
                l1 = indsym(spl, isym, lnew)
                traint2(i, j, k, l1) = traint2(i, j, k, l1) + cint2*fock_real(l, l1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and second index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = real(readint2, kind=KIND(cint2))
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(k)

            Do knew = 1, nsym(spk, isym)
                k1 = indsym(spk, isym, knew)
                traint2(i, j, k1, l) = traint2(i, j, k1, l) + cint2*fock_real(k, k1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and third index transformation  !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = real(readint2, kind=KIND(cint2))
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(j)

            Do jnew = 1, nsym(spj, isym)
                j1 = indsym(spj, isym, jnew)
                traint2(i, j1, k, l) = traint2(i, j1, k, l) + cint2*fock_real(j, j1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)
        traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and fourth index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='read')
        do
            read (unit_int2_subspace, iostat=iostat) i, j, k, l, readint2
            cint2 = real(readint2, kind=KIND(cint2))
            call check_iostat(iostat=iostat, file=trim(fname), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            isym = irpamo(i)

            Do inew = 1, nsym(spi, isym)
                i1 = indsym(spi, isym, inew)
                traint2(i1, j, k, l) = traint2(i1, j, k, l) + cint2*fock_real(i, i1)
            End do
        end do
        close (unit_int2_subspace)

#ifdef HAVE_MPI
        call allreduce_wrapper(traint2)
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2_subspace, file=trim(fname), status='old', optional_action='write')
        call write_traint2_to_disk_fourth(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        close (unit_int2_subspace)

        if (debug .and. rank == 0) print *, 'read and write file properly. filename : ', trim(fname)
        Call memminus(KIND(traint2), SIZE(traint2), 1); deallocate (traint2)
        Call memminus(KIND(indsym), SIZE(indsym), 1); deallocate (indsym)
        Call memminus(KIND(nsym), SIZE(nsym), 1); deallocate (nsym)

    end subroutine intra_3_real

    subroutine write_traint2_to_disk_fourth_complex(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
!==============================================================================================
! This is a writing subroutine for two-electron integrals
! after the fourth integral transformation.
!
! For the clarity of later processing, this subroutine will write the numbers
! within each subspace orbitals,
! instead of the numbers of the entire inactive, active, and secondary orbitals.
!
! Author: K.Noda
!==============================================================================================
        use module_global_variables, only: nprocs, rank
        implicit none
        integer                 :: n_cnt, i, j, k, l
        integer, intent(in)     :: ii, ie, ji, je, ki, ke, li, le, unit_int2_subspace
        real(8)                 :: cutoff_threshold
        complex*16, intent(in)  :: traint2(ii:ie, ji:je, ki:ke, li:le)
        integer                 :: i_tra, j_tra, k_tra, l_tra
        integer                 :: i_ini, i_end, j_ini, j_end, k_ini, k_end, l_ini, l_end

! Initialization
        i_ini = ii; i_end = ie; j_ini = ji; j_end = je; k_ini = ki; k_end = ke; l_ini = li; l_end = le

! Check where to start and end the writing
        call where_subspace_is(i_ini, i_end)
        call where_subspace_is(j_ini, j_end)
        call where_subspace_is(k_ini, k_end)
        call where_subspace_is(l_ini, l_end)

! Write the integrals to disk
        n_cnt = 0
        i_tra = ii - 1
        j_tra = ji - 1
        k_tra = ki - 1
        l_tra = li - 1
        Do l = l_ini, l_end
            Do k = k_ini, k_end
                Do j = j_ini, j_end
                    Do i = i_ini, i_end
                        !===================================================================================================
                        ! About the index of traint2
                        ! If ii is in the active subspace, ii is equal to nact + 1,
                        ! so i + i_tra is needed to get the correct index. (The same applies to inactive and secondary.)
                        ! When nact is zero, it never pass through the do loop even once, so there is no problem.
                        !===================================================================================================
                        if (ABS(traint2(i + i_tra, j + j_tra, k + k_tra, l + l_tra)) > cutoff_threshold) then
                            if (mod(n_cnt, nprocs) == rank) then ! Averaging the size of the subspace 2-integral file per a MPI process
                                write (unit_int2_subspace) i, j, k, l, traint2(i + i_tra, j + j_tra, k + k_tra, l + l_tra)
                            end if
                            n_cnt = n_cnt + 1
                        end if
                    End do
                End do
            End do
        End do
    end subroutine write_traint2_to_disk_fourth_complex

    subroutine write_traint2_to_disk_complex(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        use module_global_variables, only: nprocs, rank
        implicit none
        integer                 :: n_cnt, i, j, k, l
        integer, intent(in)     :: ii, ie, ji, je, ki, ke, li, le, unit_int2_subspace
        real(8)                 :: cutoff_threshold
        complex*16, intent(in)  :: traint2(ii:ie, ji:je, ki:ke, li:le)

        n_cnt = 0
        Do l = li, le
            Do k = ki, ke
                Do j = ji, je
                    Do i = ii, ie
                        ! Very small integrals are not written to disk
                        if (ABS(traint2(i, j, k, l)) > cutoff_threshold) then
                            if (mod(n_cnt, nprocs) == rank) then ! Averaging the size of the subspace 2-integral file per a MPI process
                                write (unit_int2_subspace) i, j, k, l, traint2(i, j, k, l)
                            end if
                            ! if traint2(i,j,k,l)>cutoff_threshold, all MPI process need to count up n_cnt!!!
                            n_cnt = n_cnt + 1
                        end if
                    End do
                End do
            End do
        End do
    end subroutine write_traint2_to_disk_complex

    subroutine write_traint2_to_disk_fourth_real(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
!==============================================================================================
! This is a writing subroutine for two-electron integrals
! after the fourth integral transformation.
!
! For the clarity of later processing, this subroutine will write the numbers
! within each subspace orbitals,
! instead of the numbers of the entire inactive, active, and secondary orbitals.
!
! Author: K.Noda
!==============================================================================================
        use module_global_variables, only: nprocs, rank
        implicit none
        integer                 :: n_cnt, i, j, k, l
        integer, intent(in)     :: ii, ie, ji, je, ki, ke, li, le, unit_int2_subspace
        real(8)                 :: cutoff_threshold
        real(8), intent(in)  :: traint2(ii:ie, ji:je, ki:ke, li:le)
        integer                 :: i_tra, j_tra, k_tra, l_tra
        integer                 :: i_ini, i_end, j_ini, j_end, k_ini, k_end, l_ini, l_end

! Initialization
        i_ini = ii; i_end = ie; j_ini = ji; j_end = je; k_ini = ki; k_end = ke; l_ini = li; l_end = le

! Check where to start and end the writing
        call where_subspace_is(i_ini, i_end)
        call where_subspace_is(j_ini, j_end)
        call where_subspace_is(k_ini, k_end)
        call where_subspace_is(l_ini, l_end)

! Write the integrals to disk
        n_cnt = 0
        i_tra = ii - 1
        j_tra = ji - 1
        k_tra = ki - 1
        l_tra = li - 1
        Do l = l_ini, l_end
            Do k = k_ini, k_end
                Do j = j_ini, j_end
                    Do i = i_ini, i_end
                        !===================================================================================================
                        ! About the index of traint2
                        ! If ii is in the active subspace, ii is equal to nact + 1,
                        ! so i + i_tra is needed to get the correct index. (The same applies to inactive and secondary.)
                        ! When nact is zero, it never pass through the do loop even once, so there is no problem.
                        !===================================================================================================
                        if (ABS(traint2(i + i_tra, j + j_tra, k + k_tra, l + l_tra)) > cutoff_threshold) then
                            if (mod(n_cnt, nprocs) == rank) then ! Averaging the size of the subspace 2-integral file per a MPI process
                                write (unit_int2_subspace) i, j, k, l, &
                                    traint2(i + i_tra, j + j_tra, k + k_tra, l + l_tra), 0.0d+00
                            end if
                            n_cnt = n_cnt + 1
                        end if
                    End do
                End do
            End do
        End do
    end subroutine write_traint2_to_disk_fourth_real

    subroutine write_traint2_to_disk_real(ii, ie, ji, je, ki, ke, li, le, traint2, cutoff_threshold, unit_int2_subspace)
        use module_global_variables, only: nprocs, rank
        implicit none
        integer                 :: n_cnt, i, j, k, l
        integer, intent(in)     :: ii, ie, ji, je, ki, ke, li, le, unit_int2_subspace
        real(8)                 :: cutoff_threshold
        real(8), intent(in)  :: traint2(ii:ie, ji:je, ki:ke, li:le)

        n_cnt = 0
        Do l = li, le
            Do k = ki, ke
                Do j = ji, je
                    Do i = ii, ie
                        ! Very small integrals are not written to disk
                        if (ABS(traint2(i, j, k, l)) > cutoff_threshold) then
                            if (mod(n_cnt, nprocs) == rank) then ! Averaging the size of the subspace 2-integral file per a MPI process
                                write (unit_int2_subspace) i, j, k, l, traint2(i, j, k, l), 0.0d+00
                            end if
                            ! if traint2(i,j,k,l)>cutoff_threshold, all MPI process need to count up n_cnt!!!
                            n_cnt = n_cnt + 1
                        end if
                    End do
                End do
            End do
        End do
    end subroutine write_traint2_to_disk_real

    subroutine where_subspace_is(ini, end)
        !=============================================================================================
        ! This subroutine returns the indices of the subspace where the integrals are to be written
        ! to disk.
        ! (e.g.) (input)  ini = global_act_start, end = global_act_end
        !        (output) ini = 1,          end = nact  => the subspace is the active space
        !=============================================================================================
        use module_global_variables, only: ninact, nact, nsec, rank, global_act_end, global_sec_end
        use module_error, only: stop_with_errorcode
        implicit none
        integer, intent(inout) :: ini, end
        if (end == 0) then
            ini = 1
            end = ninact
        elseif (ini <= ninact) then
            ini = 1
            end = ninact
        else if (ini <= global_act_end) then
            ini = 1
            end = nact
        else if (ini <= global_sec_end) then
            ini = 1
            end = nsec
        else
            if (rank == 0) print '(2(a,i0),a)', "Error: The indices of the subspace are out of range. (ini, end) = (", &
                ini, ",", end, ")"
            call stop_with_errorcode(1)
        end if
    end subroutine where_subspace_is
end module module_intra
