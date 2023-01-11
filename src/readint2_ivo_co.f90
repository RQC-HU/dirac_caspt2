! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE readint2_ivo_co(filename) ! 2 electorn integrals created by typart in utchem

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module
    use module_file_manager, only: open_unformatted_file, check_iostat

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    character*50, intent(in) :: filename

    character  :: datex*10, timex*8

    integer    :: mdcint, nkr, nuniq, nmoc, nsec_start, nsec_end, iostat
    integer    :: nz, j0, i0
    integer    :: i, j, k, l, jtr0, itr0
    integer    :: SignIJ, SignKL, itr, jtr, ltr, ktr, inz, totalint, save, count

    complex*16 :: cint2

    integer, allocatable :: indk(:), indl(:), kr(:)
    real*8, allocatable  :: rklr(:), rkli(:)
    logical              :: end_of_file

    nmoc = ninact + nact
    nsec_start = ninact + nact + 1
    nsec_end = ninact + nact + nsec

    Allocate (int2r_f1(nsec_start:nsec_end, nsec_start:nsec_end, nmoc, nmoc))
    Allocate (int2i_f1(nsec_start:nsec_end, nsec_start:nsec_end, nmoc, nmoc))
    Allocate (int2r_f2(nsec_start:nsec_end, nmoc, nmoc, nsec_start:nsec_end))
    Allocate (int2i_f2(nsec_start:nsec_end, nmoc, nmoc, nsec_start:nsec_end))
    Call memplus(KIND(int2r_f1), SIZE(int2r_f1), 1)
    Call memplus(KIND(int2i_f1), SIZE(int2i_f1), 1)
    Call memplus(KIND(int2r_f2), SIZE(int2r_f2), 1)
    Call memplus(KIND(int2i_f2), SIZE(int2i_f2), 1)

    Allocate (indk((nmo/2)**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    Allocate (indl((nmo/2)**2)); Call memplus(KIND(indl), SIZE(indl), 1)
    Allocate (rklr((nmo/2)**2)); Call memplus(KIND(rklr), SIZE(rklr), 1)
    Allocate (rkli((nmo/2)**2)); Call memplus(KIND(rkli), SIZE(rkli), 1)

!Iwamuro modify
    Allocate (kr(-nmo/2:nmo/2)); Call memplus(KIND(kr), SIZE(kr), 1)

    if (rank == 0) print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024

    nuniq = 0
    indk(:) = 0
    indl(:) = 0
    rklr(:) = 0.0d+00
    rkli(:) = 0.0d+00
    int2r_f1 = 0.0d+00
    int2i_f1 = 0.0d+00
    int2r_f2 = 0.0d+00
    int2i_f2 = 0.0d+00

    totalint = 0
    mdcint = 11
    call open_unformatted_file(unit=mdcint, file=trim(filename), status='old', optional_action="read")

    read (mdcint, iostat=iostat) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
    call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=end_of_file)
    if (end_of_file) then
        if (rank == 0) print *, "End of file reached"
        return ! End of file
    end if

    do
        read (mdcint, iostat=iostat) i, j, nz, (indk(inz), indl(inz), inz=1, nz), (rklr(inz), rkli(inz), inz=1, nz)
        call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=end_of_file)
        if (end_of_file) then
            if (rank == 0) print *, "End of file reached"
            exit
        end if

        if (i == 0) exit ! End of file

        totalint = totalint + nz

        itr = i + (-1)**(mod(i, 2) + 1)
        jtr = j + (-1)**(mod(j, 2) + 1)

        i0 = i
        itr0 = itr
        j0 = j
        jtr0 = jtr

        ivo_inz: Do inz = 1, nz

            i = i0
            itr = itr0
            j = j0
            jtr = jtr0

            k = indk(inz)
            ktr = k + (-1)**(mod(k, 2) + 1)
            l = indl(inz)
            ltr = l + (-1)**(mod(l, 2) + 1)

            If (i > nmoc .and. j > nmoc .and. k > nmoc .and. l > nmoc) cycle ivo_inz ! (33|33) is ignored
            If (i == j .and. k > l) cycle ivo_inz

            if (sp(i) == 3 .and. sp(j) == 3 .and. sp(k) < 3 .and. sp(l) == sp(k)) then !(33|11) or (33|22) type

                count = 0
                do
                    if (mod(i, 2) == 0) then
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

                    SignIJ = (-1)**mod(i + j, 2)
                    SignKL = (-1)**mod(k + l, 2)

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
                        cycle
                    else
                        cycle ivo_inz
                    end if
                end do
            elseif (sp(k) == 3 .and. sp(l) == 3 .and. sp(i) < 3 .and. sp(i) == sp(j)) then !(11|33) or (22|33) type

                count = 0
                do
                    if (mod(i, 2) == 0) then
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

                    SignIJ = (-1)**mod(i + j, 2)
                    SignKL = (-1)**mod(k + l, 2)

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
                        cycle
                    else
                        cycle ivo_inz
                    end if
                end do
            elseif (max(sp(i), sp(j)) == 3 .and. max(sp(k), sp(l)) == 3 .and. &
                 &  min(sp(i), sp(j)) == min(sp(k), sp(l))) then                !(31|31) or (32|32) series

                count = 0
                do
                    if (mod(i, 2) == 0) then
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

                    SignIJ = (-1)**mod(i + j, 2)
                    SignKL = (-1)**mod(k + l, 2)

                    if (i > j .and. k > l) then ! (31|31) or (32|32) ==> (31|13) or (32|23)

                        int2r_f2(i, j, ltr, ktr) = signKL*rklr(inz)
                        int2i_f2(i, j, ltr, ktr) = signKL*rkli(inz)

                    elseif (i > j .and. k < l) then ! (31|13) or (32|23) ==> (31|13) or (32|23)

                        int2r_f2(i, j, k, l) = rklr(inz)
                        int2i_f2(i, j, k, l) = rkli(inz)

                    elseif (i < j .and. k < l) then ! (13|13) or (23|23) ==> (31|13) or (32|23)

                        int2r_f2(jtr, itr, k, l) = signIJ*rklr(inz)
                        int2i_f2(jtr, itr, k, l) = signIJ*rkli(inz)

                    elseif (i < j .and. k > l) then ! (13|31) or (23|32) ==> (31|13) or (32|23)

                        int2r_f2(jtr, itr, ltr, ktr) = signIJ*signKL*rklr(inz)
                        int2i_f2(jtr, itr, ltr, ktr) = signIJ*signKL*rkli(inz)

                    end if

                    count = count + 1
                    cint2 = DCMPLX(rklr(inz), rkli(inz))
                    if (count == 1 .or. count == 3) then
                        Call takekr(i, j, k, l, cint2)              ! Consider Kramers pair
                        rklr(inz) = DBLE(cint2)
                        rkli(inz) = DIMAG(cint2)
                        cycle
                    elseif (count == 2) then           ! variables exchange (AA|BB) => (BB|AA)
                        save = i
                        i = k
                        k = save
                        save = j
                        j = l
                        l = save
                        cycle
                    else
                        cycle ivo_inz
                    end if
                end do
            end if

        End do ivo_inz

    end do

    close (mdcint)
    print *, nuniq, totalint, rank
    deallocate (indk); Call memminus(KIND(indk), SIZE(indk), 1)
    deallocate (indl); Call memminus(KIND(indl), SIZE(indl), 1)
    deallocate (rklr); Call memminus(KIND(rklr), SIZE(rklr), 1)
    deallocate (rkli); Call memminus(KIND(rkli), SIZE(rkli), 1)
    deallocate (kr); Call memminus(KIND(kr), SIZE(kr), 1)
#ifdef HAVE_MPI
    ! integer(4) max : 3805902864 == 2**31 - 1 >= ncount
    ! call MPI_Allreduce(MPI_IN_PLACE, int2r_f1(nsec_start, nsec_start, 1, 1), &
    !                    nsec*nsec*nmoc*nmoc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    ! call MPI_Allreduce(MPI_IN_PLACE, int2i_f1(nsec_start, nsec_start, 1, 1), &
    !                    nsec*nsec*nmoc*nmoc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (rank == 0) print *, 'Start MPI_Allreduce int2r_f1, int2i_f1, nmoc = ', nmoc
    do i = 1, nmoc
        if (rank == 0) print *, 'Start MPI_Allreduce int2r_f2, int2i_f2', i
        call MPI_Allreduce(MPI_IN_PLACE, int2r_f1(nsec_start, nsec_start, 1, i), &
                           nsec*nsec*nmoc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, int2i_f1(nsec_start, nsec_start, 1, i), &
                           nsec*nsec*nmoc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    end do
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (rank == 0) print *, "int2r_f1 and int2i_f1 are summed up, nsec = ", nsec
    ! call MPI_Allreduce(MPI_IN_PLACE, int2r_f2(nsec_start, 1, 1, nsec_start), &
    !                    nsec*nmoc*nmoc*nsec, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    ! call MPI_Allreduce(MPI_IN_PLACE, int2i_f2(nsec_start, 1, 1, nsec_start), &
    !                    nsec*nmoc*nmoc*nsec, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    do i = 1, nsec
        if (rank == 0) print *, 'Start MPI_Allreduce int2r_f2, int2i_f2', i
        call MPI_Allreduce(MPI_IN_PLACE, int2r_f2(nsec_start, 1, 1, ninact + nact + i), &
                           nsec*nmoc*nmoc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, int2i_f2(nsec_start, 1, 1, ninact + nact + i), &
                           nsec*nmoc*nmoc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    end do
    if (rank == 0) print *, 'End MPI_Allreduce int2r_f1, int2i_f1, int2r_f2, int2i_f2'
#endif
end subroutine readint2_ivo_co
