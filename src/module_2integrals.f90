module module_2integrals
    implicit none
    private
    public readint2_casci, get_2_elec_integral

contains

    complex*16 function get_2_elec_integral(integral_real_name, idx1, idx2, idx3, idx4)
        use four_caspt2_module, only: rank, int2r_f1, int2i_f1, int2r_f2, int2i_f2, inttwr, inttwi
        use module_error, only: stop_with_errorcode
        use module_realonly, only: realonly
        integer, intent(in) :: idx1, idx2, idx3, idx4
        character(8), parameter :: integarl_names(3) = ["int2r_f1", "int2r_f2", "inttwr  "]
        character(len=*), intent(in) :: integral_real_name ! Name of the integral (int2r_f1, int2r_f2, inttwr)
        real(8), pointer :: real_integrals(:, :, :, :), imag_integrals(:, :, :, :)

        ! Get the pointer to the real and imaginary part of the integrals
        select case (trim(integral_real_name))
        case (trim(integarl_names(1))) ! int2r_f1
            real_integrals => int2r_f1
            imag_integrals => int2i_f1
        case (trim(integarl_names(2))) ! int2r_f2
            real_integrals => int2r_f2
            imag_integrals => int2i_f2
        case (trim(integarl_names(3))) ! inttwr
            real_integrals => inttwr
            imag_integrals => inttwi
        case default
            if (rank == 0) then
                print *, "Error: Integral name is not correct"
                print *, "Integral name must be one of the following:"
                print *, integarl_names
            end if
            call stop_with_errorcode(1)
        end select

        if (realonly%is_realonly()) then
            get_2_elec_integral = DCMPLX(real_integrals(idx1, idx2, idx3, idx4), 0.0d+00)
        else
            get_2_elec_integral = DCMPLX(real_integrals(idx1, idx2, idx3, idx4), &
                                         imag_integrals(idx1, idx2, idx3, idx4))
        end if
    end function get_2_elec_integral

    subroutine readint2_casci(filename, nuniq)
        use module_realonly, only: realonly
        implicit none
        character*50, intent(in) :: filename
        integer, intent(out) :: nuniq
        real(8) :: dummy_real
        complex*16 :: dummy_complex
        dummy_real = 0.0d+00
        dummy_complex = (0.0d+00, 0.0d+00)
        if (realonly%is_realonly()) then ! Realonly
            call readint2_casci_realonly(filename, nuniq, dummy_real)
        else ! Complex
            call readint2_casci_complex(filename, nuniq, dummy_complex)
        end if

    end subroutine readint2_casci

    SUBROUTINE readint2_casci_realonly(filename, nuniq, dummy_real)  ! 2 electorn integrals created by typart in utchem (Realonly)

        use four_caspt2_module
        use module_file_manager
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE
        character*50, intent(in) :: filename
        real(8), intent(in) :: dummy_real

        character  :: datex*10, timex*8
        integer    :: unit_mdcint, nkr, nuniq, nmom, nmoc
        integer    :: j0, i0
        integer    :: k, l
        integer, allocatable    :: i(:), j(:), nz(:)
        integer    :: inz
        integer    :: jtr0, itr0
        integer    :: SignIJ, SignKL, itr, jtr, ltr, ktr, totalint, save, count
        complex*16 :: cint2
        integer, allocatable :: indk(:, :), indl(:, :), kr(:)
        real*8, allocatable  :: rklr(:, :)
        logical :: continue_read, is_end_of_file
        integer :: idx, read_line_len, iostat
        read_line_len = read_line_max ! Set read_line_len as parameter "read_line_max"
        continue_read = .true.
        nmoc = ninact + nact
        nmom = ninact + nact + nsec
        if (rank == 0) print *, "Enter readint2_casci subroutine (realonly)"

        Allocate (i(read_line_max)); call memplus(kind(i), size(i), 1)
        Allocate (j(read_line_max)); call memplus(kind(j), size(j), 1)
        Allocate (nz(read_line_max)); call memplus(kind(nz), size(nz), 1)
        Allocate (int2r_f1(nmoc + 1:nmoc + nsec, nmoc + 1:nmoc + nsec, nmoc, nmoc)); Call memplus(KIND(int2r_f1), SIZE(int2r_f1), 1)
        Allocate (int2r_f2(nmoc + 1:nmoc + nsec, nmoc, nmoc, nmoc + 1:nmoc + nsec)); Call memplus(KIND(int2r_f2), SIZE(int2r_f2), 1)
        Allocate (inttwr(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(inttwr), SIZE(inttwr), 1)
        Allocate (indk(read_line_max, nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
        Allocate (indl(read_line_max, nmo**2)); Call memplus(KIND(indl), SIZE(indl), 1)
        Allocate (rklr(read_line_max, nmo**2)); Call memplus(KIND(rklr), SIZE(rklr), 1)
        !Iwamuro modify
        Allocate (kr(-nmo/2:nmo/2)); Call memplus(KIND(kr), SIZE(kr), 1)

        if (rank == 0) print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024

        ! Initialize variables
        nuniq = 0
        i(:) = 0
        j(:) = 0
        indk(:, :) = 0
        indl(:, :) = 0
        rklr(:, :) = 0.0d+00
        inttwr = 0.0d+00
        int2r_f1 = 0.0d+00
        int2r_f2 = 0.0d+00

        totalint = 0
        call open_unformatted_file(unit=unit_mdcint, file=trim(filename), status='old', optional_action='read')

        read (unit_mdcint, iostat=iostat) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
        call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
        if (is_end_of_file) then
            continue_read = .false.
        end if

        if (rank == 0) then
            print *, datex, timex
            print *, 'readint2', 'nkr', nkr, 'kr(+),kr(-)', (kr(i0), kr(-1*i0), i0=1, nkr)
        end if
        do while (continue_read)
            do idx = 1, read_line_max
                read (unit_mdcint, iostat=iostat) i(idx), j(idx), nz(idx), &
                    (indk(idx, inz), indl(idx, inz), rklr(idx, inz), inz=1, nz(idx))
                call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
                if (is_end_of_file) then
                    continue_read = .false.
                    exit
                end if
            end do

            ! The length of read line is equal to min(read_line_max, idx)
            read_line_len = min(read_line_max, idx)

!$OMP parallel do private(idx,itr,jtr,i0,itr0,j0,jtr0,inz,k,ktr,l,ltr,SIGNIJ,SIGNKL,cint2,save,count) &
!$OMP & reduction(+:totalint,nuniq)
            do idx = 1, read_line_len
                if (i(idx) == 0) cycle ! Go to next idx

                totalint = totalint + nz(idx)

                itr = i(idx) + (-1)**(mod(i(idx), 2) + 1)
                jtr = j(idx) + (-1)**(mod(j(idx), 2) + 1)

                i0 = i(idx)
                itr0 = itr
                j0 = j(idx)
                jtr0 = jtr

                loop_inz: Do inz = 1, nz(idx)

                    i(idx) = i0
                    itr = itr0
                    j(idx) = j0
                    jtr = jtr0

                    k = indk(idx, inz)
                    ktr = k + (-1)**(mod(k, 2) + 1)
                    l = indl(idx, inz)
                    ltr = l + (-1)**(mod(l, 2) + 1)

                    If (i(idx) > nmoc .and. j(idx) > nmoc .and. k > nmoc .and. l > nmoc) cycle loop_inz ! (33|33) is ignored
                    If (i(idx) == j(idx) .and. k > l) cycle loop_inz

                    If (i(idx) <= nmoc .and. j(idx) <= nmoc .and. k <= nmoc .and. l <= nmoc) then
                        SignIJ = (-1)**(mod(i(idx), 2) + mod(j(idx), 2))
                        SignKL = (-1)**(mod(k, 2) + mod(l, 2))
                        nuniq = nuniq + 1
                        !=-> Original integral plus time-reversed partners
                        INTTWR(I(idx), J(idx), K, L) = rklr(idx, inz)
                        INTTWR(JTR, ITR, K, L) = rklr(idx, inz)*SignIJ
                        INTTWR(I(idx), J(idx), LTR, KTR) = rklr(idx, inz)*SignKL
                        INTTWR(JTR, ITR, LTR, KTR) = rklr(idx, inz)*SignIJ*SignKL
                        !=-> Complex conjugate plus time-reversed partners
                        INTTWR(J(idx), I(idx), L, K) = rklr(idx, inz)
                        INTTWR(ITR, JTR, L, K) = rklr(idx, inz)*SignIJ
                        INTTWR(J(idx), I(idx), KTR, LTR) = rklr(idx, inz)*SignKL
                        INTTWR(ITR, JTR, KTR, LTR) = rklr(idx, inz)*SignIJ*SignKL
                        !=-> Particle interchanged plus time-reversed partners
                        INTTWR(K, L, I(idx), J(idx)) = rklr(idx, inz)
                        INTTWR(LTR, KTR, I(idx), J(idx)) = rklr(idx, inz)*SignKL
                        INTTWR(K, L, JTR, ITR) = rklr(idx, inz)*SignIJ
                        INTTWR(LTR, KTR, JTR, ITR) = rklr(idx, inz)*SignIJ*SignKL
                        !=-> Particle interchanged and complex conjugated plus time-reversed partners
                        INTTWR(L, K, J(idx), I(idx)) = rklr(idx, inz)
                        INTTWR(KTR, LTR, J(idx), I(idx)) = rklr(idx, inz)*SignKL
                        INTTWR(L, K, ITR, JTR) = rklr(idx, inz)*SignIJ
                        INTTWR(KTR, LTR, ITR, JTR) = rklr(idx, inz)*SignIJ*SignKL

                    elseif (space_idx(i(idx)) == 3 .and. space_idx(j(idx)) == 3 .and. &
                            space_idx(k) < 3 .and. space_idx(l) == space_idx(k)) then !(33|11) or (33|22) type
                        count = 0
                        do
                            if (mod(i(idx), 2) == 0) then
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

                            SignIJ = (-1)**mod(i(idx) + j(idx), 2)
                            SignKL = (-1)**mod(k + l, 2)

                            int2r_f1(i(idx), j(idx), k, l) = rklr(idx, inz)

                            int2r_f1(jtr, itr, k, l) = SignIJ*rklr(idx, inz)

                            int2r_f1(i(idx), j(idx), ltr, ktr) = SignKL*rklr(idx, inz)

                            int2r_f1(jtr, itr, ltr, ktr) = SignIJ*SignKL*rklr(idx, inz)

                            count = count + 1
                            cint2 = DCMPLX(rklr(idx, inz), 0.0d+00)

                            if (count == 1) then
                                Call takekr(i(idx), j(idx), k, l, cint2)              ! Consider Kramers pair
                                rklr(idx, inz) = DBLE(cint2)
                                cycle ! Go to the next count loop
                            else
                                cycle loop_inz ! Go to the next inz
                            end if
                        end do
                    elseif (space_idx(k) == 3 .and. space_idx(l) == 3 .and. &
                            space_idx(i(idx)) < 3 .and. space_idx(i(idx)) == space_idx(j(idx))) then !(11|33) or (22|33) type
                        count = 0
                        do
                            if (mod(i(idx), 2) == 0) then
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

                            SignIJ = (-1)**mod(i(idx) + j(idx), 2)
                            SignKL = (-1)**mod(k + l, 2)

                            int2r_f1(k, l, i(idx), j(idx)) = rklr(idx, inz)

                            int2r_f1(k, l, jtr, itr) = SignIJ*rklr(idx, inz)

                            int2r_f1(ltr, ktr, i(idx), j(idx)) = SignKL*rklr(idx, inz)

                            int2r_f1(ltr, ktr, jtr, itr) = SignIJ*SignKL*rklr(idx, inz)

                            count = count + 1
                            cint2 = DCMPLX(rklr(idx, inz), 0.0d+00)
                            if (count == 1) then
                                Call takekr(i(idx), j(idx), k, l, cint2)              ! Consider Kramers pair
                                rklr(idx, inz) = DBLE(cint2)
                                cycle ! Go to the next count loop
                            else
                                cycle loop_inz ! Go to the next inz
                            end if
                        end do

                    elseif (max(space_idx(i(idx)), space_idx(j(idx))) == 3 .and. max(space_idx(k), space_idx(l)) == 3 .and. &
                          &  min(space_idx(i(idx)), space_idx(j(idx))) == min(space_idx(k), space_idx(l))) then                !(31|31) or (32|32) series

                        count = 0

                        do
                            if (mod(i(idx), 2) == 0) then
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

                            SignIJ = (-1)**mod(i(idx) + j(idx), 2)
                            SignKL = (-1)**mod(k + l, 2)

                            if (i(idx) > j(idx) .and. k > l) then ! (31|31) or (32|32) ==> (31|13) or (32|23)

                                int2r_f2(i(idx), j(idx), ltr, ktr) = signKL*rklr(idx, inz)

                            elseif (i(idx) > j(idx) .and. k < l) then ! (31|13) or (32|23) ==> (31|13) or (32|23)

                                int2r_f2(i(idx), j(idx), k, l) = rklr(idx, inz)

                            elseif (i(idx) < j(idx) .and. k < l) then ! (13|13) or (23|23) ==> (31|13) or (32|23)

                                int2r_f2(jtr, itr, k, l) = signIJ*rklr(idx, inz)

                            elseif (i(idx) < j(idx) .and. k > l) then ! (13|31) or (23|32) ==> (31|13) or (32|23)

                                int2r_f2(jtr, itr, ltr, ktr) = signIJ*signKL*rklr(idx, inz)

                            end if

                            count = count + 1
                            cint2 = DCMPLX(rklr(idx, inz), 0.0d+00)
                            if (count == 1 .or. count == 3) then
                                Call takekr(i(idx), j(idx), k, l, cint2)              ! Consider Kramers pair
                                rklr(idx, inz) = DBLE(cint2)
                                cycle ! Go to the next count loop
                            elseif (count == 2) then           ! variables exchange (AA|BB) => (BB|AA)
                                save = i(idx)
                                i(idx) = k
                                k = save
                                save = j(idx)
                                j(idx) = l
                                l = save
                                cycle ! Go to the next count loop
                            else
                                cycle loop_inz ! Go to the next inz
                            end if
                        end do
                    else
                    end if

                End do loop_inz
            end do
!$OMP end parallel do

            ! Initialize i and continue to read
            i(:) = 0

        end do
        if (rank == 0) print *, 'end Read mdcint normal'

        close (unit_mdcint)
#ifdef HAVE_MPI
        call reduce_wrapper(mat=nuniq, root_rank=0)
        call reduce_wrapper(mat=totalint, root_rank=0)
#endif
        if (rank == 0) then
            print *, nuniq, totalint
        end if

        if (allocated(indk)) deallocate (indk); Call memminus(KIND(indk), SIZE(indk), 1)
        if (allocated(indl)) deallocate (indl); Call memminus(KIND(indl), SIZE(indl), 1)
        if (allocated(rklr)) deallocate (rklr); Call memminus(KIND(rklr), SIZE(rklr), 1)
        if (allocated(kr)) deallocate (kr); Call memminus(KIND(kr), SIZE(kr), 1)
        if (allocated(i)) deallocate (i); Call memminus(KIND(i), SIZE(i), 1)
        if (allocated(j)) deallocate (j); Call memminus(KIND(j), SIZE(j), 1)
        if (allocated(nz)) deallocate (nz); call memminus(kind(nz), size(nz), 1)
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=inttwr)
        call allreduce_wrapper(mat=int2r_f1)
        call allreduce_wrapper(mat=int2r_f2)
        if (rank == 0) print *, 'End MPI_Allreduce inttwr, inttwi, int2r_f1, int2i_f1, int2r_f2, int2i_f2'
#endif
    end subroutine readint2_casci_realonly

    SUBROUTINE readint2_casci_complex(filename, nuniq, dummy_complex)  ! 2 electorn integrals created by typart in utchem (Complex)

        use four_caspt2_module
        use module_file_manager
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE
        character*50, intent(in) :: filename
        complex*16, intent(in) :: dummy_complex

        character  :: datex*10, timex*8
        integer    :: unit_mdcint, nkr, nuniq, nmom, nmoc
        integer    :: j0, i0
        integer    :: k, l
        integer, allocatable    :: i(:), j(:), nz(:)
        integer    :: inz
        integer    :: jtr0, itr0
        integer    :: SignIJ, SignKL, itr, jtr, ltr, ktr, totalint, save, count
        complex*16 :: cint2
        integer, allocatable :: indk(:, :), indl(:, :), kr(:)
        real*8, allocatable  :: rklr(:, :), rkli(:, :)
        logical :: continue_read, is_end_of_file
        integer :: idx, read_line_len, iostat
        read_line_len = read_line_max ! Set read_line_len as parameter "read_line_max"
        continue_read = .true.
        nmoc = ninact + nact
        nmom = ninact + nact + nsec
        if (rank == 0) print *, "Enter readint2_casci"

        Allocate (i(read_line_max)); call memplus(kind(i), size(i), 1)
        Allocate (j(read_line_max)); call memplus(kind(j), size(j), 1)
        Allocate (nz(read_line_max)); call memplus(kind(nz), size(nz), 1)
        Allocate (int2r_f1(nmoc + 1:nmoc + nsec, nmoc + 1:nmoc + nsec, nmoc, nmoc)); Call memplus(KIND(int2r_f1), SIZE(int2r_f1), 1)
        Allocate (int2i_f1(nmoc + 1:nmoc + nsec, nmoc + 1:nmoc + nsec, nmoc, nmoc)); Call memplus(KIND(int2i_f1), SIZE(int2i_f1), 1)
        Allocate (int2r_f2(nmoc + 1:nmoc + nsec, nmoc, nmoc, nmoc + 1:nmoc + nsec)); Call memplus(KIND(int2r_f2), SIZE(int2r_f2), 1)
        Allocate (int2i_f2(nmoc + 1:nmoc + nsec, nmoc, nmoc, nmoc + 1:nmoc + nsec)); Call memplus(KIND(int2i_f2), SIZE(int2i_f2), 1)
        Allocate (inttwr(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(inttwr), SIZE(inttwr), 1)
        Allocate (inttwi(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(inttwi), SIZE(inttwi), 1)
        Allocate (indk(read_line_max, nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
        Allocate (indl(read_line_max, nmo**2)); Call memplus(KIND(indl), SIZE(indl), 1)
        Allocate (rklr(read_line_max, nmo**2)); Call memplus(KIND(rklr), SIZE(rklr), 1)
        Allocate (rkli(read_line_max, nmo**2)); Call memplus(KIND(rkli), SIZE(rkli), 1)
        !Iwamuro modify
        Allocate (kr(-nmo/2:nmo/2)); Call memplus(KIND(kr), SIZE(kr), 1)

        if (rank == 0) print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024

        ! Initialize variables
        nuniq = 0
        i(:) = 0
        j(:) = 0
        indk(:, :) = 0
        indl(:, :) = 0
        rklr(:, :) = 0.0d+00
        rkli(:, :) = 0.0d+00
        inttwr = 0.0d+00
        inttwi = 0.0d+00
        int2r_f1 = 0.0d+00
        int2i_f1 = 0.0d+00
        int2r_f2 = 0.0d+00
        int2i_f2 = 0.0d+00

        totalint = 0
        call open_unformatted_file(unit=unit_mdcint, file=trim(filename), status='old', optional_action='read')

        read (unit_mdcint, iostat=iostat) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)

        call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
        if (is_end_of_file) then
            continue_read = .false.
        end if

        if (rank == 0) then
            print *, datex, timex
            print *, 'readint2', 'nkr', nkr, 'kr(+),kr(-)', (kr(i0), kr(-1*i0), i0=1, nkr)
        end if
        do while (continue_read)
            do idx = 1, read_line_max
                read (unit_mdcint, iostat=iostat) i(idx), j(idx), nz(idx), &
                    (indk(idx, inz), indl(idx, inz), rklr(idx, inz), rkli(idx, inz), inz=1, nz(idx))
                call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
                if (is_end_of_file) then
                    continue_read = .false.
                    exit
                end if
            end do

            ! The length of read line is equal to min(read_line_max, idx)
            read_line_len = min(read_line_max, idx)

!$OMP parallel do private(idx,itr,jtr,i0,itr0,j0,jtr0,inz,k,ktr,l,ltr,SIGNIJ,SIGNKL,cint2,save,count) &
!$OMP & reduction(+:totalint,nuniq)
            do idx = 1, read_line_len
                if (i(idx) == 0) cycle ! Go to next idx

                totalint = totalint + nz(idx)

                itr = i(idx) + (-1)**(mod(i(idx), 2) + 1)
                jtr = j(idx) + (-1)**(mod(j(idx), 2) + 1)

                i0 = i(idx)
                itr0 = itr
                j0 = j(idx)
                jtr0 = jtr

                loop_inz: Do inz = 1, nz(idx)

                    i(idx) = i0
                    itr = itr0
                    j(idx) = j0
                    jtr = jtr0

                    k = indk(idx, inz)
                    ktr = k + (-1)**(mod(k, 2) + 1)
                    l = indl(idx, inz)
                    ltr = l + (-1)**(mod(l, 2) + 1)

                    If (i(idx) > nmoc .and. j(idx) > nmoc .and. k > nmoc .and. l > nmoc) cycle loop_inz ! (33|33) is ignored
                    If (i(idx) == j(idx) .and. k > l) cycle loop_inz

                    If (i(idx) <= nmoc .and. j(idx) <= nmoc .and. k <= nmoc .and. l <= nmoc) then
                        SignIJ = (-1)**(mod(i(idx), 2) + mod(j(idx), 2))
                        SignKL = (-1)**(mod(k, 2) + mod(l, 2))
                        nuniq = nuniq + 1
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

                    elseif (space_idx(i(idx)) == 3 .and. space_idx(j(idx)) == 3 .and. &
                            space_idx(k) < 3 .and. space_idx(l) == space_idx(k)) then !(33|11) or (33|22) type
                        count = 0
                        do
                            if (mod(i(idx), 2) == 0) then
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

                            SignIJ = (-1)**mod(i(idx) + j(idx), 2)
                            SignKL = (-1)**mod(k + l, 2)

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
                                cycle ! Go to the next count loop
                            else
                                cycle loop_inz ! Go to the next inz
                            end if
                        end do
                    elseif (space_idx(k) == 3 .and. space_idx(l) == 3 .and. &
                            space_idx(i(idx)) < 3 .and. space_idx(i(idx)) == space_idx(j(idx))) then !(11|33) or (22|33) type
                        count = 0
                        do
                            if (mod(i(idx), 2) == 0) then
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

                            SignIJ = (-1)**mod(i(idx) + j(idx), 2)
                            SignKL = (-1)**mod(k + l, 2)

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
                                cycle ! Go to the next count loop
                            else
                                cycle loop_inz ! Go to the next inz
                            end if
                        end do

                    elseif (max(space_idx(i(idx)), space_idx(j(idx))) == 3 .and. max(space_idx(k), space_idx(l)) == 3 .and. &
                          &  min(space_idx(i(idx)), space_idx(j(idx))) == min(space_idx(k), space_idx(l))) then                !(31|31) or (32|32) series

                        count = 0

                        do
                            if (mod(i(idx), 2) == 0) then
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

                            SignIJ = (-1)**mod(i(idx) + j(idx), 2)
                            SignKL = (-1)**mod(k + l, 2)

                            if (i(idx) > j(idx) .and. k > l) then ! (31|31) or (32|32) ==> (31|13) or (32|23)

                                int2r_f2(i(idx), j(idx), ltr, ktr) = signKL*rklr(idx, inz)
                                int2i_f2(i(idx), j(idx), ltr, ktr) = signKL*rkli(idx, inz)

                            elseif (i(idx) > j(idx) .and. k < l) then ! (31|13) or (32|23) ==> (31|13) or (32|23)

                                int2r_f2(i(idx), j(idx), k, l) = rklr(idx, inz)
                                int2i_f2(i(idx), j(idx), k, l) = rkli(idx, inz)

                            elseif (i(idx) < j(idx) .and. k < l) then ! (13|13) or (23|23) ==> (31|13) or (32|23)

                                int2r_f2(jtr, itr, k, l) = signIJ*rklr(idx, inz)
                                int2i_f2(jtr, itr, k, l) = signIJ*rkli(idx, inz)

                            elseif (i(idx) < j(idx) .and. k > l) then ! (13|31) or (23|32) ==> (31|13) or (32|23)

                                int2r_f2(jtr, itr, ltr, ktr) = signIJ*signKL*rklr(idx, inz)
                                int2i_f2(jtr, itr, ltr, ktr) = signIJ*signKL*rkli(idx, inz)

                            end if

                            count = count + 1
                            cint2 = DCMPLX(rklr(idx, inz), rkli(idx, inz))
                            if (count == 1 .or. count == 3) then
                                Call takekr(i(idx), j(idx), k, l, cint2)              ! Consider Kramers pair
                                rklr(idx, inz) = DBLE(cint2)
                                rkli(idx, inz) = DIMAG(cint2)
                                cycle ! Go to the next count loop
                            elseif (count == 2) then           ! variables exchange (AA|BB) => (BB|AA)
                                save = i(idx)
                                i(idx) = k
                                k = save
                                save = j(idx)
                                j(idx) = l
                                l = save
                                cycle ! Go to the next count loop
                            else
                                cycle loop_inz ! Go to the next inz
                            end if
                        end do
                    else
                    end if

                End do loop_inz
            end do
!$OMP end parallel do

            ! Initialize i and continue to read
            i(:) = 0

        end do
        if (rank == 0) print *, 'end Read mdcint normal'

        close (unit_mdcint)
#ifdef HAVE_MPI
        call reduce_wrapper(mat=nuniq, root_rank=0)
        call reduce_wrapper(mat=totalint, root_rank=0)
#endif
        if (rank == 0) then
            print *, nuniq, totalint
        end if

        if (allocated(indk)) deallocate (indk); Call memminus(KIND(indk), SIZE(indk), 1)
        if (allocated(indl)) deallocate (indl); Call memminus(KIND(indl), SIZE(indl), 1)
        if (allocated(rklr)) deallocate (rklr); Call memminus(KIND(rklr), SIZE(rklr), 1)
        if (allocated(rkli)) deallocate (rkli); Call memminus(KIND(rkli), SIZE(rkli), 1)
        if (allocated(kr)) deallocate (kr); Call memminus(KIND(kr), SIZE(kr), 1)
        if (allocated(i)) deallocate (i); Call memminus(KIND(i), SIZE(i), 1)
        if (allocated(j)) deallocate (j); Call memminus(KIND(j), SIZE(j), 1)
        if (allocated(nz)) deallocate (nz); call memminus(kind(nz), size(nz), 1)
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=inttwr)
        call allreduce_wrapper(mat=inttwi)
        call allreduce_wrapper(mat=int2r_f1)
        call allreduce_wrapper(mat=int2i_f1)
        call allreduce_wrapper(mat=int2r_f2)
        call allreduce_wrapper(mat=int2i_f2)
        if (rank == 0) print *, 'End MPI_Allreduce inttwr, inttwi, int2r_f1, int2i_f1, int2r_f2, int2i_f2'
#endif
    end subroutine readint2_casci_complex

end module module_2integrals