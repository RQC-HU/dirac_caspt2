module module_2integrals
    use module_takekr
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
        use module_sort_swap, only: swap
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
        integer    :: i, j, nz
        integer    :: inz
        integer    :: jtr0, itr0
        integer    :: SignIJ, SignKL, itr, jtr, ltr, ktr, totalint, count
        complex*16 :: cint2
        integer, allocatable :: indk(:), indl(:), kr(:)
        real*8, allocatable  :: rklr(:)
        logical :: continue_read, is_end_of_file
        integer :: iostat
        continue_read = .true.
        nmoc = ninact + nact
        nmom = ninact + nact + nsec
        if (rank == 0) print *, "Enter readint2_casci subroutine (realonly)"

        Allocate (int2r_f1(nmoc + 1:nmoc + nsec, nmoc + 1:nmoc + nsec, nmoc, nmoc)); Call memplus(KIND(int2r_f1), SIZE(int2r_f1), 1)
        Allocate (int2r_f2(nmoc + 1:nmoc + nsec, nmoc, nmoc, nmoc + 1:nmoc + nsec)); Call memplus(KIND(int2r_f2), SIZE(int2r_f2), 1)
        Allocate (inttwr(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(inttwr), SIZE(inttwr), 1)
        Allocate (indk(nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
        Allocate (indl(nmo**2)); Call memplus(KIND(indl), SIZE(indl), 1)
        Allocate (rklr(nmo**2)); Call memplus(KIND(rklr), SIZE(rklr), 1)
        !Iwamuro modify
        Allocate (kr(-nmo/2:nmo/2)); Call memplus(KIND(kr), SIZE(kr), 1)

        if (rank == 0) print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024

        ! Initialize variables
        nuniq = 0
        i = 0
        j = 0
        indk(:) = 0
        indl(:) = 0
        rklr(:) = 0.0d+00
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
        do
            read (unit_mdcint, iostat=iostat) i, j, nz, (indk(inz), indl(inz), inz=1, nz), (rklr(inz), inz=1, nz)
            call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            if (i == 0) exit ! Last line of the file

            totalint = totalint + nz

            itr = i + (-1)**(mod(i, 2) + 1)
            jtr = j + (-1)**(mod(j, 2) + 1)

            i0 = i
            itr0 = itr
            j0 = j
            jtr0 = jtr

            loop_inz: Do inz = 1, nz

                i = i0
                itr = itr0
                j = j0
                jtr = jtr0

                k = indk(inz)
                ktr = k + (-1)**(mod(k, 2) + 1)
                l = indl(inz)
                ltr = l + (-1)**(mod(l, 2) + 1)

                If (i > nmoc .and. j > nmoc .and. k > nmoc .and. l > nmoc) cycle loop_inz ! (33|33) is ignored
                If (i == j .and. k > l) cycle loop_inz

                If (i <= nmoc .and. j <= nmoc .and. k <= nmoc .and. l <= nmoc) then
                    SignIJ = (-1)**(mod(i, 2) + mod(j, 2))
                    SignKL = (-1)**(mod(k, 2) + mod(l, 2))
                    nuniq = nuniq + 1
                    !=-> Original integral plus time-reversed partners
                    INTTWR(I, J, K, L) = rklr(inz)
                    INTTWR(JTR, ITR, K, L) = rklr(inz)*SignIJ
                    INTTWR(I, J, LTR, KTR) = rklr(inz)*SignKL
                    INTTWR(JTR, ITR, LTR, KTR) = rklr(inz)*SignIJ*SignKL
                    !=-> Complex conjugate plus time-reversed partners
                    INTTWR(J, I, L, K) = rklr(inz)
                    INTTWR(ITR, JTR, L, K) = rklr(inz)*SignIJ
                    INTTWR(J, I, KTR, LTR) = rklr(inz)*SignKL
                    INTTWR(ITR, JTR, KTR, LTR) = rklr(inz)*SignIJ*SignKL
                    !=-> Particle interchanged plus time-reversed partners
                    INTTWR(K, L, I, J) = rklr(inz)
                    INTTWR(LTR, KTR, I, J) = rklr(inz)*SignKL
                    INTTWR(K, L, JTR, ITR) = rklr(inz)*SignIJ
                    INTTWR(LTR, KTR, JTR, ITR) = rklr(inz)*SignIJ*SignKL
                    !=-> Particle interchanged and complex conjugated plus time-reversed partners
                    INTTWR(L, K, J, I) = rklr(inz)
                    INTTWR(KTR, LTR, J, I) = rklr(inz)*SignKL
                    INTTWR(L, K, ITR, JTR) = rklr(inz)*SignIJ
                    INTTWR(KTR, LTR, ITR, JTR) = rklr(inz)*SignIJ*SignKL

                elseif (space_idx(i) == 3 .and. space_idx(j) == 3 .and. &
                        space_idx(k) < 3 .and. space_idx(l) == space_idx(k)) then !(33|11) or (33|22) type
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

                        int2r_f1(jtr, itr, k, l) = SignIJ*rklr(inz)

                        int2r_f1(i, j, ltr, ktr) = SignKL*rklr(inz)

                        int2r_f1(jtr, itr, ltr, ktr) = SignIJ*SignKL*rklr(inz)

                        count = count + 1
                        cint2 = DCMPLX(rklr(inz), 0.0d+00)

                        if (count == 1) then
                            Call takekr(i, j, k, l, cint2)              ! Consider Kramers pair
                            rklr(inz) = DBLE(cint2)
                            cycle ! Go to the next count loop
                        else
                            cycle loop_inz ! Go to the next inz
                        end if
                    end do
                elseif (space_idx(k) == 3 .and. space_idx(l) == 3 .and. &
                        space_idx(i) < 3 .and. space_idx(i) == space_idx(j)) then !(11|33) or (22|33) type
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

                        int2r_f1(k, l, jtr, itr) = SignIJ*rklr(inz)

                        int2r_f1(ltr, ktr, i, j) = SignKL*rklr(inz)

                        int2r_f1(ltr, ktr, jtr, itr) = SignIJ*SignKL*rklr(inz)

                        count = count + 1
                        cint2 = DCMPLX(rklr(inz), 0.0d+00)
                        if (count == 1) then
                            Call takekr(i, j, k, l, cint2)              ! Consider Kramers pair
                            rklr(inz) = DBLE(cint2)
                            cycle ! Go to the next count loop
                        else
                            cycle loop_inz ! Go to the next inz
                        end if
                    end do

                elseif (max(space_idx(i), space_idx(j)) == 3 .and. max(space_idx(k), space_idx(l)) == 3 .and. &
                      &  min(space_idx(i), space_idx(j)) == min(space_idx(k), space_idx(l))) then                !(31|31) or (32|32) series

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

                        elseif (i > j .and. k < l) then ! (31|13) or (32|23) ==> (31|13) or (32|23)

                            int2r_f2(i, j, k, l) = rklr(inz)

                        elseif (i < j .and. k < l) then ! (13|13) or (23|23) ==> (31|13) or (32|23)

                            int2r_f2(jtr, itr, k, l) = signIJ*rklr(inz)

                        elseif (i < j .and. k > l) then ! (13|31) or (23|32) ==> (31|13) or (32|23)

                            int2r_f2(jtr, itr, ltr, ktr) = signIJ*signKL*rklr(inz)

                        end if

                        count = count + 1
                        cint2 = DCMPLX(rklr(inz), 0.0d+00)
                        if (count == 1 .or. count == 3) then
                            Call takekr(i, j, k, l, cint2)              ! Consider Kramers pair
                            rklr(inz) = DBLE(cint2)
                            cycle ! Go to the next count loop
                        elseif (count == 2) then           ! variables exchange (AA|BB) => (BB|AA)
                            call swap(i, k)
                            call swap(j, l)
                            cycle ! Go to the next count loop
                        else
                            cycle loop_inz ! Go to the next inz
                        end if
                    end do
                else
                end if

            End do loop_inz

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

        if (allocated(indk)) Call memminus(KIND(indk), SIZE(indk), 1); deallocate (indk)
        if (allocated(indl)) Call memminus(KIND(indl), SIZE(indl), 1); deallocate (indl)
        if (allocated(rklr)) Call memminus(KIND(rklr), SIZE(rklr), 1); deallocate (rklr)
        if (allocated(kr)) Call memminus(KIND(kr), SIZE(kr), 1); deallocate (kr)

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
        use module_sort_swap, only: swap
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
        integer    :: i, j, nz
        integer    :: inz
        integer    :: jtr0, itr0
        integer    :: SignIJ, SignKL, itr, jtr, ltr, ktr, totalint, count
        complex*16 :: cint2
        integer, allocatable :: indk(:), indl(:), kr(:)
        real*8, allocatable  :: rklr(:), rkli(:)
        logical :: continue_read, is_end_of_file
        integer :: read_line_len, iostat
        read_line_len = read_line_max ! Set read_line_len as parameter "read_line_max"
        continue_read = .true.
        nmoc = ninact + nact
        nmom = ninact + nact + nsec
        if (rank == 0) print *, "Enter readint2_casci"

        Allocate (int2r_f1(nmoc + 1:nmoc + nsec, nmoc + 1:nmoc + nsec, nmoc, nmoc)); Call memplus(KIND(int2r_f1), SIZE(int2r_f1), 1)
        Allocate (int2i_f1(nmoc + 1:nmoc + nsec, nmoc + 1:nmoc + nsec, nmoc, nmoc)); Call memplus(KIND(int2i_f1), SIZE(int2i_f1), 1)
        Allocate (int2r_f2(nmoc + 1:nmoc + nsec, nmoc, nmoc, nmoc + 1:nmoc + nsec)); Call memplus(KIND(int2r_f2), SIZE(int2r_f2), 1)
        Allocate (int2i_f2(nmoc + 1:nmoc + nsec, nmoc, nmoc, nmoc + 1:nmoc + nsec)); Call memplus(KIND(int2i_f2), SIZE(int2i_f2), 1)
        Allocate (inttwr(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(inttwr), SIZE(inttwr), 1)
        Allocate (inttwi(nmoc, nmoc, nmoc, nmoc)); Call memplus(KIND(inttwi), SIZE(inttwi), 1)
        Allocate (indk(nmo**2)); Call memplus(KIND(indk), SIZE(indk), 1)
        Allocate (indl(nmo**2)); Call memplus(KIND(indl), SIZE(indl), 1)
        Allocate (rklr(nmo**2)); Call memplus(KIND(rklr), SIZE(rklr), 1)
        Allocate (rkli(nmo**2)); Call memplus(KIND(rkli), SIZE(rkli), 1)
        !Iwamuro modify
        Allocate (kr(-nmo/2:nmo/2)); Call memplus(KIND(kr), SIZE(kr), 1)

        if (rank == 0) print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024

        ! Initialize variables
        nuniq = 0
        i = 0
        j = 0
        indk(:) = 0
        indl(:) = 0
        rklr(:) = 0.0d+00
        rkli(:) = 0.0d+00
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
            read (unit_mdcint, iostat=iostat) i, j, nz, &
                (indk(inz), indl(inz), inz=1, nz), (rklr(inz), rkli(inz), inz=1, nz)
            call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
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

            loop_inz: Do inz = 1, nz

                i = i0
                itr = itr0
                j = j0
                jtr = jtr0

                k = indk(inz)
                ktr = k + (-1)**(mod(k, 2) + 1)
                l = indl(inz)
                ltr = l + (-1)**(mod(l, 2) + 1)

                If (i > nmoc .and. j > nmoc .and. k > nmoc .and. l > nmoc) cycle loop_inz ! (33|33) is ignored
                If (i == j .and. k > l) cycle loop_inz

                If (i <= nmoc .and. j <= nmoc .and. k <= nmoc .and. l <= nmoc) then
                    SignIJ = (-1)**(mod(i, 2) + mod(j, 2))
                    SignKL = (-1)**(mod(k, 2) + mod(l, 2))
                    nuniq = nuniq + 1
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

                elseif (space_idx(i) == 3 .and. space_idx(j) == 3 .and. &
                        space_idx(k) < 3 .and. space_idx(l) == space_idx(k)) then !(33|11) or (33|22) type
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
                            cycle ! Go to the next count loop
                        else
                            cycle loop_inz ! Go to the next inz
                        end if
                    end do
                elseif (space_idx(k) == 3 .and. space_idx(l) == 3 .and. &
                        space_idx(i) < 3 .and. space_idx(i) == space_idx(j)) then !(11|33) or (22|33) type
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
                            cycle ! Go to the next count loop
                        else
                            cycle loop_inz ! Go to the next inz
                        end if
                    end do

                elseif (max(space_idx(i), space_idx(j)) == 3 .and. max(space_idx(k), space_idx(l)) == 3 .and. &
                      &  min(space_idx(i), space_idx(j)) == min(space_idx(k), space_idx(l))) then                !(31|31) or (32|32) series

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
                            cycle ! Go to the next count loop
                        elseif (count == 2) then           ! variables exchange (AA|BB) => (BB|AA)
                            call swap(i, k)
                            call swap(j, l)
                            cycle ! Go to the next count loop
                        else
                            cycle loop_inz ! Go to the next inz
                        end if
                    end do
                else
                end if

            End do loop_inz
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

        if (allocated(indk)) Call memminus(KIND(indk), SIZE(indk), 1); deallocate (indk)
        if (allocated(indl)) Call memminus(KIND(indl), SIZE(indl), 1); deallocate (indl)
        if (allocated(rklr)) Call memminus(KIND(rklr), SIZE(rklr), 1); deallocate (rklr)
        if (allocated(rkli)) Call memminus(KIND(rkli), SIZE(rkli), 1); deallocate (rkli)
        if (allocated(kr)) Call memminus(KIND(kr), SIZE(kr), 1); deallocate (kr)
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
