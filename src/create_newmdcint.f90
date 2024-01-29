!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Subroutine create_newmdcint ! 2 Electorn Integrals In Mdcint

! This subroutine creates UTChem type MDCINTNEW files from DIRAC type MDCINT files.
! MDCINTNEW files are used in this program as 2-electron integrals.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use module_error, only: stop_with_errorcode
    use module_file_manager
    Use module_global_variables
    use module_realonly, only: realonly
    use module_sort_swap, only: swap
    Implicit None
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    Character  :: datex*10, timex*8
    integer :: i0, inz, nnz, n, dirac23_int_amp_factor, idx_amp_factor
    integer :: ikr, jkr
    integer :: cur_i, cur_j, cur_k, cur_l
    integer :: ii, jj, kk, ll
    integer :: iikr, jjkr, kkkr, llkr, iii, jjj, kkk, lll
    integer, allocatable :: indk(:), indl(:), kr(:)
    double precision, allocatable  :: rklr(:), rkli(:)
    double precision :: cur_int2_real, cur_int2_imag
    real(8) :: cutoff
    integer :: iiit, jjjt, kkkt, lllt
    integer :: nkr, nz, file_idx, iostat
    integer :: unit_mdcint, unit_mdcintnew
    logical :: is_file_exist, is_end_of_file
    character(*), parameter :: fmt_filename = "dcas_MDCINT"
    integer :: unit_fmt_filename

    if (rank == 0) print *, 'Start create_newmdcint'
    Allocate (kr(-nmo/2:nmo/2))
    kr = 0
    ! Get datex, timex, nkr, and kr from MDCINT becasuse there is no kr information in the MDCINXXX files.
    if (rank == 0) then
        call open_unformatted_file(unit=unit_mdcint, file="MDCINT", status="old", optional_action="read")
        read (unit_mdcint) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
        close (unit_mdcint)
    end if
    Allocate (indk(nmo**2))
    Allocate (indl(nmo**2))
    Allocate (rklr(nmo**2))
    if (.not. realonly%is_realonly()) Allocate (rkli(nmo**2))

#ifdef HAVE_MPI
    ! Broadcast kr and other data that are not included in the MDCINXXX files
    call MPI_Bcast(datex, sizeof(datex), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
#ifdef DEBUG
    if (rank == 0) then
        print *, "datex broadcast"
        print *, "if ierr == 0, datex broadcast successed. ierr=", ierr
    end if
#endif
    if (ierr /= 0) call stop_with_errorcode(ierr)
    call MPI_Bcast(timex, sizeof(timex), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
#ifdef DEBUG
    if (rank == 0) then
        print *, "timex broadcast"
        print *, "if ierr == 0, timex broadcast successed. ierr=", ierr
    end if
#endif
    if (ierr /= 0) call stop_with_errorcode(ierr)
    call MPI_Bcast(nkr, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
#ifdef DEBUG
    if (rank == 0) then
        print *, "nkr broadcast"
        print *, "if ierr == 0, nkr broadcast successed. ierr=", ierr
    end if
#endif
    if (ierr /= 0) call stop_with_errorcode(ierr)
    call MPI_Bcast(kr(-nmo/2), nmo + 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
#ifdef DEBUG
    if (rank == 0) then
        print *, "kr broadcast"
        print *, "if ierr == 0, kr broadcast successed. ierr=", ierr
    end if
#endif
    if (ierr /= 0) call stop_with_errorcode(ierr)
    call MPI_Bcast(indmo_dirac_to_cas(1), nmo, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
#ifdef DEBUG
    if (rank == 0) then
        print *, "datex broadcast"
        print *, "if ierr == 0, datex broadcast successed. ierr=", ierr
    end if
#endif
    if (ierr /= 0) call stop_with_errorcode(ierr)
#endif

    cutoff = 0.25D-12
    nnz = 1

    ! Initialization
    file_idx = 0

    call open_formatted_file(unit=unit_fmt_filename, file=fmt_filename, status="replace", optional_action="write")
    write (unit_fmt_filename, *) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
    ! First, All process write header information to MDCINTNEWrank
    call get_mdcint_filename(file_idx)
    call open_unformatted_file(unit=unit_mdcintnew, file=mdcintNew, status="replace", optional_action="write")
    write (unit_mdcintnew) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)

    is_file_exist = .true.
    do while (is_file_exist) ! Continue reading 2-electron integrals until mdcint_filename doesn't exist.

        inquire (file=mdcint_filename, exist=is_file_exist) ! mdcint_filename exists?
        if (.not. is_file_exist) exit ! Exit do while loop if mdcint_filename doesn't exist.
        call open_unformatted_file(unit=unit_mdcint, file=mdcint_filename, status="old", optional_action="read")
        rewind (unit_mdcint)
        read (unit_mdcint)

        ! Continue to read 2-electron integrals until mdcint_filename reaches the end of file.
        mdcint_file_read: do
            if (realonly%is_realonly()) then
                read (unit_mdcint, iostat=iostat) ikr, jkr, nz, (indk(inz), indl(inz), inz=1, nz), (rklr(inz), inz=1, nz)
            else
                read (unit_mdcint, iostat=iostat) ikr, jkr, nz, (indk(inz), indl(inz), inz=1, nz), (rklr(inz), rkli(inz), inz=1, nz)
            end if

            call check_iostat(iostat=iostat, file=mdcint_filename, end_of_file_reached=is_end_of_file)
            if (is_end_of_file .or. ikr == 0) then
                exit mdcint_file_read ! End of file
            end if

!------------------------------!
!   Create new ikr for UTChem  !
!------------------------------!

!       new ikr = iikr
!           jkr = jjkr
!           kkr = kkkr
!           lkr = llkr

            do n = 1, 2
                if (n == 2) then ! Amplify the integral
                    ! https://gitlab.com/dirac/dirac/-/blob/0a337b78c4a8ebfb392ecf23bdef5ff471016210/utils/dirac_mointegral_export.F90#L223-227
                    ! Based on the above URL, the integral is amplified, but for some reason,
                    ! if the signs of rklr and rkli are not reversed, the result is wrong.
                    ! The cause is currently under investigation.
                    ikr = -ikr
                    jkr = -jkr
                    indk(:) = -indk(:)
                    indl(:) = -indl(:)
                    rklr(:) = -rklr(:)*sign(1, ikr*jkr*indk(:)*indl(:))
                    if (.not. realonly%is_realonly()) rkli(:) = rkli(:)*sign(1, ikr*jkr*indk(:)*indl(:))
                end if
                !   !$OMP parallel do private(iii,jjj,kkk,lll,iikr,jjkr,kkkr,llkr,iiit,jjjt,kkkt,lllt,ii,jj,kk,ll)
                dirac23_int_amp_factor = 1
                if (dirac_version >= 23 .and. abs(ikr) /= abs(jkr)) dirac23_int_amp_factor = 2  ! https://github.com/RQC-HU/dirac_caspt2/issues/139#issuecomment-1906140334
                Do inz = 1, nz

                    do idx_amp_factor = 1, dirac23_int_amp_factor
                        cur_i = ikr; cur_j = jkr; cur_k = indk(inz); cur_l = indl(inz)
                        cur_int2_real = rklr(inz)
                        if (.not. realonly%is_realonly()) cur_int2_imag = rkli(inz)
                        if (idx_amp_factor == 2) then
                            call swap(cur_i, cur_j); call swap(cur_k, cur_l)
                            if (is_type2_or_3()) then
                                ! Type2: (i j_bar|k l_bar) => (j i_bar|l k_bar) or
                                ! Type3: (i j_bar|k_bar l) => (j i_bar|l_bar k)
                                cur_i = -cur_i; cur_j = -cur_j; cur_k = -cur_k; cur_l = -cur_l
                            else
                                ! Type1 or 4: (i j|k l) => (j i|l k)*
                                if (.not. realonly%is_realonly()) cur_int2_imag = -cur_int2_imag
                            end if
                        end if

                        ! Convert Dirac(irrep order) index to CASPT2(energy order) index.
                        iii = indmo_dirac_to_cas(kr(cur_i))
                        jjj = indmo_dirac_to_cas(kr(cur_j))
                        kkk = indmo_dirac_to_cas(kr(cur_k))
                        lll = indmo_dirac_to_cas(kr(cur_l))

                        ! Create variables for Type1-4 determination.
                        ! The conversion is shown as the following example.
                        ! Before conversion: 1  2 3  4 5  6 7  8 9 10 11 12 13 14 15 16
                        ! After  conversion: 1 -1 2 -2 3 -3 4 -4 5 -5  6 -6  7 -7  8 -8
                        ! (-1)**(mod(iii, 2) + 1) returns 1 if iii is odd (Kramers+), -1 if iii is even (Kramers-).
                        ! (iii/2 + mod(iii, 2)) returns iii/2+1 if iii is odd, and iii/2 if iii is even.
                        ! Thus, it returns iii/2+1 if iii is odd and -iii/2 if iii is even, giving the desired conversion.
                        iikr = (-1)**(mod(iii, 2) + 1)*(iii/2 + mod(iii, 2))
                        jjkr = (-1)**(mod(jjj, 2) + 1)*(jjj/2 + mod(jjj, 2))
                        kkkr = (-1)**(mod(kkk, 2) + 1)*(kkk/2 + mod(kkk, 2))
                        llkr = (-1)**(mod(lll, 2) + 1)*(lll/2 + mod(lll, 2))

                        ! Since the signs of the indices of the integrals output by DIRAC and UTChem are different
                        ! (UTCHEM TYPE1:(++++), DIRAC TYPE1:(----)),
                        ! they are converted to the indices required for UTChem's two-electron integrals by taking Kramer pairs.
                        iiit = iii - (-1)**iii
                        jjjt = jjj - (-1)**jjj
                        kkkt = kkk - (-1)**kkk
                        lllt = lll - (-1)**lll

                        ! Create variables to determine the integrals which is required for UTChem type integrals.
                        ii = abs(iikr)
                        jj = abs(jjkr)
                        kk = abs(kkkr)
                        ll = abs(llkr)

                        !---------------------------
                        ! TYPE1 (++++) = (ij|kl)
                        ! TYPE2 (+-+-) = (ij~|kl~)
                        ! TYPE3 (+--+) = (ij~|k~l)
                        ! TYPE4 (+---) = (ij~|k~l~)
                        !---------------------------

                        if (should_write_2int_to_disk()) then
                            if (realonly%is_realonly()) then
                                write (unit_mdcintnew) iiit, jjjt, nnz, kkkt, lllt, cur_int2_real
                                write (unit_fmt_filename, '(4i6, e20.7)') cur_i, cur_j, cur_k, cur_l, cur_int2_real
                            else
                                write (unit_mdcintnew) iiit, jjjt, nnz, kkkt, lllt, cur_int2_real, -cur_int2_imag
                                write (unit_fmt_filename, '(4i6, 2e20.7)') cur_i, cur_j, cur_k, cur_l, cur_int2_real, -cur_int2_imag
                            end if
                        end if
                    end do
                End do
            end do
        end do mdcint_file_read
!--------------------------------- UTChem integral translation------------------------------------
!TYPE1       If( ((p10<=p20.and.p30<=p40.and.(p10<p30.or.(p10==p30.and.p20<=p40))) .or. &
!               (p10< p20.and.p40<=p30.and.(p10<p40.or.(p10==p40.and.p20<=p30))).and.(isp==isq)) &
!                .or. (isp/=isq) ) then
!TYPE2       If( ((p30<=p40.and.(p10<p30.or.(p10==p30.and.p20<=p40))).and.(isp==isq)) &
!                .or. (isp/=isq) ) then
!TYPE3       If( ((p30<=p40.and.(p10<p30.or.(p10==p30.and.p20<=p40))).and.(isp==isq)) &
!                .or. (isp/=isq) ) then
!TYPE4       If( ((isp==isq)) &
!              .or. (isp/=isq) ) then
!-------------------------------------------------------------------------------------------------

        close (unit_mdcint)
        file_idx = file_idx + 1
        call get_mdcint_filename(file_idx) ! Get the next MDCINT filename

    end do
    write (unit_mdcintnew) 0, 0, 0
    close (unit_mdcintnew)
    deallocate (indk)
    deallocate (indl)
    deallocate (rklr)
    if (allocated(rkli)) deallocate (rkli)

    if (rank == 0) print *, 'End create_newmdcint'
    deallocate (kr)
contains

    logical function is_type2_or_3()
        implicit none
        integer :: left_indices(2), right_indices(2)
        left_indices(1) = cur_i; left_indices(2) = cur_j
        right_indices(1) = cur_k; right_indices(2) = cur_l
        ! If the integral is TYPE2 or TYPE3, a negative index is included in the left_indices and right_indices.
        is_type2_or_3 = (count(left_indices < 0) == 1 .and. count(right_indices < 0) == 1)
    end function is_type2_or_3

    logical function should_write_2int_to_disk()
        implicit none
        should_write_2int_to_disk = .false. ! Do not write to disk by default

        ! If the integral is nearly zero, do not write to disk.
        if (realonly%is_realonly()) then
            if (abs(rklr(inz)) <= cutoff) then
                should_write_2int_to_disk = .false.
                return
            end if
        else
            if (abs(rklr(inz)) <= cutoff .and. abs(rkli(inz)) <= cutoff) then
                should_write_2int_to_disk = .false.
                return
            end if
        end if
        If (iikr > 0 .and. jjkr > 0 .and. kkkr > 0 .and. llkr > 0) then  !TYPE1
            if ((ii <= jj .and. kk <= ll .and. (ii < kk .or. (ii == kk .and. jj <= ll))) .or. &
                (ii <= jj .and. ll <= kk .and. (ii < ll .or. (ii == ll .and. jj <= kk)))) then
                should_write_2int_to_disk = .true.
            end if
        Else if (iikr > 0 .and. jjkr < 0 .and. kkkr > 0 .and. llkr < 0) then  !TYPE2
            if (ii <= jj .and. kk <= ll .and. (ii < kk .or. (ii == kk .and. jj <= ll))) then
                should_write_2int_to_disk = .true.
            end if
        Else if (iikr > 0 .and. jjkr < 0 .and. kkkr < 0 .and. llkr > 0) then  !TYPE3
            if (ii <= jj .and. kk <= ll .and. (ii < kk .or. (ii == kk .and. jj <= ll))) then
                should_write_2int_to_disk = .true.
            end if
        Else if (iikr > 0 .and. jjkr < 0 .and. kkkr < 0 .and. llkr < 0) then  !TYPE4
            if (ii <= jj) then
                should_write_2int_to_disk = .true.
            end if
        End if
    end function should_write_2int_to_disk
end Subroutine create_newmdcint
