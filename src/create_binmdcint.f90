!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Subroutine create_newmdcint ! 2 Electorn Integrals In Mdcint

! This subroutine creates UTChem type MDCINTNEW files from DIRAC type MDCINT files.
! MDCINTNEW files are used in this program as 2-electron integrals.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use module_file_manager
    use module_realonly, only: realonly
    Use module_global_variables
    ! use omp_lib
    Implicit None
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    Character  :: datex*10, timex*8
    integer :: i0, inz, nnz, n
    integer :: ikr, jkr
    integer :: ii, jj, kk, ll
    integer :: iikr, jjkr, kkkr, llkr, iii, jjj, kkk, lll
    integer, allocatable :: indk(:), indl(:), kr(:)
    double precision, allocatable  :: rklr(:), rkli(:)
    ! Iwamuro modify
    real(8) :: cutoff
    integer :: iiit, jjjt, kkkt, lllt
    integer :: nkr, nz, file_idx, iostat
    integer :: unit_mdcint, unit_mdcintnew
    logical :: is_file_exist, is_end_of_file

    Call timing(date1, tsec1, date0, tsec0)
    date1 = date0
    tsec1 = tsec0
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
    if (rank == 0) then
        print *, "datex broadcast"
        print *, "if ierr == 0, datex broadcast successed. ierr=", ierr
    end if
    call MPI_Bcast(timex, sizeof(timex), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then
        print *, "timex broadcast"
        print *, "if ierr == 0, timex broadcast successed. ierr=", ierr
    end if
    call MPI_Bcast(nkr, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then
        print *, "nkr broadcast"
        print *, "if ierr == 0, nkr broadcast successed. ierr=", ierr
    end if
    call MPI_Bcast(kr(-nmo/2), nmo + 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then
        print *, "kr broadcast"
        print *, "if ierr == 0, kr broadcast successed. ierr=", ierr
    end if
    call MPI_Bcast(indmo_dirac_to_cas(1), nmo, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then
        print *, "datex broadcast"
        print *, "if ierr == 0, datex broadcast successed. ierr=", ierr
    end if
#endif

    cutoff = 0.25D-12
    nnz = 1

    ! Initialization
    file_idx = 0

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
                read (unit_mdcint, iostat=iostat) ikr, jkr, nz, (indk(inz), indl(inz), inz=1, nz), &
                    (rklr(inz), inz=1, nz)
            else
                read (unit_mdcint, iostat=iostat) ikr, jkr, nz, (indk(inz), indl(inz), inz=1, nz), &
                    (rklr(inz), rkli(inz), inz=1, nz)
            end if

            call check_iostat(iostat=iostat, file=mdcint_filename, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit mdcint_file_read
            end if

!------------------------------!
!   Create new ikr for UTChem  !
!------------------------------!

!       new ikr = iikr
!           jkr = jjkr
!           kkr = kkkr
!           lkr = llkr

            if (ikr == 0) then
                if (rank == 0) print *, ikr, jkr, nz, mdcint_debug
                exit mdcint_file_read ! End of file
            end if

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
                Do inz = 1, nz

                    ! Convert Dirac(irrep order) index to CASPT2(energy order) index.
                    iii = indmo_dirac_to_cas(kr(ikr))
                    jjj = indmo_dirac_to_cas(kr(jkr))
                    kkk = indmo_dirac_to_cas(kr(indk(inz)))
                    lll = indmo_dirac_to_cas(kr(indl(inz)))

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
                            write (unit_mdcintnew) iiit, jjjt, nnz, kkkt, lllt, rklr(inz)
                        else
                            write (unit_mdcintnew) iiit, jjjt, nnz, kkkt, lllt, rklr(inz), -(rkli(inz))
                        end if
                    end if
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
    Call timing(date1, tsec1, date0, tsec0)
    date1 = date0
    tsec1 = tsec0
    deallocate (indk)
    deallocate (indl)
    deallocate (rklr)
    if (allocated(rkli)) deallocate (rkli)

    if (rank == 0) print *, 'end create_binmdcint.'
    deallocate (kr)
contains
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
