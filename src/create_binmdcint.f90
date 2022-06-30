!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Subroutine create_newmdcint ! 2 Electorn Integrals In Mdcint

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Use four_caspt2_module
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
    real    :: cutoff
    integer :: nnkr, iiit, jjjt, kkkt, lllt
    integer :: nkr, nz, file_idx
    logical :: is_file_exist

    Call timing(date1, tsec1, date0, tsec0)
    date1 = date0
    tsec1 = tsec0
    casci_mdcint_cnt = 0
    Allocate (kr(-nmo/2:nmo/2))
    kr = 0
    ! Get datex, timex, nkr, and kr from MDCINT becasuse there is no kr information in the MDCINXXX files.
    if (rank == 0) then
        ! open (8, file="debug", form="formatted", status="unknown")
        open (10, file="MDCINT", form="unformatted", status="unknown")
        read (10) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
        close (10)
    end if
    Allocate (indk(nmo**2))
    Allocate (indl(nmo**2))
    Allocate (rklr(nmo**2))
    Allocate (rkli(nmo**2))

    if (rank == 0) then ! Process limits for output
        write (*, *) "allocate successed."
    end if
#ifdef HAVE_MPI
    ! Broadcast kr and other data that are not included in the MDCINXXX files
    call MPI_Bcast(datex, sizeof(datex), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then ! Process limits for output
        write (*, *) "datex broadcast"
        write (*, *) "if ierr == 0, datex broadcast successed. ierr=", ierr
    end if
    call MPI_Bcast(timex, sizeof(timex), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then ! Process limits for output
        write (*, *) "timex broadcast"
        write (*, *) "if ierr == 0, timex broadcast successed. ierr=", ierr
    end if
    call MPI_Bcast(nkr, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then ! Process limits for output
        write (*, *) "nkr broadcast"
        write (*, *) "if ierr == 0, nkr broadcast successed. ierr=", ierr
    end if
    call MPI_Bcast(kr(-nmo/2), nmo + 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then ! Process limits for output
        write (*, *) "kr broadcast"
        write (*, *) "if ierr == 0, kr broadcast successed. ierr=", ierr
    end if
    call MPI_Bcast(indmor(1), nmo, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then ! Process limits for output
        write (*, *) "datex broadcast"
        write (*, *) "if ierr == 0, datex broadcast successed. ierr=", ierr
    end if
#endif
    nnkr = 0

    realonly = .false.
    cutoff = 0.25D-12
    nnz = 1

    ! Initialization
    file_idx = 0
    is_file_exist = .true.
    do while (is_file_exist)
        call get_mdcint_filename(file_idx)

        inquire (file=mdcint_filename, exist=is_file_exist)
        if (is_file_exist) then
            open (100, file=mdcint_filename, form='unformatted', status='unknown')
            read (100)

            read (100, ERR=200) ikr, jkr, nz, (indk(inz), indl(inz), rklr(inz), rkli(inz), inz=1, nz)
            goto 201

200         realonly = .true.
            if (rank == 0) then ! Process limits for output
                write (*, *) "realonly = ", realonly
            end if
201         close (100)

            open (100, file=mdcint_filename, form='unformatted', status='unknown')
            read (100)
        end if
        if (file_idx == 0) then
            open (200, file=mdcintNew, form='unformatted', status='replace')
            write (200) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
        end if
        if (is_file_exist) then
            nnkr = nkr
            rkli = 0.0d+00

!Iwamuro debug
            ! write(*,*) "new_ikr1", datex, timex, nkr, (kr(i0),kr(-1*i0),i0=1,nkr)
            ! write(*,*) Filename

100         if (realonly) then
                read (100, end=1000) ikr, jkr, nz, &
                    (indk(inz), indl(inz), inz=1, nz), &
                    (rklr(inz), inz=1, nz)
            else
                read (100, end=1000) ikr, jkr, nz, &
                    (indk(inz), indl(inz), inz=1, nz), &
                    (rklr(inz), rkli(inz), inz=1, nz)
            end if

!------------------------------!
!   Create new ikr for UTChem  !
!------------------------------!

!       new ikr = iikr
!           jkr = jjkr
!           kkr = kkkr
!           lkr = llkr

!        Do inz = 1,nz

            if (ikr < 0) then
                if (rank == 0) then ! Process limits for output
                    write (*, *) "ikr<0. ikr=", ikr
                end if
!        go to 100
            end if
            if (ikr == 0) then
                if (rank == 0) then ! Process limits for output
                    write (*, *) ikr, jkr, nz, mdcint_debug
                end if
                ! write (200) 0, 0, 0
                ! write(29,'(3I4)') 0, 0, 0
                ! write(30,'(3I4)') 0, 0, 0
                go to 1000
            end if

            n = 1
10          select case (n)
            case (1)
                ikr = ikr
                jkr = jkr
                indk(:) = indk(:)
                indl(:) = indl(:)
                rklr(:) = rklr(:)
                rkli(:) = rkli(:)
            case (2)
                ikr = -ikr
                jkr = -jkr
                indk(:) = -indk(:)
                indl(:) = -indl(:)
                rklr(:) = -rklr(:)*sign(1, ikr*jkr*indk(:)*indl(:))
                rkli(:) = rkli(:)*sign(1, ikr*jkr*indk(:)*indl(:))
            end select
            !   !$OMP parallel do private(iii,jjj,kkk,lll,iikr,jjkr,kkkr,llkr,iiit,jjjt,kkkt,lllt,ii,jj,kk,ll)
            Do inz = 1, nz
                ! write (*, *) "Immediately after the declaration of do."
                ! write (300, "(5I20,E32.16)") ikr, jkr, nz, indk(inz), indl(inz), rklr(inz)
                ! Debug output (if write(*,*))
                ! if (inz == 1) then
                ! write(*,*)"new_ikr2"
                ! write(*,*)"Filename:", Filename
                ! write(*,*)"inz:", inz
                ! end if

                iii = indmor(kr(ikr))
                jjj = indmor(kr(jkr))
                kkk = indmor(kr(indk(inz)))
                lll = indmor(kr(indl(inz)))

                iikr = (-1)**(mod(iii, 2) + 1)*(iii/2 + mod(iii, 2))
                jjkr = (-1)**(mod(jjj, 2) + 1)*(jjj/2 + mod(jjj, 2))
                kkkr = (-1)**(mod(kkk, 2) + 1)*(kkk/2 + mod(kkk, 2))
                llkr = (-1)**(mod(lll, 2) + 1)*(lll/2 + mod(lll, 2))

                iiit = iii - (-1)**iii
                jjjt = jjj - (-1)**jjj
                kkkt = kkk - (-1)**kkk
                lllt = lll - (-1)**lll

! Iwamuro debug
                ! if (inz == 1) then
                !     ! write(*,*) "new_ikr2", iikr, jjkr, kkkr, llkr
                ! end if
! Debug output end (if write(*,*))

!------------------------------------------------------------

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

                If (iikr > 0 .and. jjkr > 0 .and. kkkr > 0 .and. llkr > 0) then  !TYPE1
                    if ((ii <= jj .and. kk <= ll .and. (ii < kk .or. (ii == kk .and. jj <= ll))) .or. &
                        (ii <= jj .and. ll <= kk .and. (ii < ll .or. (ii == ll .and. jj <= kk)))) then
                        if (abs(rklr(inz)) > cutoff .or. &
                            abs(rkli(inz)) > cutoff) then
                            write (200) iiit, jjjt, nnz, kkkt, lllt, rklr(inz), -(rkli(inz))
                        end if
                    end if

                Else if (iikr > 0 .and. jjkr < 0 .and. kkkr > 0 .and. llkr < 0) then  !TYPE2
                    if (ii <= jj .and. kk <= ll .and. (ii < kk .or. (ii == kk .and. jj <= ll))) then
                        if (abs(rklr(inz)) > cutoff .or. &
                            abs(rkli(inz)) > cutoff) then
                            write (200) iiit, jjjt, nnz, kkkt, lllt, rklr(inz), -(rkli(inz))
                        end if
                    end if

                Else if (iikr > 0 .and. jjkr < 0 .and. kkkr < 0 .and. llkr > 0) then  !TYPE3
                    if (ii <= jj .and. kk <= ll .and. (ii < kk .or. (ii == kk .and. jj <= ll))) then
                        if (abs(rklr(inz)) > cutoff .or. &
                            abs(rkli(inz)) > cutoff) then
                            write (200) iiit, jjjt, nnz, kkkt, lllt, rklr(inz), -(rkli(inz))
                        end if
                    end if

                Else if (iikr > 0 .and. jjkr < 0 .and. kkkr < 0 .and. llkr < 0) then  !TYPE4
                    if (ii <= jj) then
                        if (abs(rklr(inz)) > cutoff .or. &
                            abs(rkli(inz)) > cutoff) then
                            write (200) iiit, jjjt, nnz, kkkt, lllt, rklr(inz), -(rkli(inz))
                        end if
                    end if
                    ! else
                    ! write (200) - iikr, -jjkr, nnz, -kkkr, -llkr, rklr(inz), -(rkli(inz))
                    ! write (300, '(a6,5I4,2E32.16)') 'else', -iikr, -jjkr, nnz, -kkkr, -llkr, rklr(inz), -(rkli(inz))
                End if
                ! write(*, *) "before end , inz :", inz
            End do

            If (n == 1) then
                n = 2
                go to 10
            else
                go to 100
            end if

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

1000        close (100)
            file_idx = file_idx + 1
        end if
    end do
    write (200) 0, 0, 0
    close (200)
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    Call timing(date1, tsec1, date0, tsec0)
    date1 = date0
    tsec1 = tsec0
    deallocate (indk)
    deallocate (indl)
    deallocate (rklr)
    deallocate (rkli)

    if (rank == 0) then ! Process limits for output
        write (*, *) 'end create_binmdcint.'
    end if
    deallocate (kr)
end Subroutine create_newmdcint
