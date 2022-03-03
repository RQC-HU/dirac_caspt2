!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Subroutine create_newmdcint ! 2 Electorn Integrals In Mdcint

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Use four_caspt2_module
    ! use omp_lib
    Implicit None
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    Character*50 :: Filename

    Character  :: datex*10, timex*8
    integer :: i0, mdcint, inz, nnz, n
    integer :: ikr, jkr
    integer :: ii, jj, kk, ll
    integer :: iikr, jjkr, kkkr, llkr, iii, jjj, kkk, lll
    integer, allocatable :: indk(:), indl(:), kr(:)
    double precision, allocatable  :: rklr(:), rkli(:)
    ! Iwamuro modify
    real    :: cutoff
    integer :: nnkr, iiit, jjjt, kkkt, lllt
    real            :: time_start, time_end
    integer         :: nkr, nz, loopcnt

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
        write (normaloutput, *) "allocate successed. rank=", rank
    end if
#ifdef HAVE_MPI
    ! Broadcast kr and other data that are not included in the MDCINXXX files
    call MPI_Bcast(datex, sizeof(datex), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then ! Process limits for output
        write (normaloutput, '(a,i4)') "datex broadcast rank=", rank
        write (normaloutput, '(a,i4,a,i4)') "if ierr == 0, datex broadcast successed. ierr=", ierr, "rank=", rank
    end if
    call MPI_Bcast(timex, sizeof(timex), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then ! Process limits for output
        write (normaloutput, '(a,i4)') "timex broadcast rank=", rank
        write (normaloutput, '(a,i4,a,i4)') "if ierr == 0, timex broadcast successed. ierr=", ierr, "rank=", rank
    end if
    call MPI_Bcast(nkr, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then ! Process limits for output
        write (normaloutput, '(a,i4)') "nkr broadcast rank=", rank
        write (normaloutput, '(a,i4,a,i4)') "if ierr == 0, nkr broadcast successed. ierr=", ierr, "rank=", rank
    end if
    call MPI_Bcast(kr(-nmo/2), nmo + 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then ! Process limits for output
        write (normaloutput, '(a,i4)') "kr broadcast rank=", rank
        write (normaloutput, '(a,i4,a,i4)') "if ierr == 0, kr broadcast successed. ierr=", ierr, "rank=", rank
    end if
    call MPI_Bcast(indmor(1), nmo, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then ! Process limits for output
        write (normaloutput, '(a,i4)') "datex broadcast rank=", rank
        write (normaloutput, '(a,i4,a,i4)') "if ierr == 0, datex broadcast successed. ierr=", ierr, "rank=", rank
    end if
#endif
    nnkr = 0

    realonly = .false.
    cutoff = 0.25D-12
    nnz = 1

    if (rank == 0) then ! Process limits for output
        write (normaloutput, '(3a,i20)') "end set ", mdcintNew, "valiables. rank=", rank
    end if
    ! mdcint=11
    ! tid = omp_get_thread_num()
    open (100, file=mdcint_filename, form='unformatted', status='unknown')
!     if (rank == 0) then
!         read (100) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
!     else
    read (100)
!     end if

    read (100, ERR=200) ikr, jkr, nz, (indk(inz), indl(inz), rklr(inz), rkli(inz), inz=1, nz)
    ! read (mdcint,ERR=200) ikr,jkr, nz, (indk(inz),indl(inz), rklr(inz),rkli(inz), inz=1,nz)
    goto 201

200 realonly = .true.
    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) "realonly = ", realonly, rank
    end if
201 close (100)

    open (100, file=mdcint_filename, form='unformatted', status='unknown')
    open (200, file=mdcintNew, form='unformatted', status='replace')
    ! open (300, file=mdcint_debug, form='formatted', status='unknown')
    ! open(30, file=mdcint_int, form='formatted', status='unknown')
    ! write (normaloutput, *) "before read day,time,kr. rank=", rank
    read (100)
    ! write (normaloutput, *) "end read day,time,kr. rank=", rank

    write (200) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
    ! write (300, *) datex, timex, nkr, (kr(i0), kr(-1*i0), i0=1, nkr)
    nnkr = nkr
    rkli = 0.0d+00

!Iwamuro debug
    ! write(*,*) "new_ikr1", datex, timex, nkr, (kr(i0),kr(-1*i0),i0=1,nkr)
    ! write(*,*) Filename
    ! write (normaloutput, *) "before read ikr,jkr,...and more. rank=", rank

100 if (realonly) then
        read (100, end=1000) ikr, jkr, nz, &
            (indk(inz), indl(inz), inz=1, nz), &
            (rklr(inz), inz=1, nz)
        ! rkli = 0.0d+00
    else
        read (100, end=1000) ikr, jkr, nz, &
            (indk(inz), indl(inz), inz=1, nz), &
            (rklr(inz), rkli(inz), inz=1, nz)
    end if
    ! write (normaloutput, '(A,4I6)') "end read ikr,jkr,...and more. rank=", rank, ikr, jkr, nz

! Debug output
    ! write(*,*) ""
    ! write(*,*) ikr,jkr, nz, &
    !                 (indk(inz),indl(inz),inz=1,nz), &
    !                 (rklr(inz),rkli(inz), inz=1,nz)

    ! write(*,*) ""
    ! loop = loop + 1
! Debug output end
!------------------------------------------------------------

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
            write (normaloutput, *) "ikr<0. rank=", rank, "ikr=", ikr
        end if
        go to 100
    end if
    if (ikr == 0) then
        if (rank == 0) then ! Process limits for output
            write (normaloutput, *) ikr, jkr, nz, mdcint_debug
        end if
        write (200) 0, 0, 0
        ! write(29,'(3I4)') 0, 0, 0
        ! write(30,'(3I4)') 0, 0, 0
        go to 1000
    end if

    n = 1
10  select case (n)
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
        ! write (normaloutput, *) "Immediately after the declaration of do. rank=", rank
        ! write (300, "(5I20,E32.16)") ikr, jkr, nz, indk(inz), indl(inz), rklr(inz)
        ! Debug output (if write(*,*))
        ! if (inz == 1) then
        ! write(*,*)"new_ikr2"
        ! write(*,*)"Filename:", Filename
        ! write(*,*)"inz:", inz
        ! end if

        ! iii = indmor(kr(ikr8))
        iii = indmor(kr(ikr))
        ! if (rank == 1 .and. inz == 1) then
        !     ! write(*,*) "iii",ikr,ikr8,iii,(-1)**(mod(iii,2)+1)*(iii/2+mod(iii,2))
        !     ! write(*,*) "kr(ikr)", kr(ikr8)
        !     ! write(*,*) "indmor(kr(ikr))", indmor(kr(ikr))
        ! end if

        ! jjj = indmor(kr(jkr8))
        jjj = indmor(kr(jkr))
        ! if (inz == 1) then
        !     ! write(*,*) "kr(jkr)", kr(jkr)
        !     ! write(*,*) "indmor(kr(jkr))", indmor(kr(jkr))
        ! end if

        ! kkk = indmor(kr(indk8(inz)))
        kkk = indmor(kr(indk(inz)))
        ! if (inz == 1) then
        !     ! write(*,*) "indk(inz)", indk(inz)
        !     ! write(*,*) "kr(indk(inz))", kr(indk(inz))
        !     ! write(*,*) "indmor(kr(indk(inz)))", indmor(kr(indk(inz)))
        ! end if

        ! lll = indmor(kr(indl8(inz)))
        lll = indmor(kr(indl(inz)))
        ! if (inz == 1) then
        !     ! write(*,*) "indl(inz)", indl(inz)
        !     ! write(*,*) "kr(indl(inz))", kr(indl(inz))
        !     ! write(*,*) "indmor(kr(indl(inz)))", indmor(kr(indl(inz)))
        ! end if

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
!                    write(28) -ikr,-jkr,nnz,-(indk(inz)),-(indl(inz)),rklr(inz),-(rkli(inz))
!                    write(28) -iikr,-jjkr,nnz,-kkkr,-llkr,rklr8(inz),-(rkli8(inz))
                    ! write (200) iiit, jjjt, nnz, kkkt, lllt, rklr8(inz), -(rkli8(inz))
                    write (200) iiit, jjjt, nnz, kkkt, lllt, rklr(inz), -(rkli(inz))
                    ! write (300, '(5I4,2E32.16)') iiit, jjjt, nnz, kkkt, lllt, rklr8(inz), -(rkli8(inz))
                    ! casci_mdcint_cnt = casci_mdcint_cnt + 1
                    ! write(29,'(5I4,2E32.16)') -iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
                    ! else
                    ! write(29,'(a6,5I4,2E32.16)')'else1',-iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
                    ! write(30,'(5I4,2E32.16)') iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
                end if
            end if

        Else if (iikr > 0 .and. jjkr < 0 .and. kkkr > 0 .and. llkr < 0) then  !TYPE2
            if (ii <= jj .and. kk <= ll .and. (ii < kk .or. (ii == kk .and. jj <= ll))) then
                if (abs(rklr(inz)) > cutoff .or. &
                    abs(rkli(inz)) > cutoff) then
!                    write(28) -ikr,-jkr,nnz,-(indk(inz)),-(indl(inz)),rklr(inz),-(rkli(inz))
!                    write(28) -iikr,-jjkr,nnz,-kkkr,-llkr,rklr8(inz),-(rkli8(inz))
                    ! write (200) iiit, jjjt, nnz, kkkt, lllt, rklr8(inz), -(rkli8(inz))
                    write (200) iiit, jjjt, nnz, kkkt, lllt, rklr(inz), -(rkli(inz))
                    ! write (300, '(5I4,2E32.16)') iiit, jjjt, nnz, kkkt, lllt, rklr8(inz), -(rkli8(inz))
                    ! casci_mdcint_cnt = casci_mdcint_cnt + 1

                    ! write(29,'(5I4,2E32.16)') -iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
                    ! write(30,'(5I4,2E32.16)') iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
                    ! else
                    !         write(29,'(a6,5I4,2E32.16)')'else2',-iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
                end if
            end if

        Else if (iikr > 0 .and. jjkr < 0 .and. kkkr < 0 .and. llkr > 0) then  !TYPE3
            if (ii <= jj .and. kk <= ll .and. (ii < kk .or. (ii == kk .and. jj <= ll))) then
                if (abs(rklr(inz)) > cutoff .or. &
                    abs(rkli(inz)) > cutoff) then
!                    write(28) -ikr,-jkr,nnz,-(indk(inz)),-(indl(inz)),rklr(inz),-(rkli(inz))
!                    write(28) -iikr,-jjkr,nnz,-kkkr,-llkr,rklr8(inz),-(rkli8(inz))
                    ! write (200) iiit, jjjt, nnz, kkkt, lllt, rklr8(inz), -(rkli8(inz))
                    write (200) iiit, jjjt, nnz, kkkt, lllt, rklr(inz), -(rkli(inz))
                    ! write (300, '(5I4,2E32.16)') iiit, jjjt, nnz, kkkt, lllt, rklr8(inz), -(rkli8(inz))
                    ! casci_mdcint_cnt = casci_mdcint_cnt + 1

                    ! write(29,'(5I4,2E32.16)') -iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
                    ! write(30,'(5I4,2E32.16)') iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
                    ! else
                    ! write(29,'(a6,5I4,2E32.16)')'else3',-iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
                end if
            end if

        Else if (iikr > 0 .and. jjkr < 0 .and. kkkr < 0 .and. llkr < 0) then  !TYPE4
            if (ii <= jj) then
                if (abs(rklr(inz)) > cutoff .or. &
                    abs(rkli(inz)) > cutoff) then
!                    write(28) -ikr,-jkr,nnz,-(indk(inz)),-(indl(inz)),rklr(inz),-(rkli(inz))
!                    write(28) -iikr,-jjkr,nnz,-kkkr,-llkr,rklr8(inz),-(rkli8(inz))
                    ! write (200) iiit, jjjt, nnz, kkkt, lllt, rklr8(inz), -(rkli8(inz))
                    write (200) iiit, jjjt, nnz, kkkt, lllt, rklr(inz), -(rkli(inz))
                    ! write (300, '(5I4,2E32.16)') iiit, jjjt, nnz, kkkt, lllt, rklr8(inz), -(rkli8(inz))
                    ! casci_mdcint_cnt = casci_mdcint_cnt + 1

                    ! write(29,'(5I4,2E32.16)') -iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
                    ! write(30,'(5I4,2E32.16)') iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
                    ! else
                    ! write(29,'(a6,5I4,2E32.16)')'else4',-iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
                end if
            end if
            ! else
            ! write (200) - iikr, -jjkr, nnz, -kkkr, -llkr, rklr(inz), -(rkli(inz))
            ! write (300, '(a6,5I4,2E32.16)') 'else', -iikr, -jjkr, nnz, -kkkr, -llkr, rklr(inz), -(rkli(inz))
        End if
        ! write(normaloutput, *) "before end , inz :", inz
    End do

    If (n == 1) then
        n = 2
        go to 10
    else
        go to 100
    end if
!   go to 100

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

1000 close (100)
    close (200)
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    Call timing(date1, tsec1, date0, tsec0)
    date1 = date0
    tsec1 = tsec0
    ! close (300)
    ! close(30)
    deallocate (indk)
    deallocate (indl)
    deallocate (rklr)
    deallocate (rkli)

    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'end create_binmdcint. rank=', rank
    end if
    ! Debug output for casci_mdcint_cnt
    ! do loopcnt = 0, nprocs - 1
    !     if (loopcnt == rank) then
    !         write (*, *) 'casci_mdcint_cnt : ', rank, casci_mdcint_cnt
    !     end if
    !     call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! end do
    deallocate (kr)
end Subroutine create_newmdcint
