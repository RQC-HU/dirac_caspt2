!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Subroutine create_newmdcint ! 2 Electorn Integrals In Mdcint
program create_newmdcint

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        Use four_caspt2_module, only: indmor
        ! use omp_lib
        Implicit None
        include 'mpif.h'

        Character*50 :: Filename

        Character  :: datex*10, timex*8

        integer :: nkr, nz
        integer :: i0, mdcint, inz, nnz
        integer :: ikr, jkr
        integer :: ii, jj, kk, ll
        integer :: iikr, jjkr, kkkr, llkr, iii, jjj, kkk, lll
        integer, allocatable :: indk(:), indl(:), kr(:)
        double precision, allocatable  :: rklr(:), rkli(:)
        ! Iwamuro modify
        real    :: cutoff
        integer :: nnkr, ikr8, jkr8, iiit, jjjt, kkkt, lllt
        integer, allocatable :: kkr(:), indk8(:), indl8(:)
        real*8, allocatable  :: rklr8(:), rkli8(:)

        ! integer         :: i, loop = 0, omp_max, tid
        Character*50    :: fileBaseName, mdcintBaseName, mdcintNew, mdcint_debug, mdcint_int, mdcintNum
        integer         :: ierr, nprocs, rank, procs = 8 !! TODO MPI PROCS を動的に?設定する
        real            :: time_start, time_end
        open(8,file="debug",form="formatted",status="unknown")
        ! omp_max = omp_get_max_threads()
        Allocate(kr(-nmo/2:nmo/2))
        kr = 0
        open(10,file="MDCINT",form="unformatted",status="unknown")
        read (10) datex,timex,nkr,(kr(i0),kr(-1*i0),i0=1,nkr)
        open(11,file="firstline_mdcint",form="formatted",status="unknown")
        write(11,*) datex,timex,nkr,(kr(i0),kr(-1*i0),i0=1,nkr)
        close(11)
        close(10)
        Allocate(indk(nmo**2))
        Allocate(indl(nmo**2))
        Allocate(rklr(nmo**2))
        Allocate(rkli(nmo**2))
        Allocate(rklr8(nmo**2))
        Allocate(rkli8(nmo**2))
        Allocate(indk8(nmo**2))
        Allocate(indl8(nmo**2))
        Allocate(kkr(-nmo/2:nmo/2))
        write(8,*)"allocate successed."
        kkr = 0
        nnkr = 0

        call MPI_INIT(ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
        call MPI_COMM_rank(MPI_COMM_WORLD, rank, ierr)
        !$omp parallel do private(i0,indk,indl,rklr,rkli,rklr8,rkli8,indk8,indl8,kkr,Filename,mdcintNew,mdcint_debug,datex,timex,nkr,ikr,jkr,nz,inz,realonly,iii,jjj,kkk,lll,iikr,jjkr,kkkr,llkr,nnkr,ikr8,jkr8,iiit,jjjt,kkkt,lllt,mdcint_int,mdcintNum)
        ! Do i=1,8 ! TODO MDCINTファイルの数に応じてloopの数を変更する
        realonly = .true.
        cutoff = 0.25D-12
        nnz = 1



        fileBaseName = "MDCINXXXX"
        if (rank == 0) then
                Filename = "MDCINT"
                mdcintNew = "MDCINTNEW"
                mdcint_debug = "MDCINT_debug"
                mdcint_int = "MDCINT_int"
        else
                mdcintBaseName = "MDCINXXXX"
                write(mdcintNum,"(I3)") rank
                Filename = TRIM(mdcintBaseName)//TRIM(ADJUSTL(mdcintNum))
                mdcintNew = "MDCINTNEW"//TRIM(ADJUSTL(mdcintNum))
                mdcint_debug = "MDCINT_debug"//TRIM(ADJUSTL(mdcintNum))
                mdcint_int = "MDCINT_int"//TRIM(ADJUSTL(mdcintNum))
        end if
        write(8,*) "set mdcint valiables"
        ! mdcint=11
        ! tid = omp_get_thread_num()
        open(rank+100, file=Filename, form ='unformatted', status='unknown')
        ! open(mdcint, file=Filename, form ='unformatted', status='unknown')
        ! open(mdcint, file="MDCINT", form ='unformatted', status='unknown')
        ! if (rank == 0) then
        ! read (rank+100) datex,timex,nkr,(kr(i0),kr(-1*i0),i0=1,nkr)
        ! read (mdcint) datex,timex,nkr,(kr(i0),kr(-1*i0),i0=1,nkr)
        ! else
        read (rank+100)
        ! read (mdcint)
        ! end if
        write(8,*) "read day,time,kr"
        read (rank+100,ERR=200) ikr,jkr, nz, (indk(inz),indl(inz), rklr(inz),rkli(inz), inz=1,nz)
        ! read (mdcint,ERR=200) ikr,jkr, nz, (indk(inz),indl(inz), rklr(inz),rkli(inz), inz=1,nz)
        goto 201

200     realonly = .true.
        write(*,*) "realonly = ", realonly
201     close(rank+100)
! 201     close(mdcint)

        ! open(mdcint, file="MDCINT", form='unformatted', status='unknown')
        ! open(28, file="MDCINTNEW", form='unformatted', status='unknown')
        ! open(29, file="MDCINT_debug", form='formatted', status='unknown')
        ! open(30, file="MDCINT_int", form='formatted', status='unknown')

        open(rank+100, file=Filename, form='unformatted', status='unknown')
        ! open(mdcint, file=Filename, form='unformatted', status='unknown')
        open(rank+200, file=mdcintNew, form='unformatted', status='unknown')
        open(rank+300, file=mdcint_debug, form='formatted', status='unknown')
        ! open(30, file=mdcint_int, form='formatted', status='unknown')
        ! if (rank == 0) then
        ! read (rank+100) datex,timex,nkr, (kr(i0),kr(-1*i0),i0=1,nkr)
        ! else
        read(rank+100)
        ! end if
        write(rank+300,*) rank,datex,timex,nkr, (kr(i0),kr(-1*i0),i0=1,nkr)
        nnkr = nkr
        kkr(:) = kr(:)

        ! write(28) datex, timex, nnkr, (kkr(i0),kkr(-1*i0),i0=1,nnkr)
        ! write(29,*) datex, timex, nnkr, (kkr(i0),kkr(-1*i0),i0=1,nnkr)
!Iwamuro debug
        ! write(*,*) "new_ikr1", datex, timex, nkr, (kr(i0),kr(-1*i0),i0=1,nkr)
        ! write(*,*) Filename

100     if (realonly) then
        read (rank+100,end=1000) ikr,jkr, nz, &
                (indk(inz),indl(inz),inz=1,nz), &
                (rklr(inz), inz=1,nz)
        rkli = 0.0d+00
        else
        read (rank+100,end=1000) ikr,jkr, nz, &
                (indk(inz),indl(inz),inz=1,nz), &
                (rklr(inz),rkli(inz), inz=1,nz)
        endif
        write(8,*) "read ikr,jkr,nz"
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

        if (ikr<0) go to 100
        if (ikr == 0) then
        write(20,*)ikr,jkr,nz,mdcint_debug
        write(rank+200) 0, 0, 0
        ! write(29,'(3I4)') 0, 0, 0
        ! write(30,'(3I4)') 0, 0, 0
        go to 1000
        endif

        ikr8 = ikr
        jkr8 = jkr
        indk8(:) = indk(:)
        indl8(:) = indl(:)
        rklr8(:) = rklr(:)
        rkli8(:) = rkli(:)
        write(8,*) "Before do"
        Do inz = 1,nz
        write(8,*) "after do"
        ! write (rank + 300, "(5I20,E32.16)") ikr, jkr, nz, indk(inz), indl(inz), rklr(inz)
!         ! Debug output (if write(*,*))
!         if (inz == 1) then
!         ! write(*,*)"new_ikr2"
!         ! write(*,*)"Filename:", Filename
!         ! write(*,*)"inz:", inz
!         endif

        ! write(8,*) "before iii set"
        ! write(8,*) "kr test", kr(ikr8)
        ! iii = indmor(kr(ikr8))
        ! write(8,*) "iii set", iii
!         if (rank == 1 .and. inz == 1) then
!         ! write(*,*) "iii",ikr,ikr8,iii,(-1)**(mod(iii,2)+1)*(iii/2+mod(iii,2))
!         ! write(*,*) "kr(ikr)", kr(ikr8)
!         ! write(*,*) "indmor(kr(ikr))", indmor(kr(ikr))
!         endif

!         jjj = indmor(kr(jkr8))
!         write(8,*) "jjj set"
!         if (inz == 1) then
!         ! write(*,*) "kr(jkr)", kr(jkr)
!         ! write(*,*) "indmor(kr(jkr))", indmor(kr(jkr))
!         endif

!         kkk = indmor(kr(indk8(inz)))
!         write(8,*) "kkk set"
!         if (inz == 1) then
!         ! write(*,*) "indk(inz)", indk(inz)
!         ! write(*,*) "kr(indk(inz))", kr(indk(inz))
!         ! write(*,*) "indmor(kr(indk(inz)))", indmor(kr(indk(inz)))
!         endif

!         lll = indmor(kr(indl8(inz)))
!         write(8,*) "lll set"
!         if (inz == 1) then
!         ! write(*,*) "indl(inz)", indl(inz)
!         ! write(*,*) "kr(indl(inz))", kr(indl(inz))
!         ! write(*,*) "indmor(kr(indl(inz)))", indmor(kr(indl(inz)))
!         endif

!         iikr = (-1)**(mod(iii,2)+1)*(iii/2+mod(iii,2))
!         jjkr = (-1)**(mod(jjj,2)+1)*(jjj/2+mod(jjj,2))
!         kkkr = (-1)**(mod(kkk,2)+1)*(kkk/2+mod(kkk,2))
!         llkr = (-1)**(mod(lll,2)+1)*(lll/2+mod(lll,2))
!         write(8,*) "iikr~llkr set"

!         iiit = iii-(-1)**iii
!         jjjt = jjj-(-1)**jjj
!         kkkt = kkk-(-1)**kkk
!         lllt = lll-(-1)**lll
!         write(8,*) "iiit~lllt set"


! ! Iwamuro debug
!         if (inz == 1) then
!         ! write(*,*) "new_ikr2", iikr, jjkr, kkkr, llkr
!         endif
! ! Debug output end (if write(*,*))

! !------------------------------------------------------------

!         ii = abs(iikr)
!         jj = abs(jjkr)
!         kk = abs(kkkr)
!         ll = abs(llkr)
!         write(8,*) "ii~ll set"

!         !---------------------------
!         ! TYPE1 (++++) = (ij|kl)
!         ! TYPE2 (+-+-) = (ij~|kl~)
!         ! TYPE3 (+--+) = (ij~|k~l)
!         ! TYPE4 (+---) = (ij~|k~l~)
!         !---------------------------

!         If(iikr>0 .and. jjkr>0 .and. kkkr>0 .and. llkr>0) then  !TYPE1
!         if( (ii<=jj .and. kk<=ll .and. (ii<kk .or. (ii==kk .and. jj<=ll))) .or. &
!                 (ii<=jj .and. ll<=kk .and. (ii<ll .or. (ii==ll .and. jj<=kk))) ) then
!                 if(abs(rklr(inz))>cutoff.or. &
!                 abs(rkli(inz))>cutoff      ) then
! !                    write(28) -ikr,-jkr,nnz,-(indk(inz)),-(indl(inz)),rklr(inz),-(rkli(inz))
! !                    write(28) -iikr,-jjkr,nnz,-kkkr,-llkr,rklr8(inz),-(rkli8(inz))
!                 write(rank+200) iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
!                 ! write(29,'(5I4,2E32.16)') -iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
!                 ! else
!                 ! write(29,'(a6,5I4,2E32.16)')'else1',-iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
!                 ! write(30,'(5I4,2E32.16)') iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
!                 endif
!         endif

!         Else if (iikr>0 .and. jjkr<0 .and. kkkr>0 .and. llkr<0) then  !TYPE2
!         if (ii<=jj .and. kk<=ll .and. (ii<kk .or. (ii==kk .and. jj<=ll))) then
!                 if(abs(rklr(inz))>cutoff.or. &
!                 abs(rkli(inz))>cutoff      ) then
! !                    write(28) -ikr,-jkr,nnz,-(indk(inz)),-(indl(inz)),rklr(inz),-(rkli(inz))
! !                    write(28) -iikr,-jjkr,nnz,-kkkr,-llkr,rklr8(inz),-(rkli8(inz))
!                 write(rank+200) iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
!                 ! write(29,'(5I4,2E32.16)') -iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
!                 ! write(30,'(5I4,2E32.16)') iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
!                 ! else
!                 !         write(29,'(a6,5I4,2E32.16)')'else2',-iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
!                 endif
!         endif

!         Else if (iikr>0 .and. jjkr<0 .and. kkkr<0 .and. llkr>0) then  !TYPE3
!         if(ii<=jj .and. kk<=ll .and. (ii<kk .or. (ii==kk .and. jj<=ll))) then
!                 if(abs(rklr(inz))>cutoff.or. &
!                 abs(rkli(inz))>cutoff      ) then
! !                    write(28) -ikr,-jkr,nnz,-(indk(inz)),-(indl(inz)),rklr(inz),-(rkli(inz))
! !                    write(28) -iikr,-jjkr,nnz,-kkkr,-llkr,rklr8(inz),-(rkli8(inz))
!                 write(rank+200) iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
!                 ! write(29,'(5I4,2E32.16)') -iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
!                 ! write(30,'(5I4,2E32.16)') iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
!                 ! else
!                 ! write(29,'(a6,5I4,2E32.16)')'else3',-iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
!                 endif
!         endif

!         Else if (iikr>0 .and. jjkr<0 .and. kkkr<0 .and. llkr<0) then  !TYPE4
!         if(ii<=jj) then
!                 if(abs(rklr(inz))>cutoff.or. &
!                 abs(rkli(inz))>cutoff      ) then
! !                    write(28) -ikr,-jkr,nnz,-(indk(inz)),-(indl(inz)),rklr(inz),-(rkli(inz))
! !                    write(28) -iikr,-jjkr,nnz,-kkkr,-llkr,rklr8(inz),-(rkli8(inz))
!                 write(rank+200) iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
!                 ! write(29,'(5I4,2E32.16)') -iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
!                 ! write(30,'(5I4,2E32.16)') iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
!                 ! else
!                 ! write(29,'(a6,5I4,2E32.16)')'else4',-iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
!                 endif
!         endif
!         ! else
!                 ! write(28,'(a6,5I4,2E32.16)')'else',-iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
!         Endif
300     Enddo
        write(8,*) "End do"


        go to 100

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

1000            close(rank+100)
                close(rank+200)
                close(rank+300)
                ! close(30)
                write(20,*) "1000 closed "//trim(Filename)
                deallocate(indk)
                deallocate(indl)
                deallocate(rklr)
                deallocate(rkli)
                deallocate(kkr)
                deallocate(indk8)
                deallocate(indl8)
                deallocate(rklr8)
                deallocate(rkli8)
                call MPI_FINALIZE(ierr)
        ! end do
        !$omp end parallel do
        deallocate(kr)
        time_end = mpi_wtime()
        open(10,file="output",form="formatted",status="unknown")
        write (10, *) "It took ", time_end - time_start, " seconds."
        close(10)
        close(8)
! end Subroutine create_newmdcint
end program create_newmdcint
