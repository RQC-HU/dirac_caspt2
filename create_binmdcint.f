!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Subroutine create_newmdcint ! 2 Electorn Integrals In Mdcint

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        Use four_caspt2_module

        Implicit None

        Character*50 :: Filename

        Character  :: datex*10, timex*8

        integer :: nkr, nz
        integer    :: i0, mdcint, inz, nnz
        integer :: ikr, jkr
        integer    :: ii, jj, kk, ll
        integer    :: iikr, jjkr, kkkr, llkr, iii, jjj, kkk, lll
        integer, allocatable :: indk(:), indl(:), kr(:)
        double precision, allocatable  :: rklr(:), rkli(:)
        ! Iwamuro modify
        real    :: cutoff
        integer    :: nnkr, ikr8, jkr8, iiit, jjjt, kkkt, lllt
        integer, allocatable :: kkr(:), indk8(:), indl8(:)
        real*8, allocatable  :: rklr8(:), rkli8(:)

        integer         :: i, digit
        Character*50     :: fileBaseName, mdcintBaseName, mdcintNew, mdcint_debug, mdcint_int, mdcintNum


        Do i=1,8 ! TODO MDCINTファイルの数に応じてloopの数を変更する
        realonly = .false.
        cutoff = 0.25D-12
        nnz = 1

        Allocate(kr(-nmo/2:nmo/2))
        Allocate(indk(nmo**2))
        Allocate(indl(nmo**2))
        Allocate(rklr(nmo**2))
        Allocate(rkli(nmo**2))
        Allocate(rklr8(nmo**2))
        Allocate(rkli8(nmo**2))
        Allocate(indk8(nmo**2))
        Allocate(indl8(nmo**2))
        Allocate(kkr(-nmo/2:nmo/2))

        kr = 0
        kkr = 0
        nnkr = 0
        fileBaseName = "MDCINXXXX"
        if (i == 1) then
                Filename = "MDCINT"
                mdcintNew = "MDCINTNEW"
                mdcint_debug = "MDCINT_debug"
                mdcint_int = "MDCINT_int"
        else
                mdcintBaseName = "MDCINXXXX"
                write(mdcintNum,"(I3)") i-1
                Filename = TRIM(mdcintBaseName)//TRIM(ADJUSTL(mdcintNum))
                mdcintNew = "MDCINTNEW"//TRIM(ADJUSTL(mdcintNum))
                mdcint_debug = "MDCINT_debug"//TRIM(ADJUSTL(mdcintNum))
                mdcint_int = "MDCINT_int"//TRIM(ADJUSTL(mdcintNum))
        end if
        mdcint=11
        open(mdcint, file=Filename, form ='unformatted', status='unknown')
        ! open(mdcint, file="MDCINT", form ='unformatted', status='unknown')

        read (mdcint) datex,timex,nkr,(kr(i0),kr(-1*i0),i0=1,nkr)

        read (mdcint,ERR=200) ikr,jkr, nz, (indk(inz),indl(inz), rklr(inz),rkli(inz), inz=1,nz)

        goto 201

200     realonly = .true.
        write(*,*) "realonly = ", realonly
201     close(mdcint)

        ! open(mdcint, file="MDCINT", form='unformatted', status='unknown')
        ! open(28, file="MDCINTNEW", form='unformatted', status='unknown')
        ! open(29, file="MDCINT_debug", form='formatted', status='unknown')
        ! open(30, file="MDCINT_int", form='formatted', status='unknown')


        open(mdcint, file=Filename, form='unformatted', status='unknown')
        open(28, file=mdcintNew, form='unformatted', status='unknown')
        open(29, file=mdcint_debug, form='formatted', status='unknown')
        open(30, file=mdcint_int, form='formatted', status='unknown')

        read (mdcint) datex,timex,nkr, (kr(i0),kr(-1*i0),i0=1,nkr)

        nnkr = nkr
        kkr(:) = kr(:)

        write(28) datex, timex, nnkr, (kkr(i0),kkr(-1*i0),i0=1,nnkr)
!        write(29,*) datex, timex, nnkr, (kkr(i0),kkr(-1*i0),i0=1,nnkr)
!Iwamuro debug
        ! write(*,*) "new_ikr1", datex, timex, nkr, (kr(i0),kr(-1*i0),i0=1,nkr)
        ! write(*,*) Filename

100     if (realonly) then
        read (mdcint) ikr,jkr, nz, &
                (indk(inz),indl(inz),inz=1,nz), &
                (rklr(inz), inz=1,nz)
        rkli = 0.0d+00
        else
        read (mdcint) ikr,jkr, nz, &
                (indk(inz),indl(inz),inz=1,nz), &
                (rklr(inz),rkli(inz), inz=1,nz)
        endif
! Debug output
        write(*,*) ""
        write(*,*) ikr,jkr, nz, &
                        (indk(inz),indl(inz),inz=1,nz), &
                        (rklr(inz),rkli(inz), inz=1,nz)

        write(*,*) ""
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
        write(28) 0, 0, 0
        write(29,'(3I4)') 0, 0, 0
        write(30,'(3I4)') 0, 0, 0
        go to 1000
        endif

        ikr8 = ikr
        jkr8 = jkr
        indk8(:) = indk(:)
        indl8(:) = indl(:)
        rklr8(:) = rklr(:)
        rkli8(:) = rkli(:)

        Do inz = 1,nz
        ! Debug output (if write(*,*))
        if (inz == 1) then
        write(*,*)"new_ikr2"
        write(*,*)"Filename:", Filename
        write(*,*)"inz:", inz
        endif
        iii = indmor(kr(ikr8))
        if (inz == 1) then
        write(*,*) "kr(ikr)", kr(ikr)
        write(*,*) "indmor(kr(ikr))", indmor(kr(ikr))
        endif

        jjj = indmor(kr(jkr8))
        if (inz == 1) then
        write(*,*) "kr(jkr)", kr(jkr)
        write(*,*) "indmor(kr(jkr))", indmor(kr(jkr))
        endif

        kkk = indmor(kr(indk8(inz)))
        if (inz == 1) then
        write(*,*) "indk(inz)", indk(inz)
        write(*,*) "kr(indk(inz))", kr(indk(inz))
        write(*,*) "indmor(kr(indk(inz)))", indmor(kr(indk(inz)))
        endif

        lll = indmor(kr(indl8(inz)))
        if (inz == 1) then
        write(*,*) "indl(inz)", indl(inz)
        write(*,*) "kr(indl(inz))", kr(indl(inz))
        write(*,*) "indmor(kr(indl(inz)))", indmor(kr(indl(inz)))
        endif

        iikr = (-1)**(mod(iii,2)+1)*(iii/2+mod(iii,2))
        jjkr = (-1)**(mod(jjj,2)+1)*(jjj/2+mod(jjj,2))
        kkkr = (-1)**(mod(kkk,2)+1)*(kkk/2+mod(kkk,2))
        llkr = (-1)**(mod(lll,2)+1)*(lll/2+mod(lll,2))

        iiit = iii-(-1)**iii
        jjjt = jjj-(-1)**jjj
        kkkt = kkk-(-1)**kkk
        lllt = lll-(-1)**lll


! Iwamuro debug
        if (inz == 1) then
        write(*,*) "new_ikr2", iikr, jjkr, kkkr, llkr
        endif
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

        If(iikr>0 .and. jjkr>0 .and. kkkr>0 .and. llkr>0) then  !TYPE1
        if( (ii<=jj .and. kk<=ll .and. (ii<kk .or. (ii==kk .and. jj<=ll))) .or. &
                (ii<=jj .and. ll<=kk .and. (ii<ll .or. (ii==ll .and. jj<=kk))) ) then
                if(abs(rklr(inz))>cutoff.or. &
                abs(rkli(inz))>cutoff      ) then
!                    write(28) -ikr,-jkr,nnz,-(indk(inz)),-(indl(inz)),rklr(inz),-(rkli(inz))
!                    write(28) -iikr,-jjkr,nnz,-kkkr,-llkr,rklr8(inz),-(rkli8(inz))
                write(28) iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
                write(29,'(5I4,2E32.16)') -iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
                write(30,'(5I4,2E32.16)') iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
                endif
        endif

        Else if (iikr>0 .and. jjkr<0 .and. kkkr>0 .and. llkr<0) then  !TYPE2
        if (ii<=jj .and. kk<=ll .and. (ii<kk .or. (ii==kk .and. jj<=ll))) then
                if(abs(rklr(inz))>cutoff.or. &
                abs(rkli(inz))>cutoff      ) then
!                    write(28) -ikr,-jkr,nnz,-(indk(inz)),-(indl(inz)),rklr(inz),-(rkli(inz))
!                    write(28) -iikr,-jjkr,nnz,-kkkr,-llkr,rklr8(inz),-(rkli8(inz))
                write(28) iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
                write(29,'(5I4,2E32.16)') -iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
                write(30,'(5I4,2E32.16)') iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
                endif
        endif

        Else if (iikr>0 .and. jjkr<0 .and. kkkr<0 .and. llkr>0) then  !TYPE3
        if(ii<=jj .and. kk<=ll .and. (ii<kk .or. (ii==kk .and. jj<=ll))) then
                if(abs(rklr(inz))>cutoff.or. &
                abs(rkli(inz))>cutoff      ) then
!                    write(28) -ikr,-jkr,nnz,-(indk(inz)),-(indl(inz)),rklr(inz),-(rkli(inz))
!                    write(28) -iikr,-jjkr,nnz,-kkkr,-llkr,rklr8(inz),-(rkli8(inz))
                write(28) iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
                write(29,'(5I4,2E32.16)') -iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
                write(30,'(5I4,2E32.16)') iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
                endif
        endif

        Else if (iikr>0 .and. jjkr<0 .and. kkkr<0 .and. llkr<0) then  !TYPE4
        if(ii<=jj) then
                if(abs(rklr(inz))>cutoff.or. &
                abs(rkli(inz))>cutoff      ) then
!                    write(28) -ikr,-jkr,nnz,-(indk(inz)),-(indl(inz)),rklr(inz),-(rkli(inz))
!                    write(28) -iikr,-jjkr,nnz,-kkkr,-llkr,rklr8(inz),-(rkli8(inz))
                write(28) iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
                write(29,'(5I4,2E32.16)') -iikr,-jjkr,nnz,-kkkr,-llkr,rklr(inz),-(rkli(inz))
                write(30,'(5I4,2E32.16)') iiit,jjjt,nnz,kkkt,lllt,rklr8(inz),-(rkli8(inz))
                endif
        endif
        Endif
300     Enddo

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

1000    close(mdcint)
        close(28)
        close(29)
        close(30)

        deallocate(kr)
        deallocate(indk)
        deallocate(indl)
        deallocate(rklr)
        deallocate(rkli)
        deallocate(kkr)
        deallocate(indk8)
        deallocate(indl8)
        deallocate(rklr8)
        deallocate(rkli8)
        end do
end Subroutine create_newmdcint
