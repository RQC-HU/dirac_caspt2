!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine create_newmdcint ! 2 Electorn Integrals In Mdcint

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        Use four_caspt2_module

        Implicit None

        Character*50 :: Filename

        Character  :: datex*10, timex*8

        integer(4) :: nkr, i0, mdcint, nz, inz,  !x, dammy
        integer(4) :: ikr, jkr, ii, jj, kk, ll
        integer(4) :: iikr, jjkr, kkkr, llkr, iii, jjj, kkk, lll
        integer(4), allocatable :: indk(:), indl(:), kr(:)
        double precision, allocatable  :: rklr(:), rkli(:)

! Iwamuro modify
        logical(4) :: realonly
        real    :: cutoff

        realonly = .false.
        cutoff = 0.25D-12

        Allocate(kr(-nmo/2:nmo/2))
        Allocate(indk(nmo**2))
        Allocate(indl(nmo**2))
        Allocate(rklr(nmo**2))
        Allocate(rkli(nmo**2))

        kr = 0

        mdcint=11
        open(mdcint, file="MDCINT_original", form ='unformatted', status='unknown')

        read (mdcint) datex,timex,nkr,(kr(i0),kr(-1*i0),i0=1,nkr)

        read (mdcint,ERR=200) ikr,jkr, nz, (indk(inz),indl(inz), rklr(inz),rkli(inz), inz=1,nz)

        goto 201

200     realonly = .true.
        write(*,*) "realonly = ", realonly
201     close(mdcint)

        open(mdcint, file="MDCINT_original", form='unformatted', status='unknown')
        open(28, file="MDCINT", form='unformatted', status='unknown')
        open(29, file="MDCINT_debug", form='formatted', status='unknown')
        open(30,file="MDCINT_int", form='formatted', status='unknown')

        read (mdcint) datex,timex,nkr, (kr(i0),kr(-1*i0),i0=1,nkr)
        write(28) datex, timex, nkr, (kr(i0),kr(-1*i0),i0=1,nkr)

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

!        write(*,*) 'ikr,jkr, nz', ikr,jkr, nz
!        write(*,*) 'indk(inz),indl(inz)', (indk(inz),indl(inz), inz=1,nz)
!        write(*,*) 'rklr(inz),rkli(inz)', (rklr(inz),rkli(inz), inz=1,nz)

!------------------------------------------------------------

!------------------------------!
!   Create new ikr for UTChem  !
!------------------------------!

!       new ikr = iikr
!           jkr = jjkr
!           kkr = kkkr
!           lkr = llkr

        iii = indmor(kr(-ikr))
        jjj = indmor(kr(-jkr))
        kkk = indmor(kr(-(indk(inz))))
        lll = indmor(kr(-(indl(inz))))

        iikr = (-1)**(mod(iii,2)+1)*(iii/2+mod(iii,2))
        jjkr = (-1)**(mod(jjj,2)+1)*(jjj/2+mod(jjj,2))
        kkkr = (-1)**(mod(kkk,2)+1)*(kkk/2+mod(kkk,2))
        llkr = (-1)**(mod(lll,2)+1)*(lll/2+mod(lll,2))

! Iwamuro debug
!        write(*,*) 'iikr, jjkr, kkkr, llkr', iikr, jjkr, kkkr, llkr

!------------------------------------------------------------

        if (iikr<0) go to 100
        if (iikr == 0) then
          write(28) 0, 0, 0
          go to 1000
        endif

        Do inz = 1,nz
           ii = abs(iikr)
           jj = abs(jjkr)
           kk = abs(kkkr)
           ll = abs(llkr))

       !---------------------------
       ! TYPE1 (++++) = (ij|kl)
       ! TYPE2 (+-+-) = (ij~|kl~)
       ! TYPE3 (+--+) = (ij~|k~l)
       ! TYPE4 (+---) = (ij~|k~l~)
       !---------------------------

        If(iikr>0 .and. jjkr>0 .and. kkkr>0 .and. llkr>0) then  !TYPE1
              if( (ii<=jj .and. kk<=ll .and. (ii<kk .or. (ii==kk .and. jj<=ll))) .or. &
                   (ii<=jj .and. ll<=kk .and. (ii<ll .or. (ii==ll .and. jj<=kk))) ) then
!                  if(ii==2 .and. jj==2 .and. kk==4 .and. ll==3) go to 300
                  if(abs(rklr(inz))>cutoff.or. &
                     abs(rkli(inz))>cutoff      ) then
                    write(28) -ikr,-jkr,1,-(indk(inz)),-(indl(inz)),rklr(inz),rkli(inz)
                    write(29,'(4I4,2E32.16)') iikr,jjkr,kkkr,llkr,rklr(inz),rkli(inz)
                    write(30,'(4I4,2E32.16)') -ikr,-jkr,-(indk(inz)),-(indl(inz)),rklr(inz),rkli(inz)
                  endif
!                  endif
              endif

        Else if (iikr>0 .and. jjkr<0 .and. kkkr>0 .and. llkr<0) then  !TYPE2
              if (ii<=jj .and. kk<=ll .and. (ii<kk .or. (ii==kk .and. jj<=ll))) then
                  if(abs(rklr(inz))>cutoff.or. &
                     abs(rkli(inz))>cutoff      ) then
                    write(28) -ikr,-jkr,1,-(indk(inz)),-(indl(inz)),rklr(inz),rkli(inz)
                    write(29,'(4I4,2E32.16)') iikr,jjkr,kkkr,llkr,rklr(inz),rkli(inz)
                    write(30,'(4I4,2E32.16)') -ikr,-jkr,-(indk(inz)),-(indl(inz)),rklr(inz),rkli(inz)
                  endif
              endif

        Else if (iikr>0 .and. jjkr<0 .and. kkkr<0 .and. llkr>0) then  !TYPE3
              if(ii<=jj .and. kk<=ll .and. (ii<kk .or. (ii==kk .and. jj<=ll))) then
                  if(abs(rklr(inz))>cutoff.or. &
                     abs(rkli(inz))>cutoff      ) then
                    write(28) -ikr,-jkr,1,-(indk(inz)),-(indl(inz)),rklr(inz),rkli(inz)
                    write(29,'(4I4,2E32.16)') iikr,jjkr,kkkr,llkr,rklr(inz),rkli(inz)
                    write(30,'(4I4,2E32.16)') -ikr,-jkr,-(indk(inz)),-(indl(inz)),rklr(inz),rkli(inz)
                  endif
              endif

        Else if (iikr>0 .and. jjkr<0 .and. kkkr<0 .and. llkr<0) then  !TYPE4
              if(ii<=jj) then
                  if(abs(rklr(inz))>cutoff.or. &
                     abs(rkli(inz))>cutoff      ) then
                    write(28) -ikr,-jkr,1,-(indk(inz)),-(indl(inz)),rklr(inz),rkli(inz)
                    write(29,'(4I4,2E32.16)') iikr,jjkr,kkkr,llkr,rklr(inz),rkli(inz)
                    write(30,'(4I4,2E32.16)') -ikr,-jkr,-(indk(inz)),-(indl(inz)),rklr(inz),rkli(inz)
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

end Subroutine create_newmdcint
