!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Subroutine create_newmdcint ! 2 Electorn Integrals In Mdcint

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        Use four_caspt2_module

        Implicit None

        Character*50 :: Filename

        Character  :: datex*10, timex*8

        integer    :: nkr, nz
        integer    :: i0, mdcint, inz, nnz , n
        integer    :: ikr, jkr
        integer    :: ii, jj, kk, ll
        integer    :: iikr, jjkr, kkkr, llkr, iii, jjj, kkk, lll
        integer, allocatable :: indk(:), indl(:), kr(:)
        double precision, allocatable  :: rklr(:), rkli(:)

! Iwamuro modify
        real    :: cutoff
        integer    :: nnkr, ikr8, jkr8, iiit, jjjt, kkkt, lllt
        integer, allocatable :: kkr(:), indk8(:), indl8(:) 
        real*8, allocatable  :: rklr8(:), rkli8(:)

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
 
        mdcint=11
        open(mdcint, file="MDCINT", form ='unformatted', status='unknown')

        read (mdcint) datex,timex,nkr,(kr(i0),kr(-1*i0),i0=1,nkr)

        read (mdcint,ERR=200) ikr,jkr, nz, (indk(inz),indl(inz), rklr(inz),rkli(inz), inz=1,nz)

        goto 201

200     realonly = .true.
        write(*,*) "realonly = ", realonly
201     close(mdcint)

!Iwamuro modify
        write(*,*) "nz=" ,nz

        open(mdcint, file="MDCINT", form='unformatted', status='unknown')
        open(28, file="MDCINTNEW", form='unformatted', status='unknown')
        open(29, file="MDCINT_debug", form='formatted', status='unknown')
        open(30, file="MDCINT_int", form='formatted', status='unknown')

        read (mdcint) datex,timex,nkr, (kr(i0),kr(-1*i0),i0=1,nkr)

        nnkr = nkr
        kkr(:) = kr(:)

        write(28) datex, timex, nnkr, (kkr(i0),kkr(-1*i0),i0=1,nnkr)
!        write(29,*) datex, timex, nnkr, (kkr(i0),kkr(-1*i0),i0=1,nnkr)
!Iwamuro debug
!        write(*,*) "new_ikr1", datex, timex, nkr, (kr(i0),kr(-1*i0),i0=1,nkr)   

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

!        write(*,*) ""  
!        write(*,*) ikr,jkr, nz, &
!                    (indk(inz),indl(inz),inz=1,nz), &
!                    (rklr(inz),rkli(inz), inz=1,nz)

!        write(*,*) ""
   
!------------------------------------------------------------
 
!------------------------------!
!   Create new ikr for UTChem  !
!------------------------------!

!       new ikr = iikr 
!           jkr = jjkr
!           kkr = kkkr
!           lkr = llkr

!        if (ikr<0) go to 100 ! 2021.10.18 Abe commented out
!
        if (ikr == 0) then
          write(28) 0, 0, 0
          write(29,'(3I4)') 0, 0, 0
          write(30,'(3I4)') 0, 0, 0
          go to 1000
        endif

!Abe        ikr8 = ikr
!Abe         jkr8 = jkr
!Abe         indk8(:) = indk(:)
!Abe         indl8(:) = indl(:)
!Abe         rklr8(:) = rklr(:)
!Abe         rkli8(:) = rkli(:)

! Abe 2021.10.18
        n = 1
 10          select case (n)
        case (1)
           ikr8 = ikr
           jkr8 = jkr
           indk8(:) = indk(:)
           indl8(:) = indl(:)
           rklr8(:) = rklr(:)
           rkli8(:) = rkli(:)
        case (2)
           ikr8 = -ikr
           jkr8 = -jkr
           indk8(:) = -indk(:)
           indl8(:) = -indl(:)
           rklr8(:) = -rklr(:)* sign(1,ikr*jkr*indk(:)*indl(:))
           rkli8(:) =  rkli(:)* sign(1,ikr*jkr*indk(:)*indl(:))
        end select
! Abe 2021.10.18

        Do inz = 1,nz

!        write(*,*)"new_ikr2"

        iii = indmor(kr(ikr8))
!        write(*,*) "kr(ikr)", kr(ikr)
!        write(*,*) "indmor(kr(ikr))", indmor(kr(ikr))

        jjj = indmor(kr(jkr8))
!        write(*,*) "kr(jkr)", kr(jkr)
!        write(*,*) "indmor(kr(jkr))", indmor(kr(jkr))

        kkk = indmor(kr(indk8(inz)))
!        write(*,*) "indk(inz)", indk(inz)
!        write(*,*) "kr(indk(inz))", kr(indk(inz))
!        write(*,*) "indmor(kr(indk(inz)))", indmor(kr(indk(inz)))

        lll = indmor(kr(indl8(inz)))
!        write(*,*) "indl(inz)", indl(inz)
!        write(*,*) "kr(indl(inz))", kr(indl(inz))
!        write(*,*) "indmor(kr(indl(inz)))", indmor(kr(indl(inz)))

        iikr = (-1)**(mod(iii,2)+1)*(iii/2+mod(iii,2))
        jjkr = (-1)**(mod(jjj,2)+1)*(jjj/2+mod(jjj,2))              
        kkkr = (-1)**(mod(kkk,2)+1)*(kkk/2+mod(kkk,2))
        llkr = (-1)**(mod(lll,2)+1)*(lll/2+mod(lll,2))

        iiit = iii-(-1)**iii
        jjjt = jjj-(-1)**jjj
        kkkt = kkk-(-1)**kkk
        lllt = lll-(-1)**lll


! Iwamuro debug
!        write(*,*) "new_ikr2", iikr, jjkr, kkkr, llkr
 
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

! Abe 2021.10.18
        If(n==1) then
           n=2
           go to 10
        else
           go to 100
        endif
!        go to 100
! Abe 2021.10.18

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

end Subroutine create_newmdcint
 