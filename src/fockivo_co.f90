! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE fockivo_co ! TO MAKE FOCK MATRIX for IVO

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE

 
        integer :: ii, jj, kk, ll, m
        integer :: j, i, k, l, i0, j0
        integer :: nint, n, nsym, isym, nv, numh
        integer :: imo, iao
        real*8 :: i2r, i2i, dr, di, nsign, thresd
        complex*16 :: cmplxint, dens
        logical   ::cutoff
        integer :: iostat
        complex*16,allocatable :: fsym(:,:), fdmmy(:,:)
        complex*16,allocatable :: coeff(:,:)
        real*8,    allocatable :: wsym(:), BUF(:), BUF1(:,:), eval(:)
        integer,   allocatable :: mosym(:)
        integer :: nnmo, nnao , IMAX
        character*150 :: line0, line1, line2, line3, line4, line5

      ! for new code of IVO
        integer :: npg, neg, npu, neu, nbasg, nbasu, nsum, A, nv0, ngu, B
        integer, allocatable :: syminfo(:), dmosym(:)
        complex*16,allocatable :: itrfmog(:,:), itrfmou(:,:)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=


!! NOW MAKE FOCK MATRIX FOR IVO
!! fij = hij + SIGUMA_k (ij|kk)-(ik|kj)} i, j run over virtual spinors k runs occupied spinors except HOMO

        f = 0.0d+00

        write(*,*) 'enter building fock matrix for IVO'


        if(nhomo==0) then
           numh = 0
           do i = 1, ninact+nact
              if(ABS(orbmo(i)-orbmo(nelec+ninact)) < 1.0d-01) then
                 numh = numh + 1
              endif
           enddo
        else
           numh = nhomo
        endif

        write(*,*)'number of degeneracy of HOMO is',numh, DBLE(numh),1.0d+00/DBLE(numh)

        do i = 1, nsec
           i0 = i + ninact  + nact
           f(i,i) = orbmo(i0)
           do j = i, nsec
              j0 = j + ninact  + nact
              do k = ninact+nact-numh+1, ninact+nact

                 if(k > ninact+nact-2 .and. mod(nelec,2)==1) then

                    f(i,j) = f(i,j) &
                        &  - 0.5d+00*DCMPLX(int2r_f1(i0,j0,k,k),int2i_f1(i0,j0,k,k))/DBLE(numh)
                    f(i,j) = f(i,j) &
                        &  + 0.5d+00*DCMPLX(int2r_f2(i0,k,k,j0),int2i_f2(i0,k,k,j0))/DBLE(numh)

                 else
                    f(i,j) = f(i,j) - DCMPLX(int2r_f1(i0,j0,k,k),int2i_f1(i0,j0,k,k))/DBLE(numh)
                    f(i,j) = f(i,j) + DCMPLX(int2r_f2(i0,k,k,j0),int2i_f2(i0,k,k,j0))/DBLE(numh)
                 endif

              enddo

           enddo
        enddo


           do i = 1, nsec
              do j = i, nsec
                f(j,i) = DCONJG(f(i,j))
                
!                write(*,*)"f(j,i)", f(j,i)
              enddo
           enddo

!           allocate(readmo   (nbas*2, nbas*2, 2))
!           allocate(itrfmo   (nbas*2, nbas,   2))
!           itrfmo = 0.0d+00

!           allocate(readmo (nbas, nbas))
!           allocate(itrfmo (nbas, nbas, 2))
!           allocate(mocr   (nbas, nbas))
!           allocate(moci   (nbas, nbas))

!           open(15,file='r4dorbcoeff',status='old',form='unformatted')
!           read(15,err=10) readmo

!           Allocate (BUF(nbas*lscom))
!           Allocate (BUF1(nbas,lscom*2))

           IMAX = nbas*lscom

           open(15,file='DFPCMO',status='old',form='formatted')

! From DIRAC dirgp.F WRIPCMO (Write DHF-coefficients and eigenvalues )

           if(dirac_version == 21) then
              read(15,'(A150)') line0
           elseif(dirac_version == 22) then
              read(15,'(A150)') line0
           endif
           read(15,'(A150)') line1
           if(dirac_version == 21) then
            read(15,*) A, B, npg, neg, nbasg, npu, neu, nbasu
            write(*,*) A, B, npg, neg, nbasg, npu, neu, nbasu
           elseif(dirac_version == 22) then
            read(15,*) A, B, npg, neg, nbasg, npu, neu, nbasu
            write(*,*) A, B, npg, neg, nbasg, npu, neu, nbasu
           else
            read(15,*) A, npg, neg, nbasg, npu, neu, nbasu
            write(*,*) A, npg, neg, nbasg, npu, neu, nbasu
           endif
           read(15,'(A150)') line2

           write(*,*)'end reading information, symmetry information and energy'

           nsum = npg + neg + npu + neu
           
           ! ngu : the sum of gerade and ungerade
           ngu = (npg+neg)*nbasg+(npu+neu)*nbasu

           allocate(itrfmog(nbasg, neg-noccg))
           allocate(itrfmou(nbasu, neu-noccu))
           allocate(eval(nsum))
           allocate(syminfo(nsum))
           Allocate(BUF(ngu))

           itrfmog = 0.0d+00
           itrfmou = 0.0d+00

        if(dirac_version == 21) then
            read(15,'(A150)') line3
        elseif(dirac_version == 22) then
            read(15,'(A150)') line3
        endif

          ! Read MO coefficient of DFPCMO  
             Do I = 1, ngu, 6
              Read(15,'(6F22.16)',ERR=22) BUF(I:I+5)
             Enddo 

             write(*,*)'end reading MO coefficient'

           ! unoccupid, gerade, electron 
             Do iao = 1, nbasg
               DO imo = 1, neg-noccg
                 itrfmog(iao,imo) = BUF((npg+noccg+(imo-1))*nbasg+iao)
                Enddo
             Enddo

           open(37,file='itrfmog_before',status='unknown',form='formatted')

           Do iao = 1,nbasg
                write(37,*) (itrfmog(iao,imo), imo = 1, neg-noccg)
             Enddo

           close(37)

           ! unoccupid, ungerade, electron
             Do iao = 1, nbasu
               DO imo = 1, neu-noccu
                 itrfmou(iao,imo) = BUF((npg+neg)*nbasg + (npu+noccu+(imo-1))*nbasu+iao)
               Enddo
             Enddo

             open(38,file='itrfmou_before',status='unknown',form='formatted')

             Do iao = 1,nbasu
                  write(38,*) (itrfmou(iao,imo), imo = 1, neu-noccu)
               Enddo
  
             close(38)


 22         continue

!            if(mod(nsum,6).eq.0) then
!               Do I = 1, nsum/6*6, 6
!                 Read(15,'(6E22.12)',ERR=23) eval(I:I+5)
!               End do
!            elseif(mod(nsum,6).ne.0) then
!               I = nsum/6*6 +1
!               Read(15,'(6E22.12)',ERR=23) eval(I:IMAX)
!            endif
            if(dirac_version == 21) then
                read(15,'(A150)') line4
            elseif(dirac_version == 22) then
                read(15,'(A150)') line4
            endif

             Read(15,*) eval

             Do i = 1, nsum
                write(*,*) "eval(i)", eval(i)
              Enddo

            write(*,*)'end reading eigenvalue'

 23        continue

            ! Read syminfo from DFPCMO
            if(dirac_version == 21) then
                read(15,'(A150)') line5
            elseif(dirac_version == 22) then
                read(15,'(A150)') line5
            endif

             Read(15,*) syminfo

            Do i = 1, nsum
              write(*,*) "syminfo(i)", syminfo(i)
            Enddo

            write(*,*)'end reading symmetry information2'

!           read(15,'(A150)',iostat=iostat) line4
!           if (iostat < 0) goto  120
!           read(15,'(A150)',iostat=iostat,ERR=24) line5
!           if (iostat < 0) goto  120
!           read(15,'(A150)',iostat=iostat,ERR=24) line6
!           if (iostat < 0) goto  120
!           read(15,'(A150)',iostat=iostat,ERR=24) line7
!           if (iostat < 0) goto  120
!           read(15,'(A150)',iostat=iostat,ERR=24) line8
!           if (iostat < 0) goto  120
!           read(15,'(A150)',iostat=iostat,ERR=24) line9
!           if (iostat < 0) goto  120
!           read(15,'(A150)',iostat=iostat,ERR=24) line10
!           if (iostat < 0) goto  120
!           read(15,'(A150)',iostat=iostat,ERR=24) line11
!           if (iostat < 0) goto  120
!           read(15,'(A150)',iostat=iostat,ERR=24) line12
!           if (iostat < 0) goto  120
!           read(15,'(A150)',iostat=iostat,ERR=24) line13
!           if (iostat < 0) goto  120
!           read(15,'(A150)',iostat=iostat,ERR=24) line14
!           if (iostat < 0) goto  120
!           read(15,'(A150)',iostat=iostat,ERR=24) line15

120        continue
24         continue

           close(15)

           open(105,file='BUF_write',status='unknown',form='formatted')

           Do I = 1,ngu, 6
              Write(105,'(6F22.16)') BUF(I:I+5)
           Enddo

           close(105)

!           open(16,file='DFPCMO_complement',status='unknown',form='formatted')

! DFPCMO complement in the symmetry of C32h

!            I = 0

! Electron solution / gerade    

!           Do I = 1,nbas/4
!             Do J = 1,lscom
!               BUF1(I,J) = BUF(J+lscom*(I-1+nbas/4))
!               BUF1(I,J+lscom) = BUF1(I,J)
!             Enddo
!             Do J = 1,lscom*2
!               BUF1(I+nbas/4,J) = BUF1(I,J)
!             Enddo
!           Enddo

! Electron solution / ungerade

!           Do I = 1+nbas/2, nbas/2+nbas/4
!             Do J = 1,lscom
!               BUF1(I,J) = BUF(J+lscom*(I-1+nbas/4))
!               BUF1(I,J+lscom) = -BUF1(I,J)
!             Enddo
!             Do J = 1,lscom*2
!               BUF1(I+nbas/4,J) = BUF1(I,J)
!             Enddo
!           Enddo

! Extract only electron solution 
!             Do I = 1,nbas
!               Write (16,'(24F10.5)') (BUF1(I,J),J=1,lscom*2)
!             Enddo

!           close(16)

!          Do iao = 1,lscom*2
!             Do imo = 1,nbas
!                mocr(iao,imo) = BUF1(imo,iao)
!             Enddo
!          Enddo

!          moci=0.0d0+00

!          open(87,file='mocr',status='unknown',form='formatted')
!          Do iao = 1,lscom*2
!             Do imo = 1,nbas
!                write(87,*) mocr(iao,imo)
!             Enddo
!          Enddo
!          close(87)

! Transform MO order into orbital energy order

!          do imo = 1, nbas
!             readmo(:,indmor(imo),2) = mocr(:,imo)
!          enddo

!          open(67,file='readmo',status='unknown',form='formatted')

!          Do iao = 1,lscom*2
!             Do imo = 1,nbas
!                write(67,'(24F10.5)') (real(readmo(iao,imo,2)),iao = 1,lscom*2)
!             Enddo
!          Enddo

!           close(67)

!          itrfmog(1:lscom*2,1:nbas) = readmo(1:lscom*2,1:nbas)

!           itrfmo(1:lscom*2,1:nbas,1:2) = readmo(1:lscom*2, nbas+1:nbas*2,1:2)

! IVO calculation 

! gerade   

           Do isym = 1, nsymrpa, 2
              nv = 0
              Do i = 1, nsec
                 i0 = i+ninact+nact
                 if(irpmo(i0)==isym) then
!                 if(irpamo(i0)==isym) then
                    nv = nv + 1
                endif
              enddo

              Allocate(mosym(nv))
              Allocate(fsym(nv,nv))

              fsym = 0.0d+00
              nv = 0
              Do i = 1, nsec
                 i0 = i+ninact+nact
                 if(irpmo(i0)==isym) then
!                 if(irpamo(i0)==isym) then
                    nv = nv + 1
                    mosym(nv) = i
                 endif
              enddo

!C32h gerade
            if (isym <= nsymrpa/2) then
             nv0 = 0
             Do i0 = npg+noccg+1,npg+neg
                if(ABS(syminfo(i0))==isym) then
!                 if(irpamo(i0)==isym) then
                   nv0 = nv0 + 1
                endif
             enddo

             Allocate(dmosym(nv0))

             nv0 = 0
             Do i0 = npg+noccg+1,npg+neg
                if(ABS(syminfo(i0))==isym) then
!                 if(irpamo(i0)==isym) then
                   nv0 = nv0 + 1
                   dmosym(nv0) = i0
                endif
             enddo

!C32h ungerade
            else
             nv0 = 0
             Do i0 = npg+neg+npu+noccu+1,npg+neg+npu+neu
                if(ABS(syminfo(i0))+ nsymrpa/2 ==isym) then
!                 if(irpamo(i0)==isym) then
                  nv0 = nv0 + 1
                endif
             enddo

             Allocate(dmosym(nv0))

             nv0 = 0
             Do i0 = npg+neg+npu+noccu+1,npg+neg+npu+neu
               if(ABS(syminfo(i0))+ nsymrpa/2 ==isym) then
!                 if(irpamo(i0)==isym) then
                  nv0 = nv0 + 1
                  dmosym(nv0) = i0
               endif
             enddo
            endif


              Do i = 1, nv
                 i0 = mosym(i)
                 Do j = i, nv
                    j0 = mosym(j)
                    fsym(i,j) = f(i0, j0)
                    fsym(j,i) = DCONJG(f(i0, j0))
!                    write(*,*)fsym(i,j)
                 enddo
              enddo

              Allocate(wsym(nv))
              wsym = 0.0d+00
              cutoff = .FALSE.
              thresd = 0.0d+00

              call cdiag (fsym, nv, nv,wsym , thresd, cutoff)

            ! Gerade  
              if (isym <= nsymrpa/2) then
                Allocate(coeff(nbasg,nv))  
                Do i = 1, nv0
                 i0 = dmosym(i)-npg-noccg
                 coeff(:,i)=itrfmog(:,i0)
                 write(*,*)
                Enddo

            ! Ungerade    
              else
                Allocate(coeff(nbasu,nv))  
                Do i = 1, nv0
                i0 = dmosym(i)-npg-neg-npu-noccu
                coeff(:,i)=itrfmou(:,i0)
                Enddo
              endif

              coeff(:,:) = MATMUL(coeff(:,:),fsym(:,:))
            
            ! Gerade  
              if (isym <= nsymrpa/2) then              
                Do i = 1, nv0
                i0 = dmosym(i)-npg-noccg
                itrfmog(:,i0) = coeff(:,i)
                Enddo

            ! Ungerade 
            else
                Do i = 1, nv0
                i0 = dmosym(i)-npg-neg-npu-noccu
                itrfmou(:,i0) = coeff(:,i)
                Enddo
              endif  

! Kramers - pairs
              
!              Do i = 1, nv
!                 i0 = mosym(i)+ncore+ninact+nact+1
!                 coeff(:,i)=itrfmo(:,i0)
!              Enddo

!              coeff(:,:) = MATMUL(coeff(:,:),DCONJG(fsym(:,:)))
              
!              Do i = 1, nv
!                 i0 = mosym(i)+ncore+ninact+nact+1
!                 itrfmo(:,i0,:) = coeff(:,i,:)
!              Enddo

! Kramers - pairs

              deallocate(coeff)

              Do i = 1, nv
                 i0 = mosym(i)
                 write(*,'(I4,F20.10)')i0,wsym(i)
              enddo

              Do i = 1, nv
                 i0 = mosym(i)
                 write(*,*)''
                 write(*,*)'new ',i0+ninact+nact,'th ms consists of '
                 Do j = 1, nv
                    j0 = mosym(j)
                    if(ABS(fsym(j,i))**2>1.0d-03) then
                       write(*,'(I4,"  Weights ",F20.10)')j0+ninact+nact,ABS(fsym(j,i))**2
                    endif
                 enddo
             enddo
             deallocate(fsym)
             deallocate(wsym)
             deallocate(mosym)
             deallocate(dmosym)
          enddo

!           open(37,file='itrfmo',status='unknown',form='formatted')

!           Do iao = 1,lscom*2
!             Do imo = 1,nbas
!                write(37,*) itrfmo(iao,imo,2)
!             Enddo
!           Enddo

!           close(37)

!          readmo(1:lscom*2, 1:nbas,2) = itrfmo(1:lscom*2,1:nbas,2)

!          open(15,file='r4dorbcoeff_ivo',status='unknown',form='unformatted')
!          write(15) readmo
!          close(15)

!          open(57,file='readmo_ivo',status='unknown',form='formatted')

!          Do iao = 1,lscom*2
!             Do imo = 1,nbas
!                write(57,'(24F10.5)') (real(readmo(iao,imo,2)), iao=1,lscom*2)
!             Enddo
!          Enddo

!           close(57)

! Transform MO order into irreducible representation order

!          Do imo = 1, nbas
!            Do iao = 1,lscom*2
!                 readmo1(indmo(imo),iao,2)=readmo(iao,imo,2)
!            Enddo
!          Enddo

! Extract only the necessary part of the electronic solution and return it to BUF.

!          open(77,file='readmo1',status='unknown',form='formatted')

!          Do iao = 1,lscom*2
!             Do imo = 1,nbas
!                write(77,'(24F10.5)') (real(readmo1(imo,iao,2)), iao=1,lscom*2)
!             Enddo
!          Enddo

!           close(77)

           ! unoccupid, gerade, electron 
          Do iao = 1, nbasg
            DO imo = 1, neg-noccg
             BUF((npg+noccg+(imo-1))*nbasg+iao) = itrfmog(iao,imo)
            Enddo
          Enddo

        ! unoccupid, ungerade, electron
          Do iao = 1, nbasu
            DO imo = 1, neu-noccu
              BUF((npg+neg)*nbasg + (npu+noccu+(imo-1))*nbasu+iao) = itrfmou(iao,imo)
            Enddo
          Enddo

!         I = 0

!         Do imo = 1, nbas/4
!           Do iao = 1, lscom
!                 I = I+1
!                 BUF(I+lscom*(imo+nbas/4-1))=real(readmo1(imo,iao,2))
!            Enddo
!          Enddo

!         I = 0

!         Do imo = 1+nbas/2, nbas/2+nbas/4
!           Do iao = 1, lscom
!                 I = I+1
!                 BUF(I+lscom*(imo+nbas/4-1))=real(readmo1(imo,iao,2))
!            Enddo
!          Enddo

! Create new DFPCMO : DFPCMONEW

          open(19,file='DFPCMONEW',status='unknown',form='formatted')
          if(dirac_version == 21) then
            write(19,'(A150)') line0
          elseif(dirac_version == 22) then
            write(19,'(A150)') line0
          endif
          write(19,'(A150)') line1
          if(dirac_version == 21) then
            write(19,'(8(X,I0))') A, B, npg, neg, nbasg, npu, neu, nbasu
          elseif(dirac_version == 22) then
            write(19,'(8(X,I0))') A, B, npg, neg, nbasg, npu, neu, nbasu
          else
            write(19,'(7(X,I0))') A, npg, neg, nbasg, npu, neu, nbasu
          endif
          write(19,'(A150)') line2
          if(dirac_version == 21) then
            write(19,'(A150)') line3
          elseif(dirac_version == 22) then
            write(19,'(A150)') line3
          endif
             Do I = 1, ngu, 6
               Write(19,'(6F22.16)') BUF(I:I+5)
             End do
          if(dirac_version == 21) then
            write(19,'(A150)') line4
          elseif(dirac_version == 22) then
                write(19,'(A150)') line4
          endif

          write(19,'(6E22.12)') eval
          if(dirac_version == 21) then
            write(19,'(A150)') line5
          elseif(dirac_version == 22) then
            write(19,'(A150)') line5
          endif
          write(19,'(66(X,I0))') (syminfo(i), i=1,nsum)

!          write(19,'(A150)') line4
!          write(19,'(A150)', ERR=25) line5
!          write(19,'(A150)', ERR=25) line6
!          write(19,'(A150)', ERR=25) line7
!          write(19,'(A150)', ERR=25) line8
!          write(19,'(A150)', ERR=25) line9
!          write(19,'(A150)', ERR=25) line10
!          write(19,'(A150)', ERR=25) line11
!          write(19,'(A150)', ERR=25) line12
!          write(19,'(A150)', ERR=25) line13
!          write(19,'(A150)', ERR=25) line14
!          write(19,'(A150)', ERR=25) line15

 25       continue

          close(19)

          goto 100

 10       write(*,*)'reading err of DFPCMO'
!        deallocate(fdmmy)
          deallocate(itrfmog, itrfmou)
          deallocate(BUF)
          deallocate(eval, syminfo)

 100      write(*,*)'fockivo_co end'
      end subroutine fockivo_co



