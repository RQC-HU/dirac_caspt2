! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE fockivo_ty(nhomo) ! TO MAKE FOCK MATRIX for IVO

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE

        integer, intent(in) :: nhomo

        integer :: ii, jj, kk, ll
        integer :: j, i, k, l, i0, j0
        integer :: nint, n, nsym, isym, nv, numh

        real*8 :: i2r, i2i, dr, di, nsign, thresd
        complex*16 :: cmplxint, dens
        logical   ::cutoff

        complex*16,allocatable :: fsym(:,:), fdmmy(:,:)
        complex*16,allocatable :: coeff(:,:,:),readmo(:,:,:)
        real*8,    allocatable :: wsym(:)
        integer,   allocatable :: mosym(:)


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
              enddo
           enddo

           allocate(readmo   (nbas*2, nbas*2, 2))
           allocate(itrfmo   (nbas*2, nbas,   2))
           itrfmo = 0.0d+00

           open(15,file='r4dorbcoeff',status='old',form='unformatted')
           read(15,err=10) readmo
           close(15)

           itrfmo(1:nbas*2,1:nbas,1:2) = readmo(1:nbas*2, nbas+1:nbas*2,1:2)

           Do isym = 1, nsymrpa, 2
              nv = 0
              Do i = 1, nsec
                 i0 = i+ninact+nact
                 if(irpmo(i0)==isym) then
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
                    nv = nv + 1
                    mosym(nv) = i
                 endif
              enddo
              ! Noda 2021/12/27 max(nv) = nsec. So the max dimention of fsym is nsec (fsym(nsec,nsec))
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

              Allocate(coeff(nbas*2,nv,2))
              ! Noda 2021/12/27 max(nv) = nsec. max : (coeff(nbas*2,nsec,2))
              Do i = 1, nv
                 i0 = mosym(i)+ncore+ninact+nact
                 coeff(:,i,:)=itrfmo(:,i0,:)
              Enddo

              coeff(:,:,1) = MATMUL(coeff(:,:,1),fsym(:,:))
              coeff(:,:,2) = MATMUL(coeff(:,:,2),fsym(:,:))

              Do i = 1, nv
                 i0 = mosym(i)+ncore+ninact+nact
                 itrfmo(:,i0,:) = coeff(:,i,:)
              Enddo

! Kramers - pairs

              Do i = 1, nv
                 i0 = mosym(i)+ncore+ninact+nact+1
                 coeff(:,i,:)=itrfmo(:,i0,:)
              Enddo

              coeff(:,:,1) = MATMUL(coeff(:,:,1),DCONJG(fsym(:,:)))
              coeff(:,:,2) = MATMUL(coeff(:,:,2),DCONJG(fsym(:,:)))

              Do i = 1, nv
                 i0 = mosym(i)+ncore+ninact+nact+1
                 itrfmo(:,i0,:) = coeff(:,i,:)
              Enddo

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
          enddo

          readmo(1:nbas*2, nbas+1:nbas*2,1:2) = itrfmo(1:nbas*2,1:nbas,1:2)

          open(15,file='r4dorbcoeff_ivo',status='unknown',form='unformatted')
          write(15) readmo
          close(15)
          goto 100

 10       write(*,*)'reading err of r4dorbcoeff'
!        deallocate(fdmmy)
          deallocate(readmo)
          deallocate(itrfmo)

 100      write(*,*)'fockivo_ty end'
          end
