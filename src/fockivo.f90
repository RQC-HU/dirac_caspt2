! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE fockivo(nhomo) ! TO MAKE FOCK MATRIX for IVO

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

       complex*16, allocatable :: fsym(:, :), fdmmy(:, :)
       complex*16, allocatable :: coeff(:, :, :), readmo(:, :, :)
       real*8, allocatable :: wsym(:)
       integer, allocatable :: mosym(:)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!! NOW MAKE FOCK MATRIX FOR IVO
!! fij = hij + SIGUMA_k (ij|kk)-(ik|kj)} i, j run over virtual spinors k runs occupied spinors except HOMO

       f = 0.0d+00

       write (*, *) 'enter building fock matrix for IVO'

!        Allocate(fdmmy(nsec,nsec))
!        fdmmy = 0.0d+00
!        do i = 1, nsec
!           i0 = i + ninact  + nact
!           fdmmy(i,i) = orbmo(i0)
!           do j = i, nsec
!              j0 = j + ninact  + nact
!
!              do k = ninact+nelec-1, ninact+nelec
!                 Call intmo2(i0,j0,k,k,cmplxint)
!                 fdmmy(i,j) = fdmmy(i,j) - 0.5d+00*cmplxint
!                 Call intmo2(i0,k,k,j0,cmplxint)
!                 fdmmy(i,j) = fdmmy(i,j) + 0.5d+00*cmplxint
!              enddo
!
!           enddo
!        enddo

       if (nhomo == 0) then
           numh = 0
           do i = 1, ninact + nact
               if (ABS(orbmo(i) - orbmo(nelec + ninact)) < 1.0d-01) then
                   numh = numh + 1
               end if
           end do
       else
           numh = nhomo
       end if

!        if(mod(nelec,2)==0) then
!           numh = numh
!        else
!           numh = numh-1
!        endif

       write (*, *) 'number of degeneracy of HOMO is', numh, DBLE(numh), 1.0d+00/DBLE(numh)

       do i = 1, nsec
           i0 = i + ninact + nact
           f(i, i) = orbmo(i0)
           do j = i, nsec
               j0 = j + ninact + nact
               do k = ninact + nact - numh + 1, ninact + nact

                   if (k > ninact + nact - 2 .and. mod(nelec, 2) == 1) then

                       f(i, j) = f(i, j) &
                           &  - 0.5d+00*DCMPLX(int2r_f1(i0, j0, k, k), int2i_f1(i0, j0, k, k))/DBLE(numh)
                       f(i, j) = f(i, j) &
                           &  + 0.5d+00*DCMPLX(int2r_f2(i0, k, k, j0), int2i_f2(i0, k, k, j0))/DBLE(numh)

                   else
                       f(i, j) = f(i, j) - DCMPLX(int2r_f1(i0, j0, k, k), int2i_f1(i0, j0, k, k))/DBLE(numh)
                       f(i, j) = f(i, j) + DCMPLX(int2r_f2(i0, k, k, j0), int2i_f2(i0, k, k, j0))/DBLE(numh)
                   end if

               end do

           end do
       end do

!        do i = 1, nsec
!           do j = i, nsec
!             if(ABS(fdmmy(i,j)-f(i,j))>1.0d-05) then
!                write(*,*)i,j,fdmmy(i,j),f(i,j),fdmmy(i,j)-f(i,j)
!             endif
!           enddo
!        enddo

       do i = 1, nsec
           do j = i, nsec
               f(j, i) = DCONJG(f(i, j))
           end do
       end do

       allocate (readmo(nbas*2, nbas*2, 2))
       allocate (itrfmo(nbas*2, nbas, 2))
       itrfmo = 0.0d+00

       open (15, file='r4dorbcoeff', status='old', form='unformatted')
       read (15, err=10) readmo
       close (15)

       itrfmo(1:nbas*2, 1:nbas, 1:2) = readmo(1:nbas*2, nbas + 1:nbas*2, 1:2)

       Do isym = 1, nsymrpa, 2
           nv = 0
           Do i = 1, nsec
               i0 = i + ninact + nact
               if (irpmo(i0) == isym) then
                   nv = nv + 1
               end if
           end do

           Allocate (mosym(nv))
           Allocate (fsym(nv, nv))
           fsym = 0.0d+00
           nv = 0
           Do i = 1, nsec
               i0 = i + ninact + nact
               if (irpmo(i0) == isym) then
                   nv = nv + 1
                   mosym(nv) = i
               end if
           end do

           Do i = 1, nv
               i0 = mosym(i)
               Do j = i, nv
                   j0 = mosym(j)
                   fsym(i, j) = f(i0, j0)
                   fsym(j, i) = DCONJG(f(i0, j0))
!                    write(*,*)fsym(i,j)
               end do
           end do
           Allocate (wsym(nv))
           wsym = 0.0d+00
           cutoff = .FALSE.
           thresd = 0.0d+00

           call cdiag(fsym, nv, nv, wsym, thresd, cutoff)

           Allocate (coeff(nbas*2, nv, 2))

           Do i = 1, nv
               i0 = mosym(i) + ncore + ninact + nact
               coeff(:, i, :) = itrfmo(:, i0, :)
           End do

           coeff(:, :, 1) = MATMUL(coeff(:, :, 1), fsym(:, :))
           coeff(:, :, 2) = MATMUL(coeff(:, :, 2), fsym(:, :))

           Do i = 1, nv
               i0 = mosym(i) + ncore + ninact + nact
               itrfmo(:, i0, :) = coeff(:, i, :)
           End do

! Kramers - pairs

           Do i = 1, nv
               i0 = mosym(i) + ncore + ninact + nact + 1
               coeff(:, i, :) = itrfmo(:, i0, :)
           End do

           coeff(:, :, 1) = MATMUL(coeff(:, :, 1), DCONJG(fsym(:, :)))
           coeff(:, :, 2) = MATMUL(coeff(:, :, 2), DCONJG(fsym(:, :)))

           Do i = 1, nv
               i0 = mosym(i) + ncore + ninact + nact + 1
               itrfmo(:, i0, :) = coeff(:, i, :)
           End do

! Kramers - pairs

           deallocate (coeff)

           Do i = 1, nv
               i0 = mosym(i)
               write (*, '(I4,F20.10)') i0, wsym(i)
           end do

           Do i = 1, nv
               i0 = mosym(i)
               write (*, *) ''
               write (*, *) 'new ', i0 + ninact + nact, 'th ms consists of '
               Do j = 1, nv
                   j0 = mosym(j)
                   if (ABS(fsym(j, i))**2 > 1.0d-03) then
                       write (*, '(I4,"  Weights ",F20.10)') j0 + ninact + nact, ABS(fsym(j, i))**2
                   end if
               end do
           end do
           deallocate (fsym)
           deallocate (wsym)
           deallocate (mosym)
       end do

       readmo(1:nbas*2, nbas + 1:nbas*2, 1:2) = itrfmo(1:nbas*2, 1:nbas, 1:2)

       open (15, file='r4dorbcoeff_ivo', status='unknown', form='unformatted')
       write (15) readmo
       close (15)
       goto 100

10     write (*, *) 'reading err of r4dorbcoeff'
!        deallocate(fdmmy)
       deallocate (readmo)
       deallocate (itrfmo)

100    write (*, *) 'fockivo end'
   end
