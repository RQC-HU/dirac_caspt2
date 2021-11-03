! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE read1mo (filename) ! one-electron MO integrals in MRCONEE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        use four_caspt2_module

        Implicit NONE

    	integer :: mrconee
     	character*50,intent(in) :: filename
!        integer :: j0, j, i, i0, i1
!        integer :: k0, l0, ii, jj, kk, ll, nmom
        integer :: j0, j, i, i0, i1
        integer :: k0, l0, ii, jj, kk, ll, nmom
!        real*8, allocatable :: roner(:,:), ronei(:,:)
        double precision, allocatable :: roner(:,:), ronei(:,:)

!       Write(UT_sys_ftmp) NMO,BREIT,ECORE
!       Write(UT_sys_ftmp) NSYMRP,(REPN(IRP),IRP=1,NSYMRP)
!       Write(UT_sys_ftmp) NSYMRPA,(REPNA(IRP),IRP=1,NSYMRPA*2)
!       Write(UT_sys_ftmp) ((MULTB(I,J),I=1,2*NSYMRPA),J=1,2*NSYMRPA)
!       Write(UT_sys_ftmp) (IRPMO(IMO),IRPAMO(IMO),ORBMO(IMO),IMO=1,NMO)
!       Write(UT_sys_ftmp) ((ONER(IMO,JMO),ONEI(IMO,JMO),JMO=1,NMO),IMO=1,NMO)

        write(*,*)'Enter read1mo'

    	mrconee=10

        realc = .true.

        Allocate ( roner (nmo,nmo)); Call memplus(KIND(roner ),SIZE(roner ),1)
        Allocate ( ronei (nmo,nmo)); Call memplus(KIND(ronei ),SIZE(ronei ),1)

    	open( mrconee, file=trim(filename), status='old', form='unformatted', err=10)
        rewind (mrconee)
        read(mrconee,err=10)
        read(mrconee,err=10)
        read(mrconee,err=10)
        read(mrconee,err=10)
        read(mrconee,err=10)
        read(mrconee,err=10) ((roner(i0,j0),ronei(i0,j0),j0=1,nmo),i0=1,nmo)

! Iwamuro modify
          do i0 = 1, nmo
          do j0 =1, nmo
!             Wrpite(*,'(2I4,2X,2F8.4)') i0, j0, RONER(i0,j0),RONEI(i0,j0)
          end do
       end do
    	close (mrconee)

        nmom = ninact + nact + nsec
        Allocate ( oner (nmom,nmom)); Call memplus(KIND(oner ),SIZE(oner ),1)
        Allocate ( onei (nmom,nmom)); Call memplus(KIND(onei ),SIZE(onei ),1)

        do i0 = 1, nmom
        do j0 = 1, nmom
           oner(i0,j0) = roner(indmo(i0),indmo(j0))
           onei(i0,j0) = ronei(indmo(i0),indmo(j0))
        enddo
        enddo

        deallocate (roner); Call memminus(KIND(roner ),SIZE(roner ),1)
        deallocate (ronei); Call memminus(KIND(ronei ),SIZE(ronei ),1)

        write(*,*) realc,'realc'
        goto 1000

 10     write(*,*) 'err 10 mo1'
        go to 1000
 11     write(*,*) 'err 11 mo1'
        go to 1000

 1000   end subroutine read1mo
