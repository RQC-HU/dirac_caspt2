! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE read1mo_co(filename) ! one-electron MO integrals in moint1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       use four_caspt2_module

       Implicit NONE

       integer :: mrconee, isp
       character*50, intent(in) :: filename
       integer :: j0, j, i, i0, i1
       integer :: k0, l0, ii, jj, kk, ll, nmom

!        real*8, allocatable :: roner(:,:,:), ronei(:,:,:)
       double precision, allocatable :: roner(:, :, :), ronei(:, :, :)

       if (rank == 0) then
           write (normaloutput, *) 'Enter read1mo_co'
       end if
       mrconee = 10

!  Write(UT_sys_ftmp) NMO,UT_molinp_atm_enm - DELETE, &
!                     BREIT,ETOTAL,scfru
!  Write(UT_sys_ftmp) NSYMRP,(REPN(IRP),IRP=1,NSYMRP)
!  Write(UT_sys_ftmp) ((UT_ptgsym_table_single(IJ,II),UT_ptgsym_table_double(IJ,II),IJ=0,NSYMRP-1),II=0,NSYMRP-1)
!  Write(UT_sys_ftmp) ((IRPMO(IMO,isp),ORBMO(IMO,isp), &
!                       UTCHEMIMO1(IMO,isp),UTCHEMIMO2(IMO,isp),IMO=1,NMO),isp=1,scfru)
!  Write(UT_sys_ftmp) (((ONE(JMO,IMO,isp),JMO=1,NMO),IMO=1,NMO),isp=1,scfru)

       realc = .true.

       Allocate (roner(nmo, nmo, scfru)); Call memplus(KIND(roner), SIZE(roner), 1)
       Allocate (ronei(nmo, nmo, scfru)); Call memplus(KIND(ronei), SIZE(ronei), 1)

       open (mrconee, file=trim(filename), status='old', form='unformatted', err=10)
       rewind (mrconee)
       read (mrconee, err=10)
       read (mrconee, err=10)
       read (mrconee, err=10)
       read (mrconee, err=10)
       read (mrconee, err=10)
       read (mrconee, err=10) (((roner(i0, j0, isp), ronei(i0, j0, isp), j0=1, nmo), i0=1, nmo), isp=1, scfru)

!       do i0 = 1, nmo
!          do j0 =1, nmo
!             Write(*,'(2I4,2X,2F8.4)') i0, j0, RONER(i0,j0,1),RONEI(i0,j0,1)
!          end do
!       end do
       close (mrconee)

       nmom = ninact + nact + nsec
       Allocate (oner(nmom, nmom)); Call memplus(KIND(oner), SIZE(oner), 1)
       Allocate (onei(nmom, nmom)); Call memplus(KIND(onei), SIZE(onei), 1)

!Iwamuro modify

       do i0 = 1, nmom
       do j0 = 1, nmom
!           oner(i0,j0) = roner(i0,j0,1)
!           onei(i0,j0) = ronei(i0,j0,1)
           oner(indmor(i0), indmor(j0)) = roner(i0, j0, 1) ! using alpha component for a while
           onei(indmor(i0), indmor(j0)) = ronei(i0, j0, 1)
       end do
       end do

       deallocate (roner); Call memminus(KIND(roner), SIZE(roner), 1)
       deallocate (ronei); Call memminus(KIND(ronei), SIZE(ronei), 1)

       if (rank == 0) then
           write (normaloutput, *) realc, 'realc'
       end if
       goto 1000

10     if (rank == 0) then
           write (normaloutput, *) 'err 10 mo1'
       end if
       go to 1000
11     if (rank == 0) then
           write (normaloutput, *) 'err 11 mo1'
       end if
       go to 1000

1000 end subroutine read1mo_co