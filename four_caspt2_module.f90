! *+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@

MODULE four_caspt2_module

! *+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@

    Implicit NONE

    real*8       :: thres, tmem
    integer      :: norb, ndet, nelec, nroot, iroot, selectroot
    integer      :: ninact, nact, nsec ! ninact: the number of inactive spinors, nact: active, nsec: secondary
    integer      :: ncore, nbas
    character    :: date*8, time*10, ptgrp*6

    integer, allocatable :: idet(:), sp(:)
    real*8, allocatable :: cir(:, :), cii(:, :), eigen(:) ! cir:CI coefficients(real), cii:CI coefficients(imaginary)
    real*8, allocatable :: int2r(:), int2i(:), int2r_f1(:, :, :, :), int2r_f2(:, :, :, :)
    real*8, allocatable :: int2i_f1(:, :, :, :), int2i_f2(:, :, :, :)
    integer, allocatable :: indtwr(:, :, :, :), indtwi(:, :, :, :)
    logical :: realc, realcvec, debug, realf, evenelec

    integer                 ::  val(8), initdate, date0, date1
    real*8                  :: totalsec, inittime, tsec0, tsec1, tsec, eshift, sumc2, sumc2local
    complex*16              ::  coeff1

! from MORCONEE

!   real*8 :: ecore, enuc
    real*8 :: enuc
    double precision :: ecore
!   integer :: nsymrp, nsymrpa, multb(128,128), multb2(128,128), nmo, scfru
!   integer :: nmo
    integer :: nsymrp, nsymrpa, multb(128, 128), multb2(128, 128), nmo, scfru
    character :: repn(64)*14, repna(64)*4, repn_ty(64)*6
!   character :: repn(50)*14, repna(50)*4, repn_ty(50)*6
!   integer, allocatable :: irpmo(:), irpamo(:), indmo(:), indmor(:)
    integer, allocatable :: irpmo(:), irpamo(:)
    integer, allocatable :: indmo(:), indmor(:)
    real*8, allocatable  :: oner(:, :), onei(:, :)
    real*8, allocatable  :: orbmo(:), orb(:)
    real*8, allocatable  :: orbmocas(:), orbcas(:)

!   integer, allocatable ::multb_s(:,:), multb_d(:,:), multb_ds(:,:) ! This is for typart
!   integer, allocatable ::MULTB_DF(:,:), MULTB_DB(:,:), MULTB_SB(:,:)
    integer, allocatable ::multb_s(:, :), multb_d(:, :), multb_ds(:, :) ! This is for typart
    integer, allocatable ::MULTB_DF(:, :), MULTB_DB(:, :), MULTB_SB(:, :)

    real*8, allocatable :: eps(:)

    complex*16, allocatable :: f(:, :), itrfmo(:, :, :)

! Iwamuro modify
!   integer :: nelecd(64)
    integer :: nelecd(64), nfsym, nz1, norbt
    logical :: spfr, sfform, realonly

! Old Dirac
!       Write(UT_sys_ftmp) NMO,BREIT,ETOTAL
!       Write(UT_sys_ftmp) NSYMRP,(REPN(IRP),IRP=1,NSYMRP)
!       Write(UT_sys_ftmp) NSYMRPA,(REPNA(IRP),IRP=1,NSYMRPA*2)
!       Write(UT_sys_ftmp) ((MULTB(I,J),I=1,2*NSYMRPA),J=1,2*NSYMRPA)
!       Write(UT_sys_ftmp) (IRPMO(IMO),IRPAMO(IMO),ORBMO(IMO),IMO=1,NMO)
!       Write(UT_sys_ftmp) ((ONER(IMO,JMO),ONEI(IMO,JMO),JMO=1,NMO),IMO=1,NMO)

! Dirac
!       Write NSP, BREIT, ECORE, NFSYM, NZ1, SPFR, NORBT
!       Write NSYMRP, (REPNT(IRP),IRP=1,NSYMRP), (NELEC(IRP),IRP=1,NSYMRP)
!       Write (spinor_info(1,irp),irp=1,nfsym), (spinor_info(2,irp),irp=1,nfsym),&
!             (spinor_info(3,irp),irp=1,nfsym), (spinor_info(4,irp),irp=1,nfsym),&
!             (spinor_info(5,irp),irp=1,nfsym)
!       Write NREP, (REPNA(IRP),IRP=1,2*NREP)
!       Write ((MULTB(I,J,1),I=1,2*NREP),J=1,2*NREP)
!       Write (IRPMO(IMO),IRPAMO(IMO),ORBMO(IMO),IMO=1,NMO)
!       Write ((ONER(IMO,JMO), ONEI(IMO,JMO), JMO=1, NMO), IMO=1, NMO)

    ! Valiables for MPI
    ! params
    ! ierr:
    integer         :: ierr, nprocs, rank

    ! Run on r4dcasci_co.f90 to get the filenames for each process
    character(50)   :: mdcint_filename, mdcintnew, mdcint_debug, mdcint_int
end MODULE four_caspt2_module
