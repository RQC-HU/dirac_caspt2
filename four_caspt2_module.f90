! *+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@

MODULE four_caspt2_module

! *+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@

    Implicit NONE

    real*8       :: thres, tmem
    integer      :: norb, ndet, nelec, nroot, iroot, selectroot ! ndet: the number of determinant
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
    integer, allocatable :: indmo(:), indmor(:) ! index of MO
    real*8, allocatable  :: oner(:, :), onei(:, :) ! one-electron integral (real,imaginal)
    real*8, allocatable  :: orbmo(:), orb(:)
    real*8, allocatable  :: orbmocas(:), orbcas(:)

!   integer, allocatable ::multb_s(:,:), multb_d(:,:), multb_ds(:,:) ! This is for typart
!   integer, allocatable ::MULTB_DF(:,:), MULTB_DB(:,:), MULTB_SB(:,:)
    integer, allocatable ::multb_s(:, :), multb_d(:, :), multb_ds(:, :) ! This is for typart
    integer, allocatable ::MULTB_DF(:, :), MULTB_DB(:, :), MULTB_SB(:, :)

    real*8, allocatable :: eps(:)

    complex*16, allocatable :: f(:, :), itrfmo(:, :, :) ! f: fock matrix

! Iwamuro modify
!   integer :: nelecd(64)
    integer :: nelecd(64), nfsym, nz1, norbt
    logical :: spfr, sfform, realonly ! realonly : If it is true, only real numbers are written in MDCINT.

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

    !! ---------------------------
    !! Valiables for MPI
    !! ---------------------------
    ! @params
    ! ierr   : Error code. if ierr is not 0, the MPI method is not working properly.
    ! nprocs : The number of processes in MPI
    ! rank   : Process number of MPI
    integer         :: ierr, nprocs, rank

    ! Run on r4dcasci_co.f90 and r4dcaspt2_tra_co.f90 to get the MDCINT filenames for each process
    character(50)   :: mdcint_filename, mdcintnew, mdcint_debug, mdcint_int

    ! Run on r4dcaspt2_tra_co.f90 to get the subspace filenames for each process
    character(50)   :: a1int, a2int, bint, c1int, c2int, c3int, d1int, d2int, d3int, eint, fint, gint, hint

    ! Test for MDCINT_READ_COUNT
    integer :: casci_mdcint_cnt, caspt2_mdcint_cnt, caspt2_mdcint_cnt2, simple_loop

    ! Unit number for normal output (caspt2.out)
    integer,parameter :: normaloutput = 3000
end MODULE four_caspt2_module
