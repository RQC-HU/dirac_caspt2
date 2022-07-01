! *+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@

MODULE four_caspt2_module

! *+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@

    Implicit NONE

    real(8)      :: thres, tmem
    integer      :: norb, ndet, iroot ! ndet: the number of determinant
    !! =================================================
    !! Valiables of Input (active.inp)
    !! =================================================
    ! ninact        : The number of inactive spinors
    ! nact          : The number of active spinors
    ! nsec          : The number of secondary spinors
    ! nelec         : The number of electrons in active space
    ! nroot         : The number of roots
    ! selectroot    : Which root do you want to obtain
    ! totsym
    ! ncore         : The number of core orbitals
    ! nbas          : Basis set
    ! eshift        : Real shift
    ! ptgrp         : Point group symmetry
    ! dirac_version : DIRAC version
    integer         :: ninact, nact, nsec, nelec
    integer         :: nroot, selectroot
    integer         :: totsym, ncore, nbas
    real(8)         :: eshift
    character       :: ptgrp*6
    character       :: calctype*5 = "casci" ! dmrg or casci(default)
    integer         :: dirac_version
    integer         :: ras1_start, ras1_size, ras2_start, ras2_size, ras3_start, ras3_size
    integer         :: ras1_max_hole, ras3_max_elec, min_hole_ras1 = 0
    logical         :: is_ras1_configured, is_ras2_configured, is_ras3_configured
    integer, allocatable :: ras1_list(:), ras2_list(:), ras3_list(:)
    integer, parameter :: max_ras_spinor_num = 200

    character       :: date*8, time*10
    integer, allocatable :: idet(:), sp(:), idetr(:)

    !! =================================================
    !! Valiables of CI
    !! =================================================
    ! cir   : Real numbers of CI coeffients.
    ! cii   : Imaginary numbers of CI coeffients.
    ! eigen : Eigen values
    real(8), allocatable :: cir(:, :), cii(:, :), eigen(:)

    ! integer, allocatable :: indtwr(:, :, :, :), indtwi(:, :, :, :)
    !! =================================================
    !! Valiables of two electron integrals
    !! =================================================
    ! inttwr    : Real numbers of two electron integrals (Directly specify four indices to get the integral value)
    ! inttwi    : Imaginary numbers of two electron integrals (Directly specify four indices to get the integral value)
    ! int2r_f1
    ! int2r_f2
    ! int2i_f1
    ! int2i_f2
    real(8), allocatable :: inttwr(:, :, :, :), inttwi(:, :, :, :)
    real(8), allocatable :: int2r_f1(:, :, :, :), int2r_f2(:, :, :, :)
    real(8), allocatable :: int2i_f1(:, :, :, :), int2i_f2(:, :, :, :)

    logical :: realc, realcvec, debug, realf, evenelec
    !! =================================================
    !! Valiables of timer
    !! =================================================
    ! val       : Initial time of CASCI or CASPT2 programs. (ref. https://docs.oracle.com/cd/E19205-01/820-1201/aetcf/index.html)
    ! initdate  : Initial date = val(3) (min:1, max:31)
    ! inittime  : Initial time (sec) = hour*60^2+min*60+sec+millisec*0.001
    ! totalsec  : Equals to inittime
    ! date1     : Start date (for elapsed time measurement)
    ! date0     : End date (for elapsed time measurement)
    ! tsec1     : Start seconds (for elapsed time measurement)
    ! tsec0     : End seconds (for elapsed time measurement)
    integer     ::  val(8), initdate, date0, date1
    real(8)     :: totalsec, inittime, tsec0, tsec1, tsec

    !! ========================================
    !! Valiables of second pertubation energy
    !! ========================================
    ! sumc2     : Second pertubation energy(total) / a.u.
    ! sumclocal : Second pertubation energy(each subspace) / a.u.
    ! coeff1    : The coefficient for solvH (sumclocal = sumclocal + abs(coeff1)^2)
    real*8      ::  sumc2, sumc2local
    complex*16  ::  coeff1
! from MORCONEE

!   real*8 :: ecore, enuc
    real*8 :: enuc
    double precision :: ecore ! core energy
!   integer :: nsymrp, nsymrpa, multb(128,128), multb2(128,128), nmo, scfru
!   integer :: nmo
    integer :: nsymrp, nsymrpa, multb(128, 128), multb2(128, 128), nmo, scfru
    character :: repn(64)*14, repna(64)*4, repn_ty(64)*6
!   character :: repn(50)*14, repna(50)*4, repn_ty(50)*6
!   integer, allocatable :: irpmo(:), irpamo(:), indmo(:), indmor(:)
    integer, allocatable :: irpmo(:), irpamo(:) ! symmetry number of the specific mo
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

    !! ===================
    !! Valiables of MPI
    !! ===================
    ! ierr                  : Error code. if ierr is not 0, the MPI method was not working properly.
    ! nprocs                : Total number of processes in MPI (Available up to 10000)
    ! rank                  : Process number of MPI (0 <= rank <= nprocs-1)
    ! mdcint_filename       : MDCINT filenames for each MPI process
    ! mdcintnew             : MDCINTNEW filenames for each MPI process (e.g. MDCINTNEW8)
    ! mdcint_debug          : MDCINT_debug filenames for each MPI process (e.g. MDCINT_debug1)
    ! mdcint_int            : MDCINT_int filenames for each MPI process (e.g. MDCINT_int)
    ! a1int, a2int          : A subspace filenames for each MPI process (e.g. A1int2)
    ! bint                  : B subspace filenames for each MPI process (e.g. Bint10)
    ! c1int, c2int, c3int   : C subspace filenames for each MPI process (e.g. C1int1)
    ! d1int, d2int, d3int   : D subspace filenames for each MPI process (e.g. D1int)
    ! eint                  : E subspace filenames for each MPI process (e.g. Eint3)
    ! fint                  : F subspace filenames for each MPI process (e.g. Fint4)
    ! gint                  : G subspace filenames for each MPI process (e.g. Gint11)
    ! hint                  : H subspace filenames for each MPI process (e.g. Hint12)
    ! normal_output  : The unit number for normal output (default unit number = 3000, default file name = "caspt2.out")
    ! read_line_max : The number of lines to be read at a time when reading two-electron integral related files.(e.g. A1int, MDCINTNEW)
    integer         :: ierr, nprocs, rank
    character(50)   :: mdcint_filename, mdcintnew, mdcint_debug, mdcint_int
    character(50)   :: a1int, a2int, bint, c1int, c2int, c3int, d1int, d2int, d3int, eint, fint, gint, hint
    integer, parameter :: normal_output = 3000, read_line_max = 1000

    ! Test for MDCINT_READ_COUNT
    integer :: casci_mdcint_cnt, caspt2_mdcint_cnt, caspt2_mdcint_cnt2, simple_loop
end MODULE four_caspt2_module
