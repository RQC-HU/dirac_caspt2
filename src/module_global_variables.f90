! *+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@

MODULE module_global_variables

! *+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@+*+@+*+@*+@+*+@+*+@
    use module_dict
    Implicit NONE

    !! =================================================
    !! Variables of Input (active.inp)
    !! =================================================
    ! ninact        : The number of inactive spinors
    ! nact          : The number of active spinors
    ! nsec          : The number of secondary spinors
    ! nelec         : The number of electrons in active space
    ! nroot         : The number of roots
    ! selectroot    : Which root do you want to obtain
    ! totsym
    ! ncore         : The number of core orbitals
    ! eshift        : Real shift
    ! dirac_version : DIRAC version
    ! mdcint_scheme : MDCINT **MOLTRA>.SCHEME value (https://diracprogram.org/doc/release-23/manual/moltra.html#scheme)
    integer         :: ninact, nact, nsec, nelec
    integer         :: nroot = 10, selectroot = 1
    integer         :: totsym, ncore
    real(8)         :: eshift = 0.0d+00 ! default: 0.0
    integer         :: dirac_version
    integer         :: mdcint_scheme = 0
    integer         :: ras1_start, ras2_start, ras3_start
    integer         :: ras1_size = 0, ras2_size = 0, ras3_size = 0
    integer         :: ras1_max_hole, ras3_max_elec, min_hole_ras1 = 0
    logical         :: skip_mdcint = .false., enable_restart = .false.
    integer, allocatable :: ras1_list(:), ras2_list(:), ras3_list(:)
    integer         :: nhomo = 0  ! Default value of nhomo is zero. If you want to specify the value, please use the input file.
    integer, parameter :: max_ras_spinor_num = 200, max_i4 = huge(0_4) ! 4byte integer max value
    ! vcut_mo_num: The number of virtual orbitals in each fermion irreducible representation
    ! occ_mo_num: The number of occupied orbitals in each fermion irreducible representation
    integer         :: vcut_mo_num(2) = 0, occ_mo_num(2) = 0
    logical         :: inversion = .false., no_inversion = .false., is_scheme_set = .false.
    logical         :: debug = .false. ! debugprint option
    integer, parameter :: default_scheme_dirac22_or_earlier = 6, default_scheme_dirac23_or_later = 4

    !! =================================================
    !! Variables of CI
    !! =================================================
    ! cir   : Real numbers of CI coeffients.
    ! cii   : Imaginary numbers of CI coeffients.
    ! eigen : Eigen values
    real(8), allocatable :: cir(:, :), cii(:, :), eigen(:)

    !! =================================================
    !! Variables of two electron integrals
    !! =================================================
    ! inttwr    : Real numbers of two electron integrals (range: (1:ninact+nact, 1:ninact+nact, 1:ninact+nact, 1:ninact+nact))
    ! inttwi    : Imaginary numbers of two electron integrals (range: (1:ninact+nact, 1:ninact+nact, 1:ninact+nact, 1:ninact+nact))
    ! int2r_f1  : Real numbers of two electron integrals (range: (ninact+nact+1:ninact+nact+nsec, ninact+nact+1:ninact+nact+nsec, 1:ninact+nact, 1:ninact+nact))
    ! int2r_f2  : Real numbers of two electron integrals (range: (ninact+nact+1:ninact+nact+nsec, 1:ninact+nact, 1:ninact+nact, ninact+nact+1:ninact+nact+nsec))
    ! int2i_f1  : Imaginary numbers of two electron integrals (range: (ninact+nact+1:ninact+nact+nsec, ninact+nact+1:ninact+nact+nsec, 1:ninact+nact, 1:ninact+nact))
    ! int2i_f2  : Imaginary numbers of two electron integrals (range: (ninact+nact+1:ninact+nact+nsec, 1:ninact+nact, 1:ninact+nact, ninact+nact+1:ninact+nact+nsec))
    real(8), allocatable, target :: inttwr(:, :, :, :), inttwi(:, :, :, :)
    real(8), allocatable, target :: int2r_f1(:, :, :, :), int2r_f2(:, :, :, :)
    real(8), allocatable, target :: int2i_f1(:, :, :, :), int2i_f2(:, :, :, :)

    !! =================================================
    !! Variables of timer
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
    !! Variables of second pertubation energy
    !! ========================================
    ! sumc2     : Second pertubation energy(total) / a.u.
    ! sumclocal : Second pertubation energy(each subspace) / a.u.
    ! coeff1    : The coefficient for solvH (sumclocal = sumclocal + abs(coeff1)^2)
    real(8)     :: sumc2 = 0.0d+00, e2all = 0.0d+00
    real(8)     :: sumc2_subspace(8) = 0.0d+00, e2_subspace(8) = 0.0d+00
    complex*16  ::  coeff1

    !! ========================================
    !! Variables of restart file
    !! ========================================
    character(len=1) :: next_subspace ! A-H is the valid subspace (if not, the program will stop with error message)

    !! ====================================================================================================================================
    !! Variables of MRCONEE (A file stores 1-electron integrals, symmetry information, multiplication table etc. that is created by DIRAC)
    !! ====================================================================================================================================
    double precision :: ecore ! core energy
    integer :: nmo, scfru
    integer, allocatable :: irpamo(:) ! symmetry number of the specific mo
    integer, allocatable :: indmo_cas_to_dirac(:) ! MO index transformation from CASPT2 to DIRAC (irrep: irreducible representation) order. Ex: indmo_cas_to_order(cas_index) = dirac_index [old name]: indmo
    integer, allocatable :: indmo_dirac_to_cas(:) ! MO index transformation from DIRAC (irrep) to CASPT2 order. Ex: indmo_dirac_to_order(dirac_index) = cas_index [old name]: indmor
    real(8), allocatable :: one_elec_int_r(:, :) ! one-electron integral in CASPT2 order (real part) [old name]: oner
    real(8), allocatable :: one_elec_int_i(:, :) ! one-electron integral in CAPST2 order (imaginally part) [old name]: onei
    real(8), allocatable :: dirac_mo_energy(:) ! MO energy (a.u.) (DIRAC order)
    real(8), allocatable :: caspt2_mo_energy(:) ! MO energy (a.u.) (CASPT2 order) [old name]: orbmo
    ! Symmetry
    integer :: nsymrpa ! number of irreducible representation (abelian group)
    character :: repna(64)*4 ! irreducible representation assignment (abelian group)
    ! Multiplication table for symmetry
    integer :: multb(128, 128), multb2(128, 128)
    integer, allocatable ::multb_s(:, :), multb_d(:, :), multb_ds(:, :)
    real(8) :: hf_energy_mrconee ! HF energy in MRCONEE

    !! ========================================
    !! Variables of IVO calculation
    !! ========================================
    complex*16, allocatable :: itrfmo(:, :)

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
    !! Variables of MPI
    !! ===================
    ! ierr                  : Error code. if ierr is not 0, the MPI method was not working properly.
    ! nprocs                : Total number of processes in MPI (Available up to 10000)
    ! rank                  : Process number of MPI (0 <= rank <= nprocs-1)
    ! mdcint_filename       : MDCINT filenames for each MPI process
    ! mdcintnew             : MDCINTNEW filenames for each MPI process (e.g. MDCINTNEW8)
    ! a1int, a2int          : A subspace filenames for each MPI process (e.g. A1int2)
    ! bint                  : B subspace filenames for each MPI process (e.g. Bint10)
    ! c1int, c2int, c3int   : C subspace filenames for each MPI process (e.g. C1int1)
    ! d1int, d2int, d3int   : D subspace filenames for each MPI process (e.g. D1int)
    ! eint                  : E subspace filenames for each MPI process (e.g. Eint3)
    ! fint                  : F subspace filenames for each MPI process (e.g. Fint4)
    ! gint                  : G subspace filenames for each MPI process (e.g. Gint11)
    ! hint                  : H subspace filenames for each MPI process (e.g. Hint12)
    integer :: ierr, nprocs, rank
    character(:), allocatable :: mdcint_filename, mdcintnew
    character(:), allocatable :: a1int, a2int, bint, c1int, c2int, c3int, d1int, d2int, d3int, eint, fint, gint, hint

    !! ================
    !! Others
    !! ================
    character       :: date*8, time*10 ! Timing information stored in MDCINT (DIRAC 2-electron integral file)
    integer, allocatable :: space_idx(:) ! Given the spinor index, return which space it belongs to. (1: inactive, 2: active, 3: secondary) [old name]: sp
    type(dict) :: dict_cas_idx ! Dictionary(Key: the number of CAS placement val: an integer value representing the CAS placement) [old name]: idet
    type(dict) :: dict_cas_idx_reverse ! Dictionary(key: an integer representing the CAS placement, val: the position in the CAS placement ordering) [old name]: idetr
    ! Global index of inactive, active, secondary spinors
    ! global_inact_start = 1, global_inact_end = ninact
    ! global_act_start = ninact + 1, global_act_end = ninact + nact
    ! global_sec_start = ninact + nact + 1, global_sec_end = ninact + nact + nsec
    integer    :: global_inact_start, global_inact_end, global_act_start, global_act_end, global_sec_start, global_sec_end
    ! Epsilon
    real(8), allocatable :: eps(:)
    ! Fock matrix
    complex*16, allocatable :: fock_cmplx(:, :)
    real(8), allocatable :: fock_real(:, :)
    real(8)      :: tmem ! Total allocated memory size of allcatable variables in Bytes
    real(8)      :: smat_lin_dep_threshold = 1.0d-08  ! threshold for removing linear dependencies in S-matrix
    real(8), parameter  :: bmat_no_cutoff = 0.0d+00   ! No threshold for B-matrix (all elements are used)
    real(8), parameter  :: global_threshold = 1.0d-15 ! Threshold for removing small elements
    integer      :: ndet, iroot ! ndet: the number of CAS configuration
    integer, parameter :: len_convert_int_to_chr = 30 ! Length of the string for converting integer to character
    integer, parameter :: cimat_key_size = 30 ! Length of the key in CIMAT file
    logical :: dirac_32bit_build = .false. ! If .true., the DIRAC is built in 32-bit mode (default: .false.(64-bit mode))
end MODULE module_global_variables
