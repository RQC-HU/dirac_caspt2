! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE readorb_enesym_co(filename) ! orbital energies in r4dmoin1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module
    use module_file_manager
    use module_sort_swap
    Implicit NONE

    integer :: mrconee_unit, IMO, IRP
    character*50, intent(in) :: filename
    integer :: i0, j0, k0, i, j, m, isym, jsym, ksym, iostat
    integer, allocatable :: dammo(:), UTCHEMIMO1(:, :), UTCHEMIMO2(:, :)
    integer, allocatable :: SD(:, :), DS(:, :)
    logical :: breit, is_end_of_file

!  Write(UT_sys_ftmp) NMO,UT_molinp_atm_enm - DELETE, &
!                     BREIT,ETOTAL,scfru
!  Write(UT_sys_ftmp) NSYMRP,(REPN(IRP),IRP=1,NSYMRP)
!  Write(UT_sys_ftmp) ((UT_ptgsym_table_single(IJ,II),UT_ptgsym_table_double(IJ,II),IJ=0,NSYMRP-1),II=0,NSYMRP-1)
!  Write(UT_sys_ftmp) ((IRPMO(IMO,isp),ORBMO(IMO,isp), &
!                       UTCHEMIMO1(IMO,isp),UTCHEMIMO2(IMO,isp),IMO=1,NMO),isp=1,scfru)
!  Write(UT_sys_ftmp) (((ONE(JMO,IMO,isp),JMO=1,NMO),IMO=1,NMO),isp=1,scfru)

    mrconee_unit = default_unit

    call open_unformatted_file(unit=mrconee_unit, file=trim(filename), status='old', optional_action='read')

    Read (mrconee_unit, iostat=iostat) NMO, BREIT, ECORE  ! NMO is nbas - ncore
    call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
    if (is_end_of_file) then
        print *, 'Error: error in reading NMO, BREIT, ECORE (end of file reached)'
        print *, 'iostat = ', iostat
        stop
    end if

    if (rank == 0) then
        print *, 'NMO, BREIT, ECORE, 1  ! NMO is nbas - ncore'
        print *, NMO, BREIT, ECORE, 1  ! NMO is nbas - ncore
    end if
!Iwamuro modify
    scfru = 1

    allocate (IRPMO(1:NMO))
    allocate (ORBMO(1:NMO))
    allocate (UTCHEMIMO1(1:NMO, 1:scfru))
    allocate (UTCHEMIMO2(1:NMO, 1:scfru))
    Call memplus(size(IRPMO), kind(IRPMO), 1)
    Call memplus(size(ORBMO), kind(ORBMO), 1)
    Call memplus(size(UTCHEMIMO1), kind(UTCHEMIMO1), 1)
    Call memplus(size(UTCHEMIMO2), kind(UTCHEMIMO2), 1)

    Read (mrconee_unit, iostat=iostat) NSYMRP, (REPN(IRP), IRP=1, NSYMRP)                         ! IRs chars
    call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
    if (is_end_of_file) then
        print *, 'Error: error in reading NSYMRP, REPN (end of file reached)'
        print *, 'iostat = ', iostat
        stop
    end if

    if (rank == 0) then
        print *, ' NSYMRP, (REPN(IRP),IRP=1,NSYMRP)                         ! IRs chars'
        print *, NSYMRP, (REPN(IRP), IRP=1, NSYMRP)                         ! IRs chars
    end if
!Iwamuro modify
    Read (mrconee_unit, iostat=iostat) nsymrpa, (repna(i0), i0=1, nsymrpa*2)
    call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
    if (is_end_of_file) then
        print *, 'Error: error in reading nsymrpa, repna (end of file reached)'
        print *, 'iostat = ', iostat
        stop
    end if

    if (rank == 0) then
        print *, nsymrpa, (repna(i0), i0=1, nsymrpa*2)
    end if
    allocate (MULTB_S(1:NSYMRPA, 1:NSYMRPA))
    allocate (MULTB_D(1:NSYMRPA, 1:NSYMRPA))  ! dagger
    allocate (MULTB_DF(1:NSYMRPA, 1:NSYMRPA)) ! forward
    allocate (MULTB_DB(1:NSYMRPA, 1:NSYMRPA)) ! backward
    allocate (MULTB_DS(1:NSYMRPA, 1:NSYMRPA))
    allocate (MULTB_SB(1:NSYMRPA, 1:NSYMRPA))
    Call memplus(size(MULTB_S), kind(MULTB_S), 1)
    Call memplus(size(MULTB_D), kind(MULTB_D), 1)
    Call memplus(size(MULTB_DS), kind(MULTB_DS), 1)
    Call memplus(size(MULTB_DF), kind(MULTB_DF), 1)
    Call memplus(size(MULTB_DB), kind(MULTB_DB), 1)
    Call memplus(size(MULTB_SB), kind(MULTB_SB), 1)
    allocate (DS(1:NSYMRPA, 1:NSYMRPA))
    allocate (SD(1:NSYMRPA, 1:NSYMRPA))

!     Read(UT_sys_ftmp) ((MULTB_S(J,I),MULTB_D(J,I),J=0,NSYMRP-1),I=0,NSYMRP-1)
!     Read(UT_sys_ftmp) ((IRPMO(IMO,isp),ORBMO(IMO,isp), &
!                         UTCHEMIMO1(IMO,isp),UTCHEMIMO2(IMO,isp), &
!                         IMO=1,NMO),isp=1,scfru)                                ! orbital energies <= used here

!    Read(mrconee_unit) ((MULTB_S(J,I),MULTB_D(J,I),J=1,NSYMRP),I=1,NSYMRP)
!    Read(mrconee_unit) ((IRPMO(IMO),ORBMO(IMO), &
!                         UTCHEMIMO1(IMO,isp),UTCHEMIMO2(IMO,isp), &
!                         IMO=1,NMO),isp=1,scfru)                                ! orbital energies <= used here

    Read (mrconee_unit, iostat=iostat) ((multb(i0, j0), i0=1, 2*nsymrpa), j0=1, 2*nsymrpa)
    call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
    if (is_end_of_file) then
        print *, 'Error: error in reading multb (end of file reached)'
        print *, 'iostat = ', iostat
        stop
    end if

!    Read(mrconee_unit) (IRPMO(IMO),ORBMO(IMO),IMO=1,NMO)                             ! orbital energies <= used here
!Iwamuro modify
!    Do IMO=1,NMO
!      Write(*,*) IRPMO(IMO),ORBMO(IMO)
!    Enddo

!    CLOSE(mrconee_unit)

!----------------------------------------------------------------------------------------

!Iwamuro modify

! create MULTB2

    Do i0 = 1, 2*nsymrpa
        Do j0 = 1, 2*nsymrpa
            k0 = MULTB(i0, j0)
            MULTB2(i0, k0) = j0
        End do
    End do

    if (rank == 0) then
        print *, 'MULTB'

        Do i0 = 1, 2*nsymrpa
            print '(400I3)', (MULTB(i0, j0), j0=1, 2*nsymrpa)
        End do

        print *, 'MULTB2'

        Do i0 = 1, 2*nsymrpa
            print '(400I3)', (MULTB2(i0, j0), j0=1, 2*nsymrpa)
        End do
        print *, 'end multb1,2'
    end if
! create MULTB_S, MULTB_D

    Do i0 = 1, nsymrpa
        Do j0 = 1, nsymrpa
            MULTB_S(i0, j0) = MULTB(i0 + nsymrpa, j0 + nsymrpa)
            !    MULTB_D(i0, j0) = MULTB2(j0, i0)
            MULTB_D(i0, j0) = MULTB2(i0, j0)
        End do
    End do

    MULTB_S = MULTB_S - nsymrpa
    MULTB_D = MULTB_D - nsymrpa

    if (rank == 0) then
        print *, 'MULTB_S'

        Do i0 = 1, nsymrpa
            print '(200I3)', (MULTB_S(i0, j0), j0=1, nsymrpa)
        End do

        print *, 'MULTB_D'

        Do i0 = 1, nsymrpa
            print '(200I3)', (MULTB_D(i0, j0), j0=1, nsymrpa)
        End do
    end if

!----------------------------------------------------------------------------------------

    deallocate (UTCHEMIMO1); Call memminus(KIND(UTCHEMIMO1), SIZE(UTCHEMIMO1), 1)
    deallocate (UTCHEMIMO2); Call memminus(KIND(UTCHEMIMO2), SIZE(UTCHEMIMO2), 1)

    Allocate (sp(1:nmo)); Call memplus(KIND(sp), SIZE(sp), 1)
    sp(1:ninact) = 1
    sp(ninact + 1:ninact + nact) = 2
    sp(ninact + nact + 1:ninact + nact + nsec) = 3
    sp(ninact + nact + nsec + 1:nmo) = 4

    if (rank == 0) then
        print *, 'moint1 is closed.'
    end if
!     irpmo(1:imo) = irpmo(1:imo) + 1       ! irrep starts from 1

! Create MULTB_DF, MULTB_SB and MULTB_DB

    If (trim(ptgrp) == 'C1') then  ! C1 symmetry
        if (rank == 0) then
            print *, 'If pgsym=c1, this route goes through.'
        end if
        NSYMRP = 1
        NSYMRPA = 1
        REPNA(1) = 'a'; REPNA(2) = 'a'

        SD(1, 1) = 1
        DS(1, 1) = 1
        MULTB_DS = 1
        irpmo = 1
    else ! trim(prgrp) /= 'C1'
        !Iwamuro modify
        ! Not C1 symmetry
        if (rank == 0) then
            print *, 'if pgsym=c1, this route does not go through. '
        end if

        Do jsym = 1, nsymrpa
            Do isym = 1, nsymrpa - 1, 2
                MULTB_DF(isym + 1, jsym) = MULTB_D(isym, jsym)
                MULTB_DF(isym, jsym) = MULTB_D(isym + 1, jsym)
            End do
        End do

        Do jsym = 1, nsymrpa
            Do isym = 1, nsymrpa
                ksym = MULTB_DF(isym, jsym)
                MULTB_DB(isym, ksym) = jsym
            End do
        End do

        Do jsym = 1, nsymrpa
            Do isym = 1, nsymrpa
                ksym = MULTB_S(isym, jsym)
                MULTB_SB(isym, ksym) = jsym
            End do
        End do

        if (rank == 0) then
            print *, 'MULTB_SB'
            Do I = 1, nsymrpa
                print '(50I3)', (MULTB_SB(I, J), J=1, NSYMRPA)
            End do

            print *, 'MULTB_DF'
            Do I = 1, nsymrpa
                print '(50I3)', (MULTB_DF(I, J), J=1, NSYMRPA)
            End do

            print *, 'MULTB_DB'
            Do I = 1, nsymrpa
                print '(50I3)', (MULTB_DB(I, J), J=1, NSYMRPA)
            End do
        end if
!     Write(*,'("UTCHEMIMO1",50I3)') (UTCHEMIMO1(IMO,1),IMO=1,nmo)
!     Write(*,'("UTCHEMIMO2",50I3)') (UTCHEMIMO2(IMO,1),IMO=1,nmo)
!     Write(*,'("IRPMO",50I3)') (IRPMO(IMO),IMO=1,nmo)
!     Write(*,'("ORBMO",20F10.4)') (ORBMO(IMO),IMO=1,nmo)

! Create MULTB_DS, MULTB_SD

        If (trim(ptgrp) == 'C32h') then
            REPNA(1) = '1e1/2g'; REPNA(2) = '2e1/2g'; REPNA(3) = '1e3/2g'; REPNA(4) = '2e3/2g'
            REPNA(5) = '1e5/2g'; REPNA(6) = '2e5/2g'; REPNA(7) = '1e7/2g'; REPNA(8) = '2e7/2g'
            REPNA(9) = '1e9/2g'; REPNA(10) = '2e9/2g'; REPNA(11) = '1e11/2g'; REPNA(12) = '2e11/2g'
            REPNA(13) = '1e13/2g'; REPNA(14) = '2e13/2g'; REPNA(15) = '1e15/2g'; REPNA(16) = '2e15/2g'
            REPNA(17) = '1e1/2u'; REPNA(18) = '2e1/2u'; REPNA(19) = '1e3/2u'; REPNA(20) = '2e3/2u'
            REPNA(21) = '1e5/2u'; REPNA(22) = '2e5/2u'; REPNA(23) = '1e7/2u'; REPNA(24) = '2e7/2u'
            REPNA(25) = '1e9/2u'; REPNA(26) = '2e9/2u'; REPNA(27) = '1e11/2u'; REPNA(28) = '2e11/2u'
            REPNA(29) = '1e13/2u'; REPNA(30) = '2e13/2u'; REPNA(31) = '1e15/2u'; REPNA(32) = '2e15/2u'

            REPNA(33) = 'ag    '; REPNA(34) = 'bg    '; REPNA(35) = '1e1g  '; REPNA(36) = '2e1g  '
            REPNA(37) = '1e2g  '; REPNA(38) = '2e2g  '; REPNA(39) = '1e3g  '; REPNA(40) = '2e3g  '
            REPNA(41) = '1e4g  '; REPNA(42) = '2e4g  '; REPNA(43) = '1e5g  '; REPNA(44) = '2e5g  '
            REPNA(45) = '1e6g  '; REPNA(46) = '2e6g  '; REPNA(47) = '1e7g  '; REPNA(48) = '2e7g  '
            REPNA(49) = 'au    '; REPNA(50) = 'bu    '; REPNA(51) = '1e1u  '; REPNA(52) = '2e1u  '
            REPNA(53) = '1e2u  '; REPNA(54) = '2e2u  '; REPNA(55) = '1e3u  '; REPNA(56) = '2e3u  '
            REPNA(57) = '1e4u  '; REPNA(58) = '2e4u  '; REPNA(59) = '1e5u  '; REPNA(60) = '2e5u  '
            REPNA(61) = '1e7u  '; REPNA(62) = '2e7u  '; REPNA(63) = '1e9u  '; REPNA(64) = '2e9u  '

            Do i = 1, nsymrpa/2
                Do j = 1, nsymrpa/2
                    SD(i, j) = MULTB(i + nsymrpa, j)
                End do
            End do

            Do i = 1, nsymrpa/2
                Do j = 1, nsymrpa/2
                    SD(i, j + nsymrpa/2) = SD(i, j) + nsymrpa/2
                End do
            End do

            Do i = 1, nsymrpa/2
                Do j = 1, nsymrpa/2
                    SD(i + nsymrpa/2, j) = SD(i, j + nsymrpa/2)
                End do
            End do

            Do i = 1, nsymrpa/2
                Do j = 1, nsymrpa/2
                    SD(i + nsymrpa/2, j + nsymrpa/2) = SD(i, j)
                End do
            End do

        Elseif (trim(ptgrp) == 'C32') then
            REPNA(1) = '1e1/2'; REPNA(2) = '2e1/2'; REPNA(3) = '1e3/2'; REPNA(4) = '2e3/2'
            REPNA(5) = '1e5/2'; REPNA(6) = '2e5/2'; REPNA(7) = '1e7/2'; REPNA(8) = '2e7/2'
            REPNA(9) = '1e9/2'; REPNA(10) = '2e9/2'; REPNA(11) = '1e11/2'; REPNA(12) = '2e11/2'
            REPNA(13) = '1e13/2'; REPNA(14) = '2e13/2'; REPNA(15) = '1e15/2'; REPNA(16) = '2e15/2'
            REPNA(17) = '1e1/2'; REPNA(18) = '2e1/2'; REPNA(19) = '1e3/2'; REPNA(20) = '2e3/2'
            REPNA(21) = '1e5/2'; REPNA(22) = '2e5/2'; REPNA(23) = '1e7/2'; REPNA(24) = '2e7/2'
            REPNA(25) = '1e9/2'; REPNA(26) = '2e9/2'; REPNA(27) = '1e11/2'; REPNA(28) = '2e11/2'
            REPNA(29) = '1e13/2'; REPNA(30) = '2e13/2'; REPNA(31) = '1e15/2'; REPNA(32) = '2e15/2'

            REPNA(33) = 'a    '; REPNA(34) = 'b    '; REPNA(35) = '1e1  '; REPNA(36) = '2e1  '
            REPNA(37) = '1e2  '; REPNA(38) = '2e2  '; REPNA(39) = '1e3  '; REPNA(40) = '2e3  '
            REPNA(41) = '1e4  '; REPNA(42) = '2e4  '; REPNA(43) = '1e5  '; REPNA(44) = '2e5  '
            REPNA(45) = '1e6  '; REPNA(46) = '2e6  '; REPNA(47) = '1e7  '; REPNA(48) = '2e7  '
            REPNA(49) = 'a    '; REPNA(50) = 'b    '; REPNA(51) = '1e1  '; REPNA(52) = '2e1  '
            REPNA(53) = '1e2  '; REPNA(54) = '2e2  '; REPNA(55) = '1e3  '; REPNA(56) = '2e3  '
            REPNA(57) = '1e4  '; REPNA(58) = '2e4  '; REPNA(59) = '1e5  '; REPNA(60) = '2e5  '
            REPNA(61) = '1e7  '; REPNA(62) = '2e7  '; REPNA(63) = '1e9  '; REPNA(64) = '2e9  '

            Do i = 1, nsymrpa
                Do j = 1, nsymrpa
                    SD(i, j) = MULTB(i + nsymrpa, j)
                End do
            End do

        Else

            Do i = 1, nsymrpa
                Do j = 1, nsymrpa
                    SD(i, j) = MULTB(i + nsymrpa, j)
                End do
            End do

        End if
    end if

    If (trim(ptgrp) /= 'C1') nsymrp = nsymrpa

    if (rank == 0) then
        print *, 'MULTB_SD'
        Do i = 1, nsymrpa
            print '(50I3)', (SD(i, j), j=1, nsymrpa)
        End do
    end if
    Do i = 1, nsymrpa
        Do j = 1, nsymrpa
            DS(i, j) = SD(j, i)
        End do
    End do

    if (rank == 0) then
        print *, 'MULTB_DS'
        Do i = 1, nsymrpa
            print '(50I3)', (DS(i, j), j=1, nsymrpa)
        End do
    end if
    MULTB_DS(:, :) = DS(:, :)

    deallocate (DS, SD)

    Allocate (irpamo(nmo)); Call memplus(KIND(irpamo), SIZE(irpamo), 1)
    Allocate (orb(nmo)); Call memplus(KIND(orb), SIZE(orb), 1)
    Allocate (indmo(nmo)); Call memplus(KIND(indmo), SIZE(indmo), 1)
    Allocate (indmor(nmo)); Call memplus(KIND(indmor), SIZE(indmor), 1)
    Allocate (dammo(nmo)); Call memplus(KIND(dammo), SIZE(dammo), 1)

!Iwamuro modify
    irpmo(:) = 0
    irpamo(:) = 0

    orbmo(:) = 0.0d+00
    orb(:) = 0.0d+00
    indmo(:) = 0

    Read (mrconee_unit, iostat=iostat) (IRPMO(IMO), IRPAMO(IMO), ORBMO(IMO), IMO=1, NMO)                             ! orbital energies <= used here
    if (iostat .ne. 0) then
        print *, 'Error in reading orbital energies'
        print *, 'iostat = ', iostat
        stop
    end if
    CLOSE (mrconee_unit)

!Iwamuro modify
    irpmo(:) = irpamo(:)

    if (rank == 0) then
        print '("irpamo ",20I2)', (irpamo(i0), i0=1, nmo)
    end if

    orb = orbmo

! orb is lower order of orbmo
    call heapSort(orb, .false.)
    allocate (sort_orb(nmo)); Call memplus(KIND(sort_orb), SIZE(sort_orb), 1)
    sort_orb = orb
! RAS sort
    if (is_ras1_configured .or. is_ras2_configured .or. is_ras3_configured) then
        call sort_list_energy_order_to_ras_order(sort_orb, orb)
    end if

!! orb_sort is lower order of orbmo

    ! sort_orb(i0) and sort_orb(i0+1) should be same orbital energy (kramers pair)
    do i0 = 1, nmo, 2
        m = 0
        do j0 = 1, nmo
            if (orbmo(j0) == sort_orb(i0)) then  ! orbmo(j0) is i0 th MO
                if (m == 0) then
                    indmo(i0) = j0
                    m = m + 1
                else
                    indmo(i0 + 1) = j0
                end if
            end if
        end do
    end do

    if (rank == 0) then
        print *, 'orb sort end'

        print *, 'i0,orb(i0),sort_orb(i0)'
        do i0 = 1, nmo
            print *, i0, orb(i0), sort_orb(i0)
        end do
    end if

    do i0 = 1, nmo
        indmor(indmo(i0)) = i0  ! i0 is energetic order, indmo(i0) is symmtric order (MRCONEE order)
    end do

    if (rank == 0) then
        do i0 = 1, nmo
            print '("indmor output",3I4)', indmor(i0), indmo(i0), i0
        end do
    end if
    orbmo = sort_orb

    dammo = irpmo

    do i0 = 1, nmo
        irpmo(i0) = dammo(indmo(i0))
        irpamo(i0) = dammo(indmo(i0))
    end do

    if (rank == 0) then
        print '("irpamo ",20I2)', (irpamo(i0), i0=1, nmo)

        print *, 'inactive'
        do i0 = 1, ninact
            print '(2I4,2X,E20.10,2X,I4)', i0, indmo(i0), orbmo(i0), irpmo(i0)
        end do

        print *, 'active'
        do i0 = ninact + 1, ninact + nact
            print '(2I4,2X,E20.10,2X,I4)', i0, indmo(i0), orbmo(i0), irpmo(i0)
        end do

        print *, 'secondary'
        do i0 = ninact + nact + 1, ninact + nact + nsec
            print '(2I4,2X,E20.10,2X,I4)', i0, indmo(i0), orbmo(i0), irpmo(i0)
        end do
    end if

    if (allocated(dammo)) deallocate (dammo); Call memminus(KIND(dammo), SIZE(dammo), 1)
    if (allocated(orb)) deallocate (orb); Call memminus(KIND(orb), SIZE(orb), 1)
    if (allocated(sort_orb)) deallocate (sort_orb); Call memminus(KIND(sort_orb), SIZE(sort_orb), 1)
contains

    subroutine sort_list_energy_order_to_ras_order(want_to_sort, original_orb_energy_order)
        use four_caspt2_module, only: ras1_list, ras2_list, ras3_list, ninact, nact, nsec, ras1_size, ras2_size, ras3_size
        implicit none
        real(8), intent(in) :: original_orb_energy_order(:)
        real(8), intent(inout) :: want_to_sort(:)
        integer :: current_spinor_idx, current_idx, idx
        integer :: ras1_current_idx, ras2_current_idx, ras3_current_idx
        if (rank == 0) print *, 'sizeofras', ras1_size, ras2_size, ras3_size
        if (ras1_size == 0 .and. ras2_size == 0 .and. ras3_size == 0) return ! Do nothing because ras is not configured
        current_spinor_idx = 1; current_idx = 1; ras1_current_idx = 1; ras2_current_idx = 1; ras3_current_idx = 1 ! Initialization
        ! Fill ninact
        do while (current_idx <= ninact)
            if (is_ras1_configured .and. ras1_list(ras1_current_idx) == current_spinor_idx) then
                if (ras1_size > ras1_current_idx) ras1_current_idx = ras1_current_idx + 1 ! Skip ras1_list(ras1_current_idx)
            elseif (is_ras2_configured .and. ras2_list(ras2_current_idx) == current_spinor_idx) then
                if (ras2_size > ras2_current_idx) ras2_current_idx = ras2_current_idx + 1 ! Skip ras2_list(ras2_current_idx)
            elseif (is_ras3_configured .and. ras3_list(ras3_current_idx) == current_spinor_idx) then
                if (ras3_size > ras3_current_idx) ras3_current_idx = ras3_current_idx + 1 ! Skip ras3_list(ras3_current_idx)
            else
                want_to_sort(current_idx) = original_orb_energy_order(current_spinor_idx)
                current_idx = current_idx + 1
            end if
            current_spinor_idx = current_spinor_idx + 1 ! Next spinor (energy order)
        end do
        ! current_idx must be ninact + 1
        if (current_idx /= ninact + 1) then
            print *, "ERROR: Sorting energy ascending order to RAS order is failed... STOP THE PROGRAM"
            print *, "ORIGINAL ENERGY ORDER LIST : ", original_orb_energy_order
            print *, "LIST OF SORTING IN PROGRESS: ", want_to_sort(1:current_idx)
            stop ! ERROR, STOP THE PROGRAM
        end if

        ! Fill active
        ! Fill ras1
        if (ras1_size > 0) then
            do idx = 1, ras1_size
                want_to_sort(current_idx + idx - 1) = original_orb_energy_order(ras1_list(idx))
            end do
            current_idx = current_idx + ras1_size
        end if
        ! Fill ras2
        if (ras2_size > 0) then
            do idx = 1, ras2_size
                want_to_sort(current_idx + idx - 1) = original_orb_energy_order(ras2_list(idx))
            end do
            current_idx = current_idx + ras2_size
        end if
        ! Fill ras3
        if (ras3_size > 0) then
            do idx = 1, ras3_size
                want_to_sort(current_idx + idx - 1) = original_orb_energy_order(ras3_list(idx))
            end do
            current_idx = current_idx + ras3_size
        end if

        ! current_idx must be ninact + nact +1
        if (current_idx /= ninact + nact + 1) then
            print *, "ERROR: Sorting energy ascending order to RAS order is failed... STOP THE PROGRAM"
            print *, "ORIGINAL ENERGY ORDER LIST : ", original_orb_energy_order
            print *, "LIST OF SORTING IN PROGRESS: ", want_to_sort(1:current_idx)
            stop ! ERROR, STOP THE PROGRAM
        end if
        ! Fill secondary
        do while (current_idx <= ninact + nact + nsec)
            if (ras1_size > 0 .and. ras1_list(ras1_current_idx) == current_spinor_idx) then
                if (ras1_size > ras1_current_idx) ras1_current_idx = ras1_current_idx + 1 ! Skip ras1_list(ras1_current_idx)
            elseif (ras2_size > 0 .and. ras2_list(ras2_current_idx) == current_spinor_idx) then
                if (ras2_size > ras2_current_idx) ras2_current_idx = ras2_current_idx + 1 ! Skip ras2_list(ras2_current_idx)
            elseif (ras3_size > 0 .and. ras3_list(ras3_current_idx) == current_spinor_idx) then
                if (ras3_size > ras3_current_idx) ras3_current_idx = ras3_current_idx + 1 ! Skip ras3_list(ras3_current_idx)
            else
                want_to_sort(current_idx) = original_orb_energy_order(current_spinor_idx)
                current_idx = current_idx + 1
            end if
            current_spinor_idx = current_spinor_idx + 1 ! Next spinor (energy order)
        end do
    end subroutine sort_list_energy_order_to_ras_order
end subroutine readorb_enesym_co
