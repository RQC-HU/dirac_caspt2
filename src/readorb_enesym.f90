! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE readorb_enesym(filename)

! This subroutine reads the molecular orbitals energies, symmetry information
! and NMO, BREIT, ECORE from MRCONEE file.
! - MO energies are sorted in ascending order.
! - The symmetry information is stored in the array IRPAMO.
! - Creates the multiplication tables for the irreducible representations.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module
    use module_error, only: stop_with_errorcode
    use module_file_manager
    use module_sort_swap
    Implicit NONE

    integer :: unit_mrconee, IMO, IRP
    integer :: nsymrp
    character :: repn(64)*14
    integer, allocatable :: irpmo(:)
    character*50, intent(in) :: filename
    integer :: i0, j0, k0, i, j, m, iostat
    integer, allocatable :: dammo(:)
    integer, allocatable :: SD(:, :)
    logical :: breit, is_end_of_file

    call open_unformatted_file(unit=unit_mrconee, file=trim(filename), status='old', optional_action='read')

    ! Read the number of molecular orbitals, Breit interaction and the core energy.
    Read (unit_mrconee, iostat=iostat) NMO, BREIT, ECORE  ! NMO is nbas - ncore
    call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
    if (is_end_of_file) then
        print *, 'Error: error in reading NMO, BREIT, ECORE (end of file reached)'
        print *, 'iostat = ', iostat
        call stop_with_errorcode(iostat)
    end if
    if (rank == 0) then
        print *, 'NMO, BREIT, ECORE, 1  ! NMO is nbas - ncore'
        print *, NMO, BREIT, ECORE, 1  ! NMO is nbas - ncore
    end if

    ! Read the number of irreducible representations and irreducible representation labels for each molecular orbital.
    Read (unit_mrconee, iostat=iostat) NSYMRP, (REPN(IRP), IRP=1, NSYMRP)
    call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
    if (is_end_of_file) then
        print *, 'Error: error in reading NSYMRP, REPN (end of file reached)'
        print *, 'iostat = ', iostat
        call stop_with_errorcode(iostat)
    end if
    if (rank == 0) then
        print *, ' NSYMRP, (REPN(IRP),IRP=1,NSYMRP)'
        print *, NSYMRP, (REPN(IRP), IRP=1, NSYMRP)
    end if

    ! Read the number of irreducible representations and irreducible representation labels for each molecular orbital. (abelian group)
    Read (unit_mrconee, iostat=iostat) nsymrpa, (repna(i0), i0=1, nsymrpa*2)
    call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
    if (is_end_of_file) then
        print *, 'Error: error in reading nsymrpa, repna (end of file reached)'
        print *, 'iostat = ', iostat
        call stop_with_errorcode(iostat)
    end if
    if (rank == 0) then
        print *, ' NSYMRPA, (REPNA(IRP),IRP=1,NSYMRPA*2)'
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
    allocate (SD(1:NSYMRPA, 1:NSYMRPA))
    Call memplus(size(SD), kind(SD), 1)

    ! Read the multiplication table for the irreducible representations.
    Read (unit_mrconee, iostat=iostat) ((multb(i0, j0), i0=1, 2*nsymrpa), j0=1, 2*nsymrpa)
    call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
    if (is_end_of_file) then
        print *, 'Error: error in reading multb (end of file reached)'
        print *, 'iostat = ', iostat
        call stop_with_errorcode(iostat)
    end if

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

    ! create MULTB_S, MULTB_D and MULTB_DS
    ! MULTB_S(i,j) = MULTB(i + nsymrpa, j + nsymrpa) - nsymrpa (the bottom right block of MULTB)
    MULTB_S(:, :) = MULTB(1 + nsymrpa:2*nsymrpa, 1 + nsymrpa:2*nsymrpa) - nsymrpa
    ! MULTB_D is the upper left block of MULTB2 - nsymrpa
    MULTB_D(:, :) = MULTB2(1:nsymrpa, 1:nsymrpa) - nsymrpa
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
    SD(:, :) = MULTB(nsymrpa + 1:2*nsymrpa, 1:nsymrpa)
    if (rank == 0) then
        print *, 'MULTB_SD'
        Do i = 1, nsymrpa
            print '(50I3)', (SD(i, j), j=1, nsymrpa)
        End do
    end if
    MULTB_DS = transpose(SD)
    if (rank == 0) then
        print *, 'MULTB_DS'
        Do i = 1, nsymrpa
            print '(50I3)', (MULTB_DS(i, j), j=1, nsymrpa)
        End do
    end if
    if (allocated(SD)) Call memminus(KIND(SD), SIZE(SD), 1); deallocate (SD)

    ! Define the space index for each molecular orbital.
    Allocate (space_idx(1:nmo)); Call memplus(KIND(space_idx), SIZE(space_idx), 1)
    space_idx(1:ninact) = 1 ! inactive = 1
    space_idx(ninact + 1:ninact + nact) = 2 ! active = 2
    space_idx(ninact + nact + 1:ninact + nact + nsec) = 3 ! secondary = 3
    space_idx(ninact + nact + nsec + 1:nmo) = 4 ! virtual = 4

    allocate (IRPMO(1:NMO)); call memplus(size(IRPMO), kind(IRPMO), 1)
    Allocate (irpamo(nmo)); Call memplus(KIND(irpamo), SIZE(irpamo), 1)
    allocate (dirac_mo_energy(1:NMO)); call memplus(size(dirac_mo_energy), kind(dirac_mo_energy), 1)
    allocate (caspt2_mo_energy(1:NMO)); call memplus(size(caspt2_mo_energy), kind(caspt2_mo_energy), 1)
    Allocate (indmo_cas_to_dirac(nmo)); Call memplus(KIND(indmo_cas_to_dirac), SIZE(indmo_cas_to_dirac), 1)
    Allocate (indmo_dirac_to_cas(nmo)); Call memplus(KIND(indmo_dirac_to_cas), SIZE(indmo_dirac_to_cas), 1)
    Allocate (dammo(nmo)); Call memplus(KIND(dammo), SIZE(dammo), 1)
    irpmo(:) = 0
    irpamo(:) = 0
    indmo_cas_to_dirac(:) = 0
    ! Read the irpmo, irpamo, dirac_mo_enegy.
    ! irpmo: irreducible representation number of each molecular orbital.
    ! irpamo: irreducible representation number of each molecular orbital. (abelian group)
    ! dirac_mo_energy: orbital energies of each molecular orbital.
    Read (unit_mrconee, iostat=iostat) (IRPMO(IMO), IRPAMO(IMO), dirac_mo_energy(IMO), IMO=1, NMO)                             ! orbital energies <= used here
    if (iostat .ne. 0) then
        print *, 'Error in reading orbital energies'
        print *, 'iostat = ', iostat
        call stop_with_errorcode(iostat)
    end if
    CLOSE (unit_mrconee)
    ! Print irpmo, irpamo
    if (rank == 0) then
        print '("irpmo ",20I3)', (irpmo(i0), i0=1, nmo)
    end if
    if (allocated(irpmo)) deallocate (irpmo); Call memminus(KIND(irpmo), SIZE(irpmo), 1)
    if (rank == 0) then
        print '("irpamo ",20I3)', (irpamo(i0), i0=1, nmo)
    end if

    ! Sort the orbital energies in ascending order.
    caspt2_mo_energy = dirac_mo_energy
    call heapSort(list=caspt2_mo_energy, is_reverse=.false.)

    ! RAS sort (if RAS is used)
    if (ras1_size /= 0 .or. ras2_size /= 0 .or. ras3_size /= 0) then
        call sort_list_from_energy_order_to_ras_order(caspt2_mo_energy)
    end if

    ! caspt2_mo_energy(i0) and caspt2_mo_energy(i0+1) should be same orbital energy (kramers pair)
    do i0 = 1, nmo, 2
        m = 0
        do j0 = 1, nmo
            if (dirac_mo_energy(j0) == caspt2_mo_energy(i0)) then  ! dirac_mo_energy(j0) is i0 th MO
                if (m == 0) then
                    indmo_cas_to_dirac(i0) = j0
                    m = m + 1
                else
                    indmo_cas_to_dirac(i0 + 1) = j0
                end if
            end if
        end do
    end do

    do i0 = 1, nmo
        indmo_dirac_to_cas(indmo_cas_to_dirac(i0)) = i0  ! i0 is energetic order, indmo_cas_to_dirac(i0) is symmtric order (MRCONEE order)
    end do

    if (rank == 0) then
        do i0 = 1, nmo
            print '("indmo_dirac_to_cas output",3I4)', indmo_dirac_to_cas(i0), indmo_cas_to_dirac(i0), i0
        end do
    end if

    ! irpamo is in MRCONEE order (DIRAC order)
    dammo = irpamo

    ! Convert irpamo and irpamo into energy order (CAS order)
    do i0 = 1, nmo
        irpamo(i0) = dammo(indmo_cas_to_dirac(i0))
        irpamo(i0) = dammo(indmo_cas_to_dirac(i0))
    end do

    if (rank == 0) then
        print '("irpamo ",20I3)', (irpamo(i0), i0=1, nmo)

        print *, 'inactive'
        do i0 = 1, ninact
            print '(2I4,2X,E20.10,2X,I4,1X,A)', i0, indmo_cas_to_dirac(i0), caspt2_mo_energy(i0), irpamo(i0), repna(irpamo(i0))
        end do

        print *, 'active'
        do i0 = ninact + 1, ninact + nact
            print '(2I4,2X,E20.10,2X,I4,1X,A)', i0, indmo_cas_to_dirac(i0), caspt2_mo_energy(i0), irpamo(i0), repna(irpamo(i0))
        end do

        print *, 'secondary'
        do i0 = ninact + nact + 1, ninact + nact + nsec
            print '(2I4,2X,E20.10,2X,I4,1X,A)', i0, indmo_cas_to_dirac(i0), caspt2_mo_energy(i0), irpamo(i0), repna(irpamo(i0))
        end do
    end if

    if (allocated(dammo)) deallocate (dammo); Call memminus(KIND(dammo), SIZE(dammo), 1)
contains

    subroutine sort_list_from_energy_order_to_ras_order(want_to_sort)
        !===========================================================================================================================
        ! This subroutine sorts the want_to_sort list form orbital energy order
        ! to RAS order(ninact => ras1 => ras2 => ras3 => secondary).
        !===========================================================================================================================
        use four_caspt2_module, only: ras1_list, ras2_list, ras3_list, ninact, nact, nsec, ras1_size, ras2_size, ras3_size
        implicit none
        real(8), intent(inout) :: want_to_sort(:)
        real(8), allocatable :: mo_energy_order(:)
        integer :: idx_energy_order, idx_ras_order, idx
        integer :: ras1_idx, ras2_idx, ras3_idx
        if (rank == 0) print *, 'sizeofras', ras1_size, ras2_size, ras3_size
        if (ras1_size == 0 .and. ras2_size == 0 .and. ras3_size == 0) return ! Do nothing because ras is not configured
        ! Initialization
        idx_energy_order = 1; idx_ras_order = 1
        ras1_idx = 1; ras2_idx = 1; ras3_idx = 1
        allocate (mo_energy_order(size(want_to_sort)))
        mo_energy_order = want_to_sort ! Save the original orbital energy order
        ! Fill ninact
        do while (idx_ras_order <= ninact)
            if (ras1_size /= 0 .and. ras1_list(ras1_idx) == idx_energy_order) then
                if (ras1_size > ras1_idx) ras1_idx = ras1_idx + 1 ! Skip ras1_list(ras1_idx)
            elseif (ras2_size /= 0 .and. ras2_list(ras2_idx) == idx_energy_order) then
                if (ras2_size > ras2_idx) ras2_idx = ras2_idx + 1 ! Skip ras2_list(ras2_idx)
            elseif (ras3_size /= 0 .and. ras3_list(ras3_idx) == idx_energy_order) then
                if (ras3_size > ras3_idx) ras3_idx = ras3_idx + 1 ! Skip ras3_list(ras3_idx)
            else
                want_to_sort(idx_ras_order) = mo_energy_order(idx_energy_order)
                idx_ras_order = idx_ras_order + 1
            end if
            idx_energy_order = idx_energy_order + 1 ! Next spinor (energy order)
        end do
        ! idx_ras_order must be ninact + 1
        if (idx_ras_order /= ninact + 1) then
            print *, "ERROR: Sorting energy ascending order to RAS order is failed... STOP THE PROGRAM"
            print *, "ORIGINAL ENERGY ORDER LIST : ", mo_energy_order
            print *, "LIST OF SORTING IN PROGRESS: ", want_to_sort(1:idx_ras_order)
            call stop_with_errorcode(1) ! ERROR, STOP THE PROGRAM
        end if

        ! Fill active
        ! Fill ras1
        if (ras1_size > 0) then
            do idx = 1, ras1_size
                want_to_sort(idx_ras_order + idx - 1) = mo_energy_order(ras1_list(idx))
            end do
            idx_ras_order = idx_ras_order + ras1_size
        end if
        ! Fill ras2
        if (ras2_size > 0) then
            do idx = 1, ras2_size
                want_to_sort(idx_ras_order + idx - 1) = mo_energy_order(ras2_list(idx))
            end do
            idx_ras_order = idx_ras_order + ras2_size
        end if
        ! Fill ras3
        if (ras3_size > 0) then
            do idx = 1, ras3_size
                want_to_sort(idx_ras_order + idx - 1) = mo_energy_order(ras3_list(idx))
            end do
            idx_ras_order = idx_ras_order + ras3_size
        end if

        ! idx_ras_order must be ninact + nact +1
        if (idx_ras_order /= ninact + nact + 1) then
            print *, "ERROR: Sorting energy ascending order to RAS order is failed... STOP THE PROGRAM"
            print *, "ORIGINAL ENERGY ORDER LIST : ", mo_energy_order
            print *, "LIST OF SORTING IN PROGRESS: ", want_to_sort(1:idx_ras_order)
            call stop_with_errorcode(1) ! ERROR, STOP THE PROGRAM
        end if
        ! Fill secondary
        do while (idx_ras_order <= ninact + nact + nsec)
            if (ras1_size > 0 .and. ras1_list(ras1_idx) == idx_energy_order) then
                if (ras1_size > ras1_idx) ras1_idx = ras1_idx + 1 ! Skip ras1_list(ras1_idx)
            elseif (ras2_size > 0 .and. ras2_list(ras2_idx) == idx_energy_order) then
                if (ras2_size > ras2_idx) ras2_idx = ras2_idx + 1 ! Skip ras2_list(ras2_idx)
            elseif (ras3_size > 0 .and. ras3_list(ras3_idx) == idx_energy_order) then
                if (ras3_size > ras3_idx) ras3_idx = ras3_idx + 1 ! Skip ras3_list(ras3_idx)
            else
                want_to_sort(idx_ras_order) = mo_energy_order(idx_energy_order)
                idx_ras_order = idx_ras_order + 1
            end if
            idx_energy_order = idx_energy_order + 1 ! Next spinor (energy order)
        end do
    end subroutine sort_list_from_energy_order_to_ras_order
end subroutine readorb_enesym