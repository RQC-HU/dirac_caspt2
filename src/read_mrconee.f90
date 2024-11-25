subroutine check_dirac_integer_size(filename)
    use, intrinsic :: iso_fortran_env, only: int32
    use module_global_variables
    use module_error, only: stop_with_errorcode
    use module_file_manager
    use module_sort_swap
    Implicit NONE

    integer :: unit_mrconee
    character(*), intent(in) :: filename
    integer :: iostat
    logical :: is_end_of_file
    integer :: nsymrp
    integer(kind=int32) :: leading_rec_marker, trailing_rec_marker, record_size, record_offset, seek_st
    nsymrp = 0
    record_offset = 0
    seek_st = 0

    ! Get an unused unit
    print *, 'filename = ', filename
    call open_unformatted_file(unit=unit_mrconee, file=trim(filename), status='old', optional_action='read')
    rewind (unit_mrconee)
    close (unit_mrconee)

    ! Open MRCONEE (access="stream")
    open (unit_mrconee, file=trim(filename), form="unformatted", status="old", access="stream")

    ! =============================================================================================
    ! SKIP: 1st line (the number of molecular orbitals, Breit interaction and the core energy...)
    ! read (unit_mrconee, iostat=iostat) NMO, BREIT, ECORE, nfsym, spinfr, norbt, hf_energy_mrconee
    ! =============================================================================================

    ! 1st line leading record_marker
    read (unit_mrconee, iostat=iostat) leading_rec_marker
    call check_iostat(iostat, trim(filename), is_end_of_file)
    record_offset = leading_rec_marker

    ! seek
    call seek(unit_mrconee, record_offset, seek_st)
    call check_fseek_status(status=seek_st)

    ! 1st line trailing record_marker
    read (unit_mrconee, iostat=iostat) trailing_rec_marker
    call check_iostat(iostat, trim(filename), is_end_of_file)
    if (leading_rec_marker /= trailing_rec_marker) then
        if (rank == 0) print *, "Values of leading and trailing record markers are not same. Exiting..."
        call stop_with_errorcode(1)
    end if

    ! =============================================================================================
    ! SKIP: 2nd line (nsymrp and irpmo)
    ! read (unit_mrconee, iostat=iostat) nsymrp, (REPN(IRP),IRP=1,NSYMRP)
    ! =============================================================================================

    ! 2nd line leading record_marker
    read (unit_mrconee, iostat=iostat) leading_rec_marker
    call check_iostat(iostat, trim(filename), is_end_of_file)
    record_offset = leading_rec_marker

    ! seek
    call seek(unit_mrconee, record_offset, seek_st)
    call check_fseek_status(status=seek_st)

    ! 2nd line trailing record_marker
    read (unit_mrconee, iostat=iostat) trailing_rec_marker
    call check_iostat(iostat, trim(filename), is_end_of_file)
    if (leading_rec_marker /= trailing_rec_marker) then
        if (rank == 0) print *, "Values of leading and trailing record markers are not same. Exiting..."
        call stop_with_errorcode(1)
    end if

    ! ========================================================================
    ! 3rd line (nsymrpa and repna)
    ! Read (unit_mrconee, iostat=iostat) nsymrpa, (repna(i0), i0=1, nsymrpa*2)
    ! ========================================================================

    ! 3rd line leading record_marker
    read (unit_mrconee, iostat=iostat) leading_rec_marker
    call check_iostat(iostat, trim(filename), is_end_of_file)

    ! nsymrpa
    Read (unit_mrconee, iostat=iostat) nsymrpa
    ! nsymrpa: 8byte integer if DIRAC is 64bit integer version.
    ! Read (unit_mrconee, iostat=iostat) nsymrpa, (repna(i0), i0=1, nsymrpa*2)
    record_size = kind(nsymrpa) + len(repna(1))*nsymrpa*2
    if (record_size /= leading_rec_marker) then
        call dirac_32bit_integer
        close (unit_mrconee)
        return
    end if

    ! seek
    record_offset = record_size - kind(nsymrpa)
    call seek(unit_mrconee, record_offset, seek_st)
    if (seek_st /= 0) then
        call dirac_32bit_integer
        close (unit_mrconee)
        return
    end if

    ! 3rd line trailing record_marker
    read (unit_mrconee, iostat=iostat) trailing_rec_marker
    call check_iostat(iostat, trim(filename), is_end_of_file)
    if (leading_rec_marker /= trailing_rec_marker) then
        call dirac_32bit_integer
        close (unit_mrconee)
        return
    end if

    ! ======================================================================================
    ! 4th line (multiplication table)
    ! multiplication table for the irreducible representations.
    ! Read (unit_mrconee, iostat=iostat) ((multb(i0, j0), i0=1, 2*nsymrpa), j0=1, 2*nsymrpa)
    ! ======================================================================================

    ! 4th line leading record_marker
    read (unit_mrconee, iostat=iostat) leading_rec_marker
    call check_iostat(iostat, trim(filename), is_end_of_file)

    ! MULTB
    ! Read (unit_mrconee, iostat=iostat) ((multb(i0, j0), i0=1, 2*nsymrpa), j0=1, 2*nsymrpa)
    record_size = kind(multb)*(2*nsymrpa)**2
    if (record_size /= leading_rec_marker) then
        call dirac_32bit_integer
        close (unit_mrconee)
        return
    end if

    ! seek
    record_offset = record_size
    call seek(unit_mrconee, record_offset, seek_st)
    if (seek_st > 0) then
        call dirac_32bit_integer
        close (unit_mrconee)
        return
    end if

    ! 4th line trailing record_marker
    read (unit_mrconee, iostat=iostat) trailing_rec_marker
    call check_iostat(iostat, trim(filename), is_end_of_file)
    if (leading_rec_marker /= trailing_rec_marker) then
        call dirac_32bit_integer
        close (unit_mrconee)
        return
    end if

    ! ===================================================================
    ! SKIP: 5th line (MO info)
    ! WRITE (LUMLF1) (IRPMO(I),IRPAMO(I),ORBMO(I),I=1,NSTRT*2),
    ! &               (IBSPI(I),I=1,2*NSTRT),(NORB(I),I=1,NFSYM),NBSYMRP
    ! ===================================================================

    ! 5th line leading record_marker
    read (unit_mrconee, iostat=iostat) leading_rec_marker
    call check_iostat(iostat, trim(filename), is_end_of_file)
    record_offset = leading_rec_marker

    ! seek
    call seek(unit_mrconee, record_offset, seek_st)
    call check_fseek_status(status=seek_st)

    ! 5th line trailing record_marker
    read (unit_mrconee, iostat=iostat) trailing_rec_marker
    call check_iostat(iostat, trim(filename), is_end_of_file)
    if (leading_rec_marker /= trailing_rec_marker) then
        call dirac_32bit_integer
        close (unit_mrconee)
        return
    end if

    ! ===================================================================
    ! SKIP: 6th line (1-elec integrals)
    ! WRITE (LUMLF1) ((FMOM(I,J,1),FMOM(I,J,2),
    !  &                 I=1,NSTRT*2),J=1,NSTRT*2)
    ! ===================================================================

    ! 6th line leading record_marker
    read (unit_mrconee, iostat=iostat) leading_rec_marker
    call check_iostat(iostat, trim(filename), is_end_of_file)
    record_offset = leading_rec_marker

    ! seek
    call seek(unit_mrconee, record_offset, seek_st)
    call check_fseek_status(status=seek_st)

    ! 6th line trailing record_marker
    read (unit_mrconee, iostat=iostat) trailing_rec_marker
    call check_iostat(iostat, trim(filename), is_end_of_file)
    if (leading_rec_marker /= trailing_rec_marker) then
        call dirac_32bit_integer
        close (unit_mrconee)
        return
    end if

    if (rank == 0) print *, "DIRAC seems to be compiled with 64-bit integer."
    close (unit_mrconee)
contains

    subroutine seek(unit, offset, seek_status)
        implicit none
        integer, intent(in) :: unit
        integer(kind=4), intent(in) :: offset
        integer(kind=4), intent(out) :: seek_status
        call fseek(unit, offset, 1, seek_status) ! 1: seek from current position
    end subroutine seek

    subroutine check_fseek_status(status)
        implicit none
        integer(kind=4) :: status

        if (status /= 0) then
            if (rank == 0) print *, "fseek failed. exiting..."
            call stop_with_errorcode(1)
        end if
    end subroutine check_fseek_status

    subroutine dirac_32bit_integer
        implicit none
        if (rank == 0) print *, "DIRAC seems to be compiled with 32-bit integer."
        dirac_32bit_build = .true.
    end subroutine dirac_32bit_integer
end subroutine check_dirac_integer_size

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE read_mrconee(filename)

! This subroutine reads DIRAC MRCONEE file.
! Data in MRCONEE
! - molecular orbitals energies (dirac_mo_energy)
! - symmetry information (irpamo, nsymrpa, repna)
! - NMO, BREIT, ECORE
! - 1 electron integrals (roner, ronei)
! After this subroutine is called, the following information is available:
! - MO energies are sorted in ascending order (caspt2_mo_energy)
! - The symmetry information is stored in the array IRPAMO.
! - Creates the multiplication tables for the irreducible representations.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use, intrinsic :: iso_fortran_env, only: int32
    use module_global_variables
    use module_error, only: stop_with_errorcode
    use module_file_manager
    use module_sort_swap
    Implicit NONE
    character(*), intent(in) :: filename

    integer :: unit_mrconee
    integer(kind=int32) :: nmo_32bit, nfsym_32bit, nz_32bit, norbt_32bit
    logical(kind=int32) :: breit_32bit, spinfr_32bit
    integer(kind=int64) :: IMO, IRP
    integer(kind=int64) :: i0, j0, k0, i, j, m
    logical(kind=int64) :: breit, spinfr
    integer(kind=int64) :: nfsym, nz, norbt
    integer :: iostat
    logical :: is_end_of_file
    call open_unformatted_file(unit=unit_mrconee, file=trim(filename), status='old', optional_action='read')

    ! Read the number of molecular orbitals, Breit interaction and the core energy and HF energy.
    ! reference: https://gitlab.com/dirac/dirac/-/blob/01878d230962146d8183020b51a97e12b080de99/src/moltra/traone.F#L815-816
    ! NSPC and NCORE2 added in DIRAC 21 and we don't use them, so skip them.
    ! NSPC and NCORE2 were added at the following commit.
    ! https://gitlab.com/dirac/dirac/-/commit/d0a1beb1fc0c23b4ce89c89c05bc0d421f71aea2
    if (dirac_32bit_build) then
        read (unit_mrconee, iostat=iostat) nmo_32bit, breit_32bit, ecore, &
            nfsym_32bit, nz_32bit, spinfr_32bit, norbt_32bit, hf_energy_mrconee
        nmo = nmo_32bit; breit = breit_32bit; nfsym = nfsym_32bit; nz = nz_32bit; spinfr = spinfr_32bit; norbt = norbt_32bit
    else
        read (unit_mrconee, iostat=iostat) NMO, BREIT, ECORE, nfsym, nz, spinfr, norbt, hf_energy_mrconee
    end if
    call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)
    if (is_end_of_file) then
        print *, 'Error: error in reading NMO, BREIT, ECORE (end of file reached)'
        print *, 'iostat = ', iostat
        call stop_with_errorcode(iostat)
    end if
    if (rank == 0) then
        print *, 'Information from MRCONEE'
        print *, 'NMO, BREIT, ECORE'
        print *, NMO, BREIT, ECORE
    end if

! Read the irreducible representation information. (nsymrpa, repna)
    call read_irreducible_representation_infomation

! Read DIRAC original multiplication table and convert it to the CASPT2 multiplication tables.
! MULTB -> MULTB2, MULTB_D, MULTB_S, MULTB_DS
    call create_multiplication_table

! Create mo_energy list and irreducible representation mapping list.
    call create_mo_irrep_conversion_list

!   1   1     -0.1336324989E+02     1   1g
!   2 232     -0.1336324989E+02     2  -1g
    if (rank == 0) then
        if (debug) print '("irpamo ",20I3)', (irpamo(i0), i0=1, nmo)
        print *, ' '
        print '(64A)', '----------------------------------------------------------------'
        print '(64A)', '        energy-order    Dirac     orbtal energy    irrep  irrep '
        print '(64A)', '           index        index         (a.u.)       index  string'
        print '(64A)', '----------------------------------------------------------------'
!        print '(64A)', ' inact      2000         4000   -0.1336324989E+02     20   -1g  '
        do i0 = 1, ninact
            print '(X,"inactive",I7,I13,E20.10,I7,4X,4A)', &
            &  i0, indmo_cas_to_dirac(i0), caspt2_mo_energy(i0), irpamo(i0), repna(irpamo(i0))
        end do
        do i0 = global_act_start, global_act_end
            print '(X,"active  ",I7,I13,E20.10,I7,4X,4A)', &
            &  i0, indmo_cas_to_dirac(i0), caspt2_mo_energy(i0), irpamo(i0), repna(irpamo(i0))
        end do
        do i0 = global_sec_start, global_sec_end
            print '(X,"secondary ",I5,I13,E20.10,I7,4X,4A)', &
            &  i0, indmo_cas_to_dirac(i0), caspt2_mo_energy(i0), irpamo(i0), repna(irpamo(i0))
        end do
        print *, ' '
    end if

! Read 1 electron integrals to the variables one_elec_int_r and one_elec_int_i
    call read_1_elec_integrals

    close (unit_mrconee)

contains

    subroutine read_irreducible_representation_infomation
        implicit none
        integer(kind=int64) :: nsymrp
        integer(kind=int32) :: nsymrp_32bit, nsymrpa_32bit
        character :: repn(64)*14
        ! Read the number of irreducible representations and irreducible representation labels for each molecular orbital.
        if (dirac_32bit_build) then
            read (unit_mrconee, iostat=iostat) nsymrp_32bit, (REPN(IRP), IRP=1, nsymrp_32bit)
            nsymrp = nsymrp_32bit
        else
            read (unit_mrconee, iostat=iostat) NSYMRP, (REPN(IRP), IRP=1, NSYMRP)
        end if
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
        if (dirac_32bit_build) then
            read (unit_mrconee, iostat=iostat) nsymrpa_32bit, (repna(i0), i0=1, nsymrpa_32bit*2)
            nsymrpa = nsymrpa_32bit
        else
            read (unit_mrconee, iostat=iostat) nsymrpa, (repna(i0), i0=1, nsymrpa*2)
        end if
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
    end subroutine read_irreducible_representation_infomation

    subroutine create_multiplication_table
        ! create multiplication table
        ! MULTB, MULTB2, MULTB_S, MULTB_D, MULTB_DS
        implicit none
        integer(kind=int32) :: multb_32bit(128, 128)
        integer(kind=int64), allocatable :: SD(:, :)

        allocate (MULTB_S(1:NSYMRPA, 1:NSYMRPA))
        allocate (MULTB_D(1:NSYMRPA, 1:NSYMRPA))  ! dagger
        Call memplus(size(MULTB_S), kind(MULTB_S), 1)
        Call memplus(size(MULTB_D), kind(MULTB_D), 1)
        allocate (SD(1:NSYMRPA, 1:NSYMRPA))
        Call memplus(size(SD), kind(SD), 1)

        ! Read the multiplication table for the irreducible representations.
        if (dirac_32bit_build) then
            read (unit_mrconee, iostat=iostat) ((multb_32bit(i0, j0), i0=1, 2*nsymrpa), j0=1, 2*nsymrpa)
            multb = multb_32bit
        else
            read (unit_mrconee, iostat=iostat) ((multb(i0, j0), i0=1, 2*nsymrpa), j0=1, 2*nsymrpa)
        end if
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
            If (debug) then
                print *, 'MULTB2'
                Do i0 = 1, 2*nsymrpa
                    print '(400I3)', (MULTB2(i0, j0), j0=1, 2*nsymrpa)
                End do
                print *, 'end multb1,2'
            end if
        end if

        ! create MULTB_S, MULTB_D and MULTB_DS
        ! MULTB_S(i,j) = MULTB(i + nsymrpa, j + nsymrpa) - nsymrpa (the bottom right block of MULTB)
        MULTB_S(:, :) = MULTB(1 + nsymrpa:2*nsymrpa, 1 + nsymrpa:2*nsymrpa) - nsymrpa
        ! MULTB_D is the upper left block of MULTB2 - nsymrpa
        MULTB_D(:, :) = MULTB2(1:nsymrpa, 1:nsymrpa) - nsymrpa
        if (debug .and. rank == 0) then
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
        if (debug .and. rank == 0) then
            print *, 'MULTB_SD'
            Do i = 1, nsymrpa
                print '(50I3)', (SD(i, j), j=1, nsymrpa)
            End do
        end if
        if (.not. allocated(MULTB_DS)) then
            allocate (MULTB_DS(1:NSYMRPA, 1:NSYMRPA)); Call memplus(size(MULTB_DS), kind(MULTB_DS), 1)
        end if
        MULTB_DS = transpose(SD)
        if (debug .and. rank == 0) then
            print *, 'MULTB_DS'
            Do i = 1, nsymrpa
                print '(50I3)', (MULTB_DS(i, j), j=1, nsymrpa)
            End do
        end if
        Call memminus(KIND(SD), SIZE(SD), 1); deallocate (SD)
    end subroutine create_multiplication_table

    subroutine create_mo_irrep_conversion_list
        implicit none
        integer(kind=int64), allocatable :: tmp_mo(:)
        integer(kind=int64), allocatable :: irpmo(:)
        integer(kind=int32), allocatable :: irpmo_32bit(:), irpamo_32bit(:)
        ! Define the space index for each molecular orbital.
        Allocate (space_idx(1:nmo)); Call memplus(KIND(space_idx), SIZE(space_idx), 1)
        space_idx(1:ninact) = 1 ! inactive = 1
        space_idx(global_act_start:global_act_end) = 2 ! active = 2
        space_idx(global_sec_start:global_sec_end) = 3 ! secondary = 3
        space_idx(global_sec_end + 1:nmo) = 4 ! virtual = 4

        allocate (irpmo(nmo)); call memplus(size(irpmo), kind(irpmo), 1)
        allocate (irpmo_32bit(nmo)); call memplus(size(irpmo_32bit), kind(irpmo_32bit), 1)
        allocate (irpamo_32bit(nmo)); call memplus(size(irpamo_32bit), kind(irpamo_32bit), 1)
        Allocate (irpamo(nmo)); Call memplus(KIND(irpamo), SIZE(irpamo), 1)
        allocate (dirac_mo_energy(nmo)); call memplus(size(dirac_mo_energy), kind(dirac_mo_energy), 1)
        irpmo(:) = 0
        irpamo(:) = 0
        ! Read the irpmo, irpamo, dirac_mo_enegy.
        ! irpmo: irreducible representation number of each molecular orbital.
        ! irpamo: irreducible representation number of each molecular orbital. (abelian group)
        ! dirac_mo_energy: orbital energies of each molecular orbital.
        if (dirac_32bit_build) then
            read (unit_mrconee, iostat=iostat) (irpmo_32bit(i0), irpamo_32bit(i0), dirac_mo_energy(i0), i0=1, nmo)
            irpmo = irpmo_32bit; irpamo = irpamo_32bit
        else
            read (unit_mrconee, iostat=iostat) (IRPMO(IMO), IRPAMO(IMO), dirac_mo_energy(IMO), IMO=1, NMO)
        end if
        if (iostat .ne. 0) then
            print *, 'Error in reading orbital energies'
            print *, 'iostat = ', iostat
            call stop_with_errorcode(iostat)
        end if
        CLOSE (unit_mrconee)

        ! Print irpmo and irpamo
        if (debug .and. rank == 0) then
            print '("irpmo ",20I3)', (irpmo(i0), i0=1, nmo)
        end if
        if (allocated(irpmo)) Call memminus(KIND(irpmo), SIZE(irpmo), 1); deallocate (irpmo)
        if (debug .and. rank == 0) then
            print '("irpamo ",20I3)', (irpamo(i0), i0=1, nmo)
        end if

        ! Sort the orbital energies in ascending order.
        allocate (caspt2_mo_energy(1:NMO)); call memplus(size(caspt2_mo_energy), kind(caspt2_mo_energy), 1)
        caspt2_mo_energy = dirac_mo_energy
        call heapSort(list=caspt2_mo_energy, is_descending_order=.false.)

        ! RAS sort (if RAS is used)
        if (ras1_size /= 0 .or. ras2_size /= 0 .or. ras3_size /= 0) then
            call sort_list_from_energy_order_to_ras_order(caspt2_mo_energy)
        end if

        ! Create indmo_cas_to_dirac and indmo_dirac_to_cas
        ! caspt2_mo_energy(i0) and caspt2_mo_energy(i0+1) should be same orbital energy (kramers pair)
        Allocate (indmo_cas_to_dirac(nmo)); Call memplus(KIND(indmo_cas_to_dirac), SIZE(indmo_cas_to_dirac), 1)
        Allocate (indmo_dirac_to_cas(nmo)); Call memplus(KIND(indmo_dirac_to_cas), SIZE(indmo_dirac_to_cas), 1)
        indmo_cas_to_dirac(:) = 0; indmo_dirac_to_cas(:) = 0
        do i0 = 1, nmo, 2
            m = 0
            do j0 = 1, nmo
                ! i0 is energetic order, j0 is symmtric order (MRCONEE order)
                if (dirac_mo_energy(j0) == caspt2_mo_energy(i0)) then  ! dirac_mo_energy(j0) is i0 th MO
                    if (m == 0) then
                        indmo_cas_to_dirac(i0) = j0
                        indmo_dirac_to_cas(j0) = i0
                        m = m + 1
                    else
                        indmo_cas_to_dirac(i0 + 1) = j0
                        indmo_dirac_to_cas(j0) = i0 + 1
                    end if
                end if
            end do
        end do

        ! irpamo is in MRCONEE order (DIRAC order)
        Allocate (tmp_mo(nmo)); Call memplus(KIND(tmp_mo), SIZE(tmp_mo), 1)
        tmp_mo = irpamo

        ! Convert irpamo and irpamo into energy order (CAS order)
        do i0 = 1, nmo
            irpamo(i0) = tmp_mo(indmo_cas_to_dirac(i0))
            irpamo(i0) = tmp_mo(indmo_cas_to_dirac(i0))
        end do
        if (allocated(tmp_mo)) Call memminus(KIND(tmp_mo), SIZE(tmp_mo), 1); deallocate (tmp_mo)
    end subroutine create_mo_irrep_conversion_list

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    SUBROUTINE read_1_elec_integrals ! one-electron MO integrals in MRCONEE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        use module_global_variables
        use module_file_manager

        Implicit NONE

        integer(kind=int64) :: isp, nmom
        double precision, allocatable :: roner(:, :, :), ronei(:, :, :)

        if (debug .and. rank == 0) then
            print *, 'Enter read_1_elec_integrals'
        end if

        scfru = 1
        Allocate (roner(nmo, nmo, scfru)); Call memplus(KIND(roner), SIZE(roner), 1)
        Allocate (ronei(nmo, nmo, scfru)); Call memplus(KIND(ronei), SIZE(ronei), 1)

        call open_unformatted_file(unit=unit_mrconee, file=trim(filename), status="old", optional_action="read")

        rewind (unit_mrconee)
        read (unit_mrconee, iostat=iostat)
        read (unit_mrconee, iostat=iostat)
        read (unit_mrconee, iostat=iostat)
        read (unit_mrconee, iostat=iostat)
        read (unit_mrconee, iostat=iostat)
        read (unit_mrconee, iostat=iostat) (((roner(i0, j0, isp), ronei(i0, j0, isp), j0=1, nmo), i0=1, nmo), isp=1, scfru)

        call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)

! Reverse the sign of ronei if DIRAC version is larger or equal to 21.
        if (dirac_version >= 21) then
            ronei(:, :, :) = -ronei(:, :, :)
        end if

        close (unit_mrconee)

        nmom = global_sec_end
        Allocate (one_elec_int_r(nmom, nmom)); Call memplus(KIND(one_elec_int_r), SIZE(one_elec_int_r), 1)
        Allocate (one_elec_int_i(nmom, nmom)); Call memplus(KIND(one_elec_int_i), SIZE(one_elec_int_i), 1)

! Store the one-electron integrals in energy order (CASPT2 order)
! one_elec_int_[r,i] are CASPT2 order
! roner, ronei are DIRAC order
        do i0 = 1, nmom
            do j0 = 1, nmom
                one_elec_int_r(indmo_dirac_to_cas(i0), indmo_dirac_to_cas(j0)) = roner(i0, j0, 1) ! using alpha component for a while
                one_elec_int_i(indmo_dirac_to_cas(i0), indmo_dirac_to_cas(j0)) = ronei(i0, j0, 1)
            end do
        end do

        Call memminus(KIND(roner), SIZE(roner), 1); deallocate (roner)
        Call memminus(KIND(ronei), SIZE(ronei), 1); deallocate (ronei)
    end subroutine read_1_elec_integrals

    subroutine sort_list_from_energy_order_to_ras_order(want_to_sort)
!===========================================================================================================================
! This subroutine sorts the want_to_sort list form orbital energy order
! to RAS order(ninact => ras1 => ras2 => ras3 => secondary).
!===========================================================================================================================
        use module_global_variables
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
! idx_ras_order must be global_act_start
        if (idx_ras_order /= global_act_start) then
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

! idx_ras_order must be global_act_end +1
        if (idx_ras_order /= global_sec_start) then
            print *, "ERROR: Sorting energy ascending order to RAS order is failed... STOP THE PROGRAM"
            print *, "ORIGINAL ENERGY ORDER LIST : ", mo_energy_order
            print *, "LIST OF SORTING IN PROGRESS: ", want_to_sort(1:idx_ras_order)
            call stop_with_errorcode(1) ! ERROR, STOP THE PROGRAM
        end if
! Fill secondary
        do while (idx_ras_order <= global_sec_end)
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

end subroutine read_mrconee
