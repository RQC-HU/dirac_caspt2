! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE casci

! CASCI calculation

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use module_dict, only: get_keys_vals, get_size
    use module_error, only: stop_with_errorcode
    use module_file_manager, only: open_unformatted_file
    use module_global_variables
    use module_realonly, only: realonly
    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer :: j0, j, i0, irec, unit_cimat
    real(8) :: cutoff_threshold

    complex*16, allocatable :: mat_complex(:, :) ! For complex
    real(8), allocatable    :: mat_real(:, :) ! For realonly
    real(8), allocatable    :: ecas(:)
    character(:), allocatable  :: filename
    character(len=len_convert_int_to_chr) :: chr_root
    integer :: datetmp0, datetmp1, dict_cas_idx_size, dict_cas_idx_reverse_size, idx
    integer, allocatable :: keys(:), vals(:), keys_rev(:), vals_rev(:)
    real(8) :: tsectmp0, tsectmp1

    Call search_cas_configuration

    ! Create a matrix for CI
    if (realonly%is_realonly()) then
        allocate (mat_real(ndet, ndet)); Call memplus(KIND(mat_real), SIZE(mat_real), 1)
        if (rank == 0) print *, "end allocate mat_real(ndet,ndet)"
        Call casmat_real(mat_real)
    else
        Allocate (mat_complex(ndet, ndet)); Call memplus(KIND(mat_complex), SIZE(mat_complex), 2)
        if (rank == 0) print *, "end allocate mat_complex(ndet,ndet)"
        Call casmat_complex(mat_complex)
    end if
    Allocate (ecas(ndet))
    ecas = 0.0d+00
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0

    ! Diagonalize the CI matrix
    if (rank == 0) then
        print *, 'Start mat diagonalization'
        print *, 'ndet before diagonalization', ndet
    end if
    cutoff_threshold = 0  ! No need to resolve linear dependence
    if (realonly%is_realonly()) then
        call rdiagx(mat_real, ndet, nroot, ecas)
    else
        Call cdiagx(mat_complex, ndet, nroot, ecas)
    end if
    if (rank == 0) print *, 'End mat diagonalization'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    ! keys and vals are used to store pairs of keys and values in dict_cas_idx
    dict_cas_idx_size = get_size(dict_cas_idx)
    allocate (keys(dict_cas_idx_size), vals(dict_cas_idx_size))
    call get_keys_vals(dict_cas_idx, keys, vals, dict_cas_idx_size)
    ! keys and vals are used to store pairs of keys and values in dict_cas_idx_reverse
    dict_cas_idx_reverse_size = get_size(dict_cas_idx_reverse)
    allocate (keys_rev(dict_cas_idx_reverse_size), vals_rev(dict_cas_idx_reverse_size))
    call get_keys_vals(dict_cas_idx_reverse, keys_rev, vals_rev, dict_cas_idx_reverse_size)
    ! Check if dict_cas_idx_size is equal to ndet
    if (dict_cas_idx_size /= ndet .or. dict_cas_idx_reverse_size /= ndet) then
        if (rank == 0) print *, 'ERROR: dict_cas_idx_size /= ndet .or. dict_cas_idx_reverse_size /= ndet. ndet =', ndet, &
            ",dict_cas_idx_size =", dict_cas_idx_size, ",dict_cas_idx_reverse_size =", dict_cas_idx_reverse_size
        call stop_with_errorcode(1)
    end if
! Print out CI matrix!
    if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        filename = 'CIMAT'
        call open_unformatted_file(unit=unit_cimat, file=filename, status='replace')
        write (unit_cimat) ndet
        write (unit_cimat) ecas(1:ndet)
        write (unit_cimat) dict_cas_idx_size ! The number of elements in dict_cas_idx
        do idx = 1, dict_cas_idx_size
            write (unit_cimat) keys(idx), vals(idx) ! Store pairs of keys and values in dict_cas_idx to the file
        end do
        write (unit_cimat) dict_cas_idx_reverse_size ! The number of elements in dict_cas_idx_reverse
        do idx = 1, dict_cas_idx_reverse_size
            write (unit_cimat) keys_rev(idx), vals_rev(idx) ! Store pairs of keys_rev and values in dict_cas_idx_reverse to the file
        end do
        close (unit_cimat)
    end if
    Allocate (cir(ndet, selectroot:selectroot)); Call memplus(KIND(cir), SIZE(cir), 1)
    Allocate (eigen(nroot)); Call memplus(KIND(eigen), SIZE(eigen), 1)

    cir(:, :) = 0.0d+00
    eigen(:) = 0.0d+00
    eigen(1:nroot) = ecas(1:nroot) + ecore
    if (rank == 0) then
        print '("CASCI ENERGY FOR ",I2," STATE")', totsym
        Do irec = 1, nroot
            write (chr_root, '(I4)') irec
            print '("CASCI Total Energy ROOT",a,F30.15," a.u.")', trim(adjustl(chr_root)), eigen(irec)
        End do
    end if
    if (realonly%is_realonly()) then
        cir(1:ndet, selectroot) = mat_real(1:ndet, selectroot)
        if (rank == 0) then
            do irec = 1, nroot
                print '("Root = ",I4)', irec
                do j = 1, ndet
                    if ((ABS(mat_real(j, irec))**2) > 1.0d-02) then
                        i0 = get_val(dict_cas_idx, j)
                        print *, (btest(i0, j0), j0=0, nact - 1)
                        print '(I4,2(3X,E14.7)," Weights ",E14.7)', &
                        & j, mat_real(j, irec), &
                        & ABS(mat_real(j, irec))**2
                    end if
                end do
            end do
        end if
    else
        Allocate (cii(ndet, selectroot:selectroot)); Call memplus(KIND(cii), SIZE(cii), 1)
        ! Print out the results
        cii(:, :) = 0.0d+00
        cir(1:ndet, selectroot) = DBLE(mat_complex(1:ndet, selectroot))
        cii(1:ndet, selectroot) = DIMAG(mat_complex(1:ndet, selectroot))
        if (rank == 0) then
            do irec = 1, nroot
                print '("Root = ",I4)', irec
                do j = 1, ndet
                    if ((ABS(mat_complex(j, irec))**2) > 1.0d-02) then
                        i0 = get_val(dict_cas_idx, j)
                        print *, (btest(i0, j0), j0=0, nact - 1)
                        print '(I4,2(3X,E14.7)," Weights ",E14.7)', &
                        & j, mat_complex(j, irec), &
                        & ABS(mat_complex(j, irec))**2
                    end if
                end do
            end do
        end if
        Call memminus(KIND(mat_complex), SIZE(mat_complex), 2); Deallocate (mat_complex)
    end if
    Deallocate (ecas)
    deallocate (keys, vals, keys_rev, vals_rev)
end subroutine casci
