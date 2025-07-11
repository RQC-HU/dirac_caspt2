! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

subroutine r4dcaspt2_tra   ! DO CASPT2 CALC WITH MO TRANSFORMATION

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use dcaspt2_restart_file, only: get_subspace_idx, is_skip_restart_file_subspace
    use module_2integrals, only: readint2_casci
    use module_dict, only: add
    use module_error, only: stop_with_errorcode
    use module_file_manager
    use module_global_variables
    use module_intra, only: intra_1, intra_2, intra_3
    use module_realonly, only: realonly
    use module_time
    use read_input_module, only: read_input
    Implicit NONE
    integer                 :: unit_input, cur_subspace_idx
    real(8)                 :: e0, e2, weight0
    character(:), allocatable       :: filename
    character(*), parameter         :: int_input_form = '(1x,a,1x,i0)'
    character(len=30)               :: real_str
    integer                 :: idx, i, nuniq

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    if (rank == 0) then
        call print_head_caspt2
        print *, ''
        print *, ' START RELATIVISIC CASPT2 PROGRAM'
        print *, ''
    end if
    tmem = 0.0d+00

    call write_allocated_memory_size
    call get_current_time(init_time); call print_time(init_time); start_time = init_time

    if (rank == 0) then
        print int_input_form, 'ninact        =', ninact
        print int_input_form, 'nact          =', nact
        print int_input_form, 'nsec          =', nsec
        print int_input_form, 'nelec         =', nelec
        print int_input_form, 'nroot         =', nroot
        print int_input_form, 'selectroot    =', selectroot
        print int_input_form, 'totsym        =', totsym
        write (real_str, '(E20.10)') eshift
        print '(1x,a,1x,a)', 'eshift        =', trim(adjustl(real_str))
        print int_input_form, 'diracver      =', dirac_version
        print int_input_form, 'scheme        =', mdcint_scheme
        if (ras1_size /= 0) print *, "RAS1          =", ras1_list
        if (ras2_size /= 0) print *, "RAS2          =", ras2_list
        if (ras3_size /= 0) print *, "RAS3          =", ras3_list
        print int_input_form, 'ras1_max_hole =', ras1_max_hole
        print int_input_form, 'ras3_max_elec =', ras3_max_elec
        print int_input_form, 'minholeras1   =', min_hole_ras1
        print int_input_form, 'minelecras3   =', min_elec_ras3
        print *, 'debugprint    =', debug
        if (enable_restart) print *, "restart       =", enable_restart
        print *, ''
    end if

    if (ninact == 0 .and. nsec == 0) then
        if (rank == 0) print *, "The CASPT2 energy cannot be defined when ninact = 0 and nsec = 0."
        stop
    end if

    if (rank == 0) print *, 'Reading MDCINT (2-e integrals)'

    ! Get MDCINTNEWx's filename and subspace filename
    call get_mdcint_filename(0)
    call get_subspace_filename

    ! READ MDCINTNEWx's file and devide into each subspace files.
    Call divide_2_elec_integral_into_subspaces(mdcintnew)

    call write_allocated_memory_size
    if (rank == 0) print *, 'nmo        =', nmo

    call read_cidata

    iroot = selectroot
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!            BUILDING  FOCK MATRIX               !
!  fij = hij + SIGUMA[<0|Ekl|0>{(ij|kl)-(il|kj)} !
!                 kl                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    ! Read UTChem type MDCINT files and expands the 2-electron integral in memory
    Call readint2_casci(mdcintnew, nuniq)
#ifdef DEBUG
    !! TEST TO CALCULATE FOCK MATRIX OF HF STATE fpq = hpq + SIGUMA_r[(pq|rr)-(pr|qr)]
    !! THIS MUST BE DIAGONAL MATRIX AND DIAGONAL ELEMENTS CORESPONDS TO SPINOR ENERGIES.
    if (realonly%is_realonly()) then
        if (.not. allocated(fock_real)) then
            allocate (fock_real(nmo, nmo)); Call memplus(KIND(fock_real), SIZE(fock_real), 1)
        end if
        fock_real(:, :) = 0.0d+00
        call fock_matrix_of_hf_real
    else
        if (.not. allocated(fock_cmplx)) then
            Allocate (fock_cmplx(nmo, nmo)); Call memplus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)
        end if
        fock_cmplx(:, :) = 0.0d+00
        call fock_matrix_of_hf_complex
    End if
#endif

    !! NOW MAKE FOCK MATRIX FOR CASCI STATE
    !! fij = hij + SIGUMA_kl[<0|Ekl|0>{(ij|kl)-(il|kj)}
    call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
    if (realonly%is_realonly()) then
        if (.not. allocated(fock_real)) then
            allocate (fock_real(nmo, nmo)); call memplus(KIND(fock_real), SIZE(fock_real), 1)
        end if
        fock_real(:, :) = 0.0d+00
        Call fockcasci_real
    else
        if (.not. allocated(fock_cmplx)) then
            allocate (fock_cmplx(nmo, nmo)); call memplus(KIND(fock_cmplx), SIZE(fock_cmplx), 2)
        end if
        fock_cmplx(:, :) = 0.0d+00
        Call fockcasci_complex
    end if

    if (rank == 0) then
        print *, 'end building fock'
    end if
    call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

#ifdef DEBUG
    if (rank == 0) call prtoutfock
#endif

    Allocate (eps(nmo)); Call memplus(KIND(eps), SIZE(eps), 1)
    eps = 0.0d+00

    ! Diagonalize the Fock matrix
    Call fockdiag

    ! Print orbital energies
#ifdef DEBUG
    if (rank == 0) then
        print *, debug, "debug"
        Do i = 1, nmo
            print *, 'eps(', i, ')=', eps(i)
        End do
    end if
#endif

    ! Print the irreducible representation that calculates energy
    if (rank == 0) then
        print *, ' '
        print *, '*******************************'
        print *, ' '
        print '("IRREP = ",6A)', repna(totsym)
        print *, ' '
        print *, '*******************************'
        print *, ' '
    end if

    ! Calculate eigenvalues of a 0th-order Hamiltonian applied to a 0th-order wave function
    Call calce0(e0)

    ! Initialize the date, time and the 2nd order energy
    e2 = 0.0d+00
    e2all = sum(e2_subspace)
    sumc2 = sum(sumc2_subspace)
    call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

    ! Subspace A
    cur_subspace_idx = get_subspace_idx("A")
    if (ninact == 0) then
        if (rank == 0) print *, "Skip the calculation of A subspace 2nd order energy &
&        because the 2nd order energy of A subspace cannot be defined when ninact = 0."
    else if (is_skip_restart_file_subspace("A")) then
        if (rank == 0) then
            print *, "Skip the calculation of A subspace 2nd order energy because of the caspt2_restart file."
            print '("e2a      = ",E25.15," a.u.")', e2_subspace(cur_subspace_idx)
            print '("sumc2,a  = ",E25.15)', sumc2_subspace(cur_subspace_idx)
        end if
    else
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, "Start calcultion of A subspace 2nd order energy"
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'

        ! Transform A subspace 2-electron integrals (active, inactive | active, active)
        Call intra_3(2, 1, 2, 2, a1int)

        if (debug .and. rank == 0) print *, 'End intra3 A1int'
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Transform A subspace 2-electron integrals (active, inactive | inactive, inactive)
        Call intra_3(2, 1, 1, 1, a2int)
        if (debug .and. rank == 0) print *, 'End intra_3 A2int'
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the A subspace 2nd order energy
        Call solve_A_subspace(e0)
        if (rank == 0) print *, " "
        if (rank == 0) print *, "End calcultion of A subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, " "
    end if

    ! Subspace B
    cur_subspace_idx = get_subspace_idx("B")
    if (ninact == 0) then
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, "Skip the calculation of B subspace 2nd order energy &
&        because the 2nd order energy of B subspace cannot be defined when ninact = 0."
    else if (is_skip_restart_file_subspace("B")) then
        if (rank == 0) then
            print *, "Skip the calculation of B subspace 2nd order energy because of the caspt2_restart file."
            print '("e2b      = ",E25.15," a.u.")', e2_subspace(cur_subspace_idx)
            print '("sumc2,b  = ",E25.15)', sumc2_subspace(cur_subspace_idx)
        end if
    else
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, "Start calcultion of B subspace 2nd order energy"
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        ! Transform B subspace 2-electron integrals (active, inactive | active, inactive)
        Call intra_2(2, 1, 2, 1, bint)
        if (debug .and. rank == 0) print *, "End intra_2 Bint"
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the B subspace 2nd order energy
        Call solve_B_subspace(e0)
        if (rank == 0) print *, " "
        if (rank == 0) print *, "End calcuation of B subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, " "
    end if

    ! Subspace C
    cur_subspace_idx = get_subspace_idx("C")
    if (nsec == 0) then
        if (rank == 0) print *, "Skip the calculation of C subspace 2nd order energy &
&        because the 2nd order energy of C subspace cannot be defined when nsec = 0."
    else if (is_skip_restart_file_subspace("C")) then
        if (rank == 0) then
            print *, "Skip the calculation of C subspace 2nd order energy because of the caspt2_restart file."
            print '("e2c      = ",E25.15," a.u.")', e2_subspace(cur_subspace_idx)
            print '("sumc2,c  = ",E25.15)', sumc2_subspace(cur_subspace_idx)
        end if
    else
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, "Start calcultion of C subspace 2nd order energy"
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        ! Transform C subspace 2-electron integrals (secondary, active | active, active)
        Call intra_3(3, 2, 2, 2, c1int)
        if (debug .and. rank == 0) print *, 'End intra_3 C1int'
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Transform C subspace 2-electron integrals (secondary, active | inactive, inactive)
        Call intra_3(3, 2, 1, 1, c2int)
        if (debug .and. rank == 0) print *, 'End intra_3 C2int'
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Transform C subspace 2-electron integrals (secondary, inactive | inactive, active)
        Call intra_1(3, 1, 1, 2, c3int)
        if (debug .and. rank == 0) print *, 'End intra_1 C3int'
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the B subspace 2nd order energy
        Call solve_C_subspace(e0)
        if (rank == 0) print *, " "
        if (rank == 0) print *, "End calcuation of C subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, " "
    end if

    ! Subspace D
    cur_subspace_idx = get_subspace_idx("D")
    if (ninact == 0 .or. nsec == 0) then
        if (rank == 0) print *, "Skip the calculation of D subspace 2nd order energy &
&        because the 2nd order energy of D subspace cannot be defined when ninact = 0 or nsec = 0."
    else if (is_skip_restart_file_subspace("D")) then
        if (rank == 0) then
            print *, "Skip the calculation of D subspace 2nd order energy because of the caspt2_restart file."
            print '("e2d      = ",E25.15," a.u.")', e2_subspace(cur_subspace_idx)
            print '("sumc2,d  = ",E25.15)', sumc2_subspace(cur_subspace_idx)
        end if
    else
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, "Start calcultion of D subspace 2nd order energy"
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        ! Transform D subspace 2-electron integrals (secondary, inactive | active, active)
        Call intra_3(3, 1, 2, 2, d1int)
        if (debug .and. rank == 0) print *, 'End intra_1 D1int'
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Transform D subspace 2-electron integrals (secondary, active | active, inactive)
        Call intra_1(3, 2, 2, 1, d2int)
        if (debug .and. rank == 0) print *, 'End intra_1 D2int'
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Transform D subspace 2-electron integrals (secondary, inactive | inactive, inactive)
        Call intra_3(3, 1, 1, 1, d3int)
        if (debug .and. rank == 0) print *, 'End intra_1 D3int'
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the D subspace 2nd order energy
        Call solve_D_subspace(e0)
        if (rank == 0) print *, " "
        if (rank == 0) print *, "End calcuation of D subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, " "
    end if

    ! Subspace E
    cur_subspace_idx = get_subspace_idx("E")
    if (ninact == 0 .or. nsec == 0) then
        if (rank == 0) print *, "Skip the calculation of E subspace 2nd order energy &
&        because the 2nd order energy of E subspace cannot be defined when ninact = 0 or nsec = 0."
    else if (is_skip_restart_file_subspace("E")) then
        if (rank == 0) then
            print *, "Skip the calculation of E subspace 2nd order energy because of the caspt2_restart file."
            print '("e2e      = ",E25.15," a.u.")', e2_subspace(cur_subspace_idx)
            print '("sumc2,e  = ",E25.15)', sumc2_subspace(cur_subspace_idx)
        end if
    else
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, "Start calcultion of E subspace 2nd order energy"
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        ! Transform E subspace 2-electron integrals (secondary, active | active, inactive)
        Call intra_1(3, 1, 2, 1, eint)
        if (debug .and. rank == 0) print *, 'End intra_1 Eint'
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the E subspace 2nd order energy
        Call solve_E_subspace(e0)
        if (rank == 0) print *, " "
        if (rank == 0) print *, "End calcuation of E subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, " "
    end if

    ! Subspace F
    cur_subspace_idx = get_subspace_idx("F")
    if (nsec == 0) then
        if (rank == 0) print *, "Skip the calculation of F subspace 2nd order energy &
&        because the 2nd order energy of F subspace cannot be defined when nsec = 0."
    else if (is_skip_restart_file_subspace("F")) then
        if (rank == 0) then
            print *, "Skip the calculation of F subspace 2nd order energy because of the caspt2_restart file."
            print '("e2f      = ",E25.15," a.u.")', e2_subspace(cur_subspace_idx)
            print '("sumc2,f  = ",E25.15)', sumc2_subspace(cur_subspace_idx)
        end if
    else
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, "Start calcultion of F subspace 2nd order energy"
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        ! Transform F subspace 2-electron integrals (secondary, active | secondary, active)
        Call intra_2(3, 2, 3, 2, fint)
        if (debug .and. rank == 0) print *, 'End intra_1 Fint'
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the F subspace 2nd order energy
        Call solve_F_subspace(e0)
        if (rank == 0) print *, " "
        if (rank == 0) print *, "End calcuation of F subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, " "
    end if

    ! Subspace G
    cur_subspace_idx = get_subspace_idx("G")
    if (ninact == 0 .or. nsec == 0) then
        if (rank == 0) print *, "Skip the calculation of G subspace 2nd order energy &
&        because the 2nd order energy of G subspace cannot be defined when ninact = 0 or nsec = 0."
    else if (is_skip_restart_file_subspace("G")) then
        if (rank == 0) then
            print *, "Skip the calculation of G subspace 2nd order energy because of the caspt2_restart file."
            print '("e2g      = ",E25.15," a.u.")', e2_subspace(cur_subspace_idx)
            print '("sumc2,g  = ",E25.15)', sumc2_subspace(cur_subspace_idx)
        end if
    else
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, "Start calcultion of G subspace 2nd order energy"
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        ! Transform G subspace 2-electron integrals (secondary, inactive | secondary, active)
        Call intra_1(3, 1, 3, 2, gint)
        if (debug .and. rank == 0) print *, 'End intra_1 Gint'
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the G subspace 2nd order energy
        Call solve_G_subspace(e0)
        if (rank == 0) print *, " "
        if (rank == 0) print *, "End calcuation of G subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, " "
    end if

    ! Subspace H
    cur_subspace_idx = get_subspace_idx("H")
    if (ninact == 0 .or. nsec == 0) then
        if (rank == 0) print *, "Skip the calculation of H subspace 2nd order energy &
&        because the 2nd order energy of H subspace cannot be defined when ninact = 0 or nsec = 0."
    else if (is_skip_restart_file_subspace("H")) then
        if (rank == 0) then
            print *, "Skip the calculation of H subspace 2nd order energy because of the caspt2_restart file."
            print '("e2h      = ",E25.15," a.u.")', e2_subspace(cur_subspace_idx)
            print '("sumc2,h  = ",E25.15)', sumc2_subspace(cur_subspace_idx)
        end if
    else
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, "Start calcultion of H subspace 2nd order energy"
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        ! Transform H subspace 2-electron integrals (secondary, inactive | secondary, inactive)
        if (debug .and. rank == 0) print *, 'Enter intra_2 Hint'
        Call intra_2(3, 1, 3, 1, hint)
        if (debug .and. rank == 0) print *, 'End intra_2 Hint'
        if (debug) call get_current_time_and_print_diff(start_time, end_time); start_time = end_time

        ! Calculate the H subspace 2nd order energy
        Call solve_H_subspace(e0)
        if (rank == 0) print *, " "
        if (rank == 0) print *, "End calcuation of H subspace 2nd order energy"
        call get_current_time_and_print_diff(start_time, end_time); start_time = end_time
        if (rank == 0) print '(64A)', '----------------------------------------------------------------'
        if (rank == 0) print *, " "
    end if

    ! Calculate and print the weight of the 0th wave function
    weight0 = 1.0d+00/(1.0d+00 + sumc2)
    ! Print out the total 2nd order energy

    if (rank == 0) then
        print '(" CASPT2 ENERGY FOR ",I0," STATE, selectroot = ",I0)', totsym, selectroot
        print '(" c^2 is                         ",F30.20)', sumc2
        print '(" weight of 0th wave function is ",F30.20)', weight0
        print '(" Total second order energy is   ",F30.20," a.u.")', e2all - eshift*sumc2
        print '(" Total energy is                ",F30.20," a.u.")', e2all + eigen(iroot) - eshift*sumc2
        print *, ' '
        print *, 'END OF RELATIVISTIC CASPT2 PROGRAM'
    end if

    if (allocated(cir)) then
        Call memminus(KIND(cir), SIZE(cir), 1); deallocate (cir)
    end if
    if (allocated(cii)) then
        Call memminus(KIND(cii), SIZE(cii), 1); deallocate (cii)
    end if
    if (allocated(eigen)) then
        Call memminus(KIND(eigen), SIZE(eigen), 1); deallocate (eigen)
    end if
    if (allocated(eps)) then
        Call memminus(KIND(eps), SIZE(eps), 1); deallocate (eps)
    end if
    if (allocated(inttwi)) then
        Call memminus(KIND(inttwi), SIZE(inttwi), 1); deallocate (inttwi)
    end if
    if (allocated(inttwr)) then
        Call memminus(KIND(inttwr), SIZE(inttwr), 1); deallocate (inttwr)
    end if
    if (allocated(int2r_f1)) then
        Call memminus(KIND(int2r_f1), SIZE(int2r_f1), 1); deallocate (int2r_f1)
    end if
    if (allocated(int2r_f2)) then
        Call memminus(KIND(int2r_f2), SIZE(int2r_f2), 1); deallocate (int2r_f2)
    end if
    if (.not. realonly%is_realonly()) then
        if (allocated(int2i_f1)) then
            Call memminus(KIND(int2i_f1), SIZE(int2i_f1), 1); deallocate (int2i_f1)
        end if
        if (allocated(int2i_f2)) then
            Call memminus(KIND(int2i_f2), SIZE(int2i_f2), 1); deallocate (int2i_f2)
        end if
    end if
    ! Print out the total time
!    call get_current_time_and_print_diff(init_time, end_time)

END subroutine r4dcaspt2_tra
