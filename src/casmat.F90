! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE casmat_complex(mat)

! Creates CASCI matrix(mat)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use, intrinsic :: iso_fortran_env, only: int64
    use module_global_variables
    use module_dict, only: exists, get_val
    use module_index_utils, only: convert_global_to_active_idx, convert_active_to_global_idx
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    complex*16, intent(out) :: mat(ndet, ndet)

    integer              :: occ, vir, indr, inds, inda, indb
    integer              :: ir, is, ia, ib, imo
    integer              :: i0, j0, k0, l0
    integer(kind=int64)  :: i, j, newcas_idx1, newcas_idx2
    integer              :: phase1, phase2
    real(8)              :: i2r, i2i
    complex*16           :: cmplxint, mat0
    integer, allocatable :: oc(:), vi(:)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    mat = 0.0d+00

    if (rank == 0) print *, 'Cas mat enter'
    Allocate (oc(nelec))
    Allocate (vi(nact - nelec))
    ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
    Do i = rank + 1, ndet, nprocs

        occ = 0
        oc = 0
        vir = 0
        vi = 0

        Do imo = 1, nact
            If (BTEST(get_val(dict_cas_idx, i), imo - 1)) then
                occ = occ + 1
                oc(occ) = imo
            Else
                vir = vir + 1
                vi(vir) = imo
            End if
        End do

!! IDENTICAL DETERMINANT => DIAGONAL TERM
        !   diagonal term is same as Hartree-Fock's expression

        Do i0 = 1, ninact
            ir = i0
            cmplxint = DCMPLX(one_elec_int_r(ir, ir), one_elec_int_i(ir, ir))
            mat(i, i) = mat(i, i) + cmplxint
        End do

        Do i0 = 1, nelec
            ir = convert_active_to_global_idx(oc(i0))
            cmplxint = DCMPLX(one_elec_int_r(ir, ir), one_elec_int_i(ir, ir))
            mat(i, i) = mat(i, i) + cmplxint
        End do

        mat0 = 0.0d+00

        Do i0 = 1, ninact + nelec

            if (i0 <= ninact) then
                ir = i0
            Else
                indr = oc(convert_global_to_active_idx(i0))
                ir = convert_active_to_global_idx(indr)
            End if
            Do j0 = i0 + 1, ninact + nelec

                if (j0 <= ninact) then
                    is = j0
                Else
                    inds = oc(convert_global_to_active_idx(j0))
                    is = convert_active_to_global_idx(inds)
                End if

                ! two electron integral : (ir, ir | is, is)
                i2r = inttwr(ir, ir, is, is)
                i2i = inttwi(ir, ir, is, is)
                cmplxint = DCMPLX(i2r, i2i)

                mat0 = mat0 + 0.5d+00*cmplxint

                ! two electron integral : (ir, is | is, ir)
                i2r = inttwr(ir, is, is, ir)
                i2i = inttwi(ir, is, is, ir)
                cmplxint = DCMPLX(i2r, i2i)

                mat0 = mat0 - 0.5d+00*cmplxint
            End do
        End do

        mat(i, i) = mat(i, i) + mat0 + DCONJG(mat0)

!! ONE SPINOR DIFFERENCE

        Do i0 = 1, nelec
            indr = oc(i0)
            ir = convert_active_to_global_idx(indr)

            Do k0 = 1, nact - nelec
                inda = vi(k0)
                ia = convert_active_to_global_idx(inda)

                Call one_e_exct(get_val(dict_cas_idx, i), inda, indr, newcas_idx1, phase1)

                if (exists(dict_cas_idx_reverse, newcas_idx1)) then
                    j = get_val(dict_cas_idx_reverse, newcas_idx1)
                Else
                    cycle ! Next k0 (Because newcas_idx1 is not in dict_cas_idx_reverse)
                End if

                If (j > i) then
                    cmplxint = DCMPLX(one_elec_int_r(ir, ia), one_elec_int_i(ir, ia))
                    mat(i, j) = mat(i, j) + cmplxint
                    Do l0 = 1, ninact
                        is = l0

                        ! two electron integral : (ir, ia | is, is)
                        i2r = inttwr(ir, ia, is, is)
                        i2i = inttwi(ir, ia, is, is)
                        cmplxint = DCMPLX(i2r, i2i)

                        mat(i, j) = mat(i, j) + cmplxint

                        ! two electron integral : (ir, is | is, ia)
                        i2r = inttwr(ir, is, is, ia)
                        i2i = inttwi(ir, is, is, ia)
                        cmplxint = DCMPLX(i2r, i2i)

                        mat(i, j) = mat(i, j) - cmplxint
                    End do      !l0

                    Do l0 = 1, nelec
                        inds = oc(l0)
                        is = convert_active_to_global_idx(inds)

                        ! two electron integral : (ir, ia | is, is)
                        i2r = inttwr(ir, ia, is, is)
                        i2i = inttwi(ir, ia, is, is)
                        cmplxint = DCMPLX(i2r, i2i)

                        mat(i, j) = mat(i, j) + cmplxint

                        ! two electron integral : (ir, is | is, ia)
                        i2r = inttwr(ir, is, is, ia)
                        i2i = inttwi(ir, is, is, ia)
                        cmplxint = DCMPLX(i2r, i2i)

                        mat(i, j) = mat(i, j) - cmplxint
                    End do

                    mat(i, j) = merge(1, -1, mod(phase1, 2)==0)*mat(i, j)
                    mat(j, i) = DCONJG(mat(i, j))

                End if
            End do
        End do
!! TWO ELECTRON DIFFERNT CASE
        Do i0 = 1, nelec
            Do j0 = i0 + 1, nelec
                indr = oc(i0)
                inds = oc(j0)
                ir = convert_active_to_global_idx(indr)
                is = convert_active_to_global_idx(inds)

                Do k0 = 1, nact - nelec
                    Do l0 = k0 + 1, nact - nelec
                        inda = vi(k0)
                        indb = vi(l0)
                        ia = convert_active_to_global_idx(inda)
                        ib = convert_active_to_global_idx(indb)

                        Call one_e_exct(get_val(dict_cas_idx, i), inda, indr, newcas_idx1, phase1)
                        Call one_e_exct(newcas_idx1, indb, inds, newcas_idx2, phase2)

                        if (exists(dict_cas_idx_reverse, newcas_idx2)) then
                            j = get_val(dict_cas_idx_reverse, newcas_idx2)
                        Else
                            cycle ! Next k0 (Because newcas_idx2 is not in dict_cas_idx_reverse)
                        End if

                        If (j > i) then
                            ! two electron integral : (ir, ia | is, ib)
                            i2r = inttwr(ir, ia, is, ib)
                            i2i = inttwi(ir, ia, is, ib)
                            cmplxint = DCMPLX(i2r, i2i)
                            mat(i, j) = cmplxint

                            ! two electron integral : (ir, ib | is, ia)
                            i2r = inttwr(ir, ib, is, ia)
                            i2i = inttwi(ir, ib, is, ia)
                            cmplxint = DCMPLX(i2r, i2i)
                            mat(i, j) = mat(i, j) - cmplxint

                            mat(i, j) = merge(1, -1, mod(phase1 + phase2, 2)==0)*mat(i, j)
                            mat(j, i) = DCONJG(mat(i, j))
                        End if
                    End do
                End do
            End do
        End do
    End do

    Deallocate (oc)
    Deallocate (vi)
#ifdef HAVE_MPI
    call allreduce_wrapper(mat=mat)
#endif
    if (rank == 0) print *, 'end casmat_complex'
end subroutine casmat_complex
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE casmat_real(mat)

    ! Creates CASCI matrix(mat)

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use, intrinsic :: iso_fortran_env, only: int64
    use module_global_variables
    use module_dict, only: exists, get_val
    use module_index_utils, only: convert_global_to_active_idx, convert_active_to_global_idx
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    real(8), intent(out) :: mat(ndet, ndet)

    integer              :: occ, vir, indr, inds, inda, indb
    integer              :: ir, is, ia, ib, imo
    integer              :: i0, j0, k0, l0
    integer(kind=int64)  :: i, j, newcas_idx1, newcas_idx2
    integer              :: phase1, phase2
    real(8)              :: mat0, i2r
    integer, allocatable :: oc(:), vi(:)

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    mat = 0.0d+00

    if (rank == 0) print *, 'Cas mat enter'
    Allocate (oc(nelec))
    Allocate (vi(nact - nelec))
    ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
    Do i = rank + 1, ndet, nprocs

        occ = 0
        oc = 0
        vir = 0
        vi = 0

        Do imo = 1, nact
            If (BTEST(get_val(dict_cas_idx, i), imo - 1)) then
                occ = occ + 1
                oc(occ) = imo
            Else
                vir = vir + 1
                vi(vir) = imo
            End if
        End do

    !! IDENTICAL DETERMINANT => DIAGONAL TERM
        !   diagonal term is same as Hartree-Fock's expression

        Do i0 = 1, ninact
            ir = i0
            mat(i, i) = mat(i, i) + one_elec_int_r(ir, ir)
        End do

        Do i0 = 1, nelec
            ir = convert_active_to_global_idx(oc(i0))
            mat(i, i) = mat(i, i) + one_elec_int_r(ir, ir)
        End do

        mat0 = 0.0d+00

        Do i0 = 1, ninact + nelec

            if (i0 <= ninact) then
                ir = i0
            Else
                indr = oc(convert_global_to_active_idx(i0))
                ir = convert_active_to_global_idx(indr)
            End if
            Do j0 = i0 + 1, ninact + nelec

                if (j0 <= ninact) then
                    is = j0
                Else
                    inds = oc(convert_global_to_active_idx(j0))
                    is = convert_active_to_global_idx(inds)
                End if

                ! two electron integral : (ir, ir | is, is)
                i2r = inttwr(ir, ir, is, is)
                mat0 = mat0 + 0.5d+00*i2r

                ! two electron integral : (ir, is | is, ir)
                i2r = inttwr(ir, is, is, ir)
                mat0 = mat0 - 0.5d+00*i2r
            End do
        End do

        mat(i, i) = mat(i, i) + 2*mat0

    !! ONE SPINOR DIFFERENCE

        Do i0 = 1, nelec
            indr = oc(i0)
            ir = convert_active_to_global_idx(indr)

            Do k0 = 1, nact - nelec
                inda = vi(k0)
                ia = convert_active_to_global_idx(inda)

                Call one_e_exct(get_val(dict_cas_idx, i), inda, indr, newcas_idx1, phase1)

                if (exists(dict_cas_idx_reverse, newcas_idx1)) then
                    j = get_val(dict_cas_idx_reverse, newcas_idx1)
                Else
                    cycle ! Next k0 (Because newcas_idx1 is not in dict_cas_idx_reverse)
                End if

                If (j > i) then
                    mat(i, j) = mat(i, j) + one_elec_int_r(ir, ia)
                    Do l0 = 1, ninact
                        is = l0

                        ! two electron integral : (ir, ia | is, is)
                        i2r = inttwr(ir, ia, is, is)
                        mat(i, j) = mat(i, j) + i2r

                        ! two electron integral : (ir, is | is, ia)
                        i2r = inttwr(ir, is, is, ia)
                        mat(i, j) = mat(i, j) - i2r
                    End do      !l0

                    Do l0 = 1, nelec
                        inds = oc(l0)
                        is = convert_active_to_global_idx(inds)

                        ! two electron integral : (ir, ia | is, is)
                        i2r = inttwr(ir, ia, is, is)
                        mat(i, j) = mat(i, j) + i2r

                        ! two electron integral : (ir, is | is, ia)
                        i2r = inttwr(ir, is, is, ia)
                        mat(i, j) = mat(i, j) - i2r
                    End do

                    mat(i, j) = merge(1, -1, mod(phase1, 2)==0)*mat(i, j)
                    mat(j, i) = mat(i, j)

                End if
            End do
        End do
    !! TWO ELECTRON DIFFERNT CASE
        Do i0 = 1, nelec
            Do j0 = i0 + 1, nelec
                indr = oc(i0)
                inds = oc(j0)
                ir = convert_active_to_global_idx(indr)
                is = convert_active_to_global_idx(inds)

                Do k0 = 1, nact - nelec
                    Do l0 = k0 + 1, nact - nelec
                        inda = vi(k0)
                        indb = vi(l0)
                        ia = convert_active_to_global_idx(inda)
                        ib = convert_active_to_global_idx(indb)

                        Call one_e_exct(get_val(dict_cas_idx, i), inda, indr, newcas_idx1, phase1)
                        Call one_e_exct(newcas_idx1, indb, inds, newcas_idx2, phase2)

                        if (exists(dict_cas_idx_reverse, newcas_idx2)) then
                            j = get_val(dict_cas_idx_reverse, newcas_idx2)
                        Else
                            cycle ! Next k0 (Because newcas_idx2 is not in dict_cas_idx_reverse)
                        End if

                        If (j > i) then

                            ! two electron integral : (ir, ia | is, ib)
                            i2r = inttwr(ir, ia, is, ib)
                            mat(i, j) = i2r

                            ! two electron integral : (ir, ib | is, ia)
                            i2r = inttwr(ir, ib, is, ia)
                            mat(i, j) = mat(i, j) - i2r

                            mat(i, j) = merge(1, -1, mod(phase1 + phase2, 2)==0)*mat(i, j)
                            mat(j, i) = mat(i, j)
                        End if
                    End do
                End do
            End do
        End do
    End do

    Deallocate (oc)
    Deallocate (vi)
#ifdef HAVE_MPI
    call allreduce_wrapper(mat=mat)
#endif
    if (rank == 0) print *, 'end casmat_real'
end subroutine casmat_real
