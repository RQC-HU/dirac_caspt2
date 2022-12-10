! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE casmat(mat)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    complex*16, intent(out) :: mat(ndet, ndet)

    integer              :: occ, vir, indr, inds, inda, indb
    integer              :: ir, is, ia, ib, imo
    integer              :: i0, j0, k0, l0, i, j, newidet1, newidet2
    integer              :: phase, phase1, phase2
    real*8               :: i2r, i2i
    complex*16           :: cmplxint, mat0
    integer, allocatable :: oc(:), vi(:)

    mat = 0.0d+00

    if (rank == 0) print *, 'Cas mat enter'
    Allocate (oc(nelec))
    Allocate (vi(nact - nelec))
    if (rank == 0) print *, 'allocated oc and vi'
    Do i = rank + 1, ndet, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)

        occ = 0
        oc = 0
        vir = 0
        vi = 0

        Do imo = 1, nact
            If (BTEST(idet(i), imo - 1)) then
                occ = occ + 1
                oc(occ) = imo
            Else
                vir = vir + 1
                vi(vir) = imo
            End if
        End do

!! IDENTICAL DETERMINANT => DIAGONAL TERM
        !   diagonal term is same as Hartree-Fock's expression

        !    cmplxint = 0.0d+00
        Do i0 = 1, ninact
            ir = i0
            cmplxint = DCMPLX(oner(ir, ir), onei(ir, ir))
            mat(i, i) = mat(i, i) + cmplxint
        End do

        Do i0 = 1, nelec
            indr = oc(i0)
            ir = indr + ninact
            cmplxint = DCMPLX(oner(ir, ir), onei(ir, ir))
            mat(i, i) = mat(i, i) + cmplxint
        End do

        mat0 = 0.0d+00

        Do i0 = 1, ninact + nelec

            if (i0 <= ninact) then
                ir = i0
            Else
                indr = i0 - ninact
                indr = oc(indr)
                ir = indr + ninact
            End if
            Do j0 = i0 + 1, ninact + nelec

                if (j0 <= ninact) then
                    is = j0
                Else
                    inds = j0 - ninact
                    inds = oc(inds)
                    is = inds + ninact
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
            ir = indr + ninact

            Do k0 = 1, nact - nelec
                inda = vi(k0)
                ia = inda + ninact

                Call one_e_exct(idet(i), inda, indr, newidet1, phase1)

                j = idetr(newidet1)

                If (j > i) then
                    cmplxint = DCMPLX(oner(ir, ia), onei(ir, ia))
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
                        is = inds + ninact

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

                    if (mod(phase1, 2) == 0) phase = 1.0d+00
                    if (mod(phase1, 2) == 1) phase = -1.0d+00

                    mat(i, j) = phase*mat(i, j)
                    mat(j, i) = DCONJG(mat(i, j))

                End if

            End do            ! k0
        End do               ! i0
!! TWO ELECTRON DIFFERNT CASE

        Do i0 = 1, nelec
            Do j0 = i0 + 1, nelec
                indr = oc(i0)
                inds = oc(j0)
                ir = indr + ninact
                is = inds + ninact

                Do k0 = 1, nact - nelec
                    Do l0 = k0 + 1, nact - nelec
                        inda = vi(k0)
                        indb = vi(l0)
                        ia = inda + ninact
                        ib = indb + ninact

                        Call one_e_exct(idet(i), inda, indr, newidet1, phase1)
                        Call one_e_exct(newidet1, indb, inds, newidet2, phase2)

                        j = idetr(newidet2)

                        If (j > i) then
                            if (mod(phase1 + phase2, 2) == 0) phase = 1.0d+00
                            if (mod(phase1 + phase2, 2) == 1) phase = -1.0d+00

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

                            mat(i, j) = phase*mat(i, j)
                            mat(j, i) = DCONJG(mat(i, j))

                        End if

                    End do      ! l0
                End do         ! k0
            End do            ! j0
        End do               ! i0

    End do                  ! i

    Deallocate (oc)
    Deallocate (vi)
    if (rank == 0) then
        print *, 'end casmat'
        print *, 'Reduce mat(:,:)'
    end if
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, mat(1, 1), ndet**2, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (rank == 0) print *, 'end allreduce mat(:,:)'
#endif
end subroutine casmat
