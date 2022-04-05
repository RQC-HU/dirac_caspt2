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
    integer              :: ir, is, ia, ib, imo, nint
    integer              :: i0, j0, k0, l0, i, j, newidet1, newidet2
    integer              :: phase, phase1, phase2
    real*8               :: nsign, i2r, i2i
    complex*16           :: cmplxint, mat0
    integer              :: nprocs_for_mat
    integer, allocatable :: ridet(:), oc(:), vi(:)

    mat = 0.0d+00

    if (rank == 0) then ! Process limits for output
        write (*, *) 'Cas mat enter'
    end if
    Allocate (oc(nelec))
    Allocate (vi(nact - nelec))
    if (rank == 0) then ! Process limits for output
        write (*, *) 'allocated oc and vi', rank
    end if
#ifdef BIG_MAT
    ! If only the master process has a matrix named mat,
    ! this subroutine will also excutes only by the master process.(So nprocs must be 1)
    nprocs_for_mat = 1
#else
    nprocs_for_mat = nprocs
#endif
    Do i = rank + 1, ndet, nprocs_for_mat ! MPI parallelization (Distributed loop: static scheduling, per nprocs)

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

!           write(*,*) 'i, idet(i)',i, idet(i)
!           write(*,*) occ, oc(1:occ)
!           write(*,*) vir, vi(1:vir)
!           write(*,*) ' '

!! IDENTICAL DETERMINANT => DIAGONAL TERM
        !   diagonal term is same as Hartree-Fock's expression

        !    cmplxint = 0.0d+00
        Do i0 = 1, ninact
            ir = i0
            cmplxint = CMPLX(oner(ir, ir), onei(ir, ir), 16)
            mat(i, i) = mat(i, i) + cmplxint
        End do

        Do i0 = 1, nelec
            indr = oc(i0)
            ir = indr + ninact
            cmplxint = CMPLX(oner(ir, ir), onei(ir, ir), 16)
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
                cmplxint = CMPLX(i2r, i2i, 16)
!                 write(*,*)ir,is,cmplxint

                mat0 = mat0 + 0.5d+00*cmplxint
!                 mat(i,i) = mat(i,i) + 0.5d+00*cmplxint

                ! two electron integral : (ir, is | is, ir)
                i2r = inttwr(ir, is, is, ir)
                i2i = inttwi(ir, is, is, ir)
                cmplxint = CMPLX(i2r, i2i, 16)
!                 write(*,*)ir,is,cmplxint

                mat0 = mat0 - 0.5d+00*cmplxint
!                 mat(i,i) = mat(i,i) - 0.5d+00*cmplxint

            End do
        End do

        mat(i, i) = mat(i, i) + mat0 + DCONJG(mat0)

!           write(*,*)'mat(',i,',',i,') = ',mat(i,i)

!! ONE SPINOR DIFFERENCE

        Do i0 = 1, nelec
            indr = oc(i0)
            ir = indr + ninact

            Do k0 = 1, nact - nelec
                inda = vi(k0)
                ia = inda + ninact

                Call one_e_exct(idet(i), inda, indr, newidet1, phase1)

                j = idetr(newidet1)

!                 write(*,*)'j=',j

                If (j > i) then
                    !    if (rank == 0) write (*, '(A,I5,A,I5)') 'Noda ijorder i:', i, ' j:', j
                    cmplxint = CMPLX(oner(ir, ia), onei(ir, ia), 16)
                    mat(i, j) = mat(i, j) + cmplxint
                    Do l0 = 1, ninact
                        is = l0

                        ! two electron integral : (ir, ia | is, is)
                        i2r = inttwr(ir, ia, is, is)
                        i2i = inttwi(ir, ia, is, is)
                        cmplxint = CMPLX(i2r, i2i, 16)

                        mat(i, j) = mat(i, j) + cmplxint

                        ! two electron integral : (ir, is | is, ia)
                        i2r = inttwr(ir, is, is, ia)
                        i2i = inttwi(ir, is, is, ia)
                        cmplxint = CMPLX(i2r, i2i, 16)

                        mat(i, j) = mat(i, j) - cmplxint
                    End do      !l0

                    Do l0 = 1, nelec
                        inds = oc(l0)
                        is = inds + ninact

                        ! two electron integral : (ir, ia | is, is)
                        i2r = inttwr(ir, ia, is, is)
                        i2i = inttwi(ir, ia, is, is)
                        cmplxint = CMPLX(i2r, i2i, 16)

                        mat(i, j) = mat(i, j) + cmplxint

                        ! two electron integral : (ir, is | is, ia)
                        i2r = inttwr(ir, is, is, ia)
                        i2i = inttwi(ir, is, is, ia)
                        cmplxint = CMPLX(i2r, i2i, 16)

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
!                 write(*,*)'ir,indr,is,inds',ir,indr,is,inds

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
                            !    if (rank == 0) write (*, '(A,I5,A,I5)') 'Noda ijorder 2diff i:', i, ' j:', j
                            if (mod(phase1 + phase2, 2) == 0) phase = 1.0d+00
                            if (mod(phase1 + phase2, 2) == 1) phase = -1.0d+00

                            ! two electron integral : (ir, ia | is, ib)
                            i2r = inttwr(ir, ia, is, ib)
                            i2i = inttwi(ir, ia, is, ib)
                            cmplxint = CMPLX(i2r, i2i, 16)

                            mat(i, j) = cmplxint

                            ! two electron integral : (ir, ib | is, ia)
                            i2r = inttwr(ir, ib, is, ia)
                            i2i = inttwi(ir, ib, is, ia)
                            cmplxint = CMPLX(i2r, i2i, 16)

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
    if (rank == 0) then ! Process limits for output
        write (*, '(A,I4)') 'end casmat', rank
        write (*, '(A,I4)') 'Reduce mat(:,:)', rank
    end if
#ifdef HAVE_MPI
#ifdef BIG_MAT
    ! If only the master rank execute casmat.f90, allreduce is not neccesarry.
#else
    call MPI_Allreduce(MPI_IN_PLACE, mat(1, 1), ndet**2, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
#endif
1000 end subroutine casmat

! ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    SUBROUTINE idetr(iidet, j)

! ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!        use four_caspt2_module

!        Implicit NONE

!        integer, intent(in)  :: iidet
!        integer, intent(out) :: j
!        integer              :: i0

!        j = 0

!        Do i0 = 1, ndet
!            If (idet(i0) == iidet) then
!                j = i0
!            End if
!        End do

!    END SUBROUTINE idetr
