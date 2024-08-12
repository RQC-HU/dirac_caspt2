SUBROUTINE solve_H_subspace(e0)

    use dcaspt2_restart_file, only: get_subspace_idx
    use module_global_variables
    use module_realonly, only: realonly
    implicit none
    real(8), intent(in) :: e0
    integer :: subspace_idx

    subspace_idx = get_subspace_idx('H')
    if (realonly%is_realonly()) then
        call solve_H_subspace_real()
    else
        call solve_H_subspace_complex()
    end if

contains
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_H_subspace_complex()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_file_manager, only: open_unformatted_file, check_iostat
        use module_global_variables
        use module_index_utils, only: convert_secondary_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        Integer                 :: ia, ib, ii, ij, syma, symb, i, j, k, l
        Integer                 :: i0, j0, tab, nab, tij, nij, iostat, unit_int2
        Integer, allocatable    :: ia0(:), ib0(:), ii0(:), ij0(:), iab(:, :), iij(:, :)
        complex*16              :: cint2
        complex*16, allocatable :: v(:, :)
        real(8)                  :: e
        logical                 :: is_end_of_file

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE H IS NOW CALCULATED
!
!     EaiEbj|0>         a > b, i > j
!
!   DRAS1 =-2   DRAS2 = 0    DRAS3 = +2
!
!
!  S(ckdl,aibj) = d(ac)d(bd)d(lj)d(ik)
!
!  (H0-E0S)ckdl,aibj = d(ac)d(bd)d(lj)d(ik)(eps(a)+eps(b)-eps(i)-eps(j)) = e(a,b,i,j)
!
!  V(aibj)   = (ai|bj) - (aj|bi)
!
! E2h = V(aibj)/e(a,b,i,j)
        if (debug .and. rank == 0) print *, 'ENTER solve H part'
        if (rank == 0) print '(10A)', '  '

        e = 0.0d+00

        i0 = 0
        Do ia = global_sec_start, global_sec_end
            Do ib = global_sec_start, ia - 1
                i0 = i0 + 1
            End do
        End do

        nab = i0

        Allocate (iab(nsec, nsec))
        Allocate (ia0(nab))
        Allocate (ib0(nab))
        iab = 0

        i0 = 0
        Do ia = 1, nsec
            Do ib = 1, ia - 1
                i0 = i0 + 1
                iab(ia, ib) = i0
                iab(ib, ia) = i0
                ia0(i0) = convert_secondary_to_global_idx(ia) ! secondary
                ib0(i0) = convert_secondary_to_global_idx(ib) ! secondary
            End do
        End do

        i0 = 0
        Do ii = 1, ninact
            Do ij = 1, ii - 1
                i0 = i0 + 1
            End do
        End do

        nij = i0
        Allocate (iij(1:ninact, 1:ninact))
        Allocate (ii0(nij))
        Allocate (ij0(nij))
        iij = 0

        i0 = 0
        Do ii = 1, ninact
            Do ij = 1, ii - 1
                i0 = i0 + 1
                iij(ii, ij) = i0
                iij(ij, ii) = i0
                ii0(i0) = ii
                ij0(i0) = ij
            End do
        End do

        Allocate (v(nab, nij))
        v = 0.0d+00

        call open_unformatted_file(unit=unit_int2, file=hint, status='old', optional_action='read')
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2
            call check_iostat(iostat=iostat, file=hint, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            if (i <= k .or. j == l) cycle ! Read the next line if i <= k or j == l

            tab = iab(i, k)
            tij = iij(j, l)

            if (i > k .and. j > l) then
                v(tab, tij) = v(tab, tij) + cint2

            elseif (i > k .and. j < l) then
                v(tab, tij) = v(tab, tij) - cint2

            elseif (i < k .and. j > l) then               ! (kl|ij)  l > j + ; l < j -
                v(tab, tij) = v(tab, tij) - cint2

            elseif (i < k .and. j < l) then
                v(tab, tij) = v(tab, tij) + cint2

            end if

        end do
        close (unit_int2)
        if (debug .and. rank == 0) print *, 'reading Hint is over'

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
#endif

        Do i0 = 1, nab
            ia = ia0(i0)
            ib = ib0(i0)

!     EaiEbj|0>         a > b, i > j

            Do j0 = 1, nij
                ii = ii0(j0)
                ij = ij0(j0)
                syma = MULTB_D(irpamo(ia), irpamo(ii))
                symb = MULTB_D(irpamo(ib), irpamo(ij))
                syma = MULTB_S(syma, symb)
                if (syma == 1) then

                    e = eps(ia) + eps(ib) - eps(ii) - eps(ij) + eshift  ! For Level Shift (2007/2/9)

                    coeff1 = v(i0, j0)/e
                    sumc2_subspace(subspace_idx) = sumc2_subspace(subspace_idx) + ABS(coeff1)**2

                    e2_subspace(subspace_idx) = e2_subspace(subspace_idx) - DBLE(DCONJG(v(i0, j0))*v(i0, j0)/e)
                end if
            End do
        End do

        if (rank == 0) then
            print '(" e2h      = ",E25.15," a.u.")', e2_subspace(subspace_idx)
            print '(" sumc2,h  = ",E25.15)', sumc2_subspace(subspace_idx)
        end if

        deallocate (v)
        deallocate (iab)
        deallocate (ia0)
        deallocate (ib0)
        deallocate (iij)
        deallocate (ii0)
        deallocate (ij0)

        if (debug .and. rank == 0) print *, 'end solve_H_subspace'
    End SUBROUTINE solve_H_subspace_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_H_subspace_real()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_file_manager, only: open_unformatted_file, check_iostat
        use module_global_variables
        use module_index_utils, only: convert_secondary_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        Integer                 :: ia, ib, ii, ij, syma, symb, i, j, k, l
        Integer                 :: i0, j0, tab, nab, tij, nij, iostat, unit_int2
        Integer, allocatable    :: ia0(:), ib0(:), ii0(:), ij0(:), iab(:, :), iij(:, :)
        real(8)              :: cint2
        real(8), allocatable :: v(:, :)
        real(8)                  :: e
        logical                 :: is_end_of_file

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE H IS NOW CALCULATED
!
!     EaiEbj|0>         a > b, i > j
!
!   DRAS1 =-2   DRAS2 = 0    DRAS3 = +2
!
!
!  S(ckdl,aibj) = d(ac)d(bd)d(lj)d(ik)
!
!  (H0-E0S)ckdl,aibj = d(ac)d(bd)d(lj)d(ik)(eps(a)+eps(b)-eps(i)-eps(j)) = e(a,b,i,j)
!
!  V(aibj)   = (ai|bj) - (aj|bi)
!
! E2h = V(aibj)/e(a,b,i,j)
        if (debug .and. rank == 0) print *, 'ENTER solve H part'
        if (rank == 0) print '(10A)', '  '

        e = 0.0d+00

        i0 = 0
        Do ia = global_sec_start, global_sec_end
            Do ib = global_sec_start, ia - 1
                i0 = i0 + 1
            End do
        End do

        nab = i0

        Allocate (iab(nsec, nsec))
        Allocate (ia0(nab))
        Allocate (ib0(nab))
        iab = 0

        i0 = 0
        Do ia = 1, nsec
            Do ib = 1, ia - 1
                i0 = i0 + 1
                iab(ia, ib) = i0
                iab(ib, ia) = i0
                ia0(i0) = convert_secondary_to_global_idx(ia) ! secondary
                ib0(i0) = convert_secondary_to_global_idx(ib) ! secondary
            End do
        End do

        i0 = 0
        Do ii = 1, ninact
            Do ij = 1, ii - 1
                i0 = i0 + 1
            End do
        End do

        nij = i0
        Allocate (iij(1:ninact, 1:ninact))
        Allocate (ii0(nij))
        Allocate (ij0(nij))
        iij = 0

        i0 = 0
        Do ii = 1, ninact
            Do ij = 1, ii - 1
                i0 = i0 + 1
                iij(ii, ij) = i0
                iij(ij, ii) = i0
                ii0(i0) = ii
                ij0(i0) = ij
            End do
        End do

        Allocate (v(nab, nij))
        v = 0.0d+00

        call open_unformatted_file(unit=unit_int2, file=hint, status='old', optional_action='read')
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2
            call check_iostat(iostat=iostat, file=hint, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            if (i <= k .or. j == l) cycle ! Read the next line if i <= k or j == l

            tab = iab(i, k)
            tij = iij(j, l)

            if (i > k .and. j > l) then
                v(tab, tij) = v(tab, tij) + cint2

            elseif (i > k .and. j < l) then
                v(tab, tij) = v(tab, tij) - cint2

            elseif (i < k .and. j > l) then               ! (kl|ij)  l > j + ; l < j -
                v(tab, tij) = v(tab, tij) - cint2

            elseif (i < k .and. j < l) then
                v(tab, tij) = v(tab, tij) + cint2

            end if

        end do
        close (unit_int2)
        if (debug .and. rank == 0) print *, 'reading Hint is over'

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
#endif

        Do i0 = 1, nab
            ia = ia0(i0)
            ib = ib0(i0)

!     EaiEbj|0>         a > b, i > j

            Do j0 = 1, nij
                ii = ii0(j0)
                ij = ij0(j0)
                syma = MULTB_D(irpamo(ia), irpamo(ii))
                symb = MULTB_D(irpamo(ib), irpamo(ij))
                syma = MULTB_S(syma, symb)
                if (syma == 1) then

                    e = eps(ia) + eps(ib) - eps(ii) - eps(ij) + eshift  ! For Level Shift (2007/2/9)

                    coeff1 = v(i0, j0)/e
                    sumc2_subspace(subspace_idx) = sumc2_subspace(subspace_idx) + ABS(coeff1)**2

                    e2_subspace(subspace_idx) = e2_subspace(subspace_idx) - DBLE(v(i0, j0)**2.0d+00/e)
                end if
            End do
        End do

        if (rank == 0) then
            print '(" e2h      = ",E25.15," a.u.")', e2_subspace(subspace_idx)
            print '(" sumc2,h  = ",E25.15)', sumc2_subspace(subspace_idx)
        end if

        deallocate (v)
        deallocate (iab)
        deallocate (ia0)
        deallocate (ib0)
        deallocate (iij)
        deallocate (ii0)
        deallocate (ij0)

        if (debug .and. rank == 0) print *, 'end solve_H_subspace'
    End SUBROUTINE solve_H_subspace_real

end subroutine solve_H_subspace
