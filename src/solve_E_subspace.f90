SUBROUTINE solve_E_subspace(e0)

    use dcaspt2_restart_file, only: get_subspace_idx
    use module_blas, only: gemv, gemm
    use module_ulambda_s_half, only: ulambda_s_half
    use module_global_variables
    use module_realonly, only: realonly
    use module_time
    implicit none
    real(8), intent(in) :: e0
    integer :: subspace_idx

    subspace_idx = get_subspace_idx('E')
    if (realonly%is_realonly()) then
        call solve_E_subspace_real()
    else
        call solve_E_subspace_complex()
    end if
    e2all = e2all + e2_subspace(subspace_idx)
    sumc2 = sumc2 + sumc2_subspace(subspace_idx)
contains
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_E_subspace_complex()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx, convert_secondary_to_global_idx

        Implicit NONE

        integer :: dimn, dimm, dammy
        real(8), allocatable  :: wsnew(:), ws(:), wb(:)
        real(8)               :: e2(2*nsymrpa), alpha, e
        complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
        complex*16, allocatable  :: bc(:, :)
        complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :), vc(:), vc1(:)
        integer :: j, i, syma, symb, isym, indt(1:nact)
        integer :: ia, it, ij, ii, ja, jt, jj, ji
        integer :: i0
        integer, allocatable     :: ia0(:), ii0(:), ij0(:), iaij(:, :, :)
        integer :: naij

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE E IS NOW CALCULATED
!
!     EtiEaj|0>
!
!   DRAS1 =-2   DRAS2 = +1   DRAS3 = +1
!
!   i > j
!
!  S(ukbl,tiaj) = d(ki) d(lj) d(ba) [d(ut) - <0|Etu|0>]  <= S(u,t)
!                                   ~~~~~~~~~~~~~~~~~~~
!  S(u,t) = d(ut) - <0|Etu|0>
!
!  B(u,t) = Siguma_w [eps(w){d(tu)<0|Eww|0>-<0|EtuEww|0>}] + S(u,t)eps(t)
!
!         = - Siguma_w [eps(w)<0|EtuEww|0>] + d(tu)e0 + S(u,t)eps(t)
!
!  alpha(i,j,a) = -eps(i) - eps(j) + eps(a) - e0
!
!  where
!
!  e0 = Siguma_w [eps(w)<0|Eww|0>] (<== calculated as e0 in calce0.f)
!
!  V(t,ija)   =[SIGUMA_p:active <0|Ept|0>{(ai|pj) - (aj|pi)}] + (aj|ti) - (ai|tj)
!
!  E2 = SIGUMA_iab, dimm |V1(t,ija)|^2|/{(alpha(ija) + wb(t)}
!
        e2 = 0.0d+00
        dimn = 0
        syma = 0
        indt = 0
        if (debug .and. rank == 0) print *, 'ENTER solve E part'
        if (rank == 0) print '(10A)', '  '
        if (rank == 0) print '(10A)', ' e2e(isym)'

! (ninact*(ninact-1))/2 means the number of (ii,ij) pairs (ii>ij)
        i0 = nsec*(ninact*(ninact - 1))/2
        naij = i0
        Allocate (iaij(nsec, ninact, ninact))
        iaij = 0
        Allocate (ia0(naij))
        Allocate (ii0(naij))
        Allocate (ij0(naij))

        i0 = 0
        Do ii = 1, ninact
            Do ij = 1, ii - 1                ! i > j
                Do ia = 1, nsec
                    i0 = i0 + 1
                    iaij(ia, ii, ij) = i0
                    iaij(ia, ij, ii) = i0
                    ia0(i0) = convert_secondary_to_global_idx(ia) ! secondary
                    ii0(i0) = ii ! inactive
                    ij0(i0) = ij ! inactive
                End do
            End do
        End do

        Allocate (v(naij, nact))
        v = 0.0d+00
        Call vEmat_complex(naij, iaij, v)

        Do isym = 1, nsymrpa

            dimn = 0
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                if (irpamo(jt) == isym) then
                    dimn = dimn + 1
                    indt(dimn) = it
                End if
            End do

            if (debug .and. rank == 0) print *, 'isym, dimn', isym, dimn
            If (dimn == 0) cycle ! Go to the next isym

            Allocate (sc(dimn, dimn))
            sc = 0.0d+00            ! sc N*N
            Call sEmat_complex(dimn, indt(1:dimn), sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            Allocate (ws(dimn))
            Allocate (sc0(dimn, dimn))
            sc0 = sc
            Call cdiag(sc, dimn, dimm, ws, smat_lin_dep_threshold)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (debug .and. rank == 0) print *, 'after E subspace S matrix cdiag, new dimension is', dimm
            If (dimm == 0) then
                deallocate (sc0)
                deallocate (sc)
                deallocate (ws)
                cycle ! Go to the next isym
            End if

            If (debug) then

                if (debug .and. rank == 0) print *, 'Check whether U*SU is diagonal'
                Call checkdgc(dimn, sc0, sc, ws)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (debug .and. rank == 0) print *, 'Check whether U*SU is diagonal END'
            End if

            Allocate (bc(dimn, dimn))                                 ! bc N*N
            bc = 0.0d+00
            Call bEmat_complex(dimn, sc0, indt(1:dimn), bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            deallocate (sc0)

            Allocate (uc(dimn, dimm))                                 ! uc N*M
            Allocate (wsnew(dimm))                                  ! wnew M
            uc(:, :) = 0.0d+00
            wsnew(:) = 0.0d+00
            Call ccutoff(sc, ws, dimn, dimm, smat_lin_dep_threshold, uc, wsnew)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            deallocate (ws)
            deallocate (sc)
            Call ulambda_s_half(uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            deallocate (wsnew)

            Allocate (bc0(dimm, dimn))                       ! bc0 M*N
            call gemm(transpose(DCONJG(uc)), bc, bc0)
            Allocate (bc1(dimm, dimm))                      ! bc1 M*M
            call gemm(bc0, uc, bc1)

            If (debug) then

                if (rank == 0) then
                    print *, 'Check whether bc1 is hermite or not'
                    Do i = 1, dimm
                        Do j = i, dimm
                            if (ABS(bc1(i, j) - DCONJG(bc1(j, i))) > 1.0d-6) then
                                print '(2I4,2E15.7)', i, j, bc1(i, j) - bc1(j, i)
                            End if
                        End do
                    End do
                    print *, 'Check whether bc1 is hermite or not END'
                end if

            End if

            deallocate (bc)
            deallocate (bc0)

            Allocate (wb(dimm))

            if (debug .and. rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'
            Allocate (bc0(dimm, dimm))
            bc0 = bc1
            Call cdiag(bc1, dimm, dammy, wb, bmat_no_cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            If (debug) then
                if (debug .and. rank == 0) print *, 'Check whether bc is really diagonalized or not'
                Call checkdgc(dimm, bc0, bc1, wb)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (debug .and. rank == 0) print *, 'Check whether bc is really diagonalized or not END'
            End if

            deallocate (bc0)

            if (debug .and. rank == 0) print *, 'bC1 matrix is diagonalized!'
            e2 = 0.0d+00

            Do i0 = 1, naij
                ja = ia0(i0)
                ji = ii0(i0)
                jj = ij0(i0)

!     EtiEaj|0>

                if (nsymrpa /= 1) then
                    syma = MULTB_D(irpamo(ja), irpamo(jj))
                    symb = MULTB_D(isym, irpamo(ji))
                    syma = MULTB_S(symb, syma)
                end if
                If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then

                    Allocate (vc(dimn))

                    Do it = 1, dimn
                        vc(it) = v(i0, indt(it))
                    End do

                    Allocate (vc1(dimm))
                    vc1 = 0.0d+00

                    vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn, 1:dimm))), vc(1:dimn))
                    Deallocate (vc)

                    alpha = +eps(ja) - eps(ji) - eps(jj) - e0 + eshift  ! For Level Shift (2007/2/9)

                    vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm, 1:dimm))), vc1(1:dimm))

                    Do j = 1, dimm
                        e = DBLE(DCONJG(vc1(j))*vc1(j)/(alpha + wb(j)))
                        sumc2_subspace(subspace_idx) = sumc2_subspace(subspace_idx) + e/(alpha + wb(j))
                        e2(isym) = e2(isym) - e
                    End do

                    deallocate (vc1)

                End if

            End do

            deallocate (uc)
            deallocate (wb)
            Deallocate (bc1)

            if (rank == 0) print '(" e2e(",I3,") = ",E25.15," a.u.")', isym, e2(isym)
            e2_subspace(subspace_idx) = e2_subspace(subspace_idx) + e2(isym)
        End do

        if (rank == 0) then
            print '(" e2e      = ",E25.15," a.u.")', e2_subspace(subspace_idx)
            print '(" sumc2,e  = ",E25.15)', sumc2_subspace(subspace_idx)
        end if

        deallocate (iaij)
        deallocate (ia0)
        deallocate (ii0)
        deallocate (ij0)
        deallocate (v)

        continue
        if (debug .and. rank == 0) print *, 'end solve_E_subspace'
    end subroutine solve_E_subspace_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE sEmat_complex(dimn, indt, sc) ! Assume C1 molecule, overlap matrix S in space E

! S(u,t) = d(ut) - <0|Etu|0>
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)      :: dimn, indt(dimn)
        complex*16, intent(out)  :: sc(dimn, dimn)

        real(8)  ::a, b

        integer :: it, iu
        integer :: i, j

        if (debug .and. rank == 0) print *, 'Start E subspace S matrix'
        sc = 0.0d+00

!$OMP parallel do schedule(dynamic,1) private(iu,j,it,a,b)
        Do i = rank + 1, dimn, nprocs
            iu = indt(i)

            Do j = i, dimn
                it = indt(j)
                a = 0.0d+0
                b = 0.0d+0

                Call dim1_density(it, iu, a, b)

                If (iu == it) then
                    sc(i, j) = 1 - DCMPLX(a, b)
                Else
                    sc(i, j) = -DCMPLX(a, b)
                End if

                sc(j, i) = DCONJG(sc(i, j))

            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=sc)
#endif
        if (debug .and. rank == 0) print *, 'E subspace S matrix is obtained normally'
    End subroutine sEmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE bEmat_complex(dimn, sc, indt, bc) ! Assume C1 molecule, overlap matrix B in space E
!
!
!  S(u,t) = d(ut) - <0|Etu|0>
!
!  B(u,t) = Siguma_w [eps(w){d(tu)<0|Eww|0>-<0|EtuEww|0>}] + S(u,t)eps(t)
!
!         = - Siguma_w [eps(w)<0|EtuEww|0>] + d(tu)e0 + S(u,t)eps(t)
!
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer :: it, iu, iw, jt, ju, jw
        integer :: i, j

        integer, intent(in) :: dimn, indt(dimn)
        complex*16, intent(in)  :: sc(dimn, dimn)
        complex*16, intent(out) :: bc(dimn, dimn)

        real(8)              :: denr, deni
        complex*16          :: den

        bc(:, :) = 0.0d+00

        if (debug .and. rank == 0) print *, 'Start E subspace B matrix'

!$OMP parallel do schedule(dynamic,1) private(iu,ju,j,it,jt,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs
            iu = indt(i)
            ju = convert_active_to_global_idx(iu)

            Do j = i, dimn
                it = indt(j)
                jt = convert_active_to_global_idx(it)

                Do iw = 1, nact
                    jw = convert_active_to_global_idx(iw)

!         = - Siguma_w [eps(w)<0|EtuEww|0>] + d(tu)e0 + S(u,t)eps(t)

                    Call dim2_density(it, iu, iw, iw, denr, deni)
                    den = DCMPLX(denr, deni)
                    bc(i, j) = bc(i, j) - den*eps(jw)

                End do

                if (it == iu) bc(i, j) = bc(i, j) + e0

                bc(i, j) = bc(i, j) + sc(i, j)*eps(jt)

                bc(j, i) = DCONJG(bc(i, j))

            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call reduce_wrapper(mat=bc, root_rank=0)
#endif
        if (debug .and. rank == 0) print *, 'E subspace B matrix is obtained normally'
    End subroutine bEmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE vEmat_complex(naij, iaij, v)
!
!  V(t,ija)   =[SIGUMA_p:active <0|Ept|0>{(ai|pj) - (aj|pi)}] + (aj|ti) - (ai|tj)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_file_manager
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)     :: naij, iaij(nsec, ninact, ninact)

        complex*16, intent(out) :: v(naij, nact)

        real(8)                  :: dr, di
        complex*16              :: cint2, dens

        integer :: i, j, k, l, taij
        integer :: it, iostat, unit_int2
        logical :: is_end_of_file

        if (debug .and. rank == 0) print *, 'Start E subspace V matrix'
        v = 0.0d+00

!  V(t,ija)   =[SIGUMA_p:active <0|Ept|0>{(ai|pj) - (aj|pi)}] - (ai|tj) + (aj|ti)   i > j

        call open_unformatted_file(unit=unit_int2, file=eint, status='old', optional_action='read') !  (31|21) stored
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2
            call check_iostat(iostat=iostat, file=eint, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            if (j == l) cycle ! Read the next 2-integral if j equal to l

            taij = iaij(i, j, l)

            if (j < l) then
                cint2 = -1.0d+00*cint2
            end if

            v(taij, k) = v(taij, k) - cint2

!$OMP parallel do schedule(dynamic,1) private(it,dr,di,dens)
            Do it = 1, nact
                Call dim1_density(it, k, dr, di)          ! k corresponds to p in above formula
                dens = DCMPLX(dr, di)
                v(taij, it) = v(taij, it) + cint2*dens
            End do
!$OMP end parallel do

            if (j < l) then
                cint2 = -1.0d+00*cint2          ! data cint2 becomes initial values!
            end if

        end do
        close (unit_int2)
        if (debug .and. rank == 0) print *, 'reading Eint2 is over'

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
#endif
        if (debug .and. rank == 0) print *, 'E subspace V matrix is obtained normally'
    end subroutine vEmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_E_subspace_real()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx, convert_secondary_to_global_idx

        Implicit NONE

        integer :: dimn, dimm, dammy
        real(8), allocatable  :: wsnew(:), ws(:), wb(:)
        real(8)               :: e2(2*nsymrpa), alpha, e
        real(8), allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
        real(8), allocatable  :: bc(:, :)
        real(8), allocatable  :: bc0(:, :), bc1(:, :), v(:, :), vc(:), vc1(:)
        integer :: j, i, syma, symb, isym, indt(1:nact)
        integer :: ia, it, ij, ii, ja, jt, jj, ji
        integer :: i0
        integer, allocatable     :: ia0(:), ii0(:), ij0(:), iaij(:, :, :)
        integer :: naij

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE E IS NOW CALCULATED
!
!     EtiEaj|0>
!
!   DRAS1 =-2   DRAS2 = +1   DRAS3 = +1
!
!   i > j
!
!  S(ukbl,tiaj) = d(ki) d(lj) d(ba) [d(ut) - <0|Etu|0>]  <= S(u,t)
!                                   ~~~~~~~~~~~~~~~~~~~
!  S(u,t) = d(ut) - <0|Etu|0>
!
!  B(u,t) = Siguma_w [eps(w){d(tu)<0|Eww|0>-<0|EtuEww|0>}] + S(u,t)eps(t)
!
!         = - Siguma_w [eps(w)<0|EtuEww|0>] + d(tu)e0 + S(u,t)eps(t)
!
!  alpha(i,j,a) = -eps(i) - eps(j) + eps(a) - e0
!
!  where
!
!  e0 = Siguma_w [eps(w)<0|Eww|0>] (<== calculated as e0 in calce0.f)
!
!  V(t,ija)   =[SIGUMA_p:active <0|Ept|0>{(ai|pj) - (aj|pi)}] + (aj|ti) - (ai|tj)
!
!  E2 = SIGUMA_iab, dimm |V1(t,ija)|^2|/{(alpha(ija) + wb(t)}
!
        e2 = 0.0d+00
        dimn = 0
        syma = 0
        indt = 0
        if (debug .and. rank == 0) print *, 'ENTER solve E part'
        if (rank == 0) print '(10A)', '  '
        if (rank == 0) print '(10A)', ' e2e(isym)'

! (ninact*(ninact-1))/2 means the number of (ii,ij) pairs (ii>ij)
        i0 = nsec*(ninact*(ninact - 1))/2
        naij = i0
        Allocate (iaij(nsec, ninact, ninact))
        iaij = 0
        Allocate (ia0(naij))
        Allocate (ii0(naij))
        Allocate (ij0(naij))

        i0 = 0
        Do ii = 1, ninact
            Do ij = 1, ii - 1                ! i > j
                Do ia = 1, nsec
                    i0 = i0 + 1
                    iaij(ia, ii, ij) = i0
                    iaij(ia, ij, ii) = i0
                    ia0(i0) = convert_secondary_to_global_idx(ia) ! secondary
                    ii0(i0) = ii ! inactive
                    ij0(i0) = ij ! inactive
                End do
            End do
        End do

        Allocate (v(naij, nact))
        v = 0.0d+00
        Call vEmat_real(naij, iaij, v)

        Do isym = 1, nsymrpa

            dimn = 0
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                if (irpamo(jt) == isym) then
                    dimn = dimn + 1
                    indt(dimn) = it
                End if
            End do

            if (debug .and. rank == 0) print *, 'isym, dimn', isym, dimn
            If (dimn == 0) cycle ! Go to the next isym

            Allocate (sc(dimn, dimn))
            sc = 0.0d+00            ! sc N*N
            Call sEmat_real(dimn, indt(1:dimn), sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            Allocate (ws(dimn))
            Allocate (sc0(dimn, dimn))
            sc0 = sc
            Call rdiag(sc, dimn, dimm, ws, smat_lin_dep_threshold)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (debug .and. rank == 0) print *, 'after E subspace S matrix rdiag, new dimension is', dimm
            If (dimm == 0) then
                deallocate (sc0)
                deallocate (sc)
                deallocate (ws)
                cycle ! Go to the next isym
            End if

            Allocate (bc(dimn, dimn))                                 ! bc N*N
            bc = 0.0d+00
            Call bEmat_real(dimn, sc0, indt(1:dimn), bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            deallocate (sc0)

            Allocate (uc(dimn, dimm))                                 ! uc N*M
            Allocate (wsnew(dimm))                                  ! wnew M
            uc(:, :) = 0.0d+00
            wsnew(:) = 0.0d+00
            Call rcutoff(sc, ws, dimn, dimm, smat_lin_dep_threshold, uc, wsnew)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            deallocate (ws)
            deallocate (sc)
            Call ulambda_s_half(uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            deallocate (wsnew)

            Allocate (bc0(dimm, dimn))                       ! bc0 M*N
            call gemm(transpose(uc), bc, bc0)
            Allocate (bc1(dimm, dimm))                      ! bc1 M*M
            call gemm(bc0, uc, bc1)

            If (debug) then

                if (rank == 0) then
                    print *, 'Check whether bc1 is hermite or not'
                    Do i = 1, dimm
                        Do j = i, dimm
                            if (ABS(bc1(i, j) - bc1(j, i)) > 1.0d-6) then
                                print '(2I4,2E15.7)', i, j, bc1(i, j) - bc1(j, i)
                            End if
                        End do
                    End do
                    print *, 'Check whether bc1 is hermite or not END'
                end if

            End if

            deallocate (bc)
            deallocate (bc0)

            Allocate (wb(dimm))

            if (debug .and. rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'
            Allocate (bc0(dimm, dimm))
            bc0 = bc1
            Call rdiag(bc1, dimm, dammy, wb, bmat_no_cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            deallocate (bc0)

            if (debug .and. rank == 0) print *, 'bC1 matrix is diagonalized!'
            e2 = 0.0d+00

            Do i0 = 1, naij
                ja = ia0(i0)
                ji = ii0(i0)
                jj = ij0(i0)

!     EtiEaj|0>

                if (nsymrpa /= 1) then
                    syma = MULTB_D(irpamo(ja), irpamo(jj))
                    symb = MULTB_D(isym, irpamo(ji))
                    syma = MULTB_S(symb, syma)
                end if
                If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then

                    Allocate (vc(dimn))

                    Do it = 1, dimn
                        vc(it) = v(i0, indt(it))
                    End do

                    Allocate (vc1(dimm))
                    vc1 = 0.0d+00

                    vc1(1:dimm) = MATMUL(TRANSPOSE(uc(1:dimn, 1:dimm)), vc(1:dimn))
                    Deallocate (vc)

                    alpha = +eps(ja) - eps(ji) - eps(jj) - e0 + eshift  ! For Level Shift (2007/2/9)

                    vc1(1:dimm) = MATMUL(TRANSPOSE(bc1(1:dimm, 1:dimm)), vc1(1:dimm))

                    Do j = 1, dimm
                        e = DBLE(vc1(j)**2.0d+00/(alpha + wb(j)))
                        sumc2_subspace(subspace_idx) = sumc2_subspace(subspace_idx) + e/(alpha + wb(j))
                        e2(isym) = e2(isym) - e
                    End do

                    deallocate (vc1)

                End if

            End do

            deallocate (uc)
            deallocate (wb)
            Deallocate (bc1)

            if (rank == 0) print '(" e2e(",I3,") = ",E25.15," a.u.")', isym, e2(isym)
            e2_subspace(subspace_idx) = e2_subspace(subspace_idx) + e2(isym)
        End do

        if (rank == 0) then
            print '(" e2e      = ",E25.15," a.u.")', e2_subspace(subspace_idx)
            print '(" sumc2,e  = ",E25.15)', sumc2_subspace(subspace_idx)
        end if

        deallocate (iaij)
        deallocate (ia0)
        deallocate (ii0)
        deallocate (ij0)
        deallocate (v)

        continue
        if (debug .and. rank == 0) print *, 'end solve_E_subspace'
    end subroutine solve_E_subspace_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE sEmat_real(dimn, indt, sc) ! Assume C1 molecule, overlap matrix S in space E

! S(u,t) = d(ut) - <0|Etu|0>
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)      :: dimn, indt(dimn)
        real(8), intent(out)  :: sc(dimn, dimn)

        real(8)  ::a, b

        integer :: it, iu
        integer :: i, j

        if (debug .and. rank == 0) print *, 'Start E subspace S matrix'
        sc = 0.0d+00

!$OMP parallel do schedule(dynamic,1) private(iu,j,it,a,b)
        Do i = rank + 1, dimn, nprocs
            iu = indt(i)

            Do j = i, dimn
                it = indt(j)
                a = 0.0d+0
                b = 0.0d+0

                Call dim1_density(it, iu, a, b)

                If (iu == it) then
                    sc(i, j) = 1 - a
                Else
                    sc(i, j) = -a
                End if

                sc(j, i) = sc(i, j)

            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=sc)
#endif
        if (debug .and. rank == 0) print *, 'E subspace S matrix is obtained normally'
    End subroutine sEmat_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE bEmat_real(dimn, sc, indt, bc) ! Assume C1 molecule, overlap matrix B in space E
!
!
!  S(u,t) = d(ut) - <0|Etu|0>
!
!  B(u,t) = Siguma_w [eps(w){d(tu)<0|Eww|0>-<0|EtuEww|0>}] + S(u,t)eps(t)
!
!         = - Siguma_w [eps(w)<0|EtuEww|0>] + d(tu)e0 + S(u,t)eps(t)
!
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer :: it, iu, iw, jt, ju, jw
        integer :: i, j

        integer, intent(in) :: dimn, indt(dimn)
        real(8), intent(in)  :: sc(dimn, dimn)
        real(8), intent(out) :: bc(dimn, dimn)

        real(8)              :: denr, deni
        real(8)          :: den

        bc(:, :) = 0.0d+00

        if (debug .and. rank == 0) print *, 'Start E subspace B matrix'

!$OMP parallel do schedule(dynamic,1) private(iu,ju,j,it,jt,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs
            iu = indt(i)
            ju = convert_active_to_global_idx(iu)

            Do j = i, dimn
                it = indt(j)
                jt = convert_active_to_global_idx(it)

                Do iw = 1, nact
                    jw = convert_active_to_global_idx(iw)

!         = - Siguma_w [eps(w)<0|EtuEww|0>] + d(tu)e0 + S(u,t)eps(t)

                    Call dim2_density(it, iu, iw, iw, denr, deni)
                    den = denr
                    bc(i, j) = bc(i, j) - den*eps(jw)

                End do

                if (it == iu) bc(i, j) = bc(i, j) + e0

                bc(i, j) = bc(i, j) + sc(i, j)*eps(jt)

                bc(j, i) = bc(i, j)

            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call reduce_wrapper(mat=bc, root_rank=0)
#endif
        if (debug .and. rank == 0) print *, 'E subspace B matrix is obtained normally'
    End subroutine bEmat_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE vEmat_real(naij, iaij, v)
!
!  V(t,ija)   =[SIGUMA_p:active <0|Ept|0>{(ai|pj) - (aj|pi)}] + (aj|ti) - (ai|tj)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_file_manager
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)     :: naij, iaij(nsec, ninact, ninact)

        real(8), intent(out) :: v(naij, nact)

        real(8)                  :: dr, di
        real(8)              :: cint2, dens

        integer :: i, j, k, l, taij
        integer :: it, iostat, unit_int2
        logical :: is_end_of_file

        if (debug .and. rank == 0) print *, 'Start E subspace V matrix'
        v = 0.0d+00

!  V(t,ija)   =[SIGUMA_p:active <0|Ept|0>{(ai|pj) - (aj|pi)}] - (ai|tj) + (aj|ti)   i > j

        call open_unformatted_file(unit=unit_int2, file=eint, status='old', optional_action='read') !  (31|21) stored
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2
            call check_iostat(iostat=iostat, file=eint, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            if (j == l) cycle ! Read the next 2-integral if j equal to l

            taij = iaij(i, j, l)

            if (j < l) then
                cint2 = -1.0d+00*cint2
            end if

            v(taij, k) = v(taij, k) - cint2

!$OMP parallel do schedule(dynamic,1) private(it,dr,di,dens)
            Do it = 1, nact
                Call dim1_density(it, k, dr, di)          ! k corresponds to p in above formula
                dens = dr
                v(taij, it) = v(taij, it) + cint2*dens
            End do
!$OMP end parallel do

            if (j < l) then
                cint2 = -1.0d+00*cint2          ! data cint2 becomes initial values!
            end if

        end do
        close (unit_int2)
        if (debug .and. rank == 0) print *, 'reading Eint2 is over'

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
#endif
        if (debug .and. rank == 0) print *, 'E subspace V matrix is obtained normally'
    end subroutine vEmat_real

end subroutine solve_E_subspace
