SUBROUTINE solve_F_subspace(e0)

    use, intrinsic :: iso_fortran_env, only: int64
    use dcaspt2_restart_file, only: get_subspace_idx
    use module_blas, only: gemv, gemm
    use module_ulambda_s_half, only: ulambda_s_half
    use module_global_variables
    use module_realonly, only: realonly
    use module_time
    implicit none
    real(8), intent(in) :: e0
    integer :: subspace_idx

    subspace_idx = get_subspace_idx('F')
    if (realonly%is_realonly()) then
        call solve_F_subspace_real()
    else
        call solve_F_subspace_complex()
    end if
    e2all = e2all + e2_subspace(subspace_idx)
    sumc2 = sumc2 + sumc2_subspace(subspace_idx)
contains
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_F_subspace_complex()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx, convert_secondary_to_global_idx

        Implicit NONE

        integer :: dimn, dimm, dammy
        integer, allocatable :: indsym(:, :)
        real(8), allocatable  :: wsnew(:), ws(:), wb(:)
        real(8)               :: e2(2*nsymrpa), alpha, e
        complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
        complex*16, allocatable  :: bc(:, :)
        complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :), vc(:), vc0(:), vc1(:)
        integer :: j, i, syma, isym, i0
        integer :: ia, it, ib, iu, ja, jt, jb, ju
        integer, allocatable     :: ia0(:), ib0(:), iab(:, :)
        integer                  :: nab

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE F IS NOW CALCULATED
!
!     EatEbu|0>
!
!   DRAS1 =0   DRAS2 = -2   DRAS3 = +2
!
!   a > b t > u ( c > d, v > x)
!
!  S(cvdx,atbu) = d(ac) d(bd)  [ <0|EvtExu|0> - d(xt)<0|Evu|0>]  <= S(vx,tu)
!                               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  S(vx, tu) = <0|EvtExu|0> - d(xt)<0|Evu|0>
!
!  B(vx, tu) = Siguma_w [eps(w){ <0|EvtExuEww|0> - d(xt)<0|EvuEww|0>}] + S(u,t){-eps(u)-eps(t)}
!
!  alpha(a, b) = + eps(a) + eps(b) - e0
!
!  where
!
!  e0 = Siguma_w [eps(w)<0|Eww|0>] (<== calculated as e0 in calce0.f)
!
!  V(tu,ab)   = SIGUMA_p,q:active <0|EtpEuq|0>(ap|bq) -  SIGUMA_p:active <0|Etp|0>(au|bp)
!
!  E2 = SIGUMA_iab,t:dimm |V1(t,ab)|^2|/{(alpha(ab) + wb(t)}
!
        e2 = 0.0d+00
        dimn = 0
        syma = 0
        if (debug .and. rank == 0) print *, 'ENTER solve F part'
        if (rank == 0) print '(10A)', '  '
        if (rank == 0) print '(10A)', ' e2f(isym)'

        i0 = 0
        Do ia = 1, nsec
            Do ib = 1, ia - 1
                i0 = i0 + 1
            End do
        End do

        nab = i0
        Allocate (iab(nsec, nsec))
        iab = 0
        Allocate (ia0(nab))
        Allocate (ib0(nab))

        i0 = 0
        Do ia = 1, nsec
            ja = convert_secondary_to_global_idx(ia)
            Do ib = 1, ia - 1
                jb = convert_secondary_to_global_idx(ib)
                i0 = i0 + 1
                iab(ia, ib) = i0
                iab(ib, ia) = i0
                ia0(i0) = convert_secondary_to_global_idx(ia) ! secondary
                ib0(i0) = convert_secondary_to_global_idx(ib) ! secondary
            End do
        End do

        Allocate (v(nab, nact, nact))
        v = 0.0d+00
        Call vFmat_complex(nab, iab, v)
        Do isym = 1, nsymrpa

            dimn = 0
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Do iu = 1, it - 1
                    ju = convert_active_to_global_idx(iu)

!     EatEbu|0>

                    if (nsymrpa /= 1) syma = MULTB_D(irpamo(ju) - (-1)**(mod(irpamo(ju), 2)), irpamo(jt))

                    if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                        dimn = dimn + 1
                    End if
                End do
            End do

            if (debug .and. rank == 0) print *, 'isym, dimn', isym, dimn
            If (dimn == 0) cycle ! Go to the next isym if dimn (dimention of matrix) is zero

            Allocate (indsym(2, dimn))

            dimn = 0
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Do iu = 1, it - 1
                    ju = convert_active_to_global_idx(iu)

!     EatEbu|0>

                    if (nsymrpa /= 1) syma = MULTB_D(irpamo(ju) - (-1)**(mod(irpamo(ju), 2)), irpamo(jt))

                    if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                        dimn = dimn + 1
                        indsym(1, dimn) = it
                        indsym(2, dimn) = iu
                    End if
                End do
            End do

            Allocate (sc(dimn, dimn))
            sc = 0.0d+00            ! sc N*N
            Call sFmat_complex(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            Allocate (ws(dimn))

            Allocate (sc0(dimn, dimn))
            sc0 = sc
            Call cdiag(sc, dimn, dimm, ws, smat_lin_dep_threshold)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (debug .and. rank == 0) print *, 'after F subspace S matrix cdiag, new dimension is', dimm
            If (dimm == 0) then
                deallocate (indsym)
                deallocate (sc0)
                deallocate (sc)
                deallocate (ws)
                cycle ! Go to the next isym if dimm (dimention of matrix) is zero
            End if

            If (debug) then

                if (debug .and. rank == 0) print *, 'Check whether U*SU is diagonal'
                Call checkdgc(dimn, sc0, sc, ws)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (debug .and. rank == 0) print *, 'Check whether U*SU is diagonal END'
            End if

            Allocate (bc(dimn, dimn))                                 ! bc N*N
            bc = 0.0d+00
            Call bFmat_complex(dimn, sc0, indsym, bc)
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

            Do i0 = 1, nab
                ja = ia0(i0)
                jb = ib0(i0)

!     EatEbu|0>
                if (nsymrpa /= 1) then
                    syma = MULTB_D(irpamo(ja), irpamo(jb) - (-1)**(mod(irpamo(jb), 2)))
                    syma = MULTB_S(syma, isym)
                end if
                If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then

                    Allocate (vc(dimn))
                    Do it = 1, dimn
                        vc(it) = v(i0, indsym(1, it), indsym(2, it))
                    End do
                    allocate (vc0(dimm)); Call memplus(KIND(vc0), SIZE(vc0), 2)
                    allocate (vc1(dimm)); Call memplus(KIND(vc1), SIZE(vc1), 2)
                    call gemv(transpose(DCONJG(uc)), vc, vc0)
                    Deallocate (vc)

                    alpha = +eps(ja) + eps(jb) - e0 + eshift  ! For Level Shift (2007/2/9)

                    call gemv(transpose(DCONJG(bc1)), vc0, vc1)

                    Do j = 1, dimm
                        e = (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                        sumc2_subspace(subspace_idx) = sumc2_subspace(subspace_idx) + e/(alpha + wb(j))
                        e2(isym) = e2(isym) - e
                    End do

                    Call memminus(KIND(vc0), SIZE(vc0), 2); Deallocate (vc0)
                    Call memminus(KIND(vc1), SIZE(vc1), 2); Deallocate (vc1)

                End if

            End do

            deallocate (indsym)
            deallocate (uc)
            deallocate (wb)
            Deallocate (bc1)

            if (rank == 0) print '(" e2f(",I3,") = ",E25.15," a.u.")', isym, e2(isym)
            e2_subspace(subspace_idx) = e2_subspace(subspace_idx) + e2(isym)
        End do

        if (rank == 0) then
            print '(" e2f      = ",E25.15," a.u.")', e2_subspace(subspace_idx)
            print '(" sumc2,f  = ",E25.15)', sumc2_subspace(subspace_idx)
        end if

        deallocate (iab)
        deallocate (ia0)
        deallocate (ib0)
        deallocate (v)

        if (debug .and. rank == 0) print *, 'end solve_F_subspace'
    end subroutine solve_F_subspace_complex
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE sFmat_complex(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space F

!  S(vx, tu) = <0|EvtExu|0> - d(xt)<0|Evu|0>
!
!     v > x, t > u

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)      :: dimn, indsym(2, dimn)
        complex*16, intent(out)  :: sc(dimn, dimn)
        real(8) :: a, b
        integer(kind=int64) :: it, iu, iv, ix
        integer :: i, j

        if (debug .and. rank == 0) print *, 'Start F subspace S matrix'
        sc = 0.0d+00

!$OMP parallel do schedule(dynamic,1) private(iv,ix,j,it,iu,a,b)
        Do i = rank + 1, dimn, nprocs
            iv = indsym(1, i)
            ix = indsym(2, i)
            Do j = i, dimn
                it = indsym(1, j)
                iu = indsym(2, j)

                a = 0.0d+0
                b = 0.0d+0
                Call dim2_density(iv, it, ix, iu, a, b)
                sc(i, j) = DCMPLX(a, b)

                if (ix == it) then
                    a = 0.0d+0
                    b = 0.0d+0
                    Call dim1_density(iv, iu, a, b)
                    sc(i, j) = sc(i, j) - DCMPLX(a, b)
                End if

                sc(j, i) = DCONJG(sc(i, j))

            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=sc)
#endif
        if (debug .and. rank == 0) print *, 'F subspace S matrix is obtained normally'
    End subroutine sFmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE bFmat_complex(dimn, sc, indsym, bc) ! Assume C1 molecule, overlap matrix B in space F
!
!  B(vx, tu) = Siguma_w [eps(w){ <0|EvtExuEww|0> - d(xt)<0|EvuEww|0>}] + S(u,t){-eps(u)-eps(t)}
!
!  v > x, t > u
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in) :: dimn, indsym(2, dimn)
        complex*16, intent(in)  :: sc(dimn, dimn)
        complex*16, intent(out) :: bc(dimn, dimn)

        real(8)    :: e, denr, deni
        complex*16 :: den

        integer(kind=int64) :: it, iu, iv, ix, iw
        integer :: jt, ju, jv, jx, jw, i, j

        bc(:, :) = 0.0d+00

        if (debug .and. rank == 0) print *, 'Start F subspace B matrix'

!$OMP parallel do schedule(dynamic,1) private(iv,jv,ix,jx,j,it,jt,iu,ju,e,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs

            iv = indsym(1, i)
            jv = convert_active_to_global_idx(iv)
            ix = indsym(2, i)
            jx = convert_active_to_global_idx(ix)

            Do j = i, dimn

                it = indsym(1, j)
                jt = convert_active_to_global_idx(it)
                iu = indsym(2, j)
                ju = convert_active_to_global_idx(iu)

!  B(vx, tu) = Siguma_w [eps(w){ <0|EvtExuEww|0> - d(xt)<0|EvuEww|0>}] + S(u,t){-eps(u)-eps(t)}

                e = -eps(ju) - eps(jt)

                Do iw = 1, nact
                    jw = convert_active_to_global_idx(iw)

                    Call dim3_density(iv, it, ix, iu, iw, iw, denr, deni)
                    den = DCMPLX(denr, deni)
                    bc(i, j) = bc(i, j) + den*eps(jw)

                    If (ix == it) then

                        Call dim2_density(iv, iu, iw, iw, denr, deni)
                        den = DCMPLX(denr, deni)
                        bc(i, j) = bc(i, j) - den*eps(jw)

                    End if

                End do

                bc(i, j) = bc(i, j) + sc(i, j)*e

                bc(j, i) = DCONJG(bc(i, j))

            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call reduce_wrapper(mat=bc, root_rank=0)
#endif
        if (debug .and. rank == 0) print *, 'F subspace B matrix is obtained normally'
    End subroutine bFmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE vFmat_complex(nab, iab, v)
!
! V(tu, ab) =  SIGUMA_p,q:active <0|EtpEuq|0>(ap|bq) -  SIGUMA_p:active <0|Etp|0>(au|bp)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_file_manager, only: open_unformatted_file, check_iostat
        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)     :: nab, iab(nsec, nsec)

        complex*16, intent(out) :: v(nab, nact, nact)

        real(8)                  :: dr, di
        complex*16              :: cint2, dens

        integer(kind=int64) :: i, j, k, l, it, iu 
        integer :: tab, i0, jt, ju, iostat, unit_int2, isym, syma
        integer :: multb_s_reverse(nsec, nsec)
        integer :: pattern_t(nact**2, nsymrpa), pattern_u(nact**2, nsymrpa), pattern_tu_count(nsymrpa)
        logical :: is_end_of_file

        if (debug .and. rank == 0) print *, 'Start F subspace V matrix'

! Initialization
        v = 0.0d+00
        multb_s_reverse(:, :) = 0
        call create_multb_s_reverse(multb_s_reverse)

! Save t,u patterns for each isym
        pattern_t(:, :) = 0
        pattern_u(:, :) = 0
        pattern_tu_count(:) = 0
        do isym = 1, nsymrpa
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Do iu = 1, it - 1
                    ju = convert_active_to_global_idx(iu)
                    if (nsymrpa /= 1) syma = MULTB_D(irpamo(ju) - (-1)**(mod(irpamo(ju), 2)), irpamo(jt))
                    if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                        pattern_tu_count(isym) = pattern_tu_count(isym) + 1
                        pattern_t(pattern_tu_count(isym), isym) = it
                        pattern_u(pattern_tu_count(isym), isym) = iu
                    End if
                End do
            End do
        end do
! V(ab,t,u) =  SIGUMA_p,q:active <0|EtpEuq|0>(ap|bq) -  SIGUMA_p:active <0|Etp|0>(au|bp)

        call open_unformatted_file(unit=unit_int2, file=fint, status='old', optional_action='read')  !  (32|32) stored  a > b
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2
            call check_iostat(iostat=iostat, file=fint, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            if (i <= k) cycle ! Read the next line if i is less than or equal to k

            tab = iab(i, k)

! V(ab,t,u) =  SIGUMA_p,q:active <0|EtpEuq|0>(ap|bq) -  SIGUMA_p:active <0|Etp|0>(au|bp)
!                                <0|EtjEul|0>(ij|kl)                             (ij|kl)
!
!                             p=j, q=l loop for t and u             u=j, p=l loop for t
!
!$OMP parallel
!$OMP do schedule(dynamic,1) private(it,iu,dr,di,dens)
            Do it = 1, nact
                Call dim1_density(it, l, dr, di)
                dens = DCMPLX(dr, di)
                v(tab, it, j) = v(tab, it, j) - cint2*dens
            End do
!$OMP end do

            isym = multb_s_reverse(i, k)
!$OMP do schedule(dynamic,1) private(it,iu,dr,di,dens)
            do i0 = 1, pattern_tu_count(isym)
                it = pattern_t(i0, isym)
                iu = pattern_u(i0, isym)

                Call dim2_density(it, j, iu, l, dr, di)
                dens = DCMPLX(dr, di)
                v(tab, it, iu) = v(tab, it, iu) + cint2*dens
            end do
!$OMP end do
!$OMP end parallel

        end do
        close (unit_int2)
        if (debug .and. rank == 0) print *, 'reading Fint2 is over'

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
#endif
        if (debug .and. rank == 0) print *, 'F subspace V matrix is obtained normally'

    end subroutine vFmat_complex
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_F_subspace_real()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx, convert_secondary_to_global_idx

        Implicit NONE

        integer :: dimn, dimm, dammy
        integer, allocatable :: indsym(:, :)
        real(8), allocatable  :: wsnew(:), ws(:), wb(:)
        real(8)               :: e2(2*nsymrpa), alpha, e
        real(8), allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
        real(8), allocatable  :: bc(:, :)
        real(8), allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :), vc(:), vc0(:), vc1(:)
        integer :: j, i, syma, isym, i0
        integer :: ia, it, ib, iu, ja, jt, jb, ju
        integer, allocatable     :: ia0(:), ib0(:), iab(:, :)
        integer                  :: nab

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE F IS NOW CALCULATED
!
!     EatEbu|0>
!
!   DRAS1 =0   DRAS2 = -2   DRAS3 = +2
!
!   a > b t > u ( c > d, v > x)
!
!  S(cvdx,atbu) = d(ac) d(bd)  [ <0|EvtExu|0> - d(xt)<0|Evu|0>]  <= S(vx,tu)
!                               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  S(vx, tu) = <0|EvtExu|0> - d(xt)<0|Evu|0>
!
!  B(vx, tu) = Siguma_w [eps(w){ <0|EvtExuEww|0> - d(xt)<0|EvuEww|0>}] + S(u,t){-eps(u)-eps(t)}
!
!  alpha(a, b) = + eps(a) + eps(b) - e0
!
!  where
!
!  e0 = Siguma_w [eps(w)<0|Eww|0>] (<== calculated as e0 in calce0.f)
!
!  V(tu,ab)   = SIGUMA_p,q:active <0|EtpEuq|0>(ap|bq) -  SIGUMA_p:active <0|Etp|0>(au|bp)
!
!  E2 = SIGUMA_iab,t:dimm |V1(t,ab)|^2|/{(alpha(ab) + wb(t)}
!
        e2 = 0.0d+00
        dimn = 0
        syma = 0
        if (debug .and. rank == 0) print *, 'ENTER solve F part'
        if (rank == 0) print '(10A)', '  '
        if (rank == 0) print '(10A)', ' e2f(isym)'

        i0 = 0
        Do ia = 1, nsec
            Do ib = 1, ia - 1
                i0 = i0 + 1
            End do
        End do

        nab = i0
        Allocate (iab(nsec, nsec))
        iab = 0
        Allocate (ia0(nab))
        Allocate (ib0(nab))

        i0 = 0
        Do ia = 1, nsec
            ja = convert_secondary_to_global_idx(ia)
            Do ib = 1, ia - 1
                jb = convert_secondary_to_global_idx(ib)
                i0 = i0 + 1
                iab(ia, ib) = i0
                iab(ib, ia) = i0
                ia0(i0) = convert_secondary_to_global_idx(ia) ! secondary
                ib0(i0) = convert_secondary_to_global_idx(ib) ! secondary
            End do
        End do

        Allocate (v(nab, nact, nact))
        v = 0.0d+00
        Call vFmat_real(nab, iab, v)
        Do isym = 1, nsymrpa

            dimn = 0
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Do iu = 1, it - 1
                    ju = convert_active_to_global_idx(iu)

!     EatEbu|0>

                    if (nsymrpa /= 1) syma = MULTB_D(irpamo(ju) - (-1)**(mod(irpamo(ju), 2)), irpamo(jt))

                    if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                        dimn = dimn + 1
                    End if
                End do
            End do

            if (debug .and. rank == 0) print *, 'isym, dimn', isym, dimn
            If (dimn == 0) cycle ! Go to the next isym if dimn (dimention of matrix) is zero

            Allocate (indsym(2, dimn))

            dimn = 0
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Do iu = 1, it - 1
                    ju = convert_active_to_global_idx(iu)

!     EatEbu|0>

                    if (nsymrpa /= 1) syma = MULTB_D(irpamo(ju) - (-1)**(mod(irpamo(ju), 2)), irpamo(jt))

                    if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                        dimn = dimn + 1
                        indsym(1, dimn) = it
                        indsym(2, dimn) = iu
                    End if
                End do
            End do

            Allocate (sc(dimn, dimn))
            sc = 0.0d+00            ! sc N*N
            Call sFmat_real(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            Allocate (ws(dimn))

            Allocate (sc0(dimn, dimn))
            sc0 = sc
            Call rdiag(sc, dimn, dimm, ws, smat_lin_dep_threshold)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (debug .and. rank == 0) print *, 'after F subspace S matrix rdiag, new dimension is', dimm
            If (dimm == 0) then
                deallocate (indsym)
                deallocate (sc0)
                deallocate (sc)
                deallocate (ws)
                cycle ! Go to the next isym if dimm (dimention of matrix) is zero
            End if

            Allocate (bc(dimn, dimn))                                 ! bc N*N
            bc = 0.0d+00
            Call bFmat_real(dimn, sc0, indsym, bc)
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

            Do i0 = 1, nab
                ja = ia0(i0)
                jb = ib0(i0)

!     EatEbu|0>
                if (nsymrpa /= 1) then
                    syma = MULTB_D(irpamo(ja), irpamo(jb) - (-1)**(mod(irpamo(jb), 2)))
                    syma = MULTB_S(syma, isym)
                end if
                If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then

                    Allocate (vc(dimn))
                    Do it = 1, dimn
                        vc(it) = v(i0, indsym(1, it), indsym(2, it))
                    End do
                    allocate (vc0(dimm)); Call memplus(KIND(vc0), SIZE(vc0), 2)
                    allocate (vc1(dimm)); Call memplus(KIND(vc1), SIZE(vc1), 2)
                    call gemv(transpose(uc), vc, vc0)
                    Deallocate (vc)

                    alpha = +eps(ja) + eps(jb) - e0 + eshift  ! For Level Shift (2007/2/9)

                    call gemv(transpose(bc1), vc0, vc1)

                    Do j = 1, dimm
                        e = (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                        sumc2_subspace(subspace_idx) = sumc2_subspace(subspace_idx) + e/(alpha + wb(j))
                        e2(isym) = e2(isym) - e
                    End do

                    Call memminus(KIND(vc0), SIZE(vc0), 2); Deallocate (vc0)
                    Call memminus(KIND(vc1), SIZE(vc1), 2); Deallocate (vc1)

                End if

            End do

            deallocate (indsym)
            deallocate (uc)
            deallocate (wb)
            Deallocate (bc1)

            if (rank == 0) print '(" e2f(",I3,") = ",E25.15," a.u.")', isym, e2(isym)
            e2_subspace(subspace_idx) = e2_subspace(subspace_idx) + e2(isym)
        End do

        if (rank == 0) then
            print '(" e2f      = ",E25.15," a.u.")', e2_subspace(subspace_idx)
            print '(" sumc2,f  = ",E25.15)', sumc2_subspace(subspace_idx)
        end if

        deallocate (iab)
        deallocate (ia0)
        deallocate (ib0)
        deallocate (v)

        if (debug .and. rank == 0) print *, 'end solve_F_subspace'
    end subroutine solve_F_subspace_real
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE sFmat_real(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space F

!  S(vx, tu) = <0|EvtExu|0> - d(xt)<0|Evu|0>
!
!     v > x, t > u

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)   :: dimn, indsym(2, dimn)
        real(8), intent(out)  :: sc(dimn, dimn)
        real(8) :: a, b
        integer(kind=int64) :: it, iu, iv, ix
        integer :: i, j

        if (debug .and. rank == 0) print *, 'Start F subspace S matrix'
        sc = 0.0d+00

!$OMP parallel do schedule(dynamic,1) private(iv,ix,j,it,iu,a,b)
        Do i = rank + 1, dimn, nprocs
            iv = indsym(1, i)
            ix = indsym(2, i)
            Do j = i, dimn
                it = indsym(1, j)
                iu = indsym(2, j)

                a = 0.0d+0
                b = 0.0d+0
                Call dim2_density(iv, it, ix, iu, a, b)
                sc(i, j) = a

                if (ix == it) then
                    a = 0.0d+0
                    b = 0.0d+0
                    Call dim1_density(iv, iu, a, b)
                    sc(i, j) = sc(i, j) - a
                End if

                sc(j, i) = sc(i, j)

            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=sc)
#endif
        if (debug .and. rank == 0) print *, 'F subspace S matrix is obtained normally'
    End subroutine sFmat_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE bFmat_real(dimn, sc, indsym, bc) ! Assume C1 molecule, overlap matrix B in space F
!
!  B(vx, tu) = Siguma_w [eps(w){ <0|EvtExuEww|0> - d(xt)<0|EvuEww|0>}] + S(u,t){-eps(u)-eps(t)}
!
!  v > x, t > u
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in) :: dimn, indsym(2, dimn)
        real(8), intent(in)  :: sc(dimn, dimn)
        real(8), intent(out) :: bc(dimn, dimn)

        real(8) :: e, denr, deni
        real(8) :: den

        integer(kind=int64) :: it, iu, iv, ix, iw
        integer :: jt, ju, jv, jx, jw, i, j

        bc(:, :) = 0.0d+00

        if (debug .and. rank == 0) print *, 'Start F subspace B matrix'

!$OMP parallel do schedule(dynamic,1) private(iv,jv,ix,jx,j,it,jt,iu,ju,e,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs

            iv = indsym(1, i)
            jv = convert_active_to_global_idx(iv)
            ix = indsym(2, i)
            jx = convert_active_to_global_idx(ix)

            Do j = i, dimn

                it = indsym(1, j)
                jt = convert_active_to_global_idx(it)
                iu = indsym(2, j)
                ju = convert_active_to_global_idx(iu)

!  B(vx, tu) = Siguma_w [eps(w){ <0|EvtExuEww|0> - d(xt)<0|EvuEww|0>}] + S(u,t){-eps(u)-eps(t)}

                e = -eps(ju) - eps(jt)

                Do iw = 1, nact
                    jw = convert_active_to_global_idx(iw)

                    Call dim3_density(iv, it, ix, iu, iw, iw, denr, deni)
                    den = denr
                    bc(i, j) = bc(i, j) + den*eps(jw)

                    If (ix == it) then

                        Call dim2_density(iv, iu, iw, iw, denr, deni)
                        den = denr
                        bc(i, j) = bc(i, j) - den*eps(jw)

                    End if

                End do

                bc(i, j) = bc(i, j) + sc(i, j)*e

                bc(j, i) = bc(i, j)

            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call reduce_wrapper(mat=bc, root_rank=0)
#endif
        if (debug .and. rank == 0) print *, 'F subspace B matrix is obtained normally'
    End subroutine bFmat_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE vFmat_real(nab, iab, v)
!
! V(tu, ab) =  SIGUMA_p,q:active <0|EtpEuq|0>(ap|bq) -  SIGUMA_p:active <0|Etp|0>(au|bp)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_file_manager, only: open_unformatted_file, check_iostat
        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)     :: nab, iab(nsec, nsec)

        real(8), intent(out) :: v(nab, nact, nact)

        real(8) :: dr, di
        real(8) :: cint2, dens

        integer(kind=int64) :: i, j, k, l, it, iu
        integer ::  tab, i0, jt, ju, iostat, unit_int2, isym, syma
        integer :: multb_s_reverse(nsec, nsec)
        integer :: pattern_t(nact**2, nsymrpa), pattern_u(nact**2, nsymrpa), pattern_tu_count(nsymrpa)
        logical :: is_end_of_file

        if (debug .and. rank == 0) print *, 'Start F subspace V matrix'

! Initialization
        v = 0.0d+00
        multb_s_reverse(:, :) = 0
        call create_multb_s_reverse(multb_s_reverse)

! Save t,u patterns for each isym
        pattern_t(:, :) = 0
        pattern_u(:, :) = 0
        pattern_tu_count(:) = 0
        do isym = 1, nsymrpa
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Do iu = 1, it - 1
                    ju = convert_active_to_global_idx(iu)
                    if (nsymrpa /= 1) syma = MULTB_D(irpamo(ju) - (-1)**(mod(irpamo(ju), 2)), irpamo(jt))
                    if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                        pattern_tu_count(isym) = pattern_tu_count(isym) + 1
                        pattern_t(pattern_tu_count(isym), isym) = it
                        pattern_u(pattern_tu_count(isym), isym) = iu
                    End if
                End do
            End do
        end do
! V(ab,t,u) =  SIGUMA_p,q:active <0|EtpEuq|0>(ap|bq) -  SIGUMA_p:active <0|Etp|0>(au|bp)

        call open_unformatted_file(unit=unit_int2, file=fint, status='old', optional_action='read')  !  (32|32) stored  a > b
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2
            call check_iostat(iostat=iostat, file=fint, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
            if (i <= k) cycle ! Read the next line if i is less than or equal to k

            tab = iab(i, k)

! V(ab,t,u) =  SIGUMA_p,q:active <0|EtpEuq|0>(ap|bq) -  SIGUMA_p:active <0|Etp|0>(au|bp)
!                                <0|EtjEul|0>(ij|kl)                             (ij|kl)
!
!                             p=j, q=l loop for t and u             u=j, p=l loop for t
!
!$OMP parallel
!$OMP do schedule(dynamic,1) private(it,iu,dr,di,dens)
            Do it = 1, nact
                Call dim1_density(it, l, dr, di)
                dens = dr
                v(tab, it, j) = v(tab, it, j) - cint2*dens
            End do
!$OMP end do

            isym = multb_s_reverse(i, k)
!$OMP do schedule(dynamic,1) private(it,iu,dr,di,dens)
            do i0 = 1, pattern_tu_count(isym)
                it = pattern_t(i0, isym)
                iu = pattern_u(i0, isym)

                Call dim2_density(it, j, iu, l, dr, di)
                dens = dr
                v(tab, it, iu) = v(tab, it, iu) + cint2*dens
            end do
!$OMP end do
!$OMP end parallel

        end do
        close (unit_int2)
        if (debug .and. rank == 0) print *, 'reading Fint2 is over'

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
#endif
        if (debug .and. rank == 0) print *, 'F subspace V matrix is obtained normally'

    end subroutine vFmat_real
    subroutine create_multb_s_reverse(multb_s_reverse)
!========================================================================================================
! This subroutine creates multb_s_reverse
!
! multb_s_reverse(i, j) returns the symmetry of MULTB_D(irpamo(ju) - (-1)**(mod(irpamo(ju), 2)), irpamo(jt))
!========================================================================================================
        use module_index_utils, only: convert_secondary_to_global_idx
        implicit none
        integer, intent(inout) :: multb_s_reverse(:, :)
        integer :: ia, ib, ja, jb
        integer :: isym, syma

        if (nsymrpa == 1) then
            multb_s_reverse(:, :) = 1
        else
            do ia = 1, nsec
                ja = convert_secondary_to_global_idx(ia)
                do ib = 1, ia - 1
                    jb = convert_secondary_to_global_idx(ib)
                    syma = MULTB_D(irpamo(ja), irpamo(jb) - (-1)**(mod(irpamo(jb), 2)))
                    do isym = 1, nsymrpa
                        if (MULTB_S(syma, isym) == 1) then
                            multb_s_reverse(ia, ib) = isym
                            exit
                        end if
                    end do
                end do
            end do
        end if
    end subroutine create_multb_s_reverse
end subroutine solve_F_subspace
