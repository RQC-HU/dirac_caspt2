SUBROUTINE solve_A_subspace(e0, e2a)

    use module_ulambda_s_half, only: ulambda_s_half
    use module_global_variables
    use module_realonly, only: realonly
    implicit none
    real(8), intent(in) :: e0
    real(8), intent(out):: e2a

    if (realonly%is_realonly()) then
        call solve_A_subspace_real()
    else
        call solve_A_subspace_complex()
    end if

contains
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_A_subspace_complex()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx

        Implicit NONE

        integer :: dimn, dimm, dammy

        integer, allocatable :: indsym(:, :)

        real(8), allocatable  :: wsnew(:), ws(:), wb(:)
        real(8)               :: e2(2*nsymrpa), alpha

        complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
        complex*16, allocatable  :: bc(:, :)
        complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :, :), vc(:), vc1(:)

        integer :: i, j, syma, symb, isym, sym1
        integer :: ix, iy, iz, ii, dimi, ixyz
        integer :: jx, jy, jz, it

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE A IS NOW CALCULATED
!
!     EtiEuv|0>               t > u
!
!   DRAS1 = -1   DRAS2 = +1   DRAS3 = 0
!
!! TABUN USO  x > y,  t > u,  y /= z , u /= v
!
!  S(xjyz,tiuv) = d(j,i)[ - <0|EzyEtxEuv|0> + d(tx)<0|EzyEuv|0>]
!
!  S(xyz,tuv) =  - <0|EzyEtxEuv|0> + d(tx)<0|EzyEuv|0>
!
!  B(xyz,tuv) = Siguma_w eps(w){-<0|EzyEtxEuvEww|0> + d(tx)<0|EzyEuvEww|0>}
!
!               + S(xyz,tuv)(eps(u)+eps(t)-eps(v))
!
!  alpha(i)       = - eps(i) - Siguma_w [eps(w)<0|Eww|0>] (<== calculated as e0 in calce0.f)
!
!  V(tuv,i)=  - SIGUMA_p,q,r:act <0|EvuEptEqr|0>(pi|qr)
!
!  + SIGUMA_p,q:act  <0|EvuEpq|0> (pq|ti)
!
!  - SIGUMA_p:act <0|EvuEpt|0>[h(pi)+ SIGUMA_k:inact{(pi|kk)-(pk|ki)}]
!
!  +  <0|Evu|0>[h(ti) + SIGUMA_k:inact{(ti|kk) - (tk|ki)}]
!
!  E2 = SIGUMA_i, dimm |Vc1(dimm,i)|^2|/{(alpha(i) + wb(dimm)}

        e2 = 0.0d+00
        e2a = 0.0d+00
        dimi = 0
        dimn = 0
        syma = 0
        if (rank == 0) print *, 'ENTER solve A part'
        Allocate (v(ninact, nact, nact, nact))
        Call memplus(KIND(v), SIZE(v), 2)
        Call vAmat_complex(v)
!         ExjEyz
        Do isym = 1, nsymrpa

            ixyz = 0
            Do ix = 1, nact
                Do iy = 1, nact
                    Do iz = 1, nact
                        jx = convert_active_to_global_idx(ix)
                        jy = convert_active_to_global_idx(iy)
                        jz = convert_active_to_global_idx(iz)
                        if (nsymrpa /= 1) then
                            syma = MULTB_D(irpamo(jx), isym)
                            symb = MULTB_D(irpamo(jy), irpamo(jz))
                            syma = MULTB_S(syma, symb)
                        end if
                        If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then
                            ixyz = ixyz + 1
                        End if
                    End do
                End do
            End do

            dimn = ixyz

            If (dimn == 0) cycle ! Go to the next isym.

            Allocate (indsym(3, dimn)); Call memplus(KIND(indsym), SIZE(indsym), 1)

            ixyz = 0

            Do ix = 1, nact
                Do iy = 1, nact
                    Do iz = 1, nact
                        jx = convert_active_to_global_idx(ix)
                        jy = convert_active_to_global_idx(iy)
                        jz = convert_active_to_global_idx(iz)
                        if (nsymrpa /= 1) then
                            syma = MULTB_D(irpamo(jx), isym)
                            symb = MULTB_D(irpamo(jy), irpamo(jz))
                            syma = MULTB_S(syma, symb)
                        end if

                        If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then
                            ixyz = ixyz + 1
                            indsym(1, ixyz) = ix
                            indsym(2, ixyz) = iy
                            indsym(3, ixyz) = iz
                        End if
                    End do
                End do
            End do

            if (rank == 0) print *, 'isym, dimn', isym, dimn
            Allocate (sc(dimn, dimn)); Call memplus(KIND(sc), SIZE(sc), 2)

            sc = 0.0d+00            ! sr N*N
            Call sAmat_complex(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Allocate (ws(dimn)); Call memplus(KIND(ws), SIZE(ws), 1)

            Allocate (sc0(dimn, dimn)); Call memplus(KIND(sc0), SIZE(sc0), 2)
            sc0 = sc
            Call cdiag(sc, dimn, dimm, ws, smat_lin_dep_threshold)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'after A subspace S matrix cdiag, new dimension is', dimm

            If (dimm == 0) then
                Call memminus(KIND(indsym), SIZE(indsym), 1); deallocate (indsym)
                Call memminus(KIND(sc0), SIZE(sc0), 2); deallocate (sc0)
                Call memminus(KIND(sc), SIZE(sc), 2); deallocate (sc)
                Call memminus(KIND(ws), SIZE(ws), 1); deallocate (ws)
                cycle ! Go to the next isym.
            End if

            If (debug) then
                if (rank == 0) print *, 'Check whether U*SU is diagonal'
                Call checkdgc(dimn, sc0, sc, ws)
                if (rank == 0) print *, 'Check whether U*SU is diagonal END'
            End if

            Allocate (bc(dimn, dimn)); Call memplus(KIND(bc), SIZE(bc), 2)   ! br N*N
            bc = 0.0d+00
            Call bAmat_complex(dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Call memminus(KIND(sc0), SIZE(sc0), 2); deallocate (sc0)

            Allocate (uc(dimn, dimm)); Call memplus(KIND(uc), SIZE(uc), 2)           ! uc N*M
            Allocate (wsnew(dimm)); Call memplus(KIND(wsnew), SIZE(wsnew), 1)     ! wnew M
            uc(:, :) = 0.0d+00
            wsnew(:) = 0.0d+00

            Call ccutoff(sc, ws, dimn, dimm, smat_lin_dep_threshold, uc, wsnew)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Call memminus(KIND(sc), SIZE(sc), 2); deallocate (sc)
            Call memminus(KIND(ws), SIZE(ws), 1); deallocate (ws)

            Call ulambda_s_half(uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Call memminus(KIND(wsnew), SIZE(wsnew), 1); deallocate (wsnew)

            Allocate (bc0(dimm, dimn)); Call memplus(KIND(bc0), SIZE(bc0), 2) ! bc0 M*N
            bc0 = 0.0d+00
            bc0 = MATMUL(TRANSPOSE(DCONJG(uc)), bc)

            Allocate (bc1(dimm, dimm)); Call memplus(KIND(bc1), SIZE(bc1), 2) ! bc1 M*M
            bc1 = 0.0d+00
            bc1 = MATMUL(bc0, uc)

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

            Call memminus(KIND(bc), SIZE(bc), 2); deallocate (bc)
            Call memminus(KIND(bc0), SIZE(bc0), 2); deallocate (bc0)

            Allocate (wb(dimm)); Call memplus(KIND(wb), SIZE(wb), 1)

            if (rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'

            Allocate (bc0(dimm, dimm)); Call memplus(KIND(bc0), SIZE(bc0), 2) ! bc0 M*M
            bc0 = bc1

            Call cdiag(bc1, dimm, dammy, wb, bmat_no_cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            If (debug) then
                if (rank == 0) print *, 'Check whether bc is really diagonalized or not'
                Call checkdgc(dimm, bc0, bc1, wb)
                if (rank == 0) print *, 'Check whether bc is really diagonalized or not END'
            End if
            Call memminus(KIND(bc0), SIZE(bc0), 2); deallocate (bc0)

            if (rank == 0) print *, 'bC1 matrix is diagonalized!'

            e2 = 0.0d+00

            Do ii = 1, ninact
                sym1 = irpamo(ii)
                if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. sym1 == isym)) then
                    Allocate (vc(dimn)); Call memplus(KIND(vc), SIZE(vc), 2)
                    Do it = 1, dimn
                        vc(it) = v(ii, indsym(1, it), indsym(2, it), indsym(3, it))
                    End do

                    Allocate (vc1(dimm)); Call memplus(KIND(vc1), SIZE(vc1), 2)
                    vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn, 1:dimm))), vc(1:dimn))
                    Call memminus(KIND(vc), SIZE(vc), 2); Deallocate (vc)

                    alpha = -eps(ii) - e0 + eshift   ! For Level Shift (2007/2/9)

                    vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm, 1:dimm))), vc1(1:dimm))

                    Do j = 1, dimm
                        sumc2local = sumc2local + (ABS(vc1(j))**2.0d+00)/((alpha + wb(j))**2.0d+00)
                        e2(isym) = e2(isym) - (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                    End do
                    Call memminus(KIND(vc1), SIZE(vc1), 2); Deallocate (vc1)
                End if
            End do
            if (rank == 0) print '("e2a(",I3,") = ",E20.10," a.u.")', isym, e2(isym)

            Call memminus(KIND(bc1), SIZE(bc1), 2); Deallocate (bc1)
            Call memminus(KIND(uc), SIZE(uc), 2); Deallocate (uc)
            Call memminus(KIND(wb), SIZE(wb), 1); Deallocate (wb)
            Call memminus(KIND(indsym), SIZE(indsym), 2); Deallocate (indsym)

            e2a = e2a + e2(isym)
        End do

        if (rank == 0) then
            print '("e2a      = ",E20.10," a.u.")', e2a

            print '("sumc2,a  = ",E20.10)', sumc2local
        end if
        sumc2 = sumc2 + sumc2local

        Call memminus(KIND(v), SIZE(v), 2); Deallocate (v)

        if (rank == 0) print *, 'end solve_A_subspace'
    end subroutine solve_A_subspace_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
    SUBROUTINE sAmat_complex(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space A
!
!  S(xyz,tuv) =  - <0|EzyEtxEuv|0> + d(tx)<0|EzyEuv|0>
!
!     x > y, t > u

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)      :: dimn, indsym(3, dimn)
        complex*16, intent(out)  :: sc(dimn, dimn)
        real(8)  ::a, b
        integer :: it, iu, iv, ix, iy, iz
        integer :: i, j

        if (rank == 0) print *, 'Start A subspace S matrix'
! Initialization
        sc = 0.0d+00

!$OMP parallel do private(i,ix,iy,iz,j,it,iu,iv,a,b)
        Do i = rank + 1, dimn, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
            ix = indsym(1, i)
            iy = indsym(2, i)
            iz = indsym(3, i)

            Do j = i, dimn
                it = indsym(1, j)
                iu = indsym(2, j)
                iv = indsym(3, j)

!  S(xyz,tuv) =  - <0|EzyEtxEuv|0> + d(tx)<0|EzyEuv|0>
                Call dim3_density(iz, iy, it, ix, iu, iv, a, b)
                sc(i, j) = sc(i, j) - DCMPLX(a, b)

                If (it == ix) then
                    Call dim2_density(iz, iy, iu, iv, a, b)
                    sc(i, j) = sc(i, j) + DCMPLX(a, b)
                End if
                sc(j, i) = DCONJG(sc(i, j))
            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=sc)
#endif
        if (rank == 0) print *, 'A subspace S matrix is obtained normally'
    End subroutine sAmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE bAmat_complex(dimn, sc, indsym, bc) ! Assume C1 molecule, overlap matrix B in space A
!
!  B(xyz,tuv) = Siguma_w eps(w){-<0|EzyEtxEuvEww|0> + d(tx)<0|EzyEuvEww|0>}
!
!               + S(xyz,tuv)(eps(u)+eps(t)-eps(v))
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer :: it, iu, iv, ix, iy, iz, iw
        integer :: jt, ju, jv, jx, jy, jz, jw
        integer :: i, j

        integer, intent(in)     :: dimn, indsym(3, dimn)
        complex*16, intent(in)  :: sc(dimn, dimn)
        complex*16, intent(out) :: bc(dimn, dimn)

        real(8)               :: e, denr, deni
        complex*16           :: den

        if (rank == 0) print *, 'Start A subspace B matrix'
        bc(:, :) = 0.0d+00
!$OMP parallel do private(ix,iy,iz,jx,jy,jz,it,iu,iv,jt,ju,jv,e,j,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
            ix = indsym(1, i)
            iy = indsym(2, i)
            iz = indsym(3, i)
            jx = convert_active_to_global_idx(ix)
            jy = convert_active_to_global_idx(iy)
            jz = convert_active_to_global_idx(iz)

            Do j = i, dimn

                it = indsym(1, j)
                iu = indsym(2, j)
                iv = indsym(3, j)
                jt = convert_active_to_global_idx(it)
                ju = convert_active_to_global_idx(iu)
                jv = convert_active_to_global_idx(iv)

!  B(xyz,tuv) = Siguma_w eps(w){-<0|EzyEtxEuvEww|0> + d(tx)<0|EzyEuvEww|0>}
!
!               + S(xyz,tuv)(eps(u)+eps(t)-eps(v))

                e = eps(ju) + eps(jt) - eps(jv)

                Do iw = 1, nact
                    jw = convert_active_to_global_idx(iw)

                    Call dim4_density(iz, iy, it, ix, iu, iv, iw, iw, denr, deni)
                    den = DCMPLX(denr, deni)
                    bc(i, j) = bc(i, j) - den*eps(jw)

                    If (it == ix) then
                        Call dim3_density(iz, iy, iu, iv, iw, iw, denr, deni)
                        den = DCMPLX(denr, deni)
                        bc(i, j) = bc(i, j) + den*eps(jw)

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

        if (rank == 0) print *, 'A subspace B matrix is obtained normally'
    End subroutine bAmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
    SUBROUTINE vAmat_complex(v)
!
! Assume C1 molecule, V=<0|H|i> matrix in space A
!
!  V(tuv,i)=  - SIGUMA_p,q,r:act <0|EvuEptEqr|0>(pi|qr)
!
!  + SIGUMA_p,q:act  <0|EvuEpq|0> (pq|ti)
!
!  - SIGUMA_p:act <0|EvuEpt|0>[h(pi)+ SIGUMA_k:inact{(pi|kk)-(pk|ki)}]
!
!  +  <0|Evu|0>[h(ti) + SIGUMA_k:inact{(ti|kk) - (tk|ki)}]
!
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

        complex*16, intent(out) :: v(ninact, nact, nact, nact)
        real(8)                  :: dr, di
        complex*16              :: cint2, d, dens1(nact, nact), effh(nact, ninact)
        complex*16              :: cint1
        logical                 :: is_end_of_file

        integer :: it, iu, iv, ii, ip
        integer :: jt, ju, jv, ji, jp
        integer :: i, j, k, l, dim(nsymrpa)
        integer :: dim2(nsymrpa), isym, i0, syma, symb, symc, iostat, unit_int2
        integer, allocatable :: indt(:, :), indu(:, :), indv(:, :)
        integer, allocatable :: ind2u(:, :), ind2v(:, :)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  V(tuv,i)=  - SIGUMA_p,q,r:act <0|EvuEptEqr|0>(pi|qr)
!
!  + SIGUMA_p,q:act  <0|EvuEpq|0> (pq|ti)
!
!  - SIGUMA_p:act <0|EvuEpt|0>[h(pi)+ SIGUMA_k:inact{(pi|kk)-(pk|ki)}]
!
!  +  <0|Evu|0>[h(ti) + SIGUMA_k:inact{(ti|kk) - (tk|ki)}]
!
!   = - SIGUMA_p,q,r:act <0|EvuEptEqr|0>(pi|qr)
!
!  + SIGUMA_p,q:act  <0|EvuEpq|0> (pq|ti)
!
!  - SIGUMA_p:act <0|EvuEpt|0>effh(pi)  +  <0|Evu|0>effh(ti)
!                 ========================================= This part is calculated after reading int2
!
!  effh is stored in memory while reading int2.
!
!  effh(p,i) = h(pi)+ SIGUMA_k:inact{(pi|kk)-(pk|ki)}
        if (rank == 0) print *, 'Start A subspace V matrix'
        v = 0.0d+00
        dens1 = 0.0d+00
        effh = 0.0d+00
        dim = 0

        Allocate (indt(nact**3, nsymrpa)); Call memplus(KIND(indt), SIZE(indt), 1)
        Allocate (indu(nact**3, nsymrpa)); Call memplus(KIND(indu), SIZE(indu), 1)
        Allocate (indv(nact**3, nsymrpa)); Call memplus(KIND(indv), SIZE(indv), 1)
        indt = 0
        indu = 0
        indv = 0
        dim = 0
!$OMP parallel do schedule(static) private(it,jt,iv,jv,iu,ju,syma,symb,symc)
        Do isym = 1, nsymrpa
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Do iv = 1, nact
                    jv = convert_active_to_global_idx(iv)
                    Do iu = 1, nact
                        ju = convert_active_to_global_idx(iu)
! EtiEuv
                        if (nsymrpa /= 1) then
                            syma = MULTB_D(irpamo(jt), isym)
                            symb = MULTB_D(irpamo(ju), irpamo(jv))
                            symc = MULTB_S(syma, symb)
                        end if
                        if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. symc == 1)) then
                            dim(isym) = dim(isym) + 1
                            indt(dim(isym), isym) = it
                            indu(dim(isym), isym) = iu
                            indv(dim(isym), isym) = iv
                        end if
                    End do
                End do
            End do
        End do
!$OMP end parallel do
        Allocate (ind2u(nact**2, nsymrpa)); Call memplus(KIND(ind2u), SIZE(ind2u), 1)
        Allocate (ind2v(nact**2, nsymrpa)); Call memplus(KIND(ind2v), SIZE(ind2v), 1)
        ind2u = 0.0d+00
        ind2v = 0.0d+00
        dim2 = 0
!$OMP parallel do schedule(static) private(iu,ju,iv,jv,syma)
        Do isym = 1, nsymrpa
            Do iu = 1, nact
                ju = convert_active_to_global_idx(iu)
                Do iv = 1, nact
                    jv = convert_active_to_global_idx(iv)
                    if (nsymrpa /= 1) syma = MULTB_D(irpamo(jv), irpamo(ju))
                    if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                        dim2(isym) = dim2(isym) + 1
                        ind2u(dim2(isym), isym) = iu
                        ind2v(dim2(isym), isym) = iv
                    end if
                End do
            End do
        End do
!$OMP end parallel do

!$OMP parallel do private(ji,it,jt,cint1)
        Do ii = rank + 1, ninact, nprocs
            ji = ii
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Call tramo1(jt, ji, cint1)
                effh(it, ii) = cint1
            End do
        End do
!$OMP end parallel do

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Two types of integrals are stored
!
!  (21|22) stored (pi|qr) ...TYPE 1
!  (21|11) stored (pi|jk) ...TYPE 2
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2, file=a1int, status='old', optional_action='read')
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
            call check_iostat(iostat=iostat, file=a1int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  V(tuv,i)=  - SIGUMA_p,q,r:act <0|EvuEptEqr|0>(pi|qr)
!
!  + SIGUMA_p,q:act  <0|EvuEpq|0> (ti|pq)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            isym = irpamo(j)
!$OMP parallel do private(it,iu,iv,jt,ju,jv,dr,di,d)
            Do i0 = 1, dim(isym)
                it = indt(i0, isym)
                iu = indu(i0, isym)
                iv = indv(i0, isym)
                jt = convert_active_to_global_idx(it)
                ju = convert_active_to_global_idx(iu)
                jv = convert_active_to_global_idx(iv)

                Call dim3_density(iv, iu, i, it, k, l, dr, di)
                d = DCMPLX(dr, di)
                v(j, it, iu, iv) = v(j, it, iu, iv) - cint2*d
            End do
!$OMP end parallel do

            isym = MULTB_D(irpamo(convert_active_to_global_idx(i)), irpamo(j))           ! j coresponds to ii, i coresponds to it

!$OMP parallel do private(iu,iv,ju,jv,dr,di,d)
            Do i0 = 1, dim2(isym)
                iu = ind2u(i0, isym)
                iv = ind2v(i0, isym)
                ju = convert_active_to_global_idx(iu)
                jv = convert_active_to_global_idx(iv)

                Call dim2_density(iv, iu, k, l, dr, di)
                d = DCMPLX(dr, di)
                v(j, i, iu, iv) = v(j, i, iu, iv) + cint2*d
            End do
!$OMP end parallel do

        end do

        close (unit_int2)
        if (rank == 0) print *, 'reading A1int2 is over'

        call open_unformatted_file(unit=unit_int2, file=a2int, status='old', optional_action='read') ! TYPE 2 integrals
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
            call check_iostat(iostat=iostat, file=a2int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  effh(p,i) = h(pi)+ SIGUMA_k:inact{(pi|kk)-(pk|ki)}
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (k == l .and. j /= k) then       ! (PI|KK) type
                effh(i, j) = effh(i, j) + cint2
            elseif (j == k .and. k /= l) then       ! (PK|KI) type
                effh(i, l) = effh(i, l) - cint2
            end if
        end do

        close (unit_int2)
        if (rank == 0) print *, 'reading A2int2 is over'

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=effh)
#endif

!  - SIGUMA_p:act <0|EvuEpt|0>effh(pi)  +  <0|Evu|0>effh(ti)

!$OMP parallel do private(ji,isym,it,iu,iv,dr,di,d,ip,jp)
        Do ii = rank + 1, ninact, nprocs
            ji = ii
            isym = irpamo(ji)

            Do i0 = 1, dim(isym)
                it = indt(i0, isym)
                iu = indu(i0, isym)
                iv = indv(i0, isym)

                Call dim1_density(iv, iu, dr, di)
                d = DCMPLX(dr, di)
                v(ii, it, iu, iv) = v(ii, it, iu, iv) + effh(it, ii)*d

                Do ip = 1, nact
                    jp = convert_active_to_global_idx(ip)

                    Call dim2_density(iv, iu, ip, it, dr, di)
                    d = DCMPLX(dr, di)
                    v(ii, it, iu, iv) = v(ii, it, iu, iv) - effh(ip, ii)*d
                End do
            End do
        End do
!$OMP end parallel do

        ! deallocate memory
        Call memminus(KIND(indt), SIZE(indt), 1); deallocate (indt)
        Call memminus(KIND(indu), SIZE(indu), 1); deallocate (indu)
        Call memminus(KIND(indv), SIZE(indv), 1); deallocate (indv)
        Call memminus(KIND(ind2u), SIZE(ind2u), 1); deallocate (ind2u)
        Call memminus(KIND(ind2v), SIZE(ind2v), 1); deallocate (ind2v)

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
#endif
        if (rank == 0) print *, 'A subspace V matrix is obtained normally'
    end subroutine vAmat_complex
    ! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_A_subspace_real()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx

        Implicit NONE

        integer :: dimn, dimm, dammy

        integer, allocatable :: indsym(:, :)

        real(8), allocatable  :: wsnew(:), ws(:), wb(:)
        real(8)               :: e2(2*nsymrpa), alpha

        real(8), allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
        real(8), allocatable  :: bc(:, :)
        real(8), allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :, :), vc(:), vc1(:)

        integer :: i, j, syma, symb, isym, sym1
        integer :: ix, iy, iz, ii, dimi, ixyz
        integer :: jx, jy, jz, it

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE A IS NOW CALCULATED
!
!     EtiEuv|0>               t > u
!
!   DRAS1 = -1   DRAS2 = +1   DRAS3 = 0
!
!! TABUN USO  x > y,  t > u,  y /= z , u /= v
!
!  S(xjyz,tiuv) = d(j,i)[ - <0|EzyEtxEuv|0> + d(tx)<0|EzyEuv|0>]
!
!  S(xyz,tuv) =  - <0|EzyEtxEuv|0> + d(tx)<0|EzyEuv|0>
!
!  B(xyz,tuv) = Siguma_w eps(w){-<0|EzyEtxEuvEww|0> + d(tx)<0|EzyEuvEww|0>}
!
!               + S(xyz,tuv)(eps(u)+eps(t)-eps(v))
!
!  alpha(i)       = - eps(i) - Siguma_w [eps(w)<0|Eww|0>] (<== calculated as e0 in calce0.f)
!
!  V(tuv,i)=  - SIGUMA_p,q,r:act <0|EvuEptEqr|0>(pi|qr)
!
!  + SIGUMA_p,q:act  <0|EvuEpq|0> (pq|ti)
!
!  - SIGUMA_p:act <0|EvuEpt|0>[h(pi)+ SIGUMA_k:inact{(pi|kk)-(pk|ki)}]
!
!  +  <0|Evu|0>[h(ti) + SIGUMA_k:inact{(ti|kk) - (tk|ki)}]
!
!  E2 = SIGUMA_i, dimm |Vc1(dimm,i)|^2|/{(alpha(i) + wb(dimm)}

        e2 = 0.0d+00
        e2a = 0.0d+00
        dimi = 0
        dimn = 0
        syma = 0
        if (rank == 0) print *, 'ENTER solve A part'
        Allocate (v(ninact, nact, nact, nact))
        Call memplus(KIND(v), SIZE(v), 2)
        Call vAmat_real(v)
!         ExjEyz
        Do isym = 1, nsymrpa

            ixyz = 0
            Do ix = 1, nact
                Do iy = 1, nact
                    Do iz = 1, nact
                        jx = convert_active_to_global_idx(ix)
                        jy = convert_active_to_global_idx(iy)
                        jz = convert_active_to_global_idx(iz)
                        if (nsymrpa /= 1) then
                            syma = MULTB_D(irpamo(jx), isym)
                            symb = MULTB_D(irpamo(jy), irpamo(jz))
                            syma = MULTB_S(syma, symb)
                        end if
                        If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then
                            ixyz = ixyz + 1
                        End if
                    End do
                End do
            End do

            dimn = ixyz

            If (dimn == 0) cycle ! Go to the next isym.

            Allocate (indsym(3, dimn)); Call memplus(KIND(indsym), SIZE(indsym), 1)

            ixyz = 0

            Do ix = 1, nact
                Do iy = 1, nact
                    Do iz = 1, nact
                        jx = convert_active_to_global_idx(ix)
                        jy = convert_active_to_global_idx(iy)
                        jz = convert_active_to_global_idx(iz)
                        if (nsymrpa /= 1) then
                            syma = MULTB_D(irpamo(jx), isym)
                            symb = MULTB_D(irpamo(jy), irpamo(jz))
                            syma = MULTB_S(syma, symb)
                        end if

                        If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then
                            ixyz = ixyz + 1
                            indsym(1, ixyz) = ix
                            indsym(2, ixyz) = iy
                            indsym(3, ixyz) = iz
                        End if
                    End do
                End do
            End do

            if (rank == 0) print *, 'isym, dimn', isym, dimn
            Allocate (sc(dimn, dimn)); Call memplus(KIND(sc), SIZE(sc), 2)

            sc = 0.0d+00            ! sr N*N
            Call sAmat_real(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Allocate (ws(dimn)); Call memplus(KIND(ws), SIZE(ws), 1)

            Allocate (sc0(dimn, dimn)); Call memplus(KIND(sc0), SIZE(sc0), 2)
            sc0 = sc
            Call rdiag(sc, dimn, dimm, ws, smat_lin_dep_threshold)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'after A subspace S matrix rdiag, new dimension is', dimm

            If (dimm == 0) then
                Call memminus(KIND(indsym), SIZE(indsym), 1); deallocate (indsym)
                Call memminus(KIND(sc0), SIZE(sc0), 2); deallocate (sc0)
                Call memminus(KIND(sc), SIZE(sc), 2); deallocate (sc)
                Call memminus(KIND(ws), SIZE(ws), 1); deallocate (ws)
                cycle ! Go to the next isym.
            End if

            Allocate (bc(dimn, dimn)); Call memplus(KIND(bc), SIZE(bc), 2)   ! br N*N
            bc = 0.0d+00
            Call bAmat_real(dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Call memminus(KIND(sc0), SIZE(sc0), 2); deallocate (sc0)

            Allocate (uc(dimn, dimm)); Call memplus(KIND(uc), SIZE(uc), 2)           ! uc N*M
            Allocate (wsnew(dimm)); Call memplus(KIND(wsnew), SIZE(wsnew), 1)     ! wnew M
            uc(:, :) = 0.0d+00
            wsnew(:) = 0.0d+00

            Call rcutoff(sc, ws, dimn, dimm, smat_lin_dep_threshold, uc, wsnew)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Call memminus(KIND(sc), SIZE(sc), 2); deallocate (sc)
            Call memminus(KIND(ws), SIZE(ws), 1); deallocate (ws)

            Call ulambda_s_half(uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Call memminus(KIND(wsnew), SIZE(wsnew), 1); deallocate (wsnew)

            Allocate (bc0(dimm, dimn)); Call memplus(KIND(bc0), SIZE(bc0), 2) ! bc0 M*N
            bc0 = 0.0d+00
            bc0 = MATMUL(TRANSPOSE(uc), bc)

            Allocate (bc1(dimm, dimm)); Call memplus(KIND(bc1), SIZE(bc1), 2) ! bc1 M*M
            bc1 = 0.0d+00
            bc1 = MATMUL(bc0, uc)

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

            Call memminus(KIND(bc), SIZE(bc), 2); deallocate (bc)
            Call memminus(KIND(bc0), SIZE(bc0), 2); deallocate (bc0)

            Allocate (wb(dimm)); Call memplus(KIND(wb), SIZE(wb), 1)

            if (rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'

            Allocate (bc0(dimm, dimm)); Call memplus(KIND(bc0), SIZE(bc0), 2) ! bc0 M*M
            bc0 = bc1

            Call rdiag(bc1, dimm, dammy, wb, bmat_no_cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            Call memminus(KIND(bc0), SIZE(bc0), 2); deallocate (bc0)

            if (rank == 0) print *, 'bC1 matrix is diagonalized!'

            e2 = 0.0d+00

            Do ii = 1, ninact
                sym1 = irpamo(ii)
                if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. sym1 == isym)) then
                    Allocate (vc(dimn)); Call memplus(KIND(vc), SIZE(vc), 2)
                    Do it = 1, dimn
                        vc(it) = v(ii, indsym(1, it), indsym(2, it), indsym(3, it))
                    End do

                    Allocate (vc1(dimm)); Call memplus(KIND(vc1), SIZE(vc1), 2)
                    vc1(1:dimm) = MATMUL(TRANSPOSE(uc(1:dimn, 1:dimm)), vc(1:dimn))
                    Call memminus(KIND(vc), SIZE(vc), 2); Deallocate (vc)

                    alpha = -eps(ii) - e0 + eshift   ! For Level Shift (2007/2/9)

                    vc1(1:dimm) = MATMUL(TRANSPOSE(bc1(1:dimm, 1:dimm)), vc1(1:dimm))

                    Do j = 1, dimm
                        sumc2local = sumc2local + (ABS(vc1(j))**2.0d+00)/((alpha + wb(j))**2.0d+00)
                        e2(isym) = e2(isym) - (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                    End do
                    Call memminus(KIND(vc1), SIZE(vc1), 2); Deallocate (vc1)
                End if
            End do
            if (rank == 0) print '("e2a(",I3,") = ",E20.10," a.u.")', isym, e2(isym)

            Call memminus(KIND(bc1), SIZE(bc1), 2); Deallocate (bc1)
            Call memminus(KIND(uc), SIZE(uc), 2); Deallocate (uc)
            Call memminus(KIND(wb), SIZE(wb), 1); Deallocate (wb)
            Call memminus(KIND(indsym), SIZE(indsym), 2); Deallocate (indsym)

            e2a = e2a + e2(isym)
        End do

        if (rank == 0) then
            print '("e2a      = ",E20.10," a.u.")', e2a

            print '("sumc2,a  = ",E20.10)', sumc2local
        end if
        sumc2 = sumc2 + sumc2local

        Call memminus(KIND(v), SIZE(v), 2); Deallocate (v)

        if (rank == 0) print *, 'end solve_A_subspace'
    end subroutine solve_A_subspace_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
    SUBROUTINE sAmat_real(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space A
!
!  S(xyz,tuv) =  - <0|EzyEtxEuv|0> + d(tx)<0|EzyEuv|0>
!
!     x > y, t > u

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)      :: dimn, indsym(3, dimn)
        real(8), intent(out)  :: sc(dimn, dimn)
        real(8)  ::a, b
        integer :: it, iu, iv, ix, iy, iz
        integer :: i, j

        if (rank == 0) print *, 'Start A subspace S matrix'
! Initialization
        sc = 0.0d+00

!$OMP parallel do private(i,ix,iy,iz,j,it,iu,iv,a,b)
        Do i = rank + 1, dimn, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
            ix = indsym(1, i)
            iy = indsym(2, i)
            iz = indsym(3, i)

            Do j = i, dimn
                it = indsym(1, j)
                iu = indsym(2, j)
                iv = indsym(3, j)

!  S(xyz,tuv) =  - <0|EzyEtxEuv|0> + d(tx)<0|EzyEuv|0>
                Call dim3_density(iz, iy, it, ix, iu, iv, a, b)
                sc(i, j) = sc(i, j) - a

                If (it == ix) then
                    Call dim2_density(iz, iy, iu, iv, a, b)
                    sc(i, j) = sc(i, j) + a
                End if
                sc(j, i) = sc(i, j)
            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=sc)
#endif
        if (rank == 0) print *, 'A subspace S matrix is obtained normally'
    End subroutine sAmat_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE bAmat_real(dimn, sc, indsym, bc) ! Assume C1 molecule, overlap matrix B in space A
!
!  B(xyz,tuv) = Siguma_w eps(w){-<0|EzyEtxEuvEww|0> + d(tx)<0|EzyEuvEww|0>}
!
!               + S(xyz,tuv)(eps(u)+eps(t)-eps(v))
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer :: it, iu, iv, ix, iy, iz, iw
        integer :: jt, ju, jv, jx, jy, jz, jw
        integer :: i, j

        integer, intent(in)     :: dimn, indsym(3, dimn)
        real(8), intent(in)  :: sc(dimn, dimn)
        real(8), intent(out) :: bc(dimn, dimn)

        real(8)               :: e, denr, deni
        real(8)           :: den

        if (rank == 0) print *, 'Start A subspace B matrix'
        bc(:, :) = 0.0d+00
!$OMP parallel do private(ix,iy,iz,jx,jy,jz,it,iu,iv,jt,ju,jv,e,j,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
            ix = indsym(1, i)
            iy = indsym(2, i)
            iz = indsym(3, i)
            jx = convert_active_to_global_idx(ix)
            jy = convert_active_to_global_idx(iy)
            jz = convert_active_to_global_idx(iz)

            Do j = i, dimn

                it = indsym(1, j)
                iu = indsym(2, j)
                iv = indsym(3, j)
                jt = convert_active_to_global_idx(it)
                ju = convert_active_to_global_idx(iu)
                jv = convert_active_to_global_idx(iv)

!  B(xyz,tuv) = Siguma_w eps(w){-<0|EzyEtxEuvEww|0> + d(tx)<0|EzyEuvEww|0>}
!
!               + S(xyz,tuv)(eps(u)+eps(t)-eps(v))

                e = eps(ju) + eps(jt) - eps(jv)

                Do iw = 1, nact
                    jw = convert_active_to_global_idx(iw)

                    Call dim4_density(iz, iy, it, ix, iu, iv, iw, iw, denr, deni)
                    den = denr
                    bc(i, j) = bc(i, j) - den*eps(jw)

                    If (it == ix) then
                        Call dim3_density(iz, iy, iu, iv, iw, iw, denr, deni)
                        den = denr
                        bc(i, j) = bc(i, j) + den*eps(jw)

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

        if (rank == 0) print *, 'A subspace B matrix is obtained normally'
    End subroutine bAmat_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
    SUBROUTINE vAmat_real(v)
!
! Assume C1 molecule, V=<0|H|i> matrix in space A
!
!  V(tuv,i)=  - SIGUMA_p,q,r:act <0|EvuEptEqr|0>(pi|qr)
!
!  + SIGUMA_p,q:act  <0|EvuEpq|0> (pq|ti)
!
!  - SIGUMA_p:act <0|EvuEpt|0>[h(pi)+ SIGUMA_k:inact{(pi|kk)-(pk|ki)}]
!
!  +  <0|Evu|0>[h(ti) + SIGUMA_k:inact{(ti|kk) - (tk|ki)}]
!
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

        real(8), intent(out) :: v(ninact, nact, nact, nact)
        real(8)                  :: dr, di
        real(8)              :: cint2, d, dens1(nact, nact), effh(nact, ninact)
        complex*16              :: cint1
        logical                 :: is_end_of_file

        integer :: it, iu, iv, ii, ip
        integer :: jt, ju, jv, ji, jp
        integer :: i, j, k, l, dim(nsymrpa)
        integer :: dim2(nsymrpa), isym, i0, syma, symb, symc, iostat, unit_int2
        integer, allocatable :: indt(:, :), indu(:, :), indv(:, :)
        integer, allocatable :: ind2u(:, :), ind2v(:, :)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  V(tuv,i)=  - SIGUMA_p,q,r:act <0|EvuEptEqr|0>(pi|qr)
!
!  + SIGUMA_p,q:act  <0|EvuEpq|0> (pq|ti)
!
!  - SIGUMA_p:act <0|EvuEpt|0>[h(pi)+ SIGUMA_k:inact{(pi|kk)-(pk|ki)}]
!
!  +  <0|Evu|0>[h(ti) + SIGUMA_k:inact{(ti|kk) - (tk|ki)}]
!
!   = - SIGUMA_p,q,r:act <0|EvuEptEqr|0>(pi|qr)
!
!  + SIGUMA_p,q:act  <0|EvuEpq|0> (pq|ti)
!
!  - SIGUMA_p:act <0|EvuEpt|0>effh(pi)  +  <0|Evu|0>effh(ti)
!                 ========================================= This part is calculated after reading int2
!
!  effh is stored in memory while reading int2.
!
!  effh(p,i) = h(pi)+ SIGUMA_k:inact{(pi|kk)-(pk|ki)}
        if (rank == 0) print *, 'Start A subspace V matrix'
        v = 0.0d+00
        dens1 = 0.0d+00
        effh = 0.0d+00
        dim = 0

        Allocate (indt(nact**3, nsymrpa)); Call memplus(KIND(indt), SIZE(indt), 1)
        Allocate (indu(nact**3, nsymrpa)); Call memplus(KIND(indu), SIZE(indu), 1)
        Allocate (indv(nact**3, nsymrpa)); Call memplus(KIND(indv), SIZE(indv), 1)
        indt = 0
        indu = 0
        indv = 0
        dim = 0
!$OMP parallel do schedule(static) private(it,jt,iv,jv,iu,ju,syma,symb,symc)
        Do isym = 1, nsymrpa
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Do iv = 1, nact
                    jv = convert_active_to_global_idx(iv)
                    Do iu = 1, nact
                        ju = convert_active_to_global_idx(iu)
! EtiEuv
                        if (nsymrpa /= 1) then
                            syma = MULTB_D(irpamo(jt), isym)
                            symb = MULTB_D(irpamo(ju), irpamo(jv))
                            symc = MULTB_S(syma, symb)
                        end if
                        if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. symc == 1)) then
                            dim(isym) = dim(isym) + 1
                            indt(dim(isym), isym) = it
                            indu(dim(isym), isym) = iu
                            indv(dim(isym), isym) = iv
                        end if
                    End do
                End do
            End do
        End do
!$OMP end parallel do
        Allocate (ind2u(nact**2, nsymrpa)); Call memplus(KIND(ind2u), SIZE(ind2u), 1)
        Allocate (ind2v(nact**2, nsymrpa)); Call memplus(KIND(ind2v), SIZE(ind2v), 1)
        ind2u = 0.0d+00
        ind2v = 0.0d+00
        dim2 = 0
!$OMP parallel do schedule(static) private(iu,ju,iv,jv,syma)
        Do isym = 1, nsymrpa
            Do iu = 1, nact
                ju = convert_active_to_global_idx(iu)
                Do iv = 1, nact
                    jv = convert_active_to_global_idx(iv)
                    if (nsymrpa /= 1) syma = MULTB_D(irpamo(jv), irpamo(ju))
                    if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                        dim2(isym) = dim2(isym) + 1
                        ind2u(dim2(isym), isym) = iu
                        ind2v(dim2(isym), isym) = iv
                    end if
                End do
            End do
        End do
!$OMP end parallel do

!$OMP parallel do private(ji,it,jt,cint1)
        Do ii = rank + 1, ninact, nprocs
            ji = ii
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Call tramo1(jt, ji, cint1)
                effh(it, ii) = real(cint1, kind=KIND(effh))
            End do
        End do
!$OMP end parallel do

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Two types of integrals are stored
!
!  (21|22) stored (pi|qr) ...TYPE 1
!  (21|11) stored (pi|jk) ...TYPE 2
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2, file=a1int, status='old', optional_action='read')
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
            call check_iostat(iostat=iostat, file=a1int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  V(tuv,i)=  - SIGUMA_p,q,r:act <0|EvuEptEqr|0>(pi|qr)
!
!  + SIGUMA_p,q:act  <0|EvuEpq|0> (ti|pq)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            isym = irpamo(j)
!$OMP parallel do private(it,iu,iv,jt,ju,jv,dr,di,d)
            Do i0 = 1, dim(isym)
                it = indt(i0, isym)
                iu = indu(i0, isym)
                iv = indv(i0, isym)
                jt = convert_active_to_global_idx(it)
                ju = convert_active_to_global_idx(iu)
                jv = convert_active_to_global_idx(iv)

                Call dim3_density(iv, iu, i, it, k, l, dr, di)
                d = dr
                v(j, it, iu, iv) = v(j, it, iu, iv) - cint2*d
            End do
!$OMP end parallel do

            isym = MULTB_D(irpamo(convert_active_to_global_idx(i)), irpamo(j))           ! j coresponds to ii, i coresponds to it

!$OMP parallel do private(iu,iv,ju,jv,dr,di,d)
            Do i0 = 1, dim2(isym)
                iu = ind2u(i0, isym)
                iv = ind2v(i0, isym)
                ju = convert_active_to_global_idx(iu)
                jv = convert_active_to_global_idx(iv)

                Call dim2_density(iv, iu, k, l, dr, di)
                d = dr
                v(j, i, iu, iv) = v(j, i, iu, iv) + cint2*d
            End do
!$OMP end parallel do

        end do

        close (unit_int2)
        if (rank == 0) print *, 'reading A1int2 is over'

        call open_unformatted_file(unit=unit_int2, file=a2int, status='old', optional_action='read') ! TYPE 2 integrals
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
            call check_iostat(iostat=iostat, file=a2int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  effh(p,i) = h(pi)+ SIGUMA_k:inact{(pi|kk)-(pk|ki)}
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (k == l .and. j /= k) then       ! (PI|KK) type
                effh(i, j) = effh(i, j) + cint2
            elseif (j == k .and. k /= l) then       ! (PK|KI) type
                effh(i, l) = effh(i, l) - cint2
            end if
        end do

        close (unit_int2)
        if (rank == 0) print *, 'reading A2int2 is over'

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=effh)
#endif

!  - SIGUMA_p:act <0|EvuEpt|0>effh(pi)  +  <0|Evu|0>effh(ti)

!$OMP parallel do private(ji,isym,it,iu,iv,dr,di,d,ip,jp)
        Do ii = rank + 1, ninact, nprocs
            ji = ii
            isym = irpamo(ji)

            Do i0 = 1, dim(isym)
                it = indt(i0, isym)
                iu = indu(i0, isym)
                iv = indv(i0, isym)

                Call dim1_density(iv, iu, dr, di)
                d = dr
                v(ii, it, iu, iv) = v(ii, it, iu, iv) + effh(it, ii)*d

                Do ip = 1, nact
                    jp = convert_active_to_global_idx(ip)

                    Call dim2_density(iv, iu, ip, it, dr, di)
                    d = dr
                    v(ii, it, iu, iv) = v(ii, it, iu, iv) - effh(ip, ii)*d
                End do
            End do
        End do
!$OMP end parallel do

        ! deallocate memory
        Call memminus(KIND(indt), SIZE(indt), 1); deallocate (indt)
        Call memminus(KIND(indu), SIZE(indu), 1); deallocate (indu)
        Call memminus(KIND(indv), SIZE(indv), 1); deallocate (indv)
        Call memminus(KIND(ind2u), SIZE(ind2u), 1); deallocate (ind2u)
        Call memminus(KIND(ind2v), SIZE(ind2v), 1); deallocate (ind2v)

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
#endif
        if (rank == 0) print *, 'A subspace V matrix is obtained normally'
    end subroutine vAmat_real
end SUBROUTINE solve_A_subspace
