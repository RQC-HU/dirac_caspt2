SUBROUTINE solve_C_subspace(e0)

    use dcaspt2_restart_file, only: get_subspace_idx
    use module_ulambda_s_half, only: ulambda_s_half
    use module_global_variables
    use module_realonly, only: realonly
    use module_time
    implicit none
    real(8), intent(in) :: e0
    integer :: subspace_idx

    subspace_idx = get_subspace_idx('C')
    if (realonly%is_realonly()) then
        call solve_C_subspace_real()
    else
        call solve_C_subspace_complex()
    end if
    e2all = e2all + e2_subspace(subspace_idx)
    sumc2 = sumc2 + sumc2_subspace(subspace_idx)
contains
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_C_subspace_complex()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx, convert_secondary_to_global_idx

        Implicit NONE

        integer :: dimn, dimm, dammy

        integer, allocatable :: indsym(:, :)

        real(8), allocatable  :: wsnew(:), ws(:), wb(:)
        real(8)               :: e2(nsymrpa), alpha

        complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
        complex*16, allocatable  :: bc(:, :)
        complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :, :), vc(:), vc1(:)

        integer :: j, i, syma, symb, isym
        integer :: ix, iy, iz, ia, dima, ixyz
        integer :: jx, jy, jz, ja, it

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE C IS NOW CALCULATED
!
!     EatEuv|0>
!
!   DRAS1 = 0   DRAS2 = -1   DRAS3 = +1
!
!
!  S(xyz,tuv) = <0|EzyExtEuv|0>
!
!  B(xyz,tuv) = Siguma_w [eps(w)<0|EzyExtEuvEww|0>+S(xyz,tuv)(eps(u)-eps(v)-eps(t))]
!
!  a(a)       = eps(a) - Siguma_w [eps(w)<0|Eww|0>] (<== calculated as e0 in calce0.f)
!
!  V(tuv,a)   = Siguma_p [h'ap - Siguma_q(aq|qp)]<0|EvuEtp|0> + Siguma_pqr<0|EvuEtrEpq|0>(ar|pq)
!
!  where h'ap = hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)]
!
! Indices are restricted as t > v, x > z
! So the dimension of (xyz) is (norb**3+norb**2)/2
!
!  E2 = SIGUMA_a, dimm |V1(dimm,a)|^2|/{(a(a) + wb(dimm)}

        e2 = 0.0d+00
        dima = 0
        dimn = 0
        syma = 0
        if (debug .and. rank == 0) print *, 'ENTER solve C part'
        if (rank == 0) print '(10A)', '  '
        if (rank == 0) print '(10A)', ' e2c(isym)'

        Allocate (v(nsec, nact, nact, nact))
        Call vCmat_complex(v)
        Do isym = 1, nsymrpa

            ixyz = 0
!     EatEuv|0>
!     EaxEyz|0>

            Do ix = 1, nact
                Do iy = 1, nact
                    Do iz = 1, nact

                        jx = convert_active_to_global_idx(ix)
                        jy = convert_active_to_global_idx(iy)
                        jz = convert_active_to_global_idx(iz)
                        if (nsymrpa /= 1) then
                            syma = MULTB_D(isym, irpamo(jx))
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

            If (dimn == 0) cycle ! Go to the next isym

            Allocate (indsym(3, dimn))
            indsym = 0
            ixyz = 0

            Do ix = 1, nact
                Do iy = 1, nact
                    Do iz = 1, nact

                        jx = convert_active_to_global_idx(ix)
                        jy = convert_active_to_global_idx(iy)
                        jz = convert_active_to_global_idx(iz)
                        if (nsymrpa /= 1) then
                            syma = MULTB_D(isym, irpamo(jx))
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

            if (debug .and. rank == 0) print *, 'isym, dimn', isym, dimn
            Allocate (sc(dimn, dimn))
            sc = 0.0d+00            ! sr N*N
            Call sCmat_complex(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            Allocate (ws(dimn))
            ws = 0.0d+00

            Allocate (sc0(dimn, dimn))
            sc0 = 0.0d+00
            sc0 = sc
            Call cdiag(sc, dimn, dimm, ws, smat_lin_dep_threshold)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (debug .and. rank == 0) print *, 'after C subspace S matrix cdiag, new dimension is', dimm
            If (dimm == 0) then
                deallocate (indsym)
                deallocate (sc0)
                deallocate (sc)
                deallocate (ws)
                cycle ! Go to the next isym
            End if

            If (debug) then
                if (debug .and. rank == 0) print *, 'Check whether U*SU is diagonal'
                Call checkdgc(dimn, sc0, sc, ws)
                if (debug .and. rank == 0) print *, 'Check whether U*SU is diagonal END'
            End if
            Allocate (bc(dimn, dimn))                                 ! br N*N
            bc = 0.0d+00
            Call bCmat_complex(dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            deallocate (sc0)

            Allocate (uc(dimn, dimm))                                 ! uc N*M
            Allocate (wsnew(dimm))                                    ! wnew M
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
            bc0 = 0.0d+00
            bc0 = MATMUL(TRANSPOSE(DCONJG(uc)), bc)
            Allocate (bc1(dimm, dimm))                      ! bc1 M*M
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

            deallocate (bc)
            deallocate (bc0)

            Allocate (wb(dimm))
            wb = 0.0d+00
            if (debug .and. rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'
            Allocate (bc0(dimm, dimm))
            bc0 = 0.0d+00
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
            Do ia = 1, nsec
                ja = convert_secondary_to_global_idx(ia)

                if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. irpamo(ja) == isym)) then

                    Allocate (vc(dimn))
                    Do it = 1, dimn
                        vc(it) = v(ia, indsym(1, it), indsym(2, it), indsym(3, it))
                    End do

                    Allocate (vc1(dimm))
                    vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn, 1:dimm))), vc(1:dimn))

                    Deallocate (vc)

                    alpha = eps(ja) - e0 + eshift   ! For Level Shift (2007/2/9)

                    vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm, 1:dimm))), vc1(1:dimm))

                    Do j = 1, dimm
                        sumc2_subspace(subspace_idx) = sumc2_subspace(subspace_idx) + &
                                                       (ABS(vc1(j))**2.0d+00)/((alpha + wb(j))**2.0d+00)
                        e2(isym) = e2(isym) - (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                    End do
                    Deallocate (vc1)

                End if

            End do

            if (rank == 0) print '(" e2c(",I3,") = ",E25.15," a.u.")', isym, e2(isym)
            deallocate (bc1)
            deallocate (indsym)
            Deallocate (uc)
            Deallocate (wb)

            e2_subspace(subspace_idx) = e2_subspace(subspace_idx) + e2(isym)
        End do

        if (rank == 0) then
            print '(" e2c      = ",E25.15," a.u.")', e2_subspace(subspace_idx)
            print '(" sumc2,c  = ",E25.15)', sumc2_subspace(subspace_idx)
        end if

        continue
        if (debug .and. rank == 0) print *, 'end solve_C_subspace'
    end subroutine solve_C_subspace_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE sCmat_complex(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space C

!  S(xyz,tuv) = <0|EzyExtEuv|0>
!     x > z, t > v

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

        if (debug .and. rank == 0) print *, 'Start C subspace S matrix'
        sc = 0.0d+00

!$OMP parallel do schedule(dynamic,1) private(ix,iy,iz,it,iu,iv,a,b)
        Do i = rank + 1, dimn, nprocs
            ix = indsym(1, i)
            iy = indsym(2, i)
            iz = indsym(3, i)

            Do j = i, dimn

                it = indsym(1, j)
                iu = indsym(2, j)
                iv = indsym(3, j)

                a = 0.0d+0
                b = 0.0d+0

                Call dim3_density(iz, iy, ix, it, iu, iv, a, b)

                sc(i, j) = DCMPLX(a, b)
                sc(j, i) = DCMPLX(a, -b)
                if (rank == 0) then
                    If (ABS(sc(i, j)) > 1.0d+00) then
                        print '(2I4,2E25.15)', i, j, sc(i, j)
                    End if
                end if
            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=sc)
#endif
        if (debug .and. rank == 0) print *, 'C subspace S matrix is obtained normally'
    End subroutine sCmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE bCmat_complex(dimn, sc, indsym, bc)
!
! Indices are restricted as t > v, x > z
! So the dimension of (xyz) is (norb**3+norb**2)/2
!
!  S(xyz,tuv) = <0|EzyExtEuv|0>
!
!  B(xyz,tuv) = Siguma_w [eps(w)<0|EzyExtEuvEww|0>+S(xyz,tuv)(eps(u)-eps(v)-eps(t))]
!
!  a(a)       = eps(a) - Siguma_w [eps(w)<0|Eww|0>]
!
! H0-ES = B-aS : a is iependent from the index of active orbital like, x, y, z, and so on
!
! Here B matrix is constructed.
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer :: it, iu, iv, ix, iy, iz, iw, i, j
        integer :: jt, ju, jv, jx, jy, jz, jw

        integer, intent(in) :: dimn, indsym(3, dimn)
        complex*16, intent(in)  :: sc(dimn, dimn)
        complex*16, intent(out) :: bc(dimn, dimn)

        real(8)              :: e, denr, deni
        complex*16          :: den

        if (debug .and. rank == 0) print *, 'Start C subspace B matrix'
        bc(:, :) = 0.0d+00

!$OMP parallel do schedule(dynamic,1) private(ix,iy,iz,jx,jy,jz,it,iu,iv,jt,ju,jv,e,j,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs
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

!   Siguma_w [eps(w)<0|EzyExtEuvEww|0>]+S(xyz,tuv)(eps(u)-eps(v)-eps(t))

                e = eps(ju) - eps(jv) - eps(jt)

                Do iw = 1, nact
                    jw = convert_active_to_global_idx(iw)
                    Call dim4_density(iz, iy, ix, it, iu, iv, iw, iw, denr, deni)
                    den = DCMPLX(denr, deni)
                    bc(i, j) = bc(i, j) + den*eps(jw)
                End do

                bc(i, j) = bc(i, j) + sc(i, j)*e

                bc(j, i) = DCONJG(bc(i, j))

            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call reduce_wrapper(mat=bc, root_rank=0)
#endif
        if (debug .and. rank == 0) print *, 'C subspace B matrix is obtained normally'
    End subroutine bCmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE vCmat_complex(v)

! Assume C1 molecule, V=<0|H|i> matrix in space C
!
!     EatEuv|0>
!
!  V(a,tuv)   = Siguma_p [h'ap - Siguma_w(aw|wp)]<0|EvuEtp|0> + Siguma_pqr<0|EvuEtrEpq|0>(ar|pq)
!
!  where h'ap = hap + Siguma_k(runs inactive (and frozen) orbital)[(ap|kk)-(ak|kp)]
!
! Indices are restricted as t > v, x > z
!
! So the dimension of (xyz) is (norb**3+norb**2)/2 ! <= C1 case
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_file_manager, only: open_unformatted_file, check_iostat
        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx, convert_secondary_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        complex*16, intent(out) :: v(nsec, nact, nact, nact)
        real(8)                  :: dr, di
        complex*16               :: cint1
        complex*16              :: cint2, d
        complex*16              :: effh(nsec, nact)
        integer :: i, j, k, l, dim(nsymrpa)
        integer :: isym, syma, symb, symc
        integer, allocatable :: indt(:, :), indu(:, :), indv(:, :)
        integer :: it, iu, iv, ia, ip
        integer :: jt, ju, jv, ja
        integer :: i0, iostat, unit_int2
        logical :: is_end_of_file
!^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
!  V(a,t,u,v)   = Siguma_p [h'ap - Siguma_w(aw|wp)]<0|EvuEtp|0> + Siguma_pqr<0|EvuEtrEpq|0>(ar|pq)
!
! All indices run active spinor space except below k(inactive).
!
!  where h'ap = hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)]
!
!
!  effh(a,p) =  hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)] - Siguma_w(aw|wp)
!
!  is calculated and stored in memory and after reading int2, take count in V(a,t,u,v)
!
!
!===============================================!
! Three types of integrals are stored Cint      !
!                                               !
! (ap|qr) = (32|22) TYPE 1 (includes (aw|wp) )  !
!                                               !
! (ap|kk) = (32|11) TYPE 2                      !
!                                               !
! (ak|kp) = (31|12) TYPE 3                      !
!===============================================!
!
!^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~

        if (debug .and. rank == 0) print *, 'Start C subspace V matrix'

        v = 0.0d+00
        effh = 0.0d+00
        dim = 0

        Allocate (indt(nact**3, nsymrpa))
        Allocate (indu(nact**3, nsymrpa))
        Allocate (indv(nact**3, nsymrpa))
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

!     EatEuv|0>
                        if (nsymrpa /= 1) then
                            syma = MULTB_D(irpamo(ju), irpamo(jv))
                            symb = MULTB_D(isym, irpamo(jt))
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

!$OMP parallel do schedule(dynamic,1) private(ia,ja,it,jt,cint1)
        Do ia = rank + 1, nsec, nprocs
            ja = convert_secondary_to_global_idx(ia)
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)

                Call tramo1(ja, jt, cint1)
                effh(ia, it) = cint1
            End do
        End do
!$OMP end parallel do
        call open_unformatted_file(unit=unit_int2, file=c1int, status='old', optional_action='read')
        do ! Read TYPE 1 integrals C1int until EOF
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
            call check_iostat(iostat=iostat, file=c1int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! + Siguma_pqr<0|EvuEtrEpq|0>(ar|pq)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            isym = irpamo(convert_secondary_to_global_idx(i))   ! i corresponds to a
!$OMP parallel do schedule(static,1) private(it,iu,iv,dr,di,d)
            Do i0 = 1, dim(isym)
                it = indt(i0, isym)
                iu = indu(i0, isym)
                iv = indv(i0, isym)

                Call dim3_density(iv, iu, it, j, k, l, dr, di)
                d = DCMPLX(dr, di)
                v(i, it, iu, iv) = v(i, it, iu, iv) + cint2*d

            End do
!$OMP end parallel do
!  effh(a,p) =  hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)] - Siguma_w(aw|wp)
!                                                           ~~~~~~~~~~~~~~~~~~~
            if (j == k) then
                effh(i, l) = effh(i, l) - cint2
            end if
        end do
        close (unit_int2)
        if (debug .and. rank == 0) print *, 'reading C1int2 is over'

        call open_unformatted_file(unit=unit_int2, file=c2int, status='old', optional_action='read')
        do ! Read TYPE 2 integrals C2int until EOF
            read (unit_int2, iostat=iostat) i, j, k, l, cint2
            call check_iostat(iostat=iostat, file=c2int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  effh(a,p) =  hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)] - Siguma_w(aw|wp)
!                                          ========
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (k == l) then
                effh(i, j) = effh(i, j) + cint2
            end if
        end do
        close (unit_int2)

        if (debug .and. rank == 0) print *, 'reading C2int2 is over'

        call open_unformatted_file(unit=unit_int2, file=c3int, status='old', optional_action='read') ! TYPE 3 integrals
        do ! Read TYPE 3 integrals C3int until EOF
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl):=> (ak|kp)
            call check_iostat(iostat=iostat, file=c3int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  effh(a,p) =  hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)] - Siguma_w(aw|wp)
!                                                  =========
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (j == k) then
                effh(i, l) = effh(i, l) - cint2
            end if

        end do
        close (unit_int2)

        if (debug .and. rank == 0) print *, 'reading C3int2 is over'

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=effh)
#endif

! Siguma_p effh(a,p)<0|EvuEtp|0>
!$OMP parallel do schedule(dynamic,1) private(ja,isym,i0,it,iu,iv,jt,ju,jv,dr,di,d)
        Do ia = rank + 1, nsec, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
            ja = convert_secondary_to_global_idx(ia)
            isym = irpamo(ja)

            Do ip = 1, nact

! Go to the next ip if the value of effh(ja,jp) is nearly zero
                if (ABS(effh(ia, ip)) < 1.0d-10) cycle

                Do i0 = 1, dim(isym)
                    it = indt(i0, isym)
                    iu = indu(i0, isym)
                    iv = indv(i0, isym)
                    jt = convert_active_to_global_idx(it)
                    ju = convert_active_to_global_idx(iu)
                    jv = convert_active_to_global_idx(iv)

                    Call dim2_density(iv, iu, it, ip, dr, di)
                    d = DCMPLX(dr, di)

                    v(ia, it, iu, iv) = v(ia, it, iu, iv) + effh(ia, ip)*d

                End do

            End do
        End do
!$OMP end parallel do

        deallocate (indt)
        deallocate (indu)
        deallocate (indv)

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
#endif
        if (debug .and. rank == 0) print *, 'C subspace V matrix is obtained normally'
    end subroutine vCmat_complex
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_C_subspace_real()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx, convert_secondary_to_global_idx

        Implicit NONE

        integer :: dimn, dimm, dammy

        integer, allocatable :: indsym(:, :)

        real(8), allocatable  :: wsnew(:), ws(:), wb(:)
        real(8)               :: e2(nsymrpa), alpha

        real(8), allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
        real(8), allocatable  :: bc(:, :)
        real(8), allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :, :), vc(:), vc1(:)

        integer :: j, i, syma, symb, isym
        integer :: ix, iy, iz, ia, dima, ixyz
        integer :: jx, jy, jz, ja, it

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE C IS NOW CALCULATED
!
!     EatEuv|0>
!
!   DRAS1 = 0   DRAS2 = -1   DRAS3 = +1
!
!
!  S(xyz,tuv) = <0|EzyExtEuv|0>
!
!  B(xyz,tuv) = Siguma_w [eps(w)<0|EzyExtEuvEww|0>+S(xyz,tuv)(eps(u)-eps(v)-eps(t))]
!
!  a(a)       = eps(a) - Siguma_w [eps(w)<0|Eww|0>] (<== calculated as e0 in calce0.f)
!
!  V(tuv,a)   = Siguma_p [h'ap - Siguma_q(aq|qp)]<0|EvuEtp|0> + Siguma_pqr<0|EvuEtrEpq|0>(ar|pq)
!
!  where h'ap = hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)]
!
! Indices are restricted as t > v, x > z
! So the dimension of (xyz) is (norb**3+norb**2)/2
!
!  E2 = SIGUMA_a, dimm |V1(dimm,a)|^2|/{(a(a) + wb(dimm)}

        e2 = 0.0d+00
        dima = 0
        dimn = 0
        syma = 0
        if (debug .and. rank == 0) print *, 'ENTER solve C part'
        if (rank == 0) print '(10A)', '  '
        if (rank == 0) print '(10A)', ' e2c(isym)'

        Allocate (v(nsec, nact, nact, nact))

        Call vCmat_real(v)
        Do isym = 1, nsymrpa

            ixyz = 0
!     EatEuv|0>
!     EaxEyz|0>

            Do ix = 1, nact
                Do iy = 1, nact
                    Do iz = 1, nact

                        jx = convert_active_to_global_idx(ix)
                        jy = convert_active_to_global_idx(iy)
                        jz = convert_active_to_global_idx(iz)
                        if (nsymrpa /= 1) then
                            syma = MULTB_D(isym, irpamo(jx))
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

            If (dimn == 0) cycle ! Go to the next isym

            Allocate (indsym(3, dimn))
            indsym = 0
            ixyz = 0

            Do ix = 1, nact
                Do iy = 1, nact
                    Do iz = 1, nact

                        jx = convert_active_to_global_idx(ix)
                        jy = convert_active_to_global_idx(iy)
                        jz = convert_active_to_global_idx(iz)
                        if (nsymrpa /= 1) then
                            syma = MULTB_D(isym, irpamo(jx))
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

            if (debug .and. rank == 0) print *, 'isym, dimn', isym, dimn
            Allocate (sc(dimn, dimn))
            sc = 0.0d+00            ! sr N*N
            Call sCmat_real(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            Allocate (ws(dimn))
            ws = 0.0d+00

            Allocate (sc0(dimn, dimn))
            sc0 = 0.0d+00
            sc0 = sc
            Call rdiag(sc, dimn, dimm, ws, smat_lin_dep_threshold)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (debug .and. rank == 0) print *, 'after C subspace S matrix rdiag, new dimension is', dimm
            If (dimm == 0) then
                deallocate (indsym)
                deallocate (sc0)
                deallocate (sc)
                deallocate (ws)
                cycle ! Go to the next isym
            End if

            Allocate (bc(dimn, dimn))                                 ! br N*N
            bc = 0.0d+00
            Call bCmat_real(dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            deallocate (sc0)

            Allocate (uc(dimn, dimm))                                 ! uc N*M
            Allocate (wsnew(dimm))                                    ! wnew M
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
            bc0 = 0.0d+00
            bc0 = MATMUL(TRANSPOSE(uc), bc)
            Allocate (bc1(dimm, dimm))                      ! bc1 M*M
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

            deallocate (bc)
            deallocate (bc0)

            Allocate (wb(dimm))
            wb = 0.0d+00
            if (debug .and. rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'
            Allocate (bc0(dimm, dimm))
            bc0 = 0.0d+00
            bc0 = bc1
            Call rdiag(bc1, dimm, dammy, wb, bmat_no_cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            deallocate (bc0)

            if (debug .and. rank == 0) print *, 'bC1 matrix is diagonalized!'
            Do ia = 1, nsec
                ja = convert_secondary_to_global_idx(ia)

                if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. irpamo(ja) == isym)) then

                    Allocate (vc(dimn))
                    Do it = 1, dimn
                        vc(it) = v(ia, indsym(1, it), indsym(2, it), indsym(3, it))
                    End do

                    Allocate (vc1(dimm))
                    vc1(1:dimm) = MATMUL(TRANSPOSE(uc(1:dimn, 1:dimm)), vc(1:dimn))

                    Deallocate (vc)

                    alpha = eps(ja) - e0 + eshift   ! For Level Shift (2007/2/9)

                    vc1(1:dimm) = MATMUL(TRANSPOSE(bc1(1:dimm, 1:dimm)), vc1(1:dimm))

                    Do j = 1, dimm
                        sumc2_subspace(subspace_idx) = sumc2_subspace(subspace_idx) + &
                                                       (ABS(vc1(j))**2.0d+00)/((alpha + wb(j))**2.0d+00)
                        e2(isym) = e2(isym) - (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                    End do
                    Deallocate (vc1)

                End if

            End do

            if (rank == 0) print '(" e2c(",I3,") = ",E25.15," a.u.")', isym, e2(isym)
            deallocate (bc1)
            deallocate (indsym)
            Deallocate (uc)
            Deallocate (wb)

            e2_subspace(subspace_idx) = e2_subspace(subspace_idx) + e2(isym)
        End do

        if (rank == 0) then
            print '(" e2c      = ",E25.15," a.u.")', e2_subspace(subspace_idx)
            print '(" sumc2,c  = ",E25.15)', sumc2_subspace(subspace_idx)
        end if

        continue
        if (debug .and. rank == 0) print *, 'end solve_C_subspace'
    end subroutine solve_C_subspace_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE sCmat_real(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space C

!  S(xyz,tuv) = <0|EzyExtEuv|0>
!     x > z, t > v

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

        if (debug .and. rank == 0) print *, 'Start C subspace S matrix'
        sc = 0.0d+00

!$OMP parallel do schedule(dynamic,1) private(ix,iy,iz,it,iu,iv,a,b)
        Do i = rank + 1, dimn, nprocs
            ix = indsym(1, i)
            iy = indsym(2, i)
            iz = indsym(3, i)

            Do j = i, dimn

                it = indsym(1, j)
                iu = indsym(2, j)
                iv = indsym(3, j)

                a = 0.0d+0
                b = 0.0d+0

                Call dim3_density(iz, iy, ix, it, iu, iv, a, b)

                sc(i, j) = a
                sc(j, i) = a
                if (rank == 0) then
                    If (ABS(sc(i, j)) > 1.0d+00) then
                        print '(2I4,2E25.15)', i, j, sc(i, j)
                    End if
                end if
            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=sc)
#endif
        if (debug .and. rank == 0) print *, 'C subspace S matrix is obtained normally'
    End subroutine sCmat_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE bCmat_real(dimn, sc, indsym, bc)
!
! Indices are restricted as t > v, x > z
! So the dimension of (xyz) is (norb**3+norb**2)/2
!
!  S(xyz,tuv) = <0|EzyExtEuv|0>
!
!  B(xyz,tuv) = Siguma_w [eps(w)<0|EzyExtEuvEww|0>+S(xyz,tuv)(eps(u)-eps(v)-eps(t))]
!
!  a(a)       = eps(a) - Siguma_w [eps(w)<0|Eww|0>]
!
! H0-ES = B-aS : a is iependent from the index of active orbital like, x, y, z, and so on
!
! Here B matrix is constructed.
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer :: it, iu, iv, ix, iy, iz, iw, i, j
        integer :: jt, ju, jv, jx, jy, jz, jw

        integer, intent(in) :: dimn, indsym(3, dimn)
        real(8), intent(in)  :: sc(dimn, dimn)
        real(8), intent(out) :: bc(dimn, dimn)

        real(8)              :: e, denr, deni
        real(8)          :: den

        if (debug .and. rank == 0) print *, 'Start C subspace B matrix'
        bc(:, :) = 0.0d+00

!$OMP parallel do schedule(dynamic,1) private(ix,iy,iz,jx,jy,jz,it,iu,iv,jt,ju,jv,e,j,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs
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

!   Siguma_w [eps(w)<0|EzyExtEuvEww|0>]+S(xyz,tuv)(eps(u)-eps(v)-eps(t))

                e = eps(ju) - eps(jv) - eps(jt)

                Do iw = 1, nact
                    jw = convert_active_to_global_idx(iw)
                    Call dim4_density(iz, iy, ix, it, iu, iv, iw, iw, denr, deni)
                    den = denr
                    bc(i, j) = bc(i, j) + den*eps(jw)
                End do

                bc(i, j) = bc(i, j) + sc(i, j)*e

                bc(j, i) = bc(i, j)

            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call reduce_wrapper(mat=bc, root_rank=0)
#endif
        if (debug .and. rank == 0) print *, 'C subspace B matrix is obtained normally'
    End subroutine bCmat_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE vCmat_real(v)

! Assume C1 molecule, V=<0|H|i> matrix in space C
!
!     EatEuv|0>
!
!  V(a,tuv)   = Siguma_p [h'ap - Siguma_w(aw|wp)]<0|EvuEtp|0> + Siguma_pqr<0|EvuEtrEpq|0>(ar|pq)
!
!  where h'ap = hap + Siguma_k(runs inactive (and frozen) orbital)[(ap|kk)-(ak|kp)]
!
! Indices are restricted as t > v, x > z
!
! So the dimension of (xyz) is (norb**3+norb**2)/2 ! <= C1 case
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_file_manager, only: open_unformatted_file, check_iostat
        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx, convert_secondary_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        real(8), intent(out) :: v(nsec, nact, nact, nact)
        real(8)                  :: dr, di
        complex*16               :: cint1
        real(8)              :: cint2, d
        real(8)              :: effh(nsec, nact)
        integer :: i, j, k, l, dim(nsymrpa)
        integer :: isym, syma, symb, symc
        integer, allocatable :: indt(:, :), indu(:, :), indv(:, :)
        integer :: it, iu, iv, ia, ip
        integer :: jt, ju, jv, ja
        integer :: i0, iostat, unit_int2
        logical :: is_end_of_file
!^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
!  V(a,t,u,v)   = Siguma_p [h'ap - Siguma_w(aw|wp)]<0|EvuEtp|0> + Siguma_pqr<0|EvuEtrEpq|0>(ar|pq)
!
! All indices run active spinor space except below k(inactive).
!
!  where h'ap = hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)]
!
!
!  effh(a,p) =  hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)] - Siguma_w(aw|wp)
!
!  is calculated and stored in memory and after reading int2, take count in V(a,t,u,v)
!
!
!===============================================!
! Three types of integrals are stored Cint      !
!                                               !
! (ap|qr) = (32|22) TYPE 1 (includes (aw|wp) )  !
!                                               !
! (ap|kk) = (32|11) TYPE 2                      !
!                                               !
! (ak|kp) = (31|12) TYPE 3                      !
!===============================================!
!
!^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~

        if (debug .and. rank == 0) print *, 'Start C subspace V matrix'

        v = 0.0d+00
        effh = 0.0d+00
        dim = 0

        Allocate (indt(nact**3, nsymrpa))
        Allocate (indu(nact**3, nsymrpa))
        Allocate (indv(nact**3, nsymrpa))
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

!     EatEuv|0>
                        if (nsymrpa /= 1) then
                            syma = MULTB_D(irpamo(ju), irpamo(jv))
                            symb = MULTB_D(isym, irpamo(jt))
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

!$OMP parallel do schedule(dynamic,1) private(ia,ja,it,jt,cint1)
        Do ia = rank + 1, nsec, nprocs
            ja = convert_secondary_to_global_idx(ia)
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)

                Call tramo1(ja, jt, cint1)
                effh(ia, it) = real(cint1, kind=KIND(effh))
            End do
        End do
!$OMP end parallel do
        call open_unformatted_file(unit=unit_int2, file=c1int, status='old', optional_action='read')
        do ! Read TYPE 1 integrals C1int until EOF
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
            call check_iostat(iostat=iostat, file=c1int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! + Siguma_pqr<0|EvuEtrEpq|0>(ar|pq)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            isym = irpamo(convert_secondary_to_global_idx(i))   ! i corresponds to a
!$OMP parallel do schedule(static,1) private(it,iu,iv,dr,di,d)
            Do i0 = 1, dim(isym)
                it = indt(i0, isym)
                iu = indu(i0, isym)
                iv = indv(i0, isym)

                Call dim3_density(iv, iu, it, j, k, l, dr, di)
                d = dr
                v(i, it, iu, iv) = v(i, it, iu, iv) + cint2*d

            End do
!$OMP end parallel do
!  effh(a,p) =  hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)] - Siguma_w(aw|wp)
!                                                           ~~~~~~~~~~~~~~~~~~~
            if (j == k) then
                effh(i, l) = effh(i, l) - cint2
            end if
        end do
        close (unit_int2)
        if (debug .and. rank == 0) print *, 'reading C1int2 is over'

        call open_unformatted_file(unit=unit_int2, file=c2int, status='old', optional_action='read')
        do ! Read TYPE 2 integrals C2int until EOF
            read (unit_int2, iostat=iostat) i, j, k, l, cint2
            call check_iostat(iostat=iostat, file=c2int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  effh(a,p) =  hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)] - Siguma_w(aw|wp)
!                                          ========
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (k == l) then
                effh(i, j) = effh(i, j) + cint2
            end if
        end do
        close (unit_int2)

        if (debug .and. rank == 0) print *, 'reading C2int2 is over'

        call open_unformatted_file(unit=unit_int2, file=c3int, status='old', optional_action='read') ! TYPE 3 integrals
        do ! Read TYPE 3 integrals C3int until EOF
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl):=> (ak|kp)
            call check_iostat(iostat=iostat, file=c3int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  effh(a,p) =  hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)] - Siguma_w(aw|wp)
!                                                  =========
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (j == k) then
                effh(i, l) = effh(i, l) - cint2
            end if

        end do
        close (unit_int2)

        if (debug .and. rank == 0) print *, 'reading C3int2 is over'

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=effh)
#endif

! Siguma_p effh(a,p)<0|EvuEtp|0>
!$OMP parallel do schedule(dynamic,1) private(ja,isym,i0,it,iu,iv,jt,ju,jv,dr,di,d)
        Do ia = rank + 1, nsec, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
            ja = convert_secondary_to_global_idx(ia)
            isym = irpamo(ja)

            Do ip = 1, nact

! Go to the next ip if the value of effh(ja,jp) is nearly zero
                if (ABS(effh(ia, ip)) < 1.0d-10) cycle

                Do i0 = 1, dim(isym)
                    it = indt(i0, isym)
                    iu = indu(i0, isym)
                    iv = indv(i0, isym)
                    jt = convert_active_to_global_idx(it)
                    ju = convert_active_to_global_idx(iu)
                    jv = convert_active_to_global_idx(iv)

                    Call dim2_density(iv, iu, it, ip, dr, di)
                    d = dr

                    v(ia, it, iu, iv) = v(ia, it, iu, iv) + effh(ia, ip)*d

                End do

            End do
        End do
!$OMP end parallel do

        deallocate (indt)
        deallocate (indu)
        deallocate (indv)

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
#endif
        if (debug .and. rank == 0) print *, 'C subspace V matrix is obtained normally'
    end subroutine vCmat_real
end subroutine solve_C_subspace
