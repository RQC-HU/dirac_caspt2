! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE solvC_ord_ty(e0, e2c)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    real*8, intent(in) :: e0
    real*8, intent(out):: e2c

    integer :: dimn, dimm, dammy

    integer, allocatable :: indsym(:, :)

    real*8, allocatable  :: wsnew(:), ws(:), wb(:)
    real*8               :: e2(nsymrp), alpha

    complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
    complex*16, allocatable  :: bc(:, :)
    complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :, :), vc(:), vc1(:)

    logical :: cutoff
    integer :: j, i, syma, symb, isym
    integer :: ix, iy, iz, ia, dima, ixyz
    integer :: jx, jy, jz, ja, it

    real*8  :: thresd
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1

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

    thresd = 1.0D-08
    thres = 1.0D-08

    e2 = 0.0d+00
    e2c = 0.0d+00
    dima = 0
    dimn = 0
    syma = 0
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0
    if (rank == 0) then
        print *, ' ENTER solv C part'
        print *, ' nsymrpa', nsymrpa
    end if
    Allocate (v(ninact + nact + 1:ninact + nact + nsec, ninact + 1:ninact + nact,  &
    &          ninact + 1:ninact + nact, ninact + 1:ninact + nact))
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) print *, 'end before v matrices'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Call vCmat_ord_ty(v)
    if (rank == 0) print *, 'come'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Do isym = 1, nsymrpa

        ixyz = 0

!     EatEuv|0>
!     EaxEyz|0>

        Do ix = 1, nact
            Do iy = 1, nact
                Do iz = 1, nact

                    jx = ix + ninact
                    jy = iy + ninact
                    jz = iz + ninact
                    syma = MULTB_D(isym, irpmo(jx))
                    symb = MULTB_D(irpmo(jy), irpmo(jz))
                    syma = MULTB_S(syma, symb)

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

                    jx = ix + ninact
                    jy = iy + ninact
                    jz = iz + ninact
                    syma = MULTB_D(isym, irpmo(jx))
                    symb = MULTB_D(irpmo(jy), irpmo(jz))
                    syma = MULTB_S(syma, symb)

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
        Allocate (sc(dimn, dimn))
        sc = 0.0d+00            ! sr N*N
        if (rank == 0) print *, 'before sCmat'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call sCmat(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (rank == 0) print *, 'sC matrix is obtained normally'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Allocate (ws(dimn))
        ws = 0.0d+00
        cutoff = .TRUE.
!           thresd = 1.0d-15

        Allocate (sc0(dimn, dimn))
        sc0 = 0.0d+00
        sc0 = sc
        if (rank == 0) print *, 'before cdiag'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call cdiag(sc, dimn, dimm, ws, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rank == 0) then
            print *, 'after sc cdiag'
            print *, 'after s cdiag, new dimension is', dimm
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        If (dimm == 0) then
            deallocate (indsym)
            deallocate (sc0)
            deallocate (sc)
            deallocate (ws)
            cycle ! Go to the next isym
        End if

        If (debug) then

            if (rank == 0) print *, 'Check whether U*SU is diagonal'

            Call checkdgc(dimn, sc0, sc, ws)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (rank == 0) print *, 'Check whether U*SU is diagonal END'
        End if
        if (rank == 0) print *, 'OK cdiag', dimn, dimm
        Allocate (bc(dimn, dimn))                                 ! br N*N
        bc = 0.0d+00
        if (rank == 0) print *, 'before bCmat'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call bCmat(dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        deallocate (sc0)

        if (rank == 0) print *, 'bC matrix is obtained normally'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Allocate (uc(dimn, dimm))                                 ! uc N*M
        Allocate (wsnew(dimm))                                    ! wnew M
        uc(:, :) = 0.0d+00
        wsnew(:) = 0.0d+00
        if (rank == 0) print *, 'before ccutoff'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call ccutoff(sc, ws, dimn, dimm, uc, wsnew)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (rank == 0) print *, 'OK ccutoff'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        deallocate (ws)
        deallocate (sc)
        if (rank == 0) print *, 'before ucramda_s_half'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call ucramda_s_half(uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        deallocate (wsnew)

        if (rank == 0) print *, 'ucrams half OK'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
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

        cutoff = .FALSE.

        Allocate (wb(dimm))
        wb = 0.0d+00
        if (rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'
        Allocate (bc0(dimm, dimm))
        bc0 = 0.0d+00
        bc0 = bc1
        if (rank == 0) print *, 'before cdiag'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call cdiag(bc1, dimm, dammy, wb, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rank == 0) print *, 'end cdiag'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        If (debug) then

            if (rank == 0) print *, 'Check whether bc is really diagonalized or not'
            Call checkdgc(dimm, bc0, bc1, wb)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'Check whether bc is really diagonalized or not END'
        End if

        deallocate (bc0)

        if (rank == 0) print *, 'bC1 matrix is diagonalized!'
        Do ia = 1, nsec
            ja = ia + ninact + nact

            if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. irpmo(ja) == isym)) then

                Allocate (vc(dimn))
                Do it = 1, dimn
                    vc(it) = v(ja, indsym(1, it) + ninact, indsym(2, it) + ninact, indsym(3, it) + ninact)
                End do

                Allocate (vc1(dimm))
                vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn, 1:dimm))), vc(1:dimn))

                Deallocate (vc)

                alpha = eps(ja) - e0 + eshift   ! For Level Shift (2007/2/9)

                vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm, 1:dimm))), vc1(1:dimm))

                Do j = 1, dimm
                    sumc2local = sumc2local + (ABS(vc1(j))**2.0d+00)/((alpha + wb(j))**2.0d+00)
                    e2(isym) = e2(isym) - (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                End do
                Deallocate (vc1)

            End if

        End do

        if (rank == 0) print '("e2c(",I3,") = ",E20.10,"a.u.")', isym, e2(isym)
        deallocate (bc1)
        deallocate (indsym)
        Deallocate (uc)
        Deallocate (wb)

        e2c = e2c + e2(isym)
        if (rank == 0) print *, 'End e2(isym) add'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
    End do                  ! isym

    if (rank == 0) then
        print '("e2c      = ",E20.10,"a.u.")', e2c
        print '("sumc2,c  = ",E20.10)', sumc2local
    end if
    sumc2 = sumc2 + sumc2local

    continue
    if (rank == 0) print *, 'end solvC_ord_ty'
end

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE sCmat(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space C

!  S(xyz,tuv) = <0|EzyExtEuv|0>
!     x > z, t > v

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer, intent(in)      :: dimn, indsym(3, dimn)
    complex*16, intent(out)  :: sc(dimn, dimn)
    real*8  ::a, b
    integer :: it, iu, iv, ix, iy, iz
    integer :: i, j

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
                    print '(2I4,2E20.10)', i, j, sc(i, j)
                End if
            end if
        End do               !j
    End do                  !i
    !$OMP end parallel do
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, sc(1, 1), dimn**2, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

End subroutine sCmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE bCmat(dimn, sc, indsym, bc)
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

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    integer :: it, iu, iv, ix, iy, iz, iw, i, j
    integer :: jt, ju, jv, jx, jy, jz, jw

    integer, intent(in) :: dimn, indsym(3, dimn)
    complex*16, intent(in)  :: sc(dimn, dimn)
    complex*16, intent(out) :: bc(dimn, dimn)

    real*8              :: e, denr, deni
    complex*16          :: den

    bc(:, :) = 0.0d+00

    if (rank == 0) print *, 'C space Bmat iroot=', iroot
    !$OMP parallel do schedule(dynamic,1) private(ix,iy,iz,jx,jy,jz,it,iu,iv,jt,ju,jv,e,j,iw,jw,denr,deni,den)
    Do i = rank + 1, dimn, nprocs
        ix = indsym(1, i)
        iy = indsym(2, i)
        iz = indsym(3, i)
        jx = ix + ninact
        jy = iy + ninact
        jz = iz + ninact

        Do j = i, dimn

            it = indsym(1, j)
            iu = indsym(2, j)
            iv = indsym(3, j)
            jt = it + ninact
            ju = iu + ninact
            jv = iv + ninact

!   Siguma_w [eps(w)<0|EzyExtEuvEww|0>]+S(xyz,tuv)(eps(u)-eps(v)-eps(t))

            e = eps(ju) - eps(jv) - eps(jt)

            Do iw = 1, nact
                jw = iw + ninact

                Call dim4_density(iz, iy, ix, it, iu, iv, iw, iw, denr, deni)
                den = DCMPLX(denr, deni)
                bc(i, j) = bc(i, j) + den*eps(jw)

            End do

            bc(i, j) = bc(i, j) + sc(i, j)*e

            bc(j, i) = DCONJG(bc(i, j))

        End do               !i
    End do                  !j
    !$OMP end parallel do
#ifdef HAVE_MPI
    if (rank == 0) then
        call MPI_Reduce(MPI_IN_PLACE, bc(1, 1), dimn**2, MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    else
        call MPI_Reduce(bc(1, 1), bc(1, 1), dimn**2, MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end if
#endif
    if (rank == 0) print *, 'bCmat is ended'
End subroutine bCmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE vCmat_ord_ty(v)

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

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    complex*16, intent(out) :: v(ninact + nact + 1:ninact + nact + nsec, ninact + 1:ninact + nact,  &
                &          ninact + 1:ninact + nact, ninact + 1:ninact + nact)
    real*8                  :: dr, di
    complex*16              :: cint1, cint2, d
    complex*16              :: effh(ninact + nact + 1:ninact + nact + nsec, ninact + 1:ninact + nact)
    integer :: i, j, k, l, dim(nsymrpa)
    integer :: isym, syma, symb, symc
    integer, allocatable :: indt(:, :), indu(:, :), indv(:, :)
    integer :: it, iu, iv, ia, ip
    integer :: jt, ju, jv, ja, jp
    integer :: i0, iostat
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1
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

    if (rank == 0) print *, 'Enter vCmat. Please ignore timer under this line.'
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0

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
            jt = it + ninact
            Do iv = 1, nact
                jv = iv + ninact
                Do iu = 1, nact
                    ju = iu + ninact

                    !     EatEuv|0>
                    !                    if((it == iv).and.(iu/=iv)) goto 100

                    syma = MULTB_D(irpmo(ju), irpmo(jv))
                    symb = MULTB_D(isym, irpmo(jt))
                    symc = MULTB_S(syma, symb)

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
    do isym = 1, nsymrpa
        if (rank == 0) print *, 'solvC: isym, dim(isym)', isym, dim(isym)
    end do

    !$OMP parallel do schedule(dynamic,1) private(ia,ja,it,jt,cint1)
    Do ia = rank + 1, nsec, nprocs
        ja = ia + ninact + nact
        Do it = 1, nact
            jt = it + ninact

            Call tramo1_ty(ja, jt, cint1)

            effh(ja, jt) = cint1
            !              write(*,'("1int  ",2I4,2E20.10)')ja,jt,effh(ja,jt)

        End do
    End do
    !$OMP end parallel do
    open (1, file=c1int, status='old', form='unformatted')

    do ! Read TYPE 1 integrals C1int until EOF
        read (1, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
        ! Exit loop if the iostat is less than 0  (End of File)
        if (iostat < 0) then
            if (rank == 0) print *, 'End of C1int'
            exit
        else if (iostat > 0) then
            ! Stop the program if the iostat is greater than 0
            stop 'Error: Error in reading C1int'
        end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! + Siguma_pqr<0|EvuEtrEpq|0>(ar|pq)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        isym = irpmo(i)   ! i corresponds to a
        !$OMP parallel do schedule(static,1) private(it,iu,iv,jt,ju,jv,dr,di,d)
        Do i0 = 1, dim(isym)
            it = indt(i0, isym)
            iu = indu(i0, isym)
            iv = indv(i0, isym)
            jt = it + ninact
            ju = iu + ninact
            jv = iv + ninact

            Call dim3_density(iv, iu, it, j - ninact, k - ninact, l - ninact, dr, di)
            d = DCMPLX(dr, di)
            v(i, jt, ju, jv) = v(i, jt, ju, jv) + cint2*d

        End do
        !$OMP end parallel do
!  effh(a,p) =  hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)] - Siguma_w(aw|wp)
!                                                           ~~~~~~~~~~~~~~~~~~~
        if (j == k) then
            effh(i, l) = effh(i, l) - cint2
        end if
    end do

    close (1)
    if (rank == 0) print *, 'reading C1int2 is over'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    open (1, file=c2int, status='old', form='unformatted') ! TYPE 2 integrals

    do ! Read TYPE 2 integrals C2int until EOF
        read (1, iostat=iostat) i, j, k, l, cint2
        ! Exit loop if the iostat is less than 0  (End of File)
        if (iostat < 0) then
            if (rank == 0) then
                print *, 'End of C2int'
            end if
            exit
        else if (iostat > 0) then
            ! Stop the program if the iostat is greater than 0
            stop 'Error: Error in reading C2int'
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

    close (1)
    if (rank == 0) print *, 'reading C2int2 is over'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0

    open (1, file=c3int, status='old', form='unformatted') ! TYPE 3 integrals

    do ! Read TYPE 3 integrals C3int until EOF
        read (1, iostat=iostat) i, j, k, l, cint2 !  (ij|kl):=> (ak|kp)
        ! Exit loop if the iostat is less than 0 (End of File)
        if (iostat < 0) then
            if (rank == 0) then
                print *, 'End of C3int'
            end if
            exit
        else if (iostat > 0) then
            ! Stop the program if the iostat is greater than 0
            stop 'Error: Error in reading C3int'
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
    close (1)
    if (rank == 0) print *, 'reading C3int2 is over'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, effh(ninact + nact + 1, ninact + 1), nsec*nact, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

! Siguma_p effh(a,p)<0|EvuEtp|0>
    !$OMP parallel do schedule(dynamic,1) private(ja,isym,jp,i0,it,iu,iv,jt,ju,jv,dr,di,d)
    Do ia = rank + 1, nsec, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        ja = ia + ninact + nact
        isym = irpmo(ja)

        Do ip = 1, nact
            jp = ip + ninact

            ! Go to the next ip if the value of effh(ja,jp) is nearly zero
            if (ABS(effh(ja, jp)) < 1.0d-10) cycle

            Do i0 = 1, dim(isym)
                it = indt(i0, isym)
                iu = indu(i0, isym)
                iv = indv(i0, isym)
                jt = it + ninact
                ju = iu + ninact
                jv = iv + ninact

                Call dim2_density(iv, iu, it, ip, dr, di)
                d = DCMPLX(dr, di)

                v(ja, jt, ju, jv) = v(ja, jt, ju, jv) + effh(ja, jp)*d

            End do            !i0

        End do               !ip
    End do                  !ia
    !$OMP end parallel do

    if (rank == 0) print *, 'vCmat_ord is ended'

    deallocate (indt)
    deallocate (indu)
    deallocate (indv)

#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, v(ninact + nact + 1, ninact + 1, ninact + 1, ninact + 1), nsec*nact**3, &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)

    if (rank == 0) print *, 'end Allreduce vCmat'
#endif
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
end subroutine vCmat_ord_ty
