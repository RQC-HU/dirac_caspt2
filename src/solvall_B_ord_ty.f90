! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE solvB_ord_ty(e0, e2b)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    real*8, intent(in) :: e0
    real*8, intent(out):: e2b

    integer :: dimn, dimm, dammy

    integer, allocatable :: indsym(:, :)

    real*8, allocatable  :: wsnew(:), ws(:), wb(:)
    real*8               :: e2(2*nsymrpa), e, alpha

    complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
    complex*16, allocatable  :: bc(:, :)
    complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :), vc(:), vc1(:)

    integer, allocatable     :: ii0(:), ij0(:), iij(:, :)
    integer                  :: nij, count

    logical :: cutoff
    integer :: j, i, syma, isym, i0
    integer :: ij, it, ii, iu, jj, jt, ji, ju

    real*8  :: thresd
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE B IS NOW CALCULATED
!
!     EtiEuj|0>
!
!   DRAS1 = -2   DRAS2 = 2   DRAS3 = 0
!
!   t > u, i > j
!
!
!  S(xkyl,tiuj) = d(ki)d(lj)S(xy,tu)
!
!  S(xy, tu) = <0|EtxEuy|0> -d(tx)<0|Euy|0> -d(uy)<0|Etx|0> -d(ty)<0|Eux|0> +d(tx)d(uy)-d(ty)d(ux)
!
!  B(xy,tu) = Siguma_w [eps(w){<0|EtxEuyEww|0>-d(tx)<0|EuyEww|0> -d(uy)<0|EtxEww|0> -d(ty)<0|EuxEww|0>]
!
!              +{d(tx)d(uy)-d(ty)d(ux)}*e0      +S(xy,tu){eps(t)+eps(u)}
!
!  a(i,j)       = -eps(i) - eps(j) - e0
!
! V(i,j) =  SIGUMA_p,q:active <0|EptEqu|0>(pi|qj)
!
!           + SIGUMA_p:active[<0|Ept|0> {(pj|ui) - (uj|pi)}  - <0|Epu|0> (ti|pj)]
!
!           + (uj|ti)  - (tj|ui)
!
!
!  E2 = SIGUMA_a,i, dimm |V1(dimm,ai)|^2|/{(alpha(ai) + wb(dimm)}

    thresd = 1.0D-08
    thres = 1.0D-08

    e2 = 0.0d+00
    e2b = 0.0d+00
    dimn = 0
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0
    if (rank == 0) then
        print *, ' ENTER solv B part'
        print *, ' nsymrpa', nsymrpa
    end if
    Allocate (iij(ninact, ninact)); Call memplus(KIND(iij), SIZE(iij), 1)
    iij = 0
    ! (ninact*(ninact-1))/2 means the number of (ii,ij) pairs (ii>ij)
    nij = (ninact*(ninact - 1))/2
    Allocate (ii0(nij)); Call memplus(KIND(ii0), SIZE(ii0), 1)
    Allocate (ij0(nij)); Call memplus(KIND(ii0), SIZE(ii0), 1)

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
    Allocate (v(nij, nact, nact))
    Call memplus(KIND(v), SIZE(v), 2)
    v = 0.0d+00
    if (rank == 0) print *, 'end before v matrices'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Call vBmat_ord_ty(nij, iij, v)
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0

!     EtiEuj|0>

    Do isym = 1, nsymrpa

        dimn = 0
        Do it = 1, nact
            jt = it + ninact
            Do iu = 1, it - 1
                ju = iu + ninact

                if (nsymrpa /= 1) syma = MULTB_D(irpmo(jt), irpmo(ju) - (-1)**(mod(irpmo(ju), 2)))

                if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                    dimn = dimn + 1
                End if
            End do               ! iu
        End do                  ! it

        if (rank == 0) print *, 'isym, dimn', isym, dimn
        If (dimn == 0) cycle ! Go to the next isym.

        Allocate (indsym(2, dimn)); Call memplus(KIND(indsym), SIZE(indsym), 1)

        dimn = 0
        Do it = 1, nact
            jt = it + ninact
            Do iu = 1, it - 1
                ju = iu + ninact

                if (nsymrpa /= 1) syma = MULTB_D(irpmo(jt), irpmo(ju) - (-1)**(mod(irpmo(ju), 2)))

                if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                    dimn = dimn + 1
                    indsym(1, dimn) = it
                    indsym(2, dimn) = iu
                End if
            End do               ! iu
        End do                  ! it

        Allocate (sc(dimn, dimn)); Call memplus(KIND(sc), SIZE(sc), 2)
        sc = 0.0d+00            ! sc N*N
        if (rank == 0) print *, 'before sBmat'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call sBmat(dimn, indsym, sc)

!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rank == 0) print *, 'sc matrix is obtained normally'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Allocate (ws(dimn)); Call memplus(KIND(ws), SIZE(ws), 1)

        cutoff = .TRUE.
        thresd = 1.0d-08

        Allocate (sc0(dimn, dimn)); Call memplus(KIND(sc0), SIZE(sc0), 2)
        sc0 = sc
        if (rank == 0) print *, 'before cdiag'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call cdiag(sc, dimn, dimm, ws, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rank == 0) print *, 'after s cdiag, new dimension is', dimm
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        If (dimm == 0) then
            deallocate (indsym); Call memminus(KIND(indsym), SIZE(indsym), 1)
            deallocate (sc0); Call memminus(KIND(sc0), SIZE(sc0), 2)
            deallocate (sc); Call memminus(KIND(sc), SIZE(sc), 2)
            deallocate (ws); Call memminus(KIND(ws), SIZE(ws), 1)
            cycle ! Go to the next isym.
        End if

        Allocate (bc(dimn, dimn)); Call memplus(KIND(bc), SIZE(bc), 2)   ! br N*N
        bc = 0.0d+00
        if (rank == 0) print *, 'before bBmat'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call bBmat(e0, dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (rank == 0) print *, 'bc matrix is obtained normally'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        If (debug) then

            if (rank == 0) print *, 'Check whether U*SU is diagonal'
            Call checkdgc(dimn, sc0, sc, ws)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'Check whether U*SU is diagonal END'
        End if

        deallocate (sc0); Call memminus(KIND(sc0), SIZE(sc0), 2)

        if (rank == 0) print *, 'OK cdiag', dimn, dimm
        Allocate (uc(dimn, dimm)); Call memplus(KIND(uc), SIZE(uc), 2)           ! uc N*M
        Allocate (wsnew(dimm)); Call memplus(KIND(wsnew), SIZE(wsnew), 1)     ! wnew M
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
        deallocate (sc); Call memminus(KIND(sc), SIZE(sc), 2)
        deallocate (ws); Call memminus(KIND(ws), SIZE(ws), 1)
        if (rank == 0) print *, 'before ucramda_s_half'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call ucramda_s_half(uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        deallocate (wsnew); Call memminus(KIND(wsnew), SIZE(wsnew), 1)
        if (rank == 0) print *, 'ucrams half OK'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Allocate (bc0(dimm, dimn)); Call memplus(KIND(bc0), SIZE(bc0), 2) ! bc0 M*N
        bc0 = 0.0d+00
        bc0 = MATMUL(TRANSPOSE(DCONJG(uc)), bc)

        Allocate (bc1(dimm, dimm)); Call memplus(KIND(bc1), SIZE(bc1), 2) ! bc1 M*M
        bc1 = 0.0d+00
        bc1 = MATMUL(bc0, uc)

        If (debug) then

            if (rank == 0) print *, 'Check whether bc1 is hermite or not'
            Do i = 1, dimm
                Do j = i, dimm
                    if (ABS(bc1(i, j) - DCONJG(bc1(j, i))) > 1.0d-6) then
                        if (rank == 0) print '(2I4,2E15.7)', i, j, bc1(i, j) - bc1(j, i)
                    End if
                End do
            End do
            if (rank == 0) print *, 'Check whether bc1 is hermite or not END'
        End if

        deallocate (bc); Call memminus(KIND(bc), SIZE(bc), 2)
        deallocate (bc0); Call memminus(KIND(bc0), SIZE(bc0), 2)

        cutoff = .FALSE.

        Allocate (wb(dimm)); Call memplus(KIND(wb), SIZE(wb), 1)

        if (rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'
        Allocate (bc0(dimm, dimm)); Call memplus(KIND(bc0), SIZE(bc0), 2) ! bc0 M*M
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

        deallocate (bc0); Call memminus(KIND(bc0), SIZE(bc0), 2)

        if (rank == 0) print *, 'bC1 matrix is diagonalized!'
        e2 = 0.0d+00

        Do i0 = 1, nij
            ji = ii0(i0)
            jj = ij0(i0)
            if (nsymrpa /= 1) then
                syma = MULTB_D(irpmo(ji) - (-1)**(mod(irpmo(ji), 2)), irpmo(jj))
                syma = MULTB_S(syma, isym)
            end if
            If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then

                Allocate (vc(dimn)); Call memplus(KIND(vc), SIZE(vc), 2)
                Do it = 1, dimn
                    vc(it) = v(i0, indsym(1, it), indsym(2, it))
                End do

                Allocate (vc1(dimm)); Call memplus(KIND(vc1), SIZE(vc1), 2)
                vc1 = 0.0d+00

                vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn, 1:dimm))), vc(1:dimn)) ! v => v~
                Deallocate (vc); Call memminus(KIND(vc), SIZE(vc), 2)

                alpha = -eps(ji) - eps(jj) - e0 + eshift   ! For Level Shift (2007/2/9)

                vc1(1:dimm) = &
                & MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm, 1:dimm))), vc1(1:dimm)) ! v~ => v~~

                Do j = 1, dimm
                    sumc2local = sumc2local + (ABS(vc1(j))**2.0d+00)/((alpha + wb(j))**2.0d+00)
                    e = (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                    e2(isym) = e2(isym) - e
                End do

                Deallocate (vc1); Call memminus(KIND(vc1), SIZE(vc1), 2)

            End if

        End do            !i0
        if (rank == 0) print '("e2b(",I3,") = ",E20.10,"a.u.")', isym, e2(isym)
        Deallocate (bc1); Call memminus(KIND(bc1), SIZE(bc1), 2)
        Deallocate (uc); Call memminus(KIND(uc), SIZE(uc), 2)
        Deallocate (wb); Call memminus(KIND(wb), SIZE(wb), 1)
        Deallocate (indsym); Call memminus(KIND(indsym), SIZE(indsym), 2)

        e2b = e2b + e2(isym)
        if (rank == 0) print *, 'End e2(isym) add'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0

    End do               ! isym

    if (rank == 0) then
        print '("e2b      = ",E20.10," a.u.")', e2b
        print '("sumc2,b  = ",E20.10)', sumc2local
    end if
    sumc2 = sumc2 + sumc2local

    deallocate (iij); Call memminus(KIND(iij), SIZE(iij), 1)
    deallocate (ii0); Call memminus(KIND(ii0), SIZE(ii0), 1)
    deallocate (ij0); Call memminus(KIND(ij0), SIZE(ij0), 1)
    deallocate (v); Call memminus(KIND(v), SIZE(v), 2)

    continue
    if (rank == 0) print *, 'end solvB_ord_ty'
end

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE sBmat(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space B

!  S(xy, tu) = <0|EtxEuy|0> -d(tx)<0|Euy|0> -d(uy)<0|Etx|0> -d(ty)<0|Eux|0> +d(tx)d(uy)-d(ty)d(ux)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer, intent(in)      :: dimn, indsym(2, dimn)
    complex*16, intent(out)  :: sc(dimn, dimn)

    real*8  :: a, b

    integer :: it, iu, iy, ix
    integer :: i, j

    sc = 0.0d+00

    !$OMP parallel do schedule(dynamic,1) private(ix,iy,j,it,iu,a,b)
    Do i = rank + 1, dimn, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)

        ix = indsym(1, i)
        iy = indsym(2, i)

        Do j = i, dimn

            it = indsym(1, j)
            iu = indsym(2, j)

!  S(xy, tu) = <0|EtxEuy|0> -d(tx)<0|Euy|0> -d(uy)<0|Etx|0> -d(ty)<0|Eux|0> +d(tx)d(uy)-d(ty)d(ux)
!                                                                                      ~~~~~~~~This term is0
            Call dim2_density(it, ix, iu, iy, a, b)
            sc(i, j) = sc(i, j) + DCMPLX(a, b)

            If (it == ix) then
                Call dim1_density(iu, iy, a, b)
                sc(i, j) = sc(i, j) - DCMPLX(a, b)
            End if

            If (iu == iy) then
                Call dim1_density(it, ix, a, b)
                sc(i, j) = sc(i, j) - DCMPLX(a, b)
            End if

            If (it == iy) then
                Call dim1_density(iu, ix, a, b)
                sc(i, j) = sc(i, j) - DCMPLX(a, b)
            End if

            If ((it == ix) .and. (iu == iy)) then
                sc(i, j) = sc(i, j) + 1.0d+00
            End if

            sc(j, i) = DCONJG(sc(i, j))

        End do               !j
    End do                  !i
    !$OMP end parallel do
#ifdef HAVE_MPI
    call allreduce_wrapper(mat=sc)
#endif

End subroutine sBmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE bBmat(e0, dimn, sc, indsym, bc) ! Assume C1 molecule, overlap matrix B in space B
!
!  B(xy,tu) = Siguma_w [eps(w){<0|EtxEuyEww|0>-d(tx)<0|EuyEww|0> -d(uy)<0|EtxEww|0> -d(ty)<0|EuxEww|0>]
!
!              +{d(tx)d(uy)-d(ty)d(ux)}*e0      +S(xy,tu){eps(t)+eps(u)}
!
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    integer, intent(in) :: dimn, indsym(2, dimn)
    complex*16, intent(in)  :: sc(dimn, dimn)
    complex*16, intent(out) :: bc(dimn, dimn)
    real*8, intent(in)      :: e0

    real*8              :: e, denr, deni
    complex*16          :: den

    integer :: it, iu, ix, iy, iw
    integer :: jt, ju, jy, jx, jw, i, j

    bc(:, :) = 0.0d+00

    if (rank == 0) print *, 'B space Bmat iroot=', iroot

    !$OMP parallel do schedule(dynamic,1) private(ix,iy,jx,jy,it,iu,jt,ju,e,j,iw,jw,denr,deni,den)
    Do i = rank + 1, dimn, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)

        ix = indsym(1, i)
        jx = ix + ninact
        iy = indsym(2, i)
        jy = iy + ninact

        Do j = i, dimn

            it = indsym(1, j)
            jt = it + ninact
            iu = indsym(2, j)
            ju = iu + ninact

!  B(xy,tu) = Siguma_w [eps(w){<0|EtxEuyEww|0>-d(tx)<0|EuyEww|0> -d(uy)<0|EtxEww|0> -d(ty)<0|EuxEww|0>]
!
!              +{d(tx)d(uy)-d(ty)d(ux)}*e0      +S(xy,tu){eps(t)+eps(u)}

            e = eps(jt) + eps(ju)

            Do iw = 1, nact
                jw = iw + ninact

                Call dim3_density(it, ix, iu, iy, iw, iw, denr, deni)
                den = DCMPLX(denr, deni)
                bc(i, j) = bc(i, j) + den*eps(jw)

                If (it == ix) then

                    Call dim2_density(iu, iy, iw, iw, denr, deni)
                    den = DCMPLX(denr, deni)
                    bc(i, j) = bc(i, j) - den*eps(jw)

                End if

                If (iu == iy) then

                    Call dim2_density(it, ix, iw, iw, denr, deni)
                    den = DCMPLX(denr, deni)
                    bc(i, j) = bc(i, j) - den*eps(jw)

                End if

                If (it == iy) then

                    Call dim2_density(iu, ix, iw, iw, denr, deni)
                    den = DCMPLX(denr, deni)
                    bc(i, j) = bc(i, j) - den*eps(jw)

                End if

            End do

!              +{d(tx)d(uy)-d(ty)d(ux)}*e0      +S(xy,tu){eps(t)+eps(u)}

            If ((it == ix) .and. (iu == iy)) then
                bc(i, j) = bc(i, j) + e0
            End if

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
    if (rank == 0) print *, 'bBmat is ended'
End subroutine bBmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE vBmat_ord_ty(nij, iij, v)
!
!
! V(i,j) =  SIGUMA_p,q:active <0|EptEqu|0>(pi|qj)
!
!           + SIGUMA_p:active[<0|Ept|0> {(pj|ui) - (uj|pi)}  - <0|Epu|0> (ti|pj)]
!
!           + (uj|ti)  - (tj|ui)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
    use module_file_manager
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer, intent(in)     :: nij, iij(ninact, ninact)
    complex*16, intent(out) :: v(nij, nact, nact)
    real*8                  :: dr, di
    complex*16              :: cint2, dens
    integer :: i, j, k, l, tij, i0
    integer :: it, iu, iostat, twoint_unit
    integer :: isym, syma, jt, ju
    integer :: multb_s_reverse(ninact, ninact), pattern_t(nact**2, nsymrpa), pattern_u(nact**2, nsymrpa), pattern_tu_count(nsymrpa)
    logical :: is_end_of_file

    v = 0.0d+00
    twoint_unit = default_unit
    multb_s_reverse(:, :) = 0
    call create_multb_s_reverse

    ! Save t,u patterns for each isym
    pattern_t(:, :) = 0
    pattern_u(:, :) = 0
    pattern_tu_count(:) = 0
    do isym = 1, nsymrpa
        Do it = 1, nact
            jt = it + ninact
            Do iu = 1, it - 1
                ju = iu + ninact

                if (nsymrpa /= 1) syma = MULTB_D(irpmo(jt), irpmo(ju) - (-1)**(mod(irpmo(ju), 2)))

                if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                    pattern_tu_count(isym) = pattern_tu_count(isym) + 1
                    pattern_t(pattern_tu_count(isym), isym) = it
                    pattern_u(pattern_tu_count(isym), isym) = iu
                End if
            End do
        End do
    end do
    call open_unformatted_file(unit=twoint_unit, file=bint, status='old', optional_action='read') !  (21|21) stored (ti|uj) i > j
    do
        read (twoint_unit, iostat=iostat) i, j, k, l, cint2                    !  (ij|kl)
        call check_iostat(iostat=iostat, file=bint, end_of_file_reached=is_end_of_file)
        if (is_end_of_file) then
            exit
        end if

        if (j <= l) cycle ! Read the next line if j <= l

!------------------------------------------------------------------------------------------------
!  i > j
!
! V(i,j) =  SIGUMA_p,q:active <0|EptEqu|0>(pi|qj)                                      ! term1
!
!           + SIGUMA_p:active[<0|Ept|0> {(ui|pj) - (pi|uj)}  - <0|Epu|0> (ti|pj)]      ! term2
!
!           + (ti|uj)  - (ui|tj)                                                       ! term3
!
!------------------------------------------------------------------------------------------------

        tij = iij(j, l)

!        write(*,'(5I4,2E20.10)')i,j,k,l,tij,cint2

        ! Term 3 !        + (ti|uj)  - (ui|tj)  (i > j)

        v(tij, i, k) = v(tij, i, k) + cint2 !  + (ti|uj)
        v(tij, k, i) = v(tij, k, i) - cint2 !  - (ui|tj)

        ! Term 2 !  + SIGUMA_p:active[<0|Ept|0> {(ui|pj) - (pi|uj)}  - <0|Epu|0> (ti|pj)]
        !                             ===========================      ================
        !                                loop for t                     loop for u(variable u is renamed to t)
        !!$OMP parallel
        !!$OMP do schedule(dynamic,1) private(dr,di,dens,iu)
        Do it = 1, nact

            Call dim1_density(k, it, dr, di)
            dens = DCMPLX(dr, di)
            v(tij, it, i) = v(tij, it, i) + cint2*dens
            v(tij, i, it) = v(tij, i, it) - cint2*dens

            Call dim1_density(i, it, dr, di)
            dens = DCMPLX(dr, di)
            v(tij, it, k) = v(tij, it, k) - cint2*dens
        end do
        !!$OMP end do
        isym = multb_s_reverse(j, l)
        !!$OMP do schedule(dynamic,1) private(it,iu,dr,di,dens)
        do i0 = 1, pattern_tu_count(isym)
            it = pattern_t(i0, isym)
            iu = pattern_u(i0, isym)

            ! Term1 !   SIGUMA_p,q:active <0|EptEqu|0>(pi|qj)                                      ! term1
            !                             ==================
            !                              loop for t and u

            ! Do iu = 1, it - 1
            Call dim2_density(i, it, k, iu, dr, di)
            dens = DCMPLX(dr, di)
            v(tij, it, iu) = v(tij, it, iu) + cint2*dens
            ! End do

        End do
        !!$OMP end do
        !!$OMP end parallel
    end do

    close (twoint_unit)
    if (rank == 0) print *, 'vBmat_ord_ty is ended'

#ifdef HAVE_MPI
    call allreduce_wrapper(mat=v)
#endif
contains
    subroutine create_multb_s_reverse
        !========================================================================================================
        ! This subroutine creates multb_s_reverse
        !
        ! multb_s_reverse(i, j) returns the symmetry of MULTB_D(irpmo(jt), irpmo(ju) - (-1)**(mod(irpmo(ju), 2)))
        !========================================================================================================
        implicit none
        integer :: ii, ij

        if (nsymrpa == 1) then
            multb_s_reverse(:, :) = 1
        else
            do ii = 1, ninact
                do ij = 1, ii - 1
                    syma = MULTB_D(irpmo(ii) - (-1)**(mod(irpmo(ii), 2)), irpmo(ij))
                    do isym = 1, nsymrpa
                        if (MULTB_S(syma, isym) == 1) then
                            multb_s_reverse(ii, ij) = isym
                            exit
                        end if
                    end do
                end do
            end do
        end if
    end subroutine create_multb_s_reverse
end subroutine vBmat_ord_ty
