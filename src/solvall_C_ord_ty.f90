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

    integer :: dimn, dimm, count, dammy

    integer, allocatable :: indsym(:, :)

    real*8, allocatable  :: sr(:, :), ur(:, :)
    real*8, allocatable  :: br(:, :), wsnew(:), ws(:), wb(:)
    real*8, allocatable  :: br0(:, :), br1(:, :)
    real*8               :: e2(nsymrp), dr, di, alpha

    complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
    complex*16, allocatable  :: bc(:, :)
    complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :, :), vc(:), vc1(:)

    logical :: cutoff
    integer :: j, i, k, syma, symb, isym, sym1, i0
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

!        thresd = thres
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
    if (rank == 0) then ! Process limits for output
        write (*, *) ' ENTER solv C part'
        write (*, *) ' nsymrpa', nsymrpa
    end if
    Allocate (v(ninact + nact + 1:ninact + nact + nsec, ninact + 1:ninact + nact,  &
    &          ninact + 1:ninact + nact, ninact + 1:ninact + nact))
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (*, *) 'end before v matrices'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Call vCmat_ord_ty(v)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'come'
    end if
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
                    !             Do iz = 1, ix-1
                    !                if((ix == iz).and.(iy/=iz)) goto 100

                    jx = ix + ninact
                    jy = iy + ninact
                    jz = iz + ninact
                    syma = MULTB_D(isym, irpmo(jx))
                    symb = MULTB_D(irpmo(jy), irpmo(jz))
                    syma = MULTB_S(syma, symb)

                    !Iwamuro modify
                    !                write(*,*)"syma1",syma

                    If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then
                        ixyz = ixyz + 1
                    End if
                    !Iwamuro modify
                    !                write(*,*)"ixyz1",ixyz

100             End do
            End do
        End do

        dimn = ixyz

        If (dimn == 0) goto 1000

        Allocate (indsym(3, dimn))
        indsym = 0
        ixyz = 0

        Do ix = 1, nact
            Do iy = 1, nact
                Do iz = 1, nact

                    !                if((ix == iz).and.(iy/=iz)) goto 200

                    jx = ix + ninact
                    jy = iy + ninact
                    jz = iz + ninact
                    syma = MULTB_D(isym, irpmo(jx))
                    symb = MULTB_D(irpmo(jy), irpmo(jz))
                    syma = MULTB_S(syma, symb)

                    !Iwamuro modify
                    !                write(*,*)"syma2",syma

                    If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then

                        ixyz = ixyz + 1
                        indsym(1, ixyz) = ix
                        indsym(2, ixyz) = iy
                        indsym(3, ixyz) = iz
                    End if

200             End do
            End do
        End do

        if (rank == 0) then ! Process limits for output
            write (*, *) 'isym, dimn', isym, dimn
        end if
        Allocate (sc(dimn, dimn))
        sc = 0.0d+00            ! sr N*N
        if (rank == 0) then ! Process limits for output
            write (*, *) 'before sCmat'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call sCmat(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (rank == 0) then ! Process limits for output
            write (*, *) 'sC matrix is obtained normally'
        end if
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
        if (rank == 0) then ! Process limits for output
            write (*, *) 'before cdiag'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call cdiag(sc, dimn, dimm, ws, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rank == 0) then ! Process limits for output
            write (*, *) 'after sc cdiag'
            write (*, *) 'after s cdiag, new dimension is', dimm
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        If (dimm == 0) then
            deallocate (indsym)
            deallocate (sc0)
            deallocate (sc)
            deallocate (ws)
            goto 1000
        End if

        If (debug) then

            if (rank == 0) then ! Process limits for output
                write (*, *) 'Check whether U*SU is diagonal'
            end if

            Call checkdgc(dimn, sc0, sc, ws)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (rank == 0) then ! Process limits for output
                write (*, *) 'Check whether U*SU is diagonal END'
            end if
        End if
        if (rank == 0) then ! Process limits for output
            write (*, *) 'OK cdiag', dimn, dimm
        end if
        Allocate (bc(dimn, dimn))                                 ! br N*N
        bc = 0.0d+00
        if (rank == 0) then ! Process limits for output
            write (*, *) 'before bCmat'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call bCmat(dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        deallocate (sc0)

        if (rank == 0) then ! Process limits for output
            write (*, *) 'bC matrix is obtained normally'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Allocate (uc(dimn, dimm))                                 ! uc N*M
        Allocate (wsnew(dimm))                                    ! wnew M
        uc(:, :) = 0.0d+00
        wsnew(:) = 0.0d+00
        if (rank == 0) then ! Process limits for output
            write (*, *) 'before ccutoff'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call ccutoff(sc, ws, dimn, dimm, uc, wsnew)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!           Do i0 = 1, dimm
!           write(*,'(E20.10)') wsnew(i0)
!           End do

        if (rank == 0) then ! Process limits for output
            write (*, *) 'OK ccutoff'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        deallocate (ws)
        deallocate (sc)
        if (rank == 0) then ! Process limits for output
            write (*, *) 'before ucramda_s_half'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call ucramda_s_half(uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        deallocate (wsnew)

        if (rank == 0) then ! Process limits for output
            write (*, *) 'ucrams half OK'
        end if
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
            if (rank == 0) then ! Process limits for output
                write (*, *) 'Check whether bc1 is hermite or not'
                Do i = 1, dimm
                    Do j = i, dimm
                        if (ABS(bc1(i, j) - DCONJG(bc1(j, i))) > 1.0d-6) then
                            write (*, '(2I4,2E15.7)') i, j, bc1(i, j) - bc1(j, i)
                        End if
                    End do
                End do
                write (*, *) 'Check whether bc1 is hermite or not END'
            end if
        End if

        deallocate (bc)
        deallocate (bc0)

        cutoff = .FALSE.

        Allocate (wb(dimm))
        wb = 0.0d+00
        if (rank == 0) then ! Process limits for output
            write (*, *) 'bC matrix is transrated to bc1(M*M matrix)!'
        end if
        Allocate (bc0(dimm, dimm))
        bc0 = 0.0d+00
        bc0 = bc1
        if (rank == 0) then ! Process limits for output
            write (*, *) 'before cdiag'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call cdiag(bc1, dimm, dammy, wb, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rank == 0) then ! Process limits for output
            write (*, *) 'end cdiag'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        If (debug) then

            if (rank == 0) then ! Process limits for output
                write (*, *) 'Check whether bc is really diagonalized or not'
            end if
            Call checkdgc(dimm, bc0, bc1, wb)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) then ! Process limits for output
                write (*, *) 'Check whether bc is really diagonalized or not END'
            end if
        End if

        deallocate (bc0)

        if (rank == 0) then ! Process limits for output
            write (*, *) 'bC1 matrix is diagonalized!'
        end if
        Do ia = 1, nsec
            ja = ia + ninact + nact

            if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. irpmo(ja) == isym)) then

                Allocate (vc(dimn))
                Do it = 1, dimn
                    vc(it) = v(ja, indsym(1, it) + ninact, indsym(2, it) + ninact, indsym(3, it) + ninact)
!                    write(*,'(4I4,2E20.10)') &
!                    & ja,indsym(1,it)+ninact,indsym(2,it)+ninact,indsym(3,it)+ninact,vc(it)
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

        if (rank == 0) then ! Process limits for output
            write (*, '("e2c(",I3,") = ",E20.10,"a.u.")') isym, e2(isym)
        end if
        deallocate (bc1)
        deallocate (indsym)
        Deallocate (uc)
        Deallocate (wb)

        e2c = e2c + e2(isym)
        if (rank == 0) then ! Process limits for output
            write (*, *) 'End e2(isym) add'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
1000 End do                  ! isym

    if (rank == 0) then ! Process limits for output
        write (*, '("e2c      = ",E20.10,"a.u.")') e2c
        write (*, '("sumc2,c  = ",E20.10)') sumc2local
    end if
    sumc2 = sumc2 + sumc2local

    continue
    if (rank == 0) then ! Process limits for output
        write (*, *) 'end solvC_ord_ty'
    end if
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
    integer :: jt, ju, jv, jx, jy, jz
    integer :: i, j
    integer :: count

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
            if (rank == 0) then ! Process limits for output
                If (ABS(sc(i, j)) > 1.0d+00) then
                    write (*, '(2I4,2E20.10)') i, j, sc(i, j)
                End if
            end if
        End do               !j
    End do                  !i
    !$OMP end parallel do
#ifdef HAVE_MPI
    call MPI_Allreduce (MPI_IN_PLACE, sc(1,1), dimn**2, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
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

    if (rank == 0) then ! Process limits for output
        write (*, *) 'C space Bmat iroot=', iroot
    end if
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
    if (rank == 0) then ! Process limits for output
        write (*, *) 'bCmat is ended'
    end if
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

    real*8                  :: dr, di, signkl
    complex*16              :: cint1, cint2, term1, d
    complex*16              :: effh(ninact + nact + 1:ninact + nact + nsec, ninact + 1:ninact + nact)

    integer :: i, j, k, l, kkr, lkr, count, dim(nsymrpa)
    integer :: isym, sym, syma, symb, symc

    integer, allocatable :: indt(:, :), indu(:, :), indv(:, :)
    logical :: test
    integer :: it, iu, iv, iw, ia, ip, iq, ir, ik
    integer :: jt, ju, jv, jw, ja, jp, jq, jr, jk
    integer :: i0
    integer :: loopcnt
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

    if (rank == 0) write(*,*) 'Enter vCmat. Please ignore timer under this line.'
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
                        !Iwamuro modify
                        ! write (*, *) it, iu, iv
                    end if
100             End do
            End do
        End do
    End do
    !$OMP end parallel do
    do isym = 1, nsymrpa
        if (rank == 0) write(*,*) 'solvC: isym, dim(isym)', isym, dim(isym)
    end do

!        Do isym = 1, nsymrpa
!           write(*,*)dim(isym),isym
!        Enddo
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
!  open (1, file='C1int', status='old', form='unformatted')
    ! open (1, file=c1int, status='old', form='formatted')
    open (1, file=c1int, status='old', form='unformatted')

! 30  read (1, '(4I4, 2e20.10)', err=10, end=20) i, j, k, l, cint2 !  (ij|kl)
30  read (1, err=10, end=20) i, j, k, l, cint2 !  (ij|kl)

!           write(*,'("TYPE 1  ",4I4,2E20.10)')i,j,k,l,cint2

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

!                    test =.FALSE.
!                    if(j==1.and.jt==3.and.ju==4.and.jv==8) test=.TRUE.

        Call dim3_density(iv, iu, it, j - ninact, k - ninact, l - ninact, dr, di)
        d = DCMPLX(dr, di)
        v(i, jt, ju, jv) = v(i, jt, ju, jv) + cint2*d

!                    if(test.and.ABS(d)>1.0d-10.and.AbS(cint2)>1.0d-10) &
!                    & write(*,'("3dim-2int2",6I4,2E20.10,4I4,2E20.10)') &
!                    &iv, iu,  i-ninact, it, k-ninact, l-ninact, d, i,j,k,l,cint2

    End do
    !$OMP end parallel do
!  effh(a,p) =  hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)] - Siguma_w(aw|wp)
!                                                           ~~~~~~~~~~~~~~~~~~~
    if (j == k) then
        effh(i, l) = effh(i, l) - cint2
    end if

    goto 30

20  close (1)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'reading C1int2 is over'
    end if
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
!  open (1, file='C2int', status='old', form='unformatted') ! TYPE 2 integrals
    ! open (1, file=c2int, status='old', form='formatted') ! TYPE 2 integrals
    open (1, file=c2int, status='old', form='unformatted') ! TYPE 2 integrals

! 300 read (1, '(4I4, 2e20.10)', err=10, end=200) i, j, k, l, cint2
300 read (1, err=10, end=200) i, j, k, l, cint2

!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  effh(a,p) =  hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)] - Siguma_w(aw|wp)
!                                          ========
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (k == l) then

        effh(i, j) = effh(i, j) + cint2

    end if

    goto 300

200 close (1)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'reading C2int2 is over'
    end if
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0

!  open (1, file='C3int', status='old', form='unformatted') ! TYPE 3 integrals
    ! open (1, file=c3int, status='old', form='formatted') ! TYPE 3 integrals
    open (1, file=c3int, status='old', form='unformatted') ! TYPE 3 integrals

! 3000 read (1, '(4I4, 2e20.10)', err=10, end=2000) i, j, k, l, cint2 !  (ij|kl):=> (ak|kp)
3000 read (1, err=10, end=2000) i, j, k, l, cint2 !  (ij|kl):=> (ak|kp)

!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  effh(a,p) =  hap + Siguma_k(is oqqupied)[(ap|kk)-(ak|kp)] - Siguma_w(aw|wp)
!                                                  =========
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (j == k) then

        effh(i, l) = effh(i, l) - cint2

    end if

    goto 3000

2000 close (1)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'reading C3int2 is over'
    end if
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    !    effh(ninact + nact + 1:ninact + nact + nsec, ninact + 1:ninact + nact)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, effh(ninact + nact + 1, ninact + 1), nsec*nact, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

! Siguma_p effh(a,p)<0|EvuEtp|0>
    !$OMP parallel do schedule(dynamic,1) private(ja,isym,jp,i0,it,iu,iv,jt,ju,jv,dr,di,d)
    Do ia = rank + 1, nsec, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        ja = ia + ninact + nact
!              write(*,'("effh  ",2I4,2E20.10)')ja, jp, effh(ja,jp)

        isym = irpmo(ja)

        Do ip = 1, nact
            jp = ip + ninact

!              write(*,'("effh  ",2I4,2E20.10)')ja,jp,effh(ja,jp)
            if (ABS(effh(ja, jp)) < 1.0d-10) goto 70

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

70      End do               !ip
    End do                  !ia
    !$OMP end parallel do

    goto 101

10  write (*, *) 'error while opening file Cint'; goto 101

101 if (rank == 0) write (*, *) 'vCmat_ord is ended'

    deallocate (indt)
    deallocate (indu)
    deallocate (indv)

!    Allocate (v(ninact + nact + 1:ninact + nact + nsec, ninact + 1:ninact + nact,  &
!    &          ninact + 1:ninact + nact, ninact + 1:ninact + nact))
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, v(ninact + nact + 1, ninact + 1, ninact + 1, ninact + 1), nsec*nact**3, &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write(*,*) 'end Allreduce vCmat'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
end subroutine vCmat_ord_ty