! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE solvA_ord_ty(e0, e2a)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    real*8, intent(in) :: e0
    real*8, intent(out):: e2a

    integer :: dimn, dimm, count, dammy

    integer, allocatable :: indsym(:, :)

    real*8, allocatable  :: sr(:, :), ur(:, :)
    real*8, allocatable  :: br(:, :), wsnew(:), ws(:), wb(:)
    real*8, allocatable  :: br0(:, :), br1(:, :)
    real*8               :: e2(2*nsymrp), alpha

    complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
    complex*16, allocatable  :: bc(:, :)
    complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :, :), vc(:), vc1(:)

    logical :: cutoff
    integer :: i, j, k, syma, symb, isym, i0, j0, sym1
    integer :: ix, iy, iz, ii, dimi, ixyz
    integer :: jx, jy, jz, ji, it

    real*8  :: thresd
    integer :: loopcnt
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1

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

!        thresd = thres
    thresd = 1.0D-08
    thres = 1.0D-08

    e2 = 0.0d+00
    e2a = 0.0d+00
    dimi = 0
    dimn = 0
    syma = 0
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0
    if (rank == 0) then ! Process limits for output
        write (*, *) ' ENTER solv A part'
        write (*, *) ' nsymrpa', nsymrpa
    end if

    Allocate (v(ninact, ninact + 1:ninact + nact, ninact + 1:ninact + nact, ninact + 1:ninact + nact))
    Call memplus(KIND(v), SIZE(v), 2)

    if (rank == 0) then ! Process limits for output
        write (*, *) 'before vAmat'
    end if
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (*, *) 'end before v matrices'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Call vAmat_ord_ty(v)
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
!         ExjEyz
! Noda : Probably need to paralellize it under this
    Do isym = 1, nsymrpa

        ixyz = 0

        Do ix = 1, nact
            Do iy = 1, nact
                Do iz = 1, nact

                    jx = ix + ninact
                    jy = iy + ninact
                    jz = iz + ninact
                    syma = MULTB_D(irpmo(jx), isym)
                    symb = MULTB_D(irpmo(jy), irpmo(jz))
                    syma = MULTB_S(syma, symb)

                    If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then
                        ixyz = ixyz + 1
                    End if

100             End do
            End do
        End do

        dimn = ixyz

        If (dimn == 0) goto 1000

        Allocate (indsym(3, dimn)); Call memplus(KIND(indsym), SIZE(indsym), 1)

        ixyz = 0

        Do ix = 1, nact
            Do iy = 1, nact
                Do iz = 1, nact

                    jx = ix + ninact
                    jy = iy + ninact
                    jz = iz + ninact
                    syma = MULTB_D(irpmo(jx), isym)
                    symb = MULTB_D(irpmo(jy), irpmo(jz))
                    syma = MULTB_S(syma, symb)

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
        Allocate (sc(dimn, dimn)); Call memplus(KIND(sc), SIZE(sc), 2)

        sc = 0.0d+00            ! sr N*N
        if (rank == 0) then ! Process limits for output
            write (*, *) 'before sAmat'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call sAmat(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rank == 0) then ! Process limits for output
            write (*, *) 'sc matrix is obtained normally'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Allocate (ws(dimn)); Call memplus(KIND(ws), SIZE(ws), 1)

        cutoff = .TRUE.
!           thresd = 1.0d-15

        Allocate (sc0(dimn, dimn)); Call memplus(KIND(sc0), SIZE(sc0), 2)
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
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0

        If (dimm == 0) then
            deallocate (indsym); Call memminus(KIND(indsym), SIZE(indsym), 1)
            deallocate (sc0); Call memminus(KIND(sc0), SIZE(sc0), 2)
            deallocate (sc); Call memminus(KIND(sc), SIZE(sc), 2)
            deallocate (ws); Call memminus(KIND(ws), SIZE(ws), 1)
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

        Allocate (bc(dimn, dimn)); Call memplus(KIND(bc), SIZE(bc), 2)   ! br N*N
        bc = 0.0d+00
        if (rank == 0) then ! Process limits for output
            write (*, *) 'before bAmat'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call bAmat(dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rank == 0) then ! Process limits for output
            write (*, *) 'bc matrix is obtained normally'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        deallocate (sc0); Call memminus(KIND(sc0), SIZE(sc0), 2)

        Allocate (uc(dimn, dimm)); Call memplus(KIND(uc), SIZE(uc), 2)           ! uc N*M
        Allocate (wsnew(dimm)); Call memplus(KIND(wsnew), SIZE(wsnew), 1)     ! wnew M
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
        if (rank == 0) then ! Process limits for output
            write (*, *) 'OK ccutoff'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        deallocate (sc); Call memminus(KIND(sc), SIZE(sc), 2)
        deallocate (ws); Call memminus(KIND(ws), SIZE(ws), 1)

        if (rank == 0) then ! Process limits for output
            write (*, *) 'before ucramda_s_half'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call ucramda_s_half(uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        deallocate (wsnew); Call memminus(KIND(wsnew), SIZE(wsnew), 1)

        if (rank == 0) then ! Process limits for output
            write (*, *) 'ucrams half OK'
        end if
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

        deallocate (bc); Call memminus(KIND(bc), SIZE(bc), 2)
        deallocate (bc0); Call memminus(KIND(bc0), SIZE(bc0), 2)

        cutoff = .FALSE.

        Allocate (wb(dimm)); Call memplus(KIND(wb), SIZE(wb), 1)

        if (rank == 0) then ! Process limits for output
            write (*, *) 'bC matrix is transrated to bc1(M*M matrix)!'
        end if

        Allocate (bc0(dimm, dimm)); Call memplus(KIND(bc0), SIZE(bc0), 2) ! bc0 M*M
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

!  If(debug) then

        if (rank == 0) then ! Process limits for output
            write (*, *) 'Check whether bc is really diagonalized or not'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0

        Call checkdgc(dimm, bc0, bc1, wb)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rank == 0) then ! Process limits for output
            write (*, *) 'Check whether bc is really diagonalized or not END'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0

!  End if
        deallocate (bc0); Call memminus(KIND(bc0), SIZE(bc0), 2)

        if (rank == 0) then ! Process limits for output
            write (*, *) 'bC1 matrix is diagonalized!'
        end if

        e2 = 0.0d+00

        Do ii = 1, ninact

            sym1 = irpmo(ii)

            if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. sym1 == isym)) then

                Allocate (vc(dimn)); Call memplus(KIND(vc), SIZE(vc), 2)
                Do it = 1, dimn
                    vc(it) = v(ii, indsym(1, it) + ninact, indsym(2, it) + ninact, indsym(3, it) + ninact)
!                    write(*,'(4I4,2E20.10)') &
!                    & ii,indsym(1,it)+ninact,indsym(2,it)+ninact,indsym(3,it)+ninact,vc(it)
                End do

                Allocate (vc1(dimm)); Call memplus(KIND(vc1), SIZE(vc1), 2)
                vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn, 1:dimm))), vc(1:dimn))

                Deallocate (vc); Call memminus(KIND(vc), SIZE(vc), 2)

                alpha = -eps(ii) - e0 + eshift   ! For Level Shift (2007/2/9)

                vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm, 1:dimm))), vc1(1:dimm))

                Do j = 1, dimm
                    sumc2local = sumc2local + (ABS(vc1(j))**2.0d+00)/((alpha + wb(j))**2.0d+00)
                    e2(isym) = e2(isym) - (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                End do
                Deallocate (vc1); Call memminus(KIND(vc1), SIZE(vc1), 2)

            End if

        End do
        if (rank == 0) then ! Process limits for output
            write (*, '("e2a(",I3,") = ",E20.10,"a.u.")') isym, e2(isym)
        end if

        Deallocate (bc1); Call memminus(KIND(bc1), SIZE(bc1), 2)
        Deallocate (uc); Call memminus(KIND(uc), SIZE(uc), 2)
        Deallocate (wb); Call memminus(KIND(wb), SIZE(wb), 1)
        Deallocate (indsym); Call memminus(KIND(indsym), SIZE(indsym), 2)

        e2a = e2a + e2(isym)
        if (rank == 0) then ! Process limits for output
            write (*, *) 'End e2(isym) add'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
1000 End do                  ! isym

    if (rank == 0) then ! Process limits for output
        write (*, '("e2a      = ",E20.10," a.u.")') e2a

        write (*, '("sumc2,a  = ",E20.10)') sumc2local
    end if
    sumc2 = sumc2 + sumc2local

    Deallocate (v); Call memminus(KIND(v), SIZE(v), 2)

    continue
    if (rank == 0) then ! Process limits for output
        write (*, *) 'end solvA_ord_ty'
    end if
end

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
SUBROUTINE sAmat(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space A
!
!  S(xyz,tuv) =  - <0|EzyEtxEuv|0> + d(tx)<0|EzyEuv|0>
!
!     x > y, t > u

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

    !$OMP parallel do private(i,ix,iy,iz,j,it,iu,iv,a,b)
    Do i = rank + 1, dimn, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
        ix = indsym(1, i)
        iy = indsym(2, i)
        iz = indsym(3, i)

        Do j = i, dimn
            it = indsym(1, j)
            iu = indsym(2, j)
            iv = indsym(3, j)

            a = 0.0d+0
            b = 0.0d+0

!  S(xyz,tuv) =  - <0|EzyEtxEuv|0> + d(tx)<0|EzyEuv|0>

            Call dim3_density(iz, iy, it, ix, iu, iv, a, b)

            sc(i, j) = sc(i, j) - DCMPLX(a, b)

            If (it == ix) then
                a = 0.0d+0
                b = 0.0d+0

                Call dim2_density(iz, iy, iu, iv, a, b)

                sc(i, j) = sc(i, j) + DCMPLX(a, b)

            End if

            sc(j, i) = DCONJG(sc(i, j))

        End do               !j
    End do                  !i
    !$OMP end parallel do
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, sc(1, 1), dimn**2, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

End subroutine sAmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE bAmat(dimn, sc, indsym, bc) ! Assume C1 molecule, overlap matrix B in space A
!
!  B(xyz,tuv) = Siguma_w eps(w){-<0|EzyEtxEuvEww|0> + d(tx)<0|EzyEuvEww|0>}
!
!               + S(xyz,tuv)(eps(u)+eps(t)-eps(v))
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer :: it, iu, iv, ix, iy, iz, iw
    integer :: jt, ju, jv, jx, jy, jz, jw
    integer :: i, j

    integer, intent(in)     :: dimn, indsym(3, dimn)
    complex*16, intent(in)  :: sc(dimn, dimn)
    complex*16, intent(out) :: bc(dimn, dimn)

    real*8               :: e, denr, deni
    complex*16           :: den
!           write(*,*)'sc0',sc(5,5)

    bc(:, :) = 0.0d+00
    if (rank == 0) then
        write (*, *) 'bAmat loop: rank, nprocs, dimn', rank, nprocs, dimn
    end if
    !$OMP parallel do private(ix,iy,iz,jx,jy,jz,it,iu,iv,jt,ju,jv,e,j,iw,jw,denr,deni,den)
    Do i = rank + 1, dimn, nprocs ! MPI parallelization (Distributed loop: static scheduling, per nprocs)
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

!  B(xyz,tuv) = Siguma_w eps(w){-<0|EzyEtxEuvEww|0> + d(tx)<0|EzyEuvEww|0>}
!
!               + S(xyz,tuv)(eps(u)+eps(t)-eps(v))

            e = eps(ju) + eps(jt) - eps(jv)

            Do iw = 1, nact
                jw = iw + ninact

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

        End do               !i
    End do                  !j
    !$OMP end parallel do

! Noda : なぜかMPI_Reduce ver.だとisym=2のsAmatで実行が止まるのでMPI_Allreduce
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, bc(1, 1), dimn**2, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
!    if (rank == 0) then
!        call MPI_Reduce(MPI_IN_PLACE, bc(1, 1), dimn**2, MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!    else
!        call MPI_Reduce(bc(1, 1), bc(1, 1), dimn**2, MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!    end if

    if (rank == 0) then ! Process limits for output
        write (*, *) 'bAmat is ended'
    end if
End subroutine bAmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
SUBROUTINE vAmat_ord_ty(v)
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

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    complex*16, intent(out) ::  &
    & v(ninact, ninact + 1:ninact + nact, ninact + 1:ninact + nact, ninact + 1:ninact + nact)

    real*8                  :: dr, di, signkl
    complex*16              :: cint2, d, dens1(nact, nact), effh(ninact + 1:ninact + nact, ninact)
    complex*16              :: cint1

    integer :: it, iu, iv, ii, ip, iq, ir, ik
    integer :: jt, ju, jv, ji, jp, jq, jr, jk
    integer :: i, j, k, l, kkr, lkr, count, dim(nsymrpa)
    integer :: dim2(nsymrpa), isym, sym, i0, syma, symb, symc

    integer, allocatable :: indt(:, :), indu(:, :), indv(:, :)
    integer, allocatable :: ind2u(:, :), ind2v(:, :)
    logical :: test
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1
    integer :: loopcnt, idx

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
    if (rank == 0) write (*, *) 'Enter vAmat. Please ignore timer under this line.'
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0
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
        ! !$OMP parallel do schedule(static,1) private(it,ij,iv,ij,iu,ju,syma,symb,symc) reduction(+:dim)
        Do it = 1, nact
            jt = it + ninact
            Do iv = 1, nact
                jv = iv + ninact
                Do iu = 1, nact
                    ju = iu + ninact

! EtiEuv
                    syma = MULTB_D(irpmo(jt), isym)
                    symb = MULTB_D(irpmo(ju), irpmo(jv))
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
        if (rank == 0) write (*, *) 'solvA: isym, dim(isym)', isym, dim(isym)
    end do
    Allocate (ind2u(nact**2, nsymrpa)); Call memplus(KIND(ind2u), SIZE(ind2u), 1)
    Allocate (ind2v(nact**2, nsymrpa)); Call memplus(KIND(ind2v), SIZE(ind2v), 1)
    ind2u = 0.0d+00
    ind2v = 0.0d+00
    dim2 = 0
    !$OMP parallel do schedule(static) private(iu,ju,iv,jv,syma)
    Do isym = 1, nsymrpa
        Do iu = 1, nact
            ju = iu + ninact
            Do iv = 1, nact
                jv = iv + ninact

                syma = MULTB_D(irpmo(jv), irpmo(ju))

                if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                    dim2(isym) = dim2(isym) + 1
                    ind2u(dim2(isym), isym) = iu
                    ind2v(dim2(isym), isym) = iv
                end if

            End do
        End do
    End do
    !$OMP end parallel do

    Do isym = 1, nsymrpa
        if (rank == 0) then ! Process limits for output
            write (*, '(2I4)') dim2(isym), isym
        end if
    End do
    !$OMP parallel do private(ji,it,jt,cint1)
    Do ii = rank + 1, ninact, nprocs
        ji = ii
        Do it = 1, nact
            jt = it + ninact

            Call tramo1_ty(jt, ji, cint1)
            effh(jt, ji) = cint1
!              if(jt==11.and.ji==1) write(*,'("eff 1int",2I4,2E20.10)') jt,ji,cint1
!              if(jt==11.and.ji==1) write(*,'("eff 1int",2E20.10)') effh(jt,ji)
        End do
    End do
    !$OMP end parallel do

!      write(*,*)'effh(11,1)',effh(11,1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Two types of integrals are stored
!
!  (21|22) stored (pi|qr) ...TYPE 1
!  (21|11) stored (pi|jk) ...TYPE 2
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!   open(1, file ='A1int', status='old', form='unformatted')
    ! open (1, file=a1int, status='old', form='formatted')
    open (1, file=a1int, status='old', form='unformatted')
    if (rank == 0) then ! Process limits for output
        write (*, *) 'open A1int'
    end if

! 30  read (1, '(4I4, 2e20.10)', err=10, end=20) i, j, k, l, cint2 !  (ij|kl)
30  read (1, err=10, end=20) i, j, k, l, cint2 !  (ij|kl)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  V(tuv,i)=  - SIGUMA_p,q,r:act <0|EvuEptEqr|0>(pi|qr)
!
!  + SIGUMA_p,q:act  <0|EvuEpq|0> (ti|pq)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           write(*,'("TYPE 1  ",4I4,2E20.10)')i,j,k,l,cint2

    isym = irpmo(j)
    !$OMP parallel do private(it,iu,iv,jt,ju,jv,dr,di,d)
    Do i0 = 1, dim(isym)
        it = indt(i0, isym)
        iu = indu(i0, isym)
        iv = indv(i0, isym)
        jt = it + ninact
        ju = iu + ninact
        jv = iv + ninact

        Call dim3_density(iv, iu, i - ninact, it, k - ninact, l - ninact, dr, di)
        d = DCMPLX(dr, di)
        v(j, jt, ju, jv) = v(j, jt, ju, jv) - cint2*d

    End do
    !$OMP end parallel do

    isym = MULTB_D(irpmo(i), irpmo(j))           ! j coresponds to ii, i coresponds to it

    !$OMP parallel do private(iu,iv,ju,jv,dr,di,d)
    Do i0 = 1, dim2(isym)
        iu = ind2u(i0, isym)
        iv = ind2v(i0, isym)
        ju = iu + ninact
        jv = iv + ninact

        Call dim2_density(iv, iu, k - ninact, l - ninact, dr, di)
        d = DCMPLX(dr, di)
        v(j, i, ju, jv) = v(j, i, ju, jv) + cint2*d
    End do
    !$OMP end parallel do

    goto 30

20  close (1)

    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0

!  open (1, file='A2int', status='old', form='unformatted') ! TYPE 2 integrals
    ! open (1, file=a2int, status='old', form='formatted') ! TYPE 2 integrals
    open (1, file=a2int, status='old', form='unformatted') ! TYPE 2 integrals

! 300 read (1, '(4I4, 2e20.10)', err=10, end=200) i, j, k, l, cint2 !  (ij|kl)
300 read (1, err=10, end=200) i, j, k, l, cint2 !  (ij|kl)
    count = 0

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  effh(p,i) = h(pi)+ SIGUMA_k:inact{(pi|kk)-(pk|ki)}
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (k == l .and. j /= k) then       ! (PI|KK) type

        effh(i, j) = effh(i, j) + cint2
!           write(*,'("A2int+",4I4,2E20.10)')i,j,k,l,cint2

    elseif (j == k .and. k /= l) then       ! (PK|KI) type

        effh(i, l) = effh(i, l) - cint2
!           write(*,'("A2int-",4I4,2E20.10)')i,j,k,l,cint2

    end if

    goto 300

200 close (1)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'reading A2int2 is over'
    end if

!    effh(ninact + 1:ninact + nact, ninact)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, effh(ninact + 1, 1), nact*ninact, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0

!  - SIGUMA_p:act <0|EvuEpt|0>effh(pi)  +  <0|Evu|0>effh(ti)

    !$OMP parallel do private(ji,isym,it,iu,iv,jt,ju,jv,dr,di,d,ip,jp)
    Do ii = rank + 1, ninact, nprocs
        ji = ii
        isym = irpmo(ji)

!           Do ip = 1, nact
!              jp = ip + ninact
!              if(ABS(effh(jp,ji)) > 1.0d-10) write(*,'("o effh ",2I4,2E20.10)')jp,ji,effh(jp,ji)
!           Enddo

        Do i0 = 1, dim(isym)
            it = indt(i0, isym)
            iu = indu(i0, isym)
            iv = indv(i0, isym)
            jt = it + ninact
            ju = iu + ninact
            jv = iv + ninact

            Call dim1_density(iv, iu, dr, di)

            d = DCMPLX(dr, di)
            v(ji, jt, ju, jv) = v(ji, jt, ju, jv) + effh(jt, ji)*d

            Do ip = 1, nact
                jp = ip + ninact

                Call dim2_density(iv, iu, ip, it, dr, di)
                d = DCMPLX(dr, di)
                v(ji, jt, ju, jv) = v(ji, jt, ju, jv) - effh(jp, ji)*d

            End do            ! ip

        End do               !i0
    End do                  !ii
    !$OMP end parallel do

    goto 100

10  write (*, *) 'error while opening file Aint'; goto 1000
100 continue

1000 if (rank == 0) write (*, *) 'vAmat_ord_ty is ended'

    deallocate (indt); Call memminus(KIND(indt), SIZE(indt), 1)
    deallocate (indu); Call memminus(KIND(indu), SIZE(indu), 1)
    deallocate (indv); Call memminus(KIND(indv), SIZE(indv), 1)
    deallocate (ind2u); Call memminus(KIND(ind2u), SIZE(ind2u), 1)
    deallocate (ind2v); Call memminus(KIND(ind2v), SIZE(ind2v), 1)

!    v(ninact, ninact + 1:ninact + nact, ninact + 1:ninact + nact, ninact + 1:ninact + nact)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, v(ninact, ninact + 1, ninact + 1, ninact + 1), ninact*nact*nact*nact, &
                       MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (*, *) 'end allreduce vAmat'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
end subroutine vAmat_ord_ty
