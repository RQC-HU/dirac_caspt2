! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE solvD_ord_ty(e0, e2d)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    real*8, intent(in) :: e0
    real*8, intent(out):: e2d
    integer :: dimn, dimm, dammy
    integer, allocatable :: indsym(:, :)
    real*8, allocatable  :: wsnew(:), ws(:), wb(:)
    real*8               :: e2(nsymrpa*2), e, alpha
    complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
    complex*16, allocatable  :: bc(:, :)
    complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :), vc(:), vc1(:)
    integer, allocatable     :: ia0(:), ii0(:), iai(:, :)
    integer                  :: nai
    logical :: cutoff
    integer :: j, i, syma, isym, i0
    integer :: ia, it, ii, iu
    integer :: ja, jt, ji, ju
    real*8  :: thresd
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE D IS NOW CALCULATED
!
!     EaiEtu|0>
!
!   DRAS1 = -1   DRAS2 = 0   DRAS3 = +1
!
!   t and u run all active spinor space independently!
!
!
!  S(bjxy,aitu) = d(ba)d(ij)<0|EyxEtu|0>
!
!  S(xy, tu) = <0|EyxEtu|0>
!
!  B(xy,tu) = Siguma_w [eps(w)<0|EyxEtuEww|0>+S(xy,tu)(eps(t)-eps(u))]
!
!  a(a,i)       = eps(a) - eps(i) - e0
!
!! dame!! V(a,i) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!! dame!!
!! dame!! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) + SIGUMA_p:active(ap|pi)}]
!
!
! V(a,i) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki) }
!
!
!  E2 = SIGUMA_a,i, dimm |V1(dimm,ai)|^2|/{(alpha(ai) + wb(dimm)}

    if (rank == 0) then ! Process limits for output
        write (*, *) ' ENTER solv D part'
        write (*, *) ' nsymrpa', nsymrpa
    end if
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0;
    thresd = 1.0D-08
    thres = 1.0D-08

    e2 = 0.0d+00
    e2d = 0.0d+00
    dimn = 0
    syma = 0

    i0 = 0
    Do ia = 1, nsec
        Do ii = 1, ninact
            i0 = i0 + 1
        End do
    End do

    nai = i0
    Allocate (iai(nsec, ninact))
    iai = 0
    Allocate (ia0(nai))
    Allocate (ii0(nai))

    i0 = 0
    Do ia = 1, nsec
        Do ii = 1, ninact
            i0 = i0 + 1
            iai(ia, ii) = i0
            ia0(i0) = ia
            ii0(i0) = ii
        End do
    End do
    Allocate (v(nai, nact, nact))
    v = 0.0d+00
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (*, *) 'end before v matrices'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Call vDmat_ord_ty(nai, iai, v)
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (*, *) 'end after vDmat'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    if (rank == 0) then ! Process limits for output
        write (*, *) 'come'
    end if

    Do isym = 1, nsymrpa

        dimn = 0
        Do it = 1, nact
            jt = it + ninact
            Do iu = 1, nact
                ju = iu + ninact

                syma = MULTB_D(irpmo(jt), irpmo(ju))

                if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                    dimn = dimn + 1
                End if
100         End do               ! iu
        End do                  ! it

        if (rank == 0) print *, 'isym, dimn', isym, dimn

        If (dimn == 0) goto 1000

        Allocate (indsym(2, dimn))
        indsym = 0
        dimn = 0

        Do it = 1, nact
            jt = it + ninact
            Do iu = 1, nact
                ju = iu + ninact

                syma = MULTB_D(irpmo(jt), irpmo(ju))

                if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                    dimn = dimn + 1
                    indsym(1, dimn) = it
                    indsym(2, dimn) = iu
                End if
200         End do               ! iu
        End do                  ! it

        Allocate (sc(dimn, dimn))
        sc = 0.0d+00            ! sc N*N
        if (rank == 0) then ! Process limits for output
            write (*, *) 'before sDmat'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call sDmat(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (rank == 0) then ! Process limits for output
            write (*, *) 'sc matrix is obtained normally'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Allocate (ws(dimn))
        ws = 0.0d+00
        cutoff = .TRUE.
        thresd = 1.0d-08
!           thresd = 1.0d-15

        Allocate (sc0(dimn, dimn))
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
            write (*, *) 'after s cdiag'
        end if

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

        Allocate (bc(dimn, dimn))                                 ! bc N*N
        bc = 0.0d+00
        if (rank == 0) then ! Process limits for output
            write (*, *) 'before bDmat'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call bDmat(dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (rank == 0) then ! Process limits for output
            write (*, *) 'bc matrix is obtained normally'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        deallocate (sc0)

        if (rank == 0) then ! Process limits for output
            write (*, *) 'OK cdiag', dimn, dimm
        end if

        Allocate (uc(dimn, dimm))                                 ! uc N*M
        Allocate (wsnew(dimm))                                  ! wnew M
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

        if (rank == 0) then ! Process limits for output
            IF (debug) then

                write (*, *) 'Check whether bc1 is hermite or not'
                Do i = 1, dimm
                    Do j = i, dimm
                        if (ABS(bc1(i, j) - DCONJG(bc1(j, i))) > 1.0d-6) then
                            write (*, '(2I4,2E15.7)') i, j, bc1(i, j) - bc1(j, i)
                        End if
                    End do
                End do
                write (*, *) 'Check whether bc1 is hermite or not END'

            End if
        end if

        deallocate (bc)
        deallocate (bc0)

        cutoff = .FALSE.

        Allocate (wb(dimm))
        wb = 0.0d+00

        if (rank == 0) then ! Process limits for output
            write (*, *) 'bC matrix is transrated to bc1(M*M matrix)!'
        end if
        Allocate (bc0(dimm, dimm))
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
        e2 = 0.0d+00
        Do i0 = 1, nai
            ja = ia0(i0) + ninact + nact
            ji = ii0(i0)

            syma = MULTB_D(irpmo(ja), irpmo(ji))
            syma = MULTB_S(syma, isym)

            If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then

                Allocate (vc(dimn))
                Do it = 1, dimn
                    vc(it) = v(i0, indsym(1, it), indsym(2, it))
                End do

                Allocate (vc1(dimm))

                vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn, 1:dimm))), vc(1:dimn))
                Deallocate (vc)

                alpha = +eps(ja) - eps(ji) - e0 + eshift  ! For Level Shift (2007/2/9)

                vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm, 1:dimm))), vc1(1:dimm))

                Do j = 1, dimm
                    sumc2local = sumc2local + (ABS(vc1(j))**2.0d+00)/((alpha + wb(j))**2.0d+00)
                    e = (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                    e2(isym) = e2(isym) - e
                End do

                deallocate (vc1)

            End if

        End do               !i0

        deallocate (indsym)
        deallocate (uc)
        deallocate (wb)
        Deallocate (bc1)

1000    if (rank == 0) write (*, '("e2d(",I3,") = ",E20.10," a.u.")') isym, e2(isym)
        e2d = e2d + e2(isym)
        if (rank == 0) then ! Process limits for output
            write (*, *) 'End e2(isym) add'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
    End do                  ! isym

    if (rank == 0) then ! Process limits for output
        write (*, '("e2d      = ",E20.10," a.u.")') e2d

        write (*, '("sumc2,d  = ",E20.10)') sumc2local
    end if
    sumc2 = sumc2 + sumc2local

    deallocate (iai)
    deallocate (ia0)
    deallocate (ii0)
    deallocate (v)

    continue
    if (rank == 0) then ! Process limits for output
        write (*, *) 'end solvD_ord_ty'
    end if
end

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE sDmat(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space D

!  S(xy, tu) = <0|EyxEtu|0>
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer, intent(in)      :: dimn, indsym(2, dimn)
    complex*16, intent(out)  :: sc(dimn, dimn)
    real*8  :: a, b
    integer :: it, iu, iy, ix
    integer :: i, j

    sc = 0.0d+00

    !$OMP parallel do schedule(dynamic,1) private(i,ix,iy,j,it,iu,a,b)
    Do i = rank + 1, dimn, nprocs

        ix = indsym(1, i)
        iy = indsym(2, i)
        Do j = i, dimn

            it = indsym(1, j)
            iu = indsym(2, j)

            a = 0.0d+0
            b = 0.0d+0
            Call dim2_density(iy, ix, it, iu, a, b)
            sc(i, j) = DCMPLX(a, b)
            sc(j, i) = DCONJG(sc(i, j))

        End do               !j
    End do                  !i
    !$OMP end parallel do
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, sc(1, 1), dimn**2, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
End subroutine sDmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE bDmat(dimn, sc, indsym, bc) ! Assume C1 molecule, overlap matrix B in space D
!
!  B(xy,tu) = Siguma_w [eps(w)<0|EyxEtuEww|0>+S(xy,tu)(eps(t)-eps(u))]
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

    real*8              :: e, denr, deni
    complex*16          :: den

    integer :: it, iu, ix, iy, iw
    integer :: jt, ju, jy, jx, jw, i, j

    bc(:, :) = 0.0d+00

    if (rank == 0) then ! Process limits for output
        write (*, *) 'F space Bmat iroot=', iroot
    end if
    !$OMP parallel do schedule(dynamic,1) private(ix,iy,jx,jy,it,iu,jt,ju,e,j,iw,jw,denr,deni,den)
    Do i = rank + 1, dimn, nprocs

        ix = indsym(1, i)
        jx = ix + ninact
        iy = indsym(2, i)
        jy = iy + ninact

        Do j = i, dimn

            it = indsym(1, j)
            jt = it + ninact
            iu = indsym(2, j)
            ju = iu + ninact

!  B(xy,tu) = Siguma_w [eps(w)<0|EyxEtuEww|0>+S(xy,tu)(eps(t)-eps(u))]

            e = eps(jt) - eps(ju)

            Do iw = 1, nact
                jw = iw + ninact

                Call dim3_density(iy, ix, it, iu, iw, iw, denr, deni)
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
        write (*, *) 'bDmat is ended'
    end if
End subroutine bDmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
SUBROUTINE vDmat_ord_ty(nai, iai, v)
!
!
! V(a,i) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer, intent(in)     :: nai, iai(nsec, ninact)
    complex*16, intent(out) :: v(nai, nact, nact)
    real*8                  :: dr, di
    complex*16              :: cint1, cint2,  d
    complex*16              :: effh(nsec, ninact)
    integer :: i, j, k, l, tai
    integer :: it, jt, ju, iu, ia, ii, ja, ji
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1

    if (rank == 0) write (*, *) 'Enter vDmat. Please ignore timer under this line.'
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0
    v = 0.0d+00
    effh = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! V(tai, jt, ju) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
    !$OMP parallel do schedule(dynamic,1) private(ia,ja,ii,ji,cint1)
    Do ia = rank + 1, nsec, nprocs
        ja = ia + ninact + nact
        Do ii = 1, ninact
            ji = ii
            Call tramo1_ty(ja, ji, cint1)
            effh(ia, ii) = cint1
        End do
    End do
    !$OMP end parallel do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! V(a,i, jt, ju) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Three types of integrals are stored
!
!  (31|22) stored (ai|pq) ...TYPE 1 D1int
!  (32|21) stored (ap|qi) ...TYPE 2 D2int
!
!  (31|11) stored (ai|jk) ...TYPE 3 D3int
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (rank == 0) write (*, *) 'before d1int'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    open (1, file=d1int, status='old', form='unformatted')

30  read (1, err=10, end=20) i, j, k, l, cint2 !  (ij|kl)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! V(a,i, jt, ju) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ja = i
    ji = j
    tai = iai(ja, ji)

    !$OMP parallel do schedule(dynamic,1) private(it,jt,iu,ju,dr,di,d)
    Do it = 1, nact
        jt = it + ninact
        Do iu = 1, nact
            ju = iu + ninact

            Call dim2_density(iu, it, k, l, dr, di)
            d = DCMPLX(dr, di)
            v(tai, it, iu) = v(tai, it, iu) + d*cint2

        End do
    End do
    !$OMP end parallel do

    goto 30
20  close (1)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'reading D1int2 is over'
    end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! V(a,i, jt, ju) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (rank == 0) write (*, *) 'before d2int'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    open (1, file=d2int, status='old', form='unformatted')

31  read (1, err=10, end=21) i, j, k, l, cint2 !  (ij|kl)
    ja = i
    ji = l
    tai = iai(ja, ji)
    !$OMP parallel do schedule(dynamic,1) private(it,jt,iu,ju,dr,di,d)
    Do it = 1, nact
        jt = it + ninact
        Do iu = 1, nact
            ju = iu + ninact

            Call dim2_density(iu, it, k, j, dr, di)
            d = DCMPLX(dr, di)

            v(tai, it, iu) = v(tai, it, iu) - d*cint2

        End do
    End do
    !$OMP end parallel do

    goto 31

21  close (1)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'reading D2int2 is over'
    end if
    if (rank == 0) write (*, *) 'before d3int'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    open (1, file=d3int, status='old', form='unformatted') ! (ai|jk) is stored

300 read (1, err=10, end=200) i, j, k, l, cint2 !  (ij|kl)

    if (j /= k .and. k == l) then !(ai|kk)

        effh(i, j) = effh(i, j) + cint2

    elseif (j == k .and. k /= l) then !(ak|ki)

        effh(i, l) = effh(i, l) - cint2

    end if

    goto 300

200 close (1)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'reading D3int2 is over'
    end if
    if (rank == 0) write (*, *) 'end d3int'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0

#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, effh(1, 1), nsec*ninact, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (*, *) 'end allreduce effh'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0

    !$OMP parallel do schedule(dynamic,1) private(ia,ja,ii,ji,tai,it,jt,iu,ju,dr,di,d)
    Do ia = rank + 1, nsec, nprocs
        ja = ia + ninact + nact
        Do ii = 1, ninact
            ji = ii
            tai = iai(ia, ii)

            Do it = 1, nact
                jt = it + ninact
                Do iu = 1, nact
                    ju = iu + ninact

                    Call dim1_density(iu, it, dr, di)

                    d = DCMPLX(dr, di)
                    v(tai, it, iu) = v(tai, it, iu) + effh(ia, ii)*d
                End do
            End do

        End do
    End do
    !$OMP end parallel do

    goto 100

10  write (*, *) 'error while opening file Dint'; goto 100

100 if (rank == 0) write (*, *) 'vDmat_ord_ty is ended'
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, v(1, 1, 1), nai*nact**2, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (*, *) 'end Allreduce vDmat'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
end subroutine vDmat_ord_ty
