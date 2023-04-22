! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE solvG_ord_ty(e0, e2g)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    real*8, intent(in) :: e0
    real*8, intent(out):: e2g
    integer :: dimn, dimm, dammy
    real*8, allocatable  :: wsnew(:), ws(:), wb(:)
    real*8               :: e2(2*nsymrpa), alpha, e
    complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
    complex*16, allocatable  :: bc(:, :)
    complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :), vc(:), vc1(:)
    logical                  :: cutoff
    integer                  :: j, i, i0, syma, symb, isym, indt(1:nact)
    integer                  :: ia, it, ib, ii, ja, jt, jb, ji
    integer, allocatable     :: ia0(:), ib0(:), ii0(:), iabi(:, :, :)
    integer                  :: nabi
    real*8  :: thresd
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE G IS NOW CALCULATED
!
!     EaiEbt|0>
!
!   DRAS1 =-1   DRAS2 = -1   DRAS3 = +2
!
!   c > d, a > b, and impose that c >= a (or a >= c)
!
!  S(cjdu,aibt) = d(ac) d(bd) d(ij) <0|Eut|0>  <= S(u,t)
!                                   ~~~~~~~~~
!  S(u,t) = <0|Eut|0>
!
!  B(u,t) = Siguma_w [eps(w)<0|EutEww|0>] + S(u,t)(-eps(t))
!
!  alpha(i,a,b) = -eps(i) + eps(a) + eps(b) - e0
!
!  where
!
!  e0 = Siguma_w [eps(w)<0|Eww|0>] (<== calculated as e0 in calce0.f)
!
!  V(t,iab)   =  [SIGUMA_p:active <0|Etp|0>{(ai|bp)-(ap|bi)}]
!
!
!  E2 = SIGUMA_iab, dimm |V1(t,iab)|^2|/{(alpha(iab) + wb(t)}
!
    if (rank == 0) then
        print *, ' ENTER solv G part'
        print *, ' nsymrpa', nsymrpa
    end if
    datetmp1 = date0; datetmp0 = date0

    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0
    thresd = 1.0D-08
    thres = 1.0D-08

    e2 = 0.0d+00
    e2g = 0.0d+00
    dimn = 0
    indt = 0

    i0 = 0
    Do ia = 1, nsec
        ja = ia + ninact + nact
        Do ib = 1, ia - 1
            jb = ib + ninact + nact
            Do ii = 1, ninact
                ji = ii
                i0 = i0 + 1
            End do
        End do
    End do

    nabi = i0
    Allocate (iabi(nsec, nsec, ninact))
    iabi = 0
    Allocate (ia0(nabi))
    Allocate (ib0(nabi))
    Allocate (ii0(nabi))

    i0 = 0
    Do ia = 1, nsec
        ja = ia + ninact + nact
        Do ib = 1, ia - 1
            jb = ib + ninact + nact
            Do ii = 1, ninact
                ji = ii
                i0 = i0 + 1
                iabi(ia, ib, ii) = i0
                iabi(ib, ia, ii) = i0
                ia0(i0) = ia + ninact + nact ! secondary
                ib0(i0) = ib + ninact + nact ! secondary
                ii0(i0) = ii
            End do
        End do
    End do

    Allocate (v(nabi, nact))
    v = 0.0d+00
    if (rank == 0) print *, 'end before v matrices'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Call vGmat_ord_ty(nabi, iabi, v)
    if (rank == 0) print *, 'end after vGmat'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0

    Do isym = 1, nsymrpa

        dimn = 0
        Do it = 1, nact
            jt = it + ninact
            if (irpamo(jt) == isym) then
                dimn = dimn + 1
                indt(dimn) = it
            End if
        End do                  ! it

        if (rank == 0) print *, 'isym, dimn', isym, dimn
        If (dimn == 0) cycle ! Go to the next isym

        Allocate (sc(dimn, dimn))
        sc = 0.0d+00            ! sc N*N
        if (rank == 0) print *, 'before sGmat'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call sGmat(dimn, indt(1:dimn), sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rank == 0) print *, 'sG matrix is obtained normally'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Allocate (ws(dimn))

        cutoff = .TRUE.
!           thresd = 1.0d-15

        Allocate (sc0(dimn, dimn))
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

        Allocate (bc(dimn, dimn))                                 ! bc N*N
        bc = 0.0d+00
        if (rank == 0) print *, 'before bGmat'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call bGmat(dimn, sc0, indt(1:dimn), bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (rank == 0) print *, 'bC matrix is obtained normally'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        deallocate (sc0)

        if (rank == 0) print *, 'OK cdiag', dimn, dimm
        Allocate (uc(dimn, dimm))                                 ! uc N*M
        Allocate (wsnew(dimm))                                  ! wnew M
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

        if (rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'
        Allocate (bc0(dimm, dimm))
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

        e2 = 0.0d+00

        Do i0 = 1, nabi
            ja = ia0(i0)
            jb = ib0(i0)
            ji = ii0(i0)

!     EaiEbt|0>
            if (nsymrpa /= 1) then
                syma = MULTB_D(irpamo(jb), isym)
                symb = MULTB_D(irpamo(ja), irpamo(ji))
                syma = MULTB_S(syma, symb)
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

                alpha = -eps(ji) + eps(ja) + eps(jb) - e0 + eshift  ! For Level Shift (2007/2/9)

                vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm, 1:dimm))), vc1(1:dimm))

                Do j = 1, dimm
                    e = (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                    sumc2local = sumc2local + e/(alpha + wb(j))
                    e2(isym) = e2(isym) - e
                End do

                deallocate (vc1)

            End if

        End do

        deallocate (uc)
        deallocate (wb)
        Deallocate (bc1)

        if (rank == 0) print '("e2g(",I3,") = ",E20.10," a.u.")', isym, e2(isym)
        e2g = e2g + e2(isym)
        if (rank == 0) print *, 'End e2(isym) add'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
    End do                  ! isym

    if (rank == 0) then
        print '("e2g      = ",E20.10," a.u.")', e2g
        print '("sumc2,g  = ",E20.10)', sumc2local
    end if
    sumc2 = sumc2 + sumc2local

    deallocate (iabi)
    deallocate (ia0)
    deallocate (ib0)
    deallocate (ii0)
    deallocate (v)

    continue
    if (rank == 0) print *, 'end solvG_ord_ty'
end

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE sGmat(dimn, indt, sc) ! Assume C1 molecule, overlap matrix S in space C

!  S(u,t) = <0|Eut|0>
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer, intent(in)      :: dimn, indt(dimn)
    complex*16, intent(out)  :: sc(dimn, dimn)
    real*8  ::a, b
    integer :: it, iu
    integer :: i, j

    sc = 0.0d+00

!    !$OMP parallel do schedule(dynamic,1) private(it,iu,a,b)
    Do i = rank + 1, dimn, nprocs
        it = indt(i)

        Do j = i, dimn
            iu = indt(j)
            a = 0.0d+0
            b = 0.0d+0

            Call dim1_density(it, iu, a, b)

            sc(i, j) = DCMPLX(a, b)
            sc(j, i) = DCMPLX(a, -b)
        End do               !j
    End do                  !i
!    !$OMP end parallel do
#ifdef HAVE_MPI
    call allreduce_wrapper(mat=sc)
#endif

End subroutine sGmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE bGmat(dimn, sc, indt, bc) ! Assume C1 molecule, overlap matrix B in space C
!
!
!  S(u,t) = <0|Eut|0>
!
!  B(u,t) = Siguma_w [eps(w)<0|EutEww|0>] + S(u,t)(-eps(t))
!
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer :: it, iu, iw, jt, ju, jw
    integer :: i, j

    integer, intent(in) :: dimn, indt(dimn)
    complex*16, intent(in)  :: sc(dimn, dimn)
    complex*16, intent(out) :: bc(dimn, dimn)

    real*8              :: denr, deni
    complex*16          :: den

    bc(:, :) = 0.0d+00

    if (rank == 0) print *, 'G space Bmat iroot=', iroot

    !  !$OMP parallel do schedule(dynamic,1) private(iu,ju,j,it,jt,iw,jw,denr,deni,den)
    Do i = rank + 1, dimn, nprocs
        iu = indt(i)
        ju = iu + ninact

        Do j = i, dimn
            it = indt(j)
            jt = it + ninact

!  B(u,t) = Siguma_w [eps(w)<0|EutEww|0>] + S(u,t)(-eps(t))

            Do iw = 1, nact
                jw = iw + ninact

                Call dim2_density(iu, it, iw, iw, denr, deni)
                den = DCMPLX(denr, deni)
                bc(i, j) = bc(i, j) + den*eps(jw)

            End do

            bc(i, j) = bc(i, j) - sc(i, j)*eps(ju)

            bc(j, i) = DCONJG(bc(i, j))

        End do               !i
    End do                  !j
    !   !$OMP end parallel do
#ifdef HAVE_MPI
    call reduce_wrapper(mat=bc, root_rank=0)
#endif
    if (rank == 0) print *, 'bGmat is ended'
End subroutine bGmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE vGmat_ord_ty(nabi, iabi, v)
!
!  V(t,iab)   =  [SIGUMA_p:active <0|Etp|0>{(ai|bp)-(ap|bi)}]
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
    use module_file_manager
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer, intent(in)     :: nabi, iabi(nsec, nsec, ninact)

    complex*16, intent(out) :: v(nabi, nact)

    real*8                  :: dr, di
    complex*16              :: cint2, dens

    integer :: i, j, k, l, tabi
    integer :: it, iostat, unit_int2
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1
    logical :: is_end_of_file

    if (rank == 0) print *, 'Enter vGmat. Please ignore timer under this line.'
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0
    v = 0.0d+00

!  V(t,iab)   =  [SIGUMA_p:active <0|Etp|0>{(ai|bp)-(ap|bi)}]       a > b

    call open_unformatted_file(unit=unit_int2, file=gint, status='old', optional_action='read') !  (31|32) stored
    do
        read (unit_int2, iostat=iostat) i, j, k, l, cint2
        call check_iostat(iostat=iostat, file=gint, end_of_file_reached=is_end_of_file)
        if (is_end_of_file) then
            exit
        end if
        if (i == k) cycle ! Go to the next line if i == k

        tabi = iabi(i, k, j)

        if (i < k) then
            cint2 = -1.0d+00*cint2
        end if

        Do it = 1, nact
            Call dim1_density(it, l, dr, di)
            dens = DCMPLX(dr, di)
            v(tabi, it) = v(tabi, it) + cint2*dens
        End do                  ! it

    end do
    close (unit_int2)

    if (rank == 0) print *, 'vGmat_ord_ty is ended'
#ifdef HAVE_MPI
    call allreduce_wrapper(mat=v)
    if (rank == 0) print *, 'end allreduce vGmat'
#endif
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
end subroutine vGmat_ord_ty
