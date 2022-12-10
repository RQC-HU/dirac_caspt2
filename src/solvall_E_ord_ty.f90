! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE solvE_ord_ty(e0, e2e)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    real*8, intent(in) :: e0
    real*8, intent(out):: e2e
    integer :: dimn, dimm, dammy
    real*8, allocatable  :: wsnew(:), ws(:), wb(:)
    real*8               :: e2(2*nsymrpa), alpha, e
    complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
    complex*16, allocatable  :: bc(:, :)
    complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :), vc(:), vc1(:)
    logical :: cutoff
    integer :: j, i, syma, symb, isym, indt(1:nact)
    integer :: ia, it, ij, ii, ja, jt, jj, ji
    integer :: i0
    integer, allocatable     :: ia0(:), ii0(:), ij0(:), iaij(:, :, :)
    integer :: naij
    real*8  :: thresd
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1

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
    thresd = 1.0D-08
    thres = 1.0D-08

    e2 = 0.0d+00
    e2e = 0.0d+00
    dimn = 0
    syma = 0
    indt = 0
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0
    if (rank == 0) then
        print *, ' ENTER solv E part'
        print *, ' nsymrpa', nsymrpa
    end if
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
                ia0(i0) = ia + ninact + nact ! secondary
                ii0(i0) = ii ! inactive
                ij0(i0) = ij ! inactive
            End do
        End do
    End do

    Allocate (v(naij, nact))
    v = 0.0d+00
    if (rank == 0) print *, 'end before v matrices'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Call vEmat_ord_ty(naij, iaij, v)
    if (rank == 0) print *, 'come'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0

    Do isym = 1, nsymrpa

        dimn = 0
        Do it = 1, nact
            jt = it + ninact
            if (irpmo(jt) == isym) then
                dimn = dimn + 1
                indt(dimn) = it
            End if
        End do                  ! it

        if (rank == 0) print *, 'isym, dimn', isym, dimn
        If (dimn == 0) cycle ! Go to the next isym

        Allocate (sc(dimn, dimn))
        sc = 0.0d+00            ! sc N*N
        if (rank == 0) print *, 'before sEmat'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call sEmat(dimn, indt(1:dimn), sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (rank == 0) print *, 'sc matrix is obtained normally'
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
        if (rank == 0) print *, 'before bEmat'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call bEmat(e0, dimn, sc0, indt(1:dimn), bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (rank == 0) print *, 'bc matrix is obtained normally'
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

        Do i0 = 1, naij
            ja = ia0(i0)
            ji = ii0(i0)
            jj = ij0(i0)

!     EtiEaj|0>

            if (nsymrpa /= 1) then
                syma = MULTB_D(irpmo(ja), irpmo(jj))
                symb = MULTB_D(isym, irpmo(ji))
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
                    sumc2local = sumc2local + e/(alpha + wb(j))
                    e2(isym) = e2(isym) - e
                End do

                deallocate (vc1)

            End if

        End do

        deallocate (uc)
        deallocate (wb)
        Deallocate (bc1)

        if (rank == 0) print '("e2e(",I3,") = ",E20.10," a.u.")', isym, e2(isym)
        e2e = e2e + e2(isym)
        if (rank == 0) print *, 'End e2(isym) add'
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
    End do                  ! isym

    if (rank == 0) then
        print '("e2e      = ",E20.10," a.u.")', e2e
        print '("sumc2,e  = ",E20.10)', sumc2local
    end if
    sumc2 = sumc2 + sumc2local

    deallocate (iaij)
    deallocate (ia0)
    deallocate (ii0)
    deallocate (ij0)
    deallocate (v)

    continue
    if (rank == 0) print *, 'end solveE_ord_ty'
end

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE sEmat(dimn, indt, sc) ! Assume C1 molecule, overlap matrix S in space E

! S(u,t) = d(ut) - <0|Etu|0>
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

        End do               !j
    End do                  !i
    !$OMP end parallel do
#ifdef HAVE_MPI
    call allreduce_wrapper(mat=sc)
#endif

End subroutine sEmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE bEmat(e0, dimn, sc, indt, bc) ! Assume C1 molecule, overlap matrix B in space E
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

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    integer :: it, iu, iw, jt, ju, jw
    integer :: i, j

    integer, intent(in) :: dimn, indt(dimn)
    complex*16, intent(in)  :: sc(dimn, dimn)
    complex*16, intent(out) :: bc(dimn, dimn)
    real*8, intent(in)      :: e0

    real*8              :: denr, deni
    complex*16          :: den

    bc(:, :) = 0.0d+00

    if (rank == 0) print *, 'E space Bmat iroot=', iroot

    !$OMP parallel do schedule(dynamic,1) private(iu,ju,j,it,jt,iw,jw,denr,deni,den)
    Do i = rank + 1, dimn, nprocs
        iu = indt(i)
        ju = iu + ninact

        Do j = i, dimn
            it = indt(j)
            jt = it + ninact

            Do iw = 1, nact
                jw = iw + ninact

!         = - Siguma_w [eps(w)<0|EtuEww|0>] + d(tu)e0 + S(u,t)eps(t)

                Call dim2_density(it, iu, iw, iw, denr, deni)
                den = DCMPLX(denr, deni)
                bc(i, j) = bc(i, j) - den*eps(jw)

            End do

            if (it == iu) bc(i, j) = bc(i, j) + e0

            bc(i, j) = bc(i, j) + sc(i, j)*eps(jt)

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
    if (rank == 0) print *, 'bEmat is ended'
End subroutine bEmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE vEmat_ord_ty(naij, iaij, v)
!
!  V(t,ija)   =[SIGUMA_p:active <0|Ept|0>{(ai|pj) - (aj|pi)}] + (aj|ti) - (ai|tj)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
    use module_file_manager
#ifdef HAVE_MPI
    use module_mpi
#endif
    Implicit NONE

    integer, intent(in)     :: naij, iaij(nsec, ninact, ninact)

    complex*16, intent(out) :: v(naij, nact)

    real*8                  :: dr, di
    complex*16              :: cint2, dens

    integer :: i, j, k, l, taij
    integer :: it, ik, iostat, twoint_unit
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1
    logical :: is_end_of_file

    if (rank == 0) print *, 'Enter vEmat. Please ignore timer under this line.'
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0
    v = 0.0d+00
    twoint_unit = default_unit

!  V(t,ija)   =[SIGUMA_p:active <0|Ept|0>{(ai|pj) - (aj|pi)}] - (ai|tj) + (aj|ti)   i > j

    call open_unformatted_file(unit=twoint_unit, file=eint, status='old', optional_action='read') !  (31|21) stored
    do
        read (twoint_unit, iostat=iostat) i, j, k, l, cint2
        call check_iostat(iostat=iostat, file=eint, end_of_file_reached=is_end_of_file)
        if (is_end_of_file) then
            exit
        end if

        if (j == l) cycle ! Read the next 2-integral if j equal to l

        taij = iaij(i, j, l)
        ik = k - ninact

        if (j < l) then
            cint2 = -1.0d+00*cint2
        end if

        v(taij, k) = v(taij, k) - cint2

        !$OMP parallel do schedule(dynamic,1) private(it,dr,di,dens)
        Do it = 1, nact
            Call dim1_density(it, k, dr, di)          ! k corresponds to p in above formula
            dens = DCMPLX(dr, di)
            v(taij, it) = v(taij, it) + cint2*dens
        End do                  ! it
        !$OMP end parallel do

        if (j < l) then
            cint2 = -1.0d+00*cint2          ! data cint2 becomes initial values!
        end if

    end do
    close (twoint_unit)

    if (rank == 0) print *, 'vEmat_ord_ty is ended'

#ifdef HAVE_MPI
    call allreduce_wrapper(mat=v)
    if (rank == 0) print *, 'end Allreduce vEmat'
#endif
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
end subroutine vEmat_ord_ty
