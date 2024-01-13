SUBROUTINE solve_G_subspace(e0, e2g)

    use module_ulambda_s_half, only: ulambda_s_half
    use module_global_variables
    use module_realonly, only: realonly
    implicit none
    real(8), intent(in) :: e0
    real(8), intent(out):: e2g

    if (realonly%is_realonly()) then
        call solve_G_subspace_real()
    else
        call solve_G_subspace_complex()
    end if

contains
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_G_subspace_complex()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx, convert_secondary_to_global_idx

        Implicit NONE
#ifdef HAVE_MPI
        include 'mpif.h'
#endif

        integer :: dimn, dimm, dammy
        real(8), allocatable  :: wsnew(:), ws(:), wb(:)
        real(8)               :: e2(2*nsymrpa), alpha, e
        complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
        complex*16, allocatable  :: bc(:, :)
        complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :), vc(:), vc1(:)
        integer                  :: j, i, i0, syma, symb, isym, indt(1:nact)
        integer                  :: ia, it, ib, ii, ja, jt, jb, ji
        integer, allocatable     :: ia0(:), ib0(:), ii0(:), iabi(:, :, :)
        integer                  :: nabi

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

        e2 = 0.0d+00
        e2g = 0.0d+00
        dimn = 0
        indt = 0

        i0 = 0
        Do ia = 1, nsec
            ja = convert_secondary_to_global_idx(ia)
            Do ib = 1, ia - 1
                jb = convert_secondary_to_global_idx(ib)
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
            ja = convert_secondary_to_global_idx(ia)
            Do ib = 1, ia - 1
                jb = convert_secondary_to_global_idx(ib)
                Do ii = 1, ninact
                    ji = ii
                    i0 = i0 + 1
                    iabi(ia, ib, ii) = i0
                    iabi(ib, ia, ii) = i0
                    ia0(i0) = convert_secondary_to_global_idx(ia) ! secondary
                    ib0(i0) = convert_secondary_to_global_idx(ib) ! secondary
                    ii0(i0) = ii
                End do
            End do
        End do

        Allocate (v(nabi, nact))
        v = 0.0d+00
        if (rank == 0) print *, 'end before v matrices'
        Call vGmat_complex(nabi, iabi, v)
        if (rank == 0) print *, 'end after vGmat'

        Do isym = 1, nsymrpa

            dimn = 0
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                if (irpamo(jt) == isym) then
                    dimn = dimn + 1
                    indt(dimn) = it
                End if
            End do

            if (rank == 0) print *, 'isym, dimn', isym, dimn
            If (dimn == 0) cycle ! Go to the next isym

            Allocate (sc(dimn, dimn))
            sc = 0.0d+00            ! sc N*N
            if (rank == 0) print *, 'before sGmat'
            Call sGmat_complex(dimn, indt(1:dimn), sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'sG matrix is obtained normally'
            Allocate (ws(dimn))
            Allocate (sc0(dimn, dimn))
            sc0 = sc
            if (rank == 0) print *, 'before cdiag'
            Call cdiag(sc, dimn, dimm, ws, smat_lin_dep_threshold)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'after s cdiag, new dimension is', dimm
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
            Call bGmat_complex(dimn, sc0, indt(1:dimn), bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (rank == 0) print *, 'bC matrix is obtained normally'
            deallocate (sc0)

            if (rank == 0) print *, 'OK cdiag', dimn, dimm
            Allocate (uc(dimn, dimm))                                 ! uc N*M
            Allocate (wsnew(dimm))                                  ! wnew M
            uc(:, :) = 0.0d+00
            wsnew(:) = 0.0d+00
            if (rank == 0) print *, 'before ccutoff'
            Call ccutoff(sc, ws, dimn, dimm, smat_lin_dep_threshold, uc, wsnew)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'OK ccutoff'
            deallocate (ws)
            deallocate (sc)
            if (rank == 0) print *, 'before ulambda_s_half'
            Call ulambda_s_half(uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            deallocate (wsnew)

            if (rank == 0) print *, 'ucrams half OK'
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

            if (rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'
            Allocate (bc0(dimm, dimm))
            bc0 = bc1
            if (rank == 0) print *, 'before cdiag'
            Call cdiag(bc1, dimm, dammy, wb, bmat_no_cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'end cdiag'
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
        End do

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
        if (rank == 0) print *, 'end solve_G_subspace'
    end subroutine solve_G_subspace_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE sGmat_complex(dimn, indt, sc) ! Assume C1 molecule, overlap matrix S in space C

!  S(u,t) = <0|Eut|0>
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)      :: dimn, indt(dimn)
        complex*16, intent(out)  :: sc(dimn, dimn)
        real(8)  ::a, b
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
            End do
        End do
!    !$OMP end parallel do
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=sc)
#endif

    End subroutine sGmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE bGmat_complex(dimn, sc, indt, bc) ! Assume C1 molecule, overlap matrix B in space C
!
!
!  S(u,t) = <0|Eut|0>
!
!  B(u,t) = Siguma_w [eps(w)<0|EutEww|0>] + S(u,t)(-eps(t))
!
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer :: it, iu, iw, jt, ju, jw
        integer :: i, j

        integer, intent(in) :: dimn, indt(dimn)
        complex*16, intent(in)  :: sc(dimn, dimn)
        complex*16, intent(out) :: bc(dimn, dimn)

        real(8)              :: denr, deni
        complex*16          :: den

        bc(:, :) = 0.0d+00

        if (rank == 0) print *, 'G space Bmat iroot=', iroot

!  !$OMP parallel do schedule(dynamic,1) private(iu,ju,j,it,jt,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs
            iu = indt(i)
            ju = convert_active_to_global_idx(iu)

            Do j = i, dimn
                it = indt(j)
                jt = convert_active_to_global_idx(it)

!  B(u,t) = Siguma_w [eps(w)<0|EutEww|0>] + S(u,t)(-eps(t))

                Do iw = 1, nact
                    jw = convert_active_to_global_idx(iw)

                    Call dim2_density(iu, it, iw, iw, denr, deni)
                    den = DCMPLX(denr, deni)
                    bc(i, j) = bc(i, j) + den*eps(jw)

                End do

                bc(i, j) = bc(i, j) - sc(i, j)*eps(ju)

                bc(j, i) = DCONJG(bc(i, j))

            End do
        End do
!   !$OMP end parallel do
#ifdef HAVE_MPI
        call reduce_wrapper(mat=bc, root_rank=0)
#endif
        if (rank == 0) print *, 'bGmat is ended'
    End subroutine bGmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE vGmat_complex(nabi, iabi, v)
!
!  V(t,iab)   =  [SIGUMA_p:active <0|Etp|0>{(ai|bp)-(ap|bi)}]
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_file_manager, only: open_unformatted_file, check_iostat
        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)     :: nabi, iabi(nsec, nsec, ninact)

        complex*16, intent(out) :: v(nabi, nact)

        real(8)                  :: dr, di
        complex*16              :: cint2, dens

        integer :: i, j, k, l, tabi
        integer :: it, iostat, unit_int2
        logical :: is_end_of_file

        if (rank == 0) print *, 'Enter vGmat. Please ignore timer under this line.'
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
            End do

        end do
        close (unit_int2)

        if (rank == 0) print *, 'vGmat_ord is ended'
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
        if (rank == 0) print *, 'end allreduce vGmat'
#endif
    end subroutine vGmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_G_subspace_real()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx, convert_secondary_to_global_idx

        Implicit NONE
#ifdef HAVE_MPI
        include 'mpif.h'
#endif

        integer :: dimn, dimm, dammy
        real(8), allocatable  :: wsnew(:), ws(:), wb(:)
        real(8)               :: e2(2*nsymrpa), alpha, e
        real(8), allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
        real(8), allocatable  :: bc(:, :)
        real(8), allocatable  :: bc0(:, :), bc1(:, :), v(:, :), vc(:), vc1(:)
        integer                  :: j, i, i0, syma, symb, isym, indt(1:nact)
        integer                  :: ia, it, ib, ii, ja, jt, jb, ji
        integer, allocatable     :: ia0(:), ib0(:), ii0(:), iabi(:, :, :)
        integer                  :: nabi

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

        e2 = 0.0d+00
        e2g = 0.0d+00
        dimn = 0
        indt = 0

        i0 = 0
        Do ia = 1, nsec
            ja = convert_secondary_to_global_idx(ia)
            Do ib = 1, ia - 1
                jb = convert_secondary_to_global_idx(ib)
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
            ja = convert_secondary_to_global_idx(ia)
            Do ib = 1, ia - 1
                jb = convert_secondary_to_global_idx(ib)
                Do ii = 1, ninact
                    ji = ii
                    i0 = i0 + 1
                    iabi(ia, ib, ii) = i0
                    iabi(ib, ia, ii) = i0
                    ia0(i0) = convert_secondary_to_global_idx(ia) ! secondary
                    ib0(i0) = convert_secondary_to_global_idx(ib) ! secondary
                    ii0(i0) = ii
                End do
            End do
        End do

        Allocate (v(nabi, nact))
        v = 0.0d+00
        if (rank == 0) print *, 'end before v matrices'
        Call vGmat_real(nabi, iabi, v)
        if (rank == 0) print *, 'end after vGmat'

        Do isym = 1, nsymrpa

            dimn = 0
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                if (irpamo(jt) == isym) then
                    dimn = dimn + 1
                    indt(dimn) = it
                End if
            End do

            if (rank == 0) print *, 'isym, dimn', isym, dimn
            If (dimn == 0) cycle ! Go to the next isym

            Allocate (sc(dimn, dimn))
            sc = 0.0d+00            ! sc N*N
            if (rank == 0) print *, 'before sGmat'
            Call sGmat_real(dimn, indt(1:dimn), sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'sG matrix is obtained normally'
            Allocate (ws(dimn))
            Allocate (sc0(dimn, dimn))
            sc0 = sc
            if (rank == 0) print *, 'before cdiag'
            Call rdiag(sc, dimn, dimm, ws, smat_lin_dep_threshold)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'after s cdiag, new dimension is', dimm
            If (dimm == 0) then
                deallocate (sc0)
                deallocate (sc)
                deallocate (ws)
                cycle ! Go to the next isym
            End if

            Allocate (bc(dimn, dimn))                                 ! bc N*N
            bc = 0.0d+00
            if (rank == 0) print *, 'before bGmat'
            Call bGmat_real(dimn, sc0, indt(1:dimn), bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (rank == 0) print *, 'bC matrix is obtained normally'
            deallocate (sc0)

            if (rank == 0) print *, 'OK cdiag', dimn, dimm
            Allocate (uc(dimn, dimm))                                 ! uc N*M
            Allocate (wsnew(dimm))                                  ! wnew M
            uc(:, :) = 0.0d+00
            wsnew(:) = 0.0d+00
            if (rank == 0) print *, 'before ccutoff'
            Call rcutoff(sc, ws, dimn, dimm, smat_lin_dep_threshold, uc, wsnew)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'OK ccutoff'
            deallocate (ws)
            deallocate (sc)
            if (rank == 0) print *, 'before ulambda_s_half'
            Call ulambda_s_half(uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            deallocate (wsnew)

            if (rank == 0) print *, 'ucrams half OK'
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

            if (rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'
            Allocate (bc0(dimm, dimm))
            bc0 = bc1
            if (rank == 0) print *, 'before cdiag'
            Call rdiag(bc1, dimm, dammy, wb, bmat_no_cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'end cdiag'

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

                    vc1(1:dimm) = MATMUL(TRANSPOSE(uc(1:dimn, 1:dimm)), vc(1:dimn))
                    Deallocate (vc)

                    alpha = -eps(ji) + eps(ja) + eps(jb) - e0 + eshift  ! For Level Shift (2007/2/9)

                    vc1(1:dimm) = MATMUL(TRANSPOSE(bc1(1:dimm, 1:dimm)), vc1(1:dimm))

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
        End do

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
        if (rank == 0) print *, 'end solve_G_subspace'
    end subroutine solve_G_subspace_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE sGmat_real(dimn, indt, sc) ! Assume C1 molecule, overlap matrix S in space C

!  S(u,t) = <0|Eut|0>
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)      :: dimn, indt(dimn)
        real(8), intent(out)  :: sc(dimn, dimn)
        real(8)  ::a, b
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

                sc(i, j) = a
                sc(j, i) = a
            End do
        End do
!    !$OMP end parallel do
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=sc)
#endif

    End subroutine sGmat_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE bGmat_real(dimn, sc, indt, bc) ! Assume C1 molecule, overlap matrix B in space C
!
!
!  S(u,t) = <0|Eut|0>
!
!  B(u,t) = Siguma_w [eps(w)<0|EutEww|0>] + S(u,t)(-eps(t))
!
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer :: it, iu, iw, jt, ju, jw
        integer :: i, j

        integer, intent(in) :: dimn, indt(dimn)
        real(8), intent(in)  :: sc(dimn, dimn)
        real(8), intent(out) :: bc(dimn, dimn)

        real(8)              :: denr, deni
        real(8)          :: den

        bc(:, :) = 0.0d+00

        if (rank == 0) print *, 'G space Bmat iroot=', iroot

!  !$OMP parallel do schedule(dynamic,1) private(iu,ju,j,it,jt,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs
            iu = indt(i)
            ju = convert_active_to_global_idx(iu)

            Do j = i, dimn
                it = indt(j)
                jt = convert_active_to_global_idx(it)

!  B(u,t) = Siguma_w [eps(w)<0|EutEww|0>] + S(u,t)(-eps(t))

                Do iw = 1, nact
                    jw = convert_active_to_global_idx(iw)

                    Call dim2_density(iu, it, iw, iw, denr, deni)
                    den = denr
                    bc(i, j) = bc(i, j) + den*eps(jw)

                End do

                bc(i, j) = bc(i, j) - sc(i, j)*eps(ju)

                bc(j, i) = bc(i, j)

            End do
        End do
!   !$OMP end parallel do
#ifdef HAVE_MPI
        call reduce_wrapper(mat=bc, root_rank=0)
#endif
        if (rank == 0) print *, 'bGmat is ended'
    End subroutine bGmat_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE vGmat_real(nabi, iabi, v)
!
!  V(t,iab)   =  [SIGUMA_p:active <0|Etp|0>{(ai|bp)-(ap|bi)}]
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_file_manager, only: open_unformatted_file, check_iostat
        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)     :: nabi, iabi(nsec, nsec, ninact)

        real(8), intent(out) :: v(nabi, nact)

        real(8)                  :: dr, di
        real(8)              :: cint2, dens

        integer :: i, j, k, l, tabi
        integer :: it, iostat, unit_int2
        logical :: is_end_of_file

        if (rank == 0) print *, 'Enter vGmat. Please ignore timer under this line.'
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
                dens = dr
                v(tabi, it) = v(tabi, it) + cint2*dens
            End do

        end do
        close (unit_int2)

        if (rank == 0) print *, 'vGmat_ord is ended'
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
        if (rank == 0) print *, 'end allreduce vGmat'
#endif
    end subroutine vGmat_real

end subroutine solve_G_subspace
