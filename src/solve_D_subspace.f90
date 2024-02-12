SUBROUTINE solve_D_subspace(e0, e2d)

    use module_global_variables
    use module_realonly, only: realonly
    use module_time
    implicit none
    real(8), intent(in) :: e0
    real(8), intent(out):: e2d

    if (realonly%is_realonly()) then
        call solve_D_subspace_real()
    else
        call solve_D_subspace_complex()
    end if

contains
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_D_subspace_complex()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx, convert_secondary_to_global_idx
        use module_ulambda_s_half, only: ulambda_s_half

        Implicit NONE

        integer :: dimn, dimm, dammy
        integer, allocatable :: indsym(:, :)
        real(8), allocatable  :: wsnew(:), ws(:), wb(:)
        real(8)               :: e2(nsymrpa*2), e, alpha
        complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
        complex*16, allocatable  :: bc(:, :)
        complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :), vc(:), vc1(:)
        integer, allocatable     :: ia0(:), ii0(:), iai(:, :)
        integer                  :: nai
        integer :: j, i, syma, isym, i0
        integer :: ia, it, ii, iu
        integer :: ja, jt, ji, ju

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

        if (rank == 0) print *, 'ENTER solve D part'

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
                ia0(i0) = convert_secondary_to_global_idx(ia) ! Secondary
                ii0(i0) = ii ! inactive
            End do
        End do
        Allocate (v(nai, nact, nact))
        v = 0.0d+00
        Call vDmat_complex(nai, iai, v)

        Do isym = 1, nsymrpa

            dimn = 0
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Do iu = 1, nact
                    ju = convert_active_to_global_idx(iu)
                    if (nsymrpa /= 1) syma = MULTB_D(irpamo(jt), irpamo(ju))
                    if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                        dimn = dimn + 1
                    End if
                End do
            End do

            if (rank == 0) print *, 'isym, dimn', isym, dimn

            If (dimn == 0) cycle ! Go to the next isym

            Allocate (indsym(2, dimn))
            indsym = 0
            dimn = 0

            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Do iu = 1, nact
                    ju = convert_active_to_global_idx(iu)

                    if (nsymrpa /= 1) syma = MULTB_D(irpamo(jt), irpamo(ju))

                    if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                        dimn = dimn + 1
                        indsym(1, dimn) = it
                        indsym(2, dimn) = iu
                    End if
                End do
            End do

            Allocate (sc(dimn, dimn))
            sc = 0.0d+00            ! sc N*N
            Call sDmat_complex(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            Allocate (ws(dimn))
            ws = 0.0d+00

            Allocate (sc0(dimn, dimn))
            sc0 = sc
            Call cdiag(sc, dimn, dimm, ws, smat_lin_dep_threshold)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'after D subspace S matrix cdiag, new dimension is', dimm

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
                if (rank == 0) print *, 'Check whether U*SU is diagonal END'
            End if

            Allocate (bc(dimn, dimn))                                 ! bc N*N
            bc = 0.0d+00
            Call bDmat_complex(dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            deallocate (sc0)

            Allocate (uc(dimn, dimm))                                 ! uc N*M
            Allocate (wsnew(dimm))                                  ! wnew M
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

            if (rank == 0) then
                IF (debug) then

                    print *, 'Check whether bc1 is hermite or not'
                    Do i = 1, dimm
                        Do j = i, dimm
                            if (ABS(bc1(i, j) - DCONJG(bc1(j, i))) > 1.0d-6) then
                                print '(2I4,2E15.7)', i, j, bc1(i, j) - bc1(j, i)
                            End if
                        End do
                    End do
                    print *, 'Check whether bc1 is hermite or not END'

                End if
            end if

            deallocate (bc)
            deallocate (bc0)

            Allocate (wb(dimm))
            wb = 0.0d+00

            if (rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'
            Allocate (bc0(dimm, dimm))
            bc0 = bc1
            Call cdiag(bc1, dimm, dammy, wb, bmat_no_cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            If (debug) then

                if (rank == 0) print *, 'Check whether bc is really diagonalized or not'
                Call checkdgc(dimm, bc0, bc1, wb)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (rank == 0) print *, 'Check whether bc is really diagonalized or not END'
            End if

            deallocate (bc0)

            if (rank == 0) print *, 'bC1 matrix is diagonalized!'
            e2 = 0.0d+00
            Do i0 = 1, nai
                ja = ia0(i0)
                ji = ii0(i0)
                if (nsymrpa /= 1) then
                    syma = MULTB_D(irpamo(ja), irpamo(ji))
                    syma = MULTB_S(syma, isym)
                end if
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

            End do

            deallocate (indsym)
            deallocate (uc)
            deallocate (wb)
            Deallocate (bc1)

            if (rank == 0) print '("e2d(",I3,") = ",E20.10," a.u.")', isym, e2(isym)
            e2d = e2d + e2(isym)
        End do

        if (rank == 0) then
            print '("e2d      = ",E20.10," a.u.")', e2d

            print '("sumc2,d  = ",E20.10)', sumc2local
        end if
        sumc2 = sumc2 + sumc2local

        deallocate (iai)
        deallocate (ia0)
        deallocate (ii0)
        deallocate (v)

        if (rank == 0) print *, 'end solve_D_subspace'
    end subroutine solve_D_subspace_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE sDmat_complex(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space D

!  S(xy, tu) = <0|EyxEtu|0>
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)      :: dimn, indsym(2, dimn)
        complex*16, intent(out)  :: sc(dimn, dimn)
        real(8)  :: a, b
        integer :: it, iu, iy, ix
        integer :: i, j

        if (rank == 0) print *, 'Start D subspace S matrix'
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

            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=sc)
#endif
        if (rank == 0) print *, 'D subspace S matrix is obtained normally'
    End subroutine sDmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE bDmat_complex(dimn, sc, indsym, bc) ! Assume C1 molecule, overlap matrix B in space D
!
!  B(xy,tu) = Siguma_w [eps(w)<0|EyxEtuEww|0>+S(xy,tu)(eps(t)-eps(u))]
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in) :: dimn, indsym(2, dimn)
        complex*16, intent(in)  :: sc(dimn, dimn)
        complex*16, intent(out) :: bc(dimn, dimn)

        real(8)              :: e, denr, deni
        complex*16          :: den

        integer :: it, iu, ix, iy, iw
        integer :: jt, ju, jy, jx, jw, i, j

        bc(:, :) = 0.0d+00

        if (rank == 0) print *, 'Start D subspace B matrix'
!$OMP parallel do schedule(dynamic,1) private(ix,iy,jx,jy,it,iu,jt,ju,e,j,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs

            ix = indsym(1, i)
            jx = convert_active_to_global_idx(ix)
            iy = indsym(2, i)
            jy = convert_active_to_global_idx(iy)

            Do j = i, dimn

                it = indsym(1, j)
                jt = convert_active_to_global_idx(it)
                iu = indsym(2, j)
                ju = convert_active_to_global_idx(iu)

!  B(xy,tu) = Siguma_w [eps(w)<0|EyxEtuEww|0>+S(xy,tu)(eps(t)-eps(u))]

                e = eps(jt) - eps(ju)

                Do iw = 1, nact
                    jw = convert_active_to_global_idx(iw)

                    Call dim3_density(iy, ix, it, iu, iw, iw, denr, deni)
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
        if (rank == 0) print *, 'D subspace B matrix is obtained normally'
    End subroutine bDmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
    SUBROUTINE vDmat_complex(nai, iai, v)
!
!
! V(a,i) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
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

        integer, intent(in)     :: nai, iai(nsec, ninact)
        complex*16, intent(out) :: v(nai, nact, nact)
        real(8)                  :: dr, di
        complex*16              :: cint1
        complex*16              :: cint2, d
        complex*16              :: effh(nsec, ninact)
        integer :: i, j, k, l, tai, iostat, unit_int2
        integer :: it, jt, ju, iu, ia, ii, ja, ji
        logical :: is_end_of_file

        if (rank == 0) print *, 'Start D subspace V matrix'
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
            ja = convert_secondary_to_global_idx(ia)
            Do ii = 1, ninact
                ji = ii
                Call tramo1(ja, ji, cint1)
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

        call open_unformatted_file(unit=unit_int2, file=d1int, status='old', optional_action='read')
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
            call check_iostat(iostat=iostat, file=d1int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

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
                jt = convert_active_to_global_idx(it)
                Do iu = 1, nact
                    ju = convert_active_to_global_idx(iu)

                    Call dim2_density(iu, it, k, l, dr, di)
                    d = DCMPLX(dr, di)
                    v(tai, it, iu) = v(tai, it, iu) + d*cint2

                End do
            End do
!$OMP end parallel do
        end do
        close (unit_int2)

        if (rank == 0) print *, 'reading D1int2 is over'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! V(a,i, jt, ju) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2, file=d2int, status='old', optional_action='read')
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
            call check_iostat(iostat=iostat, file=d2int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            ja = i
            ji = l
            tai = iai(ja, ji)
!$OMP parallel do schedule(dynamic,1) private(it,ju,dr,di,d)
            Do it = 1, nact
                Do iu = 1, nact

                    Call dim2_density(iu, it, k, j, dr, di)
                    d = DCMPLX(dr, di)

                    v(tai, it, iu) = v(tai, it, iu) - d*cint2

                End do
            End do
!$OMP end parallel do
        end do
        close (unit_int2)

        if (rank == 0) print *, 'reading D2int2 is over'

        call open_unformatted_file(unit=unit_int2, file=d3int, status='old', optional_action='read') ! (ai|jk) is stored
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
            call check_iostat(iostat=iostat, file=d3int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            if (j /= k .and. k == l) then !(ai|kk)

                effh(i, j) = effh(i, j) + cint2

            elseif (j == k .and. k /= l) then !(ak|ki)

                effh(i, l) = effh(i, l) - cint2

            end if
        end do
        close (unit_int2)

        if (rank == 0) print *, 'reading D3int2 is over'

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=effh)
#endif

!$OMP parallel do schedule(dynamic,1) private(ia,ja,ii,ji,tai,it,jt,iu,ju,dr,di,d)
        Do ia = rank + 1, nsec, nprocs
            Do ii = 1, ninact
                tai = iai(ia, ii)
                Do it = 1, nact
                    Do iu = 1, nact
                        Call dim1_density(iu, it, dr, di)

                        d = DCMPLX(dr, di)
                        v(tai, it, iu) = v(tai, it, iu) + effh(ia, ii)*d
                    End do
                End do

            End do
        End do
!$OMP end parallel do

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
#endif
        if (rank == 0) print *, 'D subspace V matrix is obtained normally'
    end subroutine vDmat_complex

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE solve_D_subspace_real()

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx, convert_secondary_to_global_idx
        use module_ulambda_s_half, only: ulambda_s_half

        Implicit NONE

        integer :: dimn, dimm, dammy
        integer, allocatable :: indsym(:, :)
        real(8), allocatable  :: wsnew(:), ws(:), wb(:)
        real(8)               :: e2(nsymrpa*2), e, alpha
        real(8), allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
        real(8), allocatable  :: bc(:, :)
        real(8), allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :), vc(:), vc1(:)
        integer, allocatable     :: ia0(:), ii0(:), iai(:, :)
        integer                  :: nai
        integer :: j, i, syma, isym, i0
        integer :: ia, it, ii, iu
        integer :: ja, jt, ji, ju

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

        if (rank == 0) print *, 'ENTER solve D part'

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
                ia0(i0) = convert_secondary_to_global_idx(ia) ! Secondary
                ii0(i0) = ii ! inactive
            End do
        End do
        Allocate (v(nai, nact, nact))
        v = 0.0d+00
        Call vDmat_real(nai, iai, v)

        Do isym = 1, nsymrpa

            dimn = 0
            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Do iu = 1, nact
                    ju = convert_active_to_global_idx(iu)
                    if (nsymrpa /= 1) syma = MULTB_D(irpamo(jt), irpamo(ju))
                    if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                        dimn = dimn + 1
                    End if
                End do
            End do

            if (rank == 0) print *, 'isym, dimn', isym, dimn

            If (dimn == 0) cycle ! Go to the next isym

            Allocate (indsym(2, dimn))
            indsym = 0
            dimn = 0

            Do it = 1, nact
                jt = convert_active_to_global_idx(it)
                Do iu = 1, nact
                    ju = convert_active_to_global_idx(iu)

                    if (nsymrpa /= 1) syma = MULTB_D(irpamo(jt), irpamo(ju))

                    if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                        dimn = dimn + 1
                        indsym(1, dimn) = it
                        indsym(2, dimn) = iu
                    End if
                End do
            End do

            Allocate (sc(dimn, dimn))
            sc = 0.0d+00            ! sc N*N
            Call sDmat_real(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            Allocate (ws(dimn))
            ws = 0.0d+00

            Allocate (sc0(dimn, dimn))
            sc0 = sc
            Call rdiag(sc, dimn, dimm, ws, smat_lin_dep_threshold)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (rank == 0) print *, 'after D subspace S matrix rdiag, new dimension is', dimm

            If (dimm == 0) then
                deallocate (indsym)
                deallocate (sc0)
                deallocate (sc)
                deallocate (ws)
                cycle ! Go to the next isym
            End if

            Allocate (bc(dimn, dimn))                                 ! bc N*N
            bc = 0.0d+00
            Call bDmat_real(dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            deallocate (sc0)

            Allocate (uc(dimn, dimm))                                 ! uc N*M
            Allocate (wsnew(dimm))                                  ! wnew M
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

            if (rank == 0) then
                IF (debug) then

                    print *, 'Check whether bc1 is hermite or not'
                    Do i = 1, dimm
                        Do j = i, dimm
                            if (ABS(bc1(i, j) - bc1(j, i)) > 1.0d-6) then
                                print '(2I4,2E15.7)', i, j, bc1(i, j) - bc1(j, i)
                            End if
                        End do
                    End do
                    print *, 'Check whether bc1 is hermite or not END'

                End if
            end if

            deallocate (bc)
            deallocate (bc0)

            Allocate (wb(dimm))
            wb = 0.0d+00

            if (rank == 0) print *, 'bC matrix is transrated to bc1(M*M matrix)!'
            Allocate (bc0(dimm, dimm))
            bc0 = bc1
            Call rdiag(bc1, dimm, dammy, wb, bmat_no_cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            deallocate (bc0)

            if (rank == 0) print *, 'bC1 matrix is diagonalized!'
            e2 = 0.0d+00
            Do i0 = 1, nai
                ja = ia0(i0)
                ji = ii0(i0)
                if (nsymrpa /= 1) then
                    syma = MULTB_D(irpamo(ja), irpamo(ji))
                    syma = MULTB_S(syma, isym)
                end if
                If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then

                    Allocate (vc(dimn))
                    Do it = 1, dimn
                        vc(it) = v(i0, indsym(1, it), indsym(2, it))
                    End do

                    Allocate (vc1(dimm))

                    vc1(1:dimm) = MATMUL(TRANSPOSE(uc(1:dimn, 1:dimm)), vc(1:dimn))
                    Deallocate (vc)

                    alpha = +eps(ja) - eps(ji) - e0 + eshift  ! For Level Shift (2007/2/9)

                    vc1(1:dimm) = MATMUL(TRANSPOSE(bc1(1:dimm, 1:dimm)), vc1(1:dimm))

                    Do j = 1, dimm
                        sumc2local = sumc2local + (ABS(vc1(j))**2.0d+00)/((alpha + wb(j))**2.0d+00)
                        e = (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                        e2(isym) = e2(isym) - e
                    End do

                    deallocate (vc1)

                End if

            End do

            deallocate (indsym)
            deallocate (uc)
            deallocate (wb)
            Deallocate (bc1)

            if (rank == 0) print '("e2d(",I3,") = ",E20.10," a.u.")', isym, e2(isym)
            e2d = e2d + e2(isym)
        End do

        if (rank == 0) then
            print '("e2d      = ",E20.10," a.u.")', e2d

            print '("sumc2,d  = ",E20.10)', sumc2local
        end if
        sumc2 = sumc2 + sumc2local

        deallocate (iai)
        deallocate (ia0)
        deallocate (ii0)
        deallocate (v)

        if (rank == 0) print *, 'end solve_D_subspace'
    end subroutine solve_D_subspace_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE sDmat_real(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space D

!  S(xy, tu) = <0|EyxEtu|0>
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in)      :: dimn, indsym(2, dimn)
        real(8), intent(out)  :: sc(dimn, dimn)
        real(8)  :: a, b
        integer :: it, iu, iy, ix
        integer :: i, j

        if (rank == 0) print *, 'Start D subspace S matrix'
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
                sc(i, j) = a
                sc(j, i) = sc(i, j)

            End do
        End do
!$OMP end parallel do
#ifdef HAVE_MPI
        call allreduce_wrapper(mat=sc)
#endif
        if (rank == 0) print *, 'D subspace S matrix is obtained normally'
    End subroutine sDmat_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    SUBROUTINE bDmat_real(dimn, sc, indsym, bc) ! Assume C1 molecule, overlap matrix B in space D
!
!  B(xy,tu) = Siguma_w [eps(w)<0|EyxEtuEww|0>+S(xy,tu)(eps(t)-eps(u))]
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        use module_global_variables
        use module_index_utils, only: convert_active_to_global_idx
#ifdef HAVE_MPI
        use module_mpi
#endif
        Implicit NONE

        integer, intent(in) :: dimn, indsym(2, dimn)
        real(8), intent(in)  :: sc(dimn, dimn)
        real(8), intent(out) :: bc(dimn, dimn)

        real(8)              :: e, denr, deni
        real(8)          :: den

        integer :: it, iu, ix, iy, iw
        integer :: jt, ju, jy, jx, jw, i, j

        bc(:, :) = 0.0d+00

        if (rank == 0) print *, 'Start D subspace B matrix'
!$OMP parallel do schedule(dynamic,1) private(ix,iy,jx,jy,it,iu,jt,ju,e,j,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs

            ix = indsym(1, i)
            jx = convert_active_to_global_idx(ix)
            iy = indsym(2, i)
            jy = convert_active_to_global_idx(iy)

            Do j = i, dimn

                it = indsym(1, j)
                jt = convert_active_to_global_idx(it)
                iu = indsym(2, j)
                ju = convert_active_to_global_idx(iu)

!  B(xy,tu) = Siguma_w [eps(w)<0|EyxEtuEww|0>+S(xy,tu)(eps(t)-eps(u))]

                e = eps(jt) - eps(ju)

                Do iw = 1, nact
                    jw = convert_active_to_global_idx(iw)

                    Call dim3_density(iy, ix, it, iu, iw, iw, denr, deni)
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
        if (rank == 0) print *, 'D subspace B matrix is obtained normally'
    End subroutine bDmat_real

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
    SUBROUTINE vDmat_real(nai, iai, v)
!
!
! V(a,i) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
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

        integer, intent(in)     :: nai, iai(nsec, ninact)
        real(8), intent(out) :: v(nai, nact, nact)
        real(8)                  :: dr, di
        complex*16              :: cint1
        real(8)              :: cint2, d
        real(8)              :: effh(nsec, ninact)
        integer :: i, j, k, l, tai, iostat, unit_int2
        integer :: it, jt, ju, iu, ia, ii, ja, ji
        logical :: is_end_of_file

        if (rank == 0) print *, 'Start D subspace V matrix'
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
            ja = convert_secondary_to_global_idx(ia)
            Do ii = 1, ninact
                ji = ii
                Call tramo1(ja, ji, cint1)
                effh(ia, ii) = real(cint1, kind=KIND(effh))
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

        call open_unformatted_file(unit=unit_int2, file=d1int, status='old', optional_action='read')
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
            call check_iostat(iostat=iostat, file=d1int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

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
                jt = convert_active_to_global_idx(it)
                Do iu = 1, nact
                    ju = convert_active_to_global_idx(iu)

                    Call dim2_density(iu, it, k, l, dr, di)
                    d = dr
                    v(tai, it, iu) = v(tai, it, iu) + d*cint2

                End do
            End do
!$OMP end parallel do
        end do
        close (unit_int2)

        if (rank == 0) print *, 'reading D1int2 is over'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! V(a,i, jt, ju) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call open_unformatted_file(unit=unit_int2, file=d2int, status='old', optional_action='read')
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
            call check_iostat(iostat=iostat, file=d2int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            ja = i
            ji = l
            tai = iai(ja, ji)
!$OMP parallel do schedule(dynamic,1) private(it,ju,dr,di,d)
            Do it = 1, nact
                Do iu = 1, nact

                    Call dim2_density(iu, it, k, j, dr, di)
                    d = dr

                    v(tai, it, iu) = v(tai, it, iu) - d*cint2

                End do
            End do
!$OMP end parallel do
        end do
        close (unit_int2)

        if (rank == 0) print *, 'reading D2int2 is over'

        call open_unformatted_file(unit=unit_int2, file=d3int, status='old', optional_action='read') ! (ai|jk) is stored
        do
            read (unit_int2, iostat=iostat) i, j, k, l, cint2 !  (ij|kl)
            call check_iostat(iostat=iostat, file=d3int, end_of_file_reached=is_end_of_file)
            if (is_end_of_file) then
                exit
            end if

            if (j /= k .and. k == l) then !(ai|kk)

                effh(i, j) = effh(i, j) + cint2

            elseif (j == k .and. k /= l) then !(ak|ki)

                effh(i, l) = effh(i, l) - cint2

            end if
        end do
        close (unit_int2)

        if (rank == 0) print *, 'reading D3int2 is over'

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=effh)
#endif

!$OMP parallel do schedule(dynamic,1) private(ia,ja,ii,ji,tai,it,jt,iu,ju,dr,di,d)
        Do ia = rank + 1, nsec, nprocs
            Do ii = 1, ninact
                tai = iai(ia, ii)
                Do it = 1, nact
                    Do iu = 1, nact
                        Call dim1_density(iu, it, dr, di)

                        d = dr
                        v(tai, it, iu) = v(tai, it, iu) + effh(ia, ii)*d
                    End do
                End do

            End do
        End do
!$OMP end parallel do

#ifdef HAVE_MPI
        call allreduce_wrapper(mat=v)
#endif
        if (rank == 0) print *, 'D subspace V matrix is obtained normally'
    end subroutine vDmat_real

end subroutine solve_D_subspace
