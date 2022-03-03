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

    integer :: dimn, dimm, count, dammy

    integer, allocatable :: indsym(:, :)

    real*8, allocatable  :: sr(:, :), ur(:, :)
    real*8, allocatable  :: br(:, :), wsnew(:), ws(:), wb(:)
    real*8, allocatable  :: br0(:, :), br1(:, :)
    real*8               :: e2(2*nsymrpa), alpha, e

    complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
    complex*16, allocatable  :: bc(:, :)
    complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :), vc(:), vc1(:)

    logical :: cutoff
    integer :: j, i, k, syma, symb, isym, indt(1:nact)
    integer :: ia, it, ij, ii, ja, jt, jj, ji

    integer :: i0
    integer, allocatable     :: ia0(:), ii0(:), ij0(:), iaij(:, :, :)
    integer                  :: naij

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
    !        thresd = thres
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
    if (rank == 0) then ! Process limits for output
        write (*, *) ' ENTER solv E part'
        write (*, *) ' nsymrpa', nsymrpa
    end if
    i0 = 0
    Do ia = 1, nsec
        Do ii = 1, ninact
            Do ij = 1, ii - 1                ! i > j
                i0 = i0 + 1
            End do
        End do
    End do

    naij = i0
    Allocate (iaij(ninact + nact + 1:ninact + nact + nsec, 1:ninact, 1:ninact))
    iaij = 0
    Allocate (ia0(naij))
    Allocate (ii0(naij))
    Allocate (ij0(naij))

    i0 = 0
    Do ia = 1, nsec
        ja = ia + ninact + nact
        Do ii = 1, ninact
            ji = ii
            Do ij = 1, ii - 1                ! i > j
                jj = ij
                i0 = i0 + 1
                iaij(ja, ji, jj) = i0
                iaij(ja, jj, ji) = i0
                ia0(i0) = ja
                ii0(i0) = ji
                ij0(i0) = jj
            End do
        End do
    End do

    Allocate (v(naij, ninact + 1:ninact + nact))
    v = 0.0d+00
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (*, *) 'end before v matrices'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Call vEmat_ord_ty(naij, iaij, v)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'come'
    end if
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

        if (rank == 0) then ! Process limits for output
            write (*, *) 'isym, dimn', isym, dimn
        end if
        If (dimn == 0) goto 1000

        Allocate (sc(dimn, dimn))
        sc = 0.0d+00            ! sc N*N
        if (rank == 0) then ! Process limits for output
            write (*, *) 'before sEmat'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call sEmat(dimn, indt(1:dimn), sc)
        !      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (rank == 0) then ! Process limits for output
            write (*, *) 'sc matrix is obtained normally'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Allocate (ws(dimn))

        cutoff = .TRUE.
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
            write (*, *) 'after s cdiag, new dimension is', dimm
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        If (dimm == 0) then
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
            write (*, *) 'before bEmat'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
        Call bEmat(e0, dimn, sc0, indt(1:dimn), bc)
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

        Do i0 = 1, naij
            ja = ia0(i0)
            ji = ii0(i0)
            jj = ij0(i0)

            !     EtiEaj|0>

            syma = MULTB_D(irpmo(ja), irpmo(jj))
            symb = MULTB_D(isym, irpmo(ji))
            syma = MULTB_S(symb, syma)

            If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then

                Allocate (vc(dimn))

                Do it = 1, dimn
                    vc(it) = v(i0, indt(it) + ninact)
                End do

                Allocate (vc1(dimm))
                vc1 = 0.0d+00

                vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn, 1:dimm))), vc(1:dimn))
                Deallocate (vc)

                alpha = +eps(ja) - eps(ji) - eps(jj) - e0 + eshift  ! For Level Shift (2007/2/9)

                vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm, 1:dimm))), vc1(1:dimm))

                Do j = 1, dimm
                    e = DCONJG(vc1(j))*vc1(j)/(alpha + wb(j))
                    sumc2local = sumc2local + e/(alpha + wb(j))
                    e2(isym) = e2(isym) - e
                End do

                deallocate (vc1)

            End if

        End do

        deallocate (uc)
        deallocate (wb)
        Deallocate (bc1)

1000    if (rank == 0) write (*, '("e2e(",I3,") = ",E20.10,"a.u.")') isym, e2(isym)
        e2e = e2e + e2(isym)
        if (rank == 0) then ! Process limits for output
            write (*, *) 'End e2(isym) add'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
    End do                  ! isym

    if (rank == 0) then ! Process limits for output
        write (*, '("e2e      = ",E20.10,"a.u.")') e2e

        write (*, '("sumc2,e  = ",E20.10)') sumc2local
    end if
    sumc2 = sumc2 + sumc2local

    deallocate (iaij)
    deallocate (ia0)
    deallocate (ii0)
    deallocate (ij0)
    deallocate (v)

    continue
    if (rank == 0) then ! Process limits for output
        write (*, *) 'end solveE_ord_ty'
    end if
end

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE sEmat(dimn, indt, sc) ! Assume C1 molecule, overlap matrix S in space E

    ! S(u,t) = d(ut) - <0|Etu|0>
    !
    ! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    ! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
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

            !              write(*,*)i,j,sc(i,j)

        End do               !j
    End do                  !i
    !$OMP end parallel do
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, sc(1, 1), dimn**2, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
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

    if (rank == 0) then ! Process limits for output
        write (*, *) 'E space Bmat iroot=', iroot
    end if

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

            !              write(*,*)'bc',i,j, bc(i,j)

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
        write (*, *) 'bEmat is ended'
    end if
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

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif

    integer, intent(in)     :: naij, iaij(ninact + nact + 1:ninact + nact + nsec, 1:ninact, 1:ninact)

    complex*16, intent(out) :: v(naij, ninact + 1:ninact + nact)

    real*8                  :: dr, di
    complex*16              :: cint2, dens

    integer :: i, j, k, l, taij
    integer :: it, jt, ik
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1

    if (rank == 0) write (*, *) 'Enter vEmat. Please ignore timer under this line.'
    datetmp1 = date0; datetmp0 = date0
    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0
    v = 0.0d+00

    !  V(t,ija)   =[SIGUMA_p:active <0|Ept|0>{(ai|pj) - (aj|pi)}] - (ai|tj) + (aj|ti)   i > j

    !   open(1, file ='Eint', status='old', form='unformatted')  !  (31|21) stored
    ! open (1, file=eint, status='old', form='formatted')  !  (31|21) stored
    open (1, file=eint, status='old', form='unformatted')  !  (31|21) stored
    ! 30  read (1, '(4I4, 2e20.10)', err=10, end=20) i, j, k, l, cint2
30  read (1, err=10, end=20) i, j, k, l, cint2

    if (j == l) goto 30

    taij = iaij(i, j, l)
    ik = k - ninact

    !        write(*,*) i,j,k,l,taij,cint2

    if (j < l) then
        cint2 = -1.0d+00*cint2
    end if

    v(taij, k) = v(taij, k) - cint2

    !$OMP parallel do schedule(dynamic,1) private(jt,dr,di,dens)
    Do it = 1, nact
        jt = ninact + it
        Call dim1_density(it, ik, dr, di)          ! ik corresponds to p in above formula
        dens = DCMPLX(dr, di)
        v(taij, jt) = v(taij, jt) + cint2*dens
    End do                  ! it
    !$OMP end parallel do

    if (j < l) then
        cint2 = -1.0d+00*cint2          ! data cint2 becomes initial values!
    end if

   !! Take Kramers conjugate !
    !
    !        Call takekr( i, j, k, l, cint2)
    !
    !        taij = iaij(i, j, l)
    !        ik = k - ninact
    !
   !!        write(*,*) i,j,k,l,taij,cint2
    !
    !        if (j < l) then
    !           cint2 = -1.0d+00*cint2
    !        endif
    !
    !        v(taij,k) = v(taij, k) - cint2
    !
    !        Do it = 1, nact
    !           jt = ninact+it
    !           Call dim1_density (it, ik, dr, di)          ! ik corresponds to p in above formula
    !           dens = DCMPLX(dr, di)
    !           v(taij,jt) = v(taij, jt) + cint2*dens
    !        End do                  ! it

    goto 30

20  close (1); goto 100

10  write (*, *) 'error while opening file Eint'; goto 100

100 if (rank == 0) write (*, *) 'vEmat_ord_ty is ended'
    !  v(naij, ninact+1:ninact+nact)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, v(1, ninact + 1), naij*nact, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (*, *) 'end Allreduce vEmat'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
end subroutine vEmat_ord_ty
