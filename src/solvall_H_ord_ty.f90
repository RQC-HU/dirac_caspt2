! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE solvH_ord_ty(e0, e2h)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    real*8, intent(in) :: e0
    real*8, intent(out):: e2h
    Integer                :: ia, ib, ii, ij, syma, symb, i, j, k, l
    Integer                :: i0, j0, tab, nab, tij, nij, iostat
    Integer, allocatable    :: ia0(:), ib0(:), ii0(:), ij0(:), iab(:, :), iij(:, :)
    Complex*16             :: cint2
    Complex*16, allocatable :: v(:, :)
    Real*8                 :: e

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE H IS NOW CALCULATED
!
!     EaiEbj|0>         a > b, i > j
!
!   DRAS1 =-2   DRAS2 = 0    DRAS3 = +2
!
!
!  S(ckdl,aibj) = d(ac)d(bd)d(lj)d(ik)
!
!  (H0-E0S)ckdl,aibj = d(ac)d(bd)d(lj)d(ik)(eps(a)+eps(b)-eps(i)-eps(j)) = e(a,b,i,j)
!
!  V(aibj)   = (ai|bj) - (aj|bi)
!
! E2h = V(aibj)/e(a,b,i,j)

!        thresd = 1.0D-08
!        thres = 1.0D-08

    e2h = 0.0d+00
    e = 0.0d+00

    i0 = 0
    Do ia = ninact + nact + 1, ninact + nact + nsec
        Do ib = ninact + nact + 1, ia - 1
            i0 = i0 + 1
        End do
    End do

    nab = i0

    Allocate (iab(nsec, nsec))
    Allocate (ia0(nab))
    Allocate (ib0(nab))
    iab = 0

    i0 = 0
    Do ia = 1, nsec
        Do ib = 1, ia - 1
            i0 = i0 + 1
            iab(ia, ib) = i0
            iab(ib, ia) = i0
            ia0(i0) = ia + ninact + nact ! secondary
            ib0(i0) = ib + ninact + nact ! secondary
        End do
    End do

    i0 = 0
    Do ii = 1, ninact
        Do ij = 1, ii - 1
            i0 = i0 + 1
        End do
    End do

    nij = i0
    Allocate (iij(1:ninact, 1:ninact))
    Allocate (ii0(nij))
    Allocate (ij0(nij))
    iij = 0

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

    Allocate (v(nab, nij))
    v = 0.0d+00

    open (1, file=hint, status='old', form='unformatted')
    do
        read (1, iostat=iostat) i, j, k, l, cint2
        ! Exit the loop if the end of the file is reached
        if (iostat < 0) then
            if (rank == 0) print *, 'End of Eint'
            exit
        elseif (iostat > 0) then
            ! If iostat is greater than 0, error detected in the input file, so exit the program
            stop 'Error: Error in reading Eint'
        end if
        if (i <= k .or. j == l) cycle ! Read the next line if i <= k or j == l

        tab = iab(i, k)
        tij = iij(j, l)

        if (i > k .and. j > l) then
            v(tab, tij) = v(tab, tij) + cint2

        elseif (i > k .and. j < l) then
            v(tab, tij) = v(tab, tij) - cint2

        elseif (i < k .and. j > l) then               ! (kl|ij)  l > j + ; l < j -
            v(tab, tij) = v(tab, tij) - cint2

        elseif (i < k .and. j < l) then
            v(tab, tij) = v(tab, tij) + cint2

        end if

    end do

    close (1)
#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, v(1, 1), nab*nij, MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) print *, 'reading int2 is over'

    Do i0 = 1, nab
        ia = ia0(i0)
        ib = ib0(i0)

!     EaiEbj|0>         a > b, i > j

        Do j0 = 1, nij
            ii = ii0(j0)
            ij = ij0(j0)
            syma = MULTB_D(irpmo(ia), irpmo(ii))
            symb = MULTB_D(irpmo(ib), irpmo(ij))
            syma = MULTB_S(syma, symb)
            if (syma == 1) then

                e = eps(ia) + eps(ib) - eps(ii) - eps(ij) + eshift  ! For Level Shift (2007/2/9)

                coeff1 = v(i0, j0)/e
                sumc2local = sumc2local + ABS(coeff1)**2

                e2h = e2h - DCONJG(v(i0, j0))*v(i0, j0)/e
            end if
        End do
    End do

    if (rank == 0) then
        print '("e2h      = ",E20.10,"a.u.")', e2h
        print '("sumc2,h  = ",E20.10)', sumc2local
    end if
    sumc2 = sumc2 + sumc2local

    deallocate (v)
    deallocate (iab)
    deallocate (ia0)
    deallocate (ib0)
    deallocate (iij)
    deallocate (ii0)
    deallocate (ij0)

    if (rank == 0) print *, 'end solvH_ord_ty'
End SUBROUTINE solvH_ord_ty
