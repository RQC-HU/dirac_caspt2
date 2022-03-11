! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE rcutoff(sr, w, dimn, dimm, ur, wnew)
! diagonalization of real symmetric matrix
!  and remove linear dependency for any S matrix

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE

    integer, intent(in)  :: dimn, dimm
    real*8, intent(in)   :: sr(dimn, dimn), w(dimn)
    real*8, intent(out)  :: ur(dimn, dimm), wnew(dimm)
    integer :: j0, j, i, i0, i1
    integer :: k0, l0, ii, jj, kk, ll

    write (*, *) 'New dimension becomes ', dimm

    j0 = 0
    do i0 = 1, dimn
        if (w(i0) >= thres) then
            j0 = j0 + 1
            ur(:, j0) = sr(:, i0)
            wnew(j0) = w(i0)
        end if
    end do

!test

    write (*, *) 'Eigenvalue and eigen vector becomes'
    do i0 = 1, dimm
        write (*, *) i0, 'th state'
        write (*, *) wnew(i0)
!           write(*,*) ur(:,i0)
    end do

1000 continue
end subroutine rcutoff

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE ccutoff(sc, w, dimn, dimm, uc, wnew)

!      Call ccutoff (sc, ws, dimn, dimm, thres, uc, wsnew)

! diagonalization of real symmetric matrix
!  and remove linear dependency for any S matrix

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE

    integer, intent(in)      :: dimn, dimm
    real*8, intent(in)       :: w(dimn)
    real*8, intent(out)      :: wnew(dimm)

    complex*16, intent(in)   :: sc(dimn, dimn)
    complex*16, intent(out)  :: uc(dimn, dimm)
    integer :: j0, j, i, i0, i1
    integer :: k0, l0, ii, jj, kk, ll
    if (rank == 0) then ! Process limits for output
        write (*, *) 'New dimension becomes ', dimm
    end if
    uc = 0.0d+00
    wnew = 0.0d+00

    j0 = 0
    do i0 = 1, dimn
        if (ABS(w(i0)) >= thres) then
            j0 = j0 + 1
            uc(:, j0) = sc(:, i0)
            wnew(j0) = w(i0)
        end if
    end do

!test

!!!        write(*,*) 'Eigenvalue and eigen vector becomes'
!!!        do i0 = 1, dimm
!!!           write(*,*) i0, 'th state'
!!!           write(*,*) wnew(i0)
!!!!           write(*,*) uc(:,i0)
!!!        end do

1000 continue
end subroutine ccutoff
