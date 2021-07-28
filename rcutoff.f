! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE rcutoff (sr, w, dimn, dimm, thres, ur, wnew)
                           ! diagonalization of real symmetric matrix
                           !  and remove linear dependency for any S matrix

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
     
        integer, intent(in)  :: dimn, dimm
        real*8, intent(in)   :: thres, sr(dimn,dimn), w(dimn)
        real*8, intent(out)  :: ur(dimn,dimm), wnew(dimm)
        integer :: j0, j, i, i0, i1
        integer :: k0, l0, ii, jj, kk, ll

        write(*,*) 'New dimension becomes ', dimm

        j0 = 0
        do i0 = 1, dimn
           if( w(i0) >= thres ) then
              j0 = j0+1
              ur(:,j0) = sr(:,i0)               
              wnew(j0) = w(i0)
           end if
        end do

!test

        write(*,*) 'Eigenvalue and eigen vector becomes'
        do i0 = 1, dimm
           write(*,*) i0, 'th state'
           write(*,*) wnew(i0)
!           write(*,*) ur(:,i0)
        end do

 1000   continue
   end subroutine rcutoff



