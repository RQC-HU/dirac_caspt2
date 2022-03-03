! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

  SUBROUTINE prtoutfock  ! TO PRINT OUT FOCK MATRIX

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
      !   integer :: i, j, ii, jj
        integer :: i, j


! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

         write(*,*)'inactive-inactive'

         do i = 1, ninact
         do j = 1, ninact
            if((i/=j).and.(ABS(f(i,j))>1.0d-10))then
               write(*,'(2I4,3E20.10)')i,j,f(i,j),oner(i,j)
            end if
         end do
         end do

         write(*,*)'inactive-active'

         do i = 1, ninact
         do j = ninact+1, ninact+nact
            if((i/=j).and.(ABS(f(i,j))>1.0d-10))then
               write(*,'(2I4,3E20.10)')i,j,f(i,j),oner(i,j)
            end if
         end do
         end do

         write(*,*)'inactive-secondary'

         do i = 1, ninact
         do j = ninact+nact+1, ninact+nact+nsec
            if((i/=j).and.(ABS(f(i,j))>1.0d-10))then
               write(*,'(2I4,3E20.10)')i,j,f(i,j),oner(i,j)
            end if
         end do
         end do

         write(*,*)'active-active'

         do i = ninact+1, ninact+nact
         do j = ninact+1, ninact+nact
            if((i/=j).and.(ABS(f(i,j))>1.0d-10))then
               write(*,'(2I4,3E20.10)')i,j,f(i,j),oner(i,j)
            end if
         end do
         end do

         write(*,*)'active-secondary'

         do i = ninact+1, ninact+nact
         do j =ninact+nact+1, ninact+nact+nsec
            if((i/=j).and.(ABS(f(i,j))>1.0d-10))then
               write(*,'(2I4,3E20.10)')i,j,f(i,j),oner(i,j)
            end if
         end do
         end do

         write(*,*)'secondary-secondary'

         do i = ninact+nact+1, ninact+nact+nsec
         do j = ninact+nact+1, ninact+nact+nsec
            if((i/=j).and.(ABS(f(i,j))>1.0d-10))then
               write(*,'(2I4,3E20.10)')i,j,f(i,j),oner(i,j)
            end if
         end do
         end do

      end SUBROUTINE prtoutfock
