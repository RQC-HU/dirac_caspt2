! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE calce0(e0)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE

        real*8,     intent(out):: e0

        integer :: i, ii
        real*8  :: dr, di

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        e0 = 0.0d+00
        dr = 0.0d+00
        di = 0.0d+00
        write(*,*)iroot,'iroot'

        Do i = 1, nact
           ii = i

           If (realcvec) then
              Call dim1_density_R &
              (ii, ii, dr)
              e0 = e0 + dr*eps(i+ninact)
           Else
              Call dim1_density &
              (ii, ii, dr, di)
              if(ABS(di) > 1.0d-10) write(*,*)'1dim density is complex! strange',i,di
              e0 = e0 + dr*eps(i+ninact)
           End if

        End do

        write(*,*)'e0 = Siguma_w(w:active) eps(w)Dww is ',e0

 1000 continue 
      write(*,*)'end'
   end subroutine calce0


