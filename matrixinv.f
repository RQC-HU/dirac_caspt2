!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine matinv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none
 
  real*8 :: a(x,x),c,dum
  integer :: x, y, z, w

!------------------------------------------

   do z=1,x
      c=a(w,w)
      a(w,w)=1
 
      do y=1,x
        a(w,z)=a(w,z)/c
      enddo
     
        do y=1,x
        if( y/=w ) then
          dum=a(y,w)
          a(y,w)=0
          do z=1,x
            a(y,z)=a(y,z)-dum*a(w,z)
          enddo
        endif
        enddo
   enddo

   end subroutine matinv
