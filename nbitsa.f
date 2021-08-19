!  =================================================
   INTEGER FUNCTION nbitsa (jdet) !result (nmb)

!  =================================================

!   use four_caspt2_module ,only : norb

     Implicit NONE
     integer :: jdet, nmb ,i 


     nmb=0
     do i = 0, 31
       if(btest(jdet,i)) nmb=nmb+1
     end do
  
     nbitsa=nmb


    end FUNCTION nbitsa


