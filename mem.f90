!  =================================================
   SUBROUTINE memplus (i, j, k)

!  =================================================
   use four_caspt2_module

     Implicit NONE
     integer :: i, j, k

      tmem = tmem + i*j*k


    end subroutine memplus


!  =================================================
   SUBROUTINE memminus (i, j, k)

!  =================================================
   use four_caspt2_module

     Implicit NONE
     integer :: i, j, k

      tmem = tmem - i*j*k

    end subroutine memminus


