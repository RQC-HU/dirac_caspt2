!  =================================================
SUBROUTINE memplus(i, j, k)

!  =================================================
    use module_global_variables

    Implicit NONE
    integer :: i, j, k

    tmem = tmem + i*j*k

end subroutine memplus

!  =================================================
SUBROUTINE memminus(i, j, k)

!  =================================================
    use module_global_variables

    Implicit NONE
    integer :: i, j, k

    tmem = tmem - i*j*k

end subroutine memminus
