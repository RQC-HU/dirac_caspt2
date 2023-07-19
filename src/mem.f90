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

subroutine write_allocated_memory_size
    use module_global_variables, only: rank, tmem
    implicit none
    real(8), parameter :: KB = 1024.0d+00, MB = KB*1024.0d+00, GB = MB*1024.0d+00

    if (rank == 0) then
        if (tmem < KB) then
            print '(A,F10.2,A)', 'Master rank allocated memory size: ', tmem, ' bytes'
        else if (tmem < MB) then
            print '(A,F10.2,A)', 'Master rank allocated memory size: ', tmem/KB, ' KB'
        else if (tmem < GB) then
            print '(A,F10.2,A)', 'Master rank allocated memory size: ', tmem/MB, ' MB'
        else
            print '(A,F10.2,A)', 'Master rank allocated memory size: ', tmem/GB, ' GB'
        end if
    end if
end subroutine write_allocated_memory_size
