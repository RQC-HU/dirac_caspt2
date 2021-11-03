program parallel
    use omp_lib
    implicit none
    integer :: i, j = 0, tid, omp_max
    omp_max = omp_get_max_threads()
    ! call omp_set_num_threads(10)
    write(20,*)"max threads : ", omp_max

    !$omp parallel do private(tid)
    do i = 1, 8
        tid = omp_get_thread_num()
        write(20,*) i, tid
    end do
    !$omp end parallel do
end program parallel
