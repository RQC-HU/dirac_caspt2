program loop_test
    implicit none
    integer :: i0, i1, ia, ja, ib, jb, ii, nsec, ninact, i2
    real :: t1, t2
    call cpu_time(t1)
    nsec = 1000
    ninact = 1100
    i0 = 0
    i1 = 0
    Do ia = 1, nsec
        ja = ia
        Do ib = 1, ia - 1
            Do ii = 1, ninact
                jb = ib
                i0 = i0 + 1
                i1 = ninact*(((ia - 1)*(ia - 2))/2 + (ib - 1)) + ii
                if (i0 /= i1) print *, ia, ib, ii, i0, i1
                ! write (*, *) ia, ib, i0!, ninact*(((ia-1)*(ia-2))/2+(ib-1)) + ii !+ (ib-1) * ninact
            end do
        End do
    End do
    i0 = 0
    i1 = 1
    Do ia = 1, nsec
        Do ib = 1, ia - 1
            i0 = i0 + 1
            i1 = ((ia - 2)*(ia - 1))/2 + ib
            if (i0 /= i1) print *, i0, i1
        End do
    End do
    i2 = ((nsec - 2)*(nsec - 1))/2 + nsec - 1
    print *, i1, i2, i0
    call cpu_time(t2)
    print *, "cpu time:", t2 - t1, "seconds."
end program loop_test
