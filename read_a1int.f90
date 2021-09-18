program reada1int
    implicit none
    INTEGER :: i, j, k, l, cnt = 0,iostat
    real*8 :: cint1, cint2
    complex*16 :: c
    open (1, file='A1int', form='unformatted', status='unknown')
    open (2, file='A1int_formatted', form='formatted', status='unknown')
    ! a1loop1: do
60  read (1, err=10, end=20, IOSTAT=iostat) i, j, k, l, c
    write (2, '(4I4,2E20.10)') i, j, k, l, c
    goto 100
10  write (*, *) 'err read A1int', cnt, iostat
100    cnt = cnt + 1
    goto 60
    ! end do a1loop1
20  write (*, *) 'end read A1int'
    close (1)
    close (2)

    open (3, file='A1int4', form='unformatted', status='unknown')
    open (4, file='A1int4_formatted', form='formatted', status='unknown')
    ! a1loop4: do
61  read (3, err=11, end=21) i, j, k, l, c
    write (4, '(4I4,2E20.10)') i, j, k, l, c
    goto 61
    ! end do a1loop4
11  write (*, *) 'err read A1int4'
21  write (*, *) 'end read A1int4'
    close (3)
    close (4)
end program reada1int
