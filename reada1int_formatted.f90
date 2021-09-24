program reada1int_formatted
    implicit none
    integer     :: i, j, k, l, cnt = 0
    ! real(8) :: c,d
    CHARACTER(100) :: str
    complex(16) :: c
    open (1, file='A1int', form='formatted', status='old')
! 60  write(*,*) 'count', cnt
60    read (1, '(4I4, 2e20.10)', err=10, end=20) i, j, k, l, c
    ! read (q, *, err=10, end=20) str
    write(*,'(4I4, 2e20.10)') i,j,k,l,c
    go to 11
10  write (*, *) ' err read', cnt
11  cnt = cnt + 1
    go to 60
20 write(*,*) 'end read'
end program reada1int_formatted
