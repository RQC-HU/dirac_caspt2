program main
    use four_caspt2_module, only: ras3_list
    use read_input_module
    implicit none
    integer :: count
    count = 1
    print *, "start"
!     open (2, file='ref.ras3')
!     do
!         read (2, *, end=20) ras3(count)
!         count = count + 1
!     end do
! 20  close (2)
!     count = count - 1
!     open (1, file='file2', form='formatted')
!     write (1, *) ras3(1:count)
!     close (1)
!     open (1, file='file2', form='formatted')
!     read (1, '(a)') chr
!     print *, chr
!     close (1)
    open (5, file='input', form='formatted')
    call ras_read(ras3_list, 3)
    print *, "end"
    close (5)
    print *, "open"
    open (2, file='result', form="formatted")
    print *, 'before write'
    write (2, *) ras3_list
    print *, 'end write'
    close (2)
end program main