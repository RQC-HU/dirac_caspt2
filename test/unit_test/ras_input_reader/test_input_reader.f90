program main
    use four_caspt2_module, only: ras3_list, max_ras3_spinor_num
    use read_input_module
    implicit none
    integer :: i, count
    character(100) :: chr
    integer :: ras3(max_ras3_spinor_num)
    ras3 = 0
    count = 1
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
    open (5, file='file2', form='formatted')
    call ras3_read
    print *, "end"
    close (5)
    print *, "open"
    open (2, file='file', form="formatted")
    print *, 'before write'
    write (2, *) ras3_list
    print *, 'end write'
    close (2)
end program main
