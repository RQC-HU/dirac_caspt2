! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE read1mo_co(filename) ! one-electron MO integrals in moint1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module
    use module_file_manager

    Implicit NONE

    integer :: mrconee, isp, nmom, iostat
    character*50, intent(in) :: filename
    logical :: is_end_of_file
    integer :: j0, i0
    double precision, allocatable :: roner(:, :, :), ronei(:, :, :)

    if (rank == 0) then
        print *, 'Enter read1mo_co'
    end if
    mrconee = 10

    realc = .true.

    Allocate (roner(nmo, nmo, scfru)); Call memplus(KIND(roner), SIZE(roner), 1)
    Allocate (ronei(nmo, nmo, scfru)); Call memplus(KIND(ronei), SIZE(ronei), 1)

    call open_unformatted_file(unit=mrconee, file=trim(filename), status="old", optional_action="read")

    rewind (mrconee)
    read (mrconee, iostat=iostat)
    read (mrconee, iostat=iostat)
    read (mrconee, iostat=iostat)
    read (mrconee, iostat=iostat)
    read (mrconee, iostat=iostat)
    read (mrconee, iostat=iostat) (((roner(i0, j0, isp), ronei(i0, j0, isp), j0=1, nmo), i0=1, nmo), isp=1, scfru)

    call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)

! Reverse the sign of ronei if DIRAC version is larger or equal to 21.
    if (dirac_version >= 21) then
        ronei(:, :, :) = -ronei(:, :, :)
    end if

    close (mrconee)

    nmom = ninact + nact + nsec
    Allocate (oner(nmom, nmom)); Call memplus(KIND(oner), SIZE(oner), 1)
    Allocate (onei(nmom, nmom)); Call memplus(KIND(onei), SIZE(onei), 1)

!Iwamuro modify

    do i0 = 1, nmom
        do j0 = 1, nmom
!           oner(i0,j0) = roner(i0,j0,1)
!           onei(i0,j0) = ronei(i0,j0,1)
            oner(indmor(i0), indmor(j0)) = roner(i0, j0, 1) ! using alpha component for a while
            onei(indmor(i0), indmor(j0)) = ronei(i0, j0, 1)
        end do
    end do

    deallocate (roner); Call memminus(KIND(roner), SIZE(roner), 1)
    deallocate (ronei); Call memminus(KIND(ronei), SIZE(ronei), 1)

    if (rank == 0) then
        print *, realc, 'realc'
    end if
end subroutine read1mo_co
