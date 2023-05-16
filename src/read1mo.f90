! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE read1mo(filename) ! one-electron MO integrals in MRCONEE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use module_global_variables
    use module_file_manager

    Implicit NONE

    integer :: unit_mrconee, isp, nmom, iostat
    character*50, intent(in) :: filename
    logical :: is_end_of_file
    integer :: j0, i0
    double precision, allocatable :: roner(:, :, :), ronei(:, :, :)

    if (rank == 0) then
        print *, 'Enter read1mo'
    end if

    scfru = 1
    Allocate (roner(nmo, nmo, scfru)); Call memplus(KIND(roner), SIZE(roner), 1)
    Allocate (ronei(nmo, nmo, scfru)); Call memplus(KIND(ronei), SIZE(ronei), 1)

    call open_unformatted_file(unit=unit_mrconee, file=trim(filename), status="old", optional_action="read")

    rewind (unit_mrconee)
    read (unit_mrconee, iostat=iostat)
    read (unit_mrconee, iostat=iostat)
    read (unit_mrconee, iostat=iostat)
    read (unit_mrconee, iostat=iostat)
    read (unit_mrconee, iostat=iostat)
    read (unit_mrconee, iostat=iostat) (((roner(i0, j0, isp), ronei(i0, j0, isp), j0=1, nmo), i0=1, nmo), isp=1, scfru)

    call check_iostat(iostat=iostat, file=trim(filename), end_of_file_reached=is_end_of_file)

    ! Reverse the sign of ronei if DIRAC version is larger or equal to 21.
    if (dirac_version >= 21) then
        ronei(:, :, :) = -ronei(:, :, :)
    end if

    close (unit_mrconee)

    nmom = ninact + nact + nsec
    Allocate (one_elec_int_r(nmom, nmom)); Call memplus(KIND(one_elec_int_r), SIZE(one_elec_int_r), 1)
    Allocate (one_elec_int_i(nmom, nmom)); Call memplus(KIND(one_elec_int_i), SIZE(one_elec_int_i), 1)

    ! Store the one-electron integrals in energy order (CASPT2 order)
    ! one_elec_int_[r,i] are CASPT2 order
    ! roner, ronei are DIRAC order
    do i0 = 1, nmom
        do j0 = 1, nmom
            one_elec_int_r(indmo_dirac_to_cas(i0), indmo_dirac_to_cas(j0)) = roner(i0, j0, 1) ! using alpha component for a while
            one_elec_int_i(indmo_dirac_to_cas(i0), indmo_dirac_to_cas(j0)) = ronei(i0, j0, 1)
        end do
    end do

    Call memminus(KIND(roner), SIZE(roner), 1); deallocate (roner)
    Call memminus(KIND(ronei), SIZE(ronei), 1); deallocate (ronei)
end subroutine read1mo
