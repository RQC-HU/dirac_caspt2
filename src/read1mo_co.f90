! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE read1mo_co(filename) ! one-electron MO integrals in moint1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE

    integer :: mrconee, isp, nmom
    character*50, intent(in) :: filename
    integer :: j0, i0
    double precision, allocatable :: roner(:, :, :), ronei(:, :, :)

    if (rank == 0) then
        write (*, *) 'Enter read1mo_co'
    end if
    mrconee = 10

    realc = .true.

    Allocate (roner(nmo, nmo, scfru)); Call memplus(KIND(roner), SIZE(roner), 1)
    Allocate (ronei(nmo, nmo, scfru)); Call memplus(KIND(ronei), SIZE(ronei), 1)

    open (mrconee, file=trim(filename), status='old', form='unformatted', err=10)
    rewind (mrconee)
    read (mrconee, err=10)
    read (mrconee, err=10)
    read (mrconee, err=10)
    read (mrconee, err=10)
    read (mrconee, err=10)
    read (mrconee, err=10) (((roner(i0, j0, isp), ronei(i0, j0, isp), j0=1, nmo), i0=1, nmo), isp=1, scfru)

! Reverse the sign of ronei if DIRAC version is larger or equal to 21.
if (dirac_version >= 21) then
    ronei(:, :, :) = -ronei(:, :, :)
end if
    ! if (rank == 0) then
    !     write (*, *) 'Noda RONER RONEI CHECK'
    !     write (*, *) "i0,j0,RONER(i0,j0,1),RONEI(i0,j0,1)"
    !     do i0 = 1, nmo
    !         do j0 = 1, nmo
    !             Write (*, '(2I4,2X,2F8.5)') i0, j0, RONER(i0, j0, 1), RONEI(i0, j0, 1)
    !         end do
    !     end do
    ! end if
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
        write (*, *) realc, 'realc'
    end if
    goto 1000

10  if (rank == 0) then
        write (*, *) 'err 10 mo1'
    end if
    go to 1000
    if (rank == 0) then
        write (*, *) 'err 11 mo1'
    end if
    go to 1000

1000 end subroutine read1mo_co
