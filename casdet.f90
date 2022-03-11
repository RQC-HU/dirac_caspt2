! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE casdet(totsym)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE

    integer, intent(in)   :: totsym

    integer :: nbitsa
    integer :: i, isym
    integer, allocatable  :: idet0(:)

    Allocate (idet0(ndet))
    idet0 = 0
    ndet = 0

    Do i = 1, 2**nact - 1
        if (nbitsa(i) == nelec) then
            Call detsym(i, isym)
!              if((nsymrpa == 1.and.((isym == totsym).or.(isym == totsym-1))).or. &
            if ((nsymrpa == 1) .or. &
                (nsymrpa /= 1 .and. (isym == totsym))) then
                ndet = ndet + 1
                idet0(ndet) = i
            End if
        End if
    End do

    Allocate (idet(ndet))
    idet(1:ndet) = idet0(1:ndet)
    write (*, *) 'totsym = ', totsym
    write (*, *) 'ndet   = ', ndet
!        write(*,*)idet(1:ndet)
    Deallocate (idet0)

1000 end subroutine casdet

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE detsym(ii, isym)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE

    integer, intent(in)  :: ii
    integer, intent(out) :: isym

    integer :: i, j, jsym

    isym = nsymrpa + 1

    Do i = 1, nact
        if (btest(ii, i - 1) .eqv. .true.) then
            j = i + ninact
            jsym = irpamo(j)
            isym = MULTB(jsym, isym)
        End if
    End do

1000 end subroutine detsym
