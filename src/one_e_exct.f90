! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine one_e_exct(iidet, creat, anhi, newidet, phase)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE

    integer, intent(in)  :: iidet, creat, anhi
    integer, intent(out) :: newidet, phase

    integer :: ia, ib, indu, indv

    ia = 0
    ib = 0
    phase = 0
    newidet = 0

!     Index for only CASCI Active space

    indu = creat
    indv = anhi

    if ((indv == indu) .and. (btest(iidet, indv - 1) .eqv. .true.)) then

        newidet = iidet
        phase = 0

    elseif ((btest(iidet, indv - 1) .eqv. .true.) .and. (btest(iidet, indu - 1) .eqv. .false.)) then

        newidet = iidet - 2**(indv - 1) + 2**(indu - 1)
!        calculation of phase

        if (indv < nact) then
            ia = ia + 2**nact - 2**indv
        end if
        if (indu < nact) then
            ib = ib + 2**nact - 2**indu
        end if
        phase = POPCNT(iand(iidet, ia)) + POPCNT(iand(iidet - 2**(indv - 1), ib))
        ! odd => (-), even => (+)

    end if

end subroutine one_e_exct
