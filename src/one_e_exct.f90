! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine one_e_exct(icas_idx, creat, anhi, newcas_idx, phase)

! This subroutine is used to calculate the phase and
! the new CAS index for the one-electron excitation

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use module_global_variables

    Implicit NONE

    integer, intent(in)  :: icas_idx, creat, anhi
    integer, intent(out) :: newcas_idx, phase

    integer :: ia, ib, indu, indv

    ia = 0
    ib = 0
    phase = 0
    newcas_idx = 0

!     Index for only CASCI Active space

    indu = creat
    indv = anhi

    if ((indv == indu) .and. (btest(icas_idx, indv - 1) .eqv. .true.)) then

        newcas_idx = icas_idx
        phase = 0

    elseif ((btest(icas_idx, indv - 1) .eqv. .true.) .and. (btest(icas_idx, indu - 1) .eqv. .false.)) then

        newcas_idx = icas_idx - 2**(indv - 1) + 2**(indu - 1)
!        calculation of phase

        if (indv < nact) then
            ia = ia + 2**nact - 2**indv
        end if
        if (indu < nact) then
            ib = ib + 2**nact - 2**indu
        end if
        phase = POPCNT(iand(icas_idx, ia)) + POPCNT(iand(icas_idx - 2**(indv - 1), ib))
        ! odd => (-), even => (+)

    end if

end subroutine one_e_exct
