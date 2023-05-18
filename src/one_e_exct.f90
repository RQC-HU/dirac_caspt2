! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine one_e_exct(icas_idx, creation_op, annihilation_op, newcas_idx, phase)

! This subroutine is used to calculate the phase and
! the new CAS index for the one-electron excitation

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use module_global_variables

    Implicit NONE

    ! icas_idx: index of the current CAS
    ! creation_op: index of creation operator (electron will be excited to this orbital)
    ! annihilation_op: index of annihilation operator (electron will be excited from this orbital)
    integer, intent(in)  :: icas_idx, creation_op, annihilation_op
    integer, intent(out) :: newcas_idx, phase
    integer :: ia, ib

    ia = 0
    ib = 0
    phase = 0
    newcas_idx = 0

!     Index for only CASCI Active space

    if ((annihilation_op == creation_op) .and. (btest(icas_idx, annihilation_op - 1) .eqv. .true.)) then

        newcas_idx = icas_idx
        phase = 0

    elseif ((btest(icas_idx, annihilation_op - 1) .eqv. .true.) .and. (btest(icas_idx, creation_op - 1) .eqv. .false.)) then

        newcas_idx = icas_idx - 2**(annihilation_op - 1) + 2**(creation_op - 1)

!        calculation of phase
        if (annihilation_op < nact) then
            ia = ia + 2**nact - 2**annihilation_op
        end if
        if (creation_op < nact) then
            ib = ib + 2**nact - 2**creation_op
        end if
        phase = POPCNT(iand(icas_idx, ia)) + POPCNT(iand(icas_idx - 2**(annihilation_op - 1), ib))
        ! odd => (-), even => (+)

    end if

end subroutine one_e_exct
