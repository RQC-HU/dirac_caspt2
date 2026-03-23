! collection of subroutines defined/implemented in upsteam DIRAC.
! don't compile this for integrated DIRAC caspt2 module.

! https://gitlab.com/dirac/dirac/-/blob/d7dccc656ba6cb0c88a7d928b555aed8c53e054e/src/dirac/dirgp.F#L4293
subroutine atomic_id(ID, KP, J, MJ, LL)
    implicit none
    integer, intent(in)  :: ID
    integer, intent(out) :: KP, J, MJ, LL
    integer             :: INDI, INDJ
    if (ID == 0) then
        KP = 0; J = 0; MJ = 0; LL = 0
        return
    end if
    INDJ = INT(SQRT(DBLE(2*ABS(ID)) + 0.25D0) + 0.4999D0)
    INDI = ABS(ID) - INDJ*(INDJ - 1)/2
    KP = INDJ*ID/ABS(ID)
    MJ = (2*INDI - 1)*(-1)**(INDI + 1)
    J = 2*INDJ - 1
    IF (KP > 0) THEN
        LL = INDJ
    ELSE
        LL = INDJ - 1
    END IF
end subroutine atomic_id
