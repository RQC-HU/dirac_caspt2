! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE calce0(e0)

! Calculates the eigenvalue e0 = <0|H_0|0> = Sum_w eps_w D_ww
! where D_ww is the 1dim density matrix
! (e0 is the equation 22, 2nd terms of the folloing paper.
!  J. Phys. Chem. 1990, 94, 5483-5488, Kerstin. Andersson,
!  Per Aake. Malmqvist, Bjoern O. Roos, Andrzej J. Sadlej, and Krzysztof. Wolinski
!  https://pubs.acs.org/doi/pdf/10.1021/j100377a012)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables
    use module_realonly, only: realonly
    use module_index_utils, only: convert_active_to_global_idx

    Implicit NONE

    real(8), intent(out):: e0

    integer :: i
    real(8)  :: dr, di

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    e0 = 0.0d+00
    dr = 0.0d+00
    di = 0.0d+00
    if (rank == 0) print *, 'selectroot = ', selectroot
    Do i = 1, nact

        If (realonly%is_realonly()) then
            Call dim1_density_R(i, i, dr)
        Else
            Call dim1_density(i, i, dr, di)
            if (ABS(di) > 1.0d-10 .and. rank == 0) print *, '1dim density is complex! strange', i, di
        End if
        e0 = e0 + dr*eps(convert_active_to_global_idx(i))

    End do

    if (rank == 0) then
        print *, 'e0 = Siguma_w(w:active) eps(w)Dww is ', e0
    end if
end subroutine calce0
