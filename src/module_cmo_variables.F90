module module_cmo_variables
    implicit none

    real(8), allocatable      :: BUF(:)  ! One dimensional array representing MO coeff. read from DFPCMO
    real(8), allocatable      :: eval(:)
    integer, allocatable :: syminfo(:), kappa(:)
    integer :: cmo_nfsym, cmo_nz
    integer :: positronic_mo(2), electronic_mo(2), basis_ao(2), basis_all(2), mo(2)
contains
    subroutine reset_cmo_variables
        implicit none

        if (allocated(BUF)) deallocate (BUF)
        if (allocated(eval)) deallocate (eval)
        if (allocated(syminfo)) deallocate (syminfo)
        if (allocated(kappa)) deallocate (kappa)
        cmo_nfsym = 0; cmo_nz = 0; positronic_mo = 0; electronic_mo = 0; basis_ao = 0; basis_all = 0; mo = 0
    end subroutine reset_cmo_variables
end module module_cmo_variables
