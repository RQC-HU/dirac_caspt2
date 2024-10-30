subroutine dcaspt2_deallocate
    use module_global_variables
    use module_realonly, only: realonly
    implicit none
    ! Deallocate memory
    if (allocated(cir)) then
        Call memminus(KIND(cir), SIZE(cir), 1); deallocate (cir)
    end if
    if (allocated(cii)) then
        Call memminus(KIND(cii), SIZE(cii), 1); deallocate (cii)
    end if
    if (allocated(eigen)) then
        Call memminus(KIND(eigen), SIZE(eigen), 1); deallocate (eigen)
    end if
    if (allocated(eps)) then
        Call memminus(KIND(eps), SIZE(eps), 1); deallocate (eps)
    end if
    if (allocated(irpamo)) then
        Call memminus(KIND(irpamo), SIZE(irpamo), 1); deallocate (irpamo)
    end if
    if (allocated(indmo_cas_to_dirac)) then
        Call memminus(KIND(indmo_cas_to_dirac), SIZE(indmo_cas_to_dirac), 1); deallocate (indmo_cas_to_dirac)
    end if
    if (allocated(indmo_dirac_to_cas)) then
        Call memminus(KIND(indmo_dirac_to_cas), SIZE(indmo_dirac_to_cas), 1); deallocate (indmo_dirac_to_cas)
    end if
    if (allocated(space_idx)) then
        Call memminus(KIND(space_idx), SIZE(space_idx), 1); deallocate (space_idx)
    end if
    if (allocated(MULTB_S)) then
        Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1); deallocate (MULTB_S)
    end if
    if (allocated(MULTB_D)) then
        Call memminus(KIND(MULTB_D), SIZE(MULTB_D), 1); deallocate (MULTB_D)
    end if
    if (allocated(MULTB_DS)) then
        Call memminus(KIND(MULTB_DS), SIZE(MULTB_DS), 1); deallocate (MULTB_DS)
    end if
    if (allocated(irpamo)) then
        Call memminus(KIND(irpamo), SIZE(irpamo), 1); deallocate (irpamo)
    end if
    if (allocated(indmo_cas_to_dirac)) then
        Call memminus(KIND(indmo_cas_to_dirac), SIZE(indmo_cas_to_dirac), 1); deallocate (indmo_cas_to_dirac)
    end if
    if (allocated(indmo_dirac_to_cas)) then
        Call memminus(KIND(indmo_dirac_to_cas), SIZE(indmo_dirac_to_cas), 1); deallocate (indmo_dirac_to_cas)
    end if
    if (allocated(one_elec_int_i)) then
        Call memminus(KIND(one_elec_int_i), SIZE(one_elec_int_i), 1); deallocate (one_elec_int_i)
    end if
    if (allocated(one_elec_int_r)) then
        Call memminus(KIND(one_elec_int_r), SIZE(one_elec_int_r), 1); deallocate (one_elec_int_r)
    end if
    if (allocated(int2r_f1)) then
        Call memminus(KIND(int2r_f1), SIZE(int2r_f1), 1); deallocate (int2r_f1)
    end if
    if (allocated(int2r_f2)) then
        Call memminus(KIND(int2r_f2), SIZE(int2r_f2), 1); deallocate (int2r_f2)
    end if
    if (.not. realonly%is_realonly()) then
        if (allocated(int2i_f1)) then
            Call memminus(KIND(int2i_f1), SIZE(int2i_f1), 1); deallocate (int2i_f1)
        end if
        if (allocated(int2i_f2)) then
            Call memminus(KIND(int2i_f2), SIZE(int2i_f2), 1); deallocate (int2i_f2)
        end if
    end if
    if (allocated(ras1_list)) then
        Call memminus(KIND(ras1_list), SIZE(ras1_list), 1); deallocate (ras1_list)
    end if
    if (allocated(ras2_list)) then
        Call memminus(KIND(ras2_list), SIZE(ras2_list), 1); deallocate (ras2_list)
    end if
    if (allocated(ras3_list)) then
        Call memminus(KIND(ras3_list), SIZE(ras3_list), 1); deallocate (ras3_list)
    end if
end subroutine dcaspt2_deallocate
