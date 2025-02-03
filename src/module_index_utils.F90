module module_index_utils
    use module_error, only: stop_with_errorcode
    implicit none
    private
    public :: convert_active_to_global_idx, convert_secondary_to_global_idx, &
              convert_global_to_active_idx, convert_global_to_secondary_idx, &
              get_mo_range, set_global_index, sign_even_ret1, sign_odd_ret1

contains

    function convert_active_to_global_idx(active_idx) result(global_idx)
        ! ====================================================================================================
        ! Converts a active index to a global index
        ! ====================================================================================================
        ! The global index means the index of sequential ordering of all (inactive + active + secondary) orbitals
        ! The active index means the index of sequential ordering of active orbitals
        use module_global_variables, only: ninact, nact, rank
        integer, intent(in) :: active_idx
        integer :: global_idx
        global_idx = 0

        if (1 <= active_idx .and. active_idx <= nact) then
            global_idx = active_idx + ninact
        else
            if (rank == 0) print '(a,i0,a)', "Error: active_idx = ", active_idx, " is not in the range of active orbitals."
            call stop_with_errorcode(1)
        end if
    end function convert_active_to_global_idx

    function convert_secondary_to_global_idx(secondary_idx) result(global_idx)
        ! ====================================================================================================
        ! Converts a secondary index to a global index
        ! ====================================================================================================
        ! The global index means the index of sequential ordering of all (inactive + active + secondary) orbitals
        ! The acsecondarytive index means the index of sequential ordering of secondary orbitals
        use module_global_variables, only: ninact, nact, nsec, rank
        integer, intent(in) :: secondary_idx
        integer :: global_idx
        global_idx = 0

        if (1 <= secondary_idx .and. secondary_idx <= nsec) then
            global_idx = secondary_idx + ninact + nact
        else
            if (rank == 0) print '(a,i0,a)', "Error: secondary_idx = ", secondary_idx, " is not in the range of secondary orbitals."
            call stop_with_errorcode(1)
        end if
    end function convert_secondary_to_global_idx

    function convert_global_to_active_idx(global_idx) result(active_idx)
        ! ====================================================================================================
        ! Converts a global index to a active index
        ! ====================================================================================================
        ! The global index means the index of sequential ordering of all (inactive + active + secondary) orbitals
        ! The active index means the index of sequential ordering of active orbitals
        use module_global_variables, only: rank, ninact, global_act_start, global_act_end
        integer, intent(in) :: global_idx
        integer :: active_idx
        active_idx = 0

        if (global_act_start <= global_idx .and. global_idx <= global_act_end) then
            active_idx = global_idx - ninact ! active_idx are numbered from 1 to nact
        else
            if (rank == 0) print '(a,i0,a)', "Error : global_idx = ", global_idx, " is not in the range of active orbitals."
            call stop_with_errorcode(2)
        end if

    end function convert_global_to_active_idx

    function convert_global_to_secondary_idx(global_idx) result(secondary_idx)
        ! ====================================================================================================
        ! Converts a global index to a secondary index
        ! ====================================================================================================
        ! The global index means the index of sequential ordering of all (secondary + active + secondary) orbitals
        ! The secondary index means the index of sequential ordering of secondary orbitals
        use module_global_variables, only: rank, ninact, nact, global_sec_start, global_sec_end
        integer, intent(in) :: global_idx
        integer :: secondary_idx
        secondary_idx = 0

        if (global_sec_start <= global_idx .and. global_idx <= global_sec_end) then
            secondary_idx = global_idx - ninact - nact ! secondary_idx are numbered from 1 to nsec
        else
            if (rank == 0) print '(a,i0,a)', "Error: global_idx = ", global_idx, " is not in the range of secondary orbitals."
            call stop_with_errorcode(3)
        end if

    end function convert_global_to_secondary_idx

    subroutine get_mo_range(global_mo_idx, mo_start, mo_end)
        use module_global_variables, only: global_inact_start, global_act_start, global_sec_start, &
                                           global_inact_end, global_act_end, global_sec_end
        implicit none
        integer, intent(in) :: global_mo_idx
        integer, intent(out) :: mo_start, mo_end

        if (global_inact_start <= global_mo_idx .and. global_mo_idx <= global_inact_end) then
            mo_start = global_inact_start
            mo_end = global_inact_end
        else if (global_act_start <= global_mo_idx .and. global_mo_idx <= global_act_end) then
            mo_start = global_act_start
            mo_end = global_act_end
        else if (global_sec_start <= global_mo_idx .and. global_mo_idx <= global_sec_end) then
            mo_start = global_sec_start
            mo_end = global_sec_end
        else
            print '(a,i0,a)', "Error: global_mo_idx = ", global_mo_idx, " is not in the range of orbitals."
            call stop_with_errorcode(4)
        end if

    end subroutine get_mo_range

    subroutine set_global_index
        ! Sets the global index of the inactive, active and secondary orbitals
        use module_global_variables, only: ninact, nact, nsec, global_inact_start, global_act_start, global_sec_start, &
                                           global_inact_end, global_act_end, global_sec_end
        implicit none
        if (ninact == 0) then
            global_inact_start = 0
        else
            global_inact_start = 1
        end if
        global_inact_end = ninact
        if (nact == 0) then
            global_act_start = ninact
        else
            global_act_start = ninact + 1
        end if
        global_act_end = ninact + nact
        if (nsec == 0) then
            global_sec_start = ninact + nact
        else
            global_sec_start = ninact + nact + 1
        end if
        global_sec_end = ninact + nact + nsec
    end subroutine set_global_index

    function sign_even_ret1(phase) result(sign)
        ! ====================================================================================================
        ! Returns the sign of a phase, given the phase number
        ! If the phase is even, then the sign is returned as +1
        ! If the phase is odd, then the sign is returned as -1
        ! ====================================================================================================
        use module_global_variables, only: rank
        integer, intent(in) :: phase
        integer :: sign
        ! phase | sign(return value)
        ! ==========================
        ! even  |  +1
        !  odd  |  -1
        sign = 0

        if (mod(phase, 2) == 0) then
            sign = 1
        else
            sign = -1
        end if

        if (sign == 0) then
            if (rank == 0) print *, "Error: sign = 0"
            call stop_with_errorcode(4)
        end if

    end function sign_even_ret1

    function sign_odd_ret1(phase) result(sign)
        ! ====================================================================================================
        ! Returns the sign of a phase, given the phase number
        ! If the phase is even, then the sign is returned as -1
        ! If the phase is odd, then the sign is returned as +1
        ! ====================================================================================================
        use module_global_variables, only: rank
        integer, intent(in) :: phase
        integer :: sign
        ! phase | sign(return value)
        ! ==========================
        ! even  |  -1
        !  odd  |  +1
        sign = 0

        if (mod(phase, 2) == 0) then
            sign = -1
        else
            sign = 1
        end if

        if (sign == 0) then
            if (rank == 0) print *, "Error: sign = 0"
            call stop_with_errorcode(4)
        end if
    end function sign_odd_ret1

end module module_index_utils
