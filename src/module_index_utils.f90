module module_index_utils
    use module_error, only: stop_with_errorcode
    implicit none
    private
    public :: convert_active_to_global_idx, convert_global_to_active_idx, convert_global_to_secondary_idx, &
              sign_even_ret1, sign_odd_ret1

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

    function convert_global_to_active_idx(global_idx) result(active_idx)
        ! ====================================================================================================
        ! Converts a global index to a active index
        ! ====================================================================================================
        ! The global index means the index of sequential ordering of all (inactive + active + secondary) orbitals
        ! The active index means the index of sequential ordering of active orbitals
        use module_global_variables, only: ninact, nact, rank
        integer, intent(in) :: global_idx
        integer :: active_idx
        active_idx = 0

        if (ninact + 1 <= global_idx .and. global_idx <= ninact + nact) then
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
        use module_global_variables, only: ninact, nact, nsec, rank
        integer, intent(in) :: global_idx
        integer :: secondary_idx
        secondary_idx = 0

        if (ninact + nact + 1 <= global_idx .and. global_idx <= ninact + nact + nsec) then
            secondary_idx = global_idx - ninact - nact ! secondary_idx are numbered from 1 to nsec
        else
            if (rank == 0) print '(a,i0,a)', "Error: global_idx = ", global_idx, " is not in the range of secondary orbitals."
            call stop_with_errorcode(3)
        end if

    end function convert_global_to_secondary_idx

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
