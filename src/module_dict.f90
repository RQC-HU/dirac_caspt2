module module_dict
! SPDX-License-Identifier: MIT
! Copyright (c) 2019 Takuma Yoshida
! Original code is written by Takuma Yoshida.
! URL: https://github.com/ysdtkm/fortran_associative_array
! Modified by: Kohei Noda

! This is a dictionary data structure implementation based on treap.
    use module_error, only: stop_with_errorcode
    implicit none
    private
    public :: dict_int, dict_chr_logical, get_val, add, exists, remove, get_keys_vals, get_size, get_kth_key

    type dict_int
        type(node_int), pointer :: root => null()
        integer :: randstate = 1231767121
    contains
        final :: destruct_dict_int
    end type dict_int

    type dict_chr_logical
        type(node_chr_logical), pointer :: root => null()
        integer :: randstate = 1231767121
    contains
        final :: destruct_dict_chr_logical
    end type dict_chr_logical

    type node_int
        type(node_int), pointer :: left => null(), right => null()
        integer, allocatable :: key
        integer :: val
        integer :: pri  ! min-heap
        integer :: cnt = 1
    end type node_int

    type node_chr_logical
        type(node_chr_logical), pointer :: left => null(), right => null()
        character(:), allocatable :: key
        logical :: val
        integer :: pri  ! min-heap
        integer :: cnt = 1
    end type node_chr_logical
    interface get_val
        module procedure get_val_int
        module procedure get_val_chr_logical
    end interface get_val
    interface add
        module procedure add_int
        module procedure add_chr_logical
    end interface add
    interface exists
        module procedure exists_int
        module procedure exists_chr_logical
    end interface exists
    interface remove
        module procedure remove_int
        module procedure remove_chr_logical
    end interface remove
    interface get_keys_vals
        module procedure get_keys_vals_int
        module procedure get_keys_vals_chr_logical
    end interface get_keys_vals
    interface get_size
        module procedure get_size_int
        module procedure get_size_chr_logical
    end interface get_size
    interface get_kth_key
        module procedure get_kth_key_int
        module procedure get_kth_key_chr_logical
    end interface get_kth_key
contains
    ! High level wrapper of dictionary data structure
    pure function xorshift(i)
        implicit none
        integer, intent(in) :: i
        integer :: xorshift
        if (i == 0) then
            xorshift = 1231767121
        else
            xorshift = i
        end if
        xorshift = ieor(xorshift, ishft(xorshift, 13))
        xorshift = ieor(xorshift, ishft(xorshift, -17))
        xorshift = ieor(xorshift, ishft(xorshift, 15))
    end function xorshift

    function get_val_int(t, key)
        implicit none
        type(dict_int), intent(in) :: t
        integer, intent(in) :: key
        type(node_int), pointer :: nd
        integer :: get_val_int
        nd => find_node_int(t%root, key)
        if (.not. associated(nd)) then
            call stop_with_errorcode(1)
        end if
        get_val_int = nd%val
    end function get_val_int

    function get_val_chr_logical(t, key)
        implicit none
        type(dict_chr_logical), intent(in) :: t
        character(*), intent(in) :: key
        type(node_chr_logical), pointer :: nd
        logical :: get_val_chr_logical
        nd => find_node_chr_logical(t%root, key)
        if (.not. associated(nd)) then
            call stop_with_errorcode(1)
        end if
        get_val_chr_logical = nd%val
    end function get_val_chr_logical

    function exists_int(t, key)
        implicit none
        type(dict_int), intent(in) :: t
        integer, intent(in) :: key
        type(node_int), pointer :: nd
        logical :: exists_int
        nd => find_node_int(t%root, key)
        exists_int = (associated(nd))
    end function exists_int

    function exists_chr_logical(t, key)
        implicit none
        type(dict_chr_logical), intent(in) :: t
        character(*), intent(in) :: key
        type(node_chr_logical), pointer :: nd
        logical :: exists_chr_logical
        nd => find_node_chr_logical(t%root, key)
        exists_chr_logical = (associated(nd))
    end function exists_chr_logical

    subroutine add_int(t, key, val)
        ! Add a new key-value pair to the dictionary.
        implicit none
        type(dict_int), intent(inout) :: t
        integer, intent(in) :: key
        integer, intent(in) :: val
        type(node_int), pointer :: nd
        nd => find_node_int(t%root, key)
        if (associated(nd)) then
            nd%val = val
        else  ! This implementation is not optimal
            t%root => insert_int(t%root, key, val, t%randstate)
            t%randstate = xorshift(t%randstate)
        end if
    end subroutine add_int

    subroutine add_chr_logical(t, key, val)
        ! Add a new key-value pair to the dictionary.
        implicit none
        type(dict_chr_logical), intent(inout) :: t
        character(*), intent(in) :: key
        logical, intent(in) :: val
        type(node_chr_logical), pointer :: nd
        nd => find_node_chr_logical(t%root, key)
        if (associated(nd)) then
            nd%val = val
        else  ! This implementation is not optimal
            t%root => insert_chr_logical(t%root, key, val, t%randstate)
            t%randstate = xorshift(t%randstate)
        end if
    end subroutine add_chr_logical

    subroutine remove_int(t, key)
        implicit none
        type(dict_int), intent(inout) :: t
        integer, intent(in) :: key
        t%root => erase_int(t%root, key)
    end subroutine remove_int

    subroutine remove_chr_logical(t, key)
        implicit none
        type(dict_chr_logical), intent(inout) :: t
        character(*), intent(in) :: key
        t%root => erase_chr_logical(t%root, key)
    end subroutine remove_chr_logical

    function get_kth_key_int(t, k)
        implicit none
        type(dict_int), intent(in) :: t
        integer, intent(in) :: k
        type(node_int), pointer :: res
        integer, allocatable :: get_kth_key_int
        if (k < 1 .or. k > my_count_int(t%root)) then
            print *, "get_kth_key failed"
            call stop_with_errorcode(1)
        else
            res => kth_node_int(t%root, k)
            get_kth_key_int = res%key
        end if
    end function get_kth_key_int

    function get_kth_key_chr_logical(t, k)
        implicit none
        type(dict_chr_logical), intent(in) :: t
        integer, intent(in) :: k
        type(node_chr_logical), pointer :: res
        character(:), allocatable :: get_kth_key_chr_logical
        if (k < 1 .or. k > my_count_chr_logical(t%root)) then
            print *, "get_kth_key failed"
            call stop_with_errorcode(1)
        else
            res => kth_node_chr_logical(t%root, k)
            get_kth_key_chr_logical = res%key
        end if
    end function get_kth_key_chr_logical

    subroutine get_keys_vals_int(t, keys, vals, n)
        implicit none
        type(dict_int), intent(in) :: t
        integer, intent(in) :: n
        integer, intent(out) :: keys(n)
        integer, intent(out) :: vals(n)
        integer :: counter
        if (my_count_int(t%root) /= n) call stop_with_errorcode(1)
        counter = 0
        call inorder_int(t%root, keys, vals, counter)
    end subroutine get_keys_vals_int

    subroutine get_keys_vals_chr_logical(t, keys, vals, n)
        implicit none
        type(dict_chr_logical), intent(in) :: t
        integer, intent(in) :: n
        character(*), intent(out) :: keys(n)
        logical, intent(out) :: vals(n)
        integer :: counter
        if (my_count_chr_logical(t%root) /= n) call stop_with_errorcode(1)
        counter = 0
        call inorder_chr_logical(t%root, keys, vals, counter)
    end subroutine get_keys_vals_chr_logical

    function get_size_int(t)
        implicit none
        type(dict_int), intent(in) :: t
        integer :: get_size_int
        get_size_int = my_count_int(t%root)
    end function get_size_int

    function get_size_chr_logical(t)
        implicit none
        type(dict_chr_logical), intent(in) :: t
        integer :: get_size_chr_logical
        get_size_chr_logical = my_count_chr_logical(t%root)
    end function get_size_chr_logical

    subroutine destruct_dict_int(t)
        implicit none
        type(dict_int), intent(inout) :: t
        call delete_all_int(t%root)
    end subroutine destruct_dict_int

    subroutine destruct_dict_chr_logical(t)
        implicit none
        type(dict_chr_logical), intent(inout) :: t
        call delete_all_chr_logical(t%root)
    end subroutine destruct_dict_chr_logical

    ! Low level data structure and operations of treap.
    ! This allows multiple nodes with a same key.

    subroutine update_int(root)
        implicit none
        type(node_int), pointer, intent(in) :: root
        root%cnt = my_count_int(root%left) + my_count_int(root%right) + 1
    end subroutine update_int

    subroutine update_chr_logical(root)
        implicit none
        type(node_chr_logical), pointer, intent(in) :: root
        root%cnt = my_count_chr_logical(root%left) + my_count_chr_logical(root%right) + 1
    end subroutine update_chr_logical

    function my_count_int(root)
        implicit none
        type(node_int), pointer, intent(in) :: root
        integer :: my_count_int
        if (associated(root)) then
            my_count_int = root%cnt
        else
            my_count_int = 0
        end if
    end function my_count_int

    function my_count_chr_logical(root)
        implicit none
        type(node_chr_logical), pointer, intent(in) :: root
        integer :: my_count_chr_logical
        if (associated(root)) then
            my_count_chr_logical = root%cnt
        else
            my_count_chr_logical = 0
        end if
    end function my_count_chr_logical

    function rotate_ccw_int(root)
        implicit none
        type(node_int), pointer, intent(in) :: root
        type(node_int), pointer :: tmp, rotate_ccw_int
        if (.not. associated(root%right)) call stop_with_errorcode(1)
        tmp => root%right
        root%right => tmp%left
        tmp%left => root
        rotate_ccw_int => tmp
        call update_int(root)
        call update_int(tmp)
    end function rotate_ccw_int

    function rotate_ccw_chr_logical(root)
        implicit none
        type(node_chr_logical), pointer, intent(in) :: root
        type(node_chr_logical), pointer :: tmp, rotate_ccw_chr_logical
        if (.not. associated(root%right)) call stop_with_errorcode(1)
        tmp => root%right
        root%right => tmp%left
        tmp%left => root
        rotate_ccw_chr_logical => tmp
        call update_chr_logical(root)
        call update_chr_logical(tmp)
    end function rotate_ccw_chr_logical

    function rotate_cw_int(root)
        implicit none
        type(node_int), pointer, intent(in) :: root
        type(node_int), pointer :: tmp, rotate_cw_int
        if (.not. associated(root%left)) call stop_with_errorcode(1)
        tmp => root%left
        root%left => tmp%right
        tmp%right => root
        rotate_cw_int => tmp
        call update_int(root)
        call update_int(tmp)
    end function rotate_cw_int

    function rotate_cw_chr_logical(root)
        implicit none
        type(node_chr_logical), pointer, intent(in) :: root
        type(node_chr_logical), pointer :: tmp, rotate_cw_chr_logical
        if (.not. associated(root%left)) call stop_with_errorcode(1)
        tmp => root%left
        root%left => tmp%right
        tmp%right => root
        rotate_cw_chr_logical => tmp
        call update_chr_logical(root)
        call update_chr_logical(tmp)
    end function rotate_cw_chr_logical

    recursive function insert_int(root, key, val, pri) result(res)
        implicit none
        type(node_int), pointer, intent(in) :: root
        integer, intent(in) :: pri
        integer, intent(in) :: key
        integer, intent(in) :: val
        type(node_int), pointer :: res

        if (.not. associated(root)) then
            allocate (res)
            res%key = key
            res%pri = pri
            res%val = val
        else
            res => root
            if (key > root%key) then
                root%right => insert_int(root%right, key, val, pri)
                call update_int(root)
                if (root%pri > root%right%pri) then
                    res => rotate_ccw_int(res)
                end if
            else
                root%left => insert_int(root%left, key, val, pri)
                call update_int(root)
                if (root%pri > root%left%pri) then
                    res => rotate_cw_int(res)
                end if
            end if
        end if
    end function insert_int

    recursive function insert_chr_logical(root, key, val, pri) result(res)
        implicit none
        type(node_chr_logical), pointer, intent(in) :: root
        integer, intent(in) :: pri
        character(*), intent(in) :: key
        logical, intent(in) :: val
        type(node_chr_logical), pointer :: res

        if (.not. associated(root)) then
            allocate (res)
            res%key = key
            res%pri = pri
            res%val = val
        else
            res => root
            if (key > root%key) then
                root%right => insert_chr_logical(root%right, key, val, pri)
                call update_chr_logical(root)
                if (root%pri > root%right%pri) then
                    res => rotate_ccw_chr_logical(res)
                end if
            else
                root%left => insert_chr_logical(root%left, key, val, pri)
                call update_chr_logical(root)
                if (root%pri > root%left%pri) then
                    res => rotate_cw_chr_logical(res)
                end if
            end if
        end if
    end function insert_chr_logical

    recursive function erase_int(root, key) result(res)
        implicit none
        type(node_int), pointer, intent(in) :: root
        integer, intent(in) :: key
        type(node_int), pointer :: res, tmp

        if (.not. associated(root)) then
            print *, "Erase failed"
            call stop_with_errorcode(1)
        end if

        if (key < root%key) then
            root%left => erase_int(root%left, key)
            res => root
        else if (key > root%key) then
            root%right => erase_int(root%right, key)
            res => root
        else
            if ((.not. associated(root%left)) .or. (.not. associated(root%right))) then
                tmp => root
                if (.not. associated(root%left)) then
                    res => root%right
                else
                    res => root%left
                end if
                deallocate (tmp)
            else
                if (root%left%pri < root%right%pri) then
                    res => rotate_ccw_int(root)
                    res%left => erase_int(res%left, key)
                else
                    res => rotate_cw_int(root)
                    res%right => erase_int(res%right, key)
                end if
            end if
        end if
        if (associated(res)) call update_int(res)
    end function erase_int

    recursive function erase_chr_logical(root, key) result(res)
        implicit none
        type(node_chr_logical), pointer, intent(in) :: root
        character(*), intent(in) :: key
        type(node_chr_logical), pointer :: res, tmp

        if (.not. associated(root)) then
            print *, "Erase failed"
            call stop_with_errorcode(1)
        end if

        if (key < root%key) then
            root%left => erase_chr_logical(root%left, key)
            res => root
        else if (key > root%key) then
            root%right => erase_chr_logical(root%right, key)
            res => root
        else
            if ((.not. associated(root%left)) .or. (.not. associated(root%right))) then
                tmp => root
                if (.not. associated(root%left)) then
                    res => root%right
                else
                    res => root%left
                end if
                deallocate (tmp)
            else
                if (root%left%pri < root%right%pri) then
                    res => rotate_ccw_chr_logical(root)
                    res%left => erase_chr_logical(res%left, key)
                else
                    res => rotate_cw_chr_logical(root)
                    res%right => erase_chr_logical(res%right, key)
                end if
            end if
        end if
        if (associated(res)) call update_chr_logical(res)
    end function erase_chr_logical

    recursive function find_node_int(root, key) result(res)
        implicit none
        type(node_int), pointer, intent(in) :: root
        integer, intent(in) :: key
        type(node_int), pointer :: res
        if (.not. associated(root)) then
            res => null()
        else if (root%key == key) then
            res => root
        else if (key < root%key) then
            res => find_node_int(root%left, key)
        else
            res => find_node_int(root%right, key)
        end if
    end function find_node_int

    recursive function find_node_chr_logical(root, key) result(res)
        implicit none
        type(node_chr_logical), pointer, intent(in) :: root
        character(*), intent(in) :: key
        type(node_chr_logical), pointer :: res
        if (.not. associated(root)) then
            res => null()
        else if (root%key == key) then
            res => root
        else if (key < root%key) then
            res => find_node_chr_logical(root%left, key)
        else
            res => find_node_chr_logical(root%right, key)
        end if
    end function find_node_chr_logical

    recursive function kth_node_int(root, k) result(res)
        implicit none
        type(node_int), pointer, intent(in) :: root
        integer, intent(in) :: k
        type(node_int), pointer :: res
        if (.not. associated(root)) then
            res => null()
        else if (k <= my_count_int(root%left)) then
            res => kth_node_int(root%left, k)
        else if (k == my_count_int(root%left) + 1) then
            res => root
        else
            res => kth_node_int(root%right, k - my_count_int(root%left) - 1)
        end if
    end function kth_node_int

    recursive function kth_node_chr_logical(root, k) result(res)
        implicit none
        type(node_chr_logical), pointer, intent(in) :: root
        integer, intent(in) :: k
        type(node_chr_logical), pointer :: res
        if (.not. associated(root)) then
            res => null()
        else if (k <= my_count_chr_logical(root%left)) then
            res => kth_node_chr_logical(root%left, k)
        else if (k == my_count_chr_logical(root%left) + 1) then
            res => root
        else
            res => kth_node_chr_logical(root%right, k - my_count_chr_logical(root%left) - 1)
        end if
    end function kth_node_chr_logical

    recursive subroutine delete_all_int(root)
        implicit none
        type(node_int), pointer, intent(inout) :: root
        if (.not. associated(root)) return

        call delete_all_int(root%left)
        call delete_all_int(root%right)
        deallocate (root)
        nullify (root)
    end subroutine delete_all_int

    recursive subroutine delete_all_chr_logical(root)
        implicit none
        type(node_chr_logical), pointer, intent(inout) :: root
        if (.not. associated(root)) return

        call delete_all_chr_logical(root%left)
        call delete_all_chr_logical(root%right)
        deallocate (root)
        nullify (root)
    end subroutine delete_all_chr_logical

    recursive subroutine inorder_int(root, keys, vals, counter)
        implicit none
        type(node_int), pointer, intent(in) :: root
        integer, intent(inout) :: keys(:)
        integer, intent(inout) :: vals(:)
        integer, intent(inout) :: counter
        if (.not. associated(root)) return

        call inorder_int(root%left, keys, vals, counter)
        counter = counter + 1
        keys(counter) = root%key
        vals(counter) = root%val
        call inorder_int(root%right, keys, vals, counter)
    end subroutine inorder_int

    recursive subroutine inorder_chr_logical(root, keys, vals, counter)
        implicit none
        type(node_chr_logical), pointer, intent(in) :: root
        character(*), intent(inout) :: keys(:)
        logical, intent(inout) :: vals(:)
        integer, intent(inout) :: counter
        if (.not. associated(root)) return

        call inorder_chr_logical(root%left, keys, vals, counter)
        counter = counter + 1
        keys(counter) = root%key
        vals(counter) = root%val
        call inorder_chr_logical(root%right, keys, vals, counter)
    end subroutine inorder_chr_logical
end module module_dict
