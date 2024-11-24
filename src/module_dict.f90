module module_dict
! SPDX-License-Identifier: MIT
! Copyright (c) 2019 Takuma Yoshida
! Original code is written by Takuma Yoshida.
! URL: https://github.com/ysdtkm/fortran_associative_array
! Modified by: Kohei Noda

! This is a dictionary data structure implementation based on treap.
    use, intrinsic :: iso_fortran_env, only: int64
    use module_error, only: stop_with_errorcode
    implicit none
    private
    public :: dict, get_val, add, exists, remove, get_keys_vals, get_size, get_kth_key

    type dict
        type(node), pointer :: root => null()
        integer(kind=int64) :: randstate = 1231767121
    contains
        final :: destruct_dict
    end type dict

    type node
        type(node), pointer :: left => null(), right => null()
        integer(kind=int64), allocatable :: key
        integer(kind=int64) :: val
        integer(kind=int64) :: pri  ! min-heap
        integer(kind=int64) :: cnt = 1
    end type node
contains
    ! High level wrapper of dictionary data structure
    pure function xorshift(i)
        implicit none
        integer(kind=int64), intent(in) :: i
        integer(kind=int64) :: xorshift
        if (i == 0) then
            xorshift = 1231767121
        else
            xorshift = i
        end if
        xorshift = ieor(xorshift, ishft(xorshift, 13))
        xorshift = ieor(xorshift, ishft(xorshift, -17))
        xorshift = ieor(xorshift, ishft(xorshift, 15))
    end function xorshift

    function get_val(t, key)
        implicit none
        type(dict), intent(in) :: t
        integer(kind=int64), intent(in) :: key
        type(node), pointer :: nd
        integer(kind=int64) :: get_val
        nd => find_node(t%root, key)
        if (.not. associated(nd)) then
            call stop_with_errorcode(1)
        end if
        get_val = nd%val
    end function get_val

    function exists(t, key)
        implicit none
        type(dict), intent(in) :: t
        integer(kind=int64), intent(in) :: key
        type(node), pointer :: nd
        logical :: exists
        nd => find_node(t%root, key)
        exists = (associated(nd))
    end function exists

    subroutine add(t, key, val)
        ! Add a new key-value pair to the dictionary.
        implicit none
        type(dict), intent(inout) :: t
        integer(kind=int64), intent(in) :: key
        integer(kind=int64), intent(in) :: val
        type(node), pointer :: nd
        nd => find_node(t%root, key)
        if (associated(nd)) then
            nd%val = val
        else  ! This implementation is not optimal
            t%root => insert(t%root, key, val, t%randstate)
            t%randstate = xorshift(t%randstate)
        end if
    end subroutine add

    subroutine remove(t, key)
        implicit none
        type(dict), intent(inout) :: t
        integer(kind=int64), intent(in) :: key
        t%root => erase(t%root, key)
    end subroutine remove

    function get_kth_key(t, k)
        implicit none
        type(dict), intent(in) :: t
        integer(kind=int64), intent(in) :: k
        type(node), pointer :: res
        integer(kind=int64), allocatable :: get_kth_key
        if (k < 1 .or. k > my_count(t%root)) then
            print *, "get_kth_key failed"
            call stop_with_errorcode(1)
        else
            res => kth_node(t%root, k)
            get_kth_key = res%key
        end if
    end function get_kth_key

    subroutine get_keys_vals(t, keys, vals, n)
        implicit none
        type(dict), intent(in) :: t
        integer(kind=int64), intent(in) :: n
        integer(kind=int64), intent(out) :: keys(n)
        integer(kind=int64), intent(out) :: vals(n)
        integer(kind=int64) :: counter
        if (my_count(t%root) /= n) call stop_with_errorcode(1)
        counter = 0
        call inorder(t%root, keys, vals, counter)
    end subroutine get_keys_vals

    function get_size(t)
        implicit none
        type(dict), intent(in) :: t
        integer(kind=int64) :: get_size
        get_size = my_count(t%root)
    end function get_size

    subroutine destruct_dict(t)
        implicit none
        type(dict), intent(inout) :: t
        call delete_all(t%root)
    end subroutine destruct_dict
    ! Low level data structure and operations of treap.
    ! This allows multiple nodes with a same key.

    subroutine update(root)
        implicit none
        type(node), pointer, intent(in) :: root
        root%cnt = my_count(root%left) + my_count(root%right) + 1
    end subroutine update

    function my_count(root)
        implicit none
        type(node), pointer, intent(in) :: root
        integer(kind=int64) :: my_count
        if (associated(root)) then
            my_count = root%cnt
        else
            my_count = 0
        end if
    end function my_count

    function rotate_ccw(root)
        implicit none
        type(node), pointer, intent(in) :: root
        type(node), pointer :: tmp, rotate_ccw
        if (.not. associated(root%right)) call stop_with_errorcode(1)
        tmp => root%right
        root%right => tmp%left
        tmp%left => root
        rotate_ccw => tmp
        call update(root)
        call update(tmp)
    end function rotate_ccw

    function rotate_cw(root)
        implicit none
        type(node), pointer, intent(in) :: root
        type(node), pointer :: tmp, rotate_cw
        if (.not. associated(root%left)) call stop_with_errorcode(1)
        tmp => root%left
        root%left => tmp%right
        tmp%right => root
        rotate_cw => tmp
        call update(root)
        call update(tmp)
    end function rotate_cw

    recursive function insert(root, key, val, pri) result(res)
        implicit none
        type(node), pointer, intent(in) :: root
        integer(kind=int64), intent(in) :: pri
        integer(kind=int64), intent(in) :: key
        integer(kind=int64), intent(in) :: val
        type(node), pointer :: res

        if (.not. associated(root)) then
            allocate (res)
            res%key = key
            res%pri = pri
            res%val = val
        else
            res => root
            if (key > root%key) then
                root%right => insert(root%right, key, val, pri)
                call update(root)
                if (root%pri > root%right%pri) then
                    res => rotate_ccw(res)
                end if
            else
                root%left => insert(root%left, key, val, pri)
                call update(root)
                if (root%pri > root%left%pri) then
                    res => rotate_cw(res)
                end if
            end if
        end if
    end function insert

    recursive function erase(root, key) result(res)
        implicit none
        type(node), pointer, intent(in) :: root
        integer(kind=int64), intent(in) :: key
        type(node), pointer :: res, tmp

        if (.not. associated(root)) then
            print *, "Erase failed"
            call stop_with_errorcode(1)
        end if

        if (key < root%key) then
            root%left => erase(root%left, key)
            res => root
        else if (key > root%key) then
            root%right => erase(root%right, key)
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
                    res => rotate_ccw(root)
                    res%left => erase(res%left, key)
                else
                    res => rotate_cw(root)
                    res%right => erase(res%right, key)
                end if
            end if
        end if
        if (associated(res)) call update(res)
    end function erase

    recursive function find_node(root, key) result(res)
        implicit none
        type(node), pointer, intent(in) :: root
        integer(kind=int64), intent(in) :: key
        type(node), pointer :: res
        if (.not. associated(root)) then
            res => null()
        else if (root%key == key) then
            res => root
        else if (key < root%key) then
            res => find_node(root%left, key)
        else
            res => find_node(root%right, key)
        end if
    end function find_node

    recursive function kth_node(root, k) result(res)
        implicit none
        type(node), pointer, intent(in) :: root
        integer(kind=int64), intent(in) :: k
        type(node), pointer :: res
        if (.not. associated(root)) then
            res => null()
        else if (k <= my_count(root%left)) then
            res => kth_node(root%left, k)
        else if (k == my_count(root%left) + 1) then
            res => root
        else
            res => kth_node(root%right, k - my_count(root%left) - 1)
        end if
    end function kth_node

    recursive subroutine delete_all(root)
        implicit none
        type(node), pointer, intent(inout) :: root
        if (.not. associated(root)) return

        call delete_all(root%left)
        call delete_all(root%right)
        deallocate (root)
        nullify (root)
    end subroutine delete_all

    recursive subroutine inorder(root, keys, vals, counter)
        implicit none
        type(node), pointer, intent(in) :: root
        integer(kind=int64), intent(inout) :: keys(:)
        integer(kind=int64), intent(inout) :: vals(:)
        integer(kind=int64), intent(inout) :: counter
        if (.not. associated(root)) return

        call inorder(root%left, keys, vals, counter)
        counter = counter + 1
        keys(counter) = root%key
        vals(counter) = root%val
        call inorder(root%right, keys, vals, counter)
    end subroutine inorder
end module module_dict
