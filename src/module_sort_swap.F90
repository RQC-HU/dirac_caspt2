module module_sort_swap
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
! module_sort_swap
! Copyright (c) by the authors of DIRAC-CASPT2.
! Author K.Noda
!
! This is a utility module that supports you to sort list or swap values or swap values in list.
!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
    implicit none
    private
    public heapSort, swap
    interface heapSort
        module procedure heapSortInt, heapSortReal
    end interface heapSort
    interface heapify
        module procedure heapifyInt, heapifyReal
    end interface heapify
    interface swap
        module procedure swapInt, swapReal, swapCmp16, swapArrayInt, swapArrayReal, swapArrayCmp16
    end interface swap
contains
    subroutine heapifyInt(array, first, last, is_descending_order)
        implicit none
        integer, INTENT(INOUT) :: array(:)
        integer, INTENT(IN) :: first, last
        logical, intent(in) :: is_descending_order
        integer ::  parent, child
        parent = first
        child = 2*parent
        if (is_descending_order) then
            do while (child <= last)
                if (child < last) then
                    if (array(child) > array(child + 1)) child = child + 1
                end if
                if (array(child) >= array(parent)) exit
                call swap(array, child, parent)
                parent = child
                child = 2*parent
            end do
        else
            do while (child <= last)
                if (child < last) then
                    if (array(child) < array(child + 1)) child = child + 1
                end if
                if (array(child) <= array(parent)) exit
                call swap(array, child, parent)
                parent = child
                child = 2*parent
            end do
        end if
    end subroutine heapifyInt
    subroutine heapSortInt(list, is_descending_order)
        implicit none
        integer, INTENT(INOUT) :: list(:)
        logical, intent(in) :: is_descending_order
        integer :: array_size, i, idx
        i = 1
        array_size = size(list)
        ! Build heap
        do idx = array_size/2, 1, -1
            call heapify(list, idx, array_size, is_descending_order)
        end do
        do idx = array_size, 2, -1
            call swap(list, 1, idx)
            call heapify(list, 1, idx - 1, is_descending_order)
        end do
    end subroutine heapSortInt

    subroutine heapifyReal(array, first, last, is_descending_order)
        implicit none
        real(8), INTENT(INOUT) :: array(:)
        integer, INTENT(IN) :: first, last
        logical, intent(in) :: is_descending_order
        integer ::  parent, child
        parent = first
        child = 2*parent
        if (is_descending_order) then
            do while (child <= last)
                if (child < last) then
                    if (array(child) > array(child + 1)) child = child + 1
                end if
                if (array(child) >= array(parent)) exit
                call swap(array, child, parent)
                parent = child
                child = 2*parent
            end do
        else
            do while (child <= last)
                if (child < last) then
                    if (array(child) < array(child + 1)) child = child + 1
                end if
                if (array(child) <= array(parent)) exit
                call swap(array, child, parent)
                parent = child
                child = 2*parent
            end do
        end if
    end subroutine heapifyReal
    subroutine heapSortReal(list, is_descending_order)
        implicit none
        real(8), INTENT(INOUT) :: list(:)
        logical, intent(in) :: is_descending_order
        integer :: array_size, i, idx
        i = 1
        array_size = size(list)
        ! Build heap
        do idx = array_size/2, 1, -1
            call heapify(list, idx, array_size, is_descending_order)
        end do
        do idx = array_size, 2, -1
            call swap(list, 1, idx)
            call heapify(list, 1, idx - 1, is_descending_order)
        end do
    end subroutine

    subroutine swapInt(a, b)
        ! Swap values between a and b
        implicit none
        integer temp
        integer, INTENT(INOUT) :: a, b
        temp = a
        a = b
        b = temp
    end subroutine swapInt
    subroutine swapReal(a, b)
        ! Swap values between a and b
        implicit none
        real(8) temp
        real(8), INTENT(INOUT) :: a, b
        temp = a
        a = b
        b = temp
    end subroutine swapReal
    subroutine swapCmp16(a, b)
        ! Swap values between a and b
        implicit none
        complex*16 temp
        complex*16, INTENT(INOUT) :: a, b
        temp = a
        a = b
        b = temp
    end subroutine swapCmp16

    subroutine swapArrayInt(array, a, b)
        ! Swap values between array(a) and array(b)
        implicit none
        integer temp
        integer, INTENT(INOUT) :: array(:)
        integer, INTENT(IN) :: a, b
        temp = array(a)
        array(a) = array(b)
        array(b) = temp
    end subroutine swapArrayInt
    subroutine swapArrayReal(array, a, b)
        ! Swap values between array(a) and array(b)
        implicit none
        real(8) temp
        real(8), INTENT(INOUT) :: array(:)
        integer, INTENT(IN) :: a, b
        temp = array(a)
        array(a) = array(b)
        array(b) = temp
    end subroutine swapArrayReal
    subroutine swapArrayCmp16(array, a, b)
        ! Swap values between array(a) and array(b)
        implicit none
        complex*16 temp
        complex*16, INTENT(INOUT) :: array(:)
        integer, INTENT(IN) :: a, b
        temp = array(a)
        array(a) = array(b)
        array(b) = temp
    end subroutine swapArrayCmp16
end module module_sort_swap
