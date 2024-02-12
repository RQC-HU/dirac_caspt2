module module_time
    use module_global_variables, only: rank
    use module_error, only: stop_with_errorcode
    implicit none
    private
    public :: print_time, print_time_diff, get_current_time, get_current_time_and_print_diff, timing
    public :: time_type
    integer, parameter :: month_days(12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/), &
                          month_days_leap_year(12) = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    type :: time_type
        ! Structure for time
        ! Same as VALUES in DATE_AND_TIME but with more readable names
        ! Values are set by get_current_time subroutine
        integer :: year
        integer :: month
        integer :: date
        integer :: timezone
        integer :: hour
        integer :: min
        integer :: sec
        integer :: msec
    end type time_type

    type(time_type), public :: init_time, start_time, end_time
contains

    ! Public subroutines

    subroutine get_current_time(time)
        ! Get the current time and set the time to the input time
        implicit none
        class(time_type), intent(inout) :: time
        integer :: current_time(8)

        ! VALUES in DATE_AND_TIME stores the current time with 8 integers array
        ! https://gcc.gnu.org/onlinedocs/gfortran/DATE_005fAND_005fTIME.html
        call DATE_AND_TIME(VALUES=current_time)
        time%year = current_time(1)
        time%month = current_time(2)
        time%date = current_time(3)
        time%timezone = current_time(4)
        time%hour = current_time(5)
        time%min = current_time(6)
        time%sec = current_time(7)
        time%msec = current_time(8)
    end subroutine get_current_time

    subroutine print_time(time)
        ! Print time in the format of "YEAR = xxxx MONTH = xx DATE = xx HOUR = xx MIN = xx SEC = xx.xxx"
        implicit none
        class(time_type), intent(in) :: time

        if (rank == 0) then
            print '("Present time is")'
            print '("YEAR  = ",i0, " MONTH = ",i0, " DATE  = ",i0 )', time%year, time%month, time%date
            print '("HOUR  = ",i0, " MIN  = ",i0, " SEC  = ",i0,".",i3.3 )', time%hour, time%min, time%sec, time%msec
        end if
    end subroutine print_time

    subroutine get_current_time_and_print_diff(time_start, current_time)
        ! Get the current time and print the time difference from the input time_start
        implicit none
        class(time_type), intent(in) :: time_start
        class(time_type), intent(out) :: current_time

        call get_current_time(current_time)
        call print_time_diff(time_start, current_time)

    end subroutine get_current_time_and_print_diff

    subroutine print_time_diff(time_start, time_end)
        ! Print time difference in the format of "computational time = xxxx day xx h xx min xx.xxx sec"
        ! If time difference is minus, print warning message and print start and end time instead of time difference
        implicit none
        class(time_type), intent(in) :: time_start, time_end
        type(time_type) :: time_diff

        time_diff = calc_time_diff(time_start, time_end)
        if (rank == 0) then
            if (is_minus_time(time_diff)) then
                print *, "Warning: Time difference is minus. Print start and end time instead of time difference..."
                print '(A,6(i0,a),i3.3)', "start time: ", time_start%year, "-", time_start%month, "-", time_start%date, " ", &
                    time_start%hour, ":", time_start%min, ":", time_start%sec, ".", time_start%msec
                print '(A,6(i0,a),i3.3)', "end time: ", time_end%year, "-", time_end%month, "-", time_end%date, " ", &
                    time_end%hour, ":", time_end%min, ":", time_end%sec, ".", time_end%msec
            else
                print '("computational time = ",i0," day ",i0," h ",i0," min ",i0,".",i3.3," sec")', &
                    time_diff%date, time_diff%hour, time_diff%min, time_diff%sec, time_diff%msec
            end if
        end if
    end subroutine

    ! Internal (private) functions

    function calc_time_diff(time_start, time_end) result(time_diff)
        implicit none
        class(time_type), intent(in) :: time_start, time_end
        integer :: start_day, end_day
        type(time_type) :: time_diff

        ! Calculate number of days
        start_day = ymd2day(time_start)
        end_day = ymd2day(time_end)

        ! Initialize time_diff (year and month are included in the number of days)
        time_diff%year = 0; time_diff%month = 0  ! These are included in time_diff%date, thus set to 0
        time_diff%date = end_day - start_day  ! year + month + date
        time_diff%hour = time_end%hour - time_start%hour
        time_diff%min = time_end%min - time_start%min
        time_diff%sec = time_end%sec - time_start%sec
        time_diff%msec = time_end%msec - time_start%msec

        ! If time_diff%anything is minus, need to decrement the upper time unit by 1 and
        ! add the upper time unit to the lower time unit to adjust the time_diff
        ! (e.g. if time_diff%sec = -34, time_diff%min -= 1, time_diff%sec += 60 because 1 min = 60 sec)

        ! msec
        if (time_diff%msec < 0) then
            time_diff%sec = time_diff%sec - 1
            time_diff%msec = time_diff%msec + 1000
        end if
        ! sec
        if (time_diff%sec < 0) then
            time_diff%min = time_diff%min - 1
            time_diff%sec = time_diff%sec + 60
        end if
        ! min
        if (time_diff%min < 0) then
            time_diff%hour = time_diff%hour - 1
            time_diff%min = time_diff%min + 60
        end if
        ! hour
        if (time_diff%hour < 0) then
            time_diff%date = time_diff%date - 1
            time_diff%hour = time_diff%hour + 24
        end if

    end function calc_time_diff

    logical function is_minus_time(time)
        implicit none
        class(time_type), intent(in) :: time

        is_minus_time = .false.
        if (time%year < 0 .or. time%month < 0 .or. time%date < 0 &
            .or. time%hour < 0 .or. time%min < 0 .or. time%sec < 0 .or. time%msec < 0) then
            is_minus_time = .true.
        end if
    end function is_minus_time

    logical function is_leap_year(year)
        implicit none
        integer, intent(in) :: year

        is_leap_year = .false.
        if (mod(year, 400) == 0) then
            is_leap_year = .true.
        else if (mod(year, 100) == 0) then
            is_leap_year = .false.
        else if (mod(year, 4) == 0) then
            is_leap_year = .true.
        end if
    end function is_leap_year

    function ymd2day(time)
        ! Returns the number of days from 0001/01/01 to the input date - 1.
        ! -1 is needed because if time%date = 2, day xxxx/xx/02 has not yet ended.
        implicit none
        class(time_type), intent(in) :: time
        integer :: ymd2day, y, m, d

        y = year_to_day(time%year)
        m = month_to_day(time%month, is_leap_year(time%year))
        d = time%date - 1

        ymd2day = year_to_day(time%year) + month_to_day(time%month, is_leap_year(time%year)) + time%date - 1
    end function ymd2day

    function year_to_day(year)
        ! Returns the number of days from 0001/01/01 to the end of the input year - 1.
        ! -1 is needed because if year = 2, year 2 has not yet ended.
        implicit none
        integer, intent(in) :: year
        integer :: year_to_day, last_year

        if (year < 1) then
            if (rank == 0) print *, "Error: invalid year number, year must be year >= 1, but year = ", year
            call stop_with_errorcode(1)
        end if

        ! 365 days in a year, but 366 days for leap year
        last_year = year - 1
        year_to_day = 365*last_year + (last_year/4 - last_year/100 + last_year/400)
    end function year_to_day

    function month_to_day(month, leap_year)
        ! Returns the number of days from January 1st to the end of the input month - 1.
        ! -1 is needed because if month = 3, March has not yet ended.
        ! (e.g. if month = 3, then returns the number of days from January 1st to February 28th or 29th)
        implicit none
        integer, intent(in) :: month
        logical, intent(in) :: leap_year
        integer :: month_to_day

        if (month < 1 .or. month > 12) then
            if (rank == 0) print *, "Error: invalid month number, month must be 1 <= month <= 12, but month = ", month
            call stop_with_errorcode(1)
        end if
        if (month == 1) then ! January
            month_to_day = 0
        else if (leap_year) then
            month_to_day = sum(month_days_leap_year(1:month - 1))
        else
            month_to_day = sum(month_days(1:month - 1))
        end if
    end function month_to_day

    SUBROUTINE timing(date0, tsec0, date, tsec)
        implicit none
        integer, intent(in)  :: date0
        real(8), intent(in)   :: tsec0
        integer, intent(inout) :: date
        real(8), intent(inout)  :: tsec
        real(8)               :: difsec, sec, resd

        integer              ::  val(8), day, hour, min

        Call DATE_AND_TIME(VALUES=val)

        tsec = val(8)*(1.0d-03) + val(7) + val(6)*(6.0d+01) + val(5)*(6.0d+01)**2
        date = val(3)

        if (date0 < val(3)) then
            tsec = tsec + (val(3) - date0)*(6.0d+01)*(6.0d+01)*(2.4d+01)
        End if

        difsec = tsec - tsec0

        if (rank == 0) then
            print '("Present time is")'
            print '("year  = ",I4, "month = ",I4, "date  = ",I4 )', val(1), val(2), val(3)
            print '(14X, I4, "h   ",I4, "min  ",I2, ".",I3, "sec  " )',&
            & val(5), val(6), val(7), val(8)
        end if

        day = int(AINT(difsec)/(3600*24))
        resd = difsec - day*3600*24

        hour = int(AINT(resd)/3600)
        resd = resd - hour*3600

        min = int(AINT(resd)/60)
        resd = resd - min*60

        sec = resd
        if (rank == 0) then
            print '("computational time = ",I3, "day",I3, "h ",I3, &
            &"min",F7.3, "sec")', day, hour, min, sec
        end if
    end subroutine timing
end module module_time
