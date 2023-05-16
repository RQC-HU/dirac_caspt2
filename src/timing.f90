!
SUBROUTINE timing(date0, tsec0, date, tsec)
!
    use module_global_variables, only: rank
    implicit none
    integer, intent(in)  :: date0
    real*8, intent(in)   :: tsec0
    integer, intent(inout) :: date
    real*8, intent(inout)  :: tsec
    real*8               :: difsec, sec, resd

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
