!
        SUBROUTINE timing(date0, tsec0, date, tsec)
!
            use four_caspt2_module, only: rank, normaloutput
            implicit none
            integer, intent(in)  :: date0
            real*8, intent(in)   :: tsec0
            integer, intent(inout) :: date
            real*8, intent(inout)  :: tsec
            real*8               :: difsec, sec, resd

            integer              ::  val(8), day, hour, min

            Call DATE_AND_TIME(VALUES=val)
!        Write(*,*)'Year = ',val(1),'Mon = ',val(2),'Date = ',val(3)
!        Write(*,*)'Hour = ',val(5),'Min = ',val(6),'Sec = ',val(7),'.',val(8)

            ! val(8): millisec, val(7): sec, val(6):min, val(5):hours
            tsec = val(8)*(1.0d-03) + val(7) + val(6)*(6.0d+01) + val(5)*(6.0d+01)**2

            if (date0 < val(3)) then
                tsec = tsec + (val(3) - date0)*(6.0d+01)*(6.0d+01)*(2.4d+01)
            End if

            difsec = tsec - tsec0

!        write(*,*)tsec, tsec0, difsec

            if (rank == 0) then ! Process limits for output
                write (normaloutput, '("Present time is")')
                write (normaloutput, '("year  = ",I4,"month = ",I4,"date  = ",I4 )') val(1), val(2), val(3)
                write (normaloutput, '(14X,I4,"h   ",I4,"min  ",I2,".",I3,"sec  " )')&
                & val(5), val(6), val(7), val(8)
            end if

!        If(AINT(difsec) > 3600*24) then
            day = AINT(difsec)/(3600*24)
            resd = difsec - day*3600*24
!        Else
!           day = 0
!           resd = difsec
!        Endif

!        If(AINT(resd) > 3600) then
            hour = AINT(resd)/3600
            resd = resd - hour*3600
!        Else
!           hour = 0
!           resd = resd
!        Endif

!        If(AINT(resd) > 60) then
            min = AINT(resd)/60
            resd = resd - min*60
!        Else
!           min  = 0
!           resd = resd
!        Endif

            sec = resd

!        write(*,'("computational time = ", F20.10,"sec")')difsec
            if (rank == 0) then ! Process limits for output
                write (normaloutput, '("computational time = ",I3,"day",I3,"h ",I3, &
                &"min",F7.3,"sec")') day, hour, min, sec
            end if
100     end subroutine timing
