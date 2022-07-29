! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE readint2_ord_co(filename) ! 2 electorn integrals created by typart in utchem

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
    character*50, intent(in) :: filename

    character  :: datex*10, timex*8

    integer :: mdcint, nkr, nmom, max1, max2, min1, min2
    integer :: nz
    integer :: i0, i, j, k, l
    integer :: SignIJ, SignKL, itr, jtr, ltr, ktr, inz, totalint

    integer, allocatable :: indk(:), indl(:), kr(:)

    real*8, allocatable :: rklr(:), rkli(:)

    !   Unit numbers for subspace files
    integer  :: unit_a1, unit_a2, unit_b, unit_c1, unit_c2, unit_c3, &
                unit_d1, unit_d2, unit_d3, unit_e, unit_f, unit_g, unit_h
    integer :: ioerr, iostat
    integer :: a1_cnt, a2_cnt, b_cnt, c1_cnt, c2_cnt, c3_cnt, d1_cnt, d2_cnt, d3_cnt, e_cnt, f_cnt, g_cnt, h_cnt
!Iwamuro modify
!        integer :: ikr, jkr, kkr, lkr
    !  Initialization of Unit numbers for subspace files
    unit_a1 = 100; unit_a2 = 200; unit_b = 300
    unit_c1 = 400; unit_c2 = 500; unit_c3 = 600
    unit_d1 = 700; unit_d2 = 800; unit_d3 = 900
    unit_e = 1000; unit_f = 1100; unit_g = 1200; unit_h = 1300
    a1_cnt = 0; a2_cnt = 0; b_cnt = 0; c1_cnt = 0; c2_cnt = 0; c3_cnt = 0
    d1_cnt = 0; d2_cnt = 0; d3_cnt = 0; e_cnt = 0; f_cnt = 0; g_cnt = 0; h_cnt = 0
    Allocate (kr(-nmo/2:nmo/2)); Call memplus(KIND(kr), SIZE(kr), 1)
    kr(:) = 0

    Allocate (indk((nmo/2)**2)); Call memplus(KIND(indk), SIZE(indk), 1)
    Allocate (indl((nmo/2)**2)); Call memplus(KIND(indl), SIZE(indl), 1)
    Allocate (rklr((nmo/2)**2)); Call memplus(KIND(rklr), SIZE(rklr), 1)
    Allocate (rkli((nmo/2)**2)); Call memplus(KIND(rkli), SIZE(rkli), 1)
    if (rank == 0) then
        print *, "enter readint2_ord_co"
        print '("Current Memory is ",F10.2,"MB")', tmem/1024/1024
    end if
    indk(:) = 0
    indl(:) = 0
    rklr(:) = 0.0d+00
    rkli(:) = 0.0d+00

    totalint = 0

    open (unit_a1, file=a1int, form='unformatted', status='replace')
    open (unit_a2, file=a2int, form='unformatted', status='replace')
    open (unit_b, file=bint, form='unformatted', status='replace')
    open (unit_c1, file=c1int, form='unformatted', status='replace')
    open (unit_c2, file=c2int, form='unformatted', status='replace')
    open (unit_c3, file=c3int, form='unformatted', status='replace')
    open (unit_d1, file=d1int, form='unformatted', status='replace')
    open (unit_d2, file=d2int, form='unformatted', status='replace')
    open (unit_d3, file=d3int, form='unformatted', status='replace')
    open (unit_e, file=eint, form='unformatted', status='replace')
    open (unit_f, file=fint, form='unformatted', status='replace')
    open (unit_g, file=gint, form='unformatted', status='replace')
    open (unit_h, file=hint, form='unformatted', status='replace')

    mdcint = 1500

    open (mdcint, file=trim(filename), form='unformatted', status='old', iostat=iostat)

    ! Check the status of the file
    if (iostat /= 0) then
        ! If iostat is not equal to 0, error detected in opening the file, so stop the program
        print *, 'ERROR: Failed to open '//trim(filename)//" , rank:", rank
        print *, 'Stop the program'
        stop
    end if

    Read (mdcint, iostat=iostat) datex, timex, nkr, &
        (kr(i0), kr(-1*i0), i0=1, nkr)

    ! Check the status of the file
    if (iostat < 0) then
        ! End of the file is reached. Return to the main program.
        print *, 'End of the file is reached  '//trim(filename)//" , rank:", rank
        print *, 'Return to the main program.'
        return
    elseif (iostat > 0) then
        ! If iostat is greater than 0, error detected in the input file, so exit the program
        print *, "Error : Error in reading file ", trim(filename)
        stop
    end if

    if (rank == 0) then
        print *, datex, timex
        print *, 'nkr', nkr, 'kr(+),kr(-)', (kr(i0), kr(-1*i0), i0=1, nkr)
    end if

    ! Continue to read the file until the end of the file is reached
    do
        read (mdcint, iostat=iostat) i, j, nz, &
            (indk(inz), indl(inz), inz=1, nz), &
            (rklr(inz), rkli(inz), inz=1, nz)
        ! Exit the loop if the end of the file is reached
        if (iostat < 0) then
            if (rank == 0) print *, 'End of '//trim(filename)
            exit
        elseif (iostat > 0) then
            ! If iostat is greater than 0, error detected in the input file, so exit the program
            print *, "Error : Error in reading file ", trim(filename)
            stop
        end if
        if (i == 0 .and. j == 0 .and. nz == 0) exit ! End of the file is reached, exit read loop

        totalint = totalint + nz

        itr = i + (-1)**(mod(i, 2) + 1)
        jtr = j + (-1)**(mod(j, 2) + 1)

        nmom = ninact + nact + nsec

        If (sp(i) == 4 .or. sp(j) == 4) cycle ! Read the next 2-integral
        If (i > ninact + nact .and. j > ninact + nact) cycle ! Read the next 2-integral

        SignIJ = (-1)**(mod(i + j, 2))

        Do inz = 1, nz

            k = indk(inz)
            ktr = k + (-1)**(mod(k, 2) + 1)
            l = indl(inz)
            ltr = l + (-1)**(mod(l, 2) + 1)

            If (sp(k) == 4 .or. sp(l) == 4) cycle ! Go to the next idz
            If (k > ninact + nact .and. l > ninact + nact) cycle ! Go to the next idz
            If (i == j .and. k > l) cycle ! Go to the next idz

            SignKL = (-1)**(mod(k + l, 2))

            max1 = max(sp(i), sp(j))
            min1 = min(sp(i), sp(j))
            max2 = max(sp(k), sp(l))
            min2 = min(sp(k), sp(l))

!===============================================================
! Integrals for A space  (pi|qr)(21|22) (pi|jk)(21|11)  type
!===============================================================

            If (max1 == 2 .and. min1 == 2 .and. max2 == 2 .and. min2 == 1) then    ! (22|21) => (21|22)

                if (k > l) then ! (22|21) => (21|22)
                    write (unit_a1) k, l, i, j, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_a1', ioerr, 'rank', rank
                    else

                        a1_cnt = a1_cnt + 1
                    end if

                else    ! (22|12) => (22|21)* => (21|22)*

                    write (unit_a1, IOSTAT=ioerr) l, k, j, i, rklr(inz), -1.0d+00*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_a1', ioerr, 'rank', rank
                    else

                        a1_cnt = a1_cnt + 1
                    end if
                end if

            elseif (max1 == 2 .and. min1 == 1 .and. max2 == 2 .and. min2 == 2) then ! (21|22) => (21|22)

                if (i > j) then ! (21|22) => (21|22)

                    write (unit_a1, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_a1', ioerr, 'rank', rank
                    else

                        a1_cnt = a1_cnt + 1
                    end if!                           write(*,'("A1int3",4I4,2E20.10)')i  ,j  ,k  ,l  ,  rklr(inz),         rkli(inz)

                else    ! (12|22) => (21|22)*

                    write (unit_a1, IOSTAT=ioerr) j, i, l, k, rklr(inz), -1.0d+00*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_a1', ioerr, 'rank', rank
                    else

                        a1_cnt = a1_cnt + 1
                    end if!                           write(*,'("A1int4",4I4,2E20.10)')j  ,i  ,l  ,k  ,  rklr(inz),-1.0d+00*rkli(inz)

                end if

            elseif (max1 == 2 .and. min1 == 1 .and. max2 == 1 .and. min2 == 1) then  ! (21|11)=>(21|11)

                if (i > j) then ! (21|11) => (21|11)

                    write (unit_a2, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_a2', ioerr, 'rank', rank
                    else

                        a2_cnt = a2_cnt + 1
                    end if!                           write(*,'("A2int1",4I4,2E20.10)')i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                else    ! (12|11) => (21|11)* => (21|11)*

                    write (unit_a2, IOSTAT=ioerr) j, i, l, k, rklr(inz), -1.0d+00*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_a2', ioerr, 'rank', rank
                    else

                        a2_cnt = a2_cnt + 1
                    end if!                           write(*,'("A2int2",4I4,2E20.10)')j  ,i  ,l  ,k  ,         rklr(inz), -1.0d+00*rkli(inz)
                end if

            elseif (max1 == 1 .and. min1 == 1 .and. max2 == 2 .and. min2 == 1) then  ! (11|21)=>(21|11)

                if (k > l) then   ! (11|21) => (21|11)

                    write (unit_a2, IOSTAT=ioerr) k, l, i, j, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_a2', ioerr, 'rank', rank
                    else

                        a2_cnt = a2_cnt + 1
                    end if!                           write(*,'("A2int3",4I4,2E20.10)')k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)

                else    ! (11|12) => (11|21)* => (21|11)*

                    write (unit_a2, IOSTAT=ioerr) l, k, j, i, rklr(inz), -1.0d+00*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_a2', ioerr, 'rank', rank
                    else

                        a2_cnt = a2_cnt + 1
                    end if!                           write(*,'("A2int4",4I4,2E20.10)')l  ,k  ,j  ,i  ,         rklr(inz), -1.0d+00*rkli(inz)

                end if
                ! end if

!=============================================
! Integrals for B space  (pi|qj) (21|21) type
!=============================================

            elseif (max1 == 2 .and. min1 == 1 .and. max2 == 2 .and. min2 == 1) then  ! (21|21)=>(21|21)

                if (i > j .and. k > l) then   ! (21|21) => (21|21)

                    write (unit_b, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_b', ioerr, 'rank', rank
                    else

                        b_cnt = b_cnt + 1
                    end if
                elseif (i < j .and. k > l) then ! (12|21) => (21|21)

                    write (unit_b, IOSTAT=ioerr) jtr, itr, k, l, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_b', ioerr, 'rank', rank
                    else

                        b_cnt = b_cnt + 1
                    end if
                elseif (i > j .and. k < l) then ! (21|12) => (21|21)

                    write (unit_b, IOSTAT=ioerr) i, j, ltr, ktr, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_b', ioerr, 'rank', rank
                    else

                        b_cnt = b_cnt + 1
                    end if
                elseif (i < j .and. k < l) then ! (12|12) => (21|21)*

                    write (unit_b, IOSTAT=ioerr) jtr, itr, ltr, ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_b', ioerr, 'rank', rank
                    else

                        b_cnt = b_cnt + 1
                    end if
                end if

!============================================================================
! Integrals for C space (ap|qr)(32|22) type C1int
!============================================================================

            elseif (max1 == 3 .and. min1 == 2 .and. max2 == 2 .and. min2 == 2) then ! (32|22)=>(32|22)

                if (i > j) then ! (32|22)=>(32|22)

                    write (unit_c1, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c1', ioerr, 'rank', rank
                    else

                        c1_cnt = c1_cnt + 1
                    end if!Iwamuro modify

                else    ! (23|22)=>(32|22)

                    write (unit_c1, IOSTAT=ioerr) jtr, itr, k, l, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c1', ioerr, 'rank', rank
                    else

                        c1_cnt = c1_cnt + 1
                    end if!Iwamuro modify
                end if

            elseif (max1 == 2 .and. min1 == 2 .and. max2 == 3 .and. min2 == 2) then ! (22|32)=>(32|22)

                if (k > l) then ! (22|32)=>(32|22)

                    write (unit_c1, IOSTAT=ioerr) k, l, i, j, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c1', ioerr, 'rank', rank
                    else

                        c1_cnt = c1_cnt + 1
                    end if!Iwamuro modify
                else    ! (22|23)=>(32|22)

                    write (unit_c1, IOSTAT=ioerr) ltr, ktr, i, j, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c1', ioerr, 'rank', rank
                    else

                        c1_cnt = c1_cnt + 1
                    end if!Iwamuro modify
                end if

!============================================================================
! Integrals for C space (ap|kk)(32|11)  type C2int
!============================================================================

            elseif (max1 == 3 .and. min1 == 2 .and. max2 == 1 .and. min2 == 1) then   ! (32|11)=>(32|11)

                if (i > j) then ! (32|11)=>(32|11)

                    write (unit_c2, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c2', ioerr, 'rank', rank
                    else

                        c2_cnt = c2_cnt + 1
                    end if
                else    ! (23|11)=>(32|11)

                    write (unit_c2, IOSTAT=ioerr) jtr, itr, k, l, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c2', ioerr, 'rank', rank
                    else

                        c2_cnt = c2_cnt + 1
                    end if
                end if

            elseif (max1 == 1 .and. min1 == 1 .and. max2 == 3 .and. min2 == 2) then   ! (32|11)=>(32|11)

                if (k > l) then ! (11|32)=>(32|11)

                    write (unit_c2, IOSTAT=ioerr) k, l, i, j, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c2', ioerr, 'rank', rank
                    else

                        c2_cnt = c2_cnt + 1
                    end if
                else    ! (11|23)=>(32|11)

                    write (unit_c2, IOSTAT=ioerr) ltr, ktr, i, j, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c2', ioerr, 'rank', rank
                    else

                        c2_cnt = c2_cnt + 1
                    end if
                end if

!============================================================================
! Integrals for C (ai|jp) (31|12)(C3int) and E space (ai|pj)(31|21) (Eint)
!============================================================================

            elseif (max1 == 3 .and. min1 == 1 .and. max2 == 2 .and. min2 == 1) then ! (31|21)=>(31|12)

                if (i > j .and. l > k) then ! (31|12)=>(31|21) For E

                    write (unit_e, IOSTAT=ioerr) i, j, ltr, ktr, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_e', ioerr, 'rank', rank
                    else

                        e_cnt = e_cnt + 1
                    end if
                    write (unit_c3, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c3', ioerr, 'rank', rank
                    else

                        c3_cnt = c3_cnt + 1
                    end if
                elseif (j > i .and. l > k) then ! (13|12)=>(31|21) For E

                    write (unit_e, IOSTAT=ioerr) jtr, itr, ltr, ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_e', ioerr, 'rank', rank
                    else

                        e_cnt = e_cnt + 1
                    end if
                    write (unit_c3, IOSTAT=ioerr) jtr, itr, k, l, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c3', ioerr, 'rank', rank
                    else

                        c3_cnt = c3_cnt + 1
                    end if
                elseif (i > j .and. k > l) then ! (31|21)=>(31|21) For E

                    write (unit_e, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_e', ioerr, 'rank', rank
                    else

                        e_cnt = e_cnt + 1
                    end if
                    write (unit_c3, IOSTAT=ioerr) i, j, ltr, ktr, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c3', ioerr, 'rank', rank
                    else

                        c3_cnt = c3_cnt + 1
                    end if
                elseif (i < j .and. k > l) then ! (13|21)=>(31|21) For E

                    write (unit_e, IOSTAT=ioerr) jtr, itr, k, l, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_e', ioerr, 'rank', rank
                    else

                        e_cnt = e_cnt + 1
                    end if
                    write (unit_c3, IOSTAT=ioerr) jtr, itr, ltr, ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c3', ioerr, 'rank', rank
                    else

                        c3_cnt = c3_cnt + 1
                    end if
                end if

            elseif (max1 == 2 .and. min1 == 1 .and. max2 == 3 .and. min2 == 1) then ! (21|31)=>(31|12)

                if (i > j .and. l > k) then ! (21|13)=>(31|21) For E

                    write (unit_e, IOSTAT=ioerr) ltr, ktr, i, j, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_e', ioerr, 'rank', rank
                    else

                        e_cnt = e_cnt + 1
                    end if
                    write (unit_c3, IOSTAT=ioerr) ltr, ktr, jtr, itr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c3', ioerr, 'rank', rank
                    else

                        c3_cnt = c3_cnt + 1
                    end if
                elseif (j > i .and. l > k) then ! (12|13)=>(31|21) For E

                    write (unit_e, IOSTAT=ioerr) ltr, ktr, jtr, itr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_e', ioerr, 'rank', rank
                    else

                        e_cnt = e_cnt + 1
                    end if
                    write (unit_c3, IOSTAT=ioerr) ltr, ktr, i, j, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c3', ioerr, 'rank', rank
                    else

                        c3_cnt = c3_cnt + 1
                    end if
                elseif (i > j .and. k > l) then ! (21|31)=>(31|21) For E

                    write (unit_e, IOSTAT=ioerr) k, l, i, j, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_e', ioerr, 'rank', rank
                    else

                        e_cnt = e_cnt + 1
                    end if
                    write (unit_c3, IOSTAT=ioerr) k, l, jtr, itr, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c3', ioerr, 'rank', rank
                    else

                        c3_cnt = c3_cnt + 1
                    end if
                elseif (i < j .and. k > l) then ! (12|31)=>(31|21) For E

                    write (unit_e, IOSTAT=ioerr) k, l, jtr, itr, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_e', ioerr, 'rank', rank
                    else

                        e_cnt = e_cnt + 1
                    end if
                    write (unit_c3, IOSTAT=ioerr) k, l, i, j, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_c3', ioerr, 'rank', rank
                    else

                        c3_cnt = c3_cnt + 1
                    end if
                end if

!============================================================================
! Integrals for D space (ai|pq)(31|22) type (D1int)
!============================================================================

            elseif (max1 == 3 .and. min1 == 1 .and. max2 == 2 .and. min2 == 2) then ! (31|22)=>(31|22)

                if (i > j) then ! (31|22)=>(31|22)

                    write (unit_d1, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then !ioerr /= 0
                        print *, 'error write unit_d1', ioerr, 'rank', rank
                    else

                        d1_cnt = d1_cnt + 1
                    end if
                else    ! (13|22)=>(31|22)

                    write (unit_d1, IOSTAT=ioerr) jtr, itr, k, l, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d1', ioerr, 'rank', rank
                    else

                        d1_cnt = d1_cnt + 1
                    end if
                end if

            elseif (max1 == 2 .and. min1 == 2 .and. max2 == 3 .and. min2 == 1) then ! (22|31)=>(31|22)

                if (k > l) then ! (22|31)=>(31|22)

                    write (unit_d1, IOSTAT=ioerr) k, l, i, j, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d1', ioerr, 'rank', rank
                    else

                        d1_cnt = d1_cnt + 1
                    end if
                else    ! (22|13)=>(31|22)

                    write (unit_d1, IOSTAT=ioerr) ltr, ktr, i, j, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d1', ioerr, 'rank', rank
                    else

                        d1_cnt = d1_cnt + 1
                    end if
                end if

!============================================================================
! Integrals for D space (ap|qi)(32|21) type (D2int)
!============================================================================

            elseif (max1 == 3 .and. min1 == 2 .and. max2 == 2 .and. min2 == 1) then ! (32|21)=>(32|21)

                if (i > j .and. k > l) then ! (32|21)=>(32|21)

                    write (unit_d2, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d2', ioerr, 'rank', rank
                    else

                        d2_cnt = d2_cnt + 1
                    end if
                elseif (i < j .and. k > l) then ! (23|21)=>(32|21)

                    write (unit_d2, IOSTAT=ioerr) jtr, itr, k, l, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d2', ioerr, 'rank', rank
                    else

                        d2_cnt = d2_cnt + 1
                    end if
                elseif (i > j .and. k < l) then ! (32|12)=>(32|21)

                    write (unit_d2, IOSTAT=ioerr) i, j, ltr, ktr, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d2', ioerr, 'rank', rank
                    else

                        d2_cnt = d2_cnt + 1
                    end if
                elseif (i < j .and. k < l) then ! (23|12)=>(32|21)

                    write (unit_d2, IOSTAT=ioerr) jtr, itr, ltr, ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d2', ioerr, 'rank', rank
                    else

                        d2_cnt = d2_cnt + 1
                    end if
                end if

            elseif (max1 == 2 .and. min1 == 1 .and. max2 == 3 .and. min2 == 2) then ! (21|32)=>(32|21)

                if (i > j .and. k > l) then ! (21|32)=>(32|21)

                    write (unit_d2, IOSTAT=ioerr) k, l, i, j, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d2', ioerr, 'rank', rank
                    else

                        d2_cnt = d2_cnt + 1
                    end if
                elseif (i < j .and. k > l) then ! (12|32)=>(32|21)

                    write (unit_d2, IOSTAT=ioerr) k, l, jtr, itr, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d2', ioerr, 'rank', rank
                    else

                        d2_cnt = d2_cnt + 1
                    end if
                elseif (i > j .and. k < l) then ! (21|23)=>(32|21)

                    write (unit_d2, IOSTAT=ioerr) ltr, ktr, i, j, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d2', ioerr, 'rank', rank
                    else

                        d2_cnt = d2_cnt + 1
                    end if
                elseif (i < j .and. k < l) then ! (12|23)=>(32|21)

                    write (unit_d2, IOSTAT=ioerr) ltr, ktr, jtr, itr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d2', ioerr, 'rank', rank
                    else

                        d2_cnt = d2_cnt + 1
                    end if
                end if

!============================================================================
! Integrals for D space (ai|jk)  (31|11) type (D3int)
!============================================================================

            elseif (max1 == 3 .and. min1 == 1 .and. max2 == 1 .and. min2 == 1) then   ! (31|11)=>(31|11)

                if (i > j) then ! (ai|jk) (31|11)=>(31|11)

                    write (unit_d3, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d3', ioerr, 'rank', rank
                    else

                        d3_cnt = d3_cnt + 1
                    end if
                else    ! (i~a~|kk) (13|11)=>(31|11)

                    write (unit_d3, IOSTAT=ioerr) jtr, itr, k, l, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d3', ioerr, 'rank', rank
                    else

                        d3_cnt = d3_cnt + 1
                    end if
                end if

            elseif (max1 == 1 .and. min1 == 1 .and. max2 == 3 .and. min2 == 1) then  ! (11|31)=>(31|11)

                if (k > l) then ! (jk|ai) (31|11)=>(31|11)

                    write (unit_d3, IOSTAT=ioerr) k, l, i, j, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d3', ioerr, 'rank', rank
                    else

                        d3_cnt = d3_cnt + 1
                    end if
                else  ! (jk|i~a~)=>( ai|kk) (11|13)=>(31|11)

                    write (unit_d3, IOSTAT=ioerr) ltr, ktr, i, j, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_d3', ioerr, 'rank', rank
                    else

                        d3_cnt = d3_cnt + 1
                    end if
                end if

!=============================================
! Integrals for F space  (ap|bq) (32|32) type
!=============================================

            elseif (max1 == 3 .and. min1 == 2 .and. max2 == 3 .and. min2 == 2) then  ! (32|32)=>(32|32)

                if (i > j .and. k > l) then   ! (32|32) => (32|32)

                    write (unit_f, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_f', ioerr, 'rank', rank
                    else

                        f_cnt = f_cnt + 1
                    end if
                elseif (i < j .and. k > l) then ! (23|32) => (32|32)

                    write (unit_f, IOSTAT=ioerr) jtr, itr, k, l, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_f', ioerr, 'rank', rank
                    else

                        f_cnt = f_cnt + 1
                    end if
                elseif (i > j .and. k < l) then ! (32|23) => (32|32)

                    write (unit_f, IOSTAT=ioerr) i, j, ltr, ktr, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_f', ioerr, 'rank', rank
                    else

                        f_cnt = f_cnt + 1
                    end if
                elseif (i < j .and. k < l) then ! (23|23) => (32|32)

                    write (unit_f, IOSTAT=ioerr) jtr, itr, ltr, ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_f', ioerr, 'rank', rank
                    else

                        f_cnt = f_cnt + 1
                    end if
                end if

!============================================================================
! G space (ai|bp)(31|32) type
!============================================================================

            elseif (max1 == 3 .and. min1 == 1 .and. max2 == 3 .and. min2 == 2) then ! (31|32)=>(31|32)

                if (i > j .and. l > k) then ! (31|23)=>(31|32)

                    write (unit_g, IOSTAT=ioerr) i, j, ltr, ktr, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_g', ioerr, 'rank', rank
                    else

                        g_cnt = g_cnt + 1
                    end if!                           write(*,'("Gint1",4I4,2E20.10)')i  ,j  ,ltr,ktr,   SignKL*rklr(inz),  SignKL*rkli(inz)

                elseif (j > i .and. l > k) then ! (13|23)=>(31|32)

                    write (unit_g, IOSTAT=ioerr) jtr, itr, ltr, ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_g', ioerr, 'rank', rank
                    else

                        g_cnt = g_cnt + 1
                    end if!                           write(*,'("Gint2",4I4,2E20.10)')jtr,itr,ltr,ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)

                elseif (i > j .and. k > l) then ! (31|32)=>(31|32)

                    write (unit_g, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_g', ioerr, 'rank', rank
                    else

                        g_cnt = g_cnt + 1
                    end if!                           write(*,'("Gint3",4I4,2E20.10)')i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                elseif (i < j .and. k > l) then ! (13|32)=>(31|32)

                    write (unit_g, IOSTAT=ioerr) jtr, itr, k, l, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_g', ioerr, 'rank', rank
                    else

                        g_cnt = g_cnt + 1
                    end if!                           write(*,'("Gint4",4I4,2E20.10)')jtr,itr,k  ,l  ,  SignIJ*rklr(inz), SignIJ*rkli(inz)

                end if

            elseif (max1 == 3 .and. min1 == 2 .and. max2 == 3 .and. min2 == 1) then ! (32|31)=>(31|32)

                if (i > j .and. l > k) then ! (32|13)=>(31|32)

                    write (unit_g, IOSTAT=ioerr) ltr, ktr, i, j, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_g', ioerr, 'rank', rank
                    else

                        g_cnt = g_cnt + 1
                    end if!                           write(*,'("Gint5",4I4,2E20.10)')ltr,ktr,i  ,j  ,        SignKL*rklr(inz),        SignKL*rkli(inz)

                elseif (j > i .and. l > k) then ! (23|13)=>(31|32)

                    write (unit_g, IOSTAT=ioerr) ltr, ktr, jtr, itr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_g', ioerr, 'rank', rank
                    else

                        g_cnt = g_cnt + 1
                    end if                   !    print '("Gint6",4I4,2E20.10)', ltr, ktr, jtr, itr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)

                elseif (i > j .and. k > l) then ! (32|31)=>(31|32)

                    write (unit_g, IOSTAT=ioerr) k, l, i, j, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_g', ioerr, 'rank', rank
                    else

                        g_cnt = g_cnt + 1
                    end if!                           write(*,'("Gint7",4I4,2E20.10)')k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)

                elseif (i < j .and. k > l) then ! (23|31)=>(31|32)

                    write (unit_g, IOSTAT=ioerr) k, l, jtr, itr, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_g', ioerr, 'rank', rank
                    else

                        g_cnt = g_cnt + 1
                    end if!                           write(*,'("Gint8",4I4,2E20.10)')k  ,l  ,jtr,itr,  SignIJ*rklr(inz),  SignIJ*rkli(inz)

                end if

!=============================================
! Integrals for H space  (ai|bj) (31|31) type
!=============================================

            elseif (max1 == 3 .and. min1 == 1 .and. max2 == 3 .and. min2 == 1) then  ! (31|31)=>(31|31)

                if (i > j .and. k > l) then   ! (31|31) => (31|31)

                    write (unit_h, IOSTAT=ioerr) i, j, k, l, rklr(inz), rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_h', ioerr, 'rank', rank
                    else

                        h_cnt = h_cnt + 1
                    end if!                           write(*,'("Hint1",4I4,2E20.10)')i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                elseif (i < j .and. k > l) then ! (13|31) => (31|31)

                    write (unit_h, IOSTAT=ioerr) jtr, itr, k, l, SignIJ*rklr(inz), SignIJ*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_h', ioerr, 'rank', rank
                    else

                        h_cnt = h_cnt + 1
                    end if!                           write(*,'("Hint2",4I4,2E20.10)')jtr,itr,k  ,l  ,  SignIJ*rklr(inz),  SignIJ*rkli(inz)

                elseif (i > j .and. k < l) then ! (31|13) => (31|31)

                    write (unit_h, IOSTAT=ioerr) i, j, ltr, ktr, SignKL*rklr(inz), SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_h', ioerr, 'rank', rank
                    else

                        h_cnt = h_cnt + 1
                    end if!                           write(*,'("Hint3",4I4,2E20.10)')i  ,j  ,ltr,ktr,  SignKL*rklr(inz),  SignKL*rkli(inz)

                elseif (i < j .and. k < l) then ! (13|13) => (31|31)

                    write (unit_h, IOSTAT=ioerr) jtr, itr, ltr, ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                    if (ioerr .ne. 0) then
                        print *, 'error write unit_h', ioerr, 'rank', rank
                    else

                        h_cnt = h_cnt + 1
                    end if!                           write(*,'("Hint4",4I4,2E20.10)')jtr,itr,ltr,ktr,SignIJ*SignKL*rklr(inz),SignIJ*SignKL*rkli(inz)
                end if
            end if

            if (abs(rkli(inz)) > thres) realc = .false.
        end do ! Next inz
    end do ! Continue to read 2-integrals

10  if (rank == 0) print *, 'error for opening mdcint 10'
    go to 100

100 continue

    close (mdcint)
    close (unit_a1)
    close (unit_a2)
    close (unit_b)
    close (unit_c1)
    close (unit_c2)
    close (unit_c3)
    close (unit_d1)
    close (unit_d2)
    close (unit_d3)
    close (unit_e)
    close (unit_f)
    close (unit_g)
    close (unit_h)
    if (allocated(indk)) deallocate (indk); Call memminus(KIND(indk), SIZE(indk), 1)
    if (allocated(indl)) deallocate (indl); Call memminus(KIND(indl), SIZE(indl), 1)
    if (allocated(rklr)) deallocate (rklr); Call memminus(KIND(rklr), SIZE(rklr), 1)
    if (allocated(rkli)) deallocate (rkli); Call memminus(KIND(rkli), SIZE(rkli), 1)
    if (allocated(kr)) deallocate (kr); Call memminus(KIND(kr), SIZE(kr), 1)
    if (rank == 0) print *, "end readint2_ord_co"

end subroutine readint2_ord_co
