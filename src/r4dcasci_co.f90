! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM r4dcasci_co   ! DO CASCI CALC IN THIS PROGRAM!

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module
    use read_input_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    integer                 :: i0, nuniq
    integer                 ::  inisym, endsym
    real(16)                :: time0, time1

!        integer                 ::  val(8), initdate, date0, date1
!        real*8                  :: totalsec, inittime, tsec0, tsec1, tsec

    logical                 :: test
    character*50            :: filename

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!        debug = .TRUE.
    debug = .FALSE.
    thres = 1.0d-15
!        thres = 0.0d+00

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   MPI initialization and get the number of MPI processes (nprocs) and own process number.
#ifdef HAVE_MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_rank(MPI_COMM_WORLD, rank, ierr)
#else
    rank = 0; nprocs = 1
#endif
    if (rank == 0) then ! Process limits for output
        ! open (normal_output, file='caspt2.out', form='formatted', status='unknown')
        write (*, '(A,I8,A,I8)') 'initialization of mpi, rank :', rank, ' nprocs :', nprocs
        write (*, *) ''
        write (*, *) ' ENTER R4DCASCI_TY PROGRAM written by M. Abe 2007.7.19'
        write (*, *) ''
    end if
    tmem = 0.0d+00

    if (rank == 0) then ! Process limits for output
        write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024

        val = 0
        Call DATE_AND_TIME(VALUES=val)

        write (*, *) 'Year = ', val(1), 'Mon = ', val(2), 'Date = ', val(3)
        write (*, *) 'Hour = ', val(5), 'Min = ', val(6), 'Sec = ', val(7), '.', val(8)

        totalsec = val(8)*(1.0d-03) + val(7) + val(6)*(6.0d+01) + val(5)*(6.0d+01)**2
        initdate = val(3)
        inittime = totalsec

        write (*, *) inittime
        ! Call timing(val(3), totalsec, date0, tsec)
    end if
    call read_input

    if (rank == 0) then ! Process limits for output
        write (*, *) 'ninact        =', ninact
        write (*, *) 'nact          =', nact
        write (*, *) 'nsec          =', nsec
        write (*, *) 'nelec         =', nelec
        write (*, *) 'nroot         =', nroot
        write (*, *) 'selectroot    =', selectroot
        write (*, *) 'totsym        =', totsym
        write (*, *) 'ncore         =', ncore
        write (*, *) 'nbas          =', nbas
        write (*, *) 'eshift        =', eshift
        write (*, *) 'ptgrp         =', ptgrp
        write (*, *) 'dirac_version =', dirac_version
        print *, "RAS1 = ", ras1_list
        print *, "RAS3 = ", ras3_list
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    filename = 'MRCONEE'

    call readorb_enesym_co(filename)
    call read1mo_co(filename)

    if (rank == 0) then ! Process limits for output
        write (*, *) 'realc', realc, ECORE, ninact, nact, nsec, nmo
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_mdcint_filename

    !Iwamuro create new ikr for dirac
    Call create_newmdcint
    if (rank == 0) then ! Process limits for output
        write (*, *) 'Before readint2_casci_co', rank
    end if
    ! if (rank == 0) then
    filename = 'MDCINTNEW'

    ! Call readint2_casci_co(filename, nuniq)
    Call readint2_casci_co(mdcintnew, nuniq)

!        Allocate(sp(1:nmo)) ;  Call memplus(KIND(sp),SIZE(sp),1)
!        sp( 1               : ninact           )    = 1
!        sp( ninact+1        : ninact+nact      )    = 2
!        sp( ninact+nact+1   : ninact+nact+nsec )    = 3
!        sp( ninact+nact+nsec: nmo              )    = 4
    if (rank == 0) then ! Process limits for output
        write (*, *) 'nmo        =', nmo
    end if
    nmo = ninact + nact + nsec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (rank == 0) then  ! Process limits for output
        write (*, *) "iwamuro modify"
    end if
    If (mod(nelec, 2) == 0) then
        inisym = nsymrp + 1
        endsym = 2*nsymrp
    Else
        inisym = 1
        endsym = nsymrp
    End if

    if (rank == 0) then ! Process limits for output
        write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024
!   Do totsym = inisym, inisym
!   Do totsym = inisym, endsym

!      totsym = 4

        write (*, *) ' '
        write (*, *) '*******************************'
        write (*, *) ' '
        write (*, *) 'IREP IS ', repna(totsym)
        write (*, *) ' '
        write (*, *) '*******************************'
        write (*, *) ' '
    end if
    realcvec = .TRUE.

    Call casci_ty

!      goto 1000

!    This is test for bug fix about realc part
    if (rank == 0) then ! Process limits for output
        write (*, *) realc, 'realc'
        write (*, *) realcvec, 'realcvec'
    end if
    test = .true.

    if (rank == 0) then ! Process limits for output
        write (*, *) realc, 'realc'
        write (*, *) realcvec, 'realcvec'
    end if
    realc = .FALSE.      !!!      realc =.TRUE.
    realcvec = .FALSE.   !!!      realcvec =.TRUE.

    if (rank == 0) then ! Process limits for output
        write (*, *) 'FOR TEST WE DO (F,F)'
        write (*, *) realc, 'realc'
        write (*, *) realcvec, 'realcvec'
    end if
!!=============================================!
!                                              !
    iroot = selectroot
!                                              !
!!=============================================!

    Call e0test_v2

    if (rank == 0) then ! Process limits for output
        write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024
    end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!            BUILDING  FOCK MATRIX               !
!  fij = hij + SIGUMA[<0|Ekl|0>{(ij|kl)-(il|kj)} !
!                 kl                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!! TEST TO CALCULATE FOCK MATRIX OF HF STATE fpq = hpq + SIGUMA_r[(pq|rr)-(pr|qr)]
!! THIS MUST BE DIAGONAL MATRIX AND DIAGONAL ELEMENTS CORESPONDS TO SPINOR ENERGIES.

    Allocate (f(nmo, nmo)); Call memplus(KIND(f), SIZE(f), 2)

!      debug = .FALSE.
    debug = .TRUE.
    If (debug) then
        f(:, :) = 0.0d+00
        if (rank == 0) then ! Process limits for output
            write (*, *) 'fockhf1_ty start'
        end if
        Call fockhf1_ty
    End if

!! NOW MAKE FOCK MATRIX FOR CASCI STATE
!! fij = hij + SIGUMA_kl[<0|Ekl|0>{(ij|kl)-(il|kj)}

    f(:, :) = 0.0d+00
    if (rank == 0) then
        write (*, *) 'before building fock'
        date1 = date0
        tsec1 = tsec0
        Call timing(date1, tsec1, date0, tsec0)
    end if
    Call fockcasci_ty
    if (rank == 0) then
        write (*, *) 'end building fock'
        date1 = date0
        tsec1 = tsec0
        Call timing(date1, tsec1, date0, tsec0)
    end if
!      debug = .TRUE.
    debug = .FALSE.
    if (rank == 0) then ! Process limits for output
        write (*, *) debug, 'debug'
    end if
    if (debug) Call prtoutfock

    Allocate (eps(nmo)); Call memplus(KIND(eps), SIZE(eps), 1)
    eps = 0.0d+00

    Call fockdiag_ty

    if (rank == 0) then ! Process limits for output
        Do i0 = 1, nmo
            write (*, *) 'eps(', i0, ')=', eps(i0)
        End do
    end if
!      Do i0 = 1, nmo/2
!         if(ABS(eps(i0*2)-eps(i0*2-1)) > 1.0d-10) then
!            write(*,*)i0*2-1,i0*2,eps(i0*2-1),eps(i0*2)
!         Endif
!      Enddo
    if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        open (5, file='EPS', form='unformatted', status='unknown')
        write (5) nmo
        write (5) eps(1:nmo)
        close (5)
    end if
    ! end if
    deallocate (sp); Call memplus(KIND(sp), SIZE(sp), 1)
    deallocate (cir); Call memminus(KIND(cir), SIZE(cir), 1)
    deallocate (cii); Call memminus(KIND(cii), SIZE(cii), 1)
    deallocate (eigen); Call memminus(KIND(eigen), SIZE(eigen), 1)
    deallocate (f); Call memminus(KIND(f), SIZE(f), 2)
    deallocate (eps); Call memminus(KIND(eps), SIZE(eps), 1)
    deallocate (idet); Call memminus(KIND(idet), SIZE(idet), 1)
    deallocate (idetr); Call memminus(KIND(idetr), SIZE(idetr), 1)
    deallocate (MULTB_S); Call memminus(KIND(MULTB_S), SIZE(MULTB_S), 1)
    deallocate (MULTB_D); Call memminus(KIND(MULTB_D), SIZE(MULTB_D), 1)
    deallocate (MULTB_DS); Call memminus(KIND(MULTB_DS), SIZE(MULTB_DS), 1)
    deallocate (MULTB_DF); Call memminus(KIND(MULTB_DF), SIZE(MULTB_DF), 1)
    deallocate (MULTB_DB); Call memminus(KIND(MULTB_DB), SIZE(MULTB_DB), 1)
    deallocate (MULTB_SB); Call memminus(KIND(MULTB_SB), SIZE(MULTB_SB), 1)

    deallocate (orb); Call memminus(KIND(orb), SIZE(orb), 1)
    deallocate (irpmo); Call memminus(KIND(irpmo), SIZE(irpmo), 1)
    deallocate (irpamo); Call memminus(KIND(irpamo), SIZE(irpamo), 1)
    deallocate (indmo); Call memminus(KIND(indmo), SIZE(indmo), 1)
    deallocate (indmor); Call memminus(KIND(indmor), SIZE(indmor), 1)
    deallocate (onei); Call memminus(KIND(onei), SIZE(onei), 1)
    ! deallocate (int2i); Call memminus(KIND(int2i), SIZE(int2i), 1)
    deallocate (inttwi); Call memminus(KIND(inttwi), SIZE(inttwi), 1)
    ! deallocate (indtwi); Call memminus(KIND(indtwi), SIZE(indtwi), 1)
    deallocate (oner); Call memminus(KIND(oner), SIZE(oner), 1)
    ! deallocate (int2r); Call memminus(KIND(int2r), SIZE(int2r), 1)
    deallocate (inttwr); Call memminus(KIND(inttwr), SIZE(inttwr), 1)
    ! deallocate (indtwr); Call memminus(KIND(indtwr), SIZE(indtwr), 1)
    deallocate (int2r_f1); Call memminus(KIND(int2r_f1), SIZE(int2r_f1), 1)
    deallocate (int2i_f1); Call memminus(KIND(int2i_f1), SIZE(int2i_f1), 1)
    deallocate (int2r_f2); Call memminus(KIND(int2r_f2), SIZE(int2r_f2), 1)
    deallocate (int2i_f2); Call memminus(KIND(int2i_f2), SIZE(int2i_f2), 1)
    if (rank == 0) then ! Process limits for output
        write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024

        Call timing(val(3), totalsec, date0, tsec0)
        write (*, *) 'End r4dcasci_ty part'
    end if
#ifdef HAVE_MPI
    call MPI_FINALIZE(ierr)
#endif
    if (rank == 0) then ! Process limits for output
        write (*, '(a,i4,a,i4)') 'fin. rank:', rank, 'nprocs:', nprocs
    end if
1000 continue
END program r4dcasci_co
