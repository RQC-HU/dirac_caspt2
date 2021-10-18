! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM r4dcasci_co   ! DO CASCI CALC IN THIS PROGRAM!

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
    include 'mpif.h'
    integer                 :: ii, jj, kk, ll, typetype, i0, j0
    integer                 ::  j, i, k, l, nuniq
    integer                 :: k0, l0, nint, n, dimn, n0, n1, nspace(3, 3)
    integer                 ::  totsym, inisym, endsym

!        integer                 ::  val(8), initdate, date0, date1
!        real*8                  :: totalsec, inittime, tsec0, tsec1, tsec

    logical                 :: test, cutoff

    real*8                  :: i2r, i2i, dr, di, nsign, e0, e2, e2all
    complex*16              ::  cmplxint, dens, trace1, trace2, dens1, dens2

    character*50            :: filename

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!        debug = .TRUE.
    debug = .FALSE.
    thres = 1.0d-15
!        thres = 0.0d+00

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   MPI initialization and get the number of MPI processes (nprocs) and own process number.
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_rank(MPI_COMM_WORLD, rank, ierr)
    if (rank == 0) then ! Process limits for output
        open (normaloutput, file='caspt2.out', form='formatted', status='unknown')
        write (normaloutput, '(A,I8,A,I8)') 'initialization of mpi, rank :', rank, ' nprocs :', nprocs
        write (normaloutput, *) ''
        write (normaloutput, *) ' ENTER R4DCASCI_TY PROGRAM written by M. Abe 2007.7.19'
        write (normaloutput, *) ''
    end if
    tmem = 0.0d+00

    if (rank == 0) then ! Process limits for output
        write (normaloutput, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024

        val = 0
        Call DATE_AND_TIME(VALUES=val)

        write (normaloutput, *) 'Year = ', val(1), 'Mon = ', val(2), 'Date = ', val(3)
        write (normaloutput, *) 'Hour = ', val(5), 'Min = ', val(6), 'Sec = ', val(7), '.', val(8)

        totalsec = val(8)*(1.0d-03) + val(7) + val(6)*(6.0d+01) + val(5)*(6.0d+01)**2
        initdate = val(3)
        inittime = totalsec

        write (normaloutput, *) inittime
        ! Call timing(val(3), totalsec, date0, tsec)
    end if
!     end if
    open (5 + rank, file='active.inp', form='formatted', status='old')
    read (5 + rank, '(I4)') ninact
    read (5 + rank, '(I4)') nact
    read (5 + rank, '(I4)') nsec
    read (5 + rank, '(I4)') nelec
    read (5 + rank, '(I4)') nroot
    read (5 + rank, '(I4)') selectroot
    read (5 + rank, '(I4)') totsym
    read (5 + rank, '(I4)') ncore
    read (5 + rank, '(I4)') nbas
    read (5 + rank, '(E8.2)') eshift
    read (5 + rank, '(A6)') ptgrp
    close (5 + rank)

    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'ninact     =', ninact
        write (normaloutput, *) 'nact       =', nact
        write (normaloutput, *) 'nsec       =', nsec
        write (normaloutput, *) 'nelec      =', nelec
        write (normaloutput, *) 'nroot      =', nroot
        write (normaloutput, *) 'selectroot =', selectroot
        write (normaloutput, *) 'totsym     =', totsym
        write (normaloutput, *) 'ncore      =', ncore
        write (normaloutput, *) 'nbas       =', nbas
        write (normaloutput, *) 'eshift     =', eshift
        write (normaloutput, *) 'ptgrp      =', ptgrp
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    filename = 'MRCONEE'

    call readorb_enesym_co(filename)
    call read1mo_co(filename)

    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'realc', realc, ECORE, ninact, nact, nsec, nmo
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_mdcint_filename

    !Iwamuro create new ikr for dirac
    Call create_newmdcint
    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'Before readint2_casci_co', rank
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
        write (normaloutput, *) 'nmo        =', nmo
    end if
    nmo = ninact + nact + nsec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (rank == 0) then  ! Process limits for output
        write (normaloutput, *) "iwamuro modify"
    end if
    If (mod(nelec, 2) == 0) then
        inisym = nsymrp + 1
        endsym = 2*nsymrp
    Else
        inisym = 1
        endsym = nsymrp
    End if

    if (rank == 0) then ! Process limits for output
        write (normaloutput, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024
!   Do totsym = inisym, inisym
!   Do totsym = inisym, endsym

!      totsym = 4

        write (normaloutput, *) ' '
        write (normaloutput, *) '*******************************'
        write (normaloutput, *) ' '
        write (normaloutput, *) 'IREP IS ', repna(totsym)
        write (normaloutput, *) ' '
        write (normaloutput, *) '*******************************'
        write (normaloutput, *) ' '
    end if
    realcvec = .TRUE.

    Call casci_ty(totsym)

!      goto 1000

!    This is test for bug fix about realc part
    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) realc, 'realc'
        write (normaloutput, *) realcvec, 'realcvec'
    end if
    test = .true.

    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) realc, 'realc'
        write (normaloutput, *) realcvec, 'realcvec'
    end if
    realc = .FALSE.      !!!      realc =.TRUE.
    realcvec = .FALSE.   !!!      realcvec =.TRUE.

    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) 'FOR TEST WE DO (F,F)'
        write (normaloutput, *) realc, 'realc'
        write (normaloutput, *) realcvec, 'realcvec'
    end if
!!=============================================!
!                                              !
    iroot = selectroot
!                                              !
!!=============================================!

    Call e0test_v2

    if (rank == 0) then ! Process limits for output
        write (normaloutput, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024
    end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!            BUILDING  FOCK MATRIX               !
!  fij = hij + SIGUMA[<0|Ekl|0>{(ij|kl)-(il|kj)} !
!                 kl                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!! TEST TO CALCULATE FOCK MATRIX OF HF STATE fpq = hpq + SIGUMA_r[(pq|rr)-(pr|qr)]
!! THIS MUST BE DIAGONAL MATRIX AND DIAGONAL ELEMENTS CORESPONDS TO SPINOR ENERGIES.

    Allocate (f(nmo, nmo)); Call memplus(KIND(f), SIZE(f), 2)

    f(:, :) = 0.0d+00

!      debug = .FALSE.
    debug = .TRUE.
    If (debug) then
        if (rank == 0) then ! Process limits for output
            write (normaloutput, *) 'fockhf1_ty start'
        end if
        Call fockhf1_ty
    End if

!! NOW MAKE FOCK MATRIX FOR CASCI STATE
!! fij = hij + SIGUMA_kl[<0|Ekl|0>{(ij|kl)-(il|kj)}

    f(:, :) = 0.0d+00
    if (rank == 0) then
        write(normaloutput,*)'before building fock'
        date1 = date0
        tsec1 = tsec0
        Call timing(date1, tsec1, date0, tsec0)
    end if
    Call fockcasci_ty
    if (rank == 0) then
        write(normaloutput,*)'end building fock'
        date1 = date0
        tsec1 = tsec0
        Call timing(date1, tsec1, date0, tsec0)
    end if
!      debug = .TRUE.
    debug = .FALSE.
    if (rank == 0) then ! Process limits for output
        write (normaloutput, *) debug, 'debug'
    end if
    if (debug) Call prtoutfock

    Allocate (eps(nmo)); Call memplus(KIND(eps), SIZE(eps), 1)
    eps = 0.0d+00

    Call fockdiag_ty

    if (rank == 0) then ! Process limits for output
        Do i0 = 1, nmo
            write (normaloutput, *) 'eps(', i0, ')=', eps(i0)
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
    deallocate (int2i); Call memminus(KIND(int2i), SIZE(int2i), 1)
    deallocate (inttwi); Call memminus(KIND(inttwi), SIZE(inttwi), 1)
    deallocate (indtwi); Call memminus(KIND(indtwi), SIZE(indtwi), 1)
    deallocate (oner); Call memminus(KIND(oner), SIZE(oner), 1)
    deallocate (int2r); Call memminus(KIND(int2r), SIZE(int2r), 1)
    deallocate (inttwr); Call memminus(KIND(inttwr), SIZE(inttwr), 1)
    deallocate (indtwr); Call memminus(KIND(indtwr), SIZE(indtwr), 1)
    deallocate (int2r_f1); Call memminus(KIND(int2r_f1), SIZE(int2r_f1), 1)
    deallocate (int2i_f1); Call memminus(KIND(int2i_f1), SIZE(int2i_f1), 1)
    deallocate (int2r_f2); Call memminus(KIND(int2r_f2), SIZE(int2r_f2), 1)
    deallocate (int2i_f2); Call memminus(KIND(int2i_f2), SIZE(int2i_f2), 1)
    if (rank == 0) then ! Process limits for output
        write (normaloutput, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024

        Call timing(val(3), totalsec, date0, tsec0)
        write (normaloutput, *) 'End r4dcasci_ty part'
    end if
    call MPI_FINALIZE(ierr)
    if (rank == 0) then ! Process limits for output
        write (normaloutput, '(a,i4,a,i4)') 'fin. rank:', rank, 'nprocs:', nprocs
        close (normaloutput)
    end if
1000 continue
END program r4dcasci_co
