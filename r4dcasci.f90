! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM r4dcasci   ! DO CASCI CALC IN THIS PROGRAM!

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
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
!
    integer :: ierr, nprocs, rank
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_rank(MPI_COMM_WORLD, rank, ierr)
    write (*, *) ''
    write (*, *) ' ENTER R4DCASCI PROGRAM written by M. Abe'
    write (*, *) ''

    tmem = 0.0d+00

    write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024

    val = 0
    Call DATE_AND_TIME(VALUES=val)
    Write (*, *) 'Year = ', val(1), 'Mon = ', val(2), 'Date = ', val(3)
    Write (*, *) 'Hour = ', val(5), 'Min = ', val(6), 'Sec = ', val(7), '.', val(8)

    totalsec = val(8)*(1.0d-03) + val(7) + val(6)*(6.0d+01) + val(5)*(6.0d+01)**2
    initdate = val(3)
    inittime = totalsec

    write (*, *) inittime

    Call timing(val(3), totalsec, date0, tsec)

    open (5, file='active.inp', form='formatted', status='old')
    read (5, '(I4)') ninact
    read (5, '(I4)') nact
    read (5, '(I4)') nsec
    read (5, '(I4)') nelec
    read (5, '(I4)') nroot
    read (5, '(I4)') selectroot
    read (5, '(I4)') totsym
    read (5, '(I4)') ncore
    read (5, '(I4)') nbas
    close (5)

    write (*, *) 'ninact     =', ninact
    write (*, *) 'nact       =', nact
    write (*, *) 'nsec       =', nsec
    write (*, *) 'nelec      =', nelec
    write (*, *) 'nroot      =', nroot
    write (*, *) 'selectroot =', selectroot
    write (*, *) 'totsym     =', totsym
    write (*, *) 'ncore      =', ncore
    write (*, *) 'nbas       =', nbas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    filename = 'MRCONEE'

    call readorb_enesym(filename)
!       call readorb_enec1 (filename)

    call read1mo(filename)

    write (*, *) 'realc', realc, ECORE, ninact, nact, nsec, nmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Iwamuro create new ikr for dirac
    Call create_newmdcint

    filename = 'MDCINTNEW'

    Call readint2_casci_co(filename, nuniq)

!        Allocate(sp(1:nmo)) ;  Call memplus(KIND(sp),SIZE(sp),1)
!        sp( 1               : ninact           )    = 1
!        sp( ninact+1        : ninact+nact      )    = 2
!        sp( ninact+nact+1   : ninact+nact+nsec )    = 3
!        sp( ninact+nact+nsec: nmo              )    = 4
!        write(*,*)'nmo        =' ,nmo

    nmo = ninact + nact + nsec

!      write(*,*)'Iwamuro debug1'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      write(*,*)'Iwamuro debug2'

    If (mod(nelec, 2) == 0) then
        inisym = nsymrpa + 1
        endsym = 2*nsymrpa
    Else
        inisym = 1
        endsym = nsymrpa
    End if

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

    realcvec = .TRUE.

    Call casci(totsym)

!      goto 1000

!    This is test for bug fix about realc part

    write (*, *) realc, 'realc'
    write (*, *) realcvec, 'realcvec'

    test = .true.

    write (*, *) realc, 'realc'
    write (*, *) realcvec, 'realcvec'

    realc = .FALSE.      !!!      realc =.TRUE.
    realcvec = .FALSE.   !!!      realcvec =.TRUE.

    write (*, *) 'FOR TEST WE DO (F,F)'
    write (*, *) realc, 'realc'
    write (*, *) realcvec, 'realcvec'

!!=============================================!
!                                              !
    iroot = selectroot
!                                              !
!!=============================================!

    Call e0test_v2

    write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!            BUILDING  FOCK MATRIX               !
!  fij = hij + SIGUMA[<0|Ekl|0>{(ij|kl)-(il|kj)} !
!                 kl                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!! TEST TO CALCULATE FOCK MATRIX OF HF STATE fpq = hpq + SIGUMA_r[(pq|rr)-(pr|qr)]
!! THIS MUST BE DIAGONAL MATRIX AND DIAGONAL ELEMENTS CORESPONDS TO SPINOR ENERGIES.

    Allocate (f(nmo, nmo)); Call memplus(KIND(f), SIZE(f), 2)

    f(:, :) = 0.0d+00

!     debug = .FALSE.
    debug = .TRUE.
    If (debug) then
        Call fockhf1
    End if

!Iwamuro modify
!      write(*,'(20I10)') (j0,j0=1,nmo)
!      Do i0 = 1, nmo
!         write(*,'(20F10.3)') (real(f(i0,j0)),j0=1,nmo)
!      End do

!! NOW MAKE FOCK MATRIX FOR CASCI STATE
!! fij = hij + SIGUMA_kl[<0|Ekl|0>{(ij|kl)-(il|kj)}

    f(:, :) = 0.0d+00

!Iwamuro modify
!      write(*,'(20I10)') (j0,j0=1,nmo)
!      Do i0 = 1, nmo
!         write(*,'(20F10.3)') (real(f(i0,j0)),j0=1,nmo)
!      End do

    Call fockcasci

!      debug = .TRUE.
    debug = .FALSE.
    write (*, *) debug, 'debug'

!Iwamuro modify
!      write(*,'(20I10)') (j0,j0=1,nmo)
!      Do i0 = 1, nmo
!         write(*,'(20F10.3)') (real(f(i0,j0)),j0=1,nmo)
!      End do

    if (debug) Call prtoutfock

    Allocate (eps(nmo)); Call memplus(KIND(eps), SIZE(eps), 1)
    eps = 0.0d+00

!Iwamuro modify
!      write(*,'(20I10)') (j0,j0=1,nmo)
!      Do i0 = 1, nmo
!         write(*,'(20F10.3)') (real(f(i0,j0)),j0=1,nmo)
!      End do

    Call fockdiag

    Do i0 = 1, nmo
        write (*, *) 'eps(', i0, ')=', eps(i0)
    End do

!      write(*,'(20I10)') (j0,j0=1,nmo)
!      Do i0 = 1, nmo
!         write(*,'(20F10.3)') (f(i0,j0),j0=1,nmo)
!      enddo

!      Do i0 = 1, nmo/2
!         if(ABS(eps(i0*2)-eps(i0*2-1)) > 1.0d-10) then
!            write(*,*)i0*2-1,i0*2,eps(i0*2-1),eps(i0*2)
!         Endif
!      Enddo

    open (5, file='EPS', form='unformatted', status='unknown')
    write (5) nmo
    write (5) eps(1:nmo)
    close (5)

    deallocate (sp); Call memplus(KIND(sp), SIZE(sp), 1)
    deallocate (cir); Call memminus(KIND(cir), SIZE(cir), 1)
    deallocate (cii); Call memminus(KIND(cii), SIZE(cii), 1)
    deallocate (eigen); Call memminus(KIND(eigen), SIZE(eigen), 1)
    deallocate (f); Call memminus(KIND(f), SIZE(f), 2)
    deallocate (eps); Call memminus(KIND(eps), SIZE(eps), 1)
    deallocate (idet); Call memminus(KIND(idet), SIZE(idet), 1)

    deallocate (orb); Call memminus(KIND(orb), SIZE(orb), 1)
    deallocate (irpmo); Call memminus(KIND(irpmo), SIZE(irpmo), 1)
    deallocate (irpamo); Call memminus(KIND(irpamo), SIZE(irpamo), 1)
    deallocate (indmo); Call memminus(KIND(indmo), SIZE(indmo), 1)
    deallocate (indmor); Call memminus(KIND(indmor), SIZE(indmor), 1)
    deallocate (onei); Call memminus(KIND(onei), SIZE(onei), 1)
    deallocate (int2i); Call memminus(KIND(int2i), SIZE(int2i), 1)
    deallocate (indtwi); Call memminus(KIND(indtwi), SIZE(indtwi), 1)
    deallocate (oner); Call memminus(KIND(oner), SIZE(oner), 1)
    deallocate (int2r); Call memminus(KIND(int2r), SIZE(int2r), 1)
    deallocate (indtwr); Call memminus(KIND(indtwr), SIZE(indtwr), 1)
    deallocate (int2r_f1); Call memminus(KIND(int2r_f1), SIZE(int2r_f1), 1)
    deallocate (int2i_f1); Call memminus(KIND(int2i_f1), SIZE(int2i_f1), 1)
    deallocate (int2r_f2); Call memminus(KIND(int2r_f2), SIZE(int2r_f2), 1)
    deallocate (int2i_f2); Call memminus(KIND(int2i_f2), SIZE(int2i_f2), 1)

    write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024

    Call timing(val(3), totalsec, date0, tsec0)
    write (*, *) 'End r4dcasci part'

1000 continue
END program r4dcasci