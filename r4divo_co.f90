! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   PROGRAM r4divo_co   ! DO IVO CALC ONLY FOR SMALL BASIS SETS

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
        integer                 :: ii, jj, kk, ll, typetype, i0, j0, nhomo
        integer                 ::  j, i, k, l, nuniq
        integer                 :: k0, l0, nint, n, dimn, n0, n1, nspace(3,3)
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

        write(*,*)''
        write(*,*)' ENTER R4DIVO PROGRAM written by M. Abe test17_ty version 2007/7/20'
        write(*,*)''

        tmem = 0.0d+00

        write(*,'("Current Memory is ",F10.2,"MB")')tmem/1024/1024

        val = 0
        Call DATE_AND_TIME (VALUES=val)
        Write(*,*)'Year = ',val(1),'Mon = ',val(2),'Date = ',val(3)
        Write(*,*)'Hour = ',val(5),'Min = ',val(6),'Sec = ',val(7),'.',val(8)

        totalsec = val(8)*(1.0d-03)+val(7)+val(6)*(6.0d+01)+val(5)*(6.0d+01)**2
        initdate = val(3)
        inittime = totalsec

        write(*,*)inittime

        Call timing(val(3), totalsec, date0, tsec)


        open(5,file='active.inp',form='formatted',status='old')
        read(5,'(I4)')ninact
        read(5,'(I4)')nact
        read(5,'(I4)')nsec
        read(5,'(I4)')nelec
        read(5,'(I4)')nroot
        read(5,'(I4)')selectroot
        read(5,'(I4)')totsym
        read(5,'(I4)')ncore
        read(5,'(I4)')nbas
        read(5,'(E8.2)')eshift
        read(5,'(A6)')ptgrp
        read(5,'(I4)')nhomo
        close(5)

        write(*,*)'ninact     =' ,ninact
        write(*,*)'nact       =' ,nact
        write(*,*)'nsec       =' ,nsec
        write(*,*)'nelec      =' ,nelec
        write(*,*)'nroot      =' ,nroot
        write(*,*)'selectroot =' ,selectroot
        write(*,*)'totsym     =' ,totsym
        write(*,*)'ncore      =' ,ncore
        write(*,*)'nbas       =' ,nbas
        write(*,*)'eshift     =' ,eshift          ! NO USE IN IVO BUT FOR CASCI AND CASPT2 IT IS USED
        write(*,*)'ptgrp      =' ,ptgrp
        write(*,*)'nhomo      =' ,nhomo             

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       filename = 'MRCONEE'
        
       call readorb_enesym_co (filename)
       call read1mo_co (filename)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      filename = 'MDCINT'

!        nmo = ninact + nact + nsec

      Call readint2_ivo_ty (filename, nuniq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(*,'("Current Memory is ",F10.2,"MB")')tmem/1024/1024

      write(*,*)' '
      write(*,*)'*******************************'
      write(*,*)' '
      write(*,*)'IREP IS ',repna(totsym)
      write(*,*)' '
      write(*,*)'*******************************'
      write(*,*)' '

      realcvec = .TRUE.


!      goto 1000
      


!    This is test for bug fix about realc part

      write(*,*)realc,'realc'
      write(*,*)realcvec,'realcvec'

      test = .true.

      write(*,*)realc,'realc'
      write(*,*)realcvec,'realcvec'

      realc =.FALSE.      !!!      realc =.TRUE.   
      realcvec =.FALSE.   !!!      realcvec =.TRUE.

      write(*,*)'FOR TEST WE DO (F,F)'
      write(*,*)realc,'realc'
      write(*,*)realcvec,'realcvec'

!!=============================================!
!                                              !  
      iroot = selectroot
!                                              !
!!=============================================!

!    Call e0test_v2

!      write(*,'("Current Memory is ",F10.2,"MB")')tmem/1024/1024


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!            BUILDING  FOCK MATRIX               !
!  fij = hij + SIGUMA[<0|Ekl|0>{(ij|kl)-(il|kj)} !
!                 kl                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!! TEST TO CALCULATE FOCK MATRIX OF HF STATE fpq = hpq + SIGUMA_r[(pq|rr)-(pr|qr)]
!! THIS MUST BE DIAGONAL MATRIX AND DIAGONAL ELEMENTS CORESPONDS TO SPINOR ENERGIES.

      Allocate(f(nsec,nsec)); Call memplus (KIND(f),SIZE(f),2) 

      f(:,:) = 0.0d+00

!! NOW MAKE FOCK MATRIX FOR IVO (only virtual spinors

!! fij = hij + SIGUMA_a(ij|aa)-(ia|aj)}

      f(:,:) = 0.0d+00
      
      Call fockivo_ty(nhomo)
!      Call fockivo(nhomo)

      deallocate ( f   )   ;   Call memminus (KIND( f    ),SIZE( f   ),2)

      deallocate ( orb    );   Call memminus (KIND( orb    ),SIZE( orb    ),1)
      deallocate ( irpmo  );   Call memminus (KIND( irpmo  ),SIZE( irpmo  ),1)
      deallocate ( irpamo );   Call memminus (KIND( irpamo ),SIZE( irpamo ),1)
      deallocate ( indmo  );   Call memminus (KIND( indmo  ),SIZE( indmo  ),1)
      deallocate (indmor  );   Call memminus (KIND(indmor  ),SIZE(indmor  ),1)
      deallocate (onei    );   Call memminus (KIND(onei    ),SIZE(onei    ),1)
!      deallocate (int2i   );   Call memminus (KIND(int2i   ),SIZE(int2i   ),1)
!      deallocate (indtwi  );   Call memminus (KIND(indtwi  ),SIZE(indtwi  ),1)
      deallocate ( oner   );   Call memminus (KIND( oner   ),SIZE( oner   ),1)
!      deallocate (int2r   );   Call memminus (KIND(int2r   ),SIZE(int2r   ),1)
!      deallocate (indtwr  );   Call memminus (KIND(indtwr  ),SIZE(indtwr  ),1)
      deallocate (int2r_f1);   Call memminus (KIND(int2r_f1),SIZE(int2r_f1),1)
      deallocate (int2i_f1);   Call memminus (KIND(int2i_f1),SIZE(int2i_f1),1)
      deallocate (int2r_f2);   Call memminus (KIND(int2r_f2),SIZE(int2r_f2),1)
      deallocate (int2i_f2);   Call memminus (KIND(int2i_f2),SIZE(int2i_f2),1)
      deallocate (MULTB_S) ;   Call memminus (KIND(MULTB_S),SIZE(MULTB_S),1)
      deallocate (MULTB_D) ;   Call memminus (KIND(MULTB_D),SIZE(MULTB_D),1)
      deallocate (MULTB_DS);   Call memminus (KIND(MULTB_DS),SIZE(MULTB_DS),1)
      deallocate (MULTB_DF);   Call memminus (KIND(MULTB_DF),SIZE(MULTB_DF),1)
      deallocate (MULTB_DB);   Call memminus (KIND(MULTB_DB),SIZE(MULTB_DB),1)
      deallocate (MULTB_SB);   Call memminus (KIND(MULTB_SB),SIZE(MULTB_SB),1)

      write(*,'("Current Memory is ",F10.2,"MB")')tmem/1024/1024


      Call timing(val(3), totalsec, date0, tsec0)
      write(*,*)'End r4divo_ty part'

 1000 continue 
      END program r4divo_co


