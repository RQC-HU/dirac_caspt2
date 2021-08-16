! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE readint2_ord_co (filename) ! 2 electorn integrals created by typart in utchem

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        use four_caspt2_module

        Implicit NONE

        character*50,intent(in) :: filename

        character  :: datex*10, timex*8

        integer :: mdcint, nkr, idum, nmom, max1, max2, min1, min2
        integer :: nz, type
        integer :: j0, i0, i1
        integer :: k0, l0, ii, jj, kk, ll, signind
        integer :: i, j, k, l
        integer :: SignIJ, SignKL, itr, jtr, ltr, ktr, inz, totalint

        integer, allocatable :: indk(:), indl(:), kr(:)

        real*8, allocatable :: rklr(:), rkli(:)

        logical :: breit

!Iwamuro modify
!        integer :: ikr, jkr, kkr, lkr

        Allocate(kr(-nmo/2:nmo/2)); Call memplus(KIND(kr),SIZE(kr),1)
        Allocate(indk((nmo/2)**2)); Call memplus(KIND(indk),SIZE(indk),1)
        Allocate(indl((nmo/2)**2)); Call memplus(KIND(indl),SIZE(indl),1)
        Allocate(rklr((nmo/2)**2)); Call memplus(KIND(rklr),SIZE(rklr),1)
        Allocate(rkli((nmo/2)**2)); Call memplus(KIND(rkli),SIZE(rkli),1)

        write(*,'("Current Memory is ",F10.2,"MB")')tmem/1024/1024

        indk(:) = 0
        indl(:) = 0
        rklr(:) = 0.0d+00
        rkli(:) = 0.0d+00
        kr(:)   = 0

        totalint = 0

        open( 11, file='A1int',form ='unformatted', status='unknown')
        open( 12, file='A2int',form ='unformatted', status='unknown')
        open( 2 , file='Bint' ,form ='unformatted', status='unknown')
        open( 31, file='C1int',form ='unformatted', status='unknown')
        open( 32, file='C2int',form ='unformatted', status='unknown')
        open( 33, file='C3int',form ='unformatted', status='unknown')
        open( 4 , file='D1int',form ='unformatted', status='unknown')
        open( 41, file='D2int',form ='unformatted', status='unknown')
        open( 42, file='D3int',form ='unformatted', status='unknown')
        open( 5 , file='Eint' ,form ='unformatted', status='unknown')
        open( 9 , file='Fint' ,form ='unformatted', status='unknown')
        open( 7 , file='Gint' ,form ='unformatted', status='unknown')
        open( 8 , file='Hint' ,form ='unformatted', status='unknown')


        mdcint=15

        open( mdcint, file=trim(filename),form ='unformatted', status='old', err=10)

        Read (mdcint,err=20,end=30) datex,timex,nkr, &
        (kr(i0),kr(-1*i0),i0=1,nkr)

        write(*,*) datex,timex
        write(*,*) 'nkr',nkr,'kr(+),kr(-)', (kr(i0),kr(-1*i0),i0=1,nkr)


 60           read (mdcint,ERR=40,END=50) i,j,nz, &
                  (indk(inz),indl(inz),inz=1,nz), &
                  (rklr(inz),rkli(inz),inz=1,nz)
!                  write(*,'(3I4)')i,j,nz
                  if (i==0 .and. j==0.and. nz==0) goto 50

                  totalint = totalint + nz

                  itr = i+(-1)**(mod(i,2)+1)
                  jtr = j+(-1)**(mod(j,2)+1)

                  nmom = ninact + nact + nsec

                  If(sp(i)==4 .or. sp(j) == 4) goto 60
                  If(i > ninact+nact .and. j > ninact+nact) goto 60

                  SignIJ = (-1)**(mod(i+j,2))

                  Do inz = 1, nz

                     k   = indk(inz)
                     ktr = k+(-1)**(mod(k,2)+1)
                     l   = indl(inz)
                     ltr = l+(-1)**(mod(l,2)+1)

!                     write(*,'("all ints",4I4,E20.10)')i,j,k,l,rklr(inz)

                     If(sp(k)==4 .or. sp(l) == 4) goto 70
                     If(k > ninact+nact .and. l > ninact+nact) goto 70
                     If(i==j .and. k > l) goto 70

                     SignKL = (-1)**(mod(k+l,2))

                     max1 = max(sp(i), sp(j))
                     min1 = min(sp(i), sp(j))
                     max2 = max(sp(k), sp(l))
                     min2 = min(sp(k), sp(l))

!===============================================================
! Integrals for A space  (pi|qr)(21|22) (pi|jk)(21|11)  type
!===============================================================

                     If(max1==2 .and. min1==2 .and. max2==2 .and. min2==1) then    ! (22|21) => (21|22)
!                        write(*,'(4I4,2E20.10)')i,j,k,l,  rklr(inz),         rkli(inz)

                        if(k > l) then ! (22|21) => (21|22)

                           write(11)k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)
!                           write(*,'("A1int1",4I4,2E20.10)')k  ,l  ,i  ,j  ,  rklr(inz),         rkli(inz)

                        else    ! (22|12) => (22|21)* => (21|22)*

                           write(11)l  ,k  ,j  ,i  ,         rklr(inz), -1.0d+00*rkli(inz)
!                           write(*,'("A1int2",4I4,2E20.10)')l  ,k  ,j  ,i  ,  rklr(inz), -1.0d+00*rkli(inz)
                        endif

                     elseif(max1==2 .and. min1==1 .and. max2==2 .and. min2==2) then ! (21|22) => (21|22)

!                        write(*,'(4I4,2E20.10)')i,j,k,l,  rklr(inz),         rkli(inz)

                        if(i > j) then ! (21|22) => (21|22)

                           write(11)i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)
!                           write(*,'("A1int3",4I4,2E20.10)')i  ,j  ,k  ,l  ,  rklr(inz),         rkli(inz)

                        else    ! (12|22) => (21|22)*

                           write(11)j  ,i  ,l  ,k  ,         rklr(inz),-1.0d+00*rkli(inz)
!                           write(*,'("A1int4",4I4,2E20.10)')j  ,i  ,l  ,k  ,  rklr(inz),-1.0d+00*rkli(inz)

                        endif

                     elseif(max1==2 .and. min1==1 .and. max2==1 .and. min2==1) then  ! (21|11)=>(21|11)
!                        write(*,'(4I4,2E20.10)')i,j,k,l,  rklr(inz),         rkli(inz)

                        if(i > j) then ! (21|11) => (21|11)

                           write(12)i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)
!                           write(*,'("A2int1",4I4,2E20.10)')i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                        else    ! (12|11) => (21|11)* => (21|11)*

                           write(12)j  ,i  ,l  ,k  ,         rklr(inz), -1.0d+00*rkli(inz)
!                           write(*,'("A2int2",4I4,2E20.10)')j  ,i  ,l  ,k  ,         rklr(inz), -1.0d+00*rkli(inz)
                        endif

                     elseif(max1==1 .and. min1==1 .and. max2==2 .and. min2==1) then  ! (11|21)=>(21|11)
!                        write(*,'(4I4,2E20.10)')i,j,k,l,  rklr(inz),         rkli(inz)

                        if(k > l) then   ! (11|21) => (21|11)

                           write(12)k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)
!                           write(*,'("A2int3",4I4,2E20.10)')k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)

                        else    ! (11|12) => (11|21)* => (21|11)*

                           write(12)l  ,k  ,j  ,i  ,         rklr(inz), -1.0d+00*rkli(inz)
!                           write(*,'("A2int4",4I4,2E20.10)')l  ,k  ,j  ,i  ,         rklr(inz), -1.0d+00*rkli(inz)

                        endif


!=============================================
! Integrals for B space  (pi|qj) (21|21) type
!=============================================


                     elseif(max1==2 .and. min1==1 .and. max2==2 .and. min2==1) then  ! (21|21)=>(21|21)

                        if(i > j .and. k > l) then   ! (21|21) => (21|21)

                           write(2)i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                        elseif(i < j .and. k > l) then ! (12|21) => (21|21)

                           write(2)jtr,itr,k  ,l  ,  SignIJ*rklr(inz),  SignIJ*rkli(inz)

                        elseif(i > j .and. k < l) then ! (21|12) => (21|21)

                           write(2)i  ,j  ,ltr,ktr,  SignKL*rklr(inz),  SignKL*rkli(inz)

                        elseif(i < j .and. k < l) then ! (12|12) => (21|21)*

                           write(2)jtr,itr,ltr,ktr,SignIJ*SignKL*rklr(inz),SignIJ*SignKL*rkli(inz)

                        endif


!============================================================================
! Integrals for C space (ap|qr)(32|22) type C1int
!============================================================================


                     elseif(max1==3 .and. min1==2 .and. max2==2 .and. min2==2) then ! (32|22)=>(32|22)

                        if(i > j) then ! (32|22)=>(32|22)

                           write(31)i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)
!Iwamuro modify
!                           write(*,'("C1int1",4I4,2E20.10)')i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                        else    ! (23|22)=>(32|22)

                           write(31)jtr,itr,k  ,l  ,  SignIJ*rklr(inz),  SignIJ*rkli(inz)
!Iwamuro modify
!                           write(*,'("C1int2",4I4,2E20.10)')jtr,itr,k  ,l  ,  SignIJ*rklr(inz),  SignIJ*rkli(inz)
                        endif

                     elseif(max1==2 .and. min1==2 .and. max2==3 .and. min2==2) then ! (22|32)=>(32|22)

                        if(k > l) then ! (22|32)=>(32|22)

                           write(31)k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)
!Iwamuro modify
!                           write(*,'("C1int3",4I4,2E20.10)')k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)
                        else    ! (22|23)=>(32|22)

                           write(31)ltr,ktr,i  ,j  ,  SignKL*rklr(inz),  SignKL*rkli(inz)
!Iwamuro modify
!                          write(*,'("C1int4",4I4,2E20.10)')ltr,ktr,i  ,j  ,  SignKL*rklr(inz),  SignKL*rkli(inz)
                        endif

!============================================================================
! Integrals for C space (ap|kk)(32|11)  type C2int
!============================================================================


                     elseif(max1==3 .and. min1==2 .and. max2==1 .and. min2==1)then   ! (32|11)=>(32|11)

                        if(i > j) then ! (32|11)=>(32|11)

                           write(32)i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                        else    ! (23|11)=>(32|11)

                           write(32)jtr,itr,k  ,l  ,  SignIJ*rklr(inz),  SignIJ*rkli(inz)

                        endif

                     elseif(max1==1 .and. min1==1 .and. max2==3 .and. min2==2)then   ! (32|11)=>(32|11)

                        if(k > l) then ! (11|32)=>(32|11)

                           write(32)k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)

                        else    ! (11|23)=>(32|11)

                           write(32)ltr,ktr,i  ,j  ,  SignKL*rklr(inz),  SignKL*rkli(inz)

                        endif

!============================================================================
! Integrals for C (ai|jp) (31|12)(C3int) and E space (ai|pj)(31|21) (Eint)
!============================================================================


                     elseif(max1==3 .and. min1==1 .and. max2==2 .and. min2==1) then ! (31|21)=>(31|12)

                        if    (i > j .and. l > k) then ! (31|12)=>(31|21) For E
                           write(5 )i  ,j  ,ltr,ktr,   SignKL*rklr(inz),  SignKL*rkli(inz)
                           write(33)i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                        elseif(j > i .and. l > k ) then ! (13|12)=>(31|21) For E
                           write(5 )jtr,itr,ltr,ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                           write(33)jtr,itr,k  ,l  ,  SignIJ*rklr(inz), SignIJ*rkli(inz)

                        elseif(i > j .and. k > l ) then ! (31|21)=>(31|21) For E
                           write(5 )i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)
                           write(33)i  ,j  ,ltr,ktr,   SignKL*rklr(inz),  SignKL*rkli(inz)

                        elseif(i < j .and. k > l ) then ! (13|21)=>(31|21) For E
                           write(5 )jtr,itr,k  ,l  ,  SignIJ*rklr(inz), SignIJ*rkli(inz)
                           write(33)jtr,itr,ltr,ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)

                        endif


                     elseif(max1==2 .and. min1==1 .and. max2==3 .and. min2==1) then ! (21|31)=>(31|12)

                        if    (i > j .and. l > k ) then ! (21|13)=>(31|21) For E
                           write(5 )ltr,ktr,i  ,j  ,        SignKL*rklr(inz),        SignKL*rkli(inz)
                           write(33)ltr,ktr,jtr,itr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)

                        elseif(j > i .and. l > k ) then ! (12|13)=>(31|21) For E
                           write(5 )ltr,ktr,jtr,itr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                           write(33)ltr,ktr,i  ,j  ,        SignKL*rklr(inz),        SignKL*rkli(inz)

                        elseif(i > j .and. k > l) then ! (21|31)=>(31|21) For E
                           write(5 )k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)
                           write(33)k  ,l  ,jtr,itr,  SignIJ*rklr(inz),  SignIJ*rkli(inz)

                        elseif(i < j .and. k > l) then ! (12|31)=>(31|21) For E
                           write(5 )k  ,l  ,jtr,itr,  SignIJ*rklr(inz),  SignIJ*rkli(inz)
                           write(33)k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)

                        endif


!============================================================================
! Integrals for D space (ai|pq)(31|22) type (D1int)
!============================================================================


                     elseif(max1==3 .and. min1==1 .and. max2==2 .and. min2==2) then ! (31|22)=>(31|22)

                        if(i > j) then ! (31|22)=>(31|22)

                           write(4)i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)


                        else    ! (13|22)=>(31|22)

                           write(4)jtr,itr,k  ,l  ,  SignIJ*rklr(inz),  SignIJ*rkli(inz)

                        endif

                     elseif(max1==2 .and. min1==2 .and. max2==3 .and. min2==1) then ! (22|31)=>(31|22)

                        if(k > l) then ! (22|31)=>(31|22)

                           write(4)k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)

                        else    ! (22|13)=>(31|22)

                           write(4)ltr,ktr,i  ,j  ,  SignKL*rklr(inz),  SignKL*rkli(inz)

                        endif


!============================================================================
! Integrals for D space (ap|qi)(32|21) type (D2int)
!============================================================================


                     elseif(max1==3 .and. min1==2 .and. max2==2 .and. min2==1) then ! (32|21)=>(32|21)

                        if(i > j .and. k > l) then ! (32|21)=>(32|21)

                           write(41)i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                        elseif(i < j .and. k > l) then ! (23|21)=>(32|21)

                           write(41)jtr,itr,k  ,l  ,  SignIJ*rklr(inz),  SignIJ*rkli(inz)

                        elseif(i > j .and. k < l) then ! (32|12)=>(32|21)

                           write(41)i  ,j  ,ltr,ktr,  SignKL*rklr(inz),  SignKL*rkli(inz)

                        elseif(i < j .and. k < l) then ! (23|12)=>(32|21)

                           write(41)jtr,itr,ltr,ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)

                        endif

                     elseif(max1==2 .and. min1==1 .and. max2==3 .and. min2==2) then ! (21|32)=>(32|21)

                        if(i > j .and. k > l) then ! (21|32)=>(32|21)

                           write(41)k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)

                        elseif(i < j .and. k > l) then ! (12|32)=>(32|21)

                           write(41)k  ,l  ,jtr,itr,  SignIJ*rklr(inz),  SignIJ*rkli(inz)

                        elseif(i > j .and. k < l) then ! (21|23)=>(32|21)

                           write(41)ltr,ktr,i  ,j  ,  SignKL*rklr(inz),  SignKL*rkli(inz)

                        elseif(i < j .and. k < l) then ! (12|23)=>(32|21)

                           write(41)ltr,ktr,jtr,itr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)

                        endif


!============================================================================
! Integrals for D space (ai|jk)  (31|11) type (D3int)
!============================================================================

                     elseif(max1==3 .and. min1==1 .and. max2==1 .and. min2==1) then   ! (31|11)=>(31|11)

                        if(i > j) then ! (ai|jk) (31|11)=>(31|11)

                           write(42) i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                        else    ! (i~a~|kk) (13|11)=>(31|11)

                           write(42)jtr,itr,k  ,l  ,  SignIJ*rklr(inz),  SignIJ*rkli(inz)

                        endif

                     elseif(max1==1 .and. min1==1 .and. max2==3 .and. min2==1) then  ! (11|31)=>(31|11)

                        if(k > l) then ! (jk|ai) (31|11)=>(31|11)

                           write(42) k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)

                        else  ! (jk|i~a~)=>( ai|kk) (11|13)=>(31|11)

                           write(42) ltr,ktr,i  ,j  , SignKL*rklr(inz),  SignKL*rkli(inz)

                        endif


!=============================================
! Integrals for F space  (ap|bq) (32|32) type
!=============================================


                     elseif(max1==3 .and. min1==2 .and. max2==3 .and. min2==2) then  ! (32|32)=>(32|32)

                        if(i > j .and. k > l) then   ! (32|32) => (32|32)

                           write(9)i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                        elseif(i < j .and. k > l) then ! (23|32) => (32|32)

                           write(9)jtr,itr,k  ,l  ,  SignIJ*rklr(inz),  SignIJ*rkli(inz)

                        elseif(i > j .and. k < l) then ! (32|23) => (32|32)

                           write(9)i  ,j  ,ltr,ktr,  SignKL*rklr(inz),  SignKL*rkli(inz)

                        elseif(i < j .and. k < l) then ! (23|23) => (32|32)

                           write(9)jtr,itr,ltr,ktr,SignIJ*SignKL*rklr(inz),SignIJ*SignKL*rkli(inz)

                        endif


!============================================================================
! G space (ai|bp)(31|32) type
!============================================================================


                     elseif(max1==3 .and. min1==1 .and. max2==3 .and. min2==2) then ! (31|32)=>(31|32)

                        if    (i > j .and. l > k) then ! (31|23)=>(31|32)
                           write(7)i  ,j  ,ltr,ktr,   SignKL*rklr(inz),  SignKL*rkli(inz)
!                           write(*,'("Gint1",4I4,2E20.10)')i  ,j  ,ltr,ktr,   SignKL*rklr(inz),  SignKL*rkli(inz)

                        elseif(j > i .and. l > k ) then ! (13|23)=>(31|32)
                           write(7)jtr,itr,ltr,ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
!                           write(*,'("Gint2",4I4,2E20.10)')jtr,itr,ltr,ktr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)

                        elseif(i > j .and. k > l ) then ! (31|32)=>(31|32)
                           write(7)i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)
!                           write(*,'("Gint3",4I4,2E20.10)')i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                        elseif(i < j .and. k > l ) then ! (13|32)=>(31|32)
                           write(7)jtr,itr,k  ,l  ,  SignIJ*rklr(inz), SignIJ*rkli(inz)
!                           write(*,'("Gint4",4I4,2E20.10)')jtr,itr,k  ,l  ,  SignIJ*rklr(inz), SignIJ*rkli(inz)

                        endif



                     elseif(max1==3 .and. min1==2 .and. max2==3 .and. min2==1) then ! (32|31)=>(31|32)

                        if    (i > j .and. l > k ) then ! (32|13)=>(31|32)
                           write(7)ltr,ktr,i  ,j  ,        SignKL*rklr(inz),        SignKL*rkli(inz)
!                           write(*,'("Gint5",4I4,2E20.10)')ltr,ktr,i  ,j  ,        SignKL*rklr(inz),        SignKL*rkli(inz)

                        elseif(j > i .and. l > k ) then ! (23|13)=>(31|32)
                           write(7)ltr,ktr,jtr,itr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)
                           write(*,'("Gint6",4I4,2E20.10)')ltr,ktr,jtr,itr, SignIJ*SignKL*rklr(inz), SignIJ*SignKL*rkli(inz)

                        elseif(i > j .and. k > l) then ! (32|31)=>(31|32)
                           write(7)k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)
!                           write(*,'("Gint7",4I4,2E20.10)')k  ,l  ,i  ,j  ,         rklr(inz),         rkli(inz)

                        elseif(i < j .and. k > l) then ! (23|31)=>(31|32)
                           write(7)k  ,l  ,jtr,itr,  SignIJ*rklr(inz),  SignIJ*rkli(inz)
!                           write(*,'("Gint8",4I4,2E20.10)')k  ,l  ,jtr,itr,  SignIJ*rklr(inz),  SignIJ*rkli(inz)

                        endif

!=============================================
! Integrals for H space  (ai|bj) (31|31) type
!=============================================


                     elseif(max1==3 .and. min1==1 .and. max2==3 .and. min2==1) then  ! (31|31)=>(31|31)

                        if(i > j .and. k > l) then   ! (31|31) => (31|31)

                           write(8)i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)
!                           write(*,'("Hint1",4I4,2E20.10)')i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)
!                           write(*,*)i  ,j  ,k  ,l  ,         rklr(inz),         rkli(inz)

                        elseif(i < j .and. k > l) then ! (13|31) => (31|31)

                           write(8)jtr,itr,k  ,l  ,  SignIJ*rklr(inz),  SignIJ*rkli(inz)
!                           write(*,'("Hint2",4I4,2E20.10)')jtr,itr,k  ,l  ,  SignIJ*rklr(inz),  SignIJ*rkli(inz)
!                           write(*,*)jtr,itr,k  ,l  ,  SignIJ*rklr(inz),  SignIJ*rkli(inz)

                        elseif(i > j .and. k < l) then ! (31|13) => (31|31)

                           write(8)i  ,j  ,ltr,ktr,  SignKL*rklr(inz),  SignKL*rkli(inz)
!                           write(*,'("Hint3",4I4,2E20.10)')i  ,j  ,ltr,ktr,  SignKL*rklr(inz),  SignKL*rkli(inz)
!                           write(*,*)i  ,j  ,ltr,ktr,  SignKL*rklr(inz),  SignKL*rkli(inz)

                        elseif(i < j .and. k < l) then ! (13|13) => (31|31)

                           write(8)jtr,itr,ltr,ktr,SignIJ*SignKL*rklr(inz),SignIJ*SignKL*rkli(inz)
!                           write(*,'("Hint4",4I4,2E20.10)')jtr,itr,ltr,ktr,SignIJ*SignKL*rklr(inz),SignIJ*SignKL*rkli(inz)
!                           write(*,*)jtr,itr,ltr,ktr,SignIJ*SignKL*rklr(inz),SignIJ*SignKL*rkli(inz)

                        endif

                     endif

!                     if(abs(rkli(inz)) > thres) realc = .false.

 70               Enddo

                  indk(:)=0
                  indl(:)=0
                  rklr = 0.0d+00
                  rkli = 0.0d+00

                  Goto 60




 10      write(*,*)'error for opening mdcint 10'
         go to 100
 20      write(*,*)'error for reading mdcint 20'
         go to 100
 30      write(*,*)'end mdcint 30'
         go to 100
 40      write(*,*)'error for reading mdcint 40'
         go to 100
 50      write(*,*)'end mdcint 50 normal'
         go to 100

 100     continue

         close (mdcint)

!         write(11) 0, 0, 0, 0, 0.0d+00, 0.0d+00
!         write(12) 0, 0, 0, 0, 0.0d+00, 0.0d+00
!         write(2 ) 0, 0, 0, 0, 0.0d+00, 0.0d+00
!         write(31) 0, 0, 0, 0, 0.0d+00, 0.0d+00
!         write(32) 0, 0, 0, 0, 0.0d+00, 0.0d+00
!         write(33) 0, 0, 0, 0, 0.0d+00, 0.0d+00
!         write(4 ) 0, 0, 0, 0, 0.0d+00, 0.0d+00
!         write(41) 0, 0, 0, 0, 0.0d+00, 0.0d+00
!         write(42) 0, 0, 0, 0, 0.0d+00, 0.0d+00
!         write(5 ) 0, 0, 0, 0, 0.0d+00, 0.0d+00
!         write(9 ) 0, 0, 0, 0, 0.0d+00, 0.0d+00
!         write(7 ) 0, 0, 0, 0, 0.0d+00, 0.0d+00
!         write(8 ) 0, 0, 0, 0, 0.0d+00, 0.0d+00

         close (11)
         close (12)
         close (2 )
         close (31)
         close (32)
         close (33)
         close (4 )
         close (41)
         close (42)
         close (5 )
         close (9 )
         close (7 )
         close (8 )

         deallocate (indk); Call memminus(KIND(indk),SIZE(indk),1)
         deallocate (indl); Call memminus(KIND(indl),SIZE(indl),1)
         deallocate (rklr); Call memminus(KIND(rklr),SIZE(rklr),1)
         deallocate (rkli); Call memminus(KIND(rkli),SIZE(rkli),1)
         deallocate (kr  ); Call memminus(KIND(kr  ),SIZE(kr  ),1)
         end subroutine readint2_ord_co
