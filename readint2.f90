! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE readint2 (filename, nuniq) ! 2 electorn integrals in MDCINT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        use four_caspt2_module

        Implicit NONE

        character*50,intent(in) :: filename

        character  :: datex*10, timex*8

        integer :: mdcint, nkr, idum, nuniq, nmom
        integer :: nz, type
        integer :: j0, i0, i1
        integer :: k0, l0, ii, jj, kk, ll, signind
        integer :: i, j, k, l, ikr, jkr, lkr, kkr
        integer :: SignIJ, SignKL, itr, jtr, ltr, ktr, inz, totalint

        integer, allocatable :: indk(:), indl(:), kr(:) 

        real*8, allocatable :: rklr(:), rkli(:), int2rs(:), int2is(:)

        logical :: breit


        Allocate(int2rs(0:nmo**4)); Call memplus(KIND(int2rs),SIZE(int2rs),1)
        Allocate(int2is(0:nmo**4)); Call memplus(KIND(int2is),SIZE(int2is),1)

        Allocate(kr(-nmo/2:nmo/2)); Call memplus(KIND(kr),SIZE(kr),1)
        Allocate(indtwr(nmo,nmo,nmo,nmo)); Call memplus(KIND(indtwr),SIZE(indtwr),1)
        Allocate(indtwi(nmo,nmo,nmo,nmo)); Call memplus(KIND(indtwi),SIZE(indtwi),1)


        kr = 0

        Allocate(indk((nmo/2)**2)); Call memplus(KIND(indk),SIZE(indk),1)
        Allocate(indl((nmo/2)**2)); Call memplus(KIND(indl),SIZE(indl),1)
        Allocate(rklr((nmo/2)**2)); Call memplus(KIND(rklr),SIZE(rklr),1)
        Allocate(rkli((nmo/2)**2)); Call memplus(KIND(rkli),SIZE(rkli),1)


        write(*,'("Current Memory is ",F10.2,"MB")')tmem/1024/1024

        nuniq = 0
        indk(:) = 0
        indl(:) = 0
        rklr(:) = 0.0d+00
        rkli(:) = 0.0d+00
        int2r(:) = 0.0d+00
        int2i(:) = 0.0d+00
        indtwr = 0
        indtwi = 0


!###########################################################
!  THIS PART IS TAKEN FROM GOSCIP MOLFDIR PROGRAM PACKAGE
!###########################################################

        totalint = 0
        mdcint=11
        open( mdcint, file=trim(filename),form ='unformatted', status='unknown', err=10)

!old        Read (mdcint,err=20,end=30) datex,timex,nkr, & 
!old        (idum,i0=1,4*nkr),(kr(i0),kr(-1*i0),i0=1,nkr)
        Read (mdcint,err=20,end=30) datex,timex,nkr, & 
        (kr(i0),kr(-1*i0),i0=1,nkr)

        write(*,*) datex,timex
        write(*,*) 'nkr',nkr,'kr(+),kr(-)', (kr(i0),kr(-1*i0),i0=1,nkr)

 60      read (mdcint,ERR=40,END=50) ikr,jkr,nz, &
                  (indk(inz),indl(inz),inz=1,nz), &
                  (rklr(inz),rkli(inz),inz=1,nz)
                  
                  if (ikr==0) goto 50

                  totalint = totalint + nz

                  i = kr(ikr)
                  itr = kr(-ikr)
                  j = kr(jkr)
                  jtr = kr(-jkr)

                  nmom = ninact + nact + nsec

!                  If(i > ninact+nact .and. itr > ninact+nact .and. &
!                  & j > ninact+nact .and. jtr > ninact+nact) goto 60


                  SignIJ = SIGN(1,ikr) * SIGN(1,jkr)

                  Do inz = 1, nz

                     kkr = indk(inz)
                     k = kr(kkr)
                     ktr = kr(-kkr)
                     lkr = indl(inz)
                     l = kr(lkr)
                     ltr = kr(-lkr)

                     If(i > ninact+nact .and. j > ninact+nact .and. &
                     &  k > ninact+nact .and. l > ninact+nact) goto 70

!                     If(i > ninact+nact .and. j > ninact+nact .and. &
!                     &  k > ninact+nact .and. l > ninact+nact) goto 70

                     SignKL = SIGN(1,kkr) * SIGN(1,lkr)
                     nuniq = nuniq + 1


!=-> Original integral plus time-reversed partners
                     INDTWR(I,J,K,L) = NUNIQ
                     INDTWR(JTR,ITR,K,L) = NUNIQ * SignIJ
                     INDTWR(I,J,LTR,KTR) = NUNIQ * SignKL
                     INDTWR(JTR,ITR,LTR,KTR) = NUNIQ * SignIJ * SignKL
                     INDTWI(I,J,K,L) = NUNIQ
                     INDTWI(JTR,ITR,K,L) = NUNIQ * SignIJ
                     INDTWI(I,J,LTR,KTR) = NUNIQ * SignKL
                     INDTWI(JTR,ITR,LTR,KTR) = NUNIQ * SignIJ * SignKL
!=-> Complex conjugate plus time-reversed partners
                     INDTWR(J,I,L,K) = NUNIQ
                     INDTWR(ITR,JTR,L,K) = NUNIQ * SignIJ
                     INDTWR(J,I,KTR,LTR) = NUNIQ * SignKL
                     INDTWR(ITR,JTR,KTR,LTR) = NUNIQ * SignIJ * SignKL
                     INDTWI(J,I,L,K) = - NUNIQ
                     INDTWI(ITR,JTR,L,K) = - NUNIQ * SignIJ
                     INDTWI(J,I,KTR,LTR) = - NUNIQ * SignKL
                     INDTWI(ITR,JTR,KTR,LTR) = - NUNIQ * SignIJ * SignKL
!=-> Particle interchanged plus time-reversed partners
                     INDTWR(K,L,I,J) = NUNIQ
                     INDTWR(LTR,KTR,I,J) = NUNIQ * SignKL
                     INDTWR(K,L,JTR,ITR) = NUNIQ * SignIJ
                     INDTWR(LTR,KTR,JTR,ITR) = NUNIQ * SignIJ * SignKL
                     INDTWI(K,L,I,J) = NUNIQ
                     INDTWI(LTR,KTR,I,J) = NUNIQ * SignKL
                     INDTWI(K,L,JTR,ITR) = NUNIQ * SignIJ
                     INDTWI(LTR,KTR,JTR,ITR) = NUNIQ * SignIJ * SignKL
!=-> Particle interchanged and complex conjugated plus time-reversed partners
                     INDTWR(L,K,J,I) = NUNIQ
                     INDTWR(KTR,LTR,J,I) = NUNIQ * SignKL
                     INDTWR(L,K,ITR,JTR) = NUNIQ * SignIJ
                     INDTWR(KTR,LTR,ITR,JTR) = NUNIQ * SignIJ * SignKL
                     INDTWI(L,K,J,I) = - NUNIQ
                     INDTWI(KTR,LTR,J,I) = - NUNIQ * SignKL
                     INDTWI(L,K,ITR,JTR) = - NUNIQ * SignIJ
                     INDTWI(KTR,LTR,ITR,JTR) = - NUNIQ * SignIJ * SignKL

                     int2rs(nuniq) = rklr(inz)
                     int2is(nuniq) = rkli(inz)

!                     If(abs(rklr(inz))>1.0d-1) write(*,*)rklr(inz),rkli(inz), &
!                     & i, j, k, l 

                     if(abs(rkli(inz)) > thres) realc = .false.

!!                     if((nuniq == 742).or.(nuniq == 2082)) then
!!                     write(*,*)int2r(nuniq)
!!                     write(*,*)int2i(nuniq)
!!      write(*,5)I,J,K,L        ,INDTWR(I,J,K,L)        ,INDTWI(I,J,K,L)  &
!!      &  ,JTR,ITR,K,L    ,INDTWR(JTR,ITR,K,L)    ,INDTWI(JTR,ITR,K,L)    &
!!      &  ,I,J,LTR,KTR    ,INDTWR(I,J,LTR,KTR)    ,INDTWI(I,J,LTR,KTR)    &
!!      &  ,JTR,ITR,LTR,KTR,INDTWR(JTR,ITR,LTR,KTR),INDTWI(JTR,ITR,LTR,KTR)
!!      
!!      write(*,5)J,I,L,K        ,INDTWR(J,I,L,K)        ,INDTWI(J,I,L,K) &       
!!      &  ,ITR,JTR,L,K    ,INDTWR(ITR,JTR,L,K)    ,INDTWI(ITR,JTR,L,K)    &       
!!      &  ,J,I,KTR,LTR    ,INDTWR(J,I,KTR,LTR)    ,INDTWI(J,I,KTR,LTR)    &       
!!      &  ,ITR,JTR,KTR,LTR,INDTWR(ITR,JTR,KTR,LTR),INDTWI(ITR,JTR,KTR,LTR)
!!      
!!      write(*,5)K,L,I,J        ,INDTWR(K,L,I,J)        ,INDTWI(K,L,I,J)  &
!!      &  ,LTR,KTR,I,J    ,INDTWR(LTR,KTR,I,J)    ,INDTWI(LTR,KTR,I,J)    &
!!      &  ,K,L,JTR,ITR    ,INDTWR(K,L,JTR,ITR)    ,INDTWI(K,L,JTR,ITR)    &
!!      &  ,LTR,KTR,JTR,ITR,INDTWR(LTR,KTR,JTR,ITR),INDTWI(LTR,KTR,JTR,ITR)
!!
!!    write(*,5)L,K,J,I        ,INDTWR(L,K,J,I)        ,INDTWI(L,K,J,I) &       
!!      &  ,KTR,LTR,J,I    ,INDTWR(KTR,LTR,J,I)    ,INDTWI(KTR,LTR,J,I)    &       
!!    &  ,L,K,ITR,JTR    ,INDTWR(L,K,ITR,JTR)    ,INDTWI(L,K,ITR,JTR)    &       
!!    &  ,KTR,LTR,ITR,JTR,INDTWR(KTR,LTR,ITR,JTR),INDTWI(KTR,LTR,ITR,JTR)
!! end if

 5    FORMAT(4(4I3,2I6))


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
         write(*,*)nuniq,totalint

         Allocate(int2r(0:nuniq)); Call memplus(KIND(int2r),SIZE(int2r),1)

         int2r(0:nuniq) = int2rs(0:nuniq)

         Deallocate(int2rs); Call memminus(KIND(int2rs),SIZE(int2rs),1)

         Allocate(int2i(0:nuniq)); Call memplus(KIND(int2i),SIZE(int2i),1)

         int2i(0:nuniq) = int2is(0:nuniq)

         Deallocate(int2is); Call memminus(KIND(int2is),SIZE(int2is),1)



         deallocate (indk); Call memminus(KIND(indk),SIZE(indk),1)
         deallocate (indl); Call memminus(KIND(indl),SIZE(indl),1)
         deallocate (rklr); Call memminus(KIND(rklr),SIZE(rklr),1)
         deallocate (rkli); Call memminus(KIND(rkli),SIZE(rkli),1)
         deallocate (kr  ); Call memminus(KIND(kr  ),SIZE(kr  ),1)

         end subroutine readint2

