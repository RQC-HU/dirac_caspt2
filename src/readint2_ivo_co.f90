! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE readint2_ivo_co (filename, nuniq) ! 2 electorn integrals created by typart in utchem

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        use four_caspt2_module

        Implicit NONE

        character*50,intent(in) :: filename

        character  :: datex*10, timex*8

        integer    :: mdcint, nkr, idum, nuniq, nmom, nmoc
        integer    :: nz, type
        integer    :: j0, i0, i1
        integer    :: k0, l0, ii, jj, kk, ll, signind
        integer    :: i, j, k, l, jtr0, itr0
        integer    :: SignIJ, SignKL, itr, jtr, ltr, ktr, inz, totalint, save, count

        complex*16 :: cint2

        integer, allocatable :: indk(:), indl(:), kr(:)
        real*8, allocatable  :: rklr(:), rkli(:), int2rs(:), int2is(:)

        logical :: breit

        nmoc = ninact + nact
        nmom = ninact + nact + nsec


        Allocate(int2r_f1(ninact+nact+1:ninact+nact+nsec,ninact+nact+1:ninact+nact+nsec,nmoc,nmoc))
        Allocate(int2i_f1(ninact+nact+1:ninact+nact+nsec,ninact+nact+1:ninact+nact+nsec,nmoc,nmoc))
        Allocate(int2r_f2(ninact+nact+1:ninact+nact+nsec,nmoc,nmoc,ninact+nact+1:ninact+nact+nsec))
        Allocate(int2i_f2(ninact+nact+1:ninact+nact+nsec,nmoc,nmoc,ninact+nact+1:ninact+nact+nsec))
        Call memplus(KIND(int2r_f1),SIZE(int2r_f1),1)
        Call memplus(KIND(int2i_f1),SIZE(int2i_f1),1)
        Call memplus(KIND(int2r_f2),SIZE(int2r_f2),1)
        Call memplus(KIND(int2i_f2),SIZE(int2i_f2),1)

        Allocate(indk((nmo/2)**2)); Call memplus(KIND(indk),SIZE(indk),1)
        Allocate(indl((nmo/2)**2)); Call memplus(KIND(indl),SIZE(indl),1)
        Allocate(rklr((nmo/2)**2)); Call memplus(KIND(rklr),SIZE(rklr),1)
        Allocate(rkli((nmo/2)**2)); Call memplus(KIND(rkli),SIZE(rkli),1)

!Iwamuro modify
        Allocate(kr(-nmo/2:nmo/2))           ; Call memplus(KIND(kr)    ,SIZE(kr)    ,1)

        write(*,'("Current Memory is ",F10.2,"MB")')tmem/1024/1024

        nuniq = 0
        indk(:) = 0
        indl(:) = 0
        rklr(:) = 0.0d+00
        rkli(:) = 0.0d+00
        int2r_f1 = 0.0d+00
        int2i_f1 = 0.0d+00
        int2r_f2 = 0.0d+00
        int2i_f2 = 0.0d+00

        totalint = 0
        mdcint=11
        open( mdcint, file=trim(filename),form ='unformatted', status='old', err=10)

        read (mdcint,err=20,end=30) datex,timex,nkr, &
            (kr(i0),kr(-1*i0),i0=1,nkr)

 60      read (mdcint,ERR=40,END=50) i,j,nz, &
                  (indk(inz),indl(inz),inz=1,nz), &
                  (rklr(inz),rkli(inz),inz=1,nz)

         if (i==0) goto 50

         totalint = totalint + nz
         
         itr = i+(-1)**(mod(i,2)+1)
         jtr = j+(-1)**(mod(j,2)+1)

         i0   = i
         itr0 = itr
         j0   = j
         jtr0 = jtr

         Do inz = 1, nz

            i   = i0
            itr = itr0
            j   = j0
            jtr = jtr0

            k   = indk(inz)
            ktr = k+(-1)**(mod(k,2)+1)
            l   = indl(inz)
            ltr = l+(-1)**(mod(l,2)+1)

            If(i > nmoc .and. j  > nmoc  .and.  k  > nmoc .and. l  > nmoc) goto 70 ! (33|33) is ignored
            If(i==j .and. k > l) goto 70

            if(sp(i)==3.and.sp(j)==3 .and. sp(k)< 3.and.sp(l)==sp(k))  then !(33|11) or (33|22) type

               count = 0

 11            if(mod(i, 2) == 0) then
                  itr = i - 1
               else
                  itr = i + 1
               endif

               if(mod(j, 2) == 0) then
                  jtr = j - 1
               else
                  jtr = j + 1
               endif

               if(mod(k, 2) == 0) then
                  ktr = k - 1
               else
                  ktr = k + 1
               endif

               if(mod(l, 2) == 0) then
                  ltr = l - 1
               else
                  ltr = l + 1
               endif

               SignIJ = (-1.0d+00)**mod(i+j,2)
               SignKL = (-1.0d+00)**mod(k+l,2)

               int2r_f1(i,j,k,l) = rklr(inz)
               int2i_f1(i,j,k,l) = rkli(inz)

               int2r_f1(jtr,itr,k,l) = SignIJ*rklr(inz)
               int2i_f1(jtr,itr,k,l) = SignIJ*rkli(inz)

               int2r_f1(i,j,ltr,ktr) = SignKL*rklr(inz)
               int2i_f1(i,j,ltr,ktr) = SignKL*rkli(inz)

               int2r_f1(jtr,itr,ltr,ktr) = SignIJ*SignKL*rklr(inz)
               int2i_f1(jtr,itr,ltr,ktr) = SignIJ*SignKL*rkli(inz)

               count = count + 1
               cint2 = DCMPLX(rklr(inz),rkli(inz))
               if(count ==1) then
                  Call takekr(i,j,k,l,cint2)              ! Consider Kramers pair
                  rklr(inz) = DBLE(cint2)
                  rkli(inz) = DIMAG(cint2)
                  goto 11
               else
                  goto 70
               endif

            elseif(sp(k)==3.and.sp(l)==3 .and. sp(i)< 3.and.sp(i)==sp(j))  then !(11|33) or (22|33) type
!               write(*,'("type 2",4I4,2E20.10)')i,j,k,l,rklr(inz),rkli(inz)

               count = 0

 21            if(mod(i, 2) == 0) then
                  itr = i - 1
               else
                  itr = i + 1
               endif

               if(mod(j, 2) == 0) then
                  jtr = j - 1
               else
                  jtr = j + 1
               endif

               if(mod(k, 2) == 0) then
                  ktr = k - 1
               else
                  ktr = k + 1
               endif

               if(mod(l, 2) == 0) then
                  ltr = l - 1
               else
                  ltr = l + 1
               endif

               SignIJ = (-1.0d+00)**mod(i+j,2)
               SignKL = (-1.0d+00)**mod(k+l,2)

               int2r_f1(k,l,i,j) = rklr(inz)
               int2i_f1(k,l,i,j) = rkli(inz)

               int2r_f1(k,l,jtr,itr) = SignIJ*rklr(inz)
               int2i_f1(k,l,jtr,itr) = SignIJ*rkli(inz)

               int2r_f1(ltr,ktr,i,j) = SignKL*rklr(inz)
               int2i_f1(ltr,ktr,i,j) = SignKL*rkli(inz)

               int2r_f1(ltr,ktr,jtr,itr) = SignIJ*SignKL*rklr(inz)
               int2i_f1(ltr,ktr,jtr,itr) = SignIJ*SignKL*rkli(inz)

               count = count + 1
               cint2 = DCMPLX(rklr(inz),rkli(inz))
               if(count ==1) then
                  Call takekr(i,j,k,l,cint2)              ! Consider Kramers pair
                  rklr(inz) = DBLE(cint2)
                  rkli(inz) = DIMAG(cint2)
                  goto 21
               else
                  goto 70
               endif

            elseif(max(sp(i),sp(j))==3.and.max(sp(k),sp(l))==3.and. &
                 &  min(sp(i),sp(j))==min(sp(k),sp(l))) then                !(31|31) or (32|32) series 

               count = 0

 12            if(mod(i, 2) == 0) then
                  itr = i - 1
               else
                  itr = i + 1
               endif

               if(mod(j, 2) == 0) then
                  jtr = j - 1
               else
                  jtr = j + 1
               endif

               if(mod(k, 2) == 0) then
                  ktr = k - 1
               else
                  ktr = k + 1
               endif

               if(mod(l, 2) == 0) then
                  ltr = l - 1
               else
                  ltr = l + 1
               endif

               SignIJ = (-1.0d+00)**mod(i+j,2)
               SignKL = (-1.0d+00)**mod(k+l,2)

               
               if(i > j .and. k > l) then ! (31|31) or (32|32) ==> (31|13) or (32|23)

                  int2r_f2(i,j,ltr,ktr) = signKL*rklr(inz)
                  int2i_f2(i,j,ltr,ktr) = signKL*rkli(inz)

!                  write(*,*)i,j,ltr,ktr,int2r_f2(i,j,ltr,ktr),int2i_f2(i,j,ltr,ktr)

               elseif(i > j .and. k < l) then ! (31|13) or (32|23) ==> (31|13) or (32|23)

                  int2r_f2(i,j,k,l) = rklr(inz)
                  int2i_f2(i,j,k,l) = rkli(inz)

!                  write(*,*)i,j,k,l,int2r_f2(i,j,k,l),int2i_f2(i,j,k,l)

               elseif(i < j .and. k < l) then ! (13|13) or (23|23) ==> (31|13) or (32|23)

                  int2r_f2(jtr,itr,k,l) = signIJ*rklr(inz)
                  int2i_f2(jtr,itr,k,l) = signIJ*rkli(inz)

!                  write(*,*)jtr,itr,k,l,int2r_f2(jtr,itr,k,l),int2i_f2(jtr,itr,k,l)

               elseif(i < j .and. k > l) then ! (13|31) or (23|32) ==> (31|13) or (32|23)

                  int2r_f2(jtr,itr,ltr,ktr) = signIJ*signKL*rklr(inz)
                  int2i_f2(jtr,itr,ltr,ktr) = signIJ*signKL*rkli(inz)

!                  write(*,*)jtr,itr,ltr,ktr,int2r_f2(jtr,itr,ltr,ktr),int2i_f2(jtr,itr,ltr,ktr)

               endif

               count = count + 1
               cint2 = DCMPLX(rklr(inz),rkli(inz))
               if(count ==1 .or. count ==3) then
                  Call takekr(i,j,k,l,cint2)              ! Consider Kramers pair
                  rklr(inz) = DBLE(cint2)
                  rkli(inz) = DIMAG(cint2)
                  goto 12
               elseif(count == 2 ) then           ! variables exchange (AA|BB) => (BB|AA)
                  save = i
                  i    = k
                  k    = save
                  save = j
                  j    = l
                  l    = save
                  goto 12
               else
                  goto 70
               endif
            endif

 70      Enddo
	
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

!         Allocate(int2r(0:nuniq)); Call memplus(KIND(int2r),SIZE(int2r),1)
!
!         int2r(0:nuniq) = int2rs(0:nuniq)
!
!         Deallocate(int2rs); Call memminus(KIND(int2rs),SIZE(int2rs),1)
!
!         Allocate(int2i(0:nuniq)); Call memplus(KIND(int2i),SIZE(int2i),1)
!
!         int2i(0:nuniq) = int2is(0:nuniq)
!
!         Deallocate(int2is); Call memminus(KIND(int2is),SIZE(int2is),1)



         deallocate (indk); Call memminus(KIND(indk),SIZE(indk),1)
         deallocate (indl); Call memminus(KIND(indl),SIZE(indl),1)
         deallocate (rklr); Call memminus(KIND(rklr),SIZE(rklr),1)
         deallocate (rkli); Call memminus(KIND(rkli),SIZE(rkli),1)
         deallocate (kr); Call memminus(KIND(kr),SIZE(kr),1)
         end subroutine readint2_ivo_co

