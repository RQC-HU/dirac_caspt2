! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE solvG_ord (e0, e2g)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE

        real*8,     intent(in) :: e0
        real*8,     intent(out):: e2g



        integer :: dimn, dimm, count, dammy

        integer, allocatable :: indsym(:,:)

        real*8, allocatable  :: sr(:,:), ur(:,:)
        real*8, allocatable  :: br(:,:), wsnew(:), ws(:), wb(:)
        real*8, allocatable  :: br0(:,:), br1(:,:)
        real*8               :: e2(2*nsymrpa), alpha, e


        complex*16, allocatable  :: sc(:,:), uc(:,:), sc0(:,:)
        complex*16, allocatable  :: bc(:,:)
        complex*16, allocatable  :: bc0(:,:), bc1(:,:), v(:,:), vc(:), vc1(:)

        logical                  :: cutoff
        integer                  :: j, i, k, i0, syma, isym, indt(1:nact)
        integer                  :: ia, it, ib, ii, ja, jt, jb, ji
        integer, allocatable     :: ia0(:), ib0(:), ii0(:), iabi(:,:,:)
        integer                  :: nabi


        real*8  :: thresd

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE G IS NOW CALCULATED
!      
!     EaiEbt|0>
!
!   DRAS1 =-1   DRAS2 = -1   DRAS3 = +2 
!
!   c > d, a > b, and impose that c >= a (or a >= c)
!
!  S(cjdu,aibt) = d(ac) d(bd) d(ij) <0|Eut|0>  <= S(u,t)
!                                   ~~~~~~~~~
!  S(u,t) = <0|Eut|0>
!
!  B(u,t) = Siguma_w [eps(w)<0|EutEww|0>] + S(u,t)(-eps(t))
!
!  alpha(i,a,b) = -eps(i) + eps(a) + eps(b) - e0
!
!  where
!
!  e0 = Siguma_w [eps(w)<0|Eww|0>] (<== calculated as e0 in calce0.f)
!
!  V(t,iab)   =  [SIGUMA_p:active <0|Etp|0>{(ai|bp)-(ap|bi)}]
!   
!
!  E2 = SIGUMA_iab, dimm |V1(t,iab)|^2|/{(alpha(iab) + wb(t)}
!
!        thresd = thres
        thresd = 1.0D-08
        thres = 1.0D-08

        e2 = 0.0d+00
        e2g= 0.0d+00
        dimn = 0
        syma = 0
        indt=0
        write(*,*)' ENTER solv G part'
        write(*,*)' nsymrpa', nsymrpa


        i0 = 0
        Do ia = 1, nsec              
           ja = ia+ninact+nact
           Do ib = 1, ia-1
              jb = ib+ninact+nact
              Do ii = 1, ninact
                 ji = ii
                 i0 = i0 + 1
              End do
           End do
        End do

        nabi = i0
        Allocate(iabi(ninact+nact+1:ninact+nact+nsec,ninact+nact+1:ninact+nact+nsec,1:ninact))
        iabi = 0
        Allocate(ia0(nabi)) 
        Allocate(ib0(nabi)) 
        Allocate(ii0(nabi)) 

        i0 = 0
        Do ia = 1, nsec              
           ja = ia+ninact+nact
           Do ib = 1, ia-1
              jb = ib+ninact+nact
              Do ii = 1, ninact
                 ji = ii
                 i0 = i0 + 1
                 iabi(ja, jb, ji) = i0
                 iabi(jb, ja, ji) = i0
                 ia0(i0) = ja
                 ib0(i0) = jb
                 ii0(i0) = ji
              End do
           End do
        End do

        Allocate(v(nabi, ninact+1:ninact+nact))
        v = 0.0d+00

        write(*,*)'come'
        Call vGmat_ord (nabi, iabi, v)


     Do isym = 1, nsymrpa
        
        dimn = 0
        Do it = 1, nact
           jt = it + ninact
           if (irpmo(jt) == isym) then
              dimn = dimn + 1
              indt(dimn) = it
           End if
        End do                  ! it

        write(*,*)'isym, dimn',isym, dimn

        If (dimn == 0) goto 1000

        Allocate(sc(dimn,dimn))
        sc = 0.0d+00            ! sc N*N

       Call sGmat (dimn, indt(1:dimn), sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           write(*,*)'sG matrix is obtained normally'

           Allocate(ws(dimn))
           
           cutoff = .TRUE.
!           thresd = 1.0d-15

           Allocate(sc0(dimn,dimn))
           sc0 = sc

       Call cdiag (sc, dimn, dimm, ws, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           write(*,*)'after s cdiag, new dimension is', dimm

           If(dimm == 0) then
              deallocate(sc0)
              deallocate(sc)
              deallocate(ws)
              goto 1000
           Endif

   If(debug) then

           write(*,*)'Check whether U*SU is diagonal'

       Call checkdgc(dimn, sc0, sc, ws)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           write(*,*)'Check whether U*SU is diagonal END'

   End if

           Allocate(bc(dimn,dimn))                                 ! bc N*N
           bc = 0.0d+00

       Call bGmat (dimn, sc0, indt(1:dimn), bc)               
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           write(*,*)'bC matrix is obtained normally'

           deallocate (sc0)

           write(*,*)'OK cdiag',dimn,dimm

           Allocate(uc(dimn,dimm))                                 ! uc N*M
           Allocate(wsnew(dimm))                                  ! wnew M
           uc(:,:) = 0.0d+00
           wsnew(:) = 0.0d+00

       Call ccutoff (sc, ws, dimn, dimm, uc, wsnew)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           write(*,*)'OK ccutoff'
           deallocate (ws)
           deallocate (sc)

       Call ucramda_s_half (uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           deallocate(wsnew)

           write(*,*)'ucrams half OK'
           Allocate(bc0(dimm, dimn))                       ! bc0 M*N
           bc0 = 0.0d+00
           bc0 = MATMUL(TRANSPOSE(DCONJG(uc)), bc)
           Allocate(bc1(dimm, dimm))                      ! bc1 M*M 
           bc1 = 0.0d+00
           bc1 = MATMUL(bc0, uc)

   If (debug) then

           write(*,*)'Check whether bc1 is hermite or not'
           Do i = 1, dimm
              Do j = i, dimm
                 if(ABS(bc1(i,j)-DCONJG(bc1(j,i))) > 1.0d-6) then
                    write(*,'(2I4,2E15.7)')i,j,bc1(i,j)-bc1(j,i)
                 End if
              End do
           End do
           write(*,*)'Check whether bc1 is hermite or not END'

   End if

           deallocate (bc)
           deallocate (bc0)

           cutoff = .FALSE.

           Allocate(wb(dimm))

           write(*,*)'bC matrix is transrated to bc1(M*M matrix)!'

           Allocate(bc0(dimm,dimm))
           bc0 = bc1

       Call cdiag(bc1, dimm, dammy, wb, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   If (debug) then

           write(*,*)'Check whether bc is really diagonalized or not'

       Call checkdgc(dimm, bc0, bc1, wb)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           write(*,*)'Check whether bc is really diagonalized or not END'

   End if

           deallocate(bc0)

           write(*,*)'bC1 matrix is diagonalized!'

           e2 = 0.0d+00



        Do i0 = 1, nabi            
           ja = ia0(i0)
           jb = ib0(i0)
           ji = ii0(i0)
                 
           syma = MULTB2(isym, nsymrpa + 1)
           syma = MULTB(irpmo(jb), syma)
           syma = MULTB2(irpmo(ji), syma)
           syma = MULTB(irpmo(ja), syma)
  
           If(nsymrpa==1.or.(nsymrpa/=1.and.(syma == nsymrpa + 1))) then

              Allocate(vc(dimn))

              Do it = 1, dimn
                 vc(it) = v(i0,indt(it)+ninact)
              Enddo

              Allocate(vc1(dimm))
              vc1 = 0.0d+00

              vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn,1:dimm))),vc(1:dimn))
              Deallocate (vc)

              alpha = -eps(ji) + eps(ja) + eps(jb) - e0 + eshift  ! For Level Shift (2007/2/9)

              vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm,1:dimm))),vc1(1:dimm))

              Do j = 1, dimm
                 e = (ABS(vc1(j))**2.0d+00)/(alpha+wb(j))
                 sumc2local = sumc2local + e/(alpha+wb(j))
                 e2(isym) = e2(isym) - e
              End do

              deallocate(vc1)

           End if
                
        End do

           

           deallocate(uc)
           deallocate(wb)
           Deallocate (bc1)

 1000      write(*,'("e2g(",I3,") = ",E20.10,"a.u.")')isym,e2(isym)
           e2g = e2g + e2(isym)

        End do                  ! isym

        write(*,'("e2g      = ",E20.10,"a.u.")')e2g

        write(*,'("sumc2,g  = ",E20.10)')sumc2local
        sumc2 = sumc2 + sumc2local


        deallocate(iabi)
        deallocate(ia0)
        deallocate(ib0)
        deallocate(ii0)
        deallocate(v)



      continue 
      write(*,*)'end solvg_ord'
   end





! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE sGmat(dimn, indt, sc) ! Assume C1 molecule, overlap matrix S in space C


!  S(u,t) = <0|Eut|0>
!    
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
     
        integer, intent(in)      :: dimn, indt(dimn)
        complex*16, intent(out)  :: sc(dimn,dimn)

        real*8  ::a,b

        integer :: it, iu
        integer :: i, j



        sc = 0.0d+00

        Do i = 1, dimn
           it = indt(i)

           Do j = i, dimn
              iu = indt(j)
              a = 0.0d+0
              b = 0.0d+0

              Call dim1_density &
              (it, iu, a, b)

              sc(i,j) =  DCMPLX(a,b)
              sc(j,i) =  DCMPLX(a,-b)
!              write(*,*)i,j,sc(i,j)
           End do               !j
        End do                  !i

        End subroutine sGmat
        


! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE bGmat (dimn, sc, indt, bc) ! Assume C1 molecule, overlap matrix B in space C
!
!
!  S(u,t) = <0|Eut|0>
!
!  B(u,t) = Siguma_w [eps(w)<0|EutEww|0>] + S(u,t)(-eps(t))
!
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
    
        integer :: it, iu, iw, jt, ju, jw
        integer :: i, j
        
        integer, intent(in) :: dimn, indt(dimn)
        complex*16, intent(in)  :: sc(dimn,dimn)
        complex*16, intent(out) :: bc(dimn,dimn)

        real*8              :: denr, deni
        complex*16          :: den


        bc(:,:) = 0.0d+00

        write(*,*)'G space Bmat iroot=',iroot


        Do i = 1, dimn
           iu = indt(i)
           ju = iu + ninact

           Do j = i, dimn
              it = indt(j)
              jt = it + ninact

!  B(u,t) = Siguma_w [eps(w)<0|EutEww|0>] + S(u,t)(-eps(t))

              Do iw = 1, nact
                 jw = iw + ninact

                 Call dim2_density &
                 (iu, it, iw, iw, denr, deni)
                 den = DCMPLX(denr, deni)
                 bc(i,j) = bc(i,j) + den*eps(jw)

              End do
   
!              bc(i, j) = bc(i, j) - sc(i, j)*eps(jt)
              bc(i, j) = bc(i, j) - sc(i, j)*eps(ju)

!              write(*,*)'bc',i,j, bc(i,j)
              bc(j, i) = DCONJG(bc(i, j))


           End do               !i
        End do                  !j


        write(*,*)'bGmat is ended'

   End subroutine bGmat



! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE   vGmat_ord (nabi, iabi, v)
!
!  V(t,iab)   =  [SIGUMA_p:active <0|Etp|0>{(ai|bp)-(ap|bi)}]
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
     
        
        integer, intent(in)     :: nabi, &

        & iabi(ninact+nact+1:ninact+nact+nsec,ninact+nact+1:ninact+nact+nsec,1:ninact)

        complex*16, intent(out) :: v(nabi, ninact+1:ninact+nact)

        real*8                  :: dr, di, signij, signkl
        complex*16              :: cint2, dens

        integer :: i, j, k, l, tabi
        integer :: it, jt, il

        v = 0.0d+00

!  V(t,iab)   =  [SIGUMA_p:active <0|Etp|0>{(ai|bp)-(ap|bi)}]       a > b

        open(1, file ='Gint', status='old', form='unformatted')  !  (31|32) stored
 30     read(1, err=10, end=20) i,j,k,l,cint2

        if(i == k) goto 30
!        write(*,*) i,j,k,l,tabi,cint2

        tabi = iabi(i, k, j)

        if (i < k) then
           cint2 = -1.0d+00*cint2
        endif

        il = l - ninact

        Do it = 1, nact
           jt = ninact+it
           Call dim1_density (it, il, dr, di)
           dens = DCMPLX(dr, di)
           v(tabi,jt) = v(tabi, jt) + cint2*dens
        End do                  ! it

        goto 30

 20     close(1) ; goto 100

 10     write(*,*) 'error while opening file Gint' ; goto 100

 100    write(*,*)'vGmat_ord is ended'

   end subroutine vGmat_ord



