! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE solvF_ord (e0, e2f)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE

        real*8,     intent(in) :: e0
        real*8,     intent(out):: e2f



        integer :: dimn, dimm, count, dammy

        integer, allocatable :: indsym(:,:)

        real*8, allocatable  :: sr(:,:), ur(:,:)
        real*8, allocatable  :: br(:,:), wsnew(:), ws(:), wb(:)
        real*8, allocatable  :: br0(:,:), br1(:,:)
        real*8               :: e2(2*nsymrpa), alpha, e


        complex*16, allocatable  :: sc(:,:), uc(:,:), sc0(:,:)
        complex*16, allocatable  :: bc(:,:)
        complex*16, allocatable  :: bc0(:,:), bc1(:,:), v(:,:,:), vc(:), vc1(:)

        logical :: cutoff
        integer :: j, i, k, syma, isym, i0, j0
        integer :: ia, it, ib, iu, ja, jt, jb, ju

        integer, allocatable     :: ia0(:), ib0(:), iab(:,:)
        integer                  :: nab


        real*8  :: thresd

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE F IS NOW CALCULATED
!      
!     EatEbu|0>
!
!   DRAS1 =0   DRAS2 = -2   DRAS3 = +2
!
!   a > b t > u ( c > d, v > x)
!
!  S(cvdx,atbu) = d(ac) d(bd)  [ <0|EvtExu|0> - d(xt)<0|Evu|0>]  <= S(vx,tu)
!                               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  S(vx, tu) = <0|EvtExu|0> - d(xt)<0|Evu|0>
!
!  B(vx, tu) = Siguma_w [eps(w){ <0|EvtExuEww|0> - d(xt)<0|EvuEww|0>}] + S(u,t){-eps(u)-eps(t)}
!
!  alpha(a, b) = + eps(a) + eps(b) - e0
!
!  where
!
!  e0 = Siguma_w [eps(w)<0|Eww|0>] (<== calculated as e0 in calce0.f)
!
!  V(tu,ab)   = SIGUMA_p,q:active <0|EtpEuq|0>(ap|bq) -  SIGUMA_p:active <0|Etp|0>(au|bp)
!
!  E2 = SIGUMA_iab,t:dimm |V1(t,ab)|^2|/{(alpha(ab) + wb(t)}
!
!        thresd = thres
        thresd = 1.0D-08
        thres = 1.0D-08

        e2 = 0.0d+00
        e2f= 0.0d+00
        dimn = 0
        syma = 0
        write(*,*)' ENTER solv F part'
        write(*,*)' nsymrpa', nsymrpa

        i0 = 0
        Do ia = 1, nsec
           Do ib = 1, ia-1
              i0 = i0 + 1
           End do
        End do

        nab = i0
        Allocate(iab(ninact+nact+1:ninact+nact+nsec,ninact+nact+1:ninact+nact+nsec))
        iab = 0
        Allocate(ia0(nab))
        Allocate(ib0(nab))

        i0 = 0
        Do ia = 1, nsec
           ja = ia+ninact+nact
           Do ib = 1, ia-1
              jb = ib+ninact+nact
              i0 = i0 + 1
              iab(ja, jb) = i0
              iab(jb, ja) = i0
              ia0(i0) = ja
              ib0(i0) = jb
           End do
        End do


       Allocate(v(nab, ninact+1:ninact+nact, ninact+1:ninact+nact))
        v = 0.0d+00

        Call vFmat_ord (nab, iab, v)
        write(*,*)'come'  

!Iwamuro Modify
     Do Isym = Nsymrpa+1, 2*Nsymrpa
!     Do Isym = 1, Nsymrpa
        
        Dimn = 0
        Do it = 1, nact
           jt = it + ninact
!           Do iu = 1, nact
           Do iu = 1, it-1
              ju = iu + ninact

              syma = nsymrpa + 1
              syma = MULTB2( irpmo(ju), syma)   ! HERE ABELIAN IS ASSUMED. AT PRESENT I'M NOT SURE
              syma = MULTB2( irpmo(jt), syma)   ! THIS TREATMENT IS CORRECT. THE ORDER OF OPERATOR
                                                 ! MAY BE INAFFECTED TO TOTAL POINT GROUP SYMMETRY
              if (nsymrpa == 1 .or. (nsymrpa /=1 .and. syma == isym)) then             
                 dimn = dimn + 1
              End if
           End do               ! iu
        End do                  ! it

        write(*,*)'isym, dimn',isym, dimn

        If(dimn == 0) goto 1000

        Allocate(indsym(2, dimn))
        
        dimn = 0
        Do it = 1, nact
           jt = it + ninact
!           Do iu = 1, nact
           Do iu = 1, it-1
              ju = iu + ninact

              syma = nsymrpa + 1
              syma = MULTB2( irpmo(ju), syma)   ! HERE ABELIAN IS ASSUMED. AT PRESENT I'M NOT SURE
              syma = MULTB2( irpmo(jt), syma)   ! THIS TREATMENT IS CORRECT. THE ORDER OF OPERATOR
                                                 ! MAY BE INAFFECTED TO TOTAL POINT GROUP SYMMETRY
              if (nsymrpa == 1 .or. (nsymrpa /=1 .and. syma == isym)) then             
                 dimn = dimn + 1
                 indsym(1,dimn) = it
                 indsym(2,dimn) = iu
              End if
           End do               ! iu
        End do                  ! it


        Allocate(sc(dimn,dimn))
        sc = 0.0d+00            ! sc N*N

       Call sFmat (dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           write(*,*)'sc matrix is obtained normally'

           Allocate(ws(dimn))
           
           cutoff = .TRUE.
!           thresd = 1.0d-15

           Allocate(sc0(dimn,dimn))
           sc0 = sc

       Call cdiag (sc, dimn, dimm, ws, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           write(*,*)'after s cdiag, new dimension is', dimm


           If(dimm == 0) then
              deallocate(indsym)
              deallocate(sc0)
              deallocate(sc)
              deallocate(ws)
              goto 1000
           Endif



   If (debug) then

           write(*,*)'Check whether U*SU is diagonal'

       Call checkdgc(dimn, sc0, sc, ws)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           write(*,*)'Check whether U*SU is diagonal END'

   End if

           Allocate(bc(dimn,dimn))                                 ! bc N*N
           bc = 0.0d+00

       Call bFmat (dimn, sc0, indsym, bc)               
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           write(*,*)'bc matrix is obtained normally'

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

    If(debug) then

           write(*,*)'Check whether bc is really diagonalized or not'

       Call checkdgc(dimm, bc0, bc1, wb)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           write(*,*)'Check whether bc is really diagonalized or not END'

    End if

           deallocate(bc0)

           write(*,*)'bC1 matrix is diagonalized!'

           e2 = 0.0d+00

        Do i0 = 1, nab
           ja = ia0(i0)
           jb = ib0(i0)

           syma = isym
           syma = MULTB (irpmo(jb), syma)
           syma = MULTB (irpmo(ja), syma)
  
           If(nsymrpa==1.or.(nsymrpa/=1.and.(syma == nsymrpa + 1))) then

              Allocate(vc(dimn))
              Do it = 1, dimn
                 vc(it) = v(i0,indsym(1,it)+ninact,indsym(2,it)+ninact)
              Enddo
              Allocate(vc1(dimm))
              vc1 = 0.0d+00

              vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn,1:dimm))),vc(1:dimn))
              Deallocate (vc)

              alpha = + eps(ja) + eps(jb) - e0 + eshift  ! For Level Shift (2007/2/9)

              vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm,1:dimm))),vc1(1:dimm))

              Do j = 1, dimm
                 e = (ABS(vc1(j))**2.0d+00)/(alpha+wb(j))
                 sumc2local = sumc2local + e/(alpha+wb(j))           
                 e2(isym) = e2(isym) - e
              End do

              deallocate(vc1)

           End if
                
        End do                  !i0


        deallocate(indsym)
        deallocate(uc)
        deallocate(wb)
        Deallocate (bc1)

 1000   write(*,'("e2f(",I3,") = ",E20.10,"a.u.")')isym,e2(isym)
        e2f = e2f + e2(isym)

        End do                  ! isym

        write(*,'("e2f      = ",E20.10,"a.u.")')e2f

        write(*,'("sumc2,f  = ",E20.10)')sumc2local
        sumc2 = sumc2 + sumc2local

        deallocate(iab)
        deallocate(ia0)
        deallocate(ib0)
        deallocate(v)

      continue 
      write(*,*)'end solve_ord'
   end





! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE sFmat(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space F


!  S(vx, tu) = <0|EvtExu|0> - d(xt)<0|Evu|0>
!
!     v > x, t > u

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
     
        integer, intent(in)      :: dimn, indsym(2, dimn)
        complex*16, intent(out)  :: sc(dimn,dimn)

        real*8  :: a,b

        integer :: it, iu, iv, ix
        integer :: i, j
        integer :: count



        sc = 0.0d+00

        Do i = 1, dimn
           iv  = indsym(1,i)
           ix  = indsym(2,i)
           Do j = i, dimn
              it  = indsym(1,j)
              iu  = indsym(2,j)

              a = 0.0d+0
              b = 0.0d+0
              Call dim2_density (iv, it, ix, iu, a,b)
              sc(i,j) =  DCMPLX(a,b)

              if(ix == it) then
                 a = 0.0d+0
                 b = 0.0d+0
                 Call dim1_density (iv, iu, a,b)
                 sc(i,j) =  sc(i,j) - DCMPLX(a,b)
              End if

              sc(j,i) =  DCONJG(sc(i,j))
              
           End do               !j
        End do                  !i

        End subroutine sFmat
        


! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE bFmat (dimn, sc, indsym, bc) ! Assume C1 molecule, overlap matrix B in space F
!
!  B(vx, tu) = Siguma_w [eps(w){ <0|EvtExuEww|0> - d(xt)<0|EvuEww|0>}] + S(u,t){-eps(u)-eps(t)}
!
!  v > x, t > u
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
     
        integer, intent(in) :: dimn, indsym(2,dimn)
        complex*16, intent(in)  :: sc(dimn,dimn)
        complex*16, intent(out) :: bc(dimn,dimn)

        real*8              :: e, denr, deni
        complex*16          :: den

        integer :: it, iu, iv, ix, iy, iz, iw
        integer :: jt, ju, jv, jx, jw, i, j
        

        bc(:,:) = 0.0d+00

        write(*,*)'F space Bmat iroot=',iroot

        Do i = 1, dimn

           iv = indsym(1,i)
           jv = iv + ninact
           ix = indsym(2,i)
           jx = ix + ninact

           Do j = i, dimn

              it = indsym(1,j)
              jt = it + ninact
              iu = indsym(2,j)
              ju = iu + ninact

!  B(vx, tu) = Siguma_w [eps(w){ <0|EvtExuEww|0> - d(xt)<0|EvuEww|0>}] + S(u,t){-eps(u)-eps(t)}

              e = -eps(ju) - eps(jt)

              Do iw = 1, nact
                 jw = iw + ninact

                 Call dim3_density &
                 (iv, it, ix, iu, iw, iw, denr, deni)
                 den = DCMPLX(denr, deni)
                 bc(i,j) = bc(i,j) + den*eps(jw)

                 If(ix == it) then

                    Call dim2_density &
                    (iv, iu, iw, iw, denr, deni)
                    den = DCMPLX(denr, deni)
                    bc(i,j) = bc(i,j) - den*eps(jw)

                 End if

              End do
   
              bc(i, j) = bc(i, j) + sc(i, j)*e

              bc(j, i) = DCONJG(bc(i, j))


           End do               !i
        End do                  !j


        write(*,*)'bFmat is ended'

   End subroutine bFmat



! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE vFmat_ord (nab, iab, v)
!
! V(tu, ab) =  SIGUMA_p,q:active <0|EtpEuq|0>(ap|bq) -  SIGUMA_p:active <0|Etp|0>(au|bp)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
     
        
        integer, intent(in)     :: nab, &

        & iab(ninact+nact+1:ninact+nact+nsec,ninact+nact+1:ninact+nact+nsec)

        complex*16, intent(out) :: v(nab,ninact+1:ninact+nact,ninact+1:ninact+nact)

        real*8                  :: dr, di
        complex*16              :: cint2, dens

        integer :: i, j, k, l, tab, ip, iq, save
        integer :: it, jt, ju, iu

        v = 0.0d+00

! V(ab,t,u) =  SIGUMA_p,q:active <0|EtpEuq|0>(ap|bq) -  SIGUMA_p:active <0|Etp|0>(au|bp)

        open(1, file ='Fint', status='old', form='unformatted')  !  (32|32) stored  a > b
 30     read(1, err=10, end=20) i,j,k,l,cint2

        if(i <= k) goto 30

        tab = iab(i, k)

!        if (i < k ) then   ! indices exchange i<=>k j<=>l
!           save = i
!           i    = k
!           k    = save
!           save = j
!           j    = l
!           l    = save
!        endif

!        write(*,'(4I4,2E20.10)') i,j,k,l,cint2

        ip = j - ninact
        iq = l - ninact

! V(ab,t,u) =  SIGUMA_p,q:active <0|EtpEuq|0>(ap|bq) -  SIGUMA_p:active <0|Etp|0>(au|bp)
!                                <0|EtjEul|0>(ij|kl)                             (ij|kl)
!
!                             p=j, q=l loop for t and u             u=j, p=l loop for t
!                
        Do it = 1, nact
           jt = it + ninact
           Do iu = 1, it-1
              ju = iu + ninact

              Call dim2_density (it, ip, iu, iq, dr, di)
              dens = DCMPLX(dr, di)
              v(tab,jt,ju) = v(tab,jt,ju) + cint2*dens
           End do               ! iu

           Call dim1_density (it, iq, dr, di)
           dens = DCMPLX(dr, di)
           v(tab, jt, j) = v(tab, jt, j) - cint2*dens

        End do                  ! ip

        goto 30

 20     close(1) ; goto 100

 10     write(*,*) 'error while opening file Fint' ; goto 100

 100    write(*,*)'vFmat_ord is ended'


   end subroutine vFmat_ord



