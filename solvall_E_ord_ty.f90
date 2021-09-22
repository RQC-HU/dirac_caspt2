! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE solvE_ord_ty (e0, e2e)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE

        real*8,     intent(in) :: e0
        real*8,     intent(out):: e2e



        integer :: dimn, dimm, count, dammy

        integer, allocatable :: indsym(:,:)

        real*8, allocatable  :: sr(:,:), ur(:,:)
        real*8, allocatable  :: br(:,:), wsnew(:), ws(:), wb(:)
        real*8, allocatable  :: br0(:,:), br1(:,:)
        real*8               :: e2(2*nsymrpa), alpha, e


        complex*16, allocatable  :: sc(:,:), uc(:,:), sc0(:,:)
        complex*16, allocatable  :: bc(:,:)
        complex*16, allocatable  :: bc0(:,:), bc1(:,:), v(:,:), vc(:), vc1(:)

        logical :: cutoff
        integer :: j, i, k, syma, symb, isym, indt(1:nact)
        integer :: ia, it, ij, ii, ja, jt, jj, ji

        integer :: i0
        integer, allocatable     :: ia0(:), ii0(:), ij0(:), iaij(:,:,:)
        integer                  :: naij

        real*8  :: thresd

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE E IS NOW CALCULATED
!
!     EtiEaj|0>
!
!   DRAS1 =-2   DRAS2 = +1   DRAS3 = +1
!
!   i > j
!
!  S(ukbl,tiaj) = d(ki) d(lj) d(ba) [d(ut) - <0|Etu|0>]  <= S(u,t)
!                                   ~~~~~~~~~~~~~~~~~~~
!  S(u,t) = d(ut) - <0|Etu|0>
!
!  B(u,t) = Siguma_w [eps(w){d(tu)<0|Eww|0>-<0|EtuEww|0>}] + S(u,t)eps(t)
!
!         = - Siguma_w [eps(w)<0|EtuEww|0>] + d(tu)e0 + S(u,t)eps(t)
!
!  alpha(i,j,a) = -eps(i) - eps(j) + eps(a) - e0
!
!  where
!
!  e0 = Siguma_w [eps(w)<0|Eww|0>] (<== calculated as e0 in calce0.f)
!
!  V(t,ija)   =[SIGUMA_p:active <0|Ept|0>{(ai|pj) - (aj|pi)}] + (aj|ti) - (ai|tj)
!
!  E2 = SIGUMA_iab, dimm |V1(t,ija)|^2|/{(alpha(ija) + wb(t)}
!
!        thresd = thres
        thresd = 1.0D-08
        thres = 1.0D-08

        e2 = 0.0d+00
        e2e= 0.0d+00
        dimn = 0
        syma = 0
        indt=0
        write(*,*)' ENTER solv E part'
        write(*,*)' nsymrpa', nsymrpa

        i0 = 0
        Do ia = 1, nsec
           Do ii = 1, ninact
              Do ij = 1, ii-1                ! i > j
                 i0 = i0 + 1
              End do
           End do
        End do

        naij = i0
        Allocate(iaij(ninact+nact+1:ninact+nact+nsec,1:ninact,1:ninact))
        iaij = 0
        Allocate(ia0(naij))
        Allocate(ii0(naij))
        Allocate(ij0(naij))

        i0 = 0
        Do ia = 1, nsec
           ja = ia+ninact+nact
           Do ii = 1, ninact
              ji = ii
              Do ij = 1, ii-1                ! i > j
                 jj = ij
                 i0 = i0 + 1
                 iaij(ja, ji, jj) = i0
                 iaij(ja, jj, ji) = i0
                 ia0(i0) = ja
                 ii0(i0) = ji
                 ij0(i0) = jj
              End do
           End do
        End do

       Allocate(v(naij, ninact+1:ninact+nact))
        v = 0.0d+00

        Call vEmat_ord_ty (naij, iaij, v)
        write(*,*)'come'



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

        If(dimn == 0) goto 1000

        Allocate(sc(dimn,dimn))
        sc = 0.0d+00            ! sc N*N

       Call sEmat (dimn, indt(1:dimn), sc)
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

       Call bEmat (e0, dimn, sc0, indt(1:dimn), bc)
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

   If(debug) then

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

        Do i0 = 1, naij
           ja = ia0(i0)
           ji = ii0(i0)
           jj = ij0(i0)

!     EtiEaj|0>

           syma = MULTB_D (irpmo(ja),irpmo(jj))
           symb = MULTB_D (isym, irpmo(ji))
           syma = MULTB_S (symb, syma)

           If(nsymrpa==1.or.(nsymrpa/=1.and.(syma == 1))) then

              Allocate(vc(dimn))

              Do it = 1, dimn
                 vc(it) = v(i0,indt(it)+ninact)
              Enddo

              Allocate(vc1(dimm))
              vc1 = 0.0d+00

              vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn,1:dimm))),vc(1:dimn))
              Deallocate (vc)

              alpha = + eps(ja) -eps(ji) - eps(jj) - e0 + eshift  ! For Level Shift (2007/2/9)

              vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm,1:dimm))),vc1(1:dimm))

              Do j = 1, dimm
                 e = DCONJG(vc1(j))*vc1(j)/(alpha+wb(j))
                 sumc2local = sumc2local + e/(alpha+wb(j))
                 e2(isym) = e2(isym) - e
              End do

              deallocate(vc1)

           End if

        End do



           deallocate(uc)
           deallocate(wb)
           Deallocate (bc1)

 1000      write(*,'("e2e(",I3,") = ",E20.10,"a.u.")')isym,e2(isym)
           e2e = e2e + e2(isym)

        End do                  ! isym

        write(*,'("e2e      = ",E20.10,"a.u.")')e2e

        write(*,'("sumc2,e  = ",E20.10)')sumc2local
        sumc2 = sumc2 + sumc2local

        deallocate(iaij)
        deallocate(ia0)
        deallocate(ii0)
        deallocate(ij0)
        deallocate(v)




      continue
      write(*,*)'end solveE_ord_ty'
   end





! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE sEmat(dimn, indt, sc) ! Assume C1 molecule, overlap matrix S in space E


! S(u,t) = d(ut) - <0|Etu|0>
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
           iu = indt(i)

           Do j = i, dimn
              it = indt(j)
              a = 0.0d+0
              b = 0.0d+0

              Call dim1_density &
              (it, iu, a, b)

              If(iu == it) then
                 sc(i,j) = 1  - DCMPLX(a,b)
              Else
                 sc(i,j) = -DCMPLX(a,b)
              Endif

              sc(j,i) =  DCONJG(sc(i,j))

!              write(*,*)i,j,sc(i,j)

           End do               !j
        End do                  !i

        End subroutine sEmat



! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE bEmat (e0, dimn, sc, indt, bc) ! Assume C1 molecule, overlap matrix B in space E
!
!
!  S(u,t) = d(ut) - <0|Etu|0>
!
!  B(u,t) = Siguma_w [eps(w){d(tu)<0|Eww|0>-<0|EtuEww|0>}] + S(u,t)eps(t)
!
!         = - Siguma_w [eps(w)<0|EtuEww|0>] + d(tu)e0 + S(u,t)eps(t)
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
        real*8, intent(in)      :: e0

        real*8              :: denr, deni
        complex*16          :: den


        bc(:,:) = 0.0d+00

        write(*,*)'E space Bmat iroot=',iroot


        Do i = 1, dimn
           iu = indt(i)
           ju = iu + ninact

           Do j = i, dimn
              it = indt(j)
              jt = it + ninact

              Do iw = 1, nact
                 jw = iw + ninact

!         = - Siguma_w [eps(w)<0|EtuEww|0>] + d(tu)e0 + S(u,t)eps(t)

                 Call dim2_density &
                 (it, iu, iw, iw, denr, deni)
                 den = DCMPLX(denr, deni)
                 bc(i,j) = bc(i,j) - den*eps(jw)

              End do

              if(it == iu) bc(i, j) = bc(i, j) + e0

              bc(i, j) = bc(i, j) + sc(i, j)*eps(jt)

!              write(*,*)'bc',i,j, bc(i,j)

              bc(j, i) = DCONJG(bc(i, j))


           End do               !i
        End do                  !j


        write(*,*)'bEmat is ended'

   End subroutine bEmat



! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE vEmat_ord_ty (naij, iaij, v)
!
!  V(t,ija)   =[SIGUMA_p:active <0|Ept|0>{(ai|pj) - (aj|pi)}] + (aj|ti) - (ai|tj)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE

        integer, intent(in)     :: naij, &

        & iaij(ninact+nact+1:ninact+nact+nsec,1:ninact,1:ninact)

        complex*16, intent(out) :: v(naij, ninact+1:ninact+nact)

        real*8                  :: dr, di
        complex*16              :: cint2, dens

        integer :: i, j, k, l, taij
        integer :: it, jt, ik

        v = 0.0d+00

!  V(t,ija)   =[SIGUMA_p:active <0|Ept|0>{(ai|pj) - (aj|pi)}] - (ai|tj) + (aj|ti)   i > j

      !   open(1, file ='Eint', status='old', form='unformatted')  !  (31|21) stored
        open(1, file =eint, status='old', form='unformatted')  !  (31|21) stored
 30          read(1, err=10, end=20) i,j,k,l,cint2

        if(j == l) goto 30

        taij = iaij(i, j, l)
        ik = k - ninact

!        write(*,*) i,j,k,l,taij,cint2

        if (j < l) then
           cint2 = -1.0d+00*cint2
        endif

        v(taij,k) = v(taij, k) - cint2

        Do it = 1, nact
           jt = ninact+it
           Call dim1_density (it, ik, dr, di)          ! ik corresponds to p in above formula
           dens = DCMPLX(dr, di)
           v(taij,jt) = v(taij, jt) + cint2*dens
        End do                  ! it

        if (j < l) then
           cint2 = -1.0d+00*cint2          ! data cint2 becomes initial values!
        endif

!! Take Kramers conjugate !
!
!        Call takekr( i, j, k, l, cint2)
!
!        taij = iaij(i, j, l)
!        ik = k - ninact
!
!!        write(*,*) i,j,k,l,taij,cint2
!
!        if (j < l) then
!           cint2 = -1.0d+00*cint2
!        endif
!
!        v(taij,k) = v(taij, k) - cint2
!
!        Do it = 1, nact
!           jt = ninact+it
!           Call dim1_density (it, ik, dr, di)          ! ik corresponds to p in above formula
!           dens = DCMPLX(dr, di)
!           v(taij,jt) = v(taij, jt) + cint2*dens
!        End do                  ! it

        goto 30

 20          close(1) ; goto 100

 10               write(*,*) 'error while opening file Eint' ; goto 100

 100                  write(*,*)'vEmat_ord_ty is ended'

   end subroutine vEmat_ord_ty
