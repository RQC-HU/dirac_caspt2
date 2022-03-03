! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE solvG_ord_ty (e0, e2g)

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
        integer                  :: j, i, k, i0, syma, symb, isym, indt(1:nact)
        integer                  :: ia, it, ib, ii, ja, jt, jb, ji
        integer, allocatable     :: ia0(:), ib0(:), ii0(:), iabi(:,:,:)
        integer                  :: nabi


        real*8  :: thresd
        integer :: datetmp0, datetmp1, tsectmp0, tsectmp1

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
        indt=0
        if (rank == 0) then ! Process limits for output
            write(*,*)' ENTER solv G part'
            write(*,*)' nsymrpa', nsymrpa
        end if
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

        if (rank == 0) then ! Process limits for output
            write(*,*)'come'
        end if
        Call vGmat_ord_ty (nabi, iabi, v)
        Call timing(date1, tsec1, date0, tsec0)
        date1 = date0
        tsec1 = tsec0
        tsectmp0=tsec0; tsectmp1=tsec1; datetmp0=date0; datetmp1=date1

     Do isym = 1, nsymrpa

        dimn = 0
        Do it = 1, nact
           jt = it + ninact
           if (irpmo(jt) == isym) then
              dimn = dimn + 1
              indt(dimn) = it
           End if
        End do                  ! it

        if (rank == 0) then ! Process limits for output
            write(*,*)'isym, dimn',isym, dimn
        end if
        If (dimn == 0) goto 1000

        Allocate(sc(dimn,dimn))
        sc = 0.0d+00            ! sc N*N
        if (rank == 0) then ! Process limits for output
            write (*, *) 'before sGmat'
        end if
        Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
        datetmp1 = datetmp0
        tsectmp1 = tsectmp0
       Call sGmat (dimn, indt(1:dimn), sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           if (rank == 0) then ! Process limits for output
               write(*,*)'sG matrix is obtained normally'
           end if
           Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
           datetmp1 = datetmp0
           tsectmp1 = tsectmp0
           Allocate(ws(dimn))

           cutoff = .TRUE.
!           thresd = 1.0d-15

           Allocate(sc0(dimn,dimn))
           sc0 = sc
           if (rank == 0) then ! Process limits for output
               write (*, *) 'before cdiag'
           end if
           Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
           datetmp1 = datetmp0
           tsectmp1 = tsectmp0
       Call cdiag (sc, dimn, dimm, ws, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           if (rank == 0) then ! Process limits for output
               write(*,*)'after s cdiag, new dimension is', dimm
           end if
           Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
           datetmp1 = datetmp0
           tsectmp1 = tsectmp0
           If(dimm == 0) then
              deallocate(sc0)
              deallocate(sc)
              deallocate(ws)
              goto 1000
           Endif

   If(debug) then

         if (rank == 0) then ! Process limits for output
            write(*,*)'Check whether U*SU is diagonal'
         end if

       Call checkdgc(dimn, sc0, sc, ws)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (rank == 0) then ! Process limits for output
            write(*,*)'Check whether U*SU is diagonal END'
         end if

   End if

           Allocate(bc(dimn,dimn))                                 ! bc N*N
           bc = 0.0d+00
           if (rank == 0) then ! Process limits for output
               write (*, *) 'before bGmat'
           end if
           Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
           datetmp1 = datetmp0
           tsectmp1 = tsectmp0
       Call bGmat (dimn, sc0, indt(1:dimn), bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           if (rank == 0) then ! Process limits for output
               write(*,*)'bC matrix is obtained normally'
           end if
           Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
           datetmp1 = datetmp0
           tsectmp1 = tsectmp0
           deallocate (sc0)

           if (rank == 0) then ! Process limits for output
               write(*,*)'OK cdiag',dimn,dimm
           end if
           Allocate(uc(dimn,dimm))                                 ! uc N*M
           Allocate(wsnew(dimm))                                  ! wnew M
           uc(:,:) = 0.0d+00
           wsnew(:) = 0.0d+00
           if (rank == 0) then ! Process limits for output
               write (*, *) 'before ccutoff'
           end if
           Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
           datetmp1 = datetmp0
           tsectmp1 = tsectmp0
       Call ccutoff (sc, ws, dimn, dimm, uc, wsnew)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           if (rank == 0) then ! Process limits for output
               write(*,*)'OK ccutoff'
           end if
           Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
           datetmp1 = datetmp0
           tsectmp1 = tsectmp0
           deallocate (ws)
           deallocate (sc)
           if (rank == 0) then ! Process limits for output
               write (*, *) 'before ucramda_s_half'
           end if
           Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
           datetmp1 = datetmp0
           tsectmp1 = tsectmp0
       Call ucramda_s_half (uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           deallocate(wsnew)

           if (rank == 0) then ! Process limits for output
               write(*, *)'ucrams half OK'
           end if
           Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
           datetmp1 = datetmp0
           tsectmp1 = tsectmp0
           Allocate(bc0(dimm, dimn))                       ! bc0 M*N
           bc0 = 0.0d+00
           bc0 = MATMUL(TRANSPOSE(DCONJG(uc)), bc)
           Allocate(bc1(dimm, dimm))                      ! bc1 M*M
           bc1 = 0.0d+00
           bc1 = MATMUL(bc0, uc)

   If (debug) then

         if (rank == 0) then ! Process limits for output
            write(*,*)'Check whether bc1 is hermite or not'
            Do i = 1, dimm
               Do j = i, dimm
                  if(ABS(bc1(i,j)-DCONJG(bc1(j,i))) > 1.0d-6) then
                     write(*,'(2I4,2E15.7)')i,j,bc1(i,j)-bc1(j,i)
                  End if
               End do
            End do
            write(*,*)'Check whether bc1 is hermite or not END'
         end if
   End if

           deallocate (bc)
           deallocate (bc0)

           cutoff = .FALSE.

           Allocate(wb(dimm))

           if (rank == 0) then ! Process limits for output
               write(*,*)'bC matrix is transrated to bc1(M*M matrix)!'
           end if
           Allocate(bc0(dimm,dimm))
           bc0 = bc1
           if (rank == 0) then ! Process limits for output
               write (*, *) 'before cdiag'
           end if
           Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
           datetmp1 = datetmp0
           tsectmp1 = tsectmp0
       Call cdiag(bc1, dimm, dammy, wb, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       if (rank == 0) then ! Process limits for output
            write (*, *) 'end cdiag'
       end if
       Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
       datetmp1 = datetmp0
       tsectmp1 = tsectmp0
   If (debug) then

       if (rank == 0) then ! Process limits for output
         write(*,*)'Check whether bc is really diagonalized or not'
       end if
       Call checkdgc(dimm, bc0, bc1, wb)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       if (rank == 0) then ! Process limits for output
         write(*,*)'Check whether bc is really diagonalized or not END'
       end if

   End if

           deallocate(bc0)

           if (rank == 0) then ! Process limits for output
               write(*,*)'bC1 matrix is diagonalized!'
           end if

           e2 = 0.0d+00



        Do i0 = 1, nabi
           ja = ia0(i0)
           jb = ib0(i0)
           ji = ii0(i0)

!     EaiEbt|0>

           syma = MULTB_D (irpmo(jb), isym)
           symb = MULTB_D (irpmo(ja), irpmo(ji))
           syma = MULTB_S (syma, symb)

           If(nsymrpa==1.or.(nsymrpa/=1.and.(syma == 1))) then

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

 1000      if (rank == 0) write(*,'("e2g(",I3,") = ",E20.10,"a.u.")')isym,e2(isym)
           e2g = e2g + e2(isym)
           if (rank == 0) then ! Process limits for output
               write (*, *) 'End e2(isym) add'
           end if
           Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
           datetmp1 = datetmp0
           tsectmp1 = tsectmp0
        End do                  ! isym

        if (rank == 0) then ! Process limits for output
            write(*,'("e2g      = ",E20.10,"a.u.")')e2g
            write(*,'("sumc2,g  = ",E20.10)')sumc2local
        end if
        sumc2 = sumc2 + sumc2local


        deallocate(iabi)
        deallocate(ia0)
        deallocate(ib0)
        deallocate(ii0)
        deallocate(v)



      continue
      if (rank == 0) then ! Process limits for output
         write(*,*)'end solvG_ord_ty'
      end if
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
        include 'mpif.h'

        integer :: it, iu, iw, jt, ju, jw
        integer :: i, j

        integer, intent(in) :: dimn, indt(dimn)
        complex*16, intent(in)  :: sc(dimn,dimn)
        complex*16, intent(out) :: bc(dimn,dimn)

        real*8              :: denr, deni
        complex*16          :: den


        bc(:,:) = 0.0d+00

        if (rank == 0) then ! Process limits for output
            write(*,*)'G space Bmat iroot=',iroot
        end if

        !$OMP parallel do private(iu,ju,j,it,jt,iw,jw,denr,deni,den)
        Do i = rank + 1, dimn, nprocs
         ! Do i = 1, dimn
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
        if (rank == 0) then
         call MPI_Reduce (MPI_IN_PLACE, bc(1,1), dimn**2, MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        else
         call MPI_Reduce (bc(1,1), bc(1,1), dimn**2, MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        end if

        if (rank == 0) then ! Process limits for output
            write(*,*)'bGmat is ended'
        end if
   End subroutine bGmat



! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE   vGmat_ord_ty (nabi, iabi, v)
!
!  V(t,iab)   =  [SIGUMA_p:active <0|Etp|0>{(ai|bp)-(ap|bi)}]
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

      use four_caspt2_module

      Implicit NONE
      include 'mpif.h'

      integer, intent(in)     :: nabi, &

      & iabi(ninact+nact+1:ninact+nact+nsec,ninact+nact+1:ninact+nact+nsec,1:ninact)

      complex*16, intent(out) :: v(nabi, ninact+1:ninact+nact)

      real*8                  :: dr, di, signij, signkl
      complex*16              :: cint2, dens

      integer :: i, j, k, l, tabi
      integer :: it, jt, il

      v = 0.0d+00

!  V(t,iab)   =  [SIGUMA_p:active <0|Etp|0>{(ai|bp)-(ap|bi)}]       a > b

      !   open(1, file ='Gint', status='old', form='unformatted')  !  (31|32) stored
        open(1, file =gint, status='old', form='formatted')  !  (31|32) stored
 30     read(1,  '(4I4, 2e20.10)', err=10, end=20) i,j,k,l,cint2

        if(i == k) goto 30
!        write(*,*) i,j,k,l,tabi,cint2

        tabi = iabi(i, k, j)

        if (i < k) then
           cint2 = -1.0d+00*cint2
        endif

        il = l - ninact
        !$OMP parallel do private(it,jt,dr,di,dens)
        Do it = 1, nact
           jt = ninact+it
           Call dim1_density (it, il, dr, di)
           dens = DCMPLX(dr, di)
           v(tabi,jt) = v(tabi, jt) + cint2*dens
        End do                  ! it

        goto 30

 20     close(1) ; goto 100

 10     write(*,*) 'error while opening file Gint' ; goto 100

 100    if(rank == 0)  write(*,*)'vGmat_ord_ty is ended'

      !   v(nabi, ninact+1:ninact+nact)
       call MPI_Allreduce(MPI_IN_PLACE, v(1, ninact + 1), nabi*nact, &
                          MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
   end subroutine vGmat_ord_ty
