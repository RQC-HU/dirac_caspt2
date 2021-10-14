! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE solvB_ord_ty(e0, e2b)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE

       real*8, intent(in) :: e0
       real*8, intent(out):: e2b

       integer :: dimn, dimm, count, dammy

       integer, allocatable :: indsym(:, :)

       real*8, allocatable  :: sr(:, :), ur(:, :)
       real*8, allocatable  :: br(:, :), wsnew(:), ws(:), wb(:)
       real*8, allocatable  :: br0(:, :), br1(:, :)
       real*8               :: e2(2*nsymrpa), e, alpha

       complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
       complex*16, allocatable  :: bc(:, :)
       complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :), vc(:), vc1(:)

       integer, allocatable     :: ii0(:), ij0(:), iij(:, :)
       integer                  :: nij

       logical :: cutoff
       integer :: j, i, k, syma, isym, i0, j0
       integer :: ij, it, ii, iu, jj, jt, ji, ju

       real*8  :: thresd

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE B IS NOW CALCULATED
!
!     EtiEuj|0>
!
!   DRAS1 = -2   DRAS2 = 2   DRAS3 = 0
!
!   t > u, i > j
!
!
!  S(xkyl,tiuj) = d(ki)d(lj)S(xy,tu)
!
!  S(xy, tu) = <0|EtxEuy|0> -d(tx)<0|Euy|0> -d(uy)<0|Etx|0> -d(ty)<0|Eux|0> +d(tx)d(uy)-d(ty)d(ux)
!
!  B(xy,tu) = Siguma_w [eps(w){<0|EtxEuyEww|0>-d(tx)<0|EuyEww|0> -d(uy)<0|EtxEww|0> -d(ty)<0|EuxEww|0>]
!
!              +{d(tx)d(uy)-d(ty)d(ux)}*e0      +S(xy,tu){eps(t)+eps(u)}
!
!  a(i,j)       = -eps(i) - eps(j) - e0
!
! V(i,j) =  SIGUMA_p,q:active <0|EptEqu|0>(pi|qj)
!
!           + SIGUMA_p:active[<0|Ept|0> {(pj|ui) - (uj|pi)}  - <0|Epu|0> (ti|pj)]
!
!           + (uj|ti)  - (tj|ui)
!
!
!  E2 = SIGUMA_a,i, dimm |V1(dimm,ai)|^2|/{(alpha(ai) + wb(dimm)}

!        thresd = thres
       thresd = 1.0D-08
       thres = 1.0D-08

       e2 = 0.0d+00
       e2b = 0.0d+00
       dimn = 0

       if (rank == 0) then ! Process limits for output
           write (normaloutput, *) ' ENTER solv B part'
           write (normaloutput, *) ' nsymrpa', nsymrpa
       end if
       i0 = 0
       Do ii = 1, ninact
           Do ij = 1, ii - 1
               i0 = i0 + 1
           End do
       End do

       nij = i0
       Allocate (iij(ninact, ninact)); Call memminus(KIND(iij), SIZE(iij), 1)
       iij = 0
       Allocate (ii0(nij)); Call memminus(KIND(ii0), SIZE(ii0), 1)
       Allocate (ij0(nij)); Call memminus(KIND(ii0), SIZE(ii0), 1)

       i0 = 0
       Do ii = 1, ninact
           Do ij = 1, ii - 1
               i0 = i0 + 1
               iij(ii, ij) = i0
               iij(ij, ii) = i0
               ii0(i0) = ii
               ij0(i0) = ij
           End do
       End do
       Allocate (v(nij, ninact + 1:ninact + nact, ninact + 1:ninact + nact))
       Call memplus(KIND(v), SIZE(v), 2)
       v = 0.0d+00

       Call vBmat_ord_ty(nij, iij, v)

!     EtiEuj|0>

       Do isym = 1, nsymrpa

           dimn = 0
           Do it = 1, nact
               jt = it + ninact
               Do iu = 1, it - 1
                   ju = iu + ninact

                   syma = MULTB_D(irpmo(jt), irpmo(ju) - (-1)**(mod(irpmo(ju), 2)))

                   if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                       dimn = dimn + 1
                   End if
               End do               ! iu
           End do                  ! it

           if (rank == 0) then ! Process limits for output
               write (normaloutput, *) 'isym, dimn', isym, dimn
           end if
           If (dimn == 0) goto 1000

           Allocate (indsym(2, dimn)); Call memplus(KIND(indsym), SIZE(indsym), 1)

           dimn = 0
           Do it = 1, nact
               jt = it + ninact
               Do iu = 1, it - 1
                   ju = iu + ninact

                   syma = MULTB_D(irpmo(jt), irpmo(ju) - (-1)**(mod(irpmo(ju), 2)))

                   if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                       dimn = dimn + 1
                       indsym(1, dimn) = it
                       indsym(2, dimn) = iu
                   End if
200            End do               ! iu
           End do                  ! it

           Allocate (sc(dimn, dimn)); Call memplus(KIND(sc), SIZE(sc), 2)
           sc = 0.0d+00            ! sc N*N

           Call sBmat(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           if (rank == 0) then ! Process limits for output
               write (normaloutput, *) 'sc matrix is obtained normally'
           end if
           Allocate (ws(dimn)); Call memplus(KIND(ws), SIZE(ws), 1)

           cutoff = .TRUE.
           thresd = 1.0d-08

           Allocate (sc0(dimn, dimn)); Call memplus(KIND(sc0), SIZE(sc0), 2)
           sc0 = sc

           Call cdiag(sc, dimn, dimm, ws, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           if (rank == 0) then ! Process limits for output
               write (normaloutput, *) 'after s cdiag, new dimension is', dimm
           end if
           If (dimm == 0) then
               deallocate (indsym); Call memminus(KIND(indsym), SIZE(indsym), 1)
               deallocate (sc0); Call memminus(KIND(sc0), SIZE(sc0), 2)
               deallocate (sc); Call memminus(KIND(sc), SIZE(sc), 2)
               deallocate (ws); Call memminus(KIND(ws), SIZE(ws), 1)
               goto 1000
           End if

           Allocate (bc(dimn, dimn)); Call memplus(KIND(bc), SIZE(bc), 2)   ! br N*N
           bc = 0.0d+00

           Call bBmat(e0, dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           if (rank == 0) then ! Process limits for output
               write (normaloutput, *) 'bc matrix is obtained normally'
           end if
           If (debug) then

               if (rank == 0) then ! Process limits for output
                   write (normaloutput, *) 'Check whether U*SU is diagonal'
               end if
               Call checkdgc(dimn, sc0, sc, ws)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               if (rank == 0) then  ! Process limits for output
                   write (normaloutput, *) 'Check whether U*SU is diagonal END'
               end if
           End if

           deallocate (sc0); Call memminus(KIND(sc0), SIZE(sc0), 2)

           if (rank == 0) then ! Process limits for output
               write (normaloutput, *) 'OK cdiag', dimn, dimm
           end if
           Allocate (uc(dimn, dimm)); Call memplus(KIND(uc), SIZE(uc), 2)           ! uc N*M
           Allocate (wsnew(dimm)); Call memplus(KIND(wsnew), SIZE(wsnew), 1)     ! wnew M
           uc(:, :) = 0.0d+00
           wsnew(:) = 0.0d+00

           Call ccutoff(sc, ws, dimn, dimm, uc, wsnew)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           if (rank == 0) then ! Process limits for output
               write (normaloutput, *) 'OK ccutoff'
           end if
           deallocate (sc); Call memminus(KIND(sc), SIZE(sc), 2)
           deallocate (ws); Call memminus(KIND(ws), SIZE(ws), 1)

           Call ucramda_s_half(uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           deallocate (wsnew); Call memminus(KIND(wsnew), SIZE(wsnew), 1)
           if (rank == 0) then ! Process limits for output
               write (normaloutput, *) 'ucrams half OK'
           end if
           Allocate (bc0(dimm, dimn)); Call memplus(KIND(bc0), SIZE(bc0), 2) ! bc0 M*N
           bc0 = 0.0d+00
           bc0 = MATMUL(TRANSPOSE(DCONJG(uc)), bc)

           Allocate (bc1(dimm, dimm)); Call memplus(KIND(bc1), SIZE(bc1), 2) ! bc1 M*M
           bc1 = 0.0d+00
           bc1 = MATMUL(bc0, uc)

           If (debug) then

               if (rank == 0) then ! Process limits for output
                   write (normaloutput, *) 'Check whether bc1 is hermite or not'
               end if
               Do i = 1, dimm
                   Do j = i, dimm
                       if (ABS(bc1(i, j) - DCONJG(bc1(j, i))) > 1.0d-6) then
                           if (rank == 0) then ! Process limits for output
                               write (normaloutput, '(2I4,2E15.7)') i, j, bc1(i, j) - bc1(j, i)
                           end if
                       End if
                   End do
               End do
               if (rank == 0) then ! Process limits for output
                   write (normaloutput, *) 'Check whether bc1 is hermite or not END'
               end if
           End if

           deallocate (bc); Call memminus(KIND(bc), SIZE(bc), 2)
           deallocate (bc0); Call memminus(KIND(bc0), SIZE(bc0), 2)

           cutoff = .FALSE.

           Allocate (wb(dimm)); Call memplus(KIND(wb), SIZE(wb), 1)

           if (rank == 0) then ! Process limits for output
               write (normaloutput, *) 'bC matrix is transrated to bc1(M*M matrix)!'
           end if
           Allocate (bc0(dimm, dimm)); Call memplus(KIND(bc0), SIZE(bc0), 2) ! bc0 M*M
           bc0 = bc1

           Call cdiag(bc1, dimm, dammy, wb, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           If (debug) then

               if (rank == 0) then ! Process limits for output
                   write (normaloutput, *) 'Check whether bc is really diagonalized or not'
               end if
               Call checkdgc(dimm, bc0, bc1, wb)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               if (rank == 0) then ! Process limits for output
                   write (normaloutput, *) 'Check whether bc is really diagonalized or not END'
               end if
           End if

           deallocate (bc0); Call memminus(KIND(bc0), SIZE(bc0), 2)

           if (rank == 0) then ! Process limits for output
               write (normaloutput, *) 'bC1 matrix is diagonalized!'
           end if
           e2 = 0.0d+00

           Do i0 = 1, nij
               ji = ii0(i0)
               jj = ij0(i0)

               syma = MULTB_D(irpmo(ji) - (-1)**(mod(irpmo(ji), 2)), irpmo(jj))
               syma = MULTB_S(syma, isym)

               If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then

                   Allocate (vc(dimn)); Call memplus(KIND(vc), SIZE(vc), 2)
                   Do it = 1, dimn
                       vc(it) = v(i0, indsym(1, it) + ninact, indsym(2, it) + ninact)
                   End do

                   Allocate (vc1(dimm)); Call memplus(KIND(vc1), SIZE(vc1), 2)
                   vc1 = 0.0d+00

                   vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn, 1:dimm))), vc(1:dimn)) ! v => v~
                   Deallocate (vc); Call memminus(KIND(vc), SIZE(vc), 2)

                   alpha = -eps(ji) - eps(jj) - e0 + eshift   ! For Level Shift (2007/2/9)

                   vc1(1:dimm) = &
                   & MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm, 1:dimm))), vc1(1:dimm)) ! v~ => v~~

                   Do j = 1, dimm
                       sumc2local = sumc2local + (ABS(vc1(j))**2.0d+00)/((alpha + wb(j))**2.0d+00)
                       e = (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                       e2(isym) = e2(isym) - e
                   End do

                   Deallocate (vc1); Call memminus(KIND(vc1), SIZE(vc1), 2)

               End if

           End do            !i0
           if (rank == 0) then ! Process limits for output
               write (normaloutput, '("e2b(",I3,") = ",E20.10,"a.u.")') isym, e2(isym)
           end if
           Deallocate (bc1); Call memminus(KIND(bc1), SIZE(bc1), 2)
           Deallocate (uc); Call memminus(KIND(uc), SIZE(uc), 2)
           Deallocate (wb); Call memminus(KIND(wb), SIZE(wb), 1)
           Deallocate (indsym); Call memminus(KIND(indsym), SIZE(indsym), 2)

1000       e2b = e2b + e2(isym)

       End do               ! isym

       if (rank == 0) then ! Process limits for output
           write (normaloutput, '("e2b      = ",E20.10,"a.u.")') e2b
           write (normaloutput, '("sumc2,b  = ",E20.10)') sumc2local
       end if
       sumc2 = sumc2 + sumc2local

       deallocate (iij); Call memminus(KIND(iij), SIZE(iij), 1)
       deallocate (ii0); Call memminus(KIND(ii0), SIZE(ii0), 1)
       deallocate (ij0); Call memminus(KIND(ij0), SIZE(ij0), 1)
       deallocate (v); Call memminus(KIND(v), SIZE(v), 2)

       continue
       if (rank == 0) then ! Process limits for output
           write (normaloutput, *) 'end solvB_ord_ty'
       end if
   end

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE sBmat(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space B

!  S(xy, tu) = <0|EtxEuy|0> -d(tx)<0|Euy|0> -d(uy)<0|Etx|0> -d(ty)<0|Eux|0> +d(tx)d(uy)-d(ty)d(ux)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE

       integer, intent(in)      :: dimn, indsym(2, dimn)
       complex*16, intent(out)  :: sc(dimn, dimn)

       real*8  :: a, b

       integer :: it, iu, iy, ix, ivx, itu
       integer :: jt, ju, jy, jx
       integer :: i, j
       integer :: count

       sc = 0.0d+00

       Do i = 1, dimn

           ix = indsym(1, i)
           iy = indsym(2, i)

           Do j = i, dimn

               it = indsym(1, j)
               iu = indsym(2, j)

!  S(xy, tu) = <0|EtxEuy|0> -d(tx)<0|Euy|0> -d(uy)<0|Etx|0> -d(ty)<0|Eux|0> +d(tx)d(uy)-d(ty)d(ux)
!                                                                                      ~~~~~~~~This term is0
               Call dim2_density(it, ix, iu, iy, a, b)
               sc(i, j) = sc(i, j) + DCMPLX(a, b)

               If (it == ix) then
                   Call dim1_density(iu, iy, a, b)
                   sc(i, j) = sc(i, j) - DCMPLX(a, b)
               End if

               If (iu == iy) then
                   Call dim1_density(it, ix, a, b)
                   sc(i, j) = sc(i, j) - DCMPLX(a, b)
               End if

               If (it == iy) then
                   Call dim1_density(iu, ix, a, b)
                   sc(i, j) = sc(i, j) - DCMPLX(a, b)
               End if

               If ((it == ix) .and. (iu == iy)) then
                   sc(i, j) = sc(i, j) + 1.0d+00
               End if

!              If((it == iy).and.(iu == ix)) then
!                 write(*,*)'it == iy).and.(iu == ix)'
!                 sc(i,j) = sc(i,j) - 1.0d+00
!              Endif

               sc(j, i) = DCONJG(sc(i, j))

           End do               !j
       End do                  !i

   End subroutine sBmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE bBmat(e0, dimn, sc, indsym, bc) ! Assume C1 molecule, overlap matrix B in space B
!
!  B(xy,tu) = Siguma_w [eps(w){<0|EtxEuyEww|0>-d(tx)<0|EuyEww|0> -d(uy)<0|EtxEww|0> -d(ty)<0|EuxEww|0>]
!
!              +{d(tx)d(uy)-d(ty)d(ux)}*e0      +S(xy,tu){eps(t)+eps(u)}
!
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE
       include 'mpif.h'

       integer, intent(in) :: dimn, indsym(2, dimn)
       complex*16, intent(in)  :: sc(dimn, dimn)
       complex*16, intent(out) :: bc(dimn, dimn)
       real*8, intent(in)      :: e0

       real*8              :: e, denr, deni
       complex*16          :: den

       integer :: it, iu, iv, ix, iy, iz, iw
       integer :: jt, ju, jy, jx, jw, i, j

       bc(:, :) = 0.0d+00

       if (rank == 0) then ! Process limits for output
           write (normaloutput, *) 'B space Bmat iroot=', iroot
       end if

       !$OMP parallel do private(ix,iy,jx,jy,it,iu,jt,ju,e,j,iw,jw,denr,deni,den)
       Do i = rank + 1, dimn, nprocs
           ! Do i = 1, dimn

           ix = indsym(1, i)
           jx = ix + ninact
           iy = indsym(2, i)
           jy = iy + ninact

           Do j = i, dimn

               it = indsym(1, j)
               jt = it + ninact
               iu = indsym(2, j)
               ju = iu + ninact

!  B(xy,tu) = Siguma_w [eps(w){<0|EtxEuyEww|0>-d(tx)<0|EuyEww|0> -d(uy)<0|EtxEww|0> -d(ty)<0|EuxEww|0>]
!
!              +{d(tx)d(uy)-d(ty)d(ux)}*e0      +S(xy,tu){eps(t)+eps(u)}

               e = eps(jt) + eps(ju)

               Do iw = 1, nact
                   jw = iw + ninact

                   Call dim3_density &
                       (it, ix, iu, iy, iw, iw, denr, deni)
                   den = DCMPLX(denr, deni)
                   bc(i, j) = bc(i, j) + den*eps(jw)

                   If (it == ix) then

                       Call dim2_density(iu, iy, iw, iw, denr, deni)
                       den = DCMPLX(denr, deni)
                       bc(i, j) = bc(i, j) - den*eps(jw)

                   End if

                   If (iu == iy) then

                       Call dim2_density(it, ix, iw, iw, denr, deni)
                       den = DCMPLX(denr, deni)
                       bc(i, j) = bc(i, j) - den*eps(jw)

                   End if

                   If (it == iy) then

                       Call dim2_density(iu, ix, iw, iw, denr, deni)
                       den = DCMPLX(denr, deni)
                       bc(i, j) = bc(i, j) - den*eps(jw)

                   End if

               End do

!              +{d(tx)d(uy)-d(ty)d(ux)}*e0      +S(xy,tu){eps(t)+eps(u)}

               If ((it == ix) .and. (iu == iy)) then
                   bc(i, j) = bc(i, j) + e0
               End if

!              If((it == iy) .and.(iu == ix)) then        ! THIS TERM IS 0
!                 bc(i, j) = bc(i, j) - e0
!              Endif

               bc(i, j) = bc(i, j) + sc(i, j)*e

               bc(j, i) = DCONJG(bc(i, j))

           End do               !i
       End do                  !j
       if (rank == 0) then
        call MPI_Reduce (MPI_IN_PLACE, bc(1,1), dimn**2, MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       else
        call MPI_Reduce (bc(1,1), bc(1,1), dimn**2, MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       end if
       if (rank == 0) then ! Process limits for output
           write (normaloutput, *) 'bBmat is ended'
       end if
   End subroutine bBmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE vBmat_ord_ty(nij, iij, v)
!
!
! V(i,j) =  SIGUMA_p,q:active <0|EptEqu|0>(pi|qj)
!
!           + SIGUMA_p:active[<0|Ept|0> {(pj|ui) - (uj|pi)}  - <0|Epu|0> (ti|pj)]
!
!           + (uj|ti)  - (tj|ui)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE
       include 'mpif.h'

       integer, intent(in)     :: nij, iij(ninact, ninact)
       complex*16, intent(out) :: v(nij, ninact + 1:ninact + nact, ninact + 1:ninact + nact)

       real*8                  :: dr, di
       complex*16              :: cint2, dens

       integer :: i, j, k, l, tij, ip, iq, save, count
       integer :: it, jt, ju, iu

       v = 0.0d+00

       !   open(1, file ='Bint', status='old', form='unformatted')  !  (21|21) stored (ti|uj) i > j
       open (1, file=bint, status='old', form='formatted')  !  (21|21) stored (ti|uj) i > j

30     read (1, '(4I4, 2e20.10)', err=10, end=20) i, j, k, l, cint2                    !  (ij|kl)

       if (j <= l) goto 30
!        write(*,'(4I4,2E20.10)')i,j,k,l,cint2

!------------------------------------------------------------------------------------------------
!  i > j
!
! V(i,j) =  SIGUMA_p,q:active <0|EptEqu|0>(pi|qj)                                      ! term1
!
!           + SIGUMA_p:active[<0|Ept|0> {(ui|pj) - (pi|uj)}  - <0|Epu|0> (ti|pj)]      ! term2
!
!           + (ti|uj)  - (ui|tj)                                                       ! term3
!
!------------------------------------------------------------------------------------------------

       tij = iij(j, l)

!        write(*,'(5I4,2E20.10)')i,j,k,l,tij,cint2

       ! Term 3 !        + (ti|uj)  - (ui|tj)  (i > j)

       v(tij, i, k) = v(tij, i, k) + cint2 !  + (ti|uj)
       v(tij, k, i) = v(tij, k, i) - cint2 !  - (ui|tj)

       ! Term 2 !  + SIGUMA_p:active[<0|Ept|0> {(ui|pj) - (pi|uj)}  - <0|Epu|0> (ti|pj)]
       !                             ===========================      ================
       !                                loop for t                     loop for u(variable u is renamed to t)
       Do it = 1, nact
           jt = it + ninact

           Call dim1_density(k - ninact, it, dr, di)
           dens = DCMPLX(dr, di)
           v(tij, jt, i) = v(tij, jt, i) + cint2*dens
           v(tij, i, jt) = v(tij, i, jt) - cint2*dens

           Call dim1_density(i - ninact, it, dr, di)
           dens = DCMPLX(dr, di)
           v(tij, jt, k) = v(tij, jt, k) - cint2*dens

           ! Term1 !   SIGUMA_p,q:active <0|EptEqu|0>(pi|qj)                                      ! term1
           !                             ==================
           !                              loop for t and u

           Do iu = 1, it - 1
               ju = iu + ninact
               Call dim2_density(i - ninact, it, k - ninact, iu, dr, di)
               dens = DCMPLX(dr, di)
               v(tij, jt, ju) = v(tij, jt, ju) + cint2*dens
           End do

       End do

       goto 30

20     close (1)
       if (rank == 0) then ! Process limits for output
           write (normaloutput, *) 'reading int2 is over'
       end if
       goto 100

10     write (*, *) 'error while opening file Bint'; goto 100

100    if (rank == 0) write (normaloutput, *) 'vBmat_ord_ty is ended'

       !  v(nij, ninact + 1:ninact + nact, ninact + 1:ninact + nact)
       call MPI_Allreduce(MPI_IN_PLACE, v(1, ninact + 1, ninact + 1), nij*nact**2, &
                          MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
   end subroutine vBmat_ord_ty
