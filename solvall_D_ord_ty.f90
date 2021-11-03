! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE solvD_ord_ty(e0, e2d)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE

       real*8, intent(in) :: e0
       real*8, intent(out):: e2d

       integer :: dimn, dimm, count, dammy

       integer, allocatable :: indsym(:, :)

       real*8, allocatable  :: sr(:, :), ur(:, :)
       real*8, allocatable  :: br(:, :), wsnew(:), ws(:), wb(:)
       real*8, allocatable  :: br0(:, :), br1(:, :)
       real*8               :: e2(nsymrpa*2), e, alpha

       complex*16, allocatable  :: sc(:, :), uc(:, :), sc0(:, :)
       complex*16, allocatable  :: bc(:, :)
       complex*16, allocatable  :: bc0(:, :), bc1(:, :), v(:, :, :), vc(:), vc1(:)

       integer, allocatable     :: ia0(:), ii0(:), iai(:, :)
       integer                  :: nai

       logical :: cutoff
       integer :: j, i, k, syma, isym, i0, j0
       integer :: ia, it, ii, iu
       integer :: ja, jt, ji, ju

       real*8  :: thresd

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!      SPACE D IS NOW CALCULATED
!
!     EaiEtu|0>
!
!   DRAS1 = -1   DRAS2 = 0   DRAS3 = +1
!
!   t and u run all active spinor space independently!
!
!
!  S(bjxy,aitu) = d(ba)d(ij)<0|EyxEtu|0>
!
!  S(xy, tu) = <0|EyxEtu|0>
!
!  B(xy,tu) = Siguma_w [eps(w)<0|EyxEtuEww|0>+S(xy,tu)(eps(t)-eps(u))]
!
!  a(a,i)       = eps(a) - eps(i) - e0
!
!! dame!! V(a,i) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!! dame!!
!! dame!! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) + SIGUMA_p:active(ap|pi)}]
!
!
! V(a,i) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki) }
!
!
!  E2 = SIGUMA_a,i, dimm |V1(dimm,ai)|^2|/{(alpha(ai) + wb(dimm)}

!        thresd = thres
       thresd = 1.0D-08
       thres = 1.0D-08

       e2 = 0.0d+00
       e2d = 0.0d+00
       dimn = 0
       syma = 0

       write (*, *) ' ENTER solv D part'
       write (*, *) ' nsymrpa', nsymrpa

       i0 = 0
       Do ia = 1, nsec
           Do ii = 1, ninact
               i0 = i0 + 1
           End do
       End do

       nai = i0
       Allocate (iai(ninact + nact + 1:ninact + nact + nsec, ninact))
       iai = 0
       Allocate (ia0(nai))
       Allocate (ii0(nai))

       i0 = 0
       Do ia = 1, nsec
           ja = ia + ninact + nact
           Do ii = 1, ninact
               i0 = i0 + 1
               iai(ja, ii) = i0
               ia0(i0) = ja
               ii0(i0) = ii
           End do
       End do
       Allocate (v(nai, ninact + 1:ninact + nact, ninact + 1:ninact + nact))
       v = 0.0d+00

       Call vDmat_ord_ty(nai, iai, v)

       write (*, *) 'come'

       Do isym = 1, nsymrpa

           dimn = 0
           Do it = 1, nact
               jt = it + ninact
               Do iu = 1, nact
                   ju = iu + ninact

                   syma = MULTB_D(irpmo(jt), irpmo(ju))

                   if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                       dimn = dimn + 1
                   End if
100            End do               ! iu
           End do                  ! it

!        write(*,*)'isym, dimn',isym, dimn

           If (dimn == 0) goto 1000

           Allocate (indsym(2, dimn))
           indsym = 0
           dimn = 0

           Do it = 1, nact
               jt = it + ninact
               Do iu = 1, nact
                   ju = iu + ninact

                   syma = MULTB_D(irpmo(jt), irpmo(ju))

                   if (nsymrpa == 1 .or. (nsymrpa /= 1 .and. syma == isym)) then
                       dimn = dimn + 1
                       indsym(1, dimn) = it
                       indsym(2, dimn) = iu
                   End if
200            End do               ! iu
           End do                  ! it

           Allocate (sc(dimn, dimn))
           sc = 0.0d+00            ! sc N*N

           Call sDmat(dimn, indsym, sc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           write (*, *) 'sc matrix is obtained normally'
           Allocate (ws(dimn))
           ws = 0.0d+00
           cutoff = .TRUE.
           thresd = 1.0d-08
!           thresd = 1.0d-15

           Allocate (sc0(dimn, dimn))
           sc0 = sc

           Call cdiag(sc, dimn, dimm, ws, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           write (*, *) 'after s cdiag'
!           Do i0 = 1, dimn
!           write(*,'(E20.10)') ws(i0)
!           End do

           If (dimm == 0) then
               deallocate (indsym)
               deallocate (sc0)
               deallocate (sc)
               deallocate (ws)
               goto 1000
           End if

           If (debug) then

               write (*, *) 'Check whether U*SU is diagonal'

               Call checkdgc(dimn, sc0, sc, ws)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               write (*, *) 'Check whether U*SU is diagonal END'

           End if

           Allocate (bc(dimn, dimn))                                 ! bc N*N
           bc = 0.0d+00

           Call bDmat(dimn, sc0, indsym, bc)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           write (*, *) 'bc matrix is obtained normally'

           deallocate (sc0)

           write (*, *) 'OK cdiag', dimn, dimm

           Allocate (uc(dimn, dimm))                                 ! uc N*M
           Allocate (wsnew(dimm))                                  ! wnew M
           uc(:, :) = 0.0d+00
           wsnew(:) = 0.0d+00

           Call ccutoff(sc, ws, dimn, dimm, uc, wsnew)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           write (*, *) 'OK ccutoff'
           deallocate (ws)
           deallocate (sc)

           Call ucramda_s_half(uc, wsnew, dimn, dimm)    ! uc N*M matrix rewritten as uramda^(-1/2)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           deallocate (wsnew)

           write (*, *) 'ucrams half OK'
           Allocate (bc0(dimm, dimn))                       ! bc0 M*N
           bc0 = 0.0d+00
           bc0 = MATMUL(TRANSPOSE(DCONJG(uc)), bc)
           Allocate (bc1(dimm, dimm))                      ! bc1 M*M
           bc1 = 0.0d+00
           bc1 = MATMUL(bc0, uc)

           IF (debug) then

               write (*, *) 'Check whether bc1 is hermite or not'
               Do i = 1, dimm
                   Do j = i, dimm
                       if (ABS(bc1(i, j) - DCONJG(bc1(j, i))) > 1.0d-6) then
                           write (*, '(2I4,2E15.7)') i, j, bc1(i, j) - bc1(j, i)
                       End if
                   End do
               End do
               write (*, *) 'Check whether bc1 is hermite or not END'

           End if

           deallocate (bc)
           deallocate (bc0)

           cutoff = .FALSE.

           Allocate (wb(dimm))
           wb = 0.0d+00

           write (*, *) 'bC matrix is transrated to bc1(M*M matrix)!'

           Allocate (bc0(dimm, dimm))
           bc0 = bc1

           Call cdiag(bc1, dimm, dammy, wb, thresd, cutoff)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           If (debug) then

               write (*, *) 'Check whether bc is really diagonalized or not'

               Call checkdgc(dimm, bc0, bc1, wb)
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               write (*, *) 'Check whether bc is really diagonalized or not END'
           End if

           deallocate (bc0)

           write (*, *) 'bC1 matrix is diagonalized!'

           e2 = 0.0d+00
           Do i0 = 1, nai
               ja = ia0(i0)
               ji = ii0(i0)

               syma = MULTB_D(irpmo(ja), irpmo(ji))
               syma = MULTB_S(syma, isym)

               If (nsymrpa == 1 .or. (nsymrpa /= 1 .and. (syma == 1))) then

                   Allocate (vc(dimn))
                   Do it = 1, dimn
                       vc(it) = v(i0, indsym(1, it) + ninact, indsym(2, it) + ninact)
                   End do

                   Allocate (vc1(dimm))

                   vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(uc(1:dimn, 1:dimm))), vc(1:dimn))
                   Deallocate (vc)

                   alpha = +eps(ja) - eps(ji) - e0 + eshift  ! For Level Shift (2007/2/9)

                   vc1(1:dimm) = MATMUL(TRANSPOSE(DCONJG(bc1(1:dimm, 1:dimm))), vc1(1:dimm))

                   Do j = 1, dimm
                       sumc2local = sumc2local + (ABS(vc1(j))**2.0d+00)/((alpha + wb(j))**2.0d+00)
                       e = (ABS(vc1(j))**2.0d+00)/(alpha + wb(j))
                       e2(isym) = e2(isym) - e
                   End do

                   deallocate (vc1)

               End if

           End do               !i0

           deallocate (indsym)
           deallocate (uc)
           deallocate (wb)
           Deallocate (bc1)

1000       write (*, '("e2d(",I3,") = ",E20.10,"a.u.",I4)') isym, e2(isym), rank
           e2d = e2d + e2(isym)

       End do                  ! isym

       write (*, '("e2d      = ",E20.10,"a.u.",I4)') e2d, rank

       write (*, '("sumc2,d  = ",E20.10)') sumc2local
       sumc2 = sumc2 + sumc2local

       deallocate (iai)
       deallocate (ia0)
       deallocate (ii0)
       deallocate (v)

       continue
       write (*, *) 'end solvD_ord'
   end

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE sDmat(dimn, indsym, sc) ! Assume C1 molecule, overlap matrix S in space D

!  S(xy, tu) = <0|EyxEtu|0>
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE

       integer, intent(in)      :: dimn, indsym(2, dimn)
       complex*16, intent(out)  :: sc(dimn, dimn)

       real*8  :: a, b

       integer :: it, iu, iy, ix, ivx, itu
       integer :: i, j
       integer :: count

       sc = 0.0d+00

       Do i = 1, dimn

           ix = indsym(1, i)
           iy = indsym(2, i)
           Do j = i, dimn

               it = indsym(1, j)
               iu = indsym(2, j)

               a = 0.0d+0
               b = 0.0d+0
               Call dim2_density(iy, ix, it, iu, a, b)
               sc(i, j) = DCMPLX(a, b)
               sc(j, i) = DCONJG(sc(i, j))

!              If(ABS(sc(i,j)) > 1.0d+00) then
!                 write(*,'(2I4,2E20.10)')i,j,sc(i,j)
!              Endif

           End do               !j
       End do                  !i

   End subroutine sDmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE bDmat(dimn, sc, indsym, bc) ! Assume C1 molecule, overlap matrix B in space D
!
!  B(xy,tu) = Siguma_w [eps(w)<0|EyxEtuEww|0>+S(xy,tu)(eps(t)-eps(u))]
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE

       integer, intent(in) :: dimn, indsym(2, dimn)
       complex*16, intent(in)  :: sc(dimn, dimn)
       complex*16, intent(out) :: bc(dimn, dimn)

       real*8              :: e, denr, deni
       complex*16          :: den

       integer :: it, iu, iv, ix, iy, iz, iw
       integer :: jt, ju, jy, jx, jw, i, j

       bc(:, :) = 0.0d+00

       write (*, *) 'F space Bmat iroot=', iroot

       Do i = 1, dimn

           ix = indsym(1, i)
           jx = ix + ninact
           iy = indsym(2, i)
           jy = iy + ninact

           Do j = i, dimn

               it = indsym(1, j)
               jt = it + ninact
               iu = indsym(2, j)
               ju = iu + ninact

!  B(xy,tu) = Siguma_w [eps(w)<0|EyxEtuEww|0>+S(xy,tu)(eps(t)-eps(u))]

               e = eps(jt) - eps(ju)

               Do iw = 1, nact
                   jw = iw + ninact

                   Call dim3_density &
                       (iy, ix, it, iu, iw, iw, denr, deni)
                   den = DCMPLX(denr, deni)
                   bc(i, j) = bc(i, j) + den*eps(jw)

               End do

               bc(i, j) = bc(i, j) + sc(i, j)*e

               bc(j, i) = DCONJG(bc(i, j))

           End do               !i
       End do                  !j

       write (*, *) 'bDmat is ended'

   End subroutine bDmat

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!
   SUBROUTINE vDmat_ord_ty(nai, iai, v)
!
!
! V(a,i) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE
       include 'mpif.h'
       integer, intent(in)     :: nai, iai(ninact + nact + 1:ninact + nact + nsec, ninact)

       complex*16, intent(out) :: v(nai, ninact + 1:ninact + nact, ninact + 1:ninact + nact)

       real*8                  :: dr, di, signkl
       complex*16              :: cint1, cint2, dens, d
       complex*16              :: effh(ninact + nact + 1:ninact + nact + nsec, ninact)

       integer :: i, j, k, l, tai, ip, iq, save, count
       integer :: it, jt, ju, iu, ia, ii, ja, ji, kkr, lkr
       logical :: test

       v = 0.0d+00
       effh = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! V(tai, jt, ju) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
       if (rank == 0) then

           Do ia = 1, nsec
               ja = ia + ninact + nact
               Do ii = 1, ninact
                   ji = ii
                   Call tramo1_ty(ja, ji, cint1)
                   effh(ja, ji) = cint1
!              if(ja==19.and.ji==1) write(*,'("effh int1 ",2I4,2E20.10)')ja,ji,cint1
               End do
           End do
       end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! V(a,i, jt, ju) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Three types of integrals are stored
!
!  (31|22) stored (ai|pq) ...TYPE 1 D1int
!  (32|21) stored (ap|qi) ...TYPE 2 D2int
!
!  (31|11) stored (ai|jk) ...TYPE 3 D3int
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       !  open (1, file='D1int', status='old', form='unformatted')
       open (1, file=d1int, status='old', form='formatted')

30     read (1, '(4I4, 2e20.10)', err=10, end=20) i, j, k, l, cint2 !  (ij|kl)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! V(a,i, jt, ju) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       ja = i
       ji = j
       tai = iai(ja, ji)
!           write(*,'("type1 (31|22)",4I4,2E20.10)')i,j,k,l,cint2

       Do it = 1, nact
           jt = it + ninact
           Do iu = 1, nact
               ju = iu + ninact

               Call dim2_density(iu, it, k - ninact, l - ninact, dr, di)
               d = DCMPLX(dr, di)
               v(tai, jt, ju) = v(tai, jt, ju) + d*cint2

           End do
       End do

       goto 30
20     close (1)
       write (*, *) 'reading D2int2 is over'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! V(a,i, jt, ju) = SIGUMA_pq:active <0|EutEpq|0>{(ai|pq) - (aq|pi)}
!
! + <0|Eut|0>[hai +{ SIGUMA_k:inactive(ai|kk) - (ak|ki)}]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       !  open (1, file='D2int', status='old', form='unformatted')
       open (1, file=d2int, status='old', form='formatted')

31     read (1, '(4I4, 2e20.10)', err=10, end=21) i, j, k, l, cint2 !  (ij|kl)
       ja = i
       ji = l
       tai = iai(ja, ji)

       Do it = 1, nact
           jt = it + ninact
           Do iu = 1, nact
               ju = iu + ninact

               Call dim2_density(iu, it, k - ninact, j - ninact, dr, di)
               d = DCMPLX(dr, di)

               v(tai, jt, ju) = v(tai, jt, ju) - d*cint2

           End do
       End do

       goto 31

21     close (1)
       write (*, *) 'reading D2int2 is over'

       !  open (1, file='D3int', status='old', form='unformatted') ! (ai|jk) is stored
       open (1, file=d3int, status='old', form='formatted') ! (ai|jk) is stored

300    read (1, '(4I4, 2e20.10)', err=10, end=200) i, j, k, l, cint2 !  (ij|kl)
!        write(*,*)'D1int', i,j,k,l ,cint2

       if (j /= k .and. k == l) then !(ai|kk)

           effh(i, j) = effh(i, j) + cint2

       elseif (j == k .and. k /= l) then !(ak|ki)

           effh(i, l) = effh(i, l) - cint2

       end if

       goto 300

200    close (1)
       write (*, *) 'reading D3int2 is over'
    !    effh(ninact + nact + 1:ninact + nact + nsec, ninact)
       if (rank == 0) then
           call MPI_Reduce(MPI_IN_PLACE, effh(ninact + nact + 1, 1), nsec*ninact, &
                           MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       else
           call MPI_Reduce(effh(ninact + nact + 1, 1), effh(ninact + nact + 1, 1), nsec*ninact, &
                           MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       end if
       if (rank /= 0) then
        effh(:,:) = 0
       end if

       Do ia = 1, nsec
           ja = ia + ninact + nact
           Do ii = 1, ninact
               ji = ii
               tai = iai(ja, ji)
!              if(ABS(effh(ja,ji)) > 1.0d-10) write(*,'("effh ",2I4,2E20.10)')ja,ji,effh(ja,ji)

               Do it = 1, nact
                   jt = it + ninact
                   Do iu = 1, nact
                       ju = iu + ninact

                       Call dim1_density(iu, it, dr, di)

                       d = DCMPLX(dr, di)
                       v(tai, jt, ju) = v(tai, jt, ju) + effh(ja, ji)*d
                   End do
               End do

           End do
       End do

       goto 100

10     write (*, *) 'error while opening file Dint'; goto 100

100    write (*, *) 'vDmat_ord_ty is ended'
       !  v(nai, ninact + 1:ninact + nact, ninact + 1:ninact + nact)
       call MPI_Allreduce(MPI_IN_PLACE, v(1, ninact + 1, ninact + 1), nai*nact**2, &
                          MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierr)
   end subroutine vDmat_ord_ty
