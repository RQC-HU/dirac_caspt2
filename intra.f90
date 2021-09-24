! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE intra_1(spi, spj, spk, spl, fname)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE

       integer, intent(in)        :: spi, spj, spk, spl
       character*5, intent(in)    :: fname

       integer, allocatable :: i(:), j(:), k(:), l(:), indsym(:, :, :), nsym(:, :)
       complex*16, allocatable :: cint2(:), traint2(:, :, :, :)

       real*8                  :: thresd
       complex*16              :: cint2_0

       integer :: i0, j0, k0, l0, i1, j1, k1, l1, inew, jnew, knew, lnew
       integer :: ii, ji, ki, li, ie, je, ke, le
       integer :: inz, nmx, ini(3), end(3), isp, isym, imo, nmaxint, ired

       logical :: is_opened
       thresd = 1.0d-15

       ini(1) = 1
       end(1) = ninact
       ini(2) = ninact + 1
       end(2) = ninact + nact
       ini(3) = ninact + nact + 1
       end(3) = ninact + nact + nsec

       nmx = max(ninact, nact, nsec)
       Allocate (indsym(3, nsymrpa, nmx)); Call memplus(KIND(indsym), SIZE(indsym), 1)
       Allocate (nsym(3, nsymrpa)); Call memplus(KIND(nsym), SIZE(nsym), 1)
       indsym = 0
       nsym = 0

       Do isym = 1, nsymrpa
           Do isp = 1, 3                           ! i th space( inact or act or sec)
               Do imo = ini(isp), end(isp)
                   if (irpamo(imo) == isym) then
                       nsym(isp, isym) = nsym(isp, isym) + 1
                       indsym(isp, isym, nsym(isp, isym)) = imo
                   end if
               End do
           End do
       End do

       ii = ini(spi)
       ie = end(spi)
       ji = ini(spj)
       je = end(spj)
       ki = ini(spk)
       ke = end(spk)
       li = ini(spl)
       le = end(spl)

       Allocate (traint2(ii:ie, ji:je, ki:ke, li:le)); Call memplus(KIND(traint2), SIZE(traint2), 2)

!debug
!        write(*,'("C1int",8I4)')ii,ie,ji,je,ki,ke,li,le

       traint2 = 0.0d+00
       write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024

!        nmaxint = AINT(3.5d+09 - tmem)/48 ! one integrals required four integer and one complex values

! Abe modified 2016. 11.11  100 GB is the max memory
       nmaxint = AINT(5.00d+10 - tmem)/48 ! one integrals required four integer and one complex values
! Abe modified 2016. 11.11

       write (*, *) 'nmaxint = ', nmaxint
       write (*, *) 'tmem, nmaxint*48', tmem, nmaxint*48
       IF (nmaxint < 0) stop

       Allocate (i(nmaxint)); Call memplus(KIND(i), SIZE(i), 1)
       Allocate (j(nmaxint)); Call memplus(KIND(j), SIZE(j), 1)
       Allocate (k(nmaxint)); Call memplus(KIND(k), SIZE(k), 1)
       Allocate (l(nmaxint)); Call memplus(KIND(l), SIZE(l), 1)

       Allocate (cint2(nmaxint)); Call memplus(KIND(cint2), SIZE(cint2), 2)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and first index transformation !
!                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')  ! no symmetry about spi,spj,spk,spl

30     Do inz = 1, nmaxint
           read (1, '(4I4,2e20.10)', err=10, end=20) i(inz), j(inz), k(inz), l(inz), cint2(inz)
       End do

       Do inz = 1, nmaxint
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(l0)

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

           Call takekr(i0, j0, k0, l0, cint2_0)
!           write(*,'("takekr")')i0,j0,k0,l0,cint2_0
           isym = irpmo(l0)

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

       End do

       goto 30

20     ired = inz
!        write(*,*)'ired',ired,i(ired),j(ired),k(ired),l(ired)
       Do inz = 1, ired - 1
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpamo(l0)
!           write(*,*)'nsym(spl,isym)',spl,isym,nsym(spl,isym)

!debug iwamuro
!           write(*,)i0,j0,k0,l0
           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

           Call takekr(i0, j0, k0, l0, cint2_0)
!           write(*,*)i0,j0,k0,l0
           isym = irpamo(l0)

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

       End do

       close (1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

       Do l0 = li, le
           Do k0 = ki, ke
               Do j0 = ji, je
                   Do i0 = ii, ie
                       if (ABS(traint2(i0, j0, k0, l0)) > thresd) then
                           write (1, '(4I4, 2e20.10)') i0, j0, k0, l0, traint2(i0, j0, k0, l0)
                       end if
                   End do
               End do
           End do
       End do

       close (1)

       traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and second index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

31     Do inz = 1, nmaxint
           read (1, '(4I4,2e20.10)', err=10, end=21) i(inz), j(inz), k(inz), l(inz), cint2(inz)
!           write(*,*) i(inz),j(inz),k(inz),l(inz),cint2(inz)
       End do

       Do inz = 1, nmaxint
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(k0)

           Do knew = 1, nsym(spk, isym)
               k1 = indsym(spk, isym, knew)
               traint2(i0, j0, k1, l0) = traint2(i0, j0, k1, l0) + cint2_0*DCONJG(f(k0, k1))
           End do

       End do

       goto 31

21     ired = inz
!        write(*,*)'ired',ired
       Do inz = 1, ired - 1
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(k0)

           Do knew = 1, nsym(spk, isym)
               k1 = indsym(spk, isym, knew)
               traint2(i0, j0, k1, l0) = traint2(i0, j0, k1, l0) + cint2_0*DCONJG(f(k0, k1))
           End do

       End do

       close (1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

       Do l0 = li, le
           Do k0 = ki, ke
               Do j0 = ji, je
                   Do i0 = ii, ie
                       if (ABS(traint2(i0, j0, k0, l0)) > thresd) then
                           write (1, '(4I4, 2e20.10)') i0, j0, k0, l0, traint2(i0, j0, k0, l0)
                       end if
                   End do
               End do
           End do
       End do

       close (1)

       traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and third index transformation  !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

32     Do inz = 1, nmaxint
           read (1, '(4I4,2e20.10)', err=10, end=22) i(inz), j(inz), k(inz), l(inz), cint2(inz)
       End do

       Do inz = 1, nmaxint
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(j0)

           Do jnew = 1, nsym(spj, isym)
               j1 = indsym(spj, isym, jnew)
               traint2(i0, j1, k0, l0) = traint2(i0, j1, k0, l0) + cint2_0*f(j0, j1)
           End do

       End do

       goto 32

22     ired = inz
!        write(*,*)'ired',ired
       Do inz = 1, ired - 1
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(j0)

           Do jnew = 1, nsym(spj, isym)
               j1 = indsym(spj, isym, jnew)
               traint2(i0, j1, k0, l0) = traint2(i0, j1, k0, l0) + cint2_0*f(j0, j1)
           End do

       End do

       close (1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

       Do l0 = li, le
           Do k0 = ki, ke
               Do j0 = ji, je
                   Do i0 = ii, ie
                       if (ABS(traint2(i0, j0, k0, l0)) > thresd) then
                           write (1, '(4I4, 2e20.10)') i0, j0, k0, l0, traint2(i0, j0, k0, l0)
                       end if
                   End do
               End do
           End do
       End do

       close (1)

       traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and fourth index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

33     Do inz = 1, nmaxint
           read (1, '(4I4,2e20.10)', err=10, end=23) i(inz), j(inz), k(inz), l(inz), cint2(inz)
       End do

       Do inz = 1, nmaxint
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(i0)

           Do inew = 1, nsym(spi, isym)
               i1 = indsym(spi, isym, inew)
               traint2(i1, j0, k0, l0) = traint2(i1, j0, k0, l0) + cint2_0*DCONJG(f(i0, i1))
           End do

       End do

       goto 33

23     ired = inz
!        write(*,*)'ired',ired
       Do inz = 1, ired - 1
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(i0)

           Do inew = 1, nsym(spi, isym)
               i1 = indsym(spi, isym, inew)
               traint2(i1, j0, k0, l0) = traint2(i1, j0, k0, l0) + cint2_0*DCONJG(f(i0, i1))
           End do

       End do

       close (1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

       Do l0 = li, le
           Do k0 = ki, ke
               Do j0 = ji, je
                   Do i0 = ii, ie
                       if (ABS(traint2(i0, j0, k0, l0)) > thresd) then
                           write (1, '(4I4, 2e20.10)') i0, j0, k0, l0, traint2(i0, j0, k0, l0)
                       end if
                   End do
               End do
           End do
       End do

       close (1)

       goto 100
!  10     write(*,*)'error opening file' ; goto 1000
10     write (*, *) 'error opening file'
       inquire (1, opened=is_opened)
       if (is_opened .eqv. .true.) then
           close (1)
       end if

100    deallocate (i); Call memminus(KIND(i), SIZE(i), 1)
       deallocate (j); Call memminus(KIND(j), SIZE(j), 1)
       deallocate (k); Call memminus(KIND(k), SIZE(k), 1)
       deallocate (l); Call memminus(KIND(l), SIZE(l), 1)

       deallocate (cint2); Call memminus(KIND(cint2), SIZE(cint2), 2)
       deallocate (traint2); Call memminus(KIND(traint2), SIZE(traint2), 2)

       deallocate (indsym); Call memminus(KIND(indsym), SIZE(indsym), 1)
       deallocate (nsym); Call memminus(KIND(nsym), SIZE(nsym), 1)

1000 end subroutine intra_1

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE intra_2(spi, spj, spk, spl, fname)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE

       integer, intent(in)        :: spi, spj, spk, spl
       character*5, intent(in)    :: fname

       integer, allocatable :: i(:), j(:), k(:), l(:), indsym(:, :, :), nsym(:, :)
       complex*16, allocatable :: cint2(:), traint2(:, :, :, :)

       real*8                  :: thresd
       complex*16              :: cint2_0

       integer :: i0, j0, k0, l0, i1, j1, k1, l1, inew, jnew, knew, lnew
       integer :: ii, ji, ki, li, ie, je, ke, le, lkr0, kkr0
       integer :: inz, nmx, ini(3), end(3), isp, isym, imo, nmaxint, ired, save

       logical :: is_opened
       thresd = 1.0d-15

       if (.not. (spi == spk .and. spj == spl)) then
           write (*, *) 'error intra_2', spi, spj, spk, spl
           stop
       end if

       ini(1) = 1
       end(1) = ninact
       ini(2) = ninact + 1
       end(2) = ninact + nact
       ini(3) = ninact + nact + 1
       end(3) = ninact + nact + nsec

       nmx = max(ninact, nact, nsec)
       Allocate (indsym(3, nsymrpa, nmx)); Call memplus(KIND(indsym), SIZE(indsym), 1)
       Allocate (nsym(3, nsymrpa)); Call memplus(KIND(nsym), SIZE(nsym), 1)
       indsym = 0
       nsym = 0

       Do isym = 1, nsymrpa
           Do isp = 1, 3                           ! i th space( inact or act or sec)
               Do imo = ini(isp), end(isp)
                   if (irpamo(imo) == isym) then
                       nsym(isp, isym) = nsym(isp, isym) + 1
                       indsym(isp, isym, nsym(isp, isym)) = imo
                   end if
               End do
           End do
       End do

       ii = ini(spi)
       ie = end(spi)
       ji = ini(spj)
       je = end(spj)
       ki = ini(spk)
       ke = end(spk)
       li = ini(spl)
       le = end(spl)

       Allocate (traint2(ii:ie, ji:je, ki:ke, li:le)); Call memplus(KIND(traint2), SIZE(traint2), 2)

!        write(*,*)ii,ie,ji,je,ki,ke,li,le

       traint2 = 0.0d+00
       write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024

!        nmaxint = AINT(3.5d+09 - tmem)/48 ! one integrals required four integer and one complex values

! Abe modified 2016. 11.11  100 GB is the max memory
       nmaxint = AINT(5.00d+10 - tmem)/48 ! one integrals required four integer and one complex values
! Abe modified 2016. 11.11

       write (*, *) 'nmaxint = ', nmaxint
       write (*, *) 'tmem, nmaxint*48', tmem, nmaxint*48
       IF (nmaxint < 0) stop

       write (*, *) 'nmaxint = ', nmaxint

       Allocate (i(nmaxint)); Call memplus(KIND(i), SIZE(i), 1)
       Allocate (j(nmaxint)); Call memplus(KIND(j), SIZE(j), 1)
       Allocate (k(nmaxint)); Call memplus(KIND(k), SIZE(k), 1)
       Allocate (l(nmaxint)); Call memplus(KIND(l), SIZE(l), 1)

       Allocate (cint2(nmaxint)); Call memplus(KIND(cint2), SIZE(cint2), 2)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and first index transformation !
!                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

30     Do inz = 1, nmaxint
           read (1, '(4I4,2e20.10)', err=10, end=20) i(inz), j(inz), k(inz), l(inz), cint2(inz)
!Iwamuro modify
!           write(*,'("cint2debugB",4I4,2E20.10)')i(inz),j(inz),k(inz),l(inz),cint2(inz)
       End do

       Do inz = 1, nmaxint
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(l0)

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

           if (i0 == k0 .and. j0 == l0) goto 50
           if (ABS(i0 - k0) == 1 .and. ABS(j0 - l0) == 1 .and. &
           & ABS(i0/2 - k0/2) == 1 .and. ABS(j0/2 - l0/2) == 1) goto 50

           i0 = k(inz)
           j0 = l(inz)
           k0 = i(inz)
           l0 = j(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(l0)

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

50         Call takekr(i0, j0, k0, l0, cint2_0)
           isym = irpmo(l0)

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

           if (i0 == k0 .and. j0 == l0) goto 51
           if (ABS(i0 - k0) == 1 .and. ABS(j0 - l0) == 1 .and. &
           & ABS(i0/2 - k0/2) == 1 .and. ABS(j0/2 - l0/2) == 1) goto 51

           save = i0
           i0 = k0
           k0 = save
           save = j0
           j0 = l0
           l0 = save

           isym = irpmo(l0)

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

51     End do

       goto 30

20     ired = inz

       Do inz = 1, ired - 1
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(l0)

!           write(*,*) 'type1',i0,j0,k0,l0,cint2_0
           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

           if (i0 == k0 .and. j0 == l0) goto 52
           if (ABS(i0 - k0) == 1 .and. ABS(j0 - l0) == 1 .and. &
           & ABS(i0/2 - k0/2) == 1 .and. ABS(j0/2 - l0/2) == 1) goto 52

           i0 = k(inz)
           j0 = l(inz)
           k0 = i(inz)
           l0 = j(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(l0)
!           write(*,*) 'type2',i0,j0,k0,l0,cint2_0

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

52         Call takekr(i0, j0, k0, l0, cint2_0)
           isym = irpmo(l0)
!           write(*,*) 'type3',i0,j0,k0,l0,cint2_0

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

           if (i0 == k0 .and. j0 == l0) goto 53
           if (ABS(i0 - k0) == 1 .and. ABS(j0 - l0) == 1 .and. &
           & ABS(i0/2 - k0/2) == 1 .and. ABS(j0/2 - l0/2) == 1) goto 53

           save = i0
           i0 = k0
           k0 = save
           save = j0
           j0 = l0
           l0 = save

           isym = irpmo(l0)
!           write(*,*) 'type4',i0,j0,k0,l0,cint2_0

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

53     End do

       close (1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

       Do l0 = li, le
           Do k0 = ki, ke
               Do j0 = ji, je
                   Do i0 = ii, ie
                       if (ABS(traint2(i0, j0, k0, l0)) > thresd) then
                           write (1, '(4I4, 2e20.10)') i0, j0, k0, l0, traint2(i0, j0, k0, l0)
                       end if
                   End do
               End do
           End do
       End do

       close (1)

       traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and second index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

31     Do inz = 1, nmaxint
           read (1, '(4I4,2e20.10)', err=10, end=21) i(inz), j(inz), k(inz), l(inz), cint2(inz)
!           write(*,*) i(inz),j(inz),k(inz),l(inz),cint2(inz)
       End do

       Do inz = 1, nmaxint
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(k0)

           Do knew = 1, nsym(spk, isym)
               k1 = indsym(spk, isym, knew)
               traint2(i0, j0, k1, l0) = traint2(i0, j0, k1, l0) + cint2_0*DCONJG(f(k0, k1))
           End do

       End do

       goto 31

21     ired = inz
!        write(*,*)'ired',ired
       Do inz = 1, ired - 1
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(k0)
!           write(*,*)i0,j0,k0,l0,isym, nsym(spk,isym)

           Do knew = 1, nsym(spk, isym)
               k1 = indsym(spk, isym, knew)
!              write(*,*)i0,j0,k1,l0,traint2(i0,j0,k1,l0)
               traint2(i0, j0, k1, l0) = traint2(i0, j0, k1, l0) + cint2_0*DCONJG(f(k0, k1))
           End do

       End do

       close (1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

       Do l0 = li, le
           Do k0 = ki, ke
               Do j0 = ji, je
                   Do i0 = ii, ie
                       if (ABS(traint2(i0, j0, k0, l0)) > thresd) then
                           write (1, '(4I4, 2e20.10)') i0, j0, k0, l0, traint2(i0, j0, k0, l0)
                       end if
                   End do
               End do
           End do
       End do

       close (1)

       traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and third index transformation  !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

32     Do inz = 1, nmaxint
           read (1, '(4I4,2e20.10)', err=10, end=22) i(inz), j(inz), k(inz), l(inz), cint2(inz)
       End do

       Do inz = 1, nmaxint
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(j0)

           Do jnew = 1, nsym(spj, isym)
               j1 = indsym(spj, isym, jnew)
               traint2(i0, j1, k0, l0) = traint2(i0, j1, k0, l0) + cint2_0*f(j0, j1)
           End do

       End do

       goto 32

22     ired = inz
!        write(*,*)'ired',ired
       Do inz = 1, ired - 1
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(j0)

           Do jnew = 1, nsym(spj, isym)
               j1 = indsym(spj, isym, jnew)
               traint2(i0, j1, k0, l0) = traint2(i0, j1, k0, l0) + cint2_0*f(j0, j1)
           End do

       End do

       close (1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

       Do l0 = li, le
           Do k0 = ki, ke
               Do j0 = ji, je
                   Do i0 = ii, ie
                       if (ABS(traint2(i0, j0, k0, l0)) > thresd) then
                           write (1, '(4I4, 2e20.10)') i0, j0, k0, l0, traint2(i0, j0, k0, l0)
                       end if
                   End do
               End do
           End do
       End do

       close (1)

       traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and fourth index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

33     Do inz = 1, nmaxint
           read (1, '(4I4,2e20.10)', err=10, end=23) i(inz), j(inz), k(inz), l(inz), cint2(inz)
       End do

       Do inz = 1, nmaxint
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(i0)

           Do inew = 1, nsym(spi, isym)
               i1 = indsym(spi, isym, inew)
               traint2(i1, j0, k0, l0) = traint2(i1, j0, k0, l0) + cint2_0*DCONJG(f(i0, i1))
           End do

       End do

       goto 33

23     ired = inz
!        write(*,*)'ired',ired
       Do inz = 1, ired - 1
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpmo(i0)

           Do inew = 1, nsym(spi, isym)
               i1 = indsym(spi, isym, inew)
!              if(i1 > k0) then
               traint2(i1, j0, k0, l0) = traint2(i1, j0, k0, l0) + cint2_0*DCONJG(f(i0, i1))
!              endif
           End do

       End do

       close (1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

       Do l0 = li, le
           Do k0 = ki, ke
               Do j0 = ji, je
                   Do i0 = ii, ie
                       if (ABS(traint2(i0, j0, k0, l0)) > thresd) then
                           write (1, '(4I4, 2e20.10)') i0, j0, k0, l0, traint2(i0, j0, k0, l0)
                       end if
                   End do
               End do
           End do
       End do

       close (1)

       goto 100
!  10     write(*,*)'error opening file' ; goto 1000
10     write (*, *) 'error opening file'
       inquire (1, opened=is_opened)
       if (is_opened .eqv. .true.) then
           close (1)
       end if

100    deallocate (i); Call memminus(KIND(i), SIZE(i), 1)
       deallocate (j); Call memminus(KIND(j), SIZE(j), 1)
       deallocate (k); Call memminus(KIND(k), SIZE(k), 1)
       deallocate (l); Call memminus(KIND(l), SIZE(l), 1)

       deallocate (cint2); Call memminus(KIND(cint2), SIZE(cint2), 2)
       deallocate (traint2); Call memminus(KIND(traint2), SIZE(traint2), 2)

       deallocate (indsym); Call memminus(KIND(indsym), SIZE(indsym), 1)
       deallocate (nsym); Call memminus(KIND(nsym), SIZE(nsym), 1)

1000 end subroutine intra_2

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE intra_3(spi, spj, spk, spl, fname)
!
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

       use four_caspt2_module

       Implicit NONE
       include 'mpif.h'
       integer, intent(in)        :: spi, spj, spk, spl
       character*50, intent(in)    :: fname

       integer, allocatable :: i(:), j(:), k(:), l(:), indsym(:, :, :), nsym(:, :)
       complex*16, allocatable :: cint2(:), traint2(:, :, :, :)

       real*8                  :: thresd
       complex*16              :: cint2_0

       integer :: i0, j0, k0, l0, i1, j1, k1, l1, inew, jnew, knew, lnew
       integer :: ii, ji, ki, li, ie, je, ke, le, lkr0, kkr0
       integer :: inz, nmx, ini(3), end(3), isp, isym, imo, nmaxint, ired, save

       logical :: is_opened
       integer :: n_cnt
       thresd = 1.0d-15

       if (.not. (spk == spl)) then
           write (*, *) 'error intra_3', spi, spj, spk, spl
           stop
       end if
       write (*, *) 'start intra_3, filename : ', trim(fname), ' rank :', rank
       ini(1) = 1
       end(1) = ninact
       ini(2) = ninact + 1
       end(2) = ninact + nact
       ini(3) = ninact + nact + 1
       end(3) = ninact + nact + nsec

       nmx = max(ninact, nact, nsec)
       Allocate (indsym(3, nsymrpa, nmx)); Call memplus(KIND(indsym), SIZE(indsym), 1)
       Allocate (nsym(3, nsymrpa)); Call memplus(KIND(nsym), SIZE(nsym), 1)
       indsym = 0
       nsym = 0

       Do isym = 1, nsymrpa
           Do isp = 1, 3                           ! i th space( inact or act or sec)
               Do imo = ini(isp), end(isp)
                   if (irpamo(imo) == isym) then
                       nsym(isp, isym) = nsym(isp, isym) + 1
                       indsym(isp, isym, nsym(isp, isym)) = imo
                   end if
               End do
           End do
       End do

       ii = ini(spi)
       ie = end(spi)
       ji = ini(spj)
       je = end(spj)
       ki = ini(spk)
       ke = end(spk)
       li = ini(spl)
       le = end(spl)

       Allocate (traint2(ii:ie, ji:je, ki:ke, li:le)); Call memplus(KIND(traint2), SIZE(traint2), 2)

!modify iwamuro
!        write(*,'("intra_3",8I4)')ii,ie,ji,je,ki,ke,li,le

       traint2 = 0.0d+00
       write (*, '("Current Memory is ",F10.2,"MB")') tmem/1024/1024

!        nmaxint = AINT(3.5d+09 - tmem)/48 ! one integrals required four integer and one complex values

! Abe modified 2016. 11.11  100 GB is the max memory
       nmaxint = AINT(5.00d+10 - tmem)/48 ! one integrals required four integer and one complex values
! Abe modified 2016. 11.11

       write (*, *) 'nmaxint = ', nmaxint
       write (*, *) 'tmem, nmaxint*48', tmem, nmaxint*48
       IF (nmaxint < 0) stop

       write (*, *) 'nmaxint = ', nmaxint

       Allocate (i(nmaxint)); Call memplus(KIND(i), SIZE(i), 1)
       Allocate (j(nmaxint)); Call memplus(KIND(j), SIZE(j), 1)
       Allocate (k(nmaxint)); Call memplus(KIND(k), SIZE(k), 1)
       Allocate (l(nmaxint)); Call memplus(KIND(l), SIZE(l), 1)

       Allocate (cint2(nmaxint)); Call memplus(KIND(cint2), SIZE(cint2), 2)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and first index transformation !
!                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')
30     Do inz = 1, nmaxint

           read (1, '(4I4,2e20.10)', err=10, end=20) i(inz), j(inz), k(inz), l(inz), cint2(inz)
!Iwamuro modify
!           write(*,'("cint2debugC",4I4,2E20.10)')i(inz),j(inz),k(inz),l(inz),cint2(inz)
       End do
       Do inz = 1, nmaxint
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpamo(l0)

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

           Call takekr(i0, j0, k0, l0, cint2_0)
           isym = irpamo(l0)
!Iwamuro modify
!           write(*,'("takekr",4I4,2E20.10)')i0,j0,k0,l0,cint2_0

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

           if (ABS(k0 - l0) == 1 .and. ABS(k0/2 - l0/2) == 1) goto 50

           i0 = i(inz)
           j0 = j(inz)
           if (mod(k(inz), 2) == 0) l0 = k(inz) - 1
           if (mod(k(inz), 2) == 1) l0 = k(inz) + 1
           if (mod(l(inz), 2) == 0) k0 = l(inz) - 1
           if (mod(l(inz), 2) == 1) k0 = l(inz) + 1
           cint2_0 = (-1.0d+00)**mod(k(inz) + l(inz), 2)*cint2(inz)
           isym = irpamo(l0)

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

           Call takekr(i0, j0, k0, l0, cint2_0)
           isym = irpamo(l0)

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

50     End do

       goto 30

20     ired = inz
       Do inz = 1, ired - 1
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpamo(l0)

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do
!Iwamuro modify
!           write(*,'("type 1",4I4,2E20.10)')i0,j0,k0,l0,cint2_0

           Call takekr(i0, j0, k0, l0, cint2_0)
           isym = irpamo(l0)

!Iwamuro modify
!           write(*,'("type 2",4I4,2E20.10)')i0,j0,k0,l0,cint2_0
           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

           if (ABS(k0 - l0) == 1 .and. ABS(k0/2 - l0/2) == 1) goto 51

           if (mod(k(inz), 2) == 0) l0 = k(inz) - 1
           if (mod(k(inz), 2) == 1) l0 = k(inz) + 1
           if (mod(l(inz), 2) == 0) k0 = l(inz) - 1
           if (mod(l(inz), 2) == 1) k0 = l(inz) + 1
           i0 = i(inz)
           j0 = j(inz)
           cint2_0 = (-1.0d+00)**mod(k(inz) + l(inz), 2)*cint2(inz)
           isym = irpamo(l0)
!Iwamuro modify
!           write(*,'("type 3",4I4,2E20.10)')i0,j0,k0,l0,cint2_0

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

           Call takekr(i0, j0, k0, l0, cint2_0)
           isym = irpamo(l0)
!Iwamuro momdify
!           write(*,'("type 4",4I4,2E20.10)')i0,j0,k0,l0,cint2_0

           Do lnew = 1, nsym(spl, isym)
               l1 = indsym(spl, isym, lnew)
               traint2(i0, j0, k0, l1) = traint2(i0, j0, k0, l1) + cint2_0*f(l0, l1)
           End do

51     End do

       close (1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
       write (*, '(a,a,10i4)') trim(fname), ' storing integrals to disk', rank, 1, li, le, ki, ke, ji, je, ii, ie

        n_cnt = 0
       Do l0 = li, le
           Do k0 = ki, ke
               Do j0 = ji, je
                   Do i0 = ii, ie
                    if (ABS(traint2(i0, j0, k0, l0)) > thresd .and. mod(n_cnt,nprocs)==rank) then
                        write (1, '(4I4, 2e20.10)') i0, j0, k0, l0, traint2(i0, j0, k0, l0)
                    end if
                    n_cnt = n_cnt + 1
                   End do
               End do
           End do
       End do

       close (1)

       traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and second index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       write (*, *) 'Read intergals  and second index transformation'
       open (1, file=trim(fname), status='old', form='formatted')

31     Do inz = 1, nmaxint
           read (1, '(4I4,2e20.10)', err=10, end=21) i(inz), j(inz), k(inz), l(inz), cint2(inz)
!Iwamuro modify
!           write(*,*) i(inz),j(inz),k(inz),l(inz),cint2(inz)
       End do

       Do inz = 1, nmaxint
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpamo(k0)

           Do knew = 1, nsym(spk, isym)
               k1 = indsym(spk, isym, knew)
               traint2(i0, j0, k1, l0) = traint2(i0, j0, k1, l0) + cint2_0*DCONJG(f(k0, k1))
           End do

       End do

       goto 31

21     ired = inz
!        write(*,*)'ired',ired
       Do inz = 1, ired - 1
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpamo(k0)
!           write(*,*)i0,j0,k0,l0,isym, nsym(spk,isym)

           Do knew = 1, nsym(spk, isym)
               k1 = indsym(spk, isym, knew)
!              write(*,*)i0,j0,k1,l0,traint2(i0,j0,k1,l0)
               traint2(i0, j0, k1, l0) = traint2(i0, j0, k1, l0) + cint2_0*DCONJG(f(k0, k1))
           End do

       End do

       close (1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       open (1, file=trim(fname), status='old', form='formatted')
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
       write (*, '(a,a,10i4)') trim(fname), ' storing integrals to disk', rank, 2, li, le, ki, ke, ji, je, ii, ie

       Do l0 = li, le
           Do k0 = ki, ke
               Do j0 = ji, je
                   Do i0 = ii, ie
                       if (ABS(traint2(i0, j0, k0, l0)) > thresd) then
                           write (1, '(4I4, 2e20.10)') i0, j0, k0, l0, traint2(i0, j0, k0, l0)
                       end if
                   End do
               End do
           End do
       End do

       close (1)

       traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and third index transformation  !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

32     Do inz = 1, nmaxint
           read (1, '(4I4,2e20.10)', err=10, end=22) i(inz), j(inz), k(inz), l(inz), cint2(inz)
       End do

       Do inz = 1, nmaxint
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpamo(j0)

           Do jnew = 1, nsym(spj, isym)
               j1 = indsym(spj, isym, jnew)
               traint2(i0, j1, k0, l0) = traint2(i0, j1, k0, l0) + cint2_0*f(j0, j1)
           End do

       End do

       goto 32

22     ired = inz
!        write(*,*)'ired',ired
       Do inz = 1, ired - 1
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpamo(j0)

           Do jnew = 1, nsym(spj, isym)
               j1 = indsym(spj, isym, jnew)
               traint2(i0, j1, k0, l0) = traint2(i0, j1, k0, l0) + cint2_0*f(j0, j1)
           End do

       End do

       close (1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
       write (*, '(a,a,10i4)') trim(fname), ' storing integrals to disk', rank, 3, li, le, ki, ke, ji, je, ii, ie

       Do l0 = li, le
           Do k0 = ki, ke
               Do j0 = ji, je
                   Do i0 = ii, ie
                       if (ABS(traint2(i0, j0, k0, l0)) > thresd) then
                           write (1, '(4I4, 2e20.10)') i0, j0, k0, l0, traint2(i0, j0, k0, l0)
                       end if
                   End do
               End do
           End do
       End do

       close (1)

       traint2 = 0.0d+00

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Read intergals  and fourth index transformation !
!                                                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')

33     Do inz = 1, nmaxint
           read (1, '(4I4,2e20.10)', err=10, end=23) i(inz), j(inz), k(inz), l(inz), cint2(inz)
       End do

       Do inz = 1, nmaxint
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpamo(i0)

           Do inew = 1, nsym(spi, isym)
               i1 = indsym(spi, isym, inew)
               traint2(i1, j0, k0, l0) = traint2(i1, j0, k0, l0) + cint2_0*DCONJG(f(i0, i1))
           End do

       End do

       goto 33

23     ired = inz
!        write(*,*)'ired',ired
       Do inz = 1, ired - 1
           i0 = i(inz)
           j0 = j(inz)
           k0 = k(inz)
           l0 = l(inz)
           cint2_0 = cint2(inz)
           isym = irpamo(i0)

           Do inew = 1, nsym(spi, isym)
               i1 = indsym(spi, isym, inew)
               traint2(i1, j0, k0, l0) = traint2(i1, j0, k0, l0) + cint2_0*DCONJG(f(i0, i1))
           End do

       End do

       close (1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Storing integrals to disk
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       open (1, file=trim(fname), status='old', form='formatted')
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
       write (*, '(a,a,10i4)') trim(fname), ' storing integrals to disk', rank, 4, li, le, ki, ke, ji, je, ii, ie

       Do l0 = li, le
           Do k0 = ki, ke
               Do j0 = ji, je
                   Do i0 = ii, ie
                       if (ABS(traint2(i0, j0, k0, l0)) > thresd) then
                           write (1, '(4I4, 2e20.10)') i0, j0, k0, l0, traint2(i0, j0, k0, l0)
                       end if
                   End do
               End do
           End do
       End do

       close (1)

       goto 100
!  10     write(*,*)'error opening file' ; goto 1000
10     write (*, *) 'error opening file', rank
       inquire (1, opened=is_opened)
       if (is_opened .eqv. .true.) then
           close (1)
       end if
       goto 101
100    call MPI_Barrier(MPI_COMM_WORLD, ierr)
       write (*, *) 'read and write file properly. filename : ', trim(fname), ' rank :', rank
101    deallocate (i); Call memminus(KIND(i), SIZE(i), 1)
       deallocate (j); Call memminus(KIND(j), SIZE(j), 1)
       deallocate (k); Call memminus(KIND(k), SIZE(k), 1)
       deallocate (l); Call memminus(KIND(l), SIZE(l), 1)

       deallocate (cint2); Call memminus(KIND(cint2), SIZE(cint2), 2)
       deallocate (traint2); Call memminus(KIND(traint2), SIZE(traint2), 2)

       deallocate (indsym); Call memminus(KIND(indsym), SIZE(indsym), 1)
       deallocate (nsym); Call memminus(KIND(nsym), SIZE(nsym), 1)

1000 end subroutine intra_3
