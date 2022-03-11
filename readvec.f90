! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE readvec(filename)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE

    integer :: mdtriv, lenrec, ios, irec, nbitsa, midet
    character*50, intent(in) :: filename
    integer :: j0, j, i, i0, i1
    integer :: k0, l0, ii, jj, kk, ll
!        real*8 :: thres

    mdtriv = 10
    eigen(:) = 0.0d+00
    cir(:, :) = 0.0d+00
    cii(:, :) = 0.0d+00

    open (mdtriv, file=trim(filename), status='old', access='direct', recl=8, err=10)
    ios = 0
    read (mdtriv, rec=1, err=11, iostat=ios) lenrec
    if (ios .ne. 0) goto 12
    close (mdtriv)

    open (mdtriv, file=filename, access='direct', recl=lenrec, err=100)
    read (mdtriv, rec=1, err=100) lenrec, nroot

    Allocate (eigen(nroot))
    read (mdtriv, rec=1, err=100) lenrec, nroot, (eigen(i0), i0=1, nroot)

    read (mdtriv, rec=2, err=200) ndet

    Allocate (idet(ndet))

    read (mdtriv, rec=2, err=200) ndet, (idet(i), i=1, ndet)

!        write(*,*)  (idet(i), i=1,ndet)
!        do i = 1, ndet
!           write(*,*)(btest(idet(i),i0), i0=0,63)
!        end do

    midet = 0

!        do i0 = 1, ndet
!          midet = max0( midet, idet(i0))
!        end do
!
!        write(*,*) midet, 'midet'
!        write(*,*)(btest(midet,i0), i0=0,63)
!
!
!        do i0 = 63, 0, -1
!            if(BTEST(midet, i0)) then
!!              write(*,*) i0
!               norb = i0 + 1
!               goto 7
!            endif
!       end do

7   nelec = POPCNT(idet(1))
    write (*, *) POPCNT(idet(1)), idet(1)
    do i0 = 1, ndet
        if (POPCNT(idet(i0)) /= nelec) then
            write (*, *) 'error about nelec', nelec, idet(i0)
        end if
    end do

    Allocate (cir(ndet, nroot))
    Allocate (cii(ndet, nroot))

    do irec = 1, nroot
        read (mdtriv, rec=irec + 2, err=300) (cir(j, irec), cii(j, irec), j=1, ndet)
    end do

!        write(*,*)'norb=' ,norb
!        write(*,*)'nelec=' ,nelec
!        write(*,*)'nroot=',nroot
!        write(*,*)'ndet=',ndet

    do i0 = 1, nroot
        write (*, *) i0, eigen(i0)
    end do

    realcvec = .true.

    write (*, *) 'j,irec, cir(j,irec), cii(j,irec)'

    do irec = 1, nroot
        do j = 1, ndet
!!            write(*,'(2I4,2(3X,E14.7))') j,irec, cir(j,irec), cii(j,irec)
            if (ABS(cii(j, irec)) > thres) then
                realcvec = .false.
            end if
        end do
    end do

    goto 1

10  write (*, *) 'err 10'
    go to 1000
11  write (*, *) 'err 11'
    go to 1000
12  write (*, *) 'err 12'
    go to 1000

100 write (*, *) 'err 100 vec come'
    go to 1000

200 write (*, *) 'err 200'
    go to 1000

300 write (*, *) 'err 300'
    go to 1000

1   close (mdtriv)
1000 end subroutine readvec
