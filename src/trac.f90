
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE traci(fa)  ! Transform CI matrix for new spinor basis

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
    real*8, intent(in)  :: fa(ninact + 1:ninact + nact, ninact + 1:ninact + nact)

    integer :: i0, j0, i, info
    integer :: ii, ok
    integer :: occ(nelec, ndet)

    integer, allocatable     :: IPIV(:)
    complex*16, Allocatable  :: ds(:, :), dsold(:, :), ci(:), work(:)
    logical     :: error

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    occ = 0
    if (rank == 0) then ! Process limits for output
        write (*, *) 'Enter TRACI'
    end if

    Do i0 = 1, ndet
        i = 0
        ok = 0
        Do j0 = 0, 31
            if (btest(idet(i0), j0)) then
                i = i + 1
                Do ii = 1, nact
                    if (ii == j0 + 1) then  ! j0+1 means occupied spinor labeled by casci
                        occ(i, i0) = ii         ! This is energetic order inside active spinor!
                        ok = ok + 1
                        goto 200
                    End if
                End do

200         end if
        End do
    End do

    Allocate (ds(ndet, ndet))

    ds = 0.0d+00

    Do i0 = 1, ndet     ! k  (old)
        Do j0 = 1, ndet  ! k~ (new)   <k|k~>

            Call dets(fa(ninact + 1:ninact + nact, ninact + 1:ninact + nact), &
                      occ(1:nelec, i0), occ(1:nelec, j0), ds(i0, j0))

        End do
    End do

! for a while !        write(*,'(/,"REAL")')
! for a while !        Do i0 = 1, ndet
! for a while !           write(*,'(5E13.5,/,5E13.5,/,5E13.5, &
! for a while !           & /,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5, &
! for a while !           & /,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5, &
! for a while !           & /,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5, &
! for a while !           & /,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5  &
! for a while !           & )') (DBLE(ds(i0,j0)),j0 = 1,ndet)
! for a while !!!!           write(*,'(5F13.5,/,5F13.5,/,5F13.5, &
! for a while !!!!           & /,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5, &
! for a while !!!!           & /,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5, &
! for a while !!!!           & /,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5, &
! for a while !!!!           & /,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5  &
! for a while !!!!           & )') (DBLE(ds(i0,j0)),j0 = 1,ndet)
! for a while !!!!!           write(*,'(5E13.5,/,5E13.5,/,5E13.5,/)') (DBLE(ds(i0,j0)),j0 = 1,ndet)
! for a while !        End do

! for a while !        write(*,'(/,"REAL")')
! for a while !        Do i0 = 1, ndet
! for a while !        Do j0 = 1, ndet
! for a while !           if(ABS(DBLE(ds(i0,j0)))> 1.0d-10) &
! for a while !           write(*,'("ds",2I4,E13.5)') i0, j0, DBLE(ds(i0,j0))
! for a while !        End do
! for a while !        End do
! for a while !
! for a while !        write(*,'(/,"REAL")')
! for a while !        Do i0 = 1, ndet
! for a while !        Do j0 = i0, ndet
! for a while !           if(DBLE(ds(i0,j0))> 1.0d-10) then
! for a while !              if( ABS(ds(i0,j0)-ds(j0,i0)) < 1.0d-5 ) then
! for a while !                 write(*,'(2I4,2E13.5)') i0, j0, DBLE(ds(i0,j0)),DBLE(ds(j0,i0))
! for a while !              End if
! for a while !           End if
! for a while !        End do
! for a while !        End do
! for a while !
! for a while !        write(*,'(/,"IMAG")')
! for a while !        Do i0 = 1, ndet
! for a while !        Do j0 = 1, ndet
! for a while !           if(ABS(DIMAG(ds(i0,j0)))> 1.0d-10) &
! for a while !           write(*,'(E13.5)') DIMAG(ds(i0,j0))
! for a while !        End do
! for a while !        End do

    if (rank == 0) then ! Process limits for output
        write (*, *) 'Obtain inverse of ds matrix'
    end if

    Allocate (IPIV(ndet))
    Allocate (dsold(ndet, ndet))

    dsold = ds

    Call ZGETRF(ndet, ndet, ds, ndet, IPIV, INFO)!      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
    if (rank == 0) then ! Process limits for output
        write (*, *) 'info', info
    end if

    Allocate (work(ndet))

    Call ZGETRI(ndet, ds, ndet, IPIV, WORK, ndet, INFO)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'info', info
    end if

! for a while !      write(*,'(/,"REAL")')
! for a while !        Do i0 = 1, ndet
! for a while !!           write(*,'(5F13.5,/,5F13.5,/,5F13.5, &
! for a while !!           & /,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5, &
! for a while !!           & /,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5, &
! for a while !!           & /,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5, &
! for a while !!           & /,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5,/,5F13.5  &
! for a while !!           & )') (DBLE(ds(i0,j0)),j0 = 1,ndet)
! for a while !
! for a while !           write(*,'(5E13.5,/,5E13.5,/,5E13.5, &
! for a while !           & /,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5, &
! for a while !           & /,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5, &
! for a while !           & /,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5, &
! for a while !           & /,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5,/,5E13.5  &
! for a while !           & )') (DBLE(ds(i0,j0)),j0 = 1,ndet)
! for a while !!           write(*,'(5E13.5,/,5E13.5,/,5E13.5,/)') (DBLE(ds(i0,j0)),j0 = 1,ndet)
! for a while !        End do
! for a while !        write(*,'(/,"IMAG")')
! for a while !        Do i0 = 1, ndet
! for a while !        Do j0 = 1, ndet
! for a while !           if(ABS(DIMAG(ds(i0,j0)))> 1.0d-10) &
! for a while !           write(*,'(E13.5)') DIMAG(ds(i0,j0))
! for a while !        End do
! for a while !        End do

    Deallocate (work)
    Deallocate (IPIV)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'Check whether inverese matrix is really so'
    end if

    error = .FALSE.

    dsold = MATMUL(ds, dsold)
    Do i0 = 1, ndet
        Do j0 = 1, ndet

            If ((i0 /= j0) .and. ABS(dsold(i0, j0)) > 1.0d-10) then
                error = .TRUE.
                if (rank == 0) then ! Process limits for output
                    write (*, '(2I4,2E13.5)') i0, j0, dsold(i0, j0)
                end if
            Elseif (i0 == j0 .and. ABS(dsold(i0, j0) - 1.0d+00) > 1.0d-10) then
                error = .TRUE.
                if (rank == 0) then ! Process limits for output
                    write (*, '(2I4,2E13.5)') i0, j0, dsold(i0, j0)
                end if
            End if

        End do
    End do

    if (rank == 0) then ! Process limits for output
        If (.not. error) write (*, *) 'Inverse matrix is obtained correclty'
    end if
    Deallocate (dsold)

!        Now ds is inverse matrix!

!        Allocate (ci(ndet))
!
!        Do i0 = 1, nroot
!           ci = 0.0d+00
!           ci = DCMPLX(cir(1:ndet,i0), cii(1:ndet,i0))
!           ci = MATMUL ( ds, ci)
!           cir(1:ndet,i0) = DBLE(ci)
!           cii(1:ndet,i0) = DIMAG(ci)
!        End do
!
!
!        Deallocate (ci)

    Allocate (ci(ndet))

    ci = 0.0d+00
    ci = DCMPLX(cir(1:ndet, selectroot), cii(1:ndet, selectroot))
    ci = MATMUL(ds, ci)
    cir(1:ndet, selectroot) = DBLE(ci)
    cii(1:ndet, selectroot) = DIMAG(ci)

    Deallocate (ci)

    Deallocate (ds)

End subroutine traci

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE tracic(fac)  ! Transform CI matrix for new spinor basis

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
#ifdef HAVE_MPI
    include 'mpif.h'
#endif
    complex*16, intent(in)  :: fac(ninact + 1:ninact + nact, ninact + 1:ninact + nact)

    integer :: i0, j0, i, info
    integer :: ii, ok
    integer :: occ(nelec, ndet)

    integer, allocatable     :: IPIV(:)
    complex*16, Allocatable  :: ds(:, :), dsold(:, :), ci(:), work(:)
    integer :: datetmp0, datetmp1
    real(8) :: tsectmp0, tsectmp1
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    occ = 0
    if (rank == 0) then ! Process limits for output
        write (*, *) 'Enter TRACI'
    end if
    datetmp1 = date0; datetmp0 = date0

    Call timing(date0, tsec0, datetmp0, tsectmp0)
    tsectmp1 = tsectmp0
    Do i0 = 1, ndet
        i = 0
        ok = 0
        Do j0 = 0, 31
            if (btest(idet(i0), j0)) then
                i = i + 1
                Do ii = 1, nact
                    if (ii == j0 + 1) then  ! j0+1 means occupied spinor labeled by casci
                        occ(i, i0) = ii         ! This is energetic order inside active spinor!
                        ok = ok + 1
                        goto 200
                    End if
                End do

200         end if
        End do
    End do
    if (rank == 0) write (*, *) 'Before allocate a matrix named ds'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Allocate (ds(ndet, ndet))

    ds = 0.0d+00
    if (rank == 0) write (*, *) 'Initialized a matrix named ds'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    ! Noda ndet^2で回っているので遅くなりそう
    Do i0 = 1, ndet     ! k  (old)
        Do j0 = 1, ndet  ! k~ (new)   <k|k~>

            Call detsc(fac(ninact + 1:ninact + nact, ninact + 1:ninact + nact), &
                       occ(1:nelec, i0), occ(1:nelec, j0), ds(i0, j0))

        End do
    End do
    if (rank == 0) write (*, *) 'End detsc'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    if (rank == 0) then ! Process limits for output
        write (*, *) 'Obtain inverse of ds matrix'
    end if
    Allocate (IPIV(ndet))
    Allocate (dsold(ndet, ndet))

    dsold = ds
    if (rank == 0) write (*, *) 'Start get LU factorization of ds'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    ! Noda Attention :: 逆行列の計算の計算量 (n*n行列の場合)
    ! LU分解+前進後退代入による逆行列の計算        n^3+n-1n3+n−1
    ! ZGETRFはある行列に対してLU分解を与える
    ! ZGETRIはLU分解したものを使って逆行列を計算する
    ! つまりZGETRF+ZGETRIの計算量はO(n^3)でcdiagと同等の計算量が必要
    Call ZGETRF(ndet, ndet, ds, ndet, IPIV, INFO)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'info', info
    end if
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (*, *) 'End get LU factorization of ds'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Allocate (work(ndet))
    if (rank == 0) write (*, *) 'Start get a inverse matrix of ds'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Call ZGETRI(ndet, ds, ndet, IPIV, WORK, ndet, INFO)
    if (rank == 0) then ! Process limits for output
        write (*, *) 'info', info
    end if
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (*, *) 'End get a inverse matrix of ds, ndet', ndet
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Deallocate (work)
    Deallocate (IPIV)
#ifdef DEBUG
    if (rank == 0) then ! Process limits for output
        write (*, *) 'Check whether inverese matrix is really so'
    end if
    error = .FALSE.

    ! Noda ndet^2で回っているので遅くなりそう
    ! dsold=AA^(-1)?
    ! AA^(-1)=E => if(i0/=j0)dsold(i0,j0)=0,else dsold(i0,j0)=1.0 ??
    ! Matmulは一般にO(N^3)の計算量なので行列の積をとるだけの操作にO(N^3)かかっている
    dsold = MATMUL(ds, dsold)
#ifdef HAVE_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) write (*, *) 'End dsold = matmul(ds, dsold) so dsold should be a identity matrix.'
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    !    do i0 = 1, ndet
    !        dsold(i0, i0) = dsold(i0, i0) - 1.0d+00
    !    end do
    !    if (maxval(abs(real(dsold))) > 1.0d-10) then
    !        if (rank == 0) then
    !            write (*, '(E13.5)') maxval(abs(real(dsold)))
    !        end if
    !    end if
    Do i0 = 1, ndet
        Do j0 = 1, ndet

            If ((i0 /= j0) .and. ABS(dsold(i0, j0)) > 1.0d-10) then
                error = .TRUE.
                if (rank == 0) then ! Process limits for output
                    write (*, '(2I4,2E13.5)') i0, j0, dsold(i0, j0)
                end if
            Elseif (i0 == j0 .and. ABS(dsold(i0, j0) - 1.0d+00) > 1.0d-10) then
                error = .TRUE.
                if (rank == 0) then ! Process limits for output
                    write (*, '(2I4,2E13.5)') i0, j0, dsold(i0, j0)
                end if
            End if

        End do
    End do

    if (rank == 0) then ! Process limits for output
        If (.not. error) write (*, *) 'Inverse matrix is obtained correclty'
    end if
#endif
    Call timing(datetmp1, tsectmp1, datetmp0, tsectmp0)
    datetmp1 = datetmp0
    tsectmp1 = tsectmp0
    Deallocate (dsold)

!        Now ds is inverse matrix!

    Allocate (ci(ndet))

    ci = 0.0d+00
    ci = DCMPLX(cir(1:ndet, selectroot), cii(1:ndet, selectroot))
    ci = MATMUL(ds, ci)
    cir(1:ndet, selectroot) = DBLE(ci(1:ndet))
    cii(1:ndet, selectroot) = DIMAG(ci(1:ndet))
    if (rank == 0) then ! Only master ranks are allowed to create files used by CASPT2 except for MDCINTNEW.
        open (5, file='NEWCICOEFF', status='unknown', form='unformatted')
        write (5) ci(1:ndet)
        close (5)
    end if

    Deallocate (ci)

    Deallocate (ds)

End subroutine tracic
