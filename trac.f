
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE traci( fa )  ! Transform CI matrix for new spinor basis

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
        real*8,      intent(in)  :: fa(ninact+1:ninact+nact,ninact+1:ninact+nact)


        integer :: i0, j0, i, info, job
        integer :: ii, jj, ok
        integer :: occ(nelec, ndet)

        integer, allocatable     :: IPIV(:)
        complex*16, Allocatable  :: ds(:,:), dsold(:,:), ci(:), work(:), z(:)
        complex*16  :: det(2)
        logical     :: error

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        occ = 0
        write(*,*)'Enter TRACI'

        Do i0 = 1, ndet
           i = 0
           ok = 0
           Do j0 = 0, 31
              if(btest(idet(i0),j0)) then
                 i = i+1
                 Do ii = 1, nact
                    if( ii == j0+1 ) then  ! j0+1 means occupied spinor labeled by casci
                       occ(i, i0) = ii         ! This is energetic order inside active spinor!
                       ok = ok + 1
                       goto 200
                    End if
                 End do
                       
 200          endif
           End do
        End do 

        Allocate(ds(ndet, ndet))

        ds = 0.0d+00

        Do i0 = 1, ndet     ! k  (old)
           Do j0 = 1, ndet  ! k~ (new)   <k|k~>

              Call dets (fa(ninact+1:ninact+nact,ninact+1:ninact+nact), &
                         occ(1:nelec,i0),occ(1:nelec,j0),ds(i0,j0))

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

        write(*,*)'Obtain inverse of ds matrix'

        Allocate (IPIV(ndet))
        Allocate(dsold(ndet, ndet))

        dsold = ds

        Call ZGETRF( ndet, ndet, ds, ndet, IPIV, INFO )!      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
        write(*,*)'info',info

           

        Allocate (work(ndet))

        Call ZGETRI( ndet, ds, ndet, IPIV, WORK, ndet, INFO )
        write(*,*)'info',info

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

        Deallocate(work)
        Deallocate(IPIV)
           write(*,*)'Check whether inverese matrix is really so'

        error = .FALSE.

        dsold = MATMUL( ds, dsold)
        Do i0 = 1, ndet
        Do j0 = 1, ndet


           If((i0 /= j0).and. ABS(dsold(i0,j0)) > 1.0d-10) then
              error = .TRUE.
              write(*,'(2I4,2E13.5)')i0,j0,dsold(i0,j0)
           Elseif(i0 == j0 .and. ABS(dsold(i0,j0)-1.0d+00) > 1.0d-10) then
              error = .TRUE.
              write(*,'(2I4,2E13.5)')i0,j0,dsold(i0,j0)
           Endif

        End do            
        End do            

        If(.not. error)  write(*,*) 'Inverse matrix is obtained correclty'

        Deallocate(dsold)

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
        ci = DCMPLX(cir(1:ndet,selectroot), cii(1:ndet,selectroot))
        ci = MATMUL ( ds, ci)
        cir(1:ndet,selectroot) = DBLE(ci)
        cii(1:ndet,selectroot) = DIMAG(ci)

        Deallocate (ci)

        Deallocate (ds)

        End subroutine traci


! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE tracic( fac )  ! Transform CI matrix for new spinor basis

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
        complex*16, intent(in)  :: fac(ninact+1:ninact+nact,ninact+1:ninact+nact)


        integer :: i0, j0, i, info, job
        integer :: ii, jj, ok
        integer :: occ(nelec, ndet)

        integer, allocatable     :: IPIV(:)
        complex*16, Allocatable  :: ds(:,:), dsold(:,:), ci(:), work(:), z(:)
        complex*16  :: det(2)
        logical     :: error

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

        occ = 0
        write(*,*)'Enter TRACI'

        Do i0 = 1, ndet
           i = 0
           ok = 0
           Do j0 = 0, 31
              if(btest(idet(i0),j0)) then
                 i = i+1
                 Do ii = 1, nact
                    if( ii == j0+1 ) then  ! j0+1 means occupied spinor labeled by casci
                       occ(i, i0) = ii         ! This is energetic order inside active spinor!
                       ok = ok + 1
                       goto 200
                    End if
                 End do
                       
 200          endif
           End do
        End do 

        Allocate(ds(ndet, ndet))

        ds = 0.0d+00

        Do i0 = 1, ndet     ! k  (old)
           Do j0 = 1, ndet  ! k~ (new)   <k|k~>

              Call detsc (fac(ninact+1:ninact+nact,ninact+1:ninact+nact), &
                         occ(1:nelec,i0),occ(1:nelec,j0),ds(i0,j0))

           End do
        End do


        write(*,*)'Obtain inverse of ds matrix'

        Allocate (IPIV(ndet))
        Allocate(dsold(ndet, ndet))

        dsold = ds

        Call ZGETRF( ndet, ndet, ds, ndet, IPIV, INFO )
        write(*,*)'info',info

        Allocate (work(ndet))

        Call ZGETRI( ndet, ds, ndet, IPIV, WORK, ndet, INFO )
        write(*,*)'info',info


        Deallocate(work)
        Deallocate(IPIV)
           write(*,*)'Check whether inverese matrix is really so'

        error = .FALSE.

        dsold = MATMUL( ds, dsold)
        Do i0 = 1, ndet
        Do j0 = 1, ndet


           If((i0 /= j0).and. ABS(dsold(i0,j0)) > 1.0d-10) then
              error = .TRUE.
              write(*,'(2I4,2E13.5)')i0,j0,dsold(i0,j0)
           Elseif(i0 == j0 .and. ABS(dsold(i0,j0)-1.0d+00) > 1.0d-10) then
              error = .TRUE.
              write(*,'(2I4,2E13.5)')i0,j0,dsold(i0,j0)
           Endif

        End do            
        End do            

        If(.not. error)  write(*,*) 'Inverse matrix is obtained correclty'

        Deallocate(dsold)

!        Now ds is inverse matrix!

        Allocate (ci(ndet))

        ci = 0.0d+00
        ci = DCMPLX(cir(1:ndet,selectroot), cii(1:ndet,selectroot))
        ci = MATMUL ( ds, ci)
        cir(1:ndet,selectroot) = DBLE(ci(1:ndet))
        cii(1:ndet,selectroot) = DIMAG(ci(1:ndet))

        open (5, file = 'NEWCICOEFF', status = 'unknown', form = 'unformatted')
        write(5) ci(1:ndet)
!        write(*,'("ci",2E20.10)') ci(1:ndet)
        close(5)


        Deallocate (ci)

        Deallocate (ds)

        End subroutine tracic


