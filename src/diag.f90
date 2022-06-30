! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE rdiag(sr, dimn, dimm, w, thresd, cutoff)
! diagonalization of real symmetric matrix
!  and remove linear dependency for any S matrix

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE

    integer, intent(in) :: dimn
    real*8, intent(in)  :: thresd
    logical, intent(in) :: cutoff

    real*8, intent(inout)  :: sr(dimn, dimn)

    integer, intent(out) :: dimm
    real*8, intent(out)  ::  w(dimn)

    integer :: info, lda, lwork
    character :: jobz*1, uplo*1
    real*8, allocatable  ::  work(:)
    integer :: j0, i0

    w(:) = 0.0d+00
    jobz = 'V' ! calculate eigenvectors
    uplo = 'U' ! calculate upper triangle matrix

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0. = indxyz
!*  LDA     (input) INTEGER
!*          The leading dimension of the array A.  LDA >= max(1,N).
!*  LWORK   (input) INTEGER
!*          The length of the array WORK.  LWORK >= max(1,3*N-1).
!*          For optimal efficiency, LWORK >= (NB+2)*N,
!*          where NB is the blocksize for DSYTRD returned by ILAENV.
!*
!*          If LWORK = -1, then a workspace query is assumed; the routine
!*          only calculates the optimal size of the WORK array, returns
!*          this value as the first entry of the WORK array, and no error
!*          message related to LWORK is issued by XERBLA.
!*

    lda = max(1, dimn)
    lwork = max(1, 3*dimn - 1)

    allocate (work(lwork))

    call dsyev(jobz, uplo, dimn, sr, lda, w, work, lwork, info)

    deallocate (work)

    if (info /= 0 .and. rank == 0) then
        write (*, *) 'error in diagonalization, info = ', info
        goto 1000
    end if

    if (cutoff) then

        if (rank == 0) then ! Process limits for output
            write (*, *) 'cut off threshold is ', thresd
        end if
        j0 = 0
        do i0 = 1, dimn
            if (w(i0) >= thresd) then
                j0 = j0 + 1
            end if
        end do

        dimm = j0

    else
        dimm = dimn
    end if

1000 continue
end subroutine rdiag

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE cdiag(c, dimn, dimm, w, thresd, cutoff)
! diagonalization of complex symmetric matrix
!  and remove linear dependency for any S matrix

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE

    integer, intent(in) :: dimn
    real*8, intent(in)  :: thresd
    logical, intent(in) :: cutoff

    complex*16, intent(inout):: c(dimn, dimn)

    integer, intent(out) :: dimm
    real*8, intent(out)  :: w(dimn)

    integer :: info, lda, lwork
    character :: jobz*1, uplo*1

    complex*16, allocatable  ::  work(:)
    real*8, allocatable      ::  rwork(:)
    integer :: j0, i0

    if (rank == 0) then ! Process limits for output
        write (*, *) 'Enter cdiagonal part'
    end if
    w(:) = 0.0d+00

    jobz = 'V' ! calculate eigenvectors
    uplo = 'U' ! calculate upper triangle matrix

!zheev!     .. Scalar Arguments ..
!zheev!      CHARACTER          JOBZ, UPLO
!zheev!      INTEGER            INFO, LDA, LWORK, N
!zheev!*     ..
!zheev!*     .. Array Arguments ..
!zheev!      DOUBLE PRECISION   RWORK( * ), W( * )
!zheev!      COMPLEX*16         A( LDA, * ), WORK( * )
!zheev!*     ..
!zheev!*
!zheev!*  Purpose
!zheev!*  =======
!zheev!*
!zheev!*  ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
!zheev!*  complex Hermitian matrix A.
!zheev!*
!zheev!*  Arguments
!zheev!*  =========
!zheev!*
!zheev!*  JOBZ    (input) CHARACTER*1
!zheev!*          = 'N':  Compute eigenvalues only;
!zheev!*          = 'V':  Compute eigenvalues and eigenvectors.
!zheev!*
!zheev!*  UPLO    (input) CHARACTER*1
!zheev!*          = 'U':  Upper triangle of A is stored;
!zheev!*          = 'L':  Lower triangle of A is stored.
!zheev!*
!zheev!*  N       (input) INTEGER
!zheev!*          The order of the matrix A.  N >= 0.
!zheev!*
!zheev!*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
!zheev!*          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!zheev!*          leading N-by-N upper triangular part of A contains the
!zheev!*          upper triangular part of the matrix A.  If UPLO = 'L',
!zheev!*          the leading N-by-N lower triangular part of A contains
!zheev!*          the lower triangular part of the matrix A.
!zheev!*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!zheev!*          orthonormal eigenvectors of the matrix A.
!zheev!*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!zheev!*          or the upper triangle (if UPLO='U') of A, including the
!zheev!*          diagonal, is destroyed.
!zheev!*
!zheev!*  LDA     (input) INTEGER
!zheev!*          The leading dimension of the array A.  LDA >= max(1,N).
!zheev!*
!zheev!*  W       (output) DOUBLE PRECISION array, dimension (N)
!zheev!*          If INFO = 0, the eigenvalues in ascending order.
!zheev!*
!zheev!*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!zheev!*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!zheev!*
!zheev!*  LWORK   (input) INTEGER                                        =================================
!zheev!*          The length of the array WORK.  LWORK >= max(1,2*N-1).  !<= 3*N-1 is better like dsyev.f
!zheev!*          For optimal efficiency, LWORK >= (NB+1)*N,             =================================
!zheev!*          where NB is the blocksize for ZHETRD returned by ILAENV.
!zheev!*
!zheev!*          If LWORK = -1, then a workspace query is assumed; the routine
!zheev!*          only calculates the optimal size of the WORK array, returns
!zheev!*          this value as the first entry of the WORK array, and no error
!zheev!*          message related to LWORK is issued by XERBLA.
!zheev!*
!zheev!*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(1, 3*N-2))
!zheev!*
!zheev!*  INFO    (output) INTEGER
!zheev!*          = 0:  successful exit
!zheev!*          < 0:  if INFO = -i, the i-th argument had an illegal value
!zheev!*          > 0:  if INFO = i, the algorithm failed to converge; i
!zheev!*                off-diagonal elements of an intermediate tridiagonal
!zheev!*                form did not converge to zero.
!zheev!*
!zheev!*  =====================================================================
!zheev!*
!zheev!*

    lda = max(1, dimn)

    lwork = max(1, 3*dimn - 1)

    allocate (work(lwork))
    allocate (rwork(3*dimn - 2))

    work = 0.0d+00
    rwork = 0.0d+00

    Call ZHEEV(JOBZ, UPLO, dimn, c, LDA, W, WORK, LWORK, RWORK, INFO)

    deallocate (work)
    deallocate (rwork)

    if (rank == 0) then ! Process limits for output
        write (*, *) 'Finish zheev info = ', info
    end if
    if (info /= 0) then
        if (rank == 0) then ! Process limits for output
            write (*, *) 'error in diagonalization, info = ', info
        end if
        goto 1000
    end if

!        Do i0 = 1, dimn
!           write(*,'(I4,E20.10)')i0,w(i0)
!        End do

    if (cutoff) then

        if (rank == 0) then ! Process limits for output
            write (*, *) 'cut off threshold is ', thresd
        end if

        j0 = 0
        do i0 = 1, dimn
            if (ABS(w(i0)) >= thresd) then
                j0 = j0 + 1
            end if
        end do

        dimm = j0

    else
        dimm = dimn
    end if

    if (rank == 0) then ! Process limits for output
        write (*, *) "end cdiag"
    end if
1000 continue
end subroutine cdiag

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE rdiag0(n, n0, n1, fa, w)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
    integer, intent(in)     ::  n, n0, n1

    real*8, intent(out)     ::  fa(n0:n1, n0:n1)
    real*8, intent(out)     ::  w(n0:n1)

    logical                 ::  cutoff
    integer                 ::  j, i, dimn, ncount(nsymrp)
    integer                 ::  ii, sym, isym
    integer                 ::  ind(n, nsymrp)

    real*8, allocatable     ::  mat(:, :), fasym(:, :)
    real*8                  ::  wsym(n)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!  DAIAGONALIZATION OF A COMPLEX HERMITIAN MATRIX
    if (rank == 0) then ! Process limits for output
        write (*, *) 'rdiag0 start'
    end if
    w = 0.0d+00
    cutoff = .FALSE.

    fa(n0:n1, n0:n1) = 0.0d+00

!        DAIGONALIZE IN EACH SYMMETRY
!        FIRST, THE ELEMENTS ARE REORDERED BY SYMMETRIC ORDER

! SET NCOUNT(SYM) : DIMENSION OF EACH SYMMETRY

    ncount = 0

    if (rank == 0) then ! Process limits for output
        write (*, *) 'nsymrp', nsymrp
    end if
    Do sym = 1, nsymrp

        Do i = n0, n1
            ii = i
            isym = irpmo(ii)
            if (isym == sym) then
                ncount(sym) = ncount(sym) + 1
                ind(ncount(sym), sym) = i
            End if
        End do

    End do

    if (rank == 0) then ! Process limits for output
        write (*, *) 'sym,ncount(sym)', (ncount(sym), sym=1, nsymrp)
    end if
    Do sym = 1, nsymrp

        Allocate (fasym(ncount(sym), ncount(sym)))
        Do j = 1, ncount(sym)
            Do i = 1, ncount(sym)
                fasym(i, j) = f(ind(i, sym), ind(j, sym))
            End do
        End do

        dimn = ncount(sym)

        Call rdiag(fasym, dimn, dimn, wsym, thres, cutoff)
!      _________________________________________________________

        Do j = 1, ncount(sym)

            w(ind(j, sym)) = wsym(j)

            Do i = 1, ncount(sym)
                fa(ind(i, sym), ind(j, sym)) = fasym(i, j)
            End do

        End do

        Deallocate (fasym)

    End do ! sym

! NOW FA BECOMES TRANSFORM MATRIX   CONJG(Fa) Fbc Fa = W <= diagonal form!

    Allocate (mat(n, n))
    mat = 0.0d+00

    mat = TRANSPOSE(fa)
    mat = MATMUL(mat, f)
    mat = MATMUL(mat, fa)

    if (rank == 0) then ! Process limits for output
        write (*, *) 'OFF DIAGONAL TERM OF U*FU'
    end if
    do i = 1, n
        do j = 1, n
            if ((i /= j) .and. (ABS(mat(i, j)) > 1.0d-10)) then
                if (rank == 0) then ! Process limits for output
                    write (*, '(2E13.5,2I3)') mat(i, j), i, j
                end if
            end if
        end do
    end do

    if (rank == 0) then ! Process limits for output
        write (*, *) 'DIAGONAL TERM OF U*FU, W AND THEIR DIFFERENCE'
        do i = 1, n
            write (*, '(4E13.5)') mat(i, i), w(i), ABS(mat(i, i) - w(i))
        end do
        write (*, '(/)')
    end if
    deallocate (mat)

    if (rank == 0) then ! Process limits for output
        write (*, *) 'rdiag0 end'
    end if
end subroutine rdiag0

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE cdiag0(n, n0, n1, fac, wc)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
    integer, intent(in)     ::  n, n0, n1

    complex*16, intent(out) ::  fac(n0:n1, n0:n1)
    real*8, intent(out)     ::  wc(n0:n1)

    logical                 ::  cutoff, fi
    integer                 ::  j, i, dimn, ncount(nsymrpa)
    integer                 ::  ii, sym, isym
    integer                 ::  ind(n, nsymrpa)

    complex*16, allocatable ::  matc(:, :), facsym(:, :), facsymo(:, :)
    real*8, allocatable      ::  wcsym(:)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!  DAIAGONALIZATION OF A COMPLEX HERMITIAN MATRIX
    if (rank == 0) then ! Process limits for output
        write (*, *) 'cdiag0 start'
        write (*, *) 'nsymrpa', nsymrpa
    end if
!         nsymrp = nsymrpa

    wc = 0.0d+00
    cutoff = .FALSE.
    fi = .FALSE.

    Do i = n0, n1
        Do j = n0, n1
            if (ABS(DIMAG(f(i, j))) > 1.0d-10) fi = .TRUE.
        End do
    End do

    if (rank == 0) then ! Process limits for output
        write (*, *) 'fi', fi
    end if

    fac(n0:n1, n0:n1) = 0.0d+00

!        DAIGONALIZE IN EACH SYMMETRY
!        FIRST, THE ELEMENTS ARE REORDERED BY SYMMETRIC ORDER

! SET NCOUNT(SYM) : DIMENSION OF EACH SYMMETRY

    ncount = 0
    ind = 0

    Do sym = 1, nsymrpa

        Do i = n0, n1
            ii = i
            isym = irpmo(ii)
            If ((nsymrpa == 1) .or. &
                (nsymrpa /= 1 .and. (isym == sym))) then
                ncount(sym) = ncount(sym) + 1
                ind(ncount(sym), sym) = i
            End if
        End do
!            write(*,*)(ind(j,sym),j=1,ncount(sym))

    End do

!         write(*,*)'sym,ncount(sym)',(ncount(sym),sym=1,nsymrpa)

    Do sym = 1, nsymrpa

        Allocate (facsym(ncount(sym), ncount(sym)))
        facsym = 0.0d+00

        Do j = 1, ncount(sym)
            Do i = j, ncount(sym)
                facsym(i, j) = f(ind(i, sym), ind(j, sym))
                facsym(j, i) = DCONJG(f(ind(i, sym), ind(j, sym))) ! HERMITE
            End do
        End do

        dimn = ncount(sym)

        Allocate (facsymo(ncount(sym), ncount(sym)))
        facsymo = facsym
        Allocate (wcsym(ncount(sym)))
        wcsym = 0.0d+00
        cutoff = .FALSE.

        Call cdiag(facsym, dimn, dimn, wcsym, thres, cutoff)
!      _________________________________________________________

        facsym = DCONJG(facsym)
        facsymo = MATMUL(TRANSPOSE(facsym), facsymo)
        facsym = DCONJG(facsym)
        facsymo = MATMUL(facsymo, facsym)

        Do i = 1, dimn
            Do j = 1, dimn
                If (i /= j .and. ABS(facsymo(i, j)) > 1.0d-10) then
                    if (rank == 0) then ! Process limits for output
                        write (*, '("sym=",3I4,2E20.10)') sym, i, j, facsymo(i, j)
                    end if
                End if
            End do
        End do

        Do i = 1, dimn
            If (ABS(facsymo(i, i) - wcsym(i)) > 1.0d-10) then
                if (rank == 0) then ! Process limits for output
                    write (*, '("sym=",2I4,3E20.10)') sym, i, facsymo(i, i), wcsym(i)
                end if
            End if
        End do

        Deallocate (facsymo)

        Do j = 1, ncount(sym)

            wc(ind(j, sym)) = wcsym(j)

            Do i = 1, ncount(sym)
                fac(ind(i, sym), ind(j, sym)) = facsym(i, j)
            End do

        End do

        Deallocate (facsym)
        Deallocate (wcsym)

    End do ! sym

! NOW FAC BECOMES TRANSFORM MATRIX   CONJG(Fac) Fbc Fac = W <= diagonal form!

    Allocate (matc(n0:n1, n0:n1))
    matc = 0.0d+00

    fac = DCONJG(fac)
    matc = TRANSPOSE(fac)
    matc = MATMUL(matc(n0:n1, n0:n1), f(n0:n1, n0:n1))
    fac = DCONJG(fac)
    matc = MATMUL(matc, fac)

    if (rank == 0) then ! Process limits for output
        write (*, *) 'OFF DIAGONAL TERM OF U*FU'
    end if
    do i = n0, n1
        do j = n0, n1
            if ((i /= j) .and. (ABS(matc(i, j)) > 1.0d-10)) then
                if (rank == 0) then ! Process limits for output
                    write (*, '(2E13.5,2I3)') matc(i, j), i, j
                end if
            end if
        end do
    end do

    if (rank == 0) then ! Process limits for output
        write (*, *) 'DIAGONAL TERM OF U*FU, W AND THEIR DIFFERENCE'
    end if
    do i = n0, n1
        if (ABS(matc(i, i) - wc(i)) > 1.0d-10) then
            if (rank == 0) then ! Process limits for output
                write (*, '(4E13.5)') matc(i, i), wc(i), ABS(matc(i, i) - wc(i))
            end if
        End if
    end do
    if (rank == 0) then ! Process limits for output
        write (*, '(/)')
    end if
    deallocate (matc)

    if (rank == 0) then ! Process limits for output
        write (*, *) 'cdiag0 end'
    end if
end subroutine cdiag0
