! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE rdiag(sr, dimn, dimm, w, cutoff_threshold)
! diagonalization of real symmetric matrix
! and remove linear dependency for any S matrix
! (if cutoff_threshold = 0.0d+00, cutoff process is not performed.)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE

    integer, intent(in) :: dimn
    real*8, intent(in)  :: cutoff_threshold

    real*8, intent(inout)  :: sr(dimn, dimn)

    integer, intent(out) :: dimm
    real*8, intent(out)  ::  w(dimn)

    integer :: info, lda, lwork
    character :: jobz*1, uplo*1
    real*8, allocatable  ::  work(:)

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

    if (info /= 0) then
        if (rank == 0) print *, 'error in diagonalization, info = ', info
        return
    end if

    ! If cutoff_threshold == 0.0d+00, cutoff process is not performed.
    if (cutoff_threshold == 0.0d+00) then
        dimm = dimn
    else
        if (rank == 0) print *, 'cut off threshold is ', cutoff_threshold
        dimm = count(w(1:dimn) >= cutoff_threshold)
    end if

end subroutine rdiag

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE cdiag(c, dimn, dimm, w, cutoff_threshold)
! diagonalization of complex symmetric matrix
! and remove linear dependency for any S matrix
! (if cutoff_threshold = 0.0d+00, cutoff process is not performed.)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE

    integer, intent(in) :: dimn
    real*8, intent(in)  :: cutoff_threshold

    complex*16, intent(inout):: c(dimn, dimn)

    integer, intent(out) :: dimm
    real*8, intent(out)  :: w(dimn)

    integer :: info, lda, lwork
    character :: jobz*1, uplo*1

    complex*16, allocatable  ::  work(:)
    real*8, allocatable      ::  rwork(:)

    if (rank == 0) print *, 'Enter cdiagonal part'
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

    if (rank == 0) print *, 'Finish zheev info = ', info
    if (info /= 0) then
        if (rank == 0) print *, 'error in diagonalization, info = ', info
        return
    end if

    ! If cutoff_threshold == 0.0d+00, cutoff process is not performed.
    if (cutoff_threshold == 0.0d+00) then
        dimm = dimn
    else
        if (rank == 0) print *, 'cut off threshold is ', cutoff_threshold
        dimm = count(w(1:dimn) >= cutoff_threshold)
    end if

    if (rank == 0) print *, "end cdiag"
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

    real(8)                 ::  cutoff_threshold
    integer                 ::  j, i, dimn, dummy, ncount(nsymrpa)
    integer                 ::  sym, isym
    integer                 ::  ind(n, nsymrpa)

    real*8, allocatable     ::  mat(:, :), fasym(:, :)
    real*8                  ::  wsym(n)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!  DAIAGONALIZATION OF A COMPLEX HERMITIAN MATRIX
    if (rank == 0) print *, 'rdiag0 start'
    w(:) = 0.0d+00

    fa(n0:n1, n0:n1) = 0.0d+00

!        DAIGONALIZE IN EACH SYMMETRY
!        FIRST, THE ELEMENTS ARE REORDERED BY SYMMETRIC ORDER

! SET NCOUNT(SYM) : DIMENSION OF EACH SYMMETRY

    ncount(:) = 0

    if (rank == 0) print *, 'nsymrpa', nsymrpa

    Do i = n0, n1
        isym = irpamo(i)
        ncount(isym) = ncount(isym) + 1
        ind(ncount(isym), isym) = i
    End do

    if (rank == 0) print *, 'sym,ncount(sym)', (ncount(sym), sym=1, nsymrpa)
    Do sym = 1, nsymrpa

        dimn = ncount(sym)
        Allocate (fasym(dimn, dimn))
        fasym(1:dimn, 1:dimn) = fock_real(ind(1:dimn, sym), ind(1:dimn, sym))

        cutoff_threshold = 0.0d+00 ! No cutoff
        Call rdiag(fasym, dimn, dummy, wsym, cutoff_threshold)
!      _________________________________________________________

        w(ind(1:dimn, sym)) = wsym(1:dimn)
        Do j = 1, dimn
            Do i = 1, dimn
                fa(ind(i, sym), ind(j, sym)) = fasym(i, j)
            End do
        End do

        Deallocate (fasym)

    End do ! sym

! NOW FA BECOMES TRANSFORM MATRIX   CONJG(Fa) Fbc Fa = W <= diagonal form!

    Allocate (mat(n0:n1, n0:n1))
    mat = 0.0d+00

    mat = TRANSPOSE(fa)
    mat = MATMUL(mat, fock_real(n0:n1, n0:n1))
    mat = MATMUL(mat, fa)

    if (rank == 0) then
        print *, 'OFF DIAGONAL TERM OF U*FU'
        do j = n0, n1
            do i = n0, n1
                if (i /= j .and. (ABS(mat(i, j)) > 1.0d-10)) then
                    print '(2E13.5,2I3)', mat(i, j), i, j
                end if
            end do
        end do

        print *, 'DIAGONAL TERM OF U*FU, W AND THEIR DIFFERENCE'
        do i = n0, n1
            print '(4E13.5)', mat(i, i), w(i), ABS(mat(i, i) - w(i))
        end do
    end if
    deallocate (mat)

    if (rank == 0) print *, 'rdiag0 end'
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

    real(8)                 ::  cutoff_threshold
    logical                 ::  fi
    integer                 ::  j, i, dimn, dummy, ncount(nsymrpa)
    integer                 ::  sym, isym
    integer                 ::  ind(n, nsymrpa)

    complex*16, allocatable ::  matc(:, :), facsym(:, :), facsymo(:, :)
    real*8, allocatable      ::  wcsym(:)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!  DAIAGONALIZATION OF A COMPLEX HERMITIAN MATRIX
    if (rank == 0) then
        print *, 'cdiag0 start'
        print *, 'nsymrpa', nsymrpa
    end if

    wc(:) = 0.0d+00
    if (count(abs(dimag(fock_cmplx(n0:n1, n0:n1))) > 1.0d-10) > 0) then
        fi = .TRUE.
    else
        fi = .FALSE.
    end if

    if (rank == 0) print *, 'fi', fi

    fac(n0:n1, n0:n1) = 0.0d+00

!        DAIGONALIZE IN EACH SYMMETRY
!        FIRST, THE ELEMENTS ARE REORDERED BY SYMMETRIC ORDER

! SET NCOUNT(SYM) : DIMENSION OF EACH SYMMETRY

    ncount(:) = 0
    ind(:, :) = 0

    do i = n0, n1
        isym = irpamo(i)
        ncount(isym) = ncount(isym) + 1
        ind(ncount(isym), isym) = i
    end do

    Do sym = 1, nsymrpa

        dimn = ncount(sym)
        Allocate (facsym(dimn, dimn))
        facsym = 0.0d+00

        Do j = 1, dimn
            Do i = j, dimn
                facsym(i, j) = fock_cmplx(ind(i, sym), ind(j, sym))
                facsym(j, i) = DCONJG(facsym(i, j)) ! HERMITE
            End do
        End do

        Allocate (facsymo(dimn, dimn))
        facsymo = facsym
        Allocate (wcsym(dimn))
        wcsym = 0.0d+00

        cutoff_threshold = 0.0d+00  ! No cutoff
        Call cdiag(facsym, dimn, dummy, wcsym, cutoff_threshold)
!      _________________________________________________________

        facsym = DCONJG(facsym)
        facsymo = MATMUL(TRANSPOSE(facsym), facsymo)
        facsym = DCONJG(facsym)
        facsymo = MATMUL(facsymo, facsym)

        ! Check facsymo
        if (rank == 0) then
            Do j = 1, dimn
                Do i = 1, dimn
                    If (i /= j .and. ABS(facsymo(i, j)) > 1.0d-10) then
                        print '("sym=",3I4,2E20.10)', sym, i, j, facsymo(i, j)
                    End if
                End do
            End do
            Do i = 1, dimn
                If (ABS(facsymo(i, i) - wcsym(i)) > 1.0d-10) then
                    print '("sym=",2I4,3E20.10)', sym, i, facsymo(i, i), wcsym(i)
                End if
            End do
        end if

        Deallocate (facsymo)

        wc(ind(1:dimn, sym)) = wcsym(1:dimn)
        Do j = 1, dimn
            Do i = 1, dimn
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
    matc = MATMUL(matc(n0:n1, n0:n1), fock_cmplx(n0:n1, n0:n1))
    fac = DCONJG(fac)
    matc = MATMUL(matc, fac)

    ! Check U*FU
    if (rank == 0) then
        print *, 'OFF DIAGONAL TERM OF U*FU'
        do j = n0, n1
            do i = n0, n1
                if ((i /= j) .and. (ABS(matc(i, j)) > 1.0d-10)) then
                    print '(2E13.5,2I3)', matc(i, j), i, j
                end if
            end do
        end do
        print *, 'DIAGONAL TERM OF U*FU, W AND THEIR DIFFERENCE'
        do i = n0, n1
            if (ABS(matc(i, i) - wc(i)) > 1.0d-10) then
                print '(4E13.5)', matc(i, i), wc(i), ABS(matc(i, i) - wc(i))
            End if
        end do
    end if

    deallocate (matc)

    if (rank == 0) print *, 'cdiag0 end'
end subroutine cdiag0
