! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE rdiag(sr, dimn, dimm, w, cutoff_threshold)
! diagonalization of real symmetric matrix
! and remove linear dependency for any S matrix
! (if cutoff_threshold = 0.0d+00, cutoff process is not performed.)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_error, only: stop_with_errorcode
    use module_global_variables

    Implicit NONE

    integer, intent(in) :: dimn
    real(8), intent(in)  :: cutoff_threshold

    real(8), intent(inout)  :: sr(dimn, dimn)

    integer, intent(out) :: dimm
    real(8), intent(out)  ::  w(dimn)

    integer :: info, lda, lwork
    character :: jobz*1, uplo*1
    real(8), allocatable  ::  work(:)

    ! Prepare for diagonalization
    w(:) = 0.0d+00
    jobz = 'V' ! calculate eigenvectors
    uplo = 'U' ! calculate upper triangle matrix
    lda = max(1, dimn)
    lwork = max(1, 3*dimn - 1)
    allocate (work(lwork))
    work(:) = 0.0d+00

    ! Symmetric matrix diagonalization by LAPACK routine DSYEV
    ! DSYEV: https://netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html
    ! sr is overwritten by eigenvectors
    ! w is eigenvalues
    ! info is the error code (info = 0 means no error)
    call dsyev(jobz, uplo, dimn, sr, lda, w, work, lwork, info)

    deallocate (work)

    ! Error check
    if (info /= 0) then
        if (rank == 0) print *, 'error in diagonalization, info = ', info
        call stop_with_errorcode(info)
    end if

    ! If cutoff_threshold == 0.0d+00, cutoff process is not performed.
    if (cutoff_threshold == 0.0d+00) then
        dimm = dimn
    else ! cutoff process is performed
        if (debug .and. rank == 0) print *, 'cut off threshold is ', cutoff_threshold
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

    use module_error, only: stop_with_errorcode
    use module_global_variables

    Implicit NONE

    integer, intent(in) :: dimn
    real(8), intent(in)  :: cutoff_threshold

    complex*16, intent(inout):: c(dimn, dimn)

    integer, intent(out) :: dimm
    real(8), intent(out)  :: w(dimn)

    integer :: info, lda, lwork
    character :: jobz*1, uplo*1

    complex*16, allocatable  ::  work(:)
    real(8), allocatable      ::  rwork(:)

    ! Prepare for diagonalization
    w(:) = 0.0d+00
    jobz = 'V' ! calculate eigenvectors
    uplo = 'U' ! calculate upper triangle matrix
    lda = max(1, dimn)
    lwork = max(1, 3*dimn - 1)
    allocate (work(lwork))
    allocate (rwork(3*dimn - 2))
    work(:) = 0.0d+00
    rwork(:) = 0.0d+00

    ! Hermite matrix diagonalization by LAPACK routine ZHEEV
    ! ZHEEV: https://netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaf23fb5b3ae38072ef4890ba43d5cfea2.html
    ! c is overwritten by eigenvectors
    ! w is eigenvalues
    ! info is the error code (info = 0 means no error)
    Call ZHEEV(JOBZ, UPLO, dimn, c, LDA, W, WORK, LWORK, RWORK, INFO)

    deallocate (work)
    deallocate (rwork)

    ! Error check
    if (info /= 0) then
        if (rank == 0) print *, 'error in diagonalization, info = ', info
        call stop_with_errorcode(info)
    end if

    ! If cutoff_threshold == 0.0d+00, cutoff process is not performed.
    if (cutoff_threshold == 0.0d+00) then
        dimm = dimn
    else ! cutoff process is performed
        if (debug .and. rank == 0) print *, 'cut off threshold is ', cutoff_threshold
        dimm = count(w(1:dimn) >= cutoff_threshold)
    end if

end subroutine cdiag

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE rdiag0(n, n0, n1, fa, w)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables

    Implicit NONE
    integer, intent(in)     ::  n, n0, n1

    real(8), intent(out)     ::  fa(n0:n1, n0:n1)
    real(8), intent(out)     ::  w(n0:n1)

    real(8)                 ::  cutoff_threshold
    integer                 ::  j, i, dimn, dummy, ncount(nsymrpa)
    integer                 ::  isym
    integer                 ::  ind(n, nsymrpa)

    real(8), allocatable     ::  mat(:, :), fasym(:, :)
    real(8)                  ::  wsym(n)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

!  DAIAGONALIZATION OF A COMPLEX HERMITIAN MATRIX
    if (debug .and. rank == 0) print *, 'rdiag0 start'
    w(:) = 0.0d+00

    fa(n0:n1, n0:n1) = 0.0d+00

!        DAIGONALIZE IN EACH SYMMETRY
!        FIRST, THE ELEMENTS ARE REORDERED BY SYMMETRIC ORDER

! SET NCOUNT(SYM) : DIMENSION OF EACH SYMMETRY

    ncount(:) = 0

    Do i = n0, n1
        isym = irpamo(i)
        ncount(isym) = ncount(isym) + 1
        ind(ncount(isym), isym) = i
    End do

    Do isym = 1, nsymrpa

        dimn = ncount(isym)
        Allocate (fasym(dimn, dimn))
        fasym(1:dimn, 1:dimn) = fock_real(ind(1:dimn, isym), ind(1:dimn, isym))

        cutoff_threshold = 0.0d+00 ! No cutoff
        Call rdiag(fasym, dimn, dummy, wsym, cutoff_threshold)

        w(ind(1:dimn, isym)) = wsym(1:dimn)
        Do j = 1, dimn
            Do i = 1, dimn
                fa(ind(i, isym), ind(j, isym)) = fasym(i, j)
            End do
        End do

        Deallocate (fasym)

    End do

! NOW FA BECOMES TRANSFORM MATRIX   CONJG(Fa) Fbc Fa = W <= diagonal form!

    Allocate (mat(n0:n1, n0:n1))
    mat = 0.0d+00

    mat = TRANSPOSE(fa)
    mat = MATMUL(mat, fock_real(n0:n1, n0:n1))
    mat = MATMUL(mat, fa)

    if (debug .and. rank == 0) then
        print *, 'OFF DIAGONAL TERM OF U*FU (print only abs(diff) > 1.0d-08)'
        do j = n0, n1
            do i = n0, n1
                if (i /= j .and. (ABS(mat(i, j)) > 1.0d-08)) then
                    print '(E13.5,2I8)', mat(i, j), i, j
                end if
            end do
        end do

        print *, 'DIAGONAL TERM OF U*FU, W AND THEIR DIFFERENCE (print only abs(diff) > 1.0d-08)'
        do i = n0, n1
            if (ABS(mat(i, i) - w(i)) > 1.0d-08) then
                print '(3E13.5)', mat(i, i), w(i), ABS(mat(i, i) - w(i))
            end if
        end do
    end if
    deallocate (mat)

    if (debug .and. rank == 0) print *, 'rdiag0 end'
end subroutine rdiag0

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE cdiag0(n, n0, n1, fac, wc)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use module_global_variables

    Implicit NONE
    integer, intent(in)     ::  n, n0, n1

    complex*16, intent(out) ::  fac(n0:n1, n0:n1)
    real(8), intent(out)     ::  wc(n0:n1)

    real(8)                 ::  cutoff_threshold
    logical                 ::  fi
    integer                 ::  j, i, dimn, dummy, ncount(nsymrpa)
    integer                 ::  sym, isym
    integer                 ::  ind(n, nsymrpa)

    complex*16, allocatable ::  matc(:, :), facsym(:, :), facsymo(:, :)
    real(8), allocatable      ::  wcsym(:)

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

    if (debug .and. rank == 0) print *, 'fi', fi

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
                        print '("sym=",3I8,2E20.10)', sym, i, j, facsymo(i, j)
                    End if
                End do
            End do
            Do i = 1, dimn
                If (ABS(facsymo(i, i) - wcsym(i)) > 1.0d-10) then
                    print '("sym=",2I8,3E20.10)', sym, i, facsymo(i, i), wcsym(i)
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
    if (debug .and. rank == 0) then
        print *, 'OFF DIAGONAL TERM OF U*FU (print only abs(diff) > 1.0d-08)'
        do j = n0, n1
            do i = n0, n1
                if ((i /= j) .and. (ABS(matc(i, j)) > 1.0d-08)) then
                    print '(2E13.5,2I8)', matc(i, j), i, j
                end if
            end do
        end do
        print *, 'DIAGONAL TERM OF U*FU, W AND THEIR DIFFERENCE (print only abs(diff) > 1.0d-08)'
        do i = n0, n1
            if (ABS(matc(i, i) - wc(i)) > 1.0d-08) then
                print '(4E13.5)', matc(i, i), wc(i), ABS(matc(i, i) - wc(i))
            End if
        end do
    end if

    deallocate (matc)

    if (debug .and. rank == 0) print *, 'cdiag0 end'
end subroutine cdiag0
