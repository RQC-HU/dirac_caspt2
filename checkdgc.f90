! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

SUBROUTINE checkdgc(n, old, tra, w)

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
    integer, intent(in)        :: n
    complex*16, intent(in)     :: old(n, n)
    complex*16, intent(in)     :: tra(n, n)
    complex*16                 :: mat(n, n)
    real*8, intent(in)         :: w(n)

    integer :: i, j

    mat = MATMUL(TRANSPOSE(DCONJG(tra)), old)
    mat = MATMUL(mat, tra)

    Do i = 1, n
        If (ABS(mat(i, i) - w(i)) > 1.0d-08) then
            if (rank == 0) then
                write (normaloutput, '(I4,4E15.7)') i, mat(i, i), mat(i, i) - w(i)
            end if
        End if
        Do j = 1, n
            If ((i /= j) .and. ABS(mat(i, j)) > 1.d-08) then
                if (rank == 0) then
                    write (normaloutput, '(2I4,2E15.7)') i, j, mat(i, j)
                end if
            End if
        End do
    End do

end subroutine checkdgc
