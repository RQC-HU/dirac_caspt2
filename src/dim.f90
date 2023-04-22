! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine dim1_density(creat1, anhi1, sr, si)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
    integer, intent(in) :: creat1, anhi1
    real*8, intent(out) :: sr, si

    integer :: j0, i, i0
    complex*16 :: cmplxcii, cmplxcij, cmplxs

    integer ::  newcas_idx, phase, phasenew

! calculation of <0|Ec1a1|0>

    sr = 0.0d+00
    si = 0.0d+00

    do i0 = 1, ndet

        i = cas_idx(i0)

        call one_e_exct(i, creat1, anhi1, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx !励起演算子がかかった後の番号
        phasenew = phase

        j0 = cas_idx_reverse(i)

        if (j0 == 0) cycle ! Next i0

        cmplxcii = DCMPLX(cir(i0, iroot), cii(i0, iroot))
        cmplxcij = DCMPLX(cir(j0, iroot), cii(j0, iroot))

        cmplxs = cmplxcii*DCONJG(cmplxcij)

        if (mod(phasenew, 2) == 0) then
            sr = sr + REAL(cmplxs, 8)
            si = si + DIMAG(cmplxs)
        else
            sr = sr - REAL(cmplxs, 8)
            si = si - DIMAG(cmplxs)
        end if

    end do

end subroutine dim1_density

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine dim1_density_nondiag(creat1, anhi1, s)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
    integer, intent(in)     :: creat1, anhi1
    complex*16, intent(out) :: s

    integer :: j0, i, i0
    complex*16 :: cmplxcii, cmplxcij, cmplxs

    integer ::  newcas_idx, phase, phasenew

    s = 0.0d+00
    do i0 = 1, ndet
        i = cas_idx(i0)
        call one_e_exct(i, creat1, anhi1, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phase
        j0 = cas_idx_reverse(i)

        if (j0 == 0) cycle ! Next i0

        cmplxcii = DCMPLX(cir(i0, iroot), cii(i0, iroot))
        cmplxcij = DCMPLX(cir(j0, iroot), cii(j0, iroot))

        cmplxs = DCONJG(cmplxcii)*DCONJG(cmplxcij)    ! This part is unique for Aperp term.

        if (mod(phasenew, 2) == 0) then
            s = s + cmplxs
        else
            s = s - cmplxs
        end if

    end do

end subroutine dim1_density_nondiag

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine dim1_density_diag(creat1, anhi1, s)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
    integer, intent(in)     :: creat1, anhi1
    complex*16, intent(out) :: s

    integer :: j0, i, i0
    complex*16 :: cmplxcii, cmplxcij, cmplxs

    integer ::  newcas_idx, phase, phasenew

    s = 0.0d+00
    do i0 = 1, ndet
        i = cas_idx(i0)
        call one_e_exct(i, creat1, anhi1, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phase
        j0 = cas_idx_reverse(i)

        if (j0 == 0) cycle ! Next i0

        cmplxcii = DCMPLX(cir(i0, iroot), cii(i0, iroot))
        cmplxcij = DCMPLX(cir(j0, iroot), cii(j0, iroot))

        cmplxs = cmplxcii*DCONJG(cmplxcij)

        if (mod(phasenew, 2) == 0) then
            s = s + cmplxs
        else
            s = s - cmplxs
        end if

    end do

end subroutine dim1_density_diag

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine dim2_density(creat1, anhi1, creat2, anhi2, sr, si)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
    integer, intent(in) :: creat2, anhi2, anhi1, creat1
    real*8, intent(out) :: sr, si
    integer :: newcas_idx, phase, phasenew
    integer :: j0, i, i0
    complex*16 :: cmplxcii, cmplxcij, cmplxs

! calculation of <0|Ec1a1Ec2a2|0>

    sr = 0.0d+00
    si = 0.0d+00

    do i0 = 1, ndet

        i = cas_idx(i0)

        call one_e_exct(i, creat2, anhi2, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phase

        call one_e_exct(i, creat1, anhi1, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phasenew + phase
        j0 = cas_idx_reverse(i)

        if (j0 == 0) cycle ! Next i0

!       Caluculation of C(i,iroot)*conjugate(C(j,iroot))

        cmplxcii = DCMPLX(cir(i0, iroot), cii(i0, iroot))
        cmplxcij = DCMPLX(cir(j0, iroot), cii(j0, iroot))

        cmplxs = cmplxcii*DCONJG(cmplxcij)

        if (mod(phasenew, 2) == 0) then
            sr = sr + DBLE(cmplxs)
            si = si + DIMAG(cmplxs)
        else
            sr = sr - DBLE(cmplxs)
            si = si - DIMAG(cmplxs)
        end if

    end do

end subroutine dim2_density

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine dim3_density(creat1, anhi1, creat2, anhi2, creat3, anhi3, sr, si)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
    integer, intent(in) :: anhi3, anhi2, anhi1, creat3, creat2, creat1
    real*8, intent(out) :: sr, si
    integer :: newcas_idx, phase, phasenew
    integer :: j0, i, i0
    complex*16 :: cmplxcii, cmplxcij, cmplxs

! calculation of <0|Ec1a1Ec2a2Ec3a3|0>

    sr = 0.0d+00
    si = 0.0d+00

    do i0 = 1, ndet

        i = cas_idx(i0)

        call one_e_exct(i, creat3, anhi3, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phase

        call one_e_exct(i, creat2, anhi2, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phasenew + phase

        call one_e_exct(i, creat1, anhi1, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phasenew + phase

        j0 = cas_idx_reverse(i)

        if (j0 == 0) cycle ! Next i0

        cmplxcii = DCMPLX(cir(i0, iroot), cii(i0, iroot))
        cmplxcij = DCMPLX(cir(j0, iroot), cii(j0, iroot))

        cmplxs = cmplxcii*DCONJG(cmplxcij)

        if (mod(phasenew, 2) == 0) then
            sr = sr + DBLE(cmplxs)
            si = si + DIMAG(cmplxs)
        else
            sr = sr - DBLE(cmplxs)
            si = si - DIMAG(cmplxs)
        end if

    end do

end subroutine dim3_density

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine dim4_density(creat1, anhi1, creat2, anhi2, creat3, anhi3, creat4, anhi4, sr, si)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
    integer, intent(in) :: anhi4, anhi3, anhi2, anhi1, creat4, creat3, creat2, creat1
    real*8, intent(out) :: sr, si
    integer :: newcas_idx, phase, phasenew
    integer :: j0, i, i0
    complex*16 :: cmplxcii, cmplxcij, cmplxs

! calculation of <0|Ec1a1Ec2a2Ec3a3Ec4a4|0>

    sr = 0.0d+00
    si = 0.0d+00

    do i0 = 1, ndet

        i = cas_idx(i0)

        call one_e_exct(i, creat4, anhi4, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phase

        call one_e_exct(i, creat3, anhi3, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phasenew + phase

        call one_e_exct(i, creat2, anhi2, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phasenew + phase

        call one_e_exct(i, creat1, anhi1, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx

        phasenew = phasenew + phase
        j0 = cas_idx_reverse(i)

        if (j0 == 0) cycle ! Next i0

        cmplxcii = DCMPLX(cir(i0, iroot), cii(i0, iroot))
        cmplxcij = DCMPLX(cir(j0, iroot), cii(j0, iroot))

        cmplxs = cmplxcii*DCONJG(cmplxcij)

        if (mod(phasenew, 2) == 0) then
            sr = sr + DBLE(cmplxs)
            si = si + DIMAG(cmplxs)
        else
            sr = sr - DBLE(cmplxs)
            si = si - DIMAG(cmplxs)
        end if

    end do

end subroutine dim4_density

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine dim1_density_R(creat1, anhi1, sr)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
    integer, intent(in) :: creat1, anhi1
    real*8, intent(out) :: sr
    integer :: j0, i, i0

    integer ::  newcas_idx, phase, phasenew

! calculation of <0|Ec1a1|0> for REAL

    sr = 0.0d+00

    do i0 = 1, ndet

        i = cas_idx(i0)

        call one_e_exct(i, creat1, anhi1, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0

        i = newcas_idx
        phasenew = phase
        j0 = cas_idx_reverse(i)

        if (j0 == 0) cycle ! Next i0

        if (mod(phasenew, 2) == 0) then
            sr = sr + cir(i0, iroot)*cir(j0, iroot)
        else
            sr = sr - cir(i0, iroot)*cir(j0, iroot)
        end if

    end do

end subroutine dim1_density_R

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine dim2_density_R(creat1, anhi1, creat2, anhi2, sr)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
    integer, intent(in) :: creat2, anhi2, anhi1, creat1
    real*8, intent(out) :: sr
    integer :: newcas_idx, phase, phasenew
    integer :: j0, i, i0

! calculation of <0|Ec1a1Ec2a2|0>

    sr = 0.0d+00

    do i0 = 1, ndet

        i = cas_idx(i0)

        call one_e_exct(i, creat2, anhi2, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phase

        call one_e_exct(i, creat1, anhi1, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phasenew + phase

        j0 = cas_idx_reverse(i)

        if (j0 == 0) cycle ! Next i0

!       Caluculation of C(i,iroot)*conjugate(C(j,iroot))

        if (mod(phasenew, 2) == 0) then
            sr = sr + cir(i0, iroot)*cir(j0, iroot)
        else
            sr = sr - cir(i0, iroot)*cir(j0, iroot)
        end if

    end do

end subroutine dim2_density_R

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine dim3_density_R(creat1, anhi1, creat2, anhi2, creat3, anhi3, sr)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
    integer, intent(in) :: anhi3, anhi2, anhi1, creat3, creat2, creat1
    real*8, intent(out) :: sr
    integer :: newcas_idx, phase, phasenew
    integer :: j0, i, i0

! calculation of <0|Ec1a1Ec2a2Ec3a3|0>

    sr = 0.0d+00

    do i0 = 1, ndet

        i = cas_idx(i0)

        call one_e_exct(i, creat3, anhi3, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phase

        call one_e_exct(i, creat2, anhi2, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phasenew + phase

        call one_e_exct(i, creat1, anhi1, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phasenew + phase

        j0 = cas_idx_reverse(i)

        if (j0 == 0) cycle ! Next i0

        if (mod(phasenew, 2) == 0) then
            sr = sr + cir(i0, iroot)*cir(j0, iroot)
        else
            sr = sr - cir(i0, iroot)*cir(j0, iroot)
        end if

    end do

end subroutine dim3_density_R

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine dim4_density_R(creat1, anhi1, creat2, anhi2, creat3, anhi3, creat4, anhi4, sr)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use four_caspt2_module

    Implicit NONE
    integer, intent(in) :: anhi4, anhi3, anhi2, anhi1, creat4, creat3, creat2, creat1
    real*8, intent(out) :: sr
    integer :: newcas_idx, phase, phasenew
    integer :: j0, i, i0

! calculation of <0|Ec1a1Ec2a2Ec3a3Ec4a4|0>

    sr = 0.0d+00

    do i0 = 1, ndet

        i = cas_idx(i0)

        call one_e_exct(i, creat4, anhi4, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phase

        call one_e_exct(i, creat3, anhi3, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phasenew + phase

        call one_e_exct(i, creat2, anhi2, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx
        phasenew = phasenew + phase

        call one_e_exct(i, creat1, anhi1, newcas_idx, phase)
        if (newcas_idx == 0) cycle ! Next i0
        i = newcas_idx

        phasenew = phasenew + phase
        j0 = cas_idx_reverse(i)

        if (j0 == 0) cycle ! Next i0

        if (mod(phasenew, 2) == 0) then
            sr = sr + cir(i0, iroot)*cir(j0, iroot)
        else
            sr = sr - cir(i0, iroot)*cir(j0, iroot)
        end if

    end do

end subroutine dim4_density_R
