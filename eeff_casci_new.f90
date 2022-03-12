! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

PROGRAM eeff_casci  ! Hyperfine coupling constant calculation for perpendicular term at CASCI level

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    use four_caspt2_module

    Implicit NONE
    integer                 :: ii, jj, iq, i, j, imo, jmo, nhomo
    logical                 :: test, cutoff
!        real*8                  ::
    complex*16              :: dens, eeff
    complex*16, allocatable  :: ci(:), eeffmo(:, :)

    character*50            :: filename

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!

    write (*, *) ''
    write (*, *) ' Eeff calculation'
    write (*, *) ' at CASCI level written by Abe in 2019'
    write (*, *) ''

    open (5, file='active.inp', form='formatted', status='old')
    read (5, '(I4)') ninact
    read (5, '(I4)') nact
    read (5, '(I4)') nsec
    read (5, '(I4)') nelec
    read (5, '(I4)') nroot
    read (5, '(I4)') selectroot
    close (5)

    nmo = ninact + nact + nsec

    write (*, *) 'ninact     =', ninact
    write (*, *) 'nact       =', nact
    write (*, *) 'nsec       =', nsec
    write (*, *) 'nelec      =', nelec
    write (*, *) 'nroot      =', nroot
    write (*, *) 'selectroot =', selectroot
    write (*, *) 'nmo        =', nmo

    filename = 'r4dmoint1relp2'

    Allocate (eeffmo(nmo, nmo))

    open (unit=12, file=trim(filename), status='old', form='unformatted')
    read (12)
    read (12) ((eeffmo(jmo, imo), jmo=1, nmo), imo=1, nmo)
    close (12)

    open (10, file='CIMAT', form='unformatted', status='old')

    read (10) ndet
    Allocate (idet(1:ndet))
    read (10) idet(1:ndet)

    close (10)

    Allocate (ci(1:ndet))
    ci = 0.0d+00

    open (10, file='NEWCICOEFF', form='unformatted', status='old')
    read (10) ci(1:ndet)
    close (10)

    Do crei = 1, nact
        Do anhj = 1, nact

            dens = 0.0d+00

            Do i0 = 1, ndet
                i = idet(i0)

                call one_e_exct(i, crei, anhj, newidet, phase)
                if (newidet == 0) goto 10
                i = newidet
                phasenew = phase
                j0 = 0

                do i1 = 1, ndet
                    j = idet(i1)
                    if (j == i) then
                        j0 = i1
                        goto 1
                    end if
                end do
1               continue

                if (j0 == 0) then
                    go to 10
                end if

                if (mod(phasenew, 2) == 0) then
                    dens = dens - ci(j)*DCONJG(ci(i))
                else
                    dens = dens - ci(j)*DCONJG(ci(i))
                end if

                ii = i + ninact
                jj = j + ninact
                if (ABS(dens) > 1.0d-15) write (*, *) 'ii,jj,dens,eeffmo(ii,jj )', ii, jj, dens, eeffmo(ii, jj)
                eeff = eeff + dens*eeffmo(ii, jj)

                iroot = selectroot
                eeff = 0.0d+00
                nhomo = nelec + ninact
                write (*, *) 'nhomo,eeffmo(nhomo,nhomo)   ', nhomo, eeffmo(nhomo, nhomo)
                write (*, *) 'nhomo,eeffmo(nhomo,nhomo+1) ', nhomo, eeffmo(nhomo, nhomo + 1)

                Do i = 1, nact
                    Do j = 1, nact
                        Call dim1_density_diag(i, j, dens)
                    End do
                End do
                write (*, *) 'eeff', eeff

                deallocate (ci)
                deallocate (idet)
                deallocate (eeffmo)

                END program eeff_casci
