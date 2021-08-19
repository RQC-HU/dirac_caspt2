! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE casci_ty(totsym)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        use four_caspt2_module

        Implicit NONE

        integer, intent(in) :: totsym
        integer :: nbitsa, comb
        integer :: j0, j, i, i0, i1
        integer :: k0, l0, ii, jj, kk, ll, irec, cimat
        real*8 :: thresd

        complex*16, allocatable :: mat(:,:)
        real*8, allocatable     :: ecas(:)
        logical                 :: cutoff
        character*20            :: filename

        ndet = comb(nact, nelec)
        write(*,*)'ndet',ndet

        Call casdet_ty(totsym)

        Allocate (mat(ndet, ndet)); Call memplus(KIND(mat),SIZE(mat),2)

        Call casmat(mat)

        Allocate (ecas(ndet))
        ecas = 0.0d+00
        thresd = 1.0d-15
        cutoff = .FALSE.

        Call cdiag(mat, ndet, ndet, ecas, thresd, cutoff)

! Print out CI matrix!

       write(*,*) 'debug1'

        cimat = 10
        filename = 'CIMAT'
        open(10, file='CIMAT', status='unknown', form='unformatted')
        write(10) ndet
        write(10) idet(1:ndet)
        write(10) ecas(1:ndet)
!        write(10) mat(1:ndet,1:ndet)
        close(10)

! Print out C1 matrix!

! Print out CI matrix!


        write(*,*) 'debug2'

        cimat = 10
        filename = 'CIMAT1'
        open(10, file='CIMAT1', status='unknown', form='unformatted')
        write(10) ndet
        write(10) idet(1:ndet)
        write(10) ecas(1:ndet)
        write(10) mat(1:ndet,1:ndet)
        close(10)

! Print out C1 matrix!

        write(*,*) 'debug3'

        Allocate (cir(ndet, selectroot:selectroot)); Call memplus(KIND(cir),SIZE(cir),1)
        Allocate (cii(ndet, selectroot:selectroot)); Call memplus(KIND(cii),SIZE(cii),1)
        Allocate (eigen(nroot))                    ; Call memplus(KIND(eigen),SIZE(eigen),1)

        eigen(:)=0.0d+00
        cir(:,:)=0.0d+00
        cii(:,:)=0.0d+00

        eigen(1:nroot) = ecas(1:nroot) + ecore
        cir(1:ndet, selectroot) =  DBLE(mat(1:ndet, selectroot))
        cii(1:ndet, selectroot) = DIMAG(mat(1:ndet, selectroot))

        Deallocate (ecas)

        write(*,*) 'debug4'

        write(*,'("CASCI ENERGY FOR ",I2," STATE")')totsym
        Do irec = 1, nroot
           write(*,'(I4,F30.15)')irec, eigen(irec)
        End do

        do j = 1, ndet
           if(ABS(DIMAG(mat(j, selectroot))) > thres ) then
              realcvec = .false.
           end if
        end do


        do irec = 1, nroot
           write(*,'("Root = ",I4)')irec
          do j = 1, ndet
             if((ABS(mat(j, irec))**2) > 1.0d-02 ) then
                i0 = idet(j)
                write(*,*)(btest(i0,j0),j0=0,nact-1)
                write(*,'(I4,2(3X,E14.7)," Weights ",E14.7)') &
                & j, mat(j,irec), &
                & ABS(mat(j, irec))**2 
             end if
          end do
        end do

        Deallocate (mat); Call memminus(KIND(mat),SIZE(mat),2)
              

 1000   end subroutine casci_ty

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   FUNCTION comb(n, m) RESULT(res)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        Implicit NONE

        integer :: n, m, i, j, res, m0

        j = 1
        
        if(n-m < m)then
           m0 = n-m
        else
           m0 = m
        endif
           
        Do i = n-m0+1, n
           j = j*i
        End do
        
        Do i = 1, m0
           j = j/i
        End do

        res = j
 1000   end function comb


