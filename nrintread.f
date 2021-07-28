! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   SUBROUTINE nrintread

!  This part is originally writen by Dr. T. Yanai as itrf code in program package UTChem.
!  Here is modified for reading non-relativistic integrals to compute four-CASPT2
!  By M. Abe.
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
!	 real*8,      intent(in)  :: 
!	 real*8,      intent(in)  :: 

        integer                  :: ndim, intindx, ncount, redund, i, ii
        integer                  :: nrint
        integer                  :: val_integer
        integer                  :: bitsize_integer
	integer, Allocatable     :: wrtidx(:,:), idx(:,:)
	real*8, Allocatable      :: val1(:)
        character*50             :: filename

	bitsize_integer = KIND(val_integer)*8

  	filename='moint2.info.aaaa'
	nrint = 11
        open (nrint, file=filename, status='old', access='sequential', form='formatted')

        read(nrint,*) ndim, intindx, ncount, redund

	close(nrint)

!      AT PRESENT RHF ORBITALS ARE ASSUMED!

	filename ='moint2.aaaa'
        open (nrint, file=trim(filename), &
        status='old', access='sequential', form='unformatted')


	Allocate(wrtidx(intindx, ndim))
	Allocate(val1(ndim))
	Allocate(idx(4, ndim))


	Do i = 1, ncount

           Read(nrint,ERR=40,END=50) wrtidx(1:intindx,1:ndim)
           Read(nrint,ERR=40,END=50) val1(1:ndim)

           Do ii=1, ndim

              Select case(intindx)
           
           Case (1)
              idx(1,ii) = IBITS(wrtidx(1,ii),bitsize_integer*3/4,bitsize_integer/4)
              idx(2,ii) = IBITS(wrtidx(1,ii),bitsize_integer*2/4,bitsize_integer/4)
              idx(3,ii) = IBITS(wrtidx(1,ii),bitsize_integer*1/4,bitsize_integer/4)
              idx(4,ii) = IBITS(wrtidx(1,ii),bitsize_integer*0/4,bitsize_integer/4)

           Case (2)
              idx(1,ii) = IBITS(wrtidx(1,ii),bitsize_integer*1/2,bitsize_integer/2)
              idx(2,ii) = IBITS(wrtidx(1,ii),bitsize_integer*0/2,bitsize_integer/2)
              idx(3,ii) = IBITS(wrtidx(2,ii),bitsize_integer*1/2,bitsize_integer/2)
              idx(4,ii) = IBITS(wrtidx(2,ii),bitsize_integer*0/2,bitsize_integer/2)

           Case (4)
              idx(1,ii) = wrtidx(1,ii)
              idx(2,ii) = wrtidx(2,ii)
              idx(3,ii) = wrtidx(3,ii)
              idx(4,ii) = wrtidx(4,ii)

           Case default
              write(*,*) "[INPUT ERROR] @Int2_idx : out of select ( 1 / 2 / 4 )"
              stop

           end Select


              write(*,*) idx(1:intindx, ii)
              write(*,*) val1(ii)
        
           End Do

	End do

        Read(nrint) wrtidx(:,1:redund)
        Read(nrint) val1(1:redund)


	Do ii=1,redund

	Select case(intindx)

        Case (1)
        idx(1,ii) = IBITS(wrtidx(1,ii),bitsize_integer*3/4,bitsize_integer/4)
        idx(2,ii) = IBITS(wrtidx(1,ii),bitsize_integer*2/4,bitsize_integer/4)
        idx(3,ii) = IBITS(wrtidx(1,ii),bitsize_integer*1/4,bitsize_integer/4)
        idx(4,ii) = IBITS(wrtidx(1,ii),bitsize_integer*0/4,bitsize_integer/4)

        Case (2)
   	idx(1,ii) = IBITS(wrtidx(1,ii),bitsize_integer*1/2,bitsize_integer/2)
        idx(2,ii) = IBITS(wrtidx(1,ii),bitsize_integer*0/2,bitsize_integer/2)
        idx(3,ii) = IBITS(wrtidx(2,ii),bitsize_integer*1/2,bitsize_integer/2)
        idx(4,ii) = IBITS(wrtidx(2,ii),bitsize_integer*0/2,bitsize_integer/2)

        Case (4)
        idx(1,ii) = wrtidx(1,ii)
        idx(2,ii) = wrtidx(2,ii)
        idx(3,ii) = wrtidx(3,ii)
        idx(4,ii) = wrtidx(4,ii)

        Case default
        write(*,*) "[INPUT ERROR] @Int2_idx : out of select ( 1 / 2 / 4 )"
        stop

        end Select

        End Do


40 continue
50 continue

end SUBROUTINE nrintread
