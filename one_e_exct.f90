! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine one_e_exct (iidet, creat, anhi, newidet, phase)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use four_caspt2_module

   Implicit NONE

      integer, intent(in)  :: iidet, creat, anhi
      integer, intent(out) :: newidet, phase

      integer :: ia, ib, nbitsa, indu, indv
      integer :: j0, j, i, i0, i1
      integer :: k0, l0, ii, jj, kk, ll



      ia=0
      ib=0
      phase=0
      newidet=0

!     Index for only CASCI Active space

!      indu = creat - ninact
!      indv = anhi  - ninact
      indu = creat
      indv = anhi 

      if((indv==indu).and.(btest(iidet,indv-1).eqv..true.)) then

         newidet = iidet
         phase = 0

      elseif((btest(iidet,indv-1).eqv..true.).and.(btest(iidet,indu-1).eqv..false.)) then

         newidet = iidet-2**(indv-1)+2**(indu-1)
!        calculation of phase
   
   	 do i1 = indv, nact-1
   	   ia = ia+2**i1          ! to make sequnece whose bits are all '0' bellow the vth bit
   	 end do
   
   	 do i1 = indu, nact-1
   	   ib = ib+2**i1          ! to make sequnece whose bits are all '0' bellow the uth bit
   	 end do
        
	 phase = nbitsa( iand(iidet,ia) ) + nbitsa( iand(iidet-2**(indv-1),ib) )
                       ! odd => (-), even => (+)

      endif

      end subroutine one_e_exct


