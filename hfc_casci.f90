! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   PROGRAM hfc_casci  ! Hyperfine coupling constant calculation for perpendicular term at CASCI level

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

   use four_caspt2_module

        Implicit NONE
        integer                 :: ii, jj, iatom, iq, i, j, imo, jmo, nhomo, i0, j0
        logical                 :: test, cutoff
!        real*8                  :: 
        complex*16              :: dens, hfc(-1:1)
        complex*16,allocatable  :: ci(:) , hfcmo (:,:,:,:), mat(:,:)

        character*50            :: filename

! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
! +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
!

      write(*,*)''
      write(*,*)' Hyperfine coupling constant calculation for perpendicular term '
      write(*,*)' at CASCI level written by Abe in 2019'
      write(*,*)''

      open(5,file='active.inp',form='formatted',status='old')
      read(5,'(I4)')ninact
      read(5,'(I4)')nact
      read(5,'(I4)')nsec
      read(5,'(I4)')nelec
      read(5,'(I4)')nroot
      read(5,'(I4)')selectroot
      close(5)

      nmo = ninact + nact + nsec

      write(*,*)'ninact     =' ,ninact
      write(*,*)'nact       =' ,nact
      write(*,*)'nsec       =' ,nsec
      write(*,*)'nelec      =' ,nelec
      write(*,*)'nroot      =' ,nroot
      write(*,*)'selectroot =' ,selectroot
      write(*,*)'nmo        =' ,nmo

      filename = 'r4dmoint1hfcpc'

      Allocate(hfcmo(nmo,nmo,-1:1,1:2)) 
      
      open(unit=12,file=trim(filename), status='old', form='unformatted')
      read(12)
      do iatom=1,2
         do iq=-1,1
            read(12)((hfcmo(jmo,imo,iq,iatom),jmo=1,nmo),imo=1,nmo)
         end do
      end do
      close(12)

      open(10,file='CIMAT1',form='unformatted',status='old')

      read(10) ndet
      Allocate(idet(1:ndet))
      Allocate(mat(ndet,ndet))
      read(10) idet(1:ndet)
      read(10)
      read(10) mat(1:ndet,1:ndet)
      close(10)

      Allocate(ci(ndet))
      ci = mat ( :, selectroot)

      Deallocate (mat)

      do j = 1, ndet
         if((ABS(ci(j))**2) > 1.0d-02 ) then
            i0 = idet(j)
            write(*,*)(btest(i0,j0),j0=0,nact-1)
            write(*,'(I4,2(3X,E14.7)," Weights ",E14.7)') &
            & j, ci(j), ABS(ci(j))**2
         end if
      end do


      Allocate(cir(1:ndet,selectroot:selectroot))
      Allocate(cii(1:ndet,selectroot:selectroot))

      cir(1:ndet,selectroot) = DBLE(ci(1:ndet))
      cii(1:ndet,selectroot) = DIMAG(ci(1:ndet))

      deallocate(ci)


      iroot = selectroot                       
      hfc = 0.0d+00
      iatom = 1
      nhomo = nelec + ninact
      write(*,*) 'nhomo,hfcmo(nhomo,nhomo,0,iatom)   ',nhomo,hfcmo(nhomo,nhomo,0,iatom )
      write(*,*) 'nhomo,hfcmo(nhomo,nhomo+1,0,iatom) ',nhomo,hfcmo(nhomo,nhomo+1,0,iatom )
      write(*,*) 'nhomo,hfcmo(nhomo,nhomo+1,1,iatom) ',nhomo,hfcmo(nhomo,nhomo+1,1,iatom )
      write(*,*) 'nhomo,hfcmo(nhomo+1,nhomo,1,iatom) ',nhomo,hfcmo(nhomo+1,nhomo,1,iatom )
      write(*,*) 'nhomo,hfcmo(nhomo,nhomo+1,-1,iatom)',nhomo,hfcmo(nhomo,nhomo+1,-1,iatom )
      write(*,*) 'nhomo,hfcmo(nhomo+1,nhomo,-1,iatom) ',nhomo,hfcmo(nhomo+1,nhomo,-1,iatom )

      iq = -1
      Do i = 1, nact
      Do j = 1, nact
         Call dim1_density_nondiag  (i, j, dens)
         ii = i + ninact
         if(mod(j,2)==0) jj = j - 1 + ninact
         if(mod(j,2)==1) jj = j + 1 + ninact
!         write(*,*) 'ii,jj,dens,hfcmo(ii,jj,iq,iatom )',ii,jj,dens,hfcmo(ii,jj,iq,iatom )
         hfc(iq) = hfc(iq) + dens*hfcmo(ii,jj,iq,iatom )
      End do
      End do
      write(*,*)'hfc,iq', hfc(iq),iq

      iq = 0
      Do i = 1, nact
      Do j = 1, nact
         Call dim1_density_diag     (i, j, dens)
         ii = i + ninact
         jj = j + ninact
         write(*,*) 'ii,jj,dens,hfcmo(ii,jj,iq,iatom )',ii,jj,dens,hfcmo(ii,jj,iq,iatom )
         hfc(iq) = hfc(iq) + dens*hfcmo(ii,jj,iq,iatom )
      End do
      End do
      write(*,*)'hfc,iq', hfc(iq),iq

      deallocate (cir)   
      deallocate (cii)   
      deallocate (idet)  
      deallocate (hfcmo)

      END program hfc_casci



