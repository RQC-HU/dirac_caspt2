! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE readorb_enesym (filename) ! orbital energies in MRCONEE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        use four_caspt2_module

        Implicit NONE

        integer :: mrconee
        character*50,intent(in) :: filename
        integer :: j0, j, i, i0, i1, m
        integer :: k0, l0, ii, jj, kk, ll, nmomax

!iwamuro modify
        integer :: DS(16,16), SD(16,16),  hnsym

        integer, allocatable :: dammo(:)

        real*8 :: w
!        logical :: breit
        logical :: breit

!       Write(UT_sys_ftmp) NMO,BREIT,ECORE
!       Write(UT_sys_ftmp) NSYMRP,(REPN(IRP),IRP=1,NSYMRP)
!       Write(UT_sys_ftmp) NSYMRPA,(REPNA(IRP),IRP=1,NSYMRPA*2)
!       Write(UT_sys_ftmp) ((MULTB(I,J),I=1,2*NSYMRPA),J=1,2*NSYMRPA)
!       Write(UT_sys_ftmp) (IRPMO(IMO),IRPAMO(IMO),ORBMO(IMO),IMO=1,NMO)
!       Write(UT_sys_ftmp) ((ONER(IMO,JMO),ONEI(IMO,JMO),JMO=1,NMO),IMO=1,NMO)
!       Write(UT_sys_ftmp) RAS(1), RAS(2), RAS(3)

        mrconee=10
        write(*,*)filename
        open( mrconee, file=trim(filename),form ='unformatted', status='old', err=10)
        write(*,*)'come1'

!        read(mrconee,err=11) nmo, breit, ecore
        read(mrconee,err=11) nmo, breit, ecore, nfsym, nz1, sfform, norbt
        write(*,*) nmo, breit, ecore, nfsym, nz1, sfform, norbt

        read(mrconee,err=12) nsymrp, (repn(i0), i0 = 1, nsymrp), (nelecd(i0), i0 = 1, nsymrp)
!        read(mrconee,err=12) nsymrp, (repn(i0), i0 = 1, nsymrp)
        write(*,*) nsymrp, (repn(i0), i0 = 1, nsymrp), (nelecd(i0), i0 = 1, nsymrp)

        read(mrconee,err=13) nsymrpa, (repna(i0), i0 = 1, nsymrpa*2)
!        write(*,*) nsymrpa, (repna(i0), i0 = 1, nsymrpa*2)

        read(mrconee,err=14) ((multb(i0,j0),i0=1,2*nsymrpa),j0=1,2*nsymrpa)
!        write(*,*) ((multb(i0,j0),i0=1,2*nsymrpa),j0=1,2*nsymrpa)

!         MULTB(1:16, 17:32) = 0
!         MULTB(17:32, 1:16) = 0

!--------------------------------------------------------------------------------------------------------------------------
!iwamuro modify c8h MULTB

     NSYMRP = 16
     NSYMRPA = 16
      REPNA(1) ='1e1/2g'; REPNA(2) ='2e1/2g'; REPNA(3) ='1e3/2g'; REPNA(4) ='2e3/2g'
      REPNA(5) ='1e5/2g'; REPNA(6) ='2e5/2g'; REPNA(7) ='1e7/2g'; REPNA(8) ='2e7/2g'
      REPNA(9) ='1e1/2u'; REPNA(10)='2e1/2u'; REPNA(11)='1e3/2u'; REPNA(12)='2e3/2u'
      REPNA(13)='1e5/2u'; REPNA(14)='2e5/2u'; REPNA(15)='1e7/2u'; REPNA(16)='2e7/2u'

      REPNA(17)='ag    '; REPNA(18)='bg    '; REPNA(19)='1e1g  '; REPNA(20)='2e1g  '
      REPNA(21)='1e2g  '; REPNA(22)='2e2g  '; REPNA(23)='1e3g  '; REPNA(24)='2e3g  '
      REPNA(25)='au    '; REPNA(26)='bu    '; REPNA(27)='1e1u  '; REPNA(28)='2e1u  '
      REPNA(29)='1e2u  '; REPNA(30)='2e2u  '; REPNA(31)='1e3u  '; REPNA(32)='2e3u  '

!indices 1-hnsym when singles 1-hnsym when doubles hnsym+1-nsymrp

!     SD( 1, 1)= 1; SD( 1, 2)= 2; SD( 1, 3)= 3; SD( 1, 4)= 4; SD( 1, 5)= 5; SD( 1, 6)= 6; SD( 1, 7)= 7; SD( 1, 8)= 8
!     SD( 2, 1)= 7; SD( 2, 2)= 8; SD( 2, 3)= 5; SD( 2, 4)= 6; SD( 2, 5)= 3; SD( 2, 6)= 4; SD( 2, 7)= 1; SD( 2, 8)= 2
!     SD( 3, 1)= 2; SD( 3, 2)= 3; SD( 3, 3)= 6; SD( 3, 4)= 1; SD( 3, 5)= 4; SD( 3, 6)= 7; SD( 3, 7)= 8; SD( 3, 8)= 5
!     SD( 4, 1)= 4; SD( 4, 2)= 1; SD( 4, 3)= 2; SD( 4, 4)= 5; SD( 4, 5)= 8; SD( 4, 6)= 3; SD( 4, 7)= 6; SD( 4, 8)= 7
!     SD( 5, 1)= 3; SD( 5, 2)= 6; SD( 5, 3)= 7; SD( 5, 4)= 2; SD( 5, 5)= 1; SD( 5, 6)= 8; SD( 5, 7)= 5; SD( 5, 8)= 4
!     SD( 6, 1)= 5; SD( 6, 2)= 4; SD( 6, 3)= 1; SD( 6, 4)= 8; SD( 6, 5)= 7; SD( 6, 6)= 2; SD( 6, 7)= 3; SD( 6, 8)= 6
!     SD( 7, 1)= 8; SD( 7, 2)= 5; SD( 7, 3)= 4; SD( 7, 4)= 7; SD( 7, 5)= 6; SD( 7, 6)= 5; SD( 7, 7)= 2; SD( 7, 8)= 3  !SD( 7, 6)= 1
!     SD( 8, 1)= 6; SD( 8, 2)= 7; SD( 8, 3)= 8; SD( 8, 4)= 3; SD( 8, 5)= 2; SD( 8, 6)= 1; SD( 8, 7)= 4; SD( 8, 8)= 1  !SD( 8, 6)= 5

! SD
!     MULTB( 17, 1)= 1; MULTB( 17, 2)= 2; MULTB( 17, 3)= 3; MULTB( 17, 4)= 4; MULTB( 17, 5)= 5; MULTB( 17, 6)= 6; MULTB( 17, 7)= 7; MULTB( 17, 8)= 8
!     MULTB( 18, 1)= 7; MULTB( 18, 2)= 8; MULTB( 18, 3)= 5; MULTB( 18, 4)= 6; MULTB( 18, 5)= 3; MULTB( 18, 6)= 4; MULTB( 18, 7)= 1; MULTB( 18, 8)= 2
!     MULTB( 19, 1)= 2; MULTB( 19, 2)= 3; MULTB( 19, 3)= 6; MULTB( 19, 4)= 1; MULTB( 19, 5)= 4; MULTB( 19, 6)= 7; MULTB( 19, 7)= 8; MULTB( 19, 8)= 5
!     MULTB( 20, 1)= 4; MULTB( 20, 2)= 1; MULTB( 20, 3)= 2; MULTB( 20, 4)= 5; MULTB( 20, 5)= 8; MULTB( 20, 6)= 3; MULTB( 20, 7)= 6; MULTB( 20, 8)= 7
!     MULTB( 21, 1)= 3; MULTB( 21, 2)= 6; MULTB( 21, 3)= 7; MULTB( 21, 4)= 2; MULTB( 21, 5)= 1; MULTB( 21, 6)= 8; MULTB( 21, 7)= 5; MULTB( 21, 8)= 4
!     MULTB( 22, 1)= 5; MULTB( 22, 2)= 4; MULTB( 22, 3)= 1; MULTB( 22, 4)= 8; MULTB( 22, 5)= 7; MULTB( 22, 6)= 2; MULTB( 22, 7)= 3; MULTB( 22, 8)= 6
!     MULTB( 23, 1)= 8; MULTB( 23, 2)= 5; MULTB( 23, 3)= 4; MULTB( 23, 4)= 7; MULTB( 23, 5)= 6; MULTB( 23, 6)= 1; MULTB( 23, 7)= 2; MULTB( 23, 8)= 3
!     MULTB( 24, 1)= 6; MULTB( 24, 2)= 7; MULTB( 24, 3)= 8; MULTB( 24, 4)= 3; MULTB( 24, 5)= 2; MULTB( 24, 6)= 5; MULTB( 24, 7)= 4; MULTB( 24, 8)= 1

!     MULTB( 17, 9)= 9; MULTB( 17, 10)= 10; MULTB( 17, 11)= 11; MULTB( 17, 12)= 12; MULTB( 17, 13)= 13; MULTB( 17, 14)= 14; MULTB( 17, 15)= 15; MULTB( 17, 16)= 16
!     MULTB( 18, 9)= 15; MULTB( 18, 10)= 16; MULTB( 18, 11)= 13; MULTB( 18, 12)= 14; MULTB( 18, 13)= 11; MULTB( 18, 14)= 12; MULTB( 18, 15)= 9; MULTB( 18, 16)= 10
!     MULTB( 19, 9)= 10; MULTB( 19, 10)= 11; MULTB( 19, 11)= 14; MULTB( 19, 12)= 9; MULTB( 19, 13)= 12; MULTB( 19, 14)= 15; MULTB( 19, 15)= 16; MULTB( 19, 16)= 13
!     MULTB( 20, 9)= 12; MULTB( 20, 10)= 9; MULTB( 20, 11)= 10; MULTB( 20, 12)= 13; MULTB( 20, 13)= 16; MULTB( 20, 14)= 11; MULTB( 20, 15)= 14; MULTB( 20, 16)= 15
!     MULTB( 21, 9)= 11; MULTB( 21, 10)= 14; MULTB( 21, 11)= 15; MULTB( 21, 12)= 10; MULTB( 21, 13)= 9; MULTB( 21, 14)= 16; MULTB( 21, 15)= 13; MULTB( 21, 16)= 12
!     MULTB( 22, 9)= 13; MULTB( 22, 10)= 12; MULTB( 22, 11)= 9; MULTB( 22, 12)= 16; MULTB( 22, 13)= 15; MULTB( 22, 14)= 10; MULTB( 22, 15)= 11; MULTB( 22, 16)= 14
!     MULTB( 23, 9)= 16; MULTB( 23, 10)= 13; MULTB( 23, 11)= 12; MULTB( 23, 12)= 15; MULTB( 23, 13)= 14; MULTB( 23, 14)= 9; MULTB( 23, 15)= 10; MULTB( 23, 16)= 11
!     MULTB( 24, 9)= 14; MULTB( 24, 10)= 15; MULTB( 24, 11)= 16; MULTB( 24, 12)= 11; MULTB( 24, 13)= 10; MULTB( 24, 14)= 13; MULTB( 24, 15)= 12; MULTB( 24, 16)= 9

!     MULTB( 25, 1)= 9; MULTB( 25, 2)= 10; MULTB( 25, 3)= 11; MULTB( 25, 4)= 12; MULTB( 25, 5)= 13; MULTB( 25, 6)= 14; MULTB( 25, 7)= 15; MULTB( 25, 8)= 16
!     MULTB( 26, 1)= 15; MULTB( 26, 2)= 16; MULTB( 26, 3)= 13; MULTB( 26, 4)= 14; MULTB( 26, 5)= 11; MULTB( 26, 6)= 12; MULTB( 26, 7)= 9; MULTB( 26, 8)= 10
!     MULTB( 27, 1)= 10; MULTB( 27, 2)= 11; MULTB( 27, 3)= 14; MULTB( 27, 4)= 9; MULTB( 27, 5)= 12; MULTB( 27, 6)= 15; MULTB( 27, 7)= 16; MULTB( 27, 8)= 13
!     MULTB( 28, 1)= 12; MULTB( 28, 2)= 9; MULTB( 28, 3)= 10; MULTB( 28, 4)= 13; MULTB( 28, 5)= 16; MULTB( 28, 6)= 11; MULTB( 28, 7)= 14; MULTB( 28, 8)= 15
!     MULTB( 29, 1)= 11; MULTB( 29, 2)= 14; MULTB( 29, 3)= 15; MULTB( 29, 4)= 10; MULTB( 29, 5)= 9; MULTB( 29, 6)= 16; MULTB( 29, 7)= 13; MULTB( 29, 8)= 12
!     MULTB( 30, 1)= 13; MULTB( 30, 2)= 12; MULTB( 30, 3)= 9; MULTB( 30, 4)= 16; MULTB( 30, 5)= 15; MULTB( 30, 6)= 10; MULTB( 30, 7)= 11; MULTB( 30, 8)= 14
!     MULTB( 31, 1)= 16; MULTB( 31, 2)= 13; MULTB( 31, 3)= 12; MULTB( 31, 4)= 15; MULTB( 31, 5)= 14; MULTB( 31, 6)= 9; MULTB( 31, 7)= 10; MULTB( 31, 8)= 11
!     MULTB( 32, 1)= 14; MULTB( 32, 2)= 15; MULTB( 32, 3)= 16; MULTB( 32, 4)= 11; MULTB( 32, 5)= 10; MULTB( 32, 6)= 13; MULTB( 32, 7)= 12; MULTB( 32, 8)= 9

!     MULTB( 25, 9)= 1; MULTB( 25, 10)= 2; MULTB( 25, 11)= 3; MULTB( 25, 12)= 4; MULTB( 25, 13)= 5; MULTB( 25, 14)= 6; MULTB( 25, 15)= 7; MULTB( 25, 16)= 8
!     MULTB( 26, 9)= 7; MULTB( 26, 10)= 8; MULTB( 26, 11)= 5; MULTB( 26, 12)= 6; MULTB( 26, 13)= 3; MULTB( 26, 14)= 4; MULTB( 26, 15)= 1; MULTB( 26, 16)= 2
!     MULTB( 27, 9)= 2; MULTB( 27, 10)= 3; MULTB( 27, 11)= 6; MULTB( 27, 12)= 1; MULTB( 27, 13)= 4; MULTB( 27, 14)= 7; MULTB( 27, 15)= 8; MULTB( 27, 16)= 5
!     MULTB( 28, 9)= 4; MULTB( 28, 10)= 1; MULTB( 28, 11)= 2; MULTB( 28, 12)= 5; MULTB( 28, 13)= 8; MULTB( 28, 14)= 3; MULTB( 28, 15)= 6; MULTB( 28, 16)= 7
!     MULTB( 29, 9)= 3; MULTB( 29, 10)= 6; MULTB( 29, 11)= 7; MULTB( 29, 12)= 2; MULTB( 29, 13)= 1; MULTB( 29, 14)= 8; MULTB( 29, 15)= 5; MULTB( 29, 16)= 4
!     MULTB( 30, 9)= 5; MULTB( 30, 10)= 4; MULTB( 30, 11)= 1; MULTB( 30, 12)= 8; MULTB( 30, 13)= 7; MULTB( 30, 14)= 2; MULTB( 30, 15)= 3; MULTB( 30, 16)= 6
!     MULTB( 31, 9)= 8; MULTB( 31, 10)= 5; MULTB( 31, 11)= 4; MULTB( 31, 12)= 7; MULTB( 31, 13)= 6; MULTB( 31, 14)= 1; MULTB( 31, 15)= 2; MULTB( 31, 16)= 3
!     MULTB( 32, 9)= 6; MULTB( 32, 10)= 7; MULTB( 32, 11)= 8; MULTB( 32, 12)= 3; MULTB( 32, 13)= 2; MULTB( 32, 14)= 5; MULTB( 32, 15)= 4; MULTB( 32, 16)= 1

!DS
!     MULTB( 1, 17)= 1; MULTB( 1, 18)= 7; MULTB( 1, 19)= 2; MULTB( 1, 20)= 4; MULTB( 1, 21)= 3; MULTB( 1, 22)= 5; MULTB( 1, 23)= 8; MULTB( 1, 24)= 6
!     MULTB( 2, 17)= 2; MULTB( 2, 18)= 8; MULTB( 2, 19)= 3; MULTB( 2, 20)= 1; MULTB( 2, 21)= 6; MULTB( 2, 22)= 4; MULTB( 2, 23)= 5; MULTB( 2, 24)= 7
!     MULTB( 3, 17)= 3; MULTB( 3, 18)= 5; MULTB( 3, 19)= 6; MULTB( 3, 20)= 2; MULTB( 3, 21)= 7; MULTB( 3, 22)= 1; MULTB( 3, 23)= 4; MULTB( 3, 24)= 8
!     MULTB( 4, 17)= 4; MULTB( 4, 18)= 6; MULTB( 4, 19)= 1; MULTB( 4, 20)= 5; MULTB( 4, 21)= 2; MULTB( 4, 22)= 8; MULTB( 4, 23)= 7; MULTB( 4, 24)= 3
!     MULTB( 5, 17)= 5; MULTB( 5, 18)= 3; MULTB( 5, 19)= 4; MULTB( 5, 20)= 8; MULTB( 5, 21)= 1; MULTB( 5, 22)= 7; MULTB( 5, 23)= 6; MULTB( 5, 24)= 2
!     MULTB( 6, 17)= 6; MULTB( 6, 18)= 4; MULTB( 6, 19)= 7; MULTB( 6, 20)= 3; MULTB( 6, 21)= 8; MULTB( 6, 22)= 2; MULTB( 6, 23)= 1; MULTB( 6, 24)= 5
!     MULTB( 7, 17)= 7; MULTB( 7, 18)= 1; MULTB( 7, 19)= 8; MULTB( 7, 20)= 6; MULTB( 7, 21)= 5; MULTB( 7, 22)= 3; MULTB( 7, 23)= 2; MULTB( 7, 24)= 4
!     MULTB( 8, 17)= 8; MULTB( 8, 18)= 2; MULTB( 8, 19)= 5; MULTB( 8, 20)= 7; MULTB( 8, 21)= 4; MULTB( 8, 22)= 6; MULTB( 8, 23)= 3; MULTB( 8, 24)= 1

!     MULTB( 1, 25)= 9; MULTB( 1, 26)= 15; MULTB( 1, 27)= 10; MULTB( 1, 28)= 12; MULTB( 1, 29)= 11; MULTB( 1, 30)= 13; MULTB( 1, 31)= 16; MULTB( 1, 32)= 14
!     MULTB( 2, 25)= 10; MULTB( 2, 26)= 16; MULTB( 2, 27)= 11; MULTB( 2, 28)= 9; MULTB( 2, 29)= 14; MULTB( 2, 30)= 12; MULTB( 2, 31)= 13; MULTB( 2, 32)= 15
!     MULTB( 3, 25)= 11; MULTB( 3, 26)= 13; MULTB( 3, 27)= 14; MULTB( 3, 28)= 10; MULTB( 3, 29)= 15; MULTB( 3, 30)= 9; MULTB( 3, 31)= 12; MULTB( 3, 32)= 16
!     MULTB( 4, 25)= 12; MULTB( 4, 26)= 14; MULTB( 4, 27)= 9; MULTB( 4, 28)= 13; MULTB( 4, 29)= 10; MULTB( 4, 30)= 16; MULTB( 4, 31)= 15; MULTB( 4, 32)= 11
!     MULTB( 5, 25)= 13; MULTB( 5, 26)= 11; MULTB( 5, 27)= 12; MULTB( 5, 28)= 16; MULTB( 5, 29)= 9; MULTB( 5, 30)= 15; MULTB( 5, 31)= 14; MULTB( 5, 32)= 10
!     MULTB( 6, 25)= 14; MULTB( 6, 26)= 12; MULTB( 6, 27)= 15; MULTB( 6, 28)= 11; MULTB( 6, 29)= 16; MULTB( 6, 30)= 10; MULTB( 6, 31)= 9; MULTB( 6, 32)= 13
!     MULTB( 7, 25)= 15; MULTB( 7, 26)= 9; MULTB( 7, 27)= 16; MULTB( 7, 28)= 14; MULTB( 7, 29)= 13; MULTB( 7, 30)= 11; MULTB( 7, 31)= 10; MULTB( 7, 32)= 12
!     MULTB( 8, 25)= 16; MULTB( 8, 26)= 10; MULTB( 8, 27)= 13; MULTB( 8, 28)= 15; MULTB( 8, 29)= 12; MULTB( 8, 30)= 14; MULTB( 8, 31)= 11; MULTB( 8, 32)= 9

!     MULTB( 9, 17)= 9; MULTB( 9, 18)= 15; MULTB( 9, 19)= 10; MULTB( 9, 20)= 12; MULTB( 9, 21)= 11; MULTB( 9, 22)= 13; MULTB( 9, 23)= 16; MULTB( 9, 24)= 14
!     MULTB( 10, 17)= 10; MULTB( 10, 18)= 16; MULTB( 10, 19)= 11; MULTB( 10, 20)= 9; MULTB( 10, 21)= 14; MULTB( 10, 22)= 12; MULTB( 10, 23)= 13; MULTB( 10, 24)= 15
!     MULTB( 11, 17)= 11; MULTB( 11, 18)= 13; MULTB( 11, 19)= 14; MULTB( 11, 20)= 10; MULTB( 11, 21)= 15; MULTB( 11, 22)= 9; MULTB( 11, 23)= 12; MULTB( 11, 24)= 16
!     MULTB( 12, 17)= 12; MULTB( 12, 18)= 14; MULTB( 12, 19)= 9; MULTB( 12, 20)= 13; MULTB( 12, 21)= 10; MULTB( 12, 22)= 16; MULTB( 12, 23)= 15; MULTB( 12, 24)= 11
!     MULTB( 13, 17)= 13; MULTB( 13, 18)= 11; MULTB( 13, 19)= 12; MULTB( 13, 20)= 16; MULTB( 13, 21)= 9; MULTB( 13, 22)= 15; MULTB( 13, 23)= 14; MULTB( 13, 24)= 10
!     MULTB( 14, 17)= 14; MULTB( 14, 18)= 12; MULTB( 14, 19)= 15; MULTB( 14, 20)= 11; MULTB( 14, 21)= 16; MULTB( 14, 22)= 10; MULTB( 14, 23)= 9; MULTB( 14, 24)= 13
!     MULTB( 15, 17)= 15; MULTB( 15, 18)= 9; MULTB( 15, 19)= 16; MULTB( 15, 20)= 14; MULTB( 15, 21)= 13; MULTB( 15, 22)= 11; MULTB( 15, 23)= 10; MULTB( 15, 24)= 12
!     MULTB( 16, 17)= 16; MULTB( 16, 18)= 10; MULTB( 16, 19)= 13; MULTB( 16, 20)= 15; MULTB( 16, 21)= 12; MULTB( 16, 22)= 14; MULTB( 16, 23)= 11; MULTB( 16, 24)= 9

!     MULTB( 9, 25)= 1; MULTB( 9, 26)= 7; MULTB( 9, 27)= 2; MULTB( 9, 28)= 4; MULTB( 9, 29)= 3; MULTB( 9, 30)= 5; MULTB( 9, 31)= 8; MULTB( 9, 32)= 6
!     MULTB( 10, 25)= 2; MULTB( 10, 26)= 8; MULTB( 10, 27)= 3; MULTB( 10, 28)= 1; MULTB( 10, 29)= 6; MULTB( 10, 30)= 4; MULTB( 10, 31)= 5; MULTB( 10, 32)= 7
!     MULTB( 11, 25)= 3; MULTB( 11, 26)= 5; MULTB( 11, 27)= 6; MULTB( 11, 28)= 2; MULTB( 11, 29)= 7; MULTB( 11, 30)= 1; MULTB( 11, 31)= 4; MULTB( 11, 32)= 8
!     MULTB( 12, 25)= 4; MULTB( 12, 26)= 6; MULTB( 12, 27)= 1; MULTB( 12, 28)= 5; MULTB( 12, 29)= 2; MULTB( 12, 30)= 8; MULTB( 12, 31)= 7; MULTB( 12, 32)= 3
!     MULTB( 13, 25)= 5; MULTB( 13, 26)= 3; MULTB( 13, 27)= 4; MULTB( 13, 28)= 8; MULTB( 13, 29)= 1; MULTB( 13, 30)= 7; MULTB( 13, 31)= 6; MULTB( 13, 32)= 2
!     MULTB( 14, 25)= 6; MULTB( 14, 26)= 4; MULTB( 14, 27)= 7; MULTB( 14, 28)= 3; MULTB( 14, 29)= 8; MULTB( 14, 30)= 2; MULTB( 14, 31)= 1; MULTB( 14, 32)= 5
!     MULTB( 15, 25)= 7; MULTB( 15, 26)= 1; MULTB( 15, 27)= 8; MULTB( 15, 28)= 6; MULTB( 15, 29)= 5; MULTB( 15, 30)= 3; MULTB( 15, 31)= 2; MULTB( 15, 32)= 4
!     MULTB( 16, 25)= 8; MULTB( 16, 26)= 2; MULTB( 16, 27)= 5; MULTB( 16, 28)= 7; MULTB( 16, 29)= 4; MULTB( 16, 30)= 6; MULTB( 16, 31)= 3; MULTB( 16, 32)= 1

!     Do i0 = 1, 16
!     Do j0 = 1, 16
!       MULTB(i0,j0)=MULTB(i0,j0)-16
!     End do
!     End do

!     Do i0 = 17, 32
!     Do j0 = 17, 32
!       MULTB(i0,j0)=MULTB(i0,j0)-16
!     End do
!     End do

!----------------------------------------------------------------------------------------------------------------------

      open(unit=20, file='multb_c8h.dat', action='read', &
       & form='formatted', status='old')

      Do i0 = 1,32
        read(20,*) (MULTB(i0,j0), j0 = 1,32)
      End do

      close(20)
!----------------------------------------------------------------------------------------------------------------------

        Allocate(sp(1:nmo))    ;  Call memplus(KIND(sp),SIZE(sp),1)
        sp( 1                  : ninact             )    = 1
        sp( ninact+1           : ninact+nact        )    = 2
        sp( ninact+nact+1      : ninact+nact+nsec   )    = 3
        sp( ninact+nact+nsec+1 : nmo                )    = 4

        Do i0 = 1, 2*nsymrpa
           Do j0 = 1, 2*nsymrpa
              k0 = MULTB(i0, j0)
              MULTB2(i0, k0) = j0
           Enddo
        End do

        write(*,*) 'MULTB'

        Do i0 = 1, 2*nsymrpa
           write(*,'(200I4)') (MULTB(i0, j0) ,j0 = 1, 2*nsymrpa)
        End do

        write(*,*) 'MULTB2'

        Do i0 = 1, 2*nsymrpa
           write(*,'(200I4)') (MULTB2(i0, j0) ,j0 = 1, 2*nsymrpa)
        End do

        Allocate ( irpmo (nmo)); Call memplus(KIND(irpmo ),SIZE(irpmo ),1)
        Allocate ( irpamo(nmo)); Call memplus(KIND(irpamo),SIZE(irpamo),1)
        Allocate ( orbmo (nmo)); Call memplus(KIND(orbmo ),SIZE(orbmo ),1)
        Allocate ( orb   (nmo)); Call memplus(KIND(orb   ),SIZE(orb   ),1)
        Allocate ( indmo (nmo)); Call memplus(KIND(indmo ),SIZE(indmo ),1)
        Allocate ( indmor(nmo)); Call memplus(KIND(indmor),SIZE(indmor),1)

        Allocate ( dammo (nmo)); Call memplus(KIND(dammo ),SIZE(dammo ),1)


        irpmo(:) = 0
        irpamo(:) = 0
        orbmo(:) = 0.0d+00
        orb(:) = 0.0d+00
        indmo(:) = 0

        read(mrconee,err=11) (irpmo(i0),irpamo(i0),orbmo(i0),i0 =1,nmo )

        close (mrconee)

!irpamo C8h

!        write(*,'("irpmo ",20I2)')(irpmo(i0),i0=1,20)
!        write(*,'("irpamo",20I2)')(irpamo(i0),i0=1,20)
!        write(*,'("orbmo",10F10.5)')(orbmo(i0),i0=1,10)

        irpmo(:) = irpamo(:)

        write(*,'("irpmo ",20I2)')(irpmo(i0),i0=1,nmo)
        write(*,'("irpamo",20I2)')(irpamo(i0),i0=1,nmo)
        write(*,'("orbmo",10F10.5)')(orbmo(i0),i0=1,nmo)

!Iwamuro modify
        Do i = 1,nmo

           If( irpmo(i) <= 8 ) then !keep irpmo
           Elseif (irpmo(i) <= 16 ) then
              goto 100 ! error
           Elseif (irpmo(i) <= 24) then
               irpmo(i) = irpmo (i) - 8
           Else
              goto 100   !error
           Endif

           If (irpmo(i) == 3) then
              irpmo(i) = 4
           Elseif (irpmo(i) == 4) then
              irpmo(i) = 3
           Elseif (irpmo(i) == 11) then
              irpmo(i) = 12
           Elseif (irpmo(i) == 12) then
              irpmo(i) = 11
           Endif

        Enddo

        write(*,*) "Modify irpmo"

        write(*,'("irpmo ",20I2)')(irpmo(i0),i0=1,nmo)

        orb = orbmo

! orb is lower order of orbmo

        do i0 = 1, nmo-1
           m = i0
           do j0 = i0+1, nmo
              if( orb(j0) < orb(m)) m = j0
           end do
           w = orb(i0) ; orb(i0) = orb(m) ; orb(m) = w
        end do

         do i0 = 1, nmo
            write(*,*)orb(i0)
         end do

         do i0 = 1, nmo
            write(*,*)orbmo(i0)
        end do

!! orb is lower order of orbmo

        do i0 = 1, nmo, 2
              m = 0
           do j0 = 1, nmo
              if (orbmo(j0)== orb(i0)) then  ! orbmo(j0) is i0 th MO
                 if( m==0) then
                    indmo(i0) = j0
                    m = m+1
                 else
                    indmo(i0+1) = j0
                 endif

              end if
           end do
        end do

        do i0 = 1,  nmo
           indmor(indmo(i0)) = i0  ! i0 is energetic order, indmo(i0) is symmtric order (MRCONEE order)
        end do

!        do i0 = 1,  nmo
!           write(*,'(2I4)')indmor(i0), indmo(i0), i0
!        end do

        orbmo = orb

        dammo = irpmo

        do i0 = 1,  nmo
           irpmo(i0) = dammo(indmo(i0))
           irpamo(i0) = dammo(indmo(i0))
        end do

           write(*,*)'inactive'
        do i0 = 1, ninact
           write(*,'(2I4,2X,E20.10,2X,I4)')i0,indmo(i0),orbmo(i0),irpmo(i0)
        end do

           write(*,*)'active'
        do i0 = ninact+1, ninact+nact
           write(*,'(2I4,2X,E20.10,2X,I4)')i0,indmo(i0),orbmo(i0),irpmo(i0)
        end do

           write(*,*)'secondary'
        do i0 = ninact+nact+1, ninact+nact+nsec
           write(*,'(2I4,2X,E20.10,2X,I4)')i0,indmo(i0),orbmo(i0),irpmo(i0)
        end do

!        do i0 = 1, nmo
!           indmo(i0)=i0
!        end do


        deallocate (dammo); Call memminus(KIND(dammo),SIZE(dammo),1)

        goto 1000

 10     write(*,*) 'err 0'
        go to 1000
 11     write(*,*) 'err 1'
        go to 1000
 12     write(*,*) 'err 2'
        go to 1000
 13     write(*,*) 'err 3'
        go to 1000
 14     write(*,*) 'err 4'
        go to 1000
 15     write(*,*) 'err 5'
        go to 1000
 100    go to 1000

 1000   end subroutine readorb_enesym