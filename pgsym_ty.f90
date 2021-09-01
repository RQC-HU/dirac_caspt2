!  =================================================

   SUBROUTINE c1sym_sd(DS) ! double-single multiplication

!  =================================================
   use four_caspt2_module

     Implicit NONE
     Integer :: DS(1,1), SD(1,1), i, j, hnsym

     NSYMRP=1
     NSYMRPA=1
     REPNA(1) ='a'; REPNA(2) ='a'

     SD( 1, 1)= 1
     DS( 1, 1)= 1

     MULTB_S = 1
     MULTB_D = 1
     MULTB_DS = 1
     irpmo=1

    end subroutine c1sym_sd

!  =================================================

   SUBROUTINE c2sym_sd(DS) ! double-single multiplication

!  =================================================
   use four_caspt2_module

     Implicit NONE
     Integer :: DS(2,2), SD(2,2), i, j, hnsym

     NSYMRP=2
     REPNA(1) ='1e1/2'; REPNA(2) ='2e1/2'
     REPNA(3) ='a    '; REPNA(4) ='b    '

!indices 1-2 when singles 1-2
!indices 1-2 when doubles 3-4

     SD( 1, 1)= 1; SD( 1, 2)= 2
     SD( 2, 1)= 2; SD( 2, 2)= 1

     DS = SD

    end subroutine c2sym_sd

!  =================================================

   SUBROUTINE c4sym_sd(DS) ! double-single multiplication

!  =================================================

   use four_caspt2_module

     Implicit NONE
     Integer :: SD(4,4), DS(4,4), i, j, hnsym

      NSYMRP=4
      REPNA(1) ='1e1/2'; REPNA(2) ='2e1/2'; REPNA(4) ='1e3/2 '; REPNA(4)       ='2e3/2';
      REPNA(5) ='a    '; REPNA(6) ='b    '; REPNA(6) ='1e    '; REPNA(8) ='2e   '

!indices 1-hnsym when singles 1-hnsym when doubles hnsym+1-nsymrp

     SD( 1, 1)= 1; SD( 1, 2)= 2;SD( 1, 3)= 3; SD( 1, 4)= 4
     SD( 2, 1)= 3; SD( 2, 2)= 4;SD( 2, 3)= 1; SD( 2, 4)= 2
     SD( 3, 1)= 4; SD( 3, 2)= 1;SD( 3, 3)= 2; SD( 3, 4)= 3
     SD( 4, 1)= 2; SD( 4, 2)= 3;SD( 4, 3)= 4; SD( 4, 4)= 1


     Do i= 1, nsymrp
     Do j= 1, nsymrp
      DS(i,j) = SD(j,i)
     End do
     End do


    end subroutine c4sym_sd

!  =================================================

   SUBROUTINE c6sym_sd(DS) ! double-single multiplication

!  =================================================

   use four_caspt2_module

     Implicit NONE
     Integer :: DS(6,6), SD(6,6), i, j, hnsym

     NSYMRP=6
     REPNA(1) ='1e1/2'; REPNA(2) ='2e1/2'; REPNA(3) ='1e3/2'; REPNA(4) ='2e3/2'; REPNA(5) ='1e5/2'; REPNA(6) ='2e5/2'
     REPNA(7) ='a   '; REPNA(8) ='b   '; REPNA(9)  ='1e1   ';         REPNA(10) = '2e1   ';REPNA(11) ='1e2  '; REPNA(12) ='2e2   '

     SD( 1, 1)= 1; SD( 1, 2)= 2;SD( 1, 3)= 3; SD( 1, 4)= 4;SD( 1, 5)= 5; SD( 1, 6)= 6
     SD( 2, 1)= 6; SD( 2, 2)= 5;SD( 2, 3)= 4; SD( 2, 4)= 3;SD( 2, 5)= 2; SD( 2, 6)= 1
     SD( 3, 1)= 2; SD( 3, 2)= 3;SD( 3, 3)= 6; SD( 3, 4)= 1;SD( 3, 5)= 4; SD( 3, 6)= 5
     SD( 4, 1)= 4; SD( 4, 2)= 1;SD( 4, 3)= 2; SD( 4, 4)= 5;SD( 4, 5)= 6; SD( 4, 6)= 3
     SD( 5, 1)= 5; SD( 5, 2)= 4;SD( 5, 3)= 1; SD( 5, 4)= 6;SD( 5, 5)= 3; SD( 5, 6)= 2
     SD( 6, 1)= 3; SD( 6, 2)= 6;SD( 6, 3)= 5; SD( 6, 4)= 2;SD( 6, 5)= 1; SD( 6, 6)= 4

     Do i= 1, nsymrp
     Do j= 1, nsymrp
      DS(i,j) = SD(j,i)
     End do
     End do


    end subroutine c6sym_sd


!  =================================================

   SUBROUTINE c8sym_sd(DS) ! double-single multiplication

!  =================================================

   use four_caspt2_module

     Implicit NONE
     Integer :: DS(8,8), SD(8,8), i, j, hnsym

     NSYMRP=8
     REPNA(1) ='1e1/2'; REPNA(2)  ='2e1/2'
     REPNA(3) ='1e3/2'; REPNA(4)  ='2e3/2' 
     REPNA(5) ='1e5/2'; REPNA(6)  ='2e5/2'
     REPNA(7) ='1e7/2'; REPNA(8)  ='2e7/2'
     
     REPNA(9)  ='a   '; REPNA(10) ='b    ' 
     REPNA(11) ='1e1 '; REPNA(12) ='2e1  '
     REPNA(13) ='1e2 '; REPNA(14) ='2e2  '
     REPNA(15) ='1e3 '; REPNA(16) ='2e3  '

     SD( 1, 1)= 1; SD( 1, 2)= 2;SD( 1, 3)= 3; SD( 1, 4)= 4
     SD( 1, 5)= 5; SD( 1, 6)= 6;SD( 1, 7)= 7; SD( 1, 8)= 8
     SD( 2, 1)= 7; SD( 2, 2)= 8;SD( 2, 3)= 5; SD( 2, 4)= 6
     SD( 2, 5)= 3; SD( 2, 6)= 4;SD( 2, 7)= 1; SD( 2, 8)= 2
     SD( 3, 1)= 2; SD( 3, 2)= 3;SD( 3, 3)= 6; SD( 3, 4)= 1
     SD( 3, 5)= 4; SD( 3, 6)= 7;SD( 3, 7)= 8; SD( 3, 8)= 5
     SD( 4, 1)= 4; SD( 4, 2)= 1;SD( 4, 3)= 2; SD( 4, 4)= 5
     SD( 4, 5)= 8; SD( 4, 6)= 3;SD( 4, 7)= 6; SD( 4, 8)= 7
     SD( 5, 1)= 3; SD( 5, 2)= 6;SD( 5, 3)= 6; SD( 5, 4)= 2
     SD( 5, 5)= 1; SD( 5, 6)= 8;SD( 5, 7)= 5; SD( 5, 8)= 4
     SD( 6, 1)= 5; SD( 6, 2)= 4;SD( 6, 3)= 1; SD( 6, 4)= 8
     SD( 6, 5)= 7; SD( 6, 6)= 2;SD( 6, 7)= 3; SD( 6, 8)= 6
     SD( 7, 1)= 8; SD( 7, 2)= 5;SD( 7, 3)= 4; SD( 7, 4)= 7
     SD( 7, 5)= 6; SD( 7, 6)= 1;SD( 7, 7)= 2; SD( 7, 8)= 3
     SD( 8, 1)= 6; SD( 8, 2)= 7;SD( 8, 3)= 8; SD( 8, 4)= 3
     SD( 8, 5)= 2; SD( 8, 6)= 5;SD( 8, 7)= 4; SD( 8, 8)= 1

     Do i= 1, nsymrp
     Do j= 1, nsymrp
      DS(i,j) = SD(j,i)
     End do
     End do


    end subroutine c8sym_sd

!  =================================================

   SUBROUTINE c2hsym_sd(DS) ! double-single multiplication

!  =================================================
   use four_caspt2_module

     Implicit NONE
     Integer :: DS(4,4), SD(4,4), i, j, hnsym

     NSYMRP=4
     REPNA(1) ='1e1/2g'; REPNA(2) ='2e1/2g'; REPNA(3) ='1e1/2u'; REPNA(4) ='2e1/2u'
     REPNA(5) ='ag    '; REPNA(6) ='bg    '; REPNA(7) ='au    '; REPNA(8) ='bu    '

!indices 1-4 when singles 1-4
!indices 1-4 when doubles 5-8

     SD( 1, 1)= 1; SD( 1, 2)= 2
     SD( 2, 1)= 2; SD( 2, 2)= 1

     hnsym= Int(nsymrp/2)
     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i,j+hnsym) = SD(i,j) + hnsym
     End do
     End do

     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i+hnsym,j) = SD(i,j+hnsym)
     End do
     End do

     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i+hnsym,j+hnsym) = SD(i,j)
     End do
     End do

     Do i= 1, nsymrp
     Do j= 1, nsymrp
      DS(i,j) = SD(j,i)
     End do
     End do

    end subroutine c2hsym_sd

!  =================================================

   SUBROUTINE c4hsym_sd (DS) ! double-single multiplication

!  =================================================
   use four_caspt2_module

     Implicit NONE
     Integer :: DS(8,8), SD(8,8), i, j, hnsym

     NSYMRP=8
      REPNA(1) ='1e1/2g'; REPNA(2) ='2e1/2g'; REPNA(3) ='1e3/2g'; REPNA(4) ='2e3/2g'
      REPNA(5) ='1e1/2u'; REPNA(6) ='2e1/2u'; REPNA(7) ='1e3/2u'; REPNA(8) ='2e3/2u'
      REPNA(9) ='ag    '; REPNA(10)='bg    '; REPNA(11)='1eg   '; REPNA(12)='2eg   '
      REPNA(13)='au    '; REPNA(14)='bu    '; REPNA(15)='1eu   '; REPNA(16)='2eu   '

!indices 1-8 when singles 1-8
!indices 1-8 when doubles 9-16

     SD( 1, 1)= 1; SD( 1, 2)= 2; SD( 1, 3)= 3; SD( 1, 4)= 4
     SD( 2, 1)= 3; SD( 2, 2)= 4; SD( 2, 3)= 1; SD( 2, 4)= 2
     SD( 3, 1)= 4; SD( 3, 2)= 1; SD( 3, 3)= 2; SD( 3, 4)= 3
     SD( 4, 1)= 2; SD( 4, 2)= 3; SD( 4, 3)= 4; SD( 4, 4)= 1

     hnsym= Int(nsymrp/2)
     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i,j+hnsym) = SD(i,j) + hnsym
     End do
     End do

     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i+hnsym,j) = SD(i,j+hnsym)
     End do
     End do

     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i+hnsym,j+hnsym) = SD(i,j)
     End do
     End do

     Do i= 1, nsymrp
     Do j= 1, nsymrp
      DS(i,j) = SD(j,i)
     End do
     End do

    end subroutine c4hsym_sd


!  =================================================

   SUBROUTINE c6hsym_sd (DS) ! double-single multiplication

!  =================================================
   use four_caspt2_module

     Implicit NONE
     Integer :: DS(12,12), SD(12,12), i, j, hnsym

     NSYMRP=12
      REPNA(1) ='1e1/2g'; REPNA(2) ='2e1/2g'; REPNA(3) ='1e3/2g'; REPNA(4) ='2e3/2g'; REPNA(5) ='1e5/2g'; REPNA(6) ='2e5/2g'
      REPNA(7) ='1e1/2u'; REPNA(8) ='2e1/2u'; REPNA(9) ='1e3/2u'; REPNA(10)='2e3/2u'; REPNA(11)='1e5/2u'; REPNA(12)='2e5/2u'
      REPNA(13)='ag    '; REPNA(14)='bg    '; REPNA(15)='1e1g  '; REPNA(16)='2e1g  '; REPNA(17)='1e2g  '; REPNA(18)='2e2g  '
      REPNA(19)='au    '; REPNA(20)='bu    '; REPNA(21)='1e1u  '; REPNA(22)='2e1u  '; REPNA(23)='1e2u  '; REPNA(24)='2e2u  '

!indices 1-hnsym when singles 1-hnsym when doubles hnsym+1-nsymrp

     SD( 1, 1)= 1; SD( 1, 2)= 2; SD( 1, 3)= 3; SD( 1, 4)= 4; SD( 1, 5)= 5; SD( 1, 6)= 6
     SD( 2, 1)= 5; SD( 2, 2)= 6; SD( 2, 3)= 4; SD( 2, 4)= 3; SD( 2, 5)= 2; SD( 2, 6)= 1
     SD( 3, 1)= 2; SD( 3, 2)= 3; SD( 3, 3)= 6; SD( 3, 4)= 1; SD( 3, 5)= 4; SD( 3, 6)= 5
     SD( 4, 1)= 4; SD( 4, 2)= 1; SD( 4, 3)= 2; SD( 4, 4)= 5; SD( 4, 5)= 6; SD( 4, 6)= 3
     SD( 5, 1)= 5; SD( 5, 2)= 4; SD( 5, 3)= 1; SD( 5, 4)= 6; SD( 5, 5)= 3; SD( 5, 6)= 2
     SD( 6, 1)= 3; SD( 6, 2)= 6; SD( 6, 3)= 5; SD( 6, 4)= 2; SD( 6, 5)= 1; SD( 6, 6)= 4


     hnsym= Int(nsymrp/2)
     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i,j+hnsym) = SD(i,j) + hnsym
     End do
     End do

     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i+hnsym,j) = SD(i,j+hnsym)
     End do
     End do

     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i+hnsym,j+hnsym) = SD(i,j)
     End do
     End do

     Do i= 1, nsymrp
     Do j= 1, nsymrp
      DS(i,j) = SD(j,i)
     End do
     End do

    end subroutine c6hsym_sd


!  =================================================

   SUBROUTINE c8hsym_sd (DS) ! double-single multiplication

!  =================================================
   use four_caspt2_module

     Implicit NONE
     Integer :: DS(16,16), SD(16,16), i, j, hnsym

     write(*,*) 'pass c8hsym'

     NSYMRP=16
      REPNA(1) ='1e1/2g'; REPNA(2) ='2e1/2g'; REPNA(3) ='1e3/2g'; REPNA(4) ='2e3/2g'
      REPNA(5) ='1e5/2g'; REPNA(6) ='2e5/2g'; REPNA(7) ='1e7/2g'; REPNA(8) ='2e7/2g'
      REPNA(9) ='1e1/2u'; REPNA(10)='2e1/2u'; REPNA(11)='1e3/2u'; REPNA(12)='2e3/2u'
      REPNA(13)='1e5/2u'; REPNA(14)='2e5/2u'; REPNA(15)='1e7/2u'; REPNA(16)='2e7/2u'

      REPNA(17)='ag    '; REPNA(18)='bg    '; REPNA(19)='1e1g  '; REPNA(20)='2e1g  '
      REPNA(21)='1e2g  '; REPNA(22)='2e2g  '; REPNA(23)='1e3g  '; REPNA(24)='2e3g  '
      REPNA(25)='au    '; REPNA(26)='bu    '; REPNA(27)='1e1u  '; REPNA(28)='2e1u  '
      REPNA(29)='1e2u  '; REPNA(30)='2e2u  '; REPNA(31)='1e3u  '; REPNA(32)='2e3u  '

!indices 1-hnsym when singles 1-hnsym when doubles hnsym+1-nsymrp

     SD( 1, 1)= 1; SD( 1, 2)= 2; SD( 1, 3)= 3; SD( 1, 4)= 4; SD( 1, 5)= 5; SD( 1, 6)= 6; SD( 1, 7)= 7; SD( 1, 8)= 8
     SD( 2, 1)= 7; SD( 2, 2)= 8; SD( 2, 3)= 5; SD( 2, 4)= 6; SD( 2, 5)= 3; SD( 2, 6)= 4; SD( 2, 7)= 1; SD( 2, 8)= 2
     SD( 3, 1)= 2; SD( 3, 2)= 3; SD( 3, 3)= 6; SD( 3, 4)= 1; SD( 3, 5)= 4; SD( 3, 6)= 7; SD( 3, 7)= 8; SD( 3, 8)= 5
     SD( 4, 1)= 4; SD( 4, 2)= 1; SD( 4, 3)= 2; SD( 4, 4)= 5; SD( 4, 5)= 8; SD( 4, 6)= 3; SD( 4, 7)= 6; SD( 4, 8)= 7
     SD( 5, 1)= 3; SD( 5, 2)= 6; SD( 5, 3)= 7; SD( 5, 4)= 2; SD( 5, 5)= 1; SD( 5, 6)= 8; SD( 5, 7)= 5; SD( 5, 8)= 4
     SD( 6, 1)= 5; SD( 6, 2)= 4; SD( 6, 3)= 1; SD( 6, 4)= 8; SD( 6, 5)= 7; SD( 6, 6)= 2; SD( 6, 7)= 3; SD( 6, 8)= 6
     SD( 7, 1)= 8; SD( 7, 2)= 5; SD( 7, 3)= 4; SD( 7, 4)= 7; SD( 7, 5)= 6; SD( 7, 6)= 5; SD( 7, 7)= 2; SD( 7, 8)= 3
     SD( 8, 1)= 6; SD( 8, 2)= 7; SD( 8, 3)= 8; SD( 8, 4)= 3; SD( 8, 5)= 2; SD( 8, 6)= 1; SD( 8, 7)= 4; SD( 8, 8)= 1

     hnsym= nsymrp/2
     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i,j+hnsym) = SD(i,j) + hnsym
     End do
     End do

     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i+hnsym,j) = SD(i,j+hnsym)
     End do
     End do

     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i+hnsym,j+hnsym) = SD(i,j)
     End do
     End do

     write(*,*)'MULTB_SD c8hsym'
     Do i = 1, nsymrp
       write(*,'(50I3)') (SD(i,j), j = 1, nsymrp)
     Enddo

     Do i= 1, nsymrp
     Do j= 1, nsymrp
      DS(i,j) = SD(j,i)
     End do
     End do

     write(*,*)'MULTB_DS c8hsym'
     Do i = 1, nsymrp
       write(*,'(50I3)') (DS(i,j), j = 1, nsymrp)
     Enddo

    end subroutine c8hsym_sd


!  =================================================

   SUBROUTINE c32hsym_sd (DS) ! double-single multiplication

!  =================================================
   use four_caspt2_module

     Implicit NONE
     Integer :: DS(32,32), SD(32,32), i, j, hnsym

     NSYMRP=32
      REPNA(1) ='1e1/2g'; REPNA(2) ='2e1/2g'; REPNA(3) ='1e3/2g'; REPNA(4) ='2e3/2g'
      REPNA(5) ='1e5/2g'; REPNA(6) ='2e5/2g'; REPNA(7) ='1e7/2g'; REPNA(8) ='2e7/2g'
      REPNA(9)='1e9/2g'; REPNA(10)='2e9/2g'; REPNA(11)='1e11/2g'; REPNA(12)='2e11/2g'
      REPNA(13)='1e13/2g'; REPNA(14)='2e13/2g'; REPNA(15)='1e15/2g'; REPNA(16)='2e15/2g'
      REPNA(17)='1e1/2u'; REPNA(18)='2e1/2u'; REPNA(19)='1e3/2u'; REPNA(20)='2e3/2u'
      REPNA(21)='1e5/2u'; REPNA(22)='2e5/2u'; REPNA(23)='1e7/2u'; REPNA(24)='2e7/2u'
      REPNA(25)='1e9/2u'; REPNA(26)='2e9/2u'; REPNA(27)='1e11/2u'; REPNA(28)='2e11/2u' 
      REPNA(29)='1e13/2u'; REPNA(30)='2e13/2u'; REPNA(31)='1e15/2u'; REPNA(32)='2e15/2u'

      REPNA(33)='ag    '; REPNA(34)='bg    '; REPNA(35)='1e1g  '; REPNA(36)='2e1g  '
      REPNA(37)='1e2g  '; REPNA(38)='2e2g  '; REPNA(39)='1e3g  '; REPNA(40)='2e3g  '
      REPNA(41)='1e4g  '; REPNA(42)='2e4g  '; REPNA(43)='1e5g  '; REPNA(44)='2e5g  '
      REPNA(45)='1e6g  '; REPNA(46)='2e6g  '; REPNA(47)='1e7g  '; REPNA(48)='2e7g  '
      REPNA(49)='au    '; REPNA(50)='bu    '; REPNA(51)='1e1u  '; REPNA(52)='2e1u  '
      REPNA(53)='1e2u  '; REPNA(54)='2e2u  '; REPNA(55)='1e3u  '; REPNA(56)='2e3u  '
      REPNA(57)='1e4u  '; REPNA(58)='2e4u  '; REPNA(59)='1e5u  '; REPNA(60)='2e5u  '
      REPNA(61)='1e7u  '; REPNA(62)='2e7u  '; REPNA(63)='1e9u  '; REPNA(64)='2e9u  '


!indices 1-hnsym when singles 1-hnsym when doubles hnsym+1-nsymrp

     SD( 1, 1)= 1; SD( 1, 2)= 2; SD( 1, 3)= 3; SD( 1, 4)= 4; SD( 1, 5)= 5; SD( 1, 6)= 6; SD( 1, 7)= 7; SD( 1, 8)= 8; SD( 1, 9)= 9; SD( 1, 10)= 10; SD( 1, 11)= 11; SD( 1, 12)= 12; SD( 1, 13)= 13; SD( 1, 14)= 14; SD( 1, 15)= 15; SD( 1, 16)= 16
     SD( 2, 1)= 3; SD( 2, 2)= 1; SD( 2, 3)= 5; SD( 2, 4)= 2; SD( 2, 5)= 7; SD( 2, 6)= 4; SD( 2, 7)= 9; SD( 2, 8)= 6; SD( 2, 9)= 11; SD( 2, 10)= 8; SD( 2, 11)= 13; SD( 2, 12)= 10; SD( 2, 13)= 15; SD( 2, 14)= 12; SD( 2, 15)= 16; SD( 2, 16)= 14
     SD( 3, 1)= 2; SD( 3, 2)= 4; SD( 3, 3)= 1; SD( 3, 4)= 6; SD( 3, 5)= 3; SD( 3, 6)= 8; SD( 3, 7)= 5; SD( 3, 8)= 10; SD( 3, 9)= 7; SD( 3, 10)= 12; SD( 3, 11)= 9; SD( 3, 12)= 14; SD( 3, 13)= 11; SD( 3, 14)= 16; SD( 3, 15)= 13; SD( 3, 16)= 15
     SD( 4, 1)= 5; SD( 4, 2)= 3; SD( 4, 3)= 7; SD( 4, 4)= 1; SD( 4, 5)= 9; SD( 4, 6)= 2; SD( 4, 7)= 11; SD( 4, 8)= 4; SD( 4, 9)= 13; SD( 4, 10)= 6; SD( 4, 11)= 15; SD( 4, 12)= 8; SD( 4, 13)= 16; SD( 4, 14)= 10; SD( 4, 15)= 14; SD( 4, 16)= 12
     SD( 5, 1)= 4; SD( 5, 2)= 6; SD( 5, 3)= 2; SD( 5, 4)= 8; SD( 5, 5)= 1; SD( 5, 6)= 10; SD( 5, 7)= 3; SD( 5, 8)= 12; SD( 5, 9)= 5; SD( 5, 10)= 14; SD( 5, 11)= 7; SD( 5, 12)= 16; SD( 5, 13)= 9; SD( 5, 14)= 15; SD( 5, 15)= 11; SD( 5, 16)= 13
     SD( 6, 1)= 7; SD( 6, 2)= 5; SD( 6, 3)= 9; SD( 6, 4)= 3; SD( 6, 5)= 11; SD( 6, 6)= 1; SD( 6, 7)= 13; SD( 6, 8)= 2; SD( 6, 9)= 15; SD( 6, 10)= 4; SD( 6, 11)= 16; SD( 6, 12)= 6; SD( 6, 13)= 14; SD( 6, 14)= 8; SD( 6, 15)= 12; SD( 6, 16)= 10
     SD( 7, 1)= 6; SD( 7, 2)= 8; SD( 7, 3)= 4; SD( 7, 4)= 10; SD( 7, 5)= 2; SD( 7, 6)= 12; SD( 7, 7)= 1; SD( 7, 8)= 14; SD( 7, 9)= 3; SD( 7, 10)= 16; SD( 7, 11)= 5; SD( 7, 12)= 15; SD( 7, 13)= 7; SD( 7, 14)= 13; SD( 7, 15)= 9; SD( 7, 16)= 11
     SD( 8, 1)= 9; SD( 8, 2)= 7; SD( 8, 3)= 11; SD( 8, 4)= 5; SD( 8, 5)= 13; SD( 8, 6)= 3; SD( 8, 7)= 15; SD( 8, 8)= 1; SD( 8, 9)= 16; SD( 8, 10)= 2; SD( 8, 11)= 14; SD( 8, 12)= 4; SD( 8, 13)= 12; SD( 8, 14)= 6; SD( 8, 15)= 10; SD( 8, 16)= 8
     SD( 9, 1)= 8; SD( 9, 2)= 10; SD( 9, 3)= 6; SD( 9, 4)= 12; SD( 9, 5)= 4; SD( 9, 6)= 14; SD( 9, 7)= 2; SD( 9, 8)= 16; SD( 9, 9)= 1; SD( 9, 10)= 15; SD( 9, 11)= 3; SD( 9, 12)= 13; SD( 9, 13)= 5; SD( 9, 14)= 11; SD( 9, 15)= 7; SD( 9, 16)= 9
     SD( 10, 1)= 11; SD( 10, 2)= 9; SD( 10, 3)= 13; SD( 10, 4)= 7; SD( 10, 5)= 15; SD( 10, 6)= 5; SD( 10, 7)= 16; SD( 10, 8)= 3; SD( 10, 9)= 14; SD( 10, 10)= 1; SD( 10, 11)= 12; SD( 10, 12)= 2; SD( 10, 13)= 10; SD( 10, 14)= 4; SD( 10, 15)= 8; SD( 10, 16)= 6
     SD( 11, 1)= 10; SD( 11, 2)= 12; SD( 11, 3)= 8; SD( 11, 4)= 14; SD( 11, 5)= 6; SD( 11, 6)= 16; SD( 11, 7)= 4; SD( 11, 8)= 15; SD( 11, 9)= 2; SD( 11, 10)= 13; SD( 11, 11)= 1; SD( 11, 12)= 11; SD( 11, 13)= 3; SD( 11, 14)= 9; SD( 11, 15)= 5; SD( 11, 16)= 7
     SD( 12, 1)= 13; SD( 12, 2)= 11; SD( 12, 3)= 15; SD( 12, 4)= 9; SD( 12, 5)= 16; SD( 12, 6)= 7; SD( 12, 7)= 14; SD( 12, 8)= 5; SD( 12, 9)= 12; SD( 12, 10)= 3; SD( 12, 11)= 10; SD( 12, 12)= 1; SD( 12, 13)= 8; SD( 12, 14)= 2; SD( 12, 15)= 6; SD( 12, 16)= 4
     SD( 13, 1)= 12; SD( 13, 2)= 14; SD( 13, 3)= 10; SD( 13, 4)= 16; SD( 13, 5)= 8; SD( 13, 6)= 15; SD( 13, 7)= 6; SD( 13, 8)= 13; SD( 13, 9)= 4; SD( 13, 10)= 11; SD( 13, 11)= 2; SD( 13, 12)= 9; SD( 13, 13)= 1; SD( 13, 14)= 7; SD( 13, 15)= 3; SD( 13, 16)= 5
     SD( 14, 1)= 15; SD( 14, 2)= 13; SD( 14, 3)= 16; SD( 14, 4)= 11; SD( 14, 5)= 14; SD( 14, 6)= 9; SD( 14, 7)= 12; SD( 14, 8)= 7; SD( 14, 9)= 10; SD( 14, 10)= 5; SD( 14, 11)= 8; SD( 14, 12)= 3; SD( 14, 13)= 6; SD( 14, 14)= 1; SD( 14, 15)= 4; SD( 14, 16)= 2
     SD( 15, 1)= 14; SD( 15, 2)= 16; SD( 15, 3)= 12; SD( 15, 4)= 15; SD( 15, 5)= 10; SD( 15, 6)= 13; SD( 15, 7)= 8; SD( 15, 8)= 11; SD( 15, 9)= 6; SD( 15, 10)= 9; SD( 15, 11)= 4; SD( 15, 12)= 7; SD( 15, 13)= 2; SD( 15, 14)= 5; SD( 15, 15)= 1; SD( 15, 16)= 3
     SD( 16, 1)= 16; SD( 16, 2)= 15; SD( 16, 3)= 14; SD( 16, 4)= 13; SD( 16, 5)= 12; SD( 16, 6)= 11; SD( 16, 7)= 10; SD( 16, 8)= 9; SD( 16, 9)= 8; SD( 16, 10)= 7; SD( 16, 11)= 6; SD( 16, 12)= 5; SD( 16, 13)= 4; SD( 16, 14)= 3; SD( 16, 15)= 2; SD( 16, 16)= 1

     hnsym= nsymrp/2
     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i,j+hnsym) = SD(i,j) + hnsym
     End do
     End do

     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i+hnsym,j) = SD(i,j+hnsym)
     End do
     End do

     Do i= 1, hnsym
     Do j= 1, hnsym
      SD(i+hnsym,j+hnsym) = SD(i,j)
     End do
     End do

     write(*,*)'MULTB_SD c32hsym'
     Do i = 1, nsymrp
       write(*,'(50I3)') (SD(i,j), j = 1, nsymrp)
     Enddo

     Do i= 1, nsymrp
     Do j= 1, nsymrp
      DS(i,j) = SD(j,i)
     End do
     End do

     write(*,*)'MULTB_DS c32hsym'
     Do i = 1, nsymrp
       write(*,'(50I3)') (DS(i,j), j = 1, nsymrp)
     Enddo

    end subroutine c32hsym_sd





