c       This is a subset of routines derived from FFTPACK5
c               for forward and backward real Fourier Transforms
c
c       file contains:
c
c     RFFT1I 
c     RFFTI1
c     RFFT1F 
c     RFFTF1 
c     XERFFT
c     R1F4KF 
c     R1F2KF 
c     R1F5KF 
c     R1FGKF 
c     R1F3KF 
c     RFFT1B 
c     RFFTB1 
c     R1F2KB 
c     R1F3KB 
c     R1F4KB 
c     R1F5KB 
c     R1FGKB 
c
c     These routines were obtained from:
c     http://www.cisl.ucar.edu/css/software/fftpack5/


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: rfft1i.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RFFT1I ( N, WSAVE, LENSAV, IER )
      INTEGER    N, LENSAV, IER
      REAL       WSAVE(LENSAV)
C
      IER = 0
C
      IF (LENSAV .LT. N + INT(LOG(REAL(N))) +4) THEN
        IER = 2
        CALL XERFFT ('RFFT1I ', 3)
      ENDIF
C
      IF (N .EQ. 1) RETURN
C
      CALL RFFTI1 (N,WSAVE(1),WSAVE(N+1))
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: rffti1.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RFFTI1 (N,WA,FAC)
      REAL       WA(N)      ,FAC(15)
      INTEGER    NTRYH(4)
      DOUBLE PRECISION TPI,ARGH,ARGLD,ARG
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
C
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      FAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         FAC(IB+2) = FAC(IB+1)
  106 CONTINUE
      FAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      FAC(1) = N
      FAC(2) = NF
      TPI = 8.D0*DATAN(1.D0)
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = FAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = FLOAT(LD)*ARGH
            FI = 0.
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = DCOS(ARG)
               WA(I) = DSIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: rfft1f.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RFFT1F ( N, INC, R, LENR, WSAVE, LENSAV,
     1                  WORK, LENWRK, IER)
      INTEGER  N, INC, LENR, LENSAV, LENWRK, IER
      REAL     R(LENR), WSAVE(LENSAV), WORK(LENWRK)
C
      IER = 0
C
      IF (LENR .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('RFFT1F ', 6)
      ELSEIF (LENSAV .LT. N + INT(LOG(REAL(N))) +4) THEN
        IER = 2
        CALL XERFFT ('RFFT1F ', 8)
      ELSEIF (LENWRK .LT. N) THEN
        IER = 3
        CALL XERFFT ('RFFT1F ', 10)
      ENDIF
C
      IF (N .EQ. 1) RETURN
C
      CALL RFFTF1 (N,INC,R,WORK,WSAVE,WSAVE(N+1))
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: rfftf1.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RFFTF1 (N,IN,C,CH,WA,FAC)
      REAL       CH(*) ,C(IN,*)  ,WA(N)   ,FAC(15)
C
      NF = FAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = FAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL R1F4KF (IDO,L1,C,IN,CH,1,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL R1F4KF (IDO,L1,CH,1,C,IN,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL R1F2KF (IDO,L1,C,IN,CH,1,WA(IW))
         GO TO 110
  103    CALL R1F2KF (IDO,L1,CH,1,C,IN,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL R1F3KF (IDO,L1,C,IN,CH,1,WA(IW),WA(IX2))
         GO TO 110
  105    CALL R1F3KF (IDO,L1,CH,1,C,IN,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
         CALL R1F5KF (IDO,L1,C,IN,CH,1,WA(IW),WA(IX2),
     1                      WA(IX3),WA(IX4))
         GO TO 110
  107    CALL R1F5KF (IDO,L1,CH,1,C,IN,WA(IW),WA(IX2),
     1                      WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL R1FGKF (IDO,IP,L1,IDL1,C,C,C,IN,CH,CH,1,WA(IW))
         NA = 1
         GO TO 110
  109    CALL R1FGKF (IDO,IP,L1,IDL1,CH,CH,CH,1,C,C,IN,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      SN = 1./N
      TSN = 2./N
      TSNM = -TSN
      MODN = MOD(N,2)
      NL = N-2
      IF(MODN .NE. 0) NL = N-1
      IF (NA .NE. 0) GO TO 120
      C(1,1) = SN*CH(1)
      DO 118 J=2,NL,2
         C(1,J) = TSN*CH(J)
         C(1,J+1) = TSNM*CH(J+1)
  118 CONTINUE
      IF(MODN .NE. 0) RETURN
      C(1,N) = SN*CH(N)
      RETURN
  120 C(1,1) = SN*C(1,1)
      DO 122 J=2,NL,2
         C(1,J) = TSN*C(1,J)
         C(1,J+1) = TSNM*C(1,J+1)
  122 CONTINUE
      IF(MODN .NE. 0) RETURN
      C(1,N) = SN*C(1,N)
      RETURN
      END



      SUBROUTINE XERFFT( SRNAME, INFO)
C
C     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
C
C     ..
C
C  Purpose
C  =======
C
C  XERFFT  is an error handler for library FFTPACK version 5.0 routines.
C  It is called by an FFTPACK 5.0 routine if an input parameter has an
C  invalid value.  A message is printed and execution stops.
C
C  Installers may consider modifying the STOP statement in order to
C  call system-specific exception-handling facilities.
C
C  Arguments
C  =========
C
C  SRNAME  (input) CHARACTER*6
C          The name of the routine which called XERFFT.
C
C  INFO    (input) INTEGER
C          When a single  invalid parameter in the parameter list of
C          the calling routine has been detected, INFO is the position
C          of that parameter.  In the case when an illegal combination
C          of LOT, JUMP, N, and INC has been detected, the calling
C          subprogram calls XERFFT with INFO = -1.
C
C =====================================================================
C
C     .. Executable Statements ..
C
      IF (INFO .GE. 1) THEN
        WRITE( *, '(A,A,A,I3,A)') ' ** On entry to ', SRNAME,
     1    ' parameter number ', INFO, ' had an illegal value'
      ELSEIF (INFO .EQ. -1) THEN
        WRITE( *, '(A,A,A,A)') ' ** On entry to ', SRNAME,
     1    ' parameters LOT, JUMP, N and INC are inconsistent'
      ELSEIF (INFO .EQ. -2) THEN
        WRITE( *, '(A,A,A,A)') ' ** On entry to ', SRNAME,
     1    ' parameter L is greater than LDIM'
      ELSEIF (INFO .EQ. -3) THEN
        WRITE( *, '(A,A,A,A)') ' ** On entry to ', SRNAME,
     1    ' parameter M is greater than MDIM'
      ELSEIF (INFO .EQ. -5) THEN
        WRITE( *, '(A,A,A,A)') ' ** Within ', SRNAME,
     1    ' input error returned by lower level routine'
      ELSEIF (INFO .EQ. -6) THEN
        WRITE( *, '(A,A,A,A)') ' ** On entry to ', SRNAME,
     1    ' parameter LDIM is less than 2*(L/2+1)'
      ENDIF
*
      STOP
*
*     End of XERFFT
*
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: r1f4kf.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE R1F4KF (IDO,L1,CC,IN1,CH,IN2,WA1,WA2,WA3)
      REAL       CC(IN1,IDO,L1,4)   ,CH(IN2,IDO,4,L1)     ,
     1           WA1(IDO)           ,WA2(IDO)     ,WA3(IDO)
C
      HSQT2=SQRT(2.)/2.
      DO 101 K=1,L1
         CH(1,1,1,K) = (CC(1,1,K,2)+CC(1,1,K,4))
     1      +(CC(1,1,K,1)+CC(1,1,K,3))
         CH(1,IDO,4,K) = (CC(1,1,K,1)+CC(1,1,K,3))
     1      -(CC(1,1,K,2)+CC(1,1,K,4))
         CH(1,IDO,2,K) = CC(1,1,K,1)-CC(1,1,K,3)
         CH(1,1,3,K) = CC(1,1,K,4)-CC(1,1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CH(1,I-1,1,K) = ((WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*
     1       CC(1,I,K,2))+(WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*
     1       CC(1,I,K,4)))+(CC(1,I-1,K,1)+(WA2(I-2)*CC(1,I-1,K,3)+
     1       WA2(I-1)*CC(1,I,K,3)))
            CH(1,IC-1,4,K) = (CC(1,I-1,K,1)+(WA2(I-2)*CC(1,I-1,K,3)+
     1       WA2(I-1)*CC(1,I,K,3)))-((WA1(I-2)*CC(1,I-1,K,2)+
     1       WA1(I-1)*CC(1,I,K,2))+(WA3(I-2)*CC(1,I-1,K,4)+
     1       WA3(I-1)*CC(1,I,K,4)))
            CH(1,I,1,K) = ((WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*
     1       CC(1,I-1,K,2))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*
     1       CC(1,I-1,K,4)))+(CC(1,I,K,1)+(WA2(I-2)*CC(1,I,K,3)-
     1       WA2(I-1)*CC(1,I-1,K,3)))
            CH(1,IC,4,K) = ((WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*
     1       CC(1,I-1,K,2))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*
     1       CC(1,I-1,K,4)))-(CC(1,I,K,1)+(WA2(I-2)*CC(1,I,K,3)-
     1       WA2(I-1)*CC(1,I-1,K,3)))
            CH(1,I-1,3,K) = ((WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*
     1       CC(1,I-1,K,2))-(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*
     1       CC(1,I-1,K,4)))+(CC(1,I-1,K,1)-(WA2(I-2)*CC(1,I-1,K,3)+
     1       WA2(I-1)*CC(1,I,K,3)))
            CH(1,IC-1,2,K) = (CC(1,I-1,K,1)-(WA2(I-2)*CC(1,I-1,K,3)+
     1       WA2(I-1)*CC(1,I,K,3)))-((WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*
     1       CC(1,I-1,K,2))-(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*
     1       CC(1,I-1,K,4)))
            CH(1,I,3,K) = ((WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*
     1       CC(1,I,K,4))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*
     1       CC(1,I,K,2)))+(CC(1,I,K,1)-(WA2(I-2)*CC(1,I,K,3)-
     1       WA2(I-1)*CC(1,I-1,K,3)))
            CH(1,IC,2,K) = ((WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*
     1       CC(1,I,K,4))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*
     1       CC(1,I,K,2)))-(CC(1,I,K,1)-(WA2(I-2)*CC(1,I,K,3)-
     1       WA2(I-1)*CC(1,I-1,K,3)))
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
            CH(1,IDO,1,K) = (HSQT2*(CC(1,IDO,K,2)-CC(1,IDO,K,4)))+
     1       CC(1,IDO,K,1)
            CH(1,IDO,3,K) = CC(1,IDO,K,1)-(HSQT2*(CC(1,IDO,K,2)-
     1       CC(1,IDO,K,4)))
            CH(1,1,2,K) = (-HSQT2*(CC(1,IDO,K,2)+CC(1,IDO,K,4)))-
     1       CC(1,IDO,K,3)
            CH(1,1,4,K) = (-HSQT2*(CC(1,IDO,K,2)+CC(1,IDO,K,4)))+
     1       CC(1,IDO,K,3)
  106 CONTINUE
  107 RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: r1f2kf.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE R1F2KF (IDO,L1,CC,IN1,CH,IN2,WA1)
      REAL       CH(IN2,IDO,2,L1) ,CC(IN1,IDO,L1,2) , WA1(IDO)
C
      DO 101 K=1,L1
         CH(1,1,1,K) = CC(1,1,K,1)+CC(1,1,K,2)
         CH(1,IDO,2,K) = CC(1,1,K,1)-CC(1,1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CH(1,I,1,K) = CC(1,I,K,1)+(WA1(I-2)*CC(1,I,K,2)-
     1       WA1(I-1)*CC(1,I-1,K,2))
            CH(1,IC,2,K) = (WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*
     1       CC(1,I-1,K,2))-CC(1,I,K,1)
            CH(1,I-1,1,K) = CC(1,I-1,K,1)+(WA1(I-2)*CC(1,I-1,K,2)+
     1       WA1(I-1)*CC(1,I,K,2))
            CH(1,IC-1,2,K) = CC(1,I-1,K,1)-(WA1(I-2)*CC(1,I-1,K,2)+
     1       WA1(I-1)*CC(1,I,K,2))
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(1,1,2,K) = -CC(1,IDO,K,2)
         CH(1,IDO,1,K) = CC(1,IDO,K,1)
  106 CONTINUE
  107 RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: r1f5kf.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE R1F5KF (IDO,L1,CC,IN1,CH,IN2,
     1                   WA1,WA2,WA3,WA4)
      REAL       CC(IN1,IDO,L1,5)    ,CH(IN2,IDO,5,L1)     ,
     1           WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)     ,WA4(IDO)
C
      ARG=2.*4.*ATAN(1.0)/5.
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2.*ARG)
      TI12=SIN(2.*ARG)
      DO 101 K=1,L1
         CH(1,1,1,K) = CC(1,1,K,1)+(CC(1,1,K,5)+CC(1,1,K,2))+
     1    (CC(1,1,K,4)+CC(1,1,K,3))
         CH(1,IDO,2,K) = CC(1,1,K,1)+TR11*(CC(1,1,K,5)+CC(1,1,K,2))+
     1    TR12*(CC(1,1,K,4)+CC(1,1,K,3))
         CH(1,1,3,K) = TI11*(CC(1,1,K,5)-CC(1,1,K,2))+TI12*
     1    (CC(1,1,K,4)-CC(1,1,K,3))
         CH(1,IDO,4,K) = CC(1,1,K,1)+TR12*(CC(1,1,K,5)+CC(1,1,K,2))+
     1    TR11*(CC(1,1,K,4)+CC(1,1,K,3))
         CH(1,1,5,K) = TI12*(CC(1,1,K,5)-CC(1,1,K,2))-TI11*
     1    (CC(1,1,K,4)-CC(1,1,K,3))
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            CH(1,I-1,1,K) = CC(1,I-1,K,1)+((WA1(I-2)*CC(1,I-1,K,2)+
     1       WA1(I-1)*CC(1,I,K,2))+(WA4(I-2)*CC(1,I-1,K,5)+WA4(I-1)*
     1       CC(1,I,K,5)))+((WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*
     1       CC(1,I,K,3))+(WA3(I-2)*CC(1,I-1,K,4)+
     1       WA3(I-1)*CC(1,I,K,4)))
            CH(1,I,1,K) = CC(1,I,K,1)+((WA1(I-2)*CC(1,I,K,2)-
     1       WA1(I-1)*CC(1,I-1,K,2))+(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*
     1       CC(1,I-1,K,5)))+((WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*
     1       CC(1,I-1,K,3))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*
     1       CC(1,I-1,K,4)))
            CH(1,I-1,3,K) = CC(1,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2)
     1       +WA4(I-2)*CC(1,I-1,K,5)+WA4(I-1)*CC(1,I,K,5))+TR12*
     1      ( WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*CC(1,I,K,3)
     1       +WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*CC(1,I,K,4))+TI11*
     1      ( WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2)
     1       -(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*CC(1,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*CC(1,I-1,K,3)
     1       -(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*CC(1,I-1,K,4)))
            CH(1,IC-1,2,K) = CC(1,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2)
     1       +WA4(I-2)*CC(1,I-1,K,5)+WA4(I-1)*CC(1,I,K,5))+TR12*
     1     ( WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*CC(1,I,K,3)
     1      +WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*CC(1,I,K,4))-(TI11*
     1      ( WA1(I-2)*CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2)
     1       -(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*CC(1,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*CC(1,I-1,K,3)
     1       -(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*CC(1,I-1,K,4))))
            CH(1,I,3,K) = (CC(1,I,K,1)+TR11*((WA1(I-2)*CC(1,I,K,2)-
     1       WA1(I-1)*CC(1,I-1,K,2))+(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*
     1       CC(1,I-1,K,5)))+TR12*((WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*
     1       CC(1,I-1,K,3))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*
     1       CC(1,I-1,K,4))))+(TI11*((WA4(I-2)*CC(1,I-1,K,5)+
     1       WA4(I-1)*CC(1,I,K,5))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*
     1       CC(1,I,K,2)))+TI12*((WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*
     1       CC(1,I,K,4))-(WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*
     1       CC(1,I,K,3))))
            CH(1,IC,2,K) = (TI11*((WA4(I-2)*CC(1,I-1,K,5)+WA4(I-1)*
     1       CC(1,I,K,5))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*
     1       CC(1,I,K,2)))+TI12*((WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*
     1       CC(1,I,K,4))-(WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*
     1       CC(1,I,K,3))))-(CC(1,I,K,1)+TR11*((WA1(I-2)*CC(1,I,K,2)-
     1       WA1(I-1)*CC(1,I-1,K,2))+(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*
     1       CC(1,I-1,K,5)))+TR12*((WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*
     1       CC(1,I-1,K,3))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*
     1       CC(1,I-1,K,4))))
            CH(1,I-1,5,K) = (CC(1,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2))+(WA4(I-2)*
     1       CC(1,I-1,K,5)+WA4(I-1)*CC(1,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(1,I-1,K,3)+WA2(I-1)*CC(1,I,K,3))+(WA3(I-2)*
     1       CC(1,I-1,K,4)+WA3(I-1)*CC(1,I,K,4))))+(TI12*((WA1(I-2)*
     1       CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2))-(WA4(I-2)*
     1       CC(1,I,K,5)-WA4(I-1)*CC(1,I-1,K,5)))-TI11*((WA2(I-2)*
     1       CC(1,I,K,3)-WA2(I-1)*CC(1,I-1,K,3))-(WA3(I-2)*
     1       CC(1,I,K,4)-WA3(I-1)*CC(1,I-1,K,4))))
            CH(1,IC-1,4,K) = (CC(1,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2))+(WA4(I-2)*
     1       CC(1,I-1,K,5)+WA4(I-1)*CC(1,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(1,I-1,K,3)+WA2(I-1)*CC(1,I,K,3))+(WA3(I-2)*
     1       CC(1,I-1,K,4)+WA3(I-1)*CC(1,I,K,4))))-(TI12*((WA1(I-2)*
     1       CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2))-(WA4(I-2)*
     1       CC(1,I,K,5)-WA4(I-1)*CC(1,I-1,K,5)))-TI11*((WA2(I-2)*
     1       CC(1,I,K,3)-WA2(I-1)*CC(1,I-1,K,3))-(WA3(I-2)*
     1       CC(1,I,K,4)-WA3(I-1)*CC(1,I-1,K,4))))
            CH(1,I,5,K) = (CC(1,I,K,1)+TR12*((WA1(I-2)*CC(1,I,K,2)-
     1       WA1(I-1)*CC(1,I-1,K,2))+(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*
     1       CC(1,I-1,K,5)))+TR11*((WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*
     1       CC(1,I-1,K,3))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*
     1       CC(1,I-1,K,4))))+(TI12*((WA4(I-2)*CC(1,I-1,K,5)+
     1       WA4(I-1)*CC(1,I,K,5))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*
     1       CC(1,I,K,2)))-TI11*((WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*
     1       CC(1,I,K,4))-(WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*
     1       CC(1,I,K,3))))
            CH(1,IC,4,K) = (TI12*((WA4(I-2)*CC(1,I-1,K,5)+WA4(I-1)*
     1       CC(1,I,K,5))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*
     1       CC(1,I,K,2)))-TI11*((WA3(I-2)*CC(1,I-1,K,4)+WA3(I-1)*
     1       CC(1,I,K,4))-(WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*
     1       CC(1,I,K,3))))-(CC(1,I,K,1)+TR12*((WA1(I-2)*CC(1,I,K,2)-
     1       WA1(I-1)*CC(1,I-1,K,2))+(WA4(I-2)*CC(1,I,K,5)-WA4(I-1)*
     1       CC(1,I-1,K,5)))+TR11*((WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*
     1       CC(1,I-1,K,3))+(WA3(I-2)*CC(1,I,K,4)-WA3(I-1)*
     1       CC(1,I-1,K,4))))
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: r1fgkf.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE R1FGKF (IDO,IP,L1,IDL1,CC,C1,C2,IN1,
     1              CH,CH2,IN2,WA)
      REAL          CH(IN2,IDO,L1,IP)   ,CC(IN1,IDO,IP,L1),
     1              C1(IN1,IDO,L1,IP)   ,C2(IN1,IDL1,IP),
     2              CH2(IN2,IDL1,IP)    ,WA(IDO)
C
      TPI=2.*4.*ATAN(1.0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         CH2(1,IK,1) = C2(1,IK,1)
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            CH(1,1,K,J) = C1(1,1,K,J)
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               CH(1,I-1,K,J) = WA(IDIJ-1)*C1(1,I-1,K,J)+WA(IDIJ)
     1           *C1(1,I,K,J)
               CH(1,I,K,J) = WA(IDIJ-1)*C1(1,I,K,J)-WA(IDIJ)
     1           *C1(1,I-1,K,J)
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               CH(1,I-1,K,J) = WA(IDIJ-1)*C1(1,I-1,K,J)+WA(IDIJ)
     1           *C1(1,I,K,J)
               CH(1,I,K,J) = WA(IDIJ-1)*C1(1,I,K,J)-WA(IDIJ)
     1           *C1(1,I-1,K,J)
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
            DO 112 I=3,IDO,2
               C1(1,I-1,K,J) = CH(1,I-1,K,J)+CH(1,I-1,K,JC)
               C1(1,I-1,K,JC) = CH(1,I,K,J)-CH(1,I,K,JC)
               C1(1,I,K,J) = CH(1,I,K,J)+CH(1,I,K,JC)
               C1(1,I,K,JC) = CH(1,I-1,K,JC)-CH(1,I-1,K,J)
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               C1(1,I-1,K,J) = CH(1,I-1,K,J)+CH(1,I-1,K,JC)
               C1(1,I-1,K,JC) = CH(1,I,K,J)-CH(1,I,K,JC)
               C1(1,I,K,J) = CH(1,I,K,J)+CH(1,I,K,JC)
               C1(1,I,K,JC) = CH(1,I-1,K,JC)-CH(1,I-1,K,J)
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         C2(1,IK,1) = CH2(1,IK,1)
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            C1(1,1,K,J) = CH(1,1,K,J)+CH(1,1,K,JC)
            C1(1,1,K,JC) = CH(1,1,K,JC)-CH(1,1,K,J)
  122    CONTINUE
  123 CONTINUE
C
      AR1 = 1.
      AI1 = 0.
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            CH2(1,IK,L) = C2(1,IK,1)+AR1*C2(1,IK,2)
            CH2(1,IK,LC) = AI1*C2(1,IK,IP)
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               CH2(1,IK,L) = CH2(1,IK,L)+AR2*C2(1,IK,J)
               CH2(1,IK,LC) = CH2(1,IK,LC)+AI2*C2(1,IK,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            CH2(1,IK,1) = CH2(1,IK,1)+C2(1,IK,J)
  128    CONTINUE
  129 CONTINUE
C
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            CC(1,I,1,K) = CH(1,I,K,1)
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            CC(1,I,1,K) = CH(1,I,K,1)
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            CC(1,IDO,J2-2,K) = CH(1,1,K,J)
            CC(1,1,J2-1,K) = CH(1,1,K,JC)
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
            DO 138 I=3,IDO,2
               IC = IDP2-I
               CC(1,I-1,J2-1,K) = CH(1,I-1,K,J)+CH(1,I-1,K,JC)
               CC(1,IC-1,J2-2,K) = CH(1,I-1,K,J)-CH(1,I-1,K,JC)
               CC(1,I,J2-1,K) = CH(1,I,K,J)+CH(1,I,K,JC)
               CC(1,IC,J2-2,K) = CH(1,I,K,JC)-CH(1,I,K,J)
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               CC(1,I-1,J2-1,K) = CH(1,I-1,K,J)+CH(1,I-1,K,JC)
               CC(1,IC-1,J2-2,K) = CH(1,I-1,K,J)-CH(1,I-1,K,JC)
               CC(1,I,J2-1,K) = CH(1,I,K,J)+CH(1,I,K,JC)
               CC(1,IC,J2-2,K) = CH(1,I,K,JC)-CH(1,I,K,J)
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: r1f3kf.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE R1F3KF (IDO,L1,CC,IN1,CH,IN2,WA1,WA2)
      REAL       CH(IN2,IDO,3,L1)  ,CC(IN1,IDO,L1,3)     ,
     1                WA1(IDO)     ,WA2(IDO)
C
      ARG=2.*4.*ATAN(1.0)/3.
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
         CH(1,1,1,K) = CC(1,1,K,1)+(CC(1,1,K,2)+CC(1,1,K,3))
         CH(1,1,3,K) = TAUI*(CC(1,1,K,3)-CC(1,1,K,2))
         CH(1,IDO,2,K) = CC(1,1,K,1)+TAUR*
     1      (CC(1,1,K,2)+CC(1,1,K,3))
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            CH(1,I-1,1,K) = CC(1,I-1,K,1)+((WA1(I-2)*CC(1,I-1,K,2)+
     1       WA1(I-1)*CC(1,I,K,2))+(WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*
     1       CC(1,I,K,3)))
            CH(1,I,1,K) = CC(1,I,K,1)+((WA1(I-2)*CC(1,I,K,2)-
     1       WA1(I-1)*CC(1,I-1,K,2))+(WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*
     1       CC(1,I-1,K,3)))
            CH(1,I-1,3,K) = (CC(1,I-1,K,1)+TAUR*((WA1(I-2)*
     1       CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2))+(WA2(I-2)*
     1       CC(1,I-1,K,3)+WA2(I-1)*CC(1,I,K,3))))+(TAUI*((WA1(I-2)*
     1       CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2))-(WA2(I-2)*
     1       CC(1,I,K,3)-WA2(I-1)*CC(1,I-1,K,3))))
            CH(1,IC-1,2,K) = (CC(1,I-1,K,1)+TAUR*((WA1(I-2)*
     1       CC(1,I-1,K,2)+WA1(I-1)*CC(1,I,K,2))+(WA2(I-2)*
     1       CC(1,I-1,K,3)+WA2(I-1)*CC(1,I,K,3))))-(TAUI*((WA1(I-2)*
     1       CC(1,I,K,2)-WA1(I-1)*CC(1,I-1,K,2))-(WA2(I-2)*
     1       CC(1,I,K,3)-WA2(I-1)*CC(1,I-1,K,3))))
            CH(1,I,3,K) = (CC(1,I,K,1)+TAUR*((WA1(I-2)*CC(1,I,K,2)-
     1       WA1(I-1)*CC(1,I-1,K,2))+(WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*
     1       CC(1,I-1,K,3))))+(TAUI*((WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*
     1       CC(1,I,K,3))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*
     1       CC(1,I,K,2))))
            CH(1,IC,2,K) = (TAUI*((WA2(I-2)*CC(1,I-1,K,3)+WA2(I-1)*
     1       CC(1,I,K,3))-(WA1(I-2)*CC(1,I-1,K,2)+WA1(I-1)*
     1       CC(1,I,K,2))))-(CC(1,I,K,1)+TAUR*((WA1(I-2)*CC(1,I,K,2)-
     1       WA1(I-1)*CC(1,I-1,K,2))+(WA2(I-2)*CC(1,I,K,3)-WA2(I-1)*
     1       CC(1,I-1,K,3))))
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: rfft1b.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RFFT1B ( N, INC, R, LENR, WSAVE, LENSAV,
     1                  WORK, LENWRK, IER)
      INTEGER  N, INC, LENR, LENSAV, LENWRK, IER
      REAL     R(LENR), WSAVE(LENSAV)     ,WORK(LENWRK)
C
      IER = 0
C
      IF (LENR .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('RFFT1B ', 6)
      ELSEIF (LENSAV .LT. N + INT(LOG(REAL(N))) +4) THEN
        IER = 2
        CALL XERFFT ('RFFT1B ', 8)
      ELSEIF (LENWRK .LT. N) THEN
        IER = 3
        CALL XERFFT ('RFFT1B ', 10)
      ENDIF
C
      IF (N .EQ. 1) RETURN
C
      CALL RFFTB1 (N,INC,R,WORK,WSAVE,WSAVE(N+1))
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: rfftb1.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RFFTB1 (N,IN,C,CH,WA,FAC)
      REAL       CH(*), C(IN,*), WA(N) ,FAC(15)
C
      NF = FAC(2)
      NA = 0
      DO 10 K1=1,NF
      IP = FAC(K1+2)
      NA = 1-NA      
      IF(IP .LE. 5) GO TO 10
      IF(K1 .EQ. NF) GO TO 10
      NA = 1-NA
   10 CONTINUE 
      HALF = .5
      HALFM = -.5
      MODN = MOD(N,2)
      NL = N-2
      IF(MODN .NE. 0) NL = N-1
      IF (NA .EQ. 0) GO TO 120
      CH(1) = C(1,1)
      CH(N) = C(1,N)
      DO 118 J=2,NL,2
         CH(J) = HALF*C(1,J)
         CH(J+1) = HALFM*C(1,J+1)
  118 CONTINUE
      GO TO 124
  120 DO 122 J=2,NL,2
         C(1,J) = HALF*C(1,J)
         C(1,J+1) = HALFM*C(1,J+1)
  122 CONTINUE
  124 L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = FAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL R1F4KB (IDO,L1,C,IN,CH,1,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL R1F4KB (IDO,L1,CH,1,C,IN,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL R1F2KB (IDO,L1,C,IN,CH,1,WA(IW))
         GO TO 105
  104    CALL R1F2KB (IDO,L1,CH,1,C,IN,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
C rav    CALL RIF3KB (IDO,L1,C,IN,CH,1,WA(IW),WA(IX2))
         CALL R1F3KB (IDO,L1,C,IN,CH,1,WA(IW),WA(IX2))
         GO TO 108
  107    CALL R1F3KB (IDO,L1,CH,1,C,IN,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
         CALL R1F5KB (IDO,L1,C,IN,CH,1,WA(IW),WA(IX2),
     1                  WA(IX3),WA(IX4))
         GO TO 111
  110    CALL R1F5KB (IDO,L1,CH,1,C,IN,WA(IW),WA(IX2),
     1                  WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
C rav    CALL RIFGKB (IDO,IP,L1,IDL1,C,C,C,IN,CH,CH,1,WA(IW))
         CALL R1FGKB (IDO,IP,L1,IDL1,C,C,C,IN,CH,CH,1,WA(IW))
         GO TO 114
C rav 113    CALL RIFGKB (IDO,IP,L1,IDL1,CH,CH,CH,1,C,C,IN,WA(IW))
  113    CALL R1FGKB (IDO,IP,L1,IDL1,CH,CH,CH,1,C,C,IN,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: r1f2kb.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE R1F2KB (IDO,L1,CC,IN1,CH,IN2,WA1)
      REAL       CC(IN1,IDO,2,L1), CH(IN2,IDO,L1,2), WA1(IDO)
C
      DO 101 K=1,L1
         CH(1,1,K,1) = CC(1,1,1,K)+CC(1,IDO,2,K)
         CH(1,1,K,2) = CC(1,1,1,K)-CC(1,IDO,2,K)
 101  CONTINUE
      IF (IDO-2) 107,105,102
 102  IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            
            CH(1,I-1,K,1) = CC(1,I-1,1,K)+CC(1,IC-1,2,K)
            CH(1,I,K,1) = CC(1,I,1,K)-CC(1,IC,2,K)
            
            CH(1,I-1,K,2) = WA1(I-2)*(CC(1,I-1,1,K)-CC(1,IC-1,2,K))
     1           -WA1(I-1)*(CC(1,I,1,K)+CC(1,IC,2,K))
            CH(1,I,K,2) = WA1(I-2)*(CC(1,I,1,K)+CC(1,IC,2,K))+WA1(I-1)
     1           *(CC(1,I-1,1,K)-CC(1,IC-1,2,K))

 103     CONTINUE
 104  CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
 105  DO 106 K=1,L1
         CH(1,IDO,K,1) = CC(1,IDO,1,K)+CC(1,IDO,1,K)
         CH(1,IDO,K,2) = -(CC(1,1,2,K)+CC(1,1,2,K))
 106  CONTINUE
 107  RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: r1f3kb.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE R1F3KB (IDO,L1,CC,IN1,CH,IN2,WA1,WA2)
      REAL       CC(IN1,IDO,3,L1)    ,CH(IN2,IDO,L1,3),
     1           WA1(IDO)   ,WA2(IDO)
C
      ARG=2.*4.*ATAN(1.0)/3.
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
         CH(1,1,K,1) = CC(1,1,1,K)+2.*CC(1,IDO,2,K)
         CH(1,1,K,2) = CC(1,1,1,K)+(2.*TAUR)*CC(1,IDO,2,K)
     1   -(2.*TAUI)*CC(1,1,3,K)
         CH(1,1,K,3) = CC(1,1,1,K)+(2.*TAUR)*CC(1,IDO,2,K)
     1   +2.*TAUI*CC(1,1,3,K)
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
        CH(1,I-1,K,1) = CC(1,I-1,1,K)+(CC(1,I-1,3,K)+CC(1,IC-1,2,K))
        CH(1,I,K,1) = CC(1,I,1,K)+(CC(1,I,3,K)-CC(1,IC,2,K))
        CH(1,I-1,K,2) = WA1(I-2)*
     1 ((CC(1,I-1,1,K)+TAUR*(CC(1,I-1,3,K)+CC(1,IC-1,2,K)))-
     * (TAUI*(CC(1,I,3,K)+CC(1,IC,2,K))))
     2                   -WA1(I-1)*
     3 ((CC(1,I,1,K)+TAUR*(CC(1,I,3,K)-CC(1,IC,2,K)))+
     * (TAUI*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))))
            CH(1,I,K,2) = WA1(I-2)*
     4 ((CC(1,I,1,K)+TAUR*(CC(1,I,3,K)-CC(1,IC,2,K)))+
     8 (TAUI*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))))
     5                  +WA1(I-1)*
     6 ((CC(1,I-1,1,K)+TAUR*(CC(1,I-1,3,K)+CC(1,IC-1,2,K)))-
     8 (TAUI*(CC(1,I,3,K)+CC(1,IC,2,K))))
              CH(1,I-1,K,3) = WA2(I-2)*
     7 ((CC(1,I-1,1,K)+TAUR*(CC(1,I-1,3,K)+CC(1,IC-1,2,K)))+
     8 (TAUI*(CC(1,I,3,K)+CC(1,IC,2,K))))
     8   -WA2(I-1)*
     9 ((CC(1,I,1,K)+TAUR*(CC(1,I,3,K)-CC(1,IC,2,K)))-
     8 (TAUI*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))))
            CH(1,I,K,3) = WA2(I-2)*
     1 ((CC(1,I,1,K)+TAUR*(CC(1,I,3,K)-CC(1,IC,2,K)))-
     8 (TAUI*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))))
     2                 +WA2(I-1)*
     3 ((CC(1,I-1,1,K)+TAUR*(CC(1,I-1,3,K)+CC(1,IC-1,2,K)))+
     8 (TAUI*(CC(1,I,3,K)+CC(1,IC,2,K))))
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: r1f4kb.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE R1F4KB (IDO,L1,CC,IN1,CH,IN2,WA1,WA2,WA3)
      REAL       CC(IN1,IDO,4,L1)  ,CH(IN2,IDO,L1,4)    ,
     1           WA1(IDO)  ,        WA2(IDO)  ,       WA3(IDO)
C
      SQRT2=SQRT(2.)
      DO 101 K=1,L1
         CH(1,1,K,3) = (CC(1,1,1,K)+CC(1,IDO,4,K))
     1   -(CC(1,IDO,2,K)+CC(1,IDO,2,K))
         CH(1,1,K,1) = (CC(1,1,1,K)+CC(1,IDO,4,K))
     1   +(CC(1,IDO,2,K)+CC(1,IDO,2,K))
         CH(1,1,K,4) = (CC(1,1,1,K)-CC(1,IDO,4,K))
     1   +(CC(1,1,3,K)+CC(1,1,3,K))
         CH(1,1,K,2) = (CC(1,1,1,K)-CC(1,IDO,4,K))
     1   -(CC(1,1,3,K)+CC(1,1,3,K))
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
        CH(1,I-1,K,1) = (CC(1,I-1,1,K)+CC(1,IC-1,4,K))
     1  +(CC(1,I-1,3,K)+CC(1,IC-1,2,K))
        CH(1,I,K,1) = (CC(1,I,1,K)-CC(1,IC,4,K))
     1  +(CC(1,I,3,K)-CC(1,IC,2,K))
        CH(1,I-1,K,2)=WA1(I-2)*((CC(1,I-1,1,K)-CC(1,IC-1,4,K))
     1  -(CC(1,I,3,K)+CC(1,IC,2,K)))-WA1(I-1)
     1  *((CC(1,I,1,K)+CC(1,IC,4,K))+(CC(1,I-1,3,K)-CC(1,IC-1,2,K)))
        CH(1,I,K,2)=WA1(I-2)*((CC(1,I,1,K)+CC(1,IC,4,K))
     1  +(CC(1,I-1,3,K)-CC(1,IC-1,2,K)))+WA1(I-1)
     1  *((CC(1,I-1,1,K)-CC(1,IC-1,4,K))-(CC(1,I,3,K)+CC(1,IC,2,K)))
        CH(1,I-1,K,3)=WA2(I-2)*((CC(1,I-1,1,K)+CC(1,IC-1,4,K))
     1  -(CC(1,I-1,3,K)+CC(1,IC-1,2,K)))-WA2(I-1)
     1  *((CC(1,I,1,K)-CC(1,IC,4,K))-(CC(1,I,3,K)-CC(1,IC,2,K)))
        CH(1,I,K,3)=WA2(I-2)*((CC(1,I,1,K)-CC(1,IC,4,K))
     1  -(CC(1,I,3,K)-CC(1,IC,2,K)))+WA2(I-1)
     1  *((CC(1,I-1,1,K)+CC(1,IC-1,4,K))-(CC(1,I-1,3,K)
     1  +CC(1,IC-1,2,K)))
        CH(1,I-1,K,4)=WA3(I-2)*((CC(1,I-1,1,K)-CC(1,IC-1,4,K))
     1  +(CC(1,I,3,K)+CC(1,IC,2,K)))-WA3(I-1)
     1 *((CC(1,I,1,K)+CC(1,IC,4,K))-(CC(1,I-1,3,K)-CC(1,IC-1,2,K)))
        CH(1,I,K,4)=WA3(I-2)*((CC(1,I,1,K)+CC(1,IC,4,K))
     1  -(CC(1,I-1,3,K)-CC(1,IC-1,2,K)))+WA3(I-1)
     1  *((CC(1,I-1,1,K)-CC(1,IC-1,4,K))+(CC(1,I,3,K)+CC(1,IC,2,K)))
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         CH(1,IDO,K,1) = (CC(1,IDO,1,K)+CC(1,IDO,3,K))
     1   +(CC(1,IDO,1,K)+CC(1,IDO,3,K))
         CH(1,IDO,K,2) = SQRT2*((CC(1,IDO,1,K)-CC(1,IDO,3,K))
     1   -(CC(1,1,2,K)+CC(1,1,4,K)))
         CH(1,IDO,K,3) = (CC(1,1,4,K)-CC(1,1,2,K))
     1   +(CC(1,1,4,K)-CC(1,1,2,K))
         CH(1,IDO,K,4) = -SQRT2*((CC(1,IDO,1,K)-CC(1,IDO,3,K))
     1   +(CC(1,1,2,K)+CC(1,1,4,K)))
  106 CONTINUE
  107 RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: r1f5kb.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE R1F5KB (IDO,L1,CC,IN1,CH,IN2,
     1       WA1,WA2,WA3,WA4)
      REAL   CC(IN1,IDO,5,L1)    ,CH(IN2,IDO,L1,5),
     1       WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)     ,WA4(IDO)
C
      ARG=2.*4.*ATAN(1.0)/5.
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2.*ARG)
      TI12=SIN(2.*ARG)
      DO 101 K=1,L1
         CH(1,1,K,1) = CC(1,1,1,K)+2.*CC(1,IDO,2,K)+2.*CC(1,IDO,4,K)
         CH(1,1,K,2) = (CC(1,1,1,K)+TR11*2.*CC(1,IDO,2,K)
     1   +TR12*2.*CC(1,IDO,4,K))-(TI11*2.*CC(1,1,3,K)
     1   +TI12*2.*CC(1,1,5,K))
         CH(1,1,K,3) = (CC(1,1,1,K)+TR12*2.*CC(1,IDO,2,K)
     1   +TR11*2.*CC(1,IDO,4,K))-(TI12*2.*CC(1,1,3,K)
     1   -TI11*2.*CC(1,1,5,K))
         CH(1,1,K,4) = (CC(1,1,1,K)+TR12*2.*CC(1,IDO,2,K)
     1   +TR11*2.*CC(1,IDO,4,K))+(TI12*2.*CC(1,1,3,K)
     1   -TI11*2.*CC(1,1,5,K))
         CH(1,1,K,5) = (CC(1,1,1,K)+TR11*2.*CC(1,IDO,2,K)
     1   +TR12*2.*CC(1,IDO,4,K))+(TI11*2.*CC(1,1,3,K)
     1   +TI12*2.*CC(1,1,5,K))
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
        CH(1,I-1,K,1) = CC(1,I-1,1,K)+(CC(1,I-1,3,K)+CC(1,IC-1,2,K))
     1  +(CC(1,I-1,5,K)+CC(1,IC-1,4,K))
        CH(1,I,K,1) = CC(1,I,1,K)+(CC(1,I,3,K)-CC(1,IC,2,K))
     1  +(CC(1,I,5,K)-CC(1,IC,4,K))
        CH(1,I-1,K,2) = WA1(I-2)*((CC(1,I-1,1,K)+TR11*
     1  (CC(1,I-1,3,K)+CC(1,IC-1,2,K))+TR12
     1  *(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))-(TI11*(CC(1,I,3,K)
     1  +CC(1,IC,2,K))+TI12*(CC(1,I,5,K)+CC(1,IC,4,K))))
     1  -WA1(I-1)*((CC(1,I,1,K)+TR11*(CC(1,I,3,K)-CC(1,IC,2,K))
     1  +TR12*(CC(1,I,5,K)-CC(1,IC,4,K)))+(TI11*(CC(1,I-1,3,K)
     1  -CC(1,IC-1,2,K))+TI12*(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))
        CH(1,I,K,2) = WA1(I-2)*((CC(1,I,1,K)+TR11*(CC(1,I,3,K)
     1  -CC(1,IC,2,K))+TR12*(CC(1,I,5,K)-CC(1,IC,4,K)))
     1  +(TI11*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))+TI12
     1  *(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))+WA1(I-1)
     1  *((CC(1,I-1,1,K)+TR11*(CC(1,I-1,3,K)
     1  +CC(1,IC-1,2,K))+TR12*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))
     1  -(TI11*(CC(1,I,3,K)+CC(1,IC,2,K))+TI12
     1  *(CC(1,I,5,K)+CC(1,IC,4,K))))
        CH(1,I-1,K,3) = WA2(I-2)
     1  *((CC(1,I-1,1,K)+TR12*(CC(1,I-1,3,K)+CC(1,IC-1,2,K))
     1  +TR11*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))-(TI12*(CC(1,I,3,K)
     1  +CC(1,IC,2,K))-TI11*(CC(1,I,5,K)+CC(1,IC,4,K))))
     1 -WA2(I-1)
     1 *((CC(1,I,1,K)+TR12*(CC(1,I,3,K)-
     1  CC(1,IC,2,K))+TR11*(CC(1,I,5,K)-CC(1,IC,4,K)))
     1  +(TI12*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))-TI11
     1  *(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))
        CH(1,I,K,3) = WA2(I-2)
     1 *((CC(1,I,1,K)+TR12*(CC(1,I,3,K)-
     1  CC(1,IC,2,K))+TR11*(CC(1,I,5,K)-CC(1,IC,4,K)))
     1  +(TI12*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))-TI11
     1  *(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))
     1  +WA2(I-1)
     1  *((CC(1,I-1,1,K)+TR12*(CC(1,I-1,3,K)+CC(1,IC-1,2,K))
     1  +TR11*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))-(TI12*(CC(1,I,3,K)
     1  +CC(1,IC,2,K))-TI11*(CC(1,I,5,K)+CC(1,IC,4,K))))
        CH(1,I-1,K,4) = WA3(I-2)
     1  *((CC(1,I-1,1,K)+TR12*(CC(1,I-1,3,K)+CC(1,IC-1,2,K))
     1  +TR11*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))+(TI12*(CC(1,I,3,K)
     1  +CC(1,IC,2,K))-TI11*(CC(1,I,5,K)+CC(1,IC,4,K))))
     1  -WA3(I-1)
     1 *((CC(1,I,1,K)+TR12*(CC(1,I,3,K)-
     1  CC(1,IC,2,K))+TR11*(CC(1,I,5,K)-CC(1,IC,4,K)))
     1  -(TI12*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))-TI11
     1  *(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))
        CH(1,I,K,4) = WA3(I-2)
     1 *((CC(1,I,1,K)+TR12*(CC(1,I,3,K)-
     1  CC(1,IC,2,K))+TR11*(CC(1,I,5,K)-CC(1,IC,4,K)))
     1  -(TI12*(CC(1,I-1,3,K)-CC(1,IC-1,2,K))-TI11
     1  *(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))
     1  +WA3(I-1)
     1  *((CC(1,I-1,1,K)+TR12*(CC(1,I-1,3,K)+CC(1,IC-1,2,K))
     1  +TR11*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))+(TI12*(CC(1,I,3,K)
     1  +CC(1,IC,2,K))-TI11*(CC(1,I,5,K)+CC(1,IC,4,K))))
        CH(1,I-1,K,5) = WA4(I-2)
     1  *((CC(1,I-1,1,K)+TR11*(CC(1,I-1,3,K)+CC(1,IC-1,2,K))
     1  +TR12*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))+(TI11*(CC(1,I,3,K)
     1  +CC(1,IC,2,K))+TI12*(CC(1,I,5,K)+CC(1,IC,4,K))))
     1  -WA4(I-1)
     1  *((CC(1,I,1,K)+TR11*(CC(1,I,3,K)-CC(1,IC,2,K))
     1  +TR12*(CC(1,I,5,K)-CC(1,IC,4,K)))-(TI11*(CC(1,I-1,3,K)
     1  -CC(1,IC-1,2,K))+TI12*(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))
        CH(1,I,K,5) = WA4(I-2)
     1  *((CC(1,I,1,K)+TR11*(CC(1,I,3,K)-CC(1,IC,2,K))
     1  +TR12*(CC(1,I,5,K)-CC(1,IC,4,K)))-(TI11*(CC(1,I-1,3,K)
     1  -CC(1,IC-1,2,K))+TI12*(CC(1,I-1,5,K)-CC(1,IC-1,4,K))))
     1  +WA4(I-1)
     1  *((CC(1,I-1,1,K)+TR11*(CC(1,I-1,3,K)+CC(1,IC-1,2,K))
     1  +TR12*(CC(1,I-1,5,K)+CC(1,IC-1,4,K)))+(TI11*(CC(1,I,3,K)
     1  +CC(1,IC,2,K))+TI12*(CC(1,I,5,K)+CC(1,IC,4,K))))
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: r1fgkb.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE R1FGKB (IDO,IP,L1,IDL1,CC,C1,C2,IN1,
     1          CH,CH2,IN2,WA)
      REAL      CH(IN2,IDO,L1,IP)    ,CC(IN1,IDO,IP,L1) ,
     1          C1(IN1,IDO,L1,IP)    ,C2(IN1,IDL1,IP),
     2          CH2(IN2,IDL1,IP)     ,WA(IDO)
C
      TPI=2.*4.*ATAN(1.0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            CH(1,I,K,1) = CC(1,I,1,K)
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            CH(1,I,K,1) = CC(1,I,1,K)
  104    CONTINUE
  105 CONTINUE
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            CH(1,1,K,J) = CC(1,IDO,J2-2,K)+CC(1,IDO,J2-2,K)
            CH(1,1,K,JC) = CC(1,1,J2-1,K)+CC(1,1,J2-1,K)
 1007       CONTINUE
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
            DO 109 I=3,IDO,2
               IC = IDP2-I
               CH(1,I-1,K,J) = CC(1,I-1,2*J-1,K)+CC(1,IC-1,2*J-2,K)
               CH(1,I-1,K,JC) = CC(1,I-1,2*J-1,K)-CC(1,IC-1,2*J-2,K)
               CH(1,I,K,J) = CC(1,I,2*J-1,K)-CC(1,IC,2*J-2,K)
               CH(1,I,K,JC) = CC(1,I,2*J-1,K)+CC(1,IC,2*J-2,K)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               CH(1,I-1,K,J) = CC(1,I-1,2*J-1,K)+CC(1,IC-1,2*J-2,K)
               CH(1,I-1,K,JC) = CC(1,I-1,2*J-1,K)-CC(1,IC-1,2*J-2,K)
               CH(1,I,K,J) = CC(1,I,2*J-1,K)-CC(1,IC,2*J-2,K)
               CH(1,I,K,JC) = CC(1,I,2*J-1,K)+CC(1,IC,2*J-2,K)
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1.
      AI1 = 0.
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            C2(1,IK,L) = CH2(1,IK,1)+AR1*CH2(1,IK,2)
            C2(1,IK,LC) = AI1*CH2(1,IK,IP)
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               C2(1,IK,L) = C2(1,IK,L)+AR2*CH2(1,IK,J)
               C2(1,IK,LC) = C2(1,IK,LC)+AI2*CH2(1,IK,JC)
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            CH2(1,IK,1) = CH2(1,IK,1)+CH2(1,IK,J)
  121    CONTINUE
  122 CONTINUE
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            CH(1,1,K,J) = C1(1,1,K,J)-C1(1,1,K,JC)
            CH(1,1,K,JC) = C1(1,1,K,J)+C1(1,1,K,JC)
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
            DO 125 I=3,IDO,2
               CH(1,I-1,K,J) = C1(1,I-1,K,J)-C1(1,I,K,JC)
               CH(1,I-1,K,JC) = C1(1,I-1,K,J)+C1(1,I,K,JC)
               CH(1,I,K,J) = C1(1,I,K,J)+C1(1,I-1,K,JC)
               CH(1,I,K,JC) = C1(1,I,K,J)-C1(1,I-1,K,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               CH(1,I-1,K,J) = C1(1,I-1,K,J)-C1(1,I,K,JC)
               CH(1,I-1,K,JC) = C1(1,I-1,K,J)+C1(1,I,K,JC)
               CH(1,I,K,J) = C1(1,I,K,J)+C1(1,I-1,K,JC)
               CH(1,I,K,JC) = C1(1,I,K,J)-C1(1,I-1,K,JC)
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         C2(1,IK,1) = CH2(1,IK,1)
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            C1(1,1,K,J) = CH(1,1,K,J)
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               C1(1,I-1,K,J) = WA(IDIJ-1)*CH(1,I-1,K,J)-WA(IDIJ)*
     1          CH(1,I,K,J)
               C1(1,I,K,J) = WA(IDIJ-1)*CH(1,I,K,J)+WA(IDIJ)*
     1          CH(1,I-1,K,J)
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               C1(1,I-1,K,J) = WA(IDIJ-1)*CH(1,I-1,K,J)-WA(IDIJ)*
     1          CH(1,I,K,J)
               C1(1,I,K,J) = WA(IDIJ-1)*CH(1,I,K,J)+WA(IDIJ)*
     1          CH(1,I-1,K,J)
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 RETURN
      END
