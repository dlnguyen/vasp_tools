C TSO573.WOFZD.FORT (SMALL PROGRAM USED SOMETIMES IN SET UP PHASE)      00230023
CC                                                                      00000010
CC                                                                      00000020
CC     EVALUATES EXP(X*X)*ERFC(X)                                       00000030
CC                                                                      00000040
CC                                                                      00000050
CC                                                                      00000060
      FUNCTION EXERFC(X)                                                00000070
      PI=3.14159265                                                     00000080
      FAC=2./SQRT(PI)                                                   00000090
      ALOW=FAC/(X+SQRT(X*X+2.))                                         00000100
      AUPPER=FAC/(X+SQRT(X*X+4./PI))                                    00000110
      XX=0.                                                             00000120
      YY=X                                                              00000130
      CALL WOFZ(XX,YY,RE,UIM,IQ)                                        00000140
      IF(IQ.NE.1) WRITE(6,902) IQ                                       00000150
      EXERFC=RE                                                         00000160
      IF(IQ.NE.1) EXERFC=0.                                             00000170
      IF(IQ.NE.1) STOP                                                  00000180
CC    WRITE(6,901) XX,YY,RE,UIM,ALOW,AUPPER                             00000190
      RETURN                                                            00000200
  901 FORMAT(2X,6(E14.6,2X))                                            00000210
  902 FORMAT(2X,' VALUE NOT IN FIRST QUADRANT  IQ=',I5)                 00000220
      END                                                               00000230
      SUBROUTINE WOFZ(X,Y,RE,IM,IQ)                                     00000240
C                                                                       00000250
C     THIS PROCEDURE EVALUATES THE REAL AND IMAGINARY PART OF THE       00000260
C     FUNCTION W(Z) = EXP(-Z**2)ERFC(-IZ) FOR ARGUMENTS Z=X+IY IN THE   00000270
C     FIRST QUADRANT OF THE COMPLEX PLANE. THE ACCURACY IS 10 DECIMAL   00000280
C     PLACES AFTER THE DECIMAL POINT, OR BETTER. FOR THE UNDERLYING     00000290
C     ANALYSIS, SEE W. GAUTSCHI, "EFFICIENT COMPUTATION OF THE COMPLEX  00000300
C     ERROR FUNCTION", SIAM J. MATH. ANAL., (1969)                      00000310
C                                                                       00000320
C     ALGORITHM 363, COMM. ACM, 12(1969), 635                           00000330
C                                                                       00000340
      COMPLEX*16 Z,WZ                                                   00000350
      REAL*8 XR,YR,DRE,DIM                                              00000360
      REAL*8 H,H2,LAMBDA,R1,R2,S,S1,S2,T1,T2,C                          00000370
      REAL*4 IM                                                         00000380
      INTEGER CAPN                                                      00000390
      LOGICAL B                                                         00000400
C                                                                       00000410
C     BECAUSE THE CALCULATION IS ALWAYS DONE IN THE FIRST QUADRANT,     00000420
C     THE OTHER QUADRANTS ARE REACHED BY USING THE SYMMETRY RELATIONS   00000430
C                                                                       00000440
      IQ=1                                                              00000450
C     IF(X.GE.0..AND.Y.GE.0.) IQ=1                                      00000460
      IF(X.LT.0..AND.Y.GE.0.) IQ=2                                      00000470
      IF(X.LT.0..AND.Y.LT.0.) IQ=3                                      00000480
      IF(X.GE.0..AND.Y.LT.0.) IQ=4                                      00000490
C                                                                       00000500
      XR=ABS(X)                                                         00000510
      YR=ABS(Y)                                                         00000520
      Z=DCMPLX(XR,YR)                                                   00000530
C                                                                       00000540
      IF(.NOT.(YR.LT.4.29.AND.XR.LT.5.33)) GOTO 20                      00000550
C                                                                       00000560
      S=(1.-YR/4.29)*DSQRT(1.-XR*XR/28.41)                              00000570
      H=1.6*S                                                           00000580
      H2=2.*H                                                           00000590
      CAPN=6+23*S                                                       00000600
      LAMBDA=H2**CAPN                                                   00000610
      NU=9+21*S                                                         00000620
      GOTO 30                                                           00000630
   20 CONTINUE                                                          00000640
      H=0.                                                              00000650
      CAPN=0                                                            00000660
      NU=8                                                              00000670
   30 CONTINUE                                                          00000680
C                                                                       00000690
C     IN THE FOLLOWING STATEMENT, LAMBDA=0 COVERS THE UNDERFLOW         00000700
C     CASE WHEN H>0 IS VERY SMALL.                                      00000710
C                                                                       00000720
      B=(H.EQ.0..OR.LAMBDA.EQ.0.)                                       00000730
      R1=0.                                                             00000740
      R2=0.                                                             00000750
      S1=0.                                                             00000760
      S2=0.                                                             00000770
      NN=NU+1                                                           00000780
      DO50 I=1,NN                                                       00000790
      N=NU-I+1                                                          00000800
      NP1=N+1                                                           00000810
      T1=YR+H+NP1*R1                                                    00000820
      T2=XR-NP1*R2                                                      00000830
      C=.5/(T1*T1+T2*T2)                                                00000840
      R1=C*T1                                                           00000850
      R2=C*T2                                                           00000860
      IF(.NOT.(H.GT.0..AND.N.LE.CAPN)) GOTO 70                          00000870
      T1=LAMBDA+S1                                                      00000880
      S1=R1*T1-R2*S2                                                    00000890
      S2=R2*T1+R1*S2                                                    00000900
      LAMBDA=LAMBDA/H2                                                  00000910
   70 CONTINUE                                                          00000920
   50 CONTINUE                                                          00000930
C                                                                       00000940
      IF(YR.NE.0.) GOTO 40                                              00000950
      DRE=DEXP(-XR*XR)                                                  00000960
      GOTO 45                                                           00000970
   40 DRE=1.12837916709551*S1                                           00000980
      IF(B)  DRE=1.12837916709551*R1                                    00000990
   45 DIM=1.12837916709551*S2                                           00001000
      IF(B)  DIM=1.12837916709551*R2                                    00001010
C                                                                       00001020
C     RETURN REAL AND IMAGINARY PART OF W(Z) IN THE QUADRANT IQ         00001030
C                                                                       00001040
      GOTO(1,2,3,4), IQ                                                 00001050
C                                                                       00001060
    1 CONTINUE                                                          00001070
      RE=SNGL(DRE)                                                      00001080
      IM=SNGL(DIM)                                                      00001090
      GOTO 5                                                            00001100
    2 CONTINUE                                                          00001110
      WZ=DCMPLX(DRE,DIM)                                                00001120
      WZ=DCONJG(WZ)                                                     00001130
CCCC  RE=SNGL(DREAL(WZ))                                                00001140
CCCC  IM=SNGL(DIMAG(WZ))                                                00001150
      GOTO 5                                                            00001160
    3 CONTINUE                                                          00001170
      WZ=DCMPLX(DRE,DIM)                                                00001180
      WZ=2.*CDEXP(-(Z*Z))-WZ                                            00001190
CCCC  RE=SNGL(DREAL(WZ))                                                00001200
CCCC  IM=SNGL(DIMAG(WZ))                                                00001210
      GOTO 5                                                            00001220
    4 CONTINUE                                                          00001230
      WZ=DCMPLX(DRE,DIM)                                                00001240
      WZ=2.*CDEXP(-(Z*Z))-WZ                                            00001250
      WZ=DCONJG(WZ)                                                     00001260
CCCC  RE=SNGL(DREAL(WZ))                                                00001270
CCCC  IM=SNGL(DIMAG(WZ))                                                00001280
C                                                                       00001290
    5 CONTINUE                                                          00001300
C                                                                       00001310
      RETURN                                                            00001320
      END                                                               00001330
      FUNCTION ERF(XX)                                                  ERF.2
      IMPLICIT REAL*8(A-H,O-Z)
C                                                                       ADDRESS.2
C     SANDIA MATHEMATICAL PROGRAM LIBRARY                               ADDRESS.3
C     APPLIED MATHEMATICS DIVISION 2613                                 ADDRESS.4
C     SANDIA LABORATORIES                                               ADDRESS.5
C     ALBUQUERQUE, NEW MEXICO  87185                                    JUN0278.1
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978                     JUN0278.2
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.1
C                    ISSUED BY SANDIA LABORATORIES                     *MAR1378.2
C  *                   A PRIME CONTRACTOR TO THE                       *MAR1378.3
C  *                UNITED STATES DEPARTMENT OF ENERGY                 *MAR1378.4
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * * *MAR1378.5
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE   *MAR1378.6
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE      *MAR1378.7
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,    *MAR1378.8
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES  *MAR1378.9
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL      *MAR1378.10
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR     *MAR1378.11
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS      *MAR1378.12
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE          *MAR1378.13
C  * OWNED RIGHTS.                                                     *MAR1378.14
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.15
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS     *MAR1378.16
C  * PART IS SAND77-1441.                                              *MAR1378.17
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.18
C                                                                       ADDRESS.9
C     WRITTEN BY J.E. VOGEL FROM APPROXIMATIONS DERIVED BY W.J. CODY .  ERF.4
C                                                                       ERF.5
C     ABSTRACT                                                          ERF.6
C                                                                       ERF.7
C          ERF(X) COMPUTES 2.0/SQRT(PI) TIMES THE INTEGRAL FROM 0 TO X  ERF.8
C          OF EXP(-X**2). THIS IS DONE USING RATIONAL APPROXIMATIONS.   ERF.9
C          ELEVEN CORRECT SIGNIFICANT FIGURES ARE PROVIDED.             ERF.10
C                                                                       ERF.11
C     DESCRIPTION OF PARAMETERS                                         ERF.12
C                                                                       ERF.13
C          X MAY BE ANY REAL VALUE                                      ERF.14
C                                                                       ERF.15
C     ERF IS DOCUMENTED COMPLETELY IN SC-M-70-275                       ERF.16
C                                                                       ERF.18
      DIMENSION P1(4),Q1(4),P2(6),Q2(6),P3(4),Q3(4)                     ERF.19
      DATA (P1(I),I=1,4)/242.6679552305318,21.97926161829415,6.996383488ERF.20
     1619136,-3.560984370181539E-02/,                                   ERF.21
     2(Q1(I),I=1,4)/215.0588758698612,91.16490540451490,15.0827976304077ERF.22
     39,1.0/,                                                           ERF.23
     4(P2(I),I=1,6)/22.898992851659,26.094746956075,14.571898596926,4.26ERF.24
     577201070898,.56437160686381,-6.0858151959688E-06/,                ERF.25
     6(Q2(I),I=1,6)/22.898985749891,51.933570687552,50.273202863803,26.2ERF.26
     788795758761,7.5688482293618,1.0/,                                 ERF.27
     8(P3(I),I=1,4)/-1.21308276389978E-2,-.1199039552681460,-.2439110294ERF.28
     988626,-3.24319519277746E-2/,                                      ERF.29
     1(Q3(I),I=1,4)/4.30026643452770E-02,.489552441961437,1.437712279371ERF.30
     218,1.0/                                                           ERF.31
      DATA SQPI/.564189583547756/                                       ERF.32
C                                                                       MESS.2
C     CALL MLMESS(10HSANDIA7.2 )                                        MESS3.1
C                                                                       MESS.4
C                                                                       ERF.33
      X=ABS(XX)                                                         ERF.35
      IF(X.GT.6.0)GO TO 320                                             ERF.36
      X2=X*X                                                            ERF.37
      IF(X.GT.4.0)GO TO 300                                             ERF.38
      IF(X.GT..46875)GO TO 200                                          ERF.39
      A= X*(P1(1)+X2*(P1(2)+X2*(P1(3)+X2*P1(4))))                       ERF.40
      A=A/(Q1(1)+X2*(Q1(2)+X2*(Q1(3)+X2*Q1(4))))                        ERF.41
      IF(XX.LT.0.)A=-A                                                  ERF.42
      ERF=A                                                             ERF.43
      GO TO 400                                                         ERF.44
 200  A=EXP(-X2)*(P2(1)+X*(P2(2)+X*(P2(3)+X*(P2(4)+X*(P2(5)+X*P2(6))))))ERF.45
      A=A/(Q2(1)+X*(Q2(2)+X*(Q2(3)+X*(Q2(4)+X*(Q2(5)+X*Q2(6))))))       ERF.46
      ERF=SIGN((1.0-A),XX)                                              ERF.47
      GO TO 400                                                         ERF.48
 300  XI2=1./X2                                                         ERF.49
      R=XI2*(P3(1)+XI2*(P3(2)+XI2*(P3(3)+XI2*P3(4))))/(Q3(1)+XI2*(Q3(2)+ERF.50
     1XI2*(Q3(3)+XI2*Q3(4))))                                           ERF.51
      A=EXP(-X2)*(SQPI+R)/X                                             ERF.52
      ERF=SIGN((1.0-A),XX)                                              ERF.53
      GO TO 400                                                         ERF.54
 320  CONTINUE                                                          ERF.55
      ERF=XX/X                                                          ERF.56
 400  RETURN                                                            ERF.57
      END                                                               ERF.58
      FUNCTION ERFC(XX)                                                 ERF.2
      IMPLICIT REAL*8(A-H,O-Z)
C                                                                       ADDRESS.2
C     SANDIA MATHEMATICAL PROGRAM LIBRARY                               ADDRESS.3
C     APPLIED MATHEMATICS DIVISION 2613                                 ADDRESS.4
C     SANDIA LABORATORIES                                               ADDRESS.5
C     ALBUQUERQUE, NEW MEXICO  87185                                    JUN0278.1
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978                     JUN0278.2
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.1
C                    ISSUED BY SANDIA LABORATORIES                     *MAR1378.2
C  *                   A PRIME CONTRACTOR TO THE                       *MAR1378.3
C  *                UNITED STATES DEPARTMENT OF ENERGY                 *MAR1378.4
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * * *MAR1378.5
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE   *MAR1378.6
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE      *MAR1378.7
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,    *MAR1378.8
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES  *MAR1378.9
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL      *MAR1378.10
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR     *MAR1378.11
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS      *MAR1378.12
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE          *MAR1378.13
C  * OWNED RIGHTS.                                                     *MAR1378.14
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.15
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS     *MAR1378.16
C  * PART IS SAND77-1441.                                              *MAR1378.17
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.18
C                                                                       ADDRESS.9
C     WRITTEN BY J.E. VOGEL FROM APPROXIMATIONS DERIVED BY W.J. CODY .  ERF.4
C                                                                       ERF.5
C     ABSTRACT                                                          ERF.6
C                                                                       ERF.7
C          ERF(X) COMPUTES 2.0/SQRT(PI) TIMES THE INTEGRAL FROM 0 TO X  ERF.8
C          OF EXP(-X**2). THIS IS DONE USING RATIONAL APPROXIMATIONS.   ERF.9
C          ELEVEN CORRECT SIGNIFICANT FIGURES ARE PROVIDED.             ERF.10
C...........   ERFC COMPUTES  (1.0 - ERF)   ...............
C                                                                       ERF.11
C     DESCRIPTION OF PARAMETERS                                         ERF.12
C                                                                       ERF.13
C          X MAY BE ANY REAL VALUE                                      ERF.14
C                                                                       ERF.15
C     ERF IS DOCUMENTED COMPLETELY IN SC-M-70-275                       ERF.16
C                                                                       ERF.18
      DIMENSION P1(4),Q1(4),P2(6),Q2(6),P3(4),Q3(4)                     ERF.19
      DATA (P1(I),I=1,4)/242.6679552305318,21.97926161829415,6.996383488ERF.20
     1619136,-3.560984370181539E-02/,                                   ERF.21
     2(Q1(I),I=1,4)/215.0588758698612,91.16490540451490,15.0827976304077ERF.22
     39,1.0/,                                                           ERF.23
     4(P2(I),I=1,6)/22.898992851659,26.094746956075,14.571898596926,4.26ERF.24
     577201070898,.56437160686381,-6.0858151959688E-06/,                ERF.25
     6(Q2(I),I=1,6)/22.898985749891,51.933570687552,50.273202863803,26.2ERF.26
     788795758761,7.5688482293618,1.0/,                                 ERF.27
     8(P3(I),I=1,4)/-1.21308276389978E-2,-.1199039552681460,-.2439110294ERF.28
     988626,-3.24319519277746E-2/,                                      ERF.29
     1(Q3(I),I=1,4)/4.30026643452770E-02,.489552441961437,1.437712279371ERF.30
     218,1.0/                                                           ERF.31
      DATA SQPI/.564189583547756/                                       ERF.32
C                                                                       MESS.2
C     CALL MLMESS(10HSANDIA7.2 )                                        MESS3.1
C                                                                       MESS.4
C                                                                       ERF.33
      X=ABS(XX)                                                         ERF.35
      IF(X.GT.6.0)GO TO 320                                             ERF.36
      X2=X*X                                                            ERF.37
      IF(X.GT.4.0)GO TO 300                                             ERF.38
      IF(X.GT..46875)GO TO 200                                          ERF.39
      A= X*(P1(1)+X2*(P1(2)+X2*(P1(3)+X2*P1(4))))                       ERF.40
      A=A/(Q1(1)+X2*(Q1(2)+X2*(Q1(3)+X2*Q1(4))))                        ERFC.41
      IF(XX.LT.0.)A=-A                                                  ERFC.42
      ERFC=A                                                             ERFC.43
      GO TO 400                                                         ERFC.44
 200  A=EXP(-X2)*(P2(1)+X*(P2(2)+X*(P2(3)+X*(P2(4)+X*(P2(5)+X*P2(6))))))ERFC.45
      A=A/(Q2(1)+X*(Q2(2)+X*(Q2(3)+X*(Q2(4)+X*(Q2(5)+X*Q2(6))))))       ERFC.46
      ERFC=SIGN((1.0-A),XX)                                              ERFC.47
      GO TO 400                                                         ERFC.48
 300  XI2=1./X2                                                         ERFC.49
      R=XI2*(P3(1)+XI2*(P3(2)+XI2*(P3(3)+XI2*P3(4))))/(Q3(1)+XI2*(Q3(2)+ERFC.50
     1XI2*(Q3(3)+XI2*Q3(4))))                                           ERFC.51
      A=EXP(-X2)*(SQPI+R)/X                                             ERFC.52
      ERFC=SIGN((1.0-A),XX)                                              ERFC.53
      GO TO 400                                                         ERFC.54
 320  CONTINUE                                                          ERFC.55
      ERFC=XX/X                                                          ERFC.56
 400  ERFC=1.D0-ERFC                                                    ERFC.57
      RETURN
      END                                                               ERFC.58
      SUBROUTINE ERRCHK(NCHARS,NARRAY)                                  EXTERRCHK.2
C                                                                       EXTERRCHK.3
C     SANDIA MATHEMATICAL PROGRAM LIBRARY                               EXTERRCHK.4
C     APPLIED MATHEMATICS DIVISION 2642                                 EXTERRCHK.5
C     SANDIA LABORATORIES                                               EXTERRCHK.6
C     ALBUQUERQUE, NEW MEXICO 87115                                     EXTERRCHK.7
C                                                                       EXTERRCHK.8
C     SIMPLIFIED VERSION FOR STAND-ALONE USE.     APRIL 1977            EXTERRCHK.9
C                                                                       EXTERRCHK.10
C     ABSTRACT                                                          EXTERRCHK.11
C         THE ROUTINES ERRCHK, ERXSET, AND ERRGET TOGETHER PROVIDE      EXTERRCHK.12
C         A UNIFORM METHOD WITH SEVERAL OPTIONS FOR THE PROCESSING      EXTERRCHK.13
C         OF DIAGNOSTICS AND WARNING MESSAGES WHICH ORIGINATE           EXTERRCHK.14
C         IN THE MATHEMATICAL PROGRAM LIBRARY ROUTINES.                 EXTERRCHK.15
C         ERRCHK IS THE CENTRAL ROUTINE, WHICH ACTUALLY PROCESSES       EXTERRCHK.16
C         MESSAGES.                                                     EXTERRCHK.17
C                                                                       EXTERRCHK.18
C     DESCRIPTION OF ARGUMENTS                                          EXTERRCHK.19
C         NCHARS - NUMBER OF CHARACTERS IN HOLLERITH MESSAGE.           EXTERRCHK.20
C                  IF NCHARS IS NEGATED, ERRCHK WILL UNCONDITIONALLY    EXTERRCHK.21
C                  PRINT THE MESSAGE AND STOP EXECUTION.  OTHERWISE,    EXTERRCHK.22
C                  THE BEHAVIOR OF ERRCHK MAY BE CONTROLLED BY          EXTERRCHK.23
C                  AN APPROPRIATE CALL TO ERXSET.                       EXTERRCHK.24
C         NARRAY - NAME OF ARRAY OR VARIABLE CONTAINING THE MESSAGE,    EXTERRCHK.25
C                  OR ELSE A LITERAL HOLLERITH CONSTANT CONTAINING      EXTERRCHK.26
C                  THE MESSAGE.  BY CONVENTION, ALL MESSAGES SHOULD     EXTERRCHK.27
C                  BEGIN WITH *IN SUBNAM, ...*, WHERE SUBNAM IS THE     EXTERRCHK.28
C                  NAME OF THE ROUTINE CALLING ERRCHK.                  EXTERRCHK.29
C                                                                       EXTERRCHK.30
C     EXAMPLES                                                          EXTERRCHK.31
C         1. TO ALLOW CONTROL BY CALLING ERXSET, USE                    EXTERRCHK.32
C            CALL ERRCHK(30,30HIN QUAD, INVALID VALUE OF ERR.)          EXTERRCHK.33
C         2. TO UNCONDITIONALLY PRINT A MESSAGE AND STOP EXECUTION, USE EXTERRCHK.34
C            CALL ERRCHK(-30,30HIN QUAD, INVALID VALUE OF ERR.)         EXTERRCHK.35
C                                                                       EXTERRCHK.36
      DIMENSION NARRAY(32)                                              EXTERRCHK.37
C                                                                       MESS.2
C     CALL MLMESS(10HSANDIA7.2 )                                        MESS3.1
C                                                                       MESS.4
C                                                                       EXTERRCHK.38
      CALL ERRGET(NF,NT)                                                EXTERRCHK.39
C     IF ERRCHK WAS CALLED WITH NEGATIVE CHARACTER COUNT, SET FATAL FLAGEXTERRCHK.40
      IF (NCHARS.LT.0) NF = -1                                          EXTERRCHK.41
C     IF MESSAGES ARE TO BE SUPPRESSED, RETURN                          EXTERRCHK.42
      IF (NF.EQ.0) RETURN                                               EXTERRCHK.43
C     IF CHARACTER COUNT IS INVALID, STOP                               EXTERRCHK.44
      IF (NCHARS.EQ.0) PRINT 5                                          EXTERRCHK.45
    5 FORMAT(/31H ERRCHK WAS CALLED INCORRECTLY.)                       EXTERRCHK.46
      IF (NCHARS.EQ.0) STOP                                             EXTERRCHK.47
C     PRINT MESSAGE                                                     EXTERRCHK.48
      CALL ERRPRT(IABS(NCHARS),NARRAY)                                  EXTERRCHK.49
C     IF LAST MESSAGE, SAY SO                                           EXTERRCHK.50
      IF (NF.EQ.1) PRINT 10                                             EXTERRCHK.51
   10 FORMAT (30H ERRCHK MESSAGE LIMIT REACHED.)                        EXTERRCHK.52
C     PRINT TRACE-BACK IF ASKED TO                                      EXTERRCHK.53
C     IF ((NT.GT.0).OR.(NF.LT.0)) CALL SYSTEM ROUTINE FOR TRACEBACK     EXTERRCHK.54
C     DECREMENT MESSAGE COUNT                                           EXTERRCHK.55
      IF (NF.GT.0) NF = NF-1                                            EXTERRCHK.56
      CALL ERXSET(NF,NT)                                                EXTERRCHK.57
C     IF ALL IS WELL, RETURN                                            EXTERRCHK.58
      IF (NF.GE.0) RETURN                                               EXTERRCHK.59
C     IF THIS MESSAGE IS SUPPRESSABLE BY AN ERXSET CALL,                EXTERRCHK.60
C     THEN EXPLAIN ERXSET USAGE.                                        EXTERRCHK.61
      IF (NCHARS.GT.0) PRINT 15                                         EXTERRCHK.62
   15 FORMAT ('  *** NOTE ***                                         EXTERRCHK.63
     1     TO MAKE THE ERROR MESSAGE PRINTED ABOVE BE NONFATAL,         EXTERRCHK.64
     2     OR TO SUPPRESS THE MESSAGE COMPLETELY,                       EXTERRCHK.65
     3     INSERT AN APPROPRIATE CALL TO ERXSET                         EXTERRCHK.66
     4     AT THE START OF YOUR PROGRAM.                                EXTERRCHK.67
     5     FOR EXAMPLE, TO PRINT UP TO 10 NONFATAL WARNING MESSAGES, USEEXTERRCHK.68
     6              CALL ERXSET(10,0) '    )                              EXTERRCHK.69
      PRINT 20                                                          EXTERRCHK.70
   20 FORMAT (/28H PROGRAM ABORT DUE TO ERROR.)                         EXTERRCHK.71
      STOP                                                              EXTERRCHK.72
      END                                                               EXTERRCHK.73
      SUBROUTINE ONECHK(NCHARS,NARRAY)                                  EXTERRCHK.74
C                                                                       EXTERRCHK.75
C     ABSTRACT                                                          EXTERRCHK.76
C         ONECHK IS A COMPANION ROUTINE OF ERRCHK.  IT IS CALLED        EXTERRCHK.77
C         JUST LIKE ERRCHK, AND MESSAGES FROM IT MAY BE SUPPRESSED      EXTERRCHK.78
C         BY AN APPROPRIATE CALL TO ERXSET.  IT DIFFERS FROM ERRCHK     EXTERRCHK.79
C         IN THAT EACH CALL TO ONECHK WILL PRODUCE NO MORE THAN ONE     EXTERRCHK.80
C         PRINTED MESSAGE, REGARDLESS OF HOW MANY TIMES THAT CALL IS    EXTERRCHK.81
C         EXECUTED, AND ONECHK NEVER TERMINATES EXECUTION.              EXTERRCHK.82
C         ITS PURPOSE IS TO PROVIDE ONE-TIME-ONLY INFORMATIVE           EXTERRCHK.83
C         DIAGNOSTICS.                                                  EXTERRCHK.84
C                                                                       EXTERRCHK.85
C     DESCRIPTION OF ARGUMENTS                                          EXTERRCHK.86
C         NCHARS - NUMBER OF CHARACTERS IN THE MESSAGE.                 EXTERRCHK.87
C                  IF NEGATED, THE MESSAGE WILL BE PRINTED (ONCE) EVEN  EXTERRCHK.88
C                  IF NFATAL HAS BEEN SET TO 0 (SEE ERXSET).            EXTERRCHK.89
C         NARRAY - SAME AS IN ERRCHK                                    EXTERRCHK.90
C                                                                       EXTERRCHK.91
      DIMENSION NARRAY(32)                                              EXTERRCHK.92
      DATA NFLAG/4H.$,*/                                                EXTERRCHK.93
      IF (NARRAY(1).EQ.NFLAG) RETURN                                    EXTERRCHK.94
      CALL ERRGET(NF,NT)                                                EXTERRCHK.95
      IF ((NF.EQ.0).AND.(NCHARS.GT.0)) RETURN                           EXTERRCHK.96
      CALL ERRPRT (59,59HTHE FOLLOWING INFORMATIVE DIAGNOSTIC WILL APPEAEXTERRCHK.97
     1R ONLY ONCE.)                                                     EXTERRCHK.98
      CALL ERRPRT(IABS(NCHARS),NARRAY)                                  EXTERRCHK.99
      IF (NF.GT.0) NF = NF-1                                            EXTERRCHK.100
      CALL ERXSET(NF,NT)                                                EXTERRCHK.101
      NARRAY(1) = NFLAG                                                 EXTERRCHK.102
      END                                                               EXTERRCHK.103
      SUBROUTINE ERRPRT(NCHARS,NARRAY)                                  EXTERRCHK.104
C                                                                       EXTERRCHK.105
C     UTILITY ROUTINE TO SIMPLY PRINT THE HOLLERITH MESSAGE IN NARRAY,  EXTERRCHK.106
C     WHOSE LENGTH IS NCHARS CHARACTERS.                                EXTERRCHK.107
C                                                                       EXTERRCHK.108
      DIMENSION NARRAY(32)                                              EXTERRCHK.109
C                                                                       EXTERRCHK.110
C     NOTE - NCH MUST BE THE NUMBER OF HOLLERITH CHARACTERS STORED      EXTERRCHK.111
C     PER WORD.  IF NCH IS CHANGED, FORMAT 1 MUST ALSO BE               EXTERRCHK.112
C     CHANGED CORRESPONDINGLY.                                          EXTERRCHK.113
C                                                                       EXTERRCHK.114
      NCH =  4                                                          EXTERRCHK.115
C     FOR LINE PRINTERS, USE                                            EXTERRCHK.116
    1 FORMAT (1X,32A4 )                                                 EXTERRCHK.117
C     FOR DATA TERMINALS, USE                                           EXTERRCHK.118
C   1 FORMAT (1X,7A10)                                                  EXTERRCHK.119
      NWORDS = (NCHARS+NCH-1)/NCH                                       EXTERRCHK.120
      PRINT 1,(NARRAY(I),I=1,NWORDS)                                    EXTERRCHK.121
      RETURN                                                            EXTERRCHK.122
      END                                                               EXTERRCHK.123
      SUBROUTINE ERXSET(NFATAL,NTRACE)                                  EXTERRCHK.124
C                                                                       EXTERRCHK.125
C     ABSTRACT                                                          EXTERRCHK.126
C         ERXSET IS A COMPANION ROUTINE TO SUBROUTINE ERRCHK.           EXTERRCHK.127
C         ERXSET ASSIGNS THE VALUES OF NFATAL AND NTRACE RESPECTIVELY   EXTERRCHK.128
C         TO NF AND NT IN COMMON BLOCK MLBLK0 THEREBY SPECIFYING THE    EXTERRCHK.129
C         STATE OF THE OPTIONS WHICH CONTROL THE EXECUTION OF ERRCHK.   EXTERRCHK.130
C                                                                       EXTERRCHK.131
C     DESCRIPTION OF ARGUMENTS                                          EXTERRCHK.132
C         BOTH ARGUMENTS ARE INPUT ARGUMENTS OF DATA TYPE INTEGER.      EXTERRCHK.133
C         NFATAL - IS A FATAL-ERROR / MESSAGE-LIMIT FLAG. A NEGATIVE    EXTERRCHK.134
C                  VALUE DENOTES THAT DETECTED DIFFICULTIES ARE TO BE   EXTERRCHK.135
C                  TREATED AS FATAL ERRORS.  NONNEGATIVE MEANS NONFATAL.EXTERRCHK.136
C                  A NONNEGATIVE VALUE IS THE MAXIMUM NUMBER OF NONFATALEXTERRCHK.137
C                  WARNING MESSAGES WHICH WILL BE PRINTED BY ERRCHK,    EXTERRCHK.138
C                  AFTER WHICH NONFATAL MESSAGES WILL NOT BE PRINTED.   EXTERRCHK.139
C                  (DEFAULT VALUE IS -1.)                               EXTERRCHK.140
C         NTRACE - .GE.1 WILL CAUSE A TRACE-BACK TO BE GIVEN,           EXTERRCHK.141
C                        IF THIS FEATURE IS IMPLEMENTED ON THIS SYSTEM. EXTERRCHK.142
C                  .LE.0 WILL SUPPRESS ANY TRACE-BACK, EXCEPT FOR       EXTERRCHK.143
C                        CASES WHEN EXECUTION IS TERMINATED.            EXTERRCHK.144
C                  (DEFAULT VALUE IS 0.)                                EXTERRCHK.145
C                                                                       EXTERRCHK.146
C         *NOTE* -- SOME CALLS TO ERRCHK WILL CAUSE UNCONDITIONAL       EXTERRCHK.147
C         TERMINATION OF EXECUTION.  ERXSET HAS NO EFFECT ON SUCH CALLS.EXTERRCHK.148
C                                                                       EXTERRCHK.149
C     EXAMPLES                                                          EXTERRCHK.150
C         1. TO PRINT UP TO 100 MESSAGES AS NONFATAL WARNINGS USE       EXTERRCHK.151
C            CALL ERXSET(100,0)                                         EXTERRCHK.152
C         2. TO SUPPRESS ALL MATHLIB WARNING MESSAGES USE               EXTERRCHK.153
C            CALL ERXSET(0,0)                                           EXTERRCHK.154
C                                                                       EXTERRCHK.155
      CALL ERSTGT(0,NFATAL,NTRACE)                                      EXTERRCHK.156
      RETURN                                                            EXTERRCHK.157
      END                                                               EXTERRCHK.158
      SUBROUTINE ERRGET(NFATAL,NTRACE)                                  EXTERRCHK.159
C                                                                       EXTERRCHK.160
C     ABSTRACT                                                          EXTERRCHK.161
C         ERRGET IS A COMPANION ROUTINE TO SUBROUTINE ERRCHK.           EXTERRCHK.162
C         ERRGET ASSIGNS TO NFATAL AND NTRACE RESPECTIVELY THE VALUES   EXTERRCHK.163
C         OF NF AND NT IN COMMON BLOCK MLBLK0 THEREBY ASCERTAINING THE  EXTERRCHK.164
C         STATE OF THE OPTIONS WHICH CONTROL THE EXECUTION OF ERRCHK.   EXTERRCHK.165
C                                                                       EXTERRCHK.166
C     DESCRIPTION OF ARGUMENTS                                          EXTERRCHK.167
C         BOTH ARGUMENTS ARE OUTPUT ARGUMENTS OF DATA TYPE INTEGER.     EXTERRCHK.168
C         NFATAL - CURRENT VALUE OF NF (SEE DESCRIPTION OF ERXSET.)     EXTERRCHK.169
C         NTRACE - CURRENT VALUE OF NT (SEE DESCRIPTION OF ERXSET.)     EXTERRCHK.170
C                                                                       EXTERRCHK.171
      CALL ERSTGT(1,NFATAL,NTRACE)                                      EXTERRCHK.172
      RETURN                                                            EXTERRCHK.173
      END                                                               EXTERRCHK.174
      SUBROUTINE ERSTGT(K,NFATAL,NTRACE)                                EXTERRCHK.175
C                                                                       EXTERRCHK.176
C     THIS ROUTINE IS A SLAVE TO ERRGET AND ERRSET WHICH KEEPS          EXTERRCHK.177
C     THE FLAGS AS LOCAL VARIABLES.                                     EXTERRCHK.178
C                                                                       EXTERRCHK.179
C     *** IF LOCAL VARIABLES ARE NOT NORMALLY RETAINED BETWEEN          EXTERRCHK.180
C     CALLS ON THIS SYSTEM, THE VARIABLES LNF AND LNT CAN BE            EXTERRCHK.181
C     PLACED IN A COMMON BLOCK AND PRESET TO THE FOLLOWING              EXTERRCHK.182
C     VALUES IN THE MAIN PROGRAM.                                       EXTERRCHK.183
C                                                                       EXTERRCHK.184
      DATA LNF/-1/,LNT/0/                                               EXTERRCHK.185
      IF (K.LE.0) LNF = NFATAL                                          EXTERRCHK.186
      IF (K.LE.0) LNT = NTRACE                                          EXTERRCHK.187
      IF (K.GT.0) NFATAL = LNF                                          EXTERRCHK.188
      IF (K.GT.0) NTRACE = LNT                                          EXTERRCHK.189
      RETURN                                                            EXTERRCHK.190
      END                                                               EXTERRCHK.191
