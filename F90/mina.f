      SUBROUTINE MINA(FN,NV,NDIV,DEL,A,GUESS,X,FOFX,IERR)               MINA.2
      IMPLICIT REAL*8 (A-H,O-Z)
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
C     ORIGINAL ROUTINE WAS H2 SAND MIN, BY Z. BEISINGER AND S. BELL     MINA.4
C     PRESENT VERSION BY R E JONES                                      MINA.5
C                                                                       MINA.6
C     ABSTRACT                                                          MINA.7
C        MINA FINDS AN APPROXIMATE MINIMUM OF A REAL FUNCTION OF        MINA.8
C        NV VARIABLES, GIVEN AN INITIAL ESTIMATE OF THE POSITION OF     MINA.9
C        THE MINIMUM AND RANGES FOR EACH OF THE VARIABLES.              MINA.10
C        MINA USES A SELECTIVE DIRECTED SEARCH OF A SURROUNDING         MINA.11
C        NV-DIMENSIONAL GRID OF POINTS TO FIND A DIRECTION IN WHICH     MINA.12
C        THE FUNCTION DECREASES.  IT THEN PROCEEDS IN THIS DIRECTION    MINA.13
C        AS FAR AS THE FUNCTION DECREASES, THEN DETERMINES A NEW        MINA.14
C        DIRECTION TO TRAVEL.  WHEN NO SUCH DIRECTION IS FOUND THE      MINA.15
C        SEARCH INCREMENT FACTOR IS DECREASED AND THE PROCESS           MINA.16
C        IS REPEATED.                                                   MINA.17
C                                                                       MINA.18
C     DESCRIPTION OF ARGUMENTS                                          MINA.19
C        THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST  MINA.20
C              A(NV,2), GUESS(NV), X(NV)                                MINA.21
C                                                                       MINA.22
C        INPUT--                                                        MINA.23
C        FN   - NAME OF FUNCTION OF NV VARIABLES TO BE MINIMIZED.       MINA.24
C               (THIS NAME MUST APPEAR IN AN EXTERNAL STATEMENT.)       MINA.25
C               FORM OF THE CALLING SEQUENCE MUST BE FUNCTION FN(X),    MINA.26
C               WHERE X IS AN ARRAY OF NV VARIABLE VALUES. THE          MINA.27
C               ORDERING OF THE VARIABLES IS ARBITRARY, EXCEPT          MINA.28
C               THAT IT MUST AGREE WITH THE ORDERING USED IN            MINA.29
C               ARRAYS A AND GUESS.                                     MINA.30
C        NV   - NUMBER OF VARIABLES.  (NV .GE. 1)                       MINA.31
C        NDIV - NUMBER OF REFINEMENTS OF THE SEARCH INCREMENTS TO USE.  MINA.32
C               AT EACH REFINEMENT, THE INCREMENT IN EACH DIMENSION     MINA.33
C               IS DIVIDED BY 10.  (USUALLY NDIV IS ABOUT 3 OR 4.)      MINA.34
C        DEL  - FRACTION OF VARIABLE RANGE (IN EACH DIMENSION) TO USE   MINA.35
C               AS THE INITIAL INCREMENT (IN THAT DIMENSION)            MINA.36
C        A    - ARRAY OF SEARCH BOUNDS, DIMENSIONED NV BY 2.            MINA.37
C               A(I,1) SHOULD BE THE LOWER BOUND OF THE I-TH VARIABLE.  MINA.38
C               A(I,2) SHOULD BE THE UPPER BOUND OF THE I-TH VARIABLE.  MINA.39
C        GUESS- ARRAY OF NV INITIAL VALUES.  GUESS(I) SHOULD BE THE     MINA.40
C               INITIAL VALUE TO USE FOR THE I-TH VARIABLE.             MINA.41
C                                                                       MINA.42
C        OUTPUT--                                                       MINA.43
C        X    - ARRAY (DIMENSIONED NV) GIVING THE VALUES OF THE         MINA.44
C               VARIABLES AT THE MINIMUM.  X(I) WILL BE THE VALUE       MINA.45
C               OF THE I-TH VARIABLE.                                   MINA.46
C        FOFX - FUNCTION VALUE AT THE MINIMUM                           MINA.47
C        IERR - A STATUS CODE                                           MINA.48
C              -NORMAL CODE                                             MINA.49
C               =1 MEANS THE SEARCH FOR A MINIMUM PROCEEDED FOR THE     MINA.50
C                  SPECIFIED NUMBER OF REFINEMENTS.                     MINA.51
C              -ABNORMAL CODES                                          MINA.52
C               =2 MEANS NV IS GREATER THAN 50                          MINA.53
C               =3 MEANS A RANGE MINIMUM IS GREATER THAN THE            MINA.54
C                  CORRESPONDING MAXIMUM                                MINA.55
C                                                                       MINA.57
      DIMENSION A(NV,2),GUESS(NV),X(NV)                                 MINA.58
      DIMENSION XNOW(50),XNEW(50),R(50)                                 MINA.59
C                                                                       MESS.2
C     CALL MLMESS(10HSANDIA7.2 )                                        MESS3.1
C                                                                       MESS.4
C                                                                       MINA.60
C     INITIALIZE                                                        MINA.61
C                                                                       MINA.62
      IERR = 1                                                          MINA.64
      IF (NV.LE.50) GO TO 2                                             MINA.65
      CALL ERRCHK(33,33HIN MINA  , NV IS GREATER THAN 50.)              MINA.66
      IERR = 2                                                          MINA.67
      RETURN                                                            MINA.68
    2 NX = NV                                                           MINA.69
      IDIV = 0                                                          MINA.70
      DO 5 I=1,NX                                                       MINA.71
      XNOW(I) = GUESS(I)                                                MINA.72
      IF (XNOW(I).LT.A(I,1)) XNOW(I) = A(I,1)                           MINA.73
      IF (XNOW(I).GT.A(I,2)) XNOW(I) = A(I,2)                           MINA.74
      IF (A(I,1)-A(I,2)) 5,5,4                                          MINA.75
    4 CALL ERRCHK(46,46HIN MINA  , RANGE MINIMUM GREATER THAN MAXIMUM.) MINA.76
      IERR = 3                                                          MINA.77
      RETURN                                                            MINA.78
    5 R(I) = A(I,2)-A(I,1)                                              MINA.79
      DELTA = DEL                                                       MINA.80
      IF (DELTA.LE.0.0) DELTA = 0.1                                     MINA.81
      FNOW = FN(XNOW)                                                   MINA.82
C                                                                       MINA.83
C     FIND NEW DIRECTION                                                MINA.84
C                                                                       MINA.85
    7 DO 8 I=1,NX                                                       MINA.86
    8 XNEW(I) = XNOW(I)                                                 MINA.87
      FOLD = FNOW                                                       MINA.88
   10 DO 40 I=1,NX                                                      MINA.89
      IF (XNOW(I).GE.A(I,2)) GO TO 20                                   MINA.90
      XNEW(I) = MIN(XNOW(I)+DELTA*R(I),A(I,2))
      FNEW = FN(XNEW)                                                   MINA.92
      IF (FNEW.LT.FNOW) GO TO 30                                        MINA.93
   20 IF (XNOW(I) .LE. A(I,1)) GO TO 25                                 MINA.94
      XNEW(I) = MAX(XNOW(I)-DELTA*R(I),A(I,1))
      FNEW = FN(XNEW)                                                   MINA.96
      IF (FNEW.LT.FNOW) GO TO 30                                        MINA.97
   25 XNEW(I) = XNOW(I)                                                 MINA.98
      GO TO 40                                                          MINA.99
   30 FNOW = FNEW                                                       MINA.100
   40 CONTINUE                                                          MINA.101
      ISTEP = 1                                                         MINA.102
C                                                                       MINA.103
C     REFINE IF NEEDED                                                  MINA.104
C                                                                       MINA.105
      IF (FNOW.LT.FOLD) GO TO 50                                        MINA.106
      IF (IDIV.GE.NDIV) GO TO 100                                       MINA.107
      DELTA = DELTA*0.2                                                 MINA.108
      IDIV = IDIV+1                                                     MINA.109
      GO TO 10                                                          MINA.110
C                                                                       MINA.111
C     TRY TO CONTINUE IN CHOSEN DIRECTION                               MINA.112
C                                                                       MINA.113
   50 ICHNG = 0                                                         MINA.114
      FAC = 1.0                                                         MINA.115
      IF ((ISTEP/10)*10.EQ.ISTEP) FAC = 2.0                             MINA.116
      DO 60 I=1,NX                                                      MINA.117
      DX = (XNEW(I)-XNOW(I))*FAC                                        MINA.118
      XNOW(I) = XNEW(I)                                                 MINA.119
      IF (DX) 52,54,56                                                  MINA.120
   52 XNEW(I) = MAX(XNOW(I)+DX,A(I,1))
      IF (XNEW(I).LT.XNOW(I)) ICHNG = 1                                 MINA.122
      GO TO 60                                                          MINA.123
   54 XNEW(I) = XNOW(I)                                                 MINA.124
      GO TO 60                                                          MINA.125
   56 XNEW(I) = MIN(XNOW(I)+DX,A(I,2))
      IF (XNEW(I).GT.XNOW(I)) ICHNG = 1                                 MINA.127
   60 CONTINUE                                                          MINA.128
      IF (ICHNG.EQ.0) GO TO 7                                           MINA.129
      FNEW = FN(XNEW)                                                   MINA.130
      IF (FNEW.GE.FNOW) GO TO 7                                         MINA.131
      FNOW = FNEW                                                       MINA.132
      ISTEP = ISTEP+1                                                   MINA.133
      GO TO 50                                                          MINA.134
C                                                                       MINA.135
C     RETURN ANSWERS                                                    MINA.136
C                                                                       MINA.137
  100 FOFX = FOLD                                                       MINA.138
      DO 110 I=1,NX                                                     MINA.139
  110 X(I) = XNOW(I)                                                    MINA.140
      RETURN                                                            MINA.141
      END                                                               MINA.142
      SUBROUTINE SIMIN (F,K,EPS,ANS,S,NEV,ICONT,Y)                      SIMIN.2
      IMPLICIT REAL*8 (A-H,O-Z)
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
C     ORIGINAL ROUTINE BY L F SHAMPINE, AS DESCRIBED IN REF.1 BELOW.    SIMIN.4
C     PREPARATION FOR MATH LIBRARY BY R E JONES.                        SIMIN.5
C                                                                       SIMIN.6
C     ABSTRACT                                                          SIMIN.7
C        SIMIN FINDS AN APPROXIMATE MINIMUM OF A REAL FUNCTION OF K     SIMIN.8
C        VARIABLES, GIVEN AN INITIAL ESTIMATE OF THE POSITION OF THE    SIMIN.9
C        MINIMUM. THE SIMPLEX METHOD IS USED. SEE REFERENCE 1 BELOW     SIMIN.10
C        FOR A FULL EXPLANATION OF THIS METHOD.  BRIEFLY, A SET OF      SIMIN.11
C        K+1 POINTS IN K-DIMENSIONAL SPACE IS CALLED A SIMPLEX.         SIMIN.12
C        THE MINIMIZATION PROCESS ITERATES BY REPLACING THE POINT       SIMIN.13
C        WITH THE LARGEST FUNCTION VALUE BY A NEW POINT WITH A          SIMIN.14
C        SMALLER FUNCTION VALUE.  ITERATION CONTINUES UNTIL ALL THE     SIMIN.15
C        POINTS CLUSTER SUFFICIENTLY CLOSE TO A MINIMUM.                SIMIN.16
C                                                                       SIMIN.17
C     REFERENCES                                                        SIMIN.18
C      1. L F SHAMPINE, A ROUTINE FOR UNCONSTRAINED OPTIMIZATION,       SIMIN.19
C         SC-TM-72130  OR  SC-RR-720657                                 SIMIN.20
C      2. J A NELDER AND R MEAD, A SIMPLEX METHOD FOR FUNCTION          SIMIN.21
C         MINIMIZATION, COMPUTER JOURNAL, 7(1965) 308-313               SIMIN.22
C                                                                       SIMIN.23
C     DESCRIPTION OF PARAMETERS                                         SIMIN.24
C      --INPUT--                                                        SIMIN.25
C        F  - NAME OF FUNCTION OF K VARIABLES TO BE MINIMIZED.          SIMIN.26
C             (THIS NAME MUST APPEAR IN AN EXTERNAL STATEMENT.)         SIMIN.27
C             FORM OF THE CALLING SEQUENCE MUST BE FUNCTION F(X),       SIMIN.28
C             WHERE X IS AN ARRAY OF K VARIABLES.                       SIMIN.29
C        K  - THE NUMBER OF VARIABLES.  K MUST BE AT LEAST 2.           SIMIN.30
C             NORMALLY K SHOULD BE LESS THAN ABOUT 10, AS SIMIN         SIMIN.31
C             BECOMES LESS EFFECTIVE FOR LARGER VALUES OF K.            SIMIN.32
C        EPS- THE CONVERGENCE CRITERION.  LET YAVG BE THE AVERAGE       SIMIN.33
C             VALUE OF THE FUNCTION F AT THE K+1 POINTS OF THE          SIMIN.34
C             SIMPLEX, AND LET R BE THEIR STANDARD ERROR.  (THAT IS,    SIMIN.35
C             THE ROOT-MEAN-SQUARE OF THE SET OF VALUES (Y(I)-YAVG),    SIMIN.36
C             WHERE Y(I) IS THE FUNCTION VALUE AT THE I-TH POINT OF     SIMIN.37
C             THE SIMPLEX.)  THEN--                                     SIMIN.38
C             IF EPS.GT.0, CONVERGENCE IS OBTAINED IF  R.LE.EPS.        SIMIN.39
C             IF EPS.LT.0, CONVERGENCE IS IF  R.LE.ABS(EPS*YAVG).       SIMIN.40
C             IF EPS=0, THE PROCESS WILL NOT CONVERGE BUT INSTEAD WILL  SIMIN.41
C             QUIT WHEN NEV FUNCTION EVALUATIONS HAVE BEEN USED.        SIMIN.42
C        ANS- AN ARRAY OF LENGTH K CONTAINING A GUESS FOR THE LOCATION  SIMIN.43
C             OF A MINIMUM OF F.                                        SIMIN.44
C        S  - A SCALE PARAMETER, WHICH MAY BE A SIMPLE VARIABLE OR AN   SIMIN.45
C             ARRAY OF LENGTH K.  USE OF AN ARRAY IS SIGNALLED BY       SIMIN.46
C             SETTING S(1) NEGATIVE.                                    SIMIN.47
C             -SIMPLE VARIABLE CASE.  HERE S IS THE LENGTH OF EACH      SIMIN.48
C             SIDE OF THE INITIAL SIMPLEX.  THUS, THE INITIAL SEARCH    SIMIN.49
C             RANGE IS THE SAME FOR ALL THE VARIABLES.                  SIMIN.50
C             -ARRAY CASE.  HERE THE LENGTH OF SIDE I OF THE INITIAL    SIMIN.51
C             SIMPLEX IS ABS(S(I)).  THUS, THE INITIAL SEARCH RANGE     SIMIN.52
C             MAY BE DIFFERENT FOR DIFFERENT VARIABLES.                 SIMIN.53
C             NOTE-- THE VALUE(S) USED FOR S ARE NOT VERY CRITICAL.     SIMIN.54
C             ANY REASONABLE GUESS SHOULD DO O.K.                       SIMIN.55
C        NEV- THE MAXIMUM NUMBER OF FUNCTION EVALUATIONS TO BE USED.    SIMIN.56
C             (THE ACTUAL NUMBER USED MAY EXCEED THIS SLIGHTLY SO THE   SIMIN.57
C             LAST SEARCH ITERATION MAY BE COMPLETED.)                  SIMIN.58
C        ICONT - ICONT SHOULD BE ZERO ON ANY CALL TO SIMIN WHICH        SIMIN.59
C             IS NOT A CONTINUATION OF A PREVIOUS CALL.                 SIMIN.60
C             IF ICONT=1 THE PROBLEM WILL BE CONTINUED.  IN THIS        SIMIN.61
C             CASE THE WORK ARRAY Y MUST BE THE SAME ARRAY THAT WAS     SIMIN.62
C             USED IN THE CALL THAT IS BEING CONTINUED (AND THE VALUES  SIMIN.63
C             IN IT MUST BE UNCHANGED).  THE REASON FOR THIS IS THAT    SIMIN.64
C             IF ICONT=1 THEN THE ARGUMENT S IS IGNORED AND THE SIMPLEX SIMIN.65
C             AND RELATED FUNCTION VALUES THAT WERE STORED IN ARRAY Y   SIMIN.66
C             DURING A PREVIOUS EXECUTION ARE USED TO CONTINUE THAT     SIMIN.67
C             PREVIOUS PROBLEM.                                         SIMIN.68
C        Y  - A WORK ARRAY CONTAINING AT LEAST K*K + 5*K + 1 WORDS.     SIMIN.69
C             IF ICONT=1 THIS MUST BE THE SAME ARRAY USED IN THE CALL   SIMIN.70
C             THAT IS BEING CONTINUED.                                  SIMIN.71
C      --OUTPUT--                                                       SIMIN.72
C        ANS- ANS WILL CONTAIN THE LOCATION OF THE POINT WITH THE       SIMIN.73
C             SMALLEST VALUE OF THE FUNCTION THAT WAS FOUND.            SIMIN.74
C        S  - IN THE SIMPLE VARIABLE CASE S WILL BE RETURNED AS THE     SIMIN.75
C             AVERAGE DISTANCE FROM THE VERTICES TO THE CENTROID OF     SIMIN.76
C             THE SIMPLEX.                                              SIMIN.77
C             IN THE ARRAY CASE S(I) WILL BE RETURNED AS THE AVERAGE    SIMIN.78
C             DISTANCE IN THE I-TH DIMENSION OF VERTICES FROM           SIMIN.79
C             THE CENTROID.  (S(1) WILL BE NEGATED.)                    SIMIN.80
C             NOTE-- THE VALUE(S) RETURNED IN S ARE USEFUL FOR          SIMIN.81
C             ASSESSING THE FLATNESS OF THE FUNCTION NEAR THE           SIMIN.82
C             MINIMUM.  THE LARGER THE VALUE OF S (FOR A GIVEN          SIMIN.83
C             VALUE OF EPS), THE FLATTER THE FUNCTION.                  SIMIN.84
C        NEV- NEV WILL BE THE COUNT OF THE ACTUAL NUMBER OF FUNCTION    SIMIN.85
C             EVALUATIONS USED.                                         SIMIN.86
C        Y  - WILL CONTAIN ALL DATA NEEDED TO CONTINUE THE MINIMIZATION SIMIN.87
C             SEARCH EFFICIENTLY IN A SUBSEQUENT CALL.                  SIMIN.88
C             NOTE -- THE FIRST K+1 ELEMENTS OF Y WILL CONTAIN THE      SIMIN.89
C             FUNCTION VALUES AT THE K+1 POINTS OF THE LATEST SIMPLEX.  SIMIN.90
C             THE NEXT K*(K+1) ELEMENTS OF Y WILL BE THE K+1 POINTS     SIMIN.91
C             OF THE SIMPLEX (IN EXACT CORRESPONDENSE TO THE ARRAY      SIMIN.92
C             P DISCUSSED IN REFERENCE 1 ABOVE).  THE REMAINING 3*K     SIMIN.93
C             WORDS ARE TEMPORARY WORKING STORAGE ONLY.                 SIMIN.94
C                                                                       SIMIN.96
      DIMENSION ANS(K),S(K),Y(1)                                        SIMIN.97
      EXTERNAL F                                                        SIMIN.98
C                                                                       MESS.2
C     CALL MLMESS(10HSANDIA7.2 )                                        MESS3.1
C                                                                       MESS.4
C                                                                       SIMIN.99
      IF(K.GE.2 .AND. S(1).NE.0.) GO TO 10                              SIMIN.101
      CALL ERRCHK(39,39HIN SIMIN , S(1)=0. OR K IS LESS THAN 2.)        SIMIN.102
      RETURN                                                            SIMIN.103
   10 IF (K.GT.50) CALL ERRCHK(31,31HIN SIMIN , K IS LARGER THAN 50.)   SIMIN.104
C                                                                       SIMIN.105
      IP = K+2                                                          SIMIN.106
      IC = IP+K*(K+1)                                                   SIMIN.107
      IR = IC+K                                                         SIMIN.108
      IRR = IR+K                                                        SIMIN.109
      CALL SIMINA(F,K,EPS,ANS,S,NEV,ICONT,Y,Y(IP),Y(IC),Y(IR),Y(IRR))   SIMIN.110
      RETURN                                                            SIMIN.111
      END                                                               SIMIN.112
      SUBROUTINE SIMINA(F,K,EPS,ANS,S,NEV,ICONT,Y,P,PC,PR,PRR)          SIMIN.113
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ANS(K),S(K),Y(3),P(K,3),PC(K),PR(K),PRR(K)              SIMIN.114
      DATA ALPHA,BETA,GAMMA/1.0,0.5,2.0/                                SIMIN.115
C                                                                       SIMIN.116
C     SIMINA IS A CORE MINIMIZATION ROUTINE CALLED ONLY BY SIMIN.       SIMIN.117
C                                                                       SIMIN.118
C     INITIALIZATION                                                    SIMIN.119
C                                                                       SIMIN.120
      KOUNT = 0                                                         SIMIN.122
      KK = K                                                            SIMIN.123
      IF (KK.LE.1) GO TO 99                                             SIMIN.124
      ONEK = 1.0/FLOAT(KK)                                              SIMIN.125
      KP1 = KK+1                                                        SIMIN.126
      ONEKP1 = 1.0/FLOAT(KP1)                                           SIMIN.127
      TOL = FLOAT(KP1)*EPS**2                                           SIMIN.128
C                                                                       SIMIN.129
C     INITIAL SIMPLEX                                                   SIMIN.130
C                                                                       SIMIN.131
      IF (ICONT.GE.1) GO TO 10                                          SIMIN.132
      IF (S(1)) 4,99,1                                                  SIMIN.133
    1 SKP1 = S(1)*ONEKP1                                                SIMIN.134
      DO 2 I=1,KP1                                                      SIMIN.135
      DO 2 J=1,KK                                                       SIMIN.136
    2 P(J,I) = ANS(J) - SKP1                                            SIMIN.137
      DO 3 J=1,KK                                                       SIMIN.138
    3 P(J,J+1) = P(J,J+1) + S(1)                                        SIMIN.139
      GO TO 7                                                           SIMIN.140
C                                                                       SIMIN.141
    4 DO 5 I=1,KP1                                                      SIMIN.142
      DO 5 J=1,KK                                                       SIMIN.143
    5 P(J,I) = ANS(J) - ABS(S(J))*ONEKP1                                SIMIN.144
      DO 6 J=1,KK                                                       SIMIN.145
    6 P(J,J+1) = P(J,J+1) + ABS(S(J))                                   SIMIN.146
C                                                                       SIMIN.147
C     FUNCTION VALUES FOR INITIAL SIMPLEX                               SIMIN.148
C                                                                       SIMIN.149
    7 I1 = 1                                                            SIMIN.150
      DO 8 I=1,KP1                                                      SIMIN.151
      Y(I) = F(P(1,I))                                                  SIMIN.152
      IF (Y(I).GT.Y(I1)) I1 = I                                         SIMIN.153
    8 CONTINUE                                                          SIMIN.154
      YANS = F(ANS)                                                     SIMIN.155
      KOUNT = KP1+1                                                     SIMIN.156
      IF (YANS.GE.Y(I1)) GO TO 10                                       SIMIN.157
      Y(I1) = YANS                                                      SIMIN.158
      DO 9 J=1,KK                                                       SIMIN.159
    9 P(J,I1) = ANS(J)                                                  SIMIN.160
C                                                                       SIMIN.161
C     RE-START / NEXT ITERATION                                         SIMIN.162
C     IF K.LT.0 VALUES IN THE P AND Y ARRAYS (AND ONLY THESE VALUES)    SIMIN.163
C     WILL NOT HAVE BEEN DEFINED IN THIS CALL.  THIS IS NON-ANSI USAGE. SIMIN.164
C                                                                       SIMIN.165
C     FIRST FIND LARGEST, SECOND LARGEST, AND SMALLEST FUNCTION VALUES. SIMIN.166
C                                                                       SIMIN.167
   10 I1 = 1                                                            SIMIN.168
      IL = 1                                                            SIMIN.169
      DO 12 I=2,KP1                                                     SIMIN.170
      IF (Y(I).LT.Y(IL)) IL = I                                         SIMIN.171
      IF (Y(I).GT.Y(I1)) I1 = I                                         SIMIN.172
   12 CONTINUE                                                          SIMIN.173
      I2 = IL                                                           SIMIN.174
      DO 13 I=1,KP1                                                     SIMIN.175
      IF (I.EQ.I1) GO TO 13                                             SIMIN.176
      IF (Y(I).GT.Y(I2)) I2 = I                                         SIMIN.177
   13 CONTINUE                                                          SIMIN.178
C                                                                       SIMIN.179
C     COMPUTE CENTROID, LEAVING OUT P(*,I1)                             SIMIN.180
C                                                                       SIMIN.181
      DO 15 J=1,KK                                                      SIMIN.182
      SUM = 0.0                                                         SIMIN.183
      DO 14 I=1,KP1                                                     SIMIN.184
      IF (I.EQ.I1) GO TO 14                                             SIMIN.185
      SUM = SUM + P(J,I)                                                SIMIN.186
   14 CONTINUE                                                          SIMIN.187
   15 PC(J) = SUM*ONEK                                                  SIMIN.188
C                                                                       SIMIN.189
C     FORM REFLECTED POINT AND TEST                                     SIMIN.190
C                                                                       SIMIN.191
      DO 20 J=1,KK                                                      SIMIN.192
   20 PR(J) = PC(J) + ALPHA*(PC(J)-P(J,I1))                             SIMIN.193
      YR = F(PR)                                                        SIMIN.194
      KOUNT = KOUNT+1                                                   SIMIN.195
      IF (YR.LT.Y(IL)) GO TO 30                                         SIMIN.196
      IF (YR.GE.Y(I2)) GO TO 40                                         SIMIN.197
C                                                                       SIMIN.198
C     ACCEPT REFLECTED POINT                                            SIMIN.199
C                                                                       SIMIN.200
   21 Y(I1) = YR                                                        SIMIN.201
      DO 22 J=1,KK                                                      SIMIN.202
   22 P(J,I1) = PR(J)                                                   SIMIN.203
      GO TO 60                                                          SIMIN.204
C                                                                       SIMIN.205
C     EXPAND IN FAVORABLE DIRECTION AND TEST                            SIMIN.206
C                                                                       SIMIN.207
   30 DO 31 J=1,KK                                                      SIMIN.208
   31 PRR(J) = PR(J) + GAMMA*(PR(J)-PC(J))                              SIMIN.209
      YRR = F(PRR)                                                      SIMIN.210
      KOUNT = KOUNT+1                                                   SIMIN.211
      IF (YRR.GE.YR) GO TO 21                                           SIMIN.212
C                                                                       SIMIN.213
C     ACCEPT EXPANDED POINT                                             SIMIN.214
C                                                                       SIMIN.215
      Y(I1) = YRR                                                       SIMIN.216
      DO 32 J=1,KK                                                      SIMIN.217
   32 P(J,I1) = PRR(J)                                                  SIMIN.218
      GO TO 60                                                          SIMIN.219
C                                                                       SIMIN.220
C     DECIDE WHETHER TO ACCEPT REFLECTED POINT.                         SIMIN.221
C                                                                       SIMIN.222
   40 IF (YR.GE.Y(I1)) GO TO 42                                         SIMIN.223
      Y(I1) = YR                                                        SIMIN.224
      DO 41 J=1,KK                                                      SIMIN.225
   41 P(J,I1) = PR(J)                                                   SIMIN.226
C                                                                       SIMIN.227
C     TRY CONTRACTION.                                                  SIMIN.228
C                                                                       SIMIN.229
   42 DO 43 J=1,KK                                                      SIMIN.230
   43 PR(J) = PC(J) + BETA*(P(J,I1)-PC(J))                              SIMIN.231
      YCT = F(PR)                                                       SIMIN.232
      KOUNT = KOUNT+1                                                   SIMIN.233
      IF (YCT.GT.Y(I1)) GO TO 50                                        SIMIN.234
      Y(I1) = YCT                                                       SIMIN.235
      DO 44 J=1,KK                                                      SIMIN.236
   44 P(J,I1) = PR(J)                                                   SIMIN.237
      GO TO 60                                                          SIMIN.238
C                                                                       SIMIN.239
C     ALL EFFORTS FAILED.  SHRINK THE SIMPLEX ABOUT BEST POINT.         SIMIN.240
C                                                                       SIMIN.241
   50 DO 52 I=1,KP1                                                     SIMIN.242
      IF (I.EQ.IL) GO TO 52                                             SIMIN.243
      DO 51 J=1,KK                                                      SIMIN.244
   51 P(J,I) = 0.5*(P(J,I)+P(J,IL))                                     SIMIN.245
      Y(I) = F(P(1,I))                                                  SIMIN.246
   52 CONTINUE                                                          SIMIN.247
      KOUNT = KOUNT+KP1                                                 SIMIN.248
C                                                                       SIMIN.249
C     CHECK FOR CONVERGENCE                                             SIMIN.250
C                                                                       SIMIN.251
   60 IF (KOUNT.GE.NEV) GO TO 65                                        SIMIN.252
      IF (EPS.EQ.0.0) GO TO 10                                          SIMIN.253
      SUM = 0.0                                                         SIMIN.254
      DO 61 I=1,KP1                                                     SIMIN.255
   61 SUM = SUM + Y(I)                                                  SIMIN.256
      YAVG = SUM*ONEKP1                                                 SIMIN.257
      SUM = 0.0                                                         SIMIN.258
      DO 62 I=1,KP1                                                     SIMIN.259
   62 SUM = SUM + (Y(I)-YAVG)**2                                        SIMIN.260
      IF (EPS) 64,63,63                                                 SIMIN.261
   63 IF (SUM-TOL) 65,65,10                                             SIMIN.262
   64 IF (SUM-TOL*ABS(YAVG)) 65,65,10                                   SIMIN.263
C                                                                       SIMIN.264
C     CONVERGENCE OBTAINED.                                             SIMIN.265
C                                                                       SIMIN.266
C     COMPUTE CENTROID                                                  SIMIN.267
C                                                                       SIMIN.268
   65 DO 68 J=1,KK                                                      SIMIN.269
      SUM = 0.0                                                         SIMIN.270
      DO 67 I=1,KP1                                                     SIMIN.271
   67 SUM = SUM+P(J,I)                                                  SIMIN.272
   68 PC(J) = SUM*ONEKP1                                                SIMIN.273
      IF (S(1)) 73,69,69                                                SIMIN.274
C                                                                       SIMIN.275
C     COMPUTE S(1) AS AVERAGE DISTANCE OF VERTICES FROM CENTROID.       SIMIN.276
C                                                                       SIMIN.277
   69 DIST = 0.0                                                        SIMIN.278
      DO 71 I=1,KP1                                                     SIMIN.279
      SUM = 0.0                                                         SIMIN.280
      DO 70 J=1,KK                                                      SIMIN.281
   70 SUM = SUM + (P(J,I)-PC(J))**2                                     SIMIN.282
   71 DIST = DIST + SQRT(SUM)                                           SIMIN.283
      S(1) = DIST*ONEKP1                                                SIMIN.284
      GO TO 80                                                          SIMIN.285
C                                                                       SIMIN.286
C     COMPUTE S(J) AS AVERAGE DISTANCE IN J-TH DIMENSION OF             SIMIN.287
C     VERTICES FROM THE CENTROID.                                       SIMIN.288
C                                                                       SIMIN.289
   73 DO 75 J=1,KK                                                      SIMIN.290
      SUM = 0.0                                                         SIMIN.291
      DO 74 I=1,KP1                                                     SIMIN.292
   74 SUM = SUM + ABS(P(J,I)-PC(J))                                     SIMIN.293
   75 S(J) = SUM*ONEKP1                                                 SIMIN.294
      S(1) = -S(1)                                                      SIMIN.295
C                                                                       SIMIN.296
C     RETURN P(*,IL) AS ANSWER                                          SIMIN.297
C                                                                       SIMIN.298
   80 IL = 1                                                            SIMIN.299
      DO 82 I=2,KP1                                                     SIMIN.300
      IF (Y(I).LT.Y(IL)) IL = I                                         SIMIN.301
   82 CONTINUE                                                          SIMIN.302
      DO 84 J=1,KK                                                      SIMIN.303
   84 ANS(J) = P(J,IL)                                                  SIMIN.304
      NEV = KOUNT                                                       SIMIN.305
      RETURN                                                            SIMIN.306
C                                                                       SIMIN.307
C     ERROR MESSAGE                                                     SIMIN.308
C                                                                       SIMIN.309
   99 CALL ERRCHK(39,39HIN SIMINA, S(1)=0. OR K IS LESS THAN 2.)        SIMIN.310
      RETURN                                                            SIMIN.311
      END                                                               SIMIN.312
