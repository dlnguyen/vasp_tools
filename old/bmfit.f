c
c     Use the universal bonding curve to fit the data calculated by VASP
c     Ref. PRB 37, 6632 (1988)
c     modified by T. C. Leung  April 14, 2007
c     Unit   energy : eV , vol : ( A0 **3 )
c
      PROGRAM MINFI2
      IMPLICIT REAL*8 (A-H,O-Z)
      character*64 temp1      
      DIMENSION A(3,2),GUESS(3),PAR1(3),aa(3,3)
      DIMENSION TITLE(10)
      COMMON/FIT/ IFT,NP,V(20),E(20)
      DIMENSION B(20),EFIT(20)
      EXTERNAL FE
      open(11,file='POSCAR')
      read(11,*) temp1
      read(11,*) acel
      do i=1,3
      read(11,*) (aa(i,j),j=1,3)
      enddo
      vol=aa(1,1)*aa(2,2)*aa(3,2) + aa(1,2)*aa(2,3)*aa(3,1)
     &     + aa(1,3)*aa(2,1)*aa(3,2) - aa(1,3)*aa(2,2)*aa(3,1)
     &     - aa(1,2)*aa(2,1)*aa(3,3) - aa(1,1)*aa(2,3)*aa(3,2)

      open(1,file='curve.dat')
      open(2,file='point.dat')
c     READ(5,101) TITLE
  101 FORMAT(10A8)
      READ(5,*) (GUESS(I),I=1,3)
      do i = 1,3
      READ(5,*) (A(I,J),J=1,2)
      enddo
   19 FORMAT(2F10.5)
      READ(5,*) DEL
      READ(5,*) NDIV
   12 FORMAT(3F10.5)
   14 FORMAT(2I5)
      PRINT 23,(GUESS(I),I=1,3)
   23 FORMAT(5X,'INITIAL GUESS ',3F10.5)
      PRINT 25,DEL,NDIV
   25 FORMAT(5X,F10.5,5X,I5)
      READ(5,*) NP,IFT,NPLOT
      PRINT 11,NP,IFT
   11 FORMAT(1H ,'# OF DATA POINTS',I5,5X,'IFT =',I5)
    1 FORMAT(16I5)
      do i= 1,np
      READ(5,*) B(I),E(I)
      enddo
    2 FORMAT(2F10.5)
      READ(5,*) shift
      PRINT 3,(B(I),E(I),I=1,NP)
    3 FORMAT(1H ,2F12.5)
      TPI=2.*ACOS(-1.)
      
      DO 10 I=1,NP
      V(I)=B(I)
      E(I)=E(I)-shift
   10 CONTINUE
      A(1,1)=A(1,1)-shift
      A(1,2)=A(1,2)-shift
      GUESS(1)=GUESS(1)-shift
      REWIND 1
c      WRITE(2,1) NP
      WRITE(2,2) (B(I),E(I)-shift,I=1,NP)
  100 CONTINUE
      CALL MINA(FE,3,NDIV,DEL,A,GUESS,PAR1,EMIN,IERR)
      DO 20 I=1,NP
      EFIT(I)=FN(IFT,V(I),PAR1)
   20 CONTINUE
    5 FORMAT(1H ,1P8E12.5)
c     PRINT 101,TITLE
    4 FORMAT(1H ,8F15.8)
      DO 40 I=1,NP
      DE=EFIT(I)-E(I)
      PRINT 4,B(I),V(I),E(I),EFIT(I),DE
   40 CONTINUE
      BM=PAR1(3)*PAR1(2)*2.*147.
C     FAC=1.594*SQRT(0.75)
C     AC0=PAR1(3)/FAC
C     AC0=AC0**(1./3.)
C     AC1=AC0*.529
C     PRINT 5,PAR(1),PAR(2),PAR(3),BM,AC0,AC1
      PRINT 78,IERR
   78 FORMAT(5X,' IERR=  ',I5)
c
c     Unit   energy : Ry , vol : ( au **3 )           
c
      scal = 1.0
c
c     Unit   energy : eV , vol : ( A0 **3 )   (in VASP)
c
c     the unit of bulk modulii is 10^11 N/m^2
c
      vol_conv = 0.529**3
      ene_conv = 13.6
      scal = ene_conv/vol_conv
c
c     ***************************************
c
      BM=BM/scal
      PRINT 5,PAR1(1),PAR1(2),PAR1(3),BM,EMIN
      EB=-PAR1(1)
      AM=(1.5*PAR1(3)/TPI)**(1./3.)
      XL=SQRT(0.5*EB/PAR1(2))/(2.*TPI*AM*AM)
      BP=1.+(AM+AM)/(3.*XL)
      PRINT 5,EB,AM,XL,BP
      alat = (par1(3)/vol)**(1.0/3.0)
      write(6,1997) alat,par1(3),bm,emin
1997  format(/,' lattice const =',f10.5/' volume        =',f10.5/
     &   ' bulk modulus  =',f10.5/' err           =',f10.5)
      IF(NPLOT.LE.0) NPLOT=50
      VMIN=100000.
      VMAX=0.
      DO 180 I=1,NP
      IF(V(I).GT.VMAX) VMAX=V(I)
      IF(V(I).LT.VMIN) VMIN=V(I)
  180 CONTINUE
      VMAX=VMAX*1.05
      VMIN=VMIN*0.95
      DIV=(VMAX-VMIN)/FLOAT(NPLOT-1)
      XP=VMIN
      DO 181 I=1,NPLOT
      YP=FN(IFT,XP,PAR1)
      WRITE(1,2) XP,YP-shift
      XP=XP+DIV
  181 CONTINUE
      STOP
      END
      FUNCTION FE(GUESS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION GUESS(3),EE(20)
      COMMON/FIT/ IFT,NP,V(20),E(20)
      F=0.
      DO 20 I=1,NP
      EE(I)=FN(IFT,V(I),GUESS)
      F=F+(E(I)-EE(I))**2
   20 CONTINUE
      F=F/FLOAT(NP)
      FE=SQRT(F)
      RETURN
      END
      FUNCTION FN(IFT,V,P)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P(3)
      FPI=4.*ACOS(-1.)
      IF(IFT.NE.0) GO TO 10
C *** PARABOLA ***
      FN=P(1)+P(2)*(V-P(3))**2
      RETURN
   10 IF(IFT.NE.1) GO TO 20
C *** UNIVERSAL BONDING CURVE ***
      X=(3.*V/FPI)**(1./3.)
      X0=(3.*P(3)/FPI)**(1./3.)
      XL=SQRT(-0.5*P(1)/P(2))/(FPI*X0*X0)
      XX=(X-X0)/XL
      FN=P(1)*(1+XX)*EXP(-XX)
      RETURN
   20 STOP
      END
