      SUBROUTINE DEFALT(Z,P,T)

!     THIS SUBROUTINE INTERPOLATES PROFILES FROM THE 6 BUILT-IN MODEL
!     ATMOSPHERIC TO ALTITUDE Z.  THE JUNIT INDICES FROM /CARD1B/
!     INDICATE WHICH PROFILE SHOULD BE USED FOR EACH INTERPOLATION:

!                   JUNIT     MODEL ATMOSPHERE
!                     1          TROPICAL
!                     2          MID-LATITUDE SUMMER
!                     3          MID-LATITUDE WINTER
!                     4          HIGH-LAT SUMMER
!                     5          HIGH-LAT WINTER
!                     6          US STANDARD

!     DECLARE ARGUMENTS:
!       Z        INPUT ALTITUDE [KM].
!       P        OUTPUT PRESSURE FROM MODEL ATMOSPHERE.
!       T        OUTPUT TEMPERATURE FROM MODEL ATMOSPHERE.
      REAL Z,P,T

!     INCLUDE PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INTEGER JUNITP,JUNITT,JUNIT,JLOW
      REAL WMOL
      COMMON/CARD1B/JUNITP,JUNITT,JUNIT(13),WMOL(12),JLOW
      REAL ALT,PMLOG,TMATM,AMOL
      COMMON/MLATM/ALT(50),PMLOG(50,6),TMATM(50,6),AMOL(50,6,8)
      REAL AMOLX
      COMMON/MLATMX/AMOLX(NLAYX,NMOLX)
      REAL TRAC
      COMMON/TRAC/TRAC(50,21)
      REAL CO2RAT
      COMMON/CO2MIX/CO2RAT
      INTEGER JUNITX
      REAL WMOLX
      COMMON/CRD1BX/JUNITX,WMOLX(NMOLX)
      real sn2o,scfc11,scfc12
      common/larry/sn2o,scfc11,scfc12

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /MLATM/,/TRAC/,/MLATMX/
      EXTERNAL MLATMB,XMLATM

!     DECLARE FUNCTIONS:
      REAL EXPINT

!     DECLARE LOCAL VARIABLES:
      INTEGER I0,I1,I2,I3,K,KM7
      LOGICAL LINIT
      REAL Z0,Z1,Z2,Z3,FACTOR,X0,X1,X2,X3

!     DECLARE FUNCTIONS
      REAL FLG4PT

!     DATA:
!       LBOUND   NUMBER OF LAYER BOUNDARIES IN MODEL ATMOSPHERES.
      INTEGER LBOUND
      DATA LBOUND/50/

!     BEGIN CALCULATIONS:
      IF(Z.LE.ALT(2))THEN

!         INPUT ALTITUDE Z IS BELOW SECOND ALTITUDE.
          I1=1
          I2=2
      ELSEIF(Z.GE.ALT(LBOUND-1))THEN

!         INPUT ALTITUDE Z IS ABOVE SECOND TO LAST ALTITUDE.
          I1=LBOUND-1
          I2=LBOUND
      ELSE

!         INTERMEDIATE ALTITUDE.
          I2=3
          DO 10 I3=4,LBOUND
              IF(Z.LE.ALT(I2))GOTO20
              I2=I3
   10     CONTINUE
   20     CONTINUE
          I1=I2-1
          I0=I2-2

!         LINEAR AND LAGRANGE 4-POINT INTERPOLATION COEFFICIENTS:
          Z0=ALT(I0)
          Z1=ALT(I1)
          Z2=ALT(I2)
          Z3=ALT(I3)
          LINIT=.TRUE.

!         TEST PRESSURE INPUT FLAG.
          IF(JUNITP.LE.6)THEN
              X0=PMLOG(I0,JUNITP)
              X1=PMLOG(I1,JUNITP)
              X2=PMLOG(I2,JUNITP)
              X3=PMLOG(I3,JUNITP)
              P=FLG4PT(Z,Z0,Z1,Z2,Z3,X0,X1,X2,X3,LINIT)
              LINIT=.FALSE.
              P=EXP(P)
          ENDIF

!         TEST TEMPERATURE INPUT FLAG.
          IF(JUNITT.LE.6)THEN
              X0=TMATM(I0,JUNITT)
              X1=TMATM(I1,JUNITT)
              X2=TMATM(I2,JUNITT)
              X3=TMATM(I3,JUNITT)
              T=FLG4PT(Z,Z0,Z1,Z2,Z3,X0,X1,X2,X3,LINIT)
              LINIT=.FALSE.
          ENDIF

!         MOLECULAR DENSITIES VARYING WITH MODEL ATMOSPHERE:
          DO 30 K=1,7
              IF(JUNIT(K).LE.6)THEN
                  X0=AMOL(I0,JUNIT(K),K)
                  X1=AMOL(I1,JUNIT(K),K)
                  X2=AMOL(I2,JUNIT(K),K)
                  X3=AMOL(I3,JUNIT(K),K)
                  WMOL(K)=FLG4PT(Z,Z0,Z1,Z2,Z3,X0,X1,X2,X3,LINIT)
                  LINIT=.FALSE.

!                 ADJUST CO2 MIXING RATIO BASED ON CARD1A INPUT.
                  IF(K.EQ.2)WMOL(2)=CO2RAT*WMOL(2)
                  JUNIT(K)=10
              ENDIF
   30     CONTINUE

!         MOLECULAR DENSITIES CONSTANT WITH MODEL ATMOSPHERE:
          DO 40 K=8,NMOL
              IF(JUNIT(K).LE.6)THEN
                  KM7=K-7
                  X0=TRAC(I0,KM7)
                  X1=TRAC(I1,KM7)
                  X2=TRAC(I2,KM7)
                  X3=TRAC(I3,KM7)
                  WMOL(K)=FLG4PT(Z,Z0,Z1,Z2,Z3,X0,X1,X2,X3,LINIT)
                  LINIT=.FALSE.
                  JUNIT(K)=10
              ENDIF
   40     CONTINUE
          wmol(4)=sn2o*wmol(4)
          WMOL(12)=1000.*WMOL(12)

!         MOLECULAR DENSITIES FOR CFC SPECIES:
          IF(JUNITX.LE.6)THEN
              DO 50 K=1,NMOLX
                  X0=AMOLX(I0,K)
                  X1=AMOLX(I1,K)
                  X2=AMOLX(I2,K)
                  X3=AMOLX(I3,K)
                  WMOLX(K)=FLG4PT(Z,Z0,Z1,Z2,Z3,X0,X1,X2,X3,LINIT)
                  LINIT=.FALSE.
   50         CONTINUE
              wmolx(1)=scfc11*wmolx(1)
              wmolx(2)=scfc12*wmolx(2)
          ENDIF
          RETURN
      ENDIF

!     USE 2-POINT INTERPOLATION/EXTRAPOLATION NEAR ALTITUDE END POINTS:
      FACTOR=(Z-ALT(I1))/(ALT(I2)-ALT(I1))
      IF(JUNITP.LE.6)P=EXP(PMLOG(I1,JUNITP)                             &
     &  +FACTOR*(PMLOG(I2,JUNITP)-PMLOG(I1,JUNITP)))
      IF(JUNITT.LE.6)T=TMATM(I1,JUNITT)                                 &
     &  +FACTOR*(TMATM(I2,JUNITT)-TMATM(I1,JUNITT))

!     MOLECULAR DENSITIES VARYING WITH MODEL ATMOSPHERE:
      DO 60 K=1,7
          IF(JUNIT(K).LE.6)THEN
              WMOL(K)=                                                  &
     &          EXPINT(AMOL(I1,JUNIT(K),K),AMOL(I2,JUNIT(K),K),FACTOR)
              IF(K.EQ.2)WMOL(2)=CO2RAT*WMOL(2)
              JUNIT(K)=10
          ENDIF
   60 CONTINUE

!     MOLECULAR DENSITIES CONSTANT WITH MODEL ATMOSPHERE:
      DO 70 K=8,NMOL
          IF(JUNIT(K).LE.6)THEN
              WMOL(K)=EXPINT(TRAC(I1,K-7),TRAC(I2,K-7),FACTOR)
              JUNIT(K)=10
          ENDIF
   70 CONTINUE
      wmol(4)=sn2o*wmol(4)
      WMOL(12)=1000.*WMOL(12)

!     MOLECULAR DENSITIES FOR CFC SPECIES:
      IF(JUNITX.LE.6)THEN
          DO 80 K=1,NMOLX
              WMOLX(K)=EXPINT(AMOLX(I1,K),AMOLX(I2,K),FACTOR)
   80     CONTINUE
          wmolx(1)=scfc11*wmolx(1)
          wmolx(2)=scfc12*wmolx(2)
      ENDIF

!     INTERPOLATIONS COMPLETE.
      RETURN
      END
      REAL FUNCTION FLG4PT(Z,Z1,Z2,Z3,Z4,F1,F2,F3,F4,LINIT)

!     INCLUDED FILES:
      INCLUDE 'ERROR.h'

!     FUNCTION LAG4PT DOES A 4 POINT LAGRANGE INTERPOLATION FOR
!     Z BETWEEN Z2 AND Z3, INCLUSIVE.  THE ABSCISSAE MUST BE
!     MONOTIC (Z1<Z2<Z3<Z4 OR Z1>Z2>Z3>Z4) AND IF F(Z) HAS AN
!     EXTREMUM BETWEEN Z2 AND Z3, THEN A LINEAR INTERPOLATION
!     IS USED INSTEAD.  THE 4-POINT INTERPOLATION FORMULA IS:

!     F(Z)  =  COEF1 P1     COEF2 P2     COEF3 P3     COEF4 P4

!     WHERE

!                 F1              F2              F3              F4
!          P1 = ------ ;   P2 = ------ ;   P4 = ------ ;   P4 = ------
!               DEMON1          DENOM2          DENOM3          DENOM4

!          COEF1  =  (Z - Z2) (Z - Z3) (Z - Z4)

!          COEF2  =  (Z - Z1) (Z - Z3) (Z - Z4)

!          COEF3  =  (Z - Z1) (Z - Z2) (Z - Z4)

!          COEF4  =  (Z - Z1) (Z - Z2) (Z - Z3)

!          DENOM1  =   (Z1 - Z2) (Z1 - Z3) (Z1 - Z4)

!          DENOM2  =   (Z2 - Z1) (Z2 - Z3) (Z2 - Z4)

!          DENOM3  =   (Z3 - Z1) (Z3 - Z2) (Z3 - Z4)

!          DENOM4  =   (Z4 - Z1) (Z4 - Z2) (Z4 - Z3)

!     THE EXTREMA IN F(Z) ARE FOUND BY SOLVING

!           d F(Z)        2
!           ------  =  A Z   -  2 B Z  +  C  =  0
!             dZ

!     WHERE

!        A  =  3 (P1 + P2 + P3 + P4)

!        B  =  (Z2 + Z3 + Z4) P1  +  (Z1 + Z3 + Z4) P2

!           +  (Z1 + Z2 + Z4) P3  +  (Z1 + Z2 + Z3) P4

!        C  =  (Z2 Z3 + Z2 Z4 + Z3 Z4) P1  +  (Z1 Z3 + Z1 Z4 + Z3 Z4) P2

!           +  (Z1 Z2 + Z1 Z4 + Z2 Z4) P3  +  (Z1 Z2 + Z1 Z3 + Z2 Z3) P4

!     THE DISCRIMINANT OF THE QUADRATIC EQUATION FOR Z
!     IS EXPANDED IN QUADRATIC TERMS OF P:

!                           _           _
!          2               \           \
!         B  - 4 A C   =    |  Pi       |       Pj  COEFij
!                          /_          /_
!                           i      j=i and j>i

!     WHERE

!                     1           2           2             2
!          COEFii  =  -  { (Zj-Zk)   + (Zj-Zl)   +   (Zk-Zl)  }
!                     2

!                              2
!          COEFij  =  2 (Zk-Zl)  +  (Zi-Zk) (Zj-Zl)  +  (Zi-Zl) (Zj-Zk)

!     DECLARE INPUTS
!       Z        ABSCISSA OF DESIRED ORDINATE.
!       Z1       FIRST ABSCISSA
!       Z2       SECOND ABSCISSA
!       Z3       THIRD ABSCISSA
!       Z4       FOURTH ABSCISSA
!       F1       FIRST ORDINATE
!       F2       SECOND ORDINATE
!       F3       THIRD ORDINATE
!       F4       FOURTH ORDINATE
!       LINIT    FLAG, FALSE IF ABSCISSAE ARE UNCHANGED.
      REAL Z,Z1,Z2,Z3,Z4,F1,F2,F3,F4
      LOGICAL LINIT

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     SAVED VARIABLES.
      REAL BCOEF1,BCOEF2,BCOEF3,BCOEF4,COEF1,COEF2,COEF3,COEF4,FACTOR,  &
     &  DENOM1,DENOM2,DENOM3,DENOM4,COEF11,COEF22,COEF33,COEF44,COEF12, &
     &  COEF13,COEF14,COEF23,COEF24,COEF34,CCOEF1,CCOEF2,CCOEF3,CCOEF4
      SAVE BCOEF1,BCOEF2,BCOEF3,BCOEF4,COEF1,COEF2,COEF3,COEF4,FACTOR,  &
     &  DENOM1,DENOM2,DENOM3,DENOM4,COEF11,COEF22,COEF33,COEF44,COEF12, &
     &  COEF13,COEF14,COEF23,COEF24,COEF34,CCOEF1,CCOEF2,CCOEF3,CCOEF4

!     LOCAL VARIABLES.
      REAL A,B,ZEXTRM,DSCRIM,ROOT,ZDIFF1,ZDIFF2,ZDIFF3,ZDIFF4,          &
     &  ZDIF12,ZDIF13,ZDIF14,ZDIF23,ZDIF24,ZDIF34,P1,P2,P3,P4,          &
     &  ZD12SQ,ZD13SQ,ZD14SQ,ZD23SQ,ZD24SQ,ZD34SQ,STORE

!     DATA:
!       SMALL    SMALL NUMBER TOLERANCE
      REAL SMALL
      DATA SMALL/1.E-30/

!     NEW ABSCISSAE CHECK
      IF(LINIT)THEN

!         CHECK INPUTS
          IF(Z2.EQ.Z3)GOTO10
          IF(Z2.LT.Z3)THEN
              IF(Z1.GE.Z2 .OR. Z2.GT.Z .OR. Z.GT.Z3 .OR. Z3.GE.Z4)GOTO10
          ELSE
              IF(Z1.LE.Z2 .OR. Z2.LT.Z .OR. Z.LT.Z3 .OR. Z3.LE.Z4)GOTO10
          ENDIF

!         DEFINE CONSTANTS DEPENDENT ON Z, Z1, Z2, Z3 AND Z4 ONLY.
          ZDIFF1=Z-Z1
          ZDIFF2=Z-Z2
          ZDIFF3=Z-Z3
          ZDIFF4=Z-Z4
          COEF1=ZDIFF2*ZDIFF3*ZDIFF4
          COEF2=ZDIFF1*ZDIFF3*ZDIFF4
          COEF3=ZDIFF1*ZDIFF2*ZDIFF4
          COEF4=ZDIFF1*ZDIFF2*ZDIFF3
          ZDIF12=Z1-Z2
          ZDIF13=Z1-Z3
          ZDIF14=Z1-Z4
          ZDIF23=Z2-Z3
          ZDIF24=Z2-Z4
          ZDIF34=Z3-Z4
          DENOM1= ZDIF12*ZDIF13*ZDIF14
          DENOM2=-ZDIF12*ZDIF23*ZDIF24
          DENOM3= ZDIF13*ZDIF23*ZDIF34
          DENOM4=-ZDIF14*ZDIF24*ZDIF34

!         CONSTANTS USED IN EXTREMUM CHECK
          BCOEF1=Z2+Z3+Z4
          BCOEF2=Z3+Z4+Z1
          BCOEF3=Z4+Z1+Z2
          BCOEF4=Z1+Z2+Z3
          ZD12SQ=ZDIF12**2
          ZD13SQ=ZDIF13**2
          ZD14SQ=ZDIF14**2
          ZD23SQ=ZDIF23**2
          ZD24SQ=ZDIF24**2
          ZD34SQ=ZDIF34**2
          COEF11=.5*(ZD23SQ+ZD24SQ+ZD34SQ)
          COEF22=.5*(ZD13SQ+ZD14SQ+ZD34SQ)
          COEF33=.5*(ZD12SQ+ZD14SQ+ZD24SQ)
          COEF44=.5*(ZD12SQ+ZD13SQ+ZD23SQ)
          STORE=ZDIF13*ZDIF24+ZDIF14*ZDIF23
          COEF12=2*ZD34SQ+STORE
          COEF34=2*ZD12SQ+STORE
          STORE=ZDIF12*ZDIF34-ZDIF14*ZDIF23
          COEF13=2*ZD24SQ+STORE
          COEF24=2*ZD13SQ+STORE
          STORE=ZDIF12*ZDIF34+ZDIF13*ZDIF24
          COEF14=2*ZD23SQ-STORE
          COEF23=2*ZD14SQ-STORE
          CCOEF1=Z2*Z3+Z2*Z4+Z3*Z4
          CCOEF2=Z1*Z3+Z1*Z4+Z3*Z4
          CCOEF3=Z1*Z2+Z1*Z4+Z2*Z4
          CCOEF4=Z1*Z2+Z1*Z3+Z2*Z3

!         LINEAR INTERPOLATION FACTOR
          FACTOR=(Z-Z2)/(Z3-Z2)
      ENDIF

!     DETERMINE LOCATION OF EXTREMA
      P1=F1/DENOM1
      P2=F2/DENOM2
      P3=F3/DENOM3
      P4=F4/DENOM4
      A=3*(P1+P2+P3+P4)
      B=P1*BCOEF1+P2*BCOEF2+P3*BCOEF3+P4*BCOEF4
      IF(ABS(A).LT.SMALL)THEN
          IF(ABS(B).GT.SMALL)THEN

!             ONE EXTREMUM
              ZEXTRM=.5*(CCOEF1*P1+CCOEF2*P2+CCOEF3*P3+CCOEF4*P4)/B
              IF((ZEXTRM.GT.Z2 .AND. ZEXTRM.LT.Z3) .OR.                 &
     &           (ZEXTRM.LT.Z2 .AND. ZEXTRM.GT.Z3))THEN

!                 EXTREMUM BETWEEN Z2 AND Z3.  USE LINEAR INTERPOLATION.
                  FLG4PT=F2+FACTOR*(F3-F2)
                  RETURN
              ENDIF
          ENDIF
      ELSE
          DSCRIM=P1*(P1*COEF11+P2*COEF12+P3*COEF13+P4*COEF14)           &
     &          +P2*(P2*COEF22+P3*COEF23+P4*COEF24)                     &
     &          +P3*(P3*COEF33+P4*COEF34)+P4**2*COEF44
          IF(DSCRIM.GE.0.)THEN

!             TWO EXTREMUM
              ROOT=SQRT(DSCRIM)
              ZEXTRM=(B+ROOT)/A
              IF((ZEXTRM.GT.Z2 .AND. ZEXTRM.LT.Z3) .OR.                 &
     &           (ZEXTRM.LT.Z2 .AND. ZEXTRM.GT.Z3))THEN

!                 EXTREMUM BETWEEN Z2 AND Z3.  USE LINEAR INTERPOLATION.
                  FLG4PT=F2+FACTOR*(F3-F2)
                  RETURN
              ENDIF
              ZEXTRM=(B-ROOT)/A
              IF((ZEXTRM.GT.Z2 .AND. ZEXTRM.LT.Z3) .OR.                 &
     &           (ZEXTRM.LT.Z2 .AND. ZEXTRM.GT.Z3))THEN

!                 EXTREMUM BETWEEN Z2 AND Z3.  USE LINEAR INTERPOLATION.
                  FLG4PT=F2+FACTOR*(F3-F2)
                  RETURN
              ENDIF
          ENDIF
      ENDIF

!     NO EXTREMA BETWEEN Z2 AND Z3.  USE 4-POINT LAGRANGE INTERPOLATION.
      FLG4PT=COEF1*P1+COEF2*P2+COEF3*P3+COEF4*P4
      RETURN

!     INCORRECT INPUT ERROR.
   10 CONTINUE
      WRITE(IPR,'(/A,/(18X,A,F12.4))')                                  &
     &  ' ERROR in FLG4PT:  ABSCISSAE OUT OF ORDER.',                   &
     &  ' Z1 =',Z1,' Z2 =',Z2,' Z  =',Z,' Z3 =',Z3,' Z4 =',Z4
      IF(LJMASS)CALL WRTBUF(FATAL)
      STOP ' ERROR in FLG4PT:  ABSCISSAE OUT OF ORDER.'
      END
