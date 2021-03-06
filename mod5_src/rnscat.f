      SUBROUTINE RNSCAT(V,R,TT,PHASE,DIST,CSSA,ASYMR)
      INTEGER PHASE,DIST

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD
      REAL SC(3,4)

!          ARGUMENTS:

!          F = FREQUENCY (GHZ)
!          R = RAINFALL RATE (MM/HR)
!          TCEL = TEMPERATURE (DEGREES CELSIUS)
!          PHASE = PHASE PARAMETER (1=WATER, 2=ICE)
!          DIST = DROP SIZE DISTRIBUTION PARAMETER
!                     (1=MARSHALL-PALMER, 2=BEST)

!          RESULTS:

!          SC(1) = ABSORPTION COEFFICIENT (1/KM)
!          SC(2) = EXTINCTION COEFFICIENT (1/KM)
!          SC(I),I=3,NSC = LEGENDRE COEFFICIENTS #I-3  (NSC=10)
!          ERR = ERROR RETURN CODE: 0=NO ERROR, 1=BAD FREQUENCY,
!                2=BAD RAINFALL RATE, 3=BAD TEMPERATURE,
!                4=BAD PHASE PARAMETER, 5=BAD DROP SIZE DISTRIBUTION

!          THE INTERNAL DATA:
      REAL FR(9),TEMP(3)

!          FR(I),I=1,NF = TABULATED FREQUENCIES (GHZ)  (NF=9)
!          TEMP = TABULATED TEMPERATURES

!          THE BLOCK-DATA SECTION
      DATA NF/9/,NSC/4/,MAXI/3/
      DATA TK/273.15/,CMT0/1.0/,C7500/0.5/,G0/0.0/,G7500/0.85/
      DATA TEMP/-10.,0.,10./
      DATA FR/19.35,37.,50.3,89.5,100.,118.,130.,183.,231./

!      THIS SUBROUTINE REQUIRES FREQUENCIES IN GHZ
      NOPR = 0
!C    IF(IK .EQ. 1) NOPR = 1
!C    IF(IENT .GT.1) NOPR = 0
      F= V *29.97925
      FSAV=F
      RSAV=R
      INTKEY=0
!      CONVERT TEMP TO DEGREES CELSIUS
      TCEL=TT-TK
      TSAV=TCEL
!      FREQ RANGE OF DATA 19.35-231 GHZ IF LESS THAN 19.35
!      SET UP PARAMETERS FOR INTERPOLATION
      IF(F.LT.FR(1)) THEN
        FL=0.0
        FM=FR(1)
        INTKEY=1
        IF(NOPR .GT. 0) WRITE (IPR,801)
      END IF
!      IF MORE THAN 231 GHZ SET UP PARAMETERS FOR EXTRAPOLATION
       IF(F.GT.FR(NF)) THEN
         FL=FR(NF)
         FM=7500.
         INTKEY=2
         IF(NOPR .GT. 0) WRITE (IPR,801)
      END IF
!      TEMP RANGE OF DATA IS -10 TO +10 DEGREES CELCIUS
!      IF BELOW OR ABOVE EXTREME SET AND DO CALCULATIONS AT EXTREME
      IF (TCEL.LT.TEMP(1)) THEN
        TCEL=TEMP(1)
        IF(NOPR .GT. 0) WRITE (IPR,802)
      END IF

      IF (TCEL.GT.TEMP(3)) THEN
        TCEL=TEMP(3)
        IF(NOPR .GT. 0) WRITE (IPR,802)
      END IF

!      RAIN RATE OF DATA IS FOR 0-50 MM/HR
!      IF GT 50 TREAT CALCULATIONS AS IF 50 MM/HR WAS INPUT
      IF(R.GT.50) THEN
        R=50.
        IF(NOPR .GT. 0) WRITE (IPR,803)
      END IF

      KI=1
!             FIGURE OUT THE SECOND INDEX
   10 J=PHASE+2*DIST

!             GET THE TEMPERATURE INTERPOLATION PARAMETER ST
!             IF NEEDED AND AMEND THE SECOND INDEX
      IF(TCEL.LT.TEMP(2))THEN
          IF(TCEL.LE.TEMP(1))THEN
              ST=0.
          ELSE
              ST=(TCEL-TEMP(1))/(TEMP(2)-TEMP(1))
          ENDIF
      ELSE
          IF(TCEL.GE.TEMP(3))THEN
              ST=1.
          ELSE
              ST=(TCEL-TEMP(2))/(TEMP(3)-TEMP(2))
          ENDIF
      ENDIF

!             FIGURE OUT THE THIRD INDEX AND THE FREQUENCY INTERPOLATION
!             PARAMETER SF
      CALL BS(K,F,FR,NF,SF)

!             INITIALIZE SC
      DO 11 I=1,NSC
      SC(KI,I)=0.
   11 CONTINUE
      SC(KI,3)=1.

!             NOW DO THE CALCULATIONS

!             THE WATER CONTENT IS
      IF(DIST.EQ.1) THEN
      WC=.0889*R**.84
      ELSE
      WC=.067*R**.846
      END IF

!             FOR A TEMPERATURE DEPENDENT CASE, I.E.
      IF(J.LT.3) THEN
      S1=(1.-SF)*(1.-ST)
      S2=(1.-SF)*ST
      S3=SF*(1.-ST)
      S4=SF*ST
      DO 14 I=1,MAXI
      IF(I.LE.2) THEN
      ISC=I
      ELSE
      ISC=I+1
      END IF
      SC(KI,ISC)=S1*TAB(I,J,K,WC)+S2*TAB(I,J+1,K,WC)+                   &
     &             S3*TAB(I,J,K+1,WC)+S4*TAB(I,J+1,K+1,WC)
   14 CONTINUE

!             FOR A TEMPERATURE INDEPENDENT CASE
      ELSE
      S1=1.-SF
      S2=SF
      DO 17 I=1,MAXI
      IF(I.LE.2) THEN
      ISC=I
      ELSE
      ISC=I+1
      END IF
      SC(KI,ISC)=S1*TAB(I,J,K,WC)+S2*TAB(I,J,K+1,WC)
   17 CONTINUE
      END IF
      F=FSAV
      IF(INTKEY.EQ.3) GO TO 20
      IF(INTKEY.EQ.4) GO TO 30
      IF(INTKEY.EQ.0) THEN
        CSSA=SC(KI,1)/SC(KI,2)
        IF(CSSA.GT.1.0) CSSA=1.0
        ASYMR=SC(KI,4)/3.0
        F=FSAV
        R=RSAV
        TCEL=TSAV
      RETURN
      END IF
      IF(INTKEY.EQ.1) THEN
        INTKEY=3
        F=FM
        KI=2
      END IF
      IF(INTKEY.EQ.2) THEN
        INTKEY=4
        F=FL
        KI=3
      END IF
      GO TO 10
   20 CONTINUE
      FDIF=FM-F
      FTOT=FM-FL
      CM=SC(KI,1)/SC(KI,2)
      IF(CM.GT.1.0) CM=1.0
      CL=CMT0
      AM=SC(KI,4)/3.0
      AL=G0
      GO TO 40
   30 CONTINUE
      FDIF=FM-F
      FTOT=FM-FL
      CM=C7500
      CL=SC(KI,1)/SC(KI,2)
      IF(CL.GT.1.0) CL=1.0
      AM=G7500
      AL=SC(KI,4)/3.0
   40 CTOT=CM-CL
      CAMT=FDIF*CTOT/FTOT
      CSSA=CM-CAMT
      ATOT=AM-AL
      AAMT=FDIF*ATOT/FTOT
      ASYMR=AM-AAMT
      F=FSAV
      R=RSAV
      TCEL=TSAV
      RETURN
801   FORMAT(2X,'***  THE ASYMMETRY PARAMETER DUE TO RAIN IS BASED ON', &
     & 'DATA BETWEEN 19 AND 231 GHZ',                                   &
     & /2X,'***  EXTRAPOLATION IS USED FOR FREQUENCIES LOWER AND',      &
     & 'HIGHER THAN THIS RANGE')
802   FORMAT(2X,'***  TEMPERATURE RANGE OF DATA IS -10 TO +10 ',        &
     &'DEGREES CELSIUS',/2X,'***  BEYOND THESE VALUES IT IS ',          &
     &'TREATED AS IF AT THE EXTREMES')
803   FORMAT(2X,'***  RAIN RATES BETWEEN 0 AND 50 MM/HR ARE',           &
     &'WITHIN THIS DATA RANGE',/2X,'***  ABOVE THAT THE ASYMMETRY',     &
     &' PARAMETER IS CALCULATED FOR 50 MM/HR')
      END
