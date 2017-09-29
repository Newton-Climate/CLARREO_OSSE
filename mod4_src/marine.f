      SUBROUTINE MARINE(VIS,MODEL,RHH,WS,WH,ICSTL,EXTC,ABSC,NL)
      REAL A1,A2,A3,AREL,EXT55,F,PISC,QE1,QE2,QE3,RATIO,                &
     &     RH,RHH,T1AV,T1QABS,T1QEXT,T1XV
      REAL T2AV,T2QABS,T2QEXT,T2XV,T3AV,T3QABS,T3QEXT,T3XV,TOTAL,VIS,WH,&
     &     WS,T1QEXTJ,T1QEXTJM1,T2QEXTJ,T2QEXTJM1,T3QEXTJ,T3QEXTJM1,R
      INTEGER I,ICSTL,J,MODEL,NL,JM1,IWAVLO,IWAVHI
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!        THIS SUBROUTINE DETERMINES AEROSOL EXT + ABS COEFFICIENTS
!          FOR THE NAVY MARITIME MODEL
!            CODED BY STU GATHMAN                  -  NRL

!        INPUTS-
!        WSS = CURRENT WIND SPEED (M/S)
!        WHH = 24 HOUR AVERAGE WIND SPEED (M/S)
!        RHH = RELATIVE HUMIDITY (PERCENTAGE)
!        VIS = METEOROLOGICAL RANGE (KM)
!        ICTL = AIR MASS CHARACTER  1 = OPEN OCEAN
!                      10 = STRONG CONTINENTAL INFLUENCE
!        MODEL = MODEL ATMOSPHERE

!        OUTPUTS-
!        EXTC = EXTINCTION COEFFICIENT (KM-1)
!        ABSC = ABSORPTION COEFFICIENT (KM-1)
      INCLUDE 'IFIL.h'

!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/MARAER/T1QEXT(NWAVLN,4),T2QEXT(NWAVLN,4),T3QEXT(NWAVLN,4), &
     &  T1QABS(NWAVLN,4),T2QABS(NWAVLN,4),T3QABS(NWAVLN,4),AREL(4)
      REAL WSPD(8),RHD(8),EXTC(MAER,MXWVLN),ABSC(MAER,MXWVLN)
      REAL VX0
      COMMON/EXTWAV/VX0(NWAVLN)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /EXTWAV/,/MARAER/
      EXTERNAL DEVCBD,EXTDTA,MARDTA
      DATA WSPD/6.9, 4.1, 4.1, 10.29, 6.69, 12.35, 7.2, 6.9/
      DATA RHD/80., 75.63, 76.2, 77.13, 75.24, 80.53, 45.89, 80./
      PISC = PI/1000.0
      IF(.NOT.LJMASS) WRITE(IPR,'(/A)')' MARINE AEROSOL MODEL USED'

!     CHECK LIMITS OF MODEL VALIDITY:
      IF(RHH.GT.0.)THEN
          RH=RHH
      ELSE
          RH=RHD(MODEL+1)
      ENDIF
      IF(RH.GT.98.)THEN
          RH=98.
      ELSEIF(RH.LT.50.)THEN
          RH=50.
      ENDIF
      IF(ICSTL.LT.1 .OR. ICSTL.GT.10) ICSTL = 3

!     FIND SIZE DISTRIBUTION PARAMETERS FROM METEOROLOGY INPUT:
      IF(WH.GT.20.)THEN
          WH=20.
      ELSEIF(WH.LE.0.)THEN
          IF(.NOT.LJMASS)WRITE(IPR,'(/A)')                              &
     &      ' WS NOT SPECIFIED; A DEFAULT VALUE IS USED.'
          WH=WSPD(MODEL+1)
      ENDIF
      IF(WS.GT.20.)THEN
          WS=20.
      ELSEIF(WS.LE.0.)THEN
          IF(.NOT.LJMASS)WRITE(IPR,'(/A)')                              &
     &      ' WH NOT SPECIFIED; A DEFAULT VALUE IS USED.'
          WS=WH
      ENDIF
      IF(.NOT.LJMASS) WRITE(IPR,'(3(/9X,A,F9.2,A),/9X,A,I3)')           &
     &  'WIND SPEED =',WS,' M/SEC',                                     &
     &  'WIND SPEED (24 HR AVERAGE) =',WH,' M/SEC',                     &
     &  'RELATIVE HUMIDITY =',RH,' PERCENT',                            &
     &  'AIRMASS CHARACTER =',ICSTL

!        F IS A RELATIVE HUMIDITY DEPENDENT GROWTH CORRECTION
!        TO THE ATTENUATION COEFFICIENT.

      F=((2.-RH/100.)/(6.*(1.-RH/100.)))**0.33333
      A1=2000.0*ICSTL*ICSTL
      A2 = AMAX1(5.866*(WH-2.2), 0.5)
      A3 = 10**(0.06*WS-2.8)

!     FIND EXTINCTION AT 0.55 MICRONS AND NORMALIZE TO 1.

!     INTERPOLATE FOR RELATIVE HUMIDITY


      IF(RH.LE.AREL(2))THEN
          RATIO=(RH-AREL(1))/(AREL(2)-AREL(1))
          JM1=1
          J=2
      ELSEIF(RH.LE.AREL(3))THEN
          RATIO=(RH-AREL(2))/(AREL(3)-AREL(2))
          JM1=2
          J=3
      ELSE
          RATIO=(RH-AREL(3))/(AREL(4)-AREL(3))
          JM1=3
          J=4
      ENDIF

      CALL HUNT(VX0,NWAVLN,0.55,IWAVLO)
      IWAVHI=IWAVLO+1
      R=(0.55-VX0(IWAVLO))/(VX0(IWAVHI)-VX0(IWAVLO))
      T1QEXTJ=(T1QEXT(IWAVHI,J)-T1QEXT(IWAVLO,J))*R+T1QEXT(IWAVLO,J)
      T2QEXTJ=(T2QEXT(IWAVHI,J)-T2QEXT(IWAVLO,J))*R+T2QEXT(IWAVLO,J)
      T3QEXTJ=(T3QEXT(IWAVHI,J)-T3QEXT(IWAVLO,J))*R+T3QEXT(IWAVLO,J)

      T1QEXTJM1=(T1QEXT(IWAVHI,JM1)-T1QEXT(IWAVLO,JM1))*R +             &
     &     T1QEXT(IWAVLO,JM1)
      T2QEXTJM1=(T2QEXT(IWAVHI,JM1)-T2QEXT(IWAVLO,JM1))*R +             &
     &     T2QEXT(IWAVLO,JM1)
      T3QEXTJM1=(T3QEXT(IWAVHI,JM1)-T3QEXT(IWAVLO,JM1))*R +             &
     &     T3QEXT(IWAVLO,JM1)

      QE1=T1QEXTJM1+(T1QEXTJ-T1QEXTJM1)*RATIO
      QE2=T2QEXTJM1+(T2QEXTJ-T2QEXTJM1)*RATIO
      QE3=T3QEXTJM1+(T3QEXTJ-T3QEXTJM1)*RATIO
      TOTAL = A1*10.**QE1 + A2*10.**QE2 + A3*10.**QE3
      EXT55=PISC*TOTAL/F

!     IF METEOROLOGICAL RANGE NOT SPECIFIED,FIND FROM METEOR DATA

      IF(VIS.LE.0.) VIS=3.912/(EXT55+0.01159)
      A1=A1/TOTAL
      A2=A2/TOTAL
      A3=A3/TOTAL

!     CALCULATE NORMALIZED ATTENUATION COEFFICIENTS

      DO I=1,NWAVLN
          T1XV = T1QEXT(I,JM1) + (T1QEXT(I,J) - T1QEXT(I,JM1))*RATIO
          T2XV = T2QEXT(I,JM1) + (T2QEXT(I,J) - T2QEXT(I,JM1))*RATIO
          T3XV = T3QEXT(I,JM1) + (T3QEXT(I,J) - T3QEXT(I,JM1))*RATIO
          T1AV = T1QABS(I,JM1) + (T1QABS(I,J) - T1QABS(I,JM1))*RATIO
          T2AV = T2QABS(I,JM1) + (T2QABS(I,J) - T2QABS(I,JM1))*RATIO
          T3AV = T3QABS(I,JM1) + (T3QABS(I,J) - T3QABS(I,JM1))*RATIO
          EXTC(NL,I)=A1*10**T1XV+A2*10**T2XV+A3*10**T3XV
          ABSC(NL,I)=A1*10**T1AV+A2*10**T2AV+A3*10**T3AV
      ENDDO
      IF(.NOT.LJMASS)WRITE(IPR,'(/9X,A,F10.2,A)')                       &
     &  ' METEOROLOGICAL RANGE =',VIS,' KM.'
      RETURN
      END
