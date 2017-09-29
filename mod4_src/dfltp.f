      SUBROUTINE DFLTP(L,GNDALT)

!     THIS SUBROUTINE INTERPOLATES PROFILES FROM THE 6 BUILT-IN MODEL
!     ATMOSPHERES TO PRESSURE.  THE JUNIT INDICES FROM /CARD1B/
!     INDICATE WHICH PROFILE SHOULD BE USED FOR EACH INTERPOLATION:

!                   JUNIT     MODEL ATMOSPHERE
!                     1          TROPICAL
!                     2          MID-LATITUDE SUMMER
!                     3          MID-LATITUDE WINTER
!                     4          HIGH-LAT SUMMER
!                     5          HIGH-LAT WINTER
!                     6          US STANDARD

!     IF JUNITP IS 10 OR MORE, THE ALTITUDE IS DETERMINED
!     FROM THE HYDROSTATIC EQUATION.  FOR THE FIRST LEVEL,
!     THE ALTITUDE IS SET TO GNDALT.

!     ARGUMENTS:
!       L        LAYER INDEX.
!       GNDALT   GROUND ALTITUDE [KM].
      INTEGER L
      REAL GNDALT

!     PARAMETERS:
!       LBOUND   NUMBER OF LEVELS IN MODEL ATMOSPHERE PROFILES.
!       LBNDM1   NUMBER OF LAYERS IN MODEL ATMOSPHERE PROFILES.
      INTEGER LBOUND,LBNDM1
      PARAMETER(LBOUND=50,LBNDM1=LBOUND-1)
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /CARD1B/
      INTEGER JUNITP,JUNITT,JUNIT,JLOW
      REAL WMOL
      COMMON/CARD1B/JUNITP,JUNITT,JUNIT(13),WMOL(12),JLOW

!     /MLATM/
      REAL ALT,PMLOG,TMATM,AMOL
      COMMON/MLATM/ALT(LBOUND),PMLOG(LBOUND,6),                         &
     &  TMATM(LBOUND,6),AMOL(LBOUND,6,8)

!     /MLATMX/
      REAL AMOLX
      COMMON/MLATMX/AMOLX(NLAYX,NMOLX)

!     /TRAC/
      REAL TRAC
      COMMON/TRAC/TRAC(LBOUND,21)

!     /CO2MIX/
      REAL CO2RAT
      COMMON/CO2MIX/CO2RAT

!     /CRD1BX/
      INTEGER JUNITX
      REAL WMOLX
      COMMON/CRD1BX/JUNITX,WMOLX(NMOLX)

!     /MDATA/
      REAL P,T,WH,WCO2,WO,WN2O,WCO,WCH4,WO2
      COMMON/MDATA/P(LAYDIM),T(LAYDIM),WH(LAYDIM),WCO2(LAYDIM),         &
     &  WO(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),WCH4(LAYDIM),WO2(LAYDIM)

!     /MODEL/
!       ZM       PROFILE BOUNDARY ALTITUDES [KM].
!       PM       PROFILE BOUNDARY PRESSURES [MBAR].
!       TM       PROFILE BOUNDARY TEMPERATURES [K].
!       RFNDX    PROFILE BOUNDARY REFRACTIVITIES.
!       DENSTY   PROFILE BOUNDARY DENSITIES [UNITS DEPEND ON SPECIES].
!       LRHSET   FLAG, TRUE IF RELATIVE HUMIDITY CANNOT BE SCALED.
      REAL ZM,PM,TM,RFNDX,DENSTY
      LOGICAL LRHSET
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    &
     &  RFNDX(LAYDIM),DENSTY(MEXT,LAYDIM),LRHSET(LAYDIM)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /MLATM/,/TRAC/,/MLATMX/
      EXTERNAL MLATMB,XMLATM,DEVCBD

!     FUNCTIONS:
!       HYSTAT   INTEGRATES HYDROSTATIC EQUATION TO FIND ALTITUDE.
      REAL FLG4PT,HYSTAT

!     LOCAL VARIABLES:
      INTEGER I0,I1,I2,I3,J0,J1,J2,J3,K,KM7,LM1
      LOGICAL LINIT
      REAL PLOG,PBNDLO,PBNDHI,P0,P1,P2,P3,Y0,Y1,Y2,Y3,P0X,P1X,P2X,P3X

!     DEFINE LOG(PRESSURE/MBAR)
      IF(P(L).LE.0.)THEN
          WRITE(IPR,'(2A)')' Error:  Pressure must be speified',        &
     &      ' for each atmospheric level when MODEL is 8.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP ' Input pressure must be positive when MODEL is 8.'
      ELSEIF(L.GT.1)THEN
          LM1=L-1
          IF(P(L).GT.P(LM1))THEN
              WRITE(IPR,'(2A)')                                         &
     &          ' Error:  Pressure inversions are not allowed.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP ' Pressure inversions are not allowed.'
          ENDIF
      ENDIF
      PLOG=LOG(P(L))

!     TEMPERATURE:
      LINIT=.TRUE.
      IF(JUNITT.LE.6)THEN

!         TEMPERATURE DETERMINED FROM MODEL ATMOSPHERE:
          PBNDLO=PMLOG(2,JUNITT)
          PBNDHI=PMLOG(LBNDM1,JUNITT)
          IF(PLOG.GE.PBNDLO)THEN

!             INPUT PRESSURE IS ABOVE PRESSURE OF
!             SECOND ATMOSPHERIC LEVEL:
              T(L)=TMATM(2,JUNITT)
              T(L)=T(L)+(TMATM(1,JUNITT)-T(L))                          &
     &          *(PLOG-PBNDLO)/(PMLOG(1,JUNITT)-PBNDLO)
          ELSEIF(PLOG.LE.PBNDHI)THEN

!             INPUT PRESSURE IS BELOW ABOVE PRESSURE OF
!             SECOND-TO-LAST ATMOSPHERIC LEVEL:
              T(L)=TMATM(LBNDM1,JUNITT)
              T(L)=T(L)+(TMATM(LBOUND,JUNITT)-T(L))                     &
     &          *(PLOG-PBNDHI)/(PMLOG(LBOUND,JUNITT)-PBNDHI)
          ELSE

!             INTERMEDIATE PRESSURE:
              DO 10 I2=2,LBNDM1
                  IF(PLOG.GE.PMLOG(I2,JUNITT))GOTO20
   10         CONTINUE
   20         CONTINUE
              I0=I2-2
              I1=I2-1
              I3=I2+1
              P0=PMLOG(I0,JUNITT)
              P1=PMLOG(I1,JUNITT)
              P2=PMLOG(I2,JUNITT)
              P3=PMLOG(I3,JUNITT)
              Y0=TMATM(I0,JUNITT)
              Y1=TMATM(I1,JUNITT)
              Y2=TMATM(I2,JUNITT)
              Y3=TMATM(I3,JUNITT)
              T(L)=FLG4PT(PLOG,P0,P1,P2,P3,Y0,Y1,Y2,Y3,LINIT)
              IF(JUNITP.EQ.JUNITT)LINIT=.FALSE.
          ENDIF
      ENDIF

!     ALTITUDE:
      IF(JUNITP.GE.1 .AND. JUNITP.LE.6)THEN

!         ALTITUDE DETERMINED FROM MODEL ATMOSPHERE:
          PBNDLO=PMLOG(2,JUNITP)
          PBNDHI=PMLOG(LBNDM1,JUNITP)
          IF(PLOG.GE.PBNDLO)THEN

!             INPUT PRESSURE IS ABOVE PRESSURE OF
!             SECOND ATMOSPHERIC LEVEL:
              ZM(L)=ALT(2)+(ALT(1)-ALT(2))                              &
     &          *(PLOG-PBNDLO)/(PMLOG(1,JUNITP)-PBNDLO)
          ELSEIF(PLOG.LE.PBNDHI)THEN

!             INPUT PRESSURE IS BELOW ABOVE PRESSURE OF
!             SECOND-TO-LAST ATMOSPHERIC LEVEL:
              ZM(L)=ALT(LBNDM1)+(ALT(LBOUND)-ALT(LBNDM1))               &
     &          *(PLOG-PBNDHI)/(PMLOG(LBOUND,JUNITP)-PBNDHI)
          ELSE
              IF(LINIT)THEN

!                 P AND T FROM DIFFERENT MODEL ATMOSPHERES:
                  DO 30 I2=2,LBNDM1
                      IF(PLOG.GE.PMLOG(I2,JUNITP))GOTO40
   30             CONTINUE
   40             CONTINUE
                  I0=I2-2
                  I1=I2-1
                  I3=I2+1
                  P0=PMLOG(I0,JUNITP)
                  P1=PMLOG(I1,JUNITP)
                  P2=PMLOG(I2,JUNITP)
                  P3=PMLOG(I3,JUNITP)
              ENDIF
              Y0=ALT(I0)
              Y1=ALT(I1)
              Y2=ALT(I2)
              Y3=ALT(I3)
              ZM(L)=FLG4PT(PLOG,P0,P1,P2,P3,Y0,Y1,Y2,Y3,LINIT)
              LINIT=.FALSE.
          ENDIF
      ELSEIF(L.EQ.1)THEN

!         SET ALTITUDE AT P(1) TO INPUT GNDALT:
          ZM(1)=GNDALT
      ELSE

!         INTEGRATE THE HYDROSTATIC EQUATION ASSUMING
!         TEMPERATURE VARIES LINEARLY WITH LOG(PRESSURE):
          ZM(L)=HYSTAT(ZM(LM1),T(LM1),P(LM1),T(L),P(L))
      ENDIF

!     MOLECULAR DENSITIES VARYING WITH MODEL ATMOSPHERE:
      DO 70 K=1,7
          IF(JUNIT(K).LE.6)THEN
              PBNDLO=PMLOG(2,JUNIT(K))
              PBNDHI=PMLOG(LBNDM1,JUNIT(K))
              IF(PLOG.GE.PBNDLO)THEN

!                 INPUT PRESSURE IS ABOVE PRESSURE OF
!                 SECOND ATMOSPHERIC LEVEL:
                  WMOL(K)=AMOL(2,JUNIT(K),K)
                  WMOL(K)=WMOL(K)*(AMOL(1,JUNIT(K),K)/WMOL(K))          &
     &              **((PLOG-PBNDLO)/(PMLOG(1,JUNIT(K))-PBNDLO))
              ELSEIF(PLOG.LE.PBNDHI)THEN

!                 INPUT PRESSURE IS BELOW ABOVE PRESSURE OF
!                 SECOND-TO-LAST ATMOSPHERIC LEVEL:
                  WMOL(K)=AMOL(LBNDM1,JUNIT(K),K)
                  WMOL(K)=WMOL(K)*(AMOL(LBOUND,JUNIT(K),K)/WMOL(K))     &
     &              **((PLOG-PBNDHI)/(PMLOG(LBOUND,JUNIT(K))-PBNDHI))
              ELSEIF(JUNITP.EQ.JUNIT(K))THEN

!                 P AND MOL FROM THE SAME MODEL ATMOSPHERE:
                  Y0=AMOL(I0,JUNITP,K)
                  Y1=AMOL(I1,JUNITP,K)
                  Y2=AMOL(I2,JUNITP,K)
                  Y3=AMOL(I3,JUNITP,K)
                  WMOL(K)=FLG4PT(PLOG,P0,P1,P2,P3,Y0,Y1,Y2,Y3,LINIT)
                  LINIT=.FALSE.
              ELSE

!                 P AND MOL FROM DIFFERENT MODEL ATMOSPHERES:
                  DO 50 J2=2,LBNDM1
                      IF(PLOG.GE.PMLOG(J2,JUNIT(K)))GOTO60
   50             CONTINUE
   60             CONTINUE
                  J0=J2-2
                  J1=J2-1
                  J3=J2+1
                  P0X=PMLOG(J0,JUNIT(K))
                  P1X=PMLOG(J1,JUNIT(K))
                  P2X=PMLOG(J2,JUNIT(K))
                  P3X=PMLOG(J3,JUNIT(K))
                  Y0=AMOL(J0,JUNIT(K),K)
                  Y1=AMOL(J1,JUNIT(K),K)
                  Y2=AMOL(J2,JUNIT(K),K)
                  Y3=AMOL(J3,JUNIT(K),K)
                  LINIT=.TRUE.
                  WMOL(K)=FLG4PT(PLOG,P0X,P1X,P2X,P3X,Y0,Y1,Y2,Y3,LINIT)
              ENDIF
              JUNIT(K)=10

!             ADJUST CO2 MIXING RATIO BASED ON CARD1A INPUT.
              IF(K.EQ.2)WMOL(2)=CO2RAT*WMOL(2)
          ENDIF
   70 CONTINUE

!     MOLECULES CONSTANT (IN ALTITUDE) WITH MODEL ATMOSPHERE:
      DO 100 K=8,NMOL
          IF(JUNIT(K).LE.6)THEN
              KM7=K-7
              PBNDLO=PMLOG(2,JUNIT(K))
              PBNDHI=PMLOG(LBNDM1,JUNIT(K))
              IF(PLOG.GE.PBNDLO)THEN

!                 INPUT PRESSURE IS ABOVE PRESSURE OF
!                 SECOND ATMOSPHERIC LEVEL:
                  WMOL(K)=TRAC(2,KM7)
                  WMOL(K)=WMOL(K)*(TRAC(1,KM7)/WMOL(K))                 &
     &              **((PLOG-PBNDLO)/(PMLOG(1,JUNIT(K))-PBNDLO))
              ELSEIF(PLOG.LE.PBNDHI)THEN

!                 INPUT PRESSURE IS BELOW ABOVE PRESSURE OF
!                 SECOND-TO-LAST ATMOSPHERIC LEVEL:
                  WMOL(K)=TRAC(LBNDM1,KM7)
                  WMOL(K)=WMOL(K)*(TRAC(LBOUND,KM7)/WMOL(K))            &
     &              **((PLOG-PBNDHI)/(PMLOG(LBOUND,JUNIT(K))-PBNDHI))
              ELSEIF(JUNITP.EQ.JUNIT(K))THEN

!                 P AND MOL FROM THE SAME MODEL ATMOSPHERE:
                  Y0=TRAC(I0,KM7)
                  Y1=TRAC(I1,KM7)
                  Y2=TRAC(I2,KM7)
                  Y3=TRAC(I3,KM7)
                  WMOL(K)=FLG4PT(PLOG,P0,P1,P2,P3,Y0,Y1,Y2,Y3,LINIT)
                  LINIT=.FALSE.
              ELSE

!                 P AND MOL FROM DIFFERENT MODEL ATMOSPHERES:
                  DO 80 J2=2,LBNDM1
                      IF(PLOG.GE.PMLOG(J2,JUNIT(K)))GOTO90
   80             CONTINUE
   90             CONTINUE
                  J0=J2-2
                  J1=J2-1
                  J3=J2+1
                  P0X=PMLOG(J0,JUNIT(K))
                  P1X=PMLOG(J1,JUNIT(K))
                  P2X=PMLOG(J2,JUNIT(K))
                  P3X=PMLOG(J3,JUNIT(K))
                  Y0=TRAC(J0,KM7)
                  Y1=TRAC(J1,KM7)
                  Y2=TRAC(J2,KM7)
                  Y3=TRAC(J3,KM7)
                  LINIT=.TRUE.
                  WMOL(K)=FLG4PT(PLOG,P0X,P1X,P2X,P3X,Y0,Y1,Y2,Y3,LINIT)
              ENDIF
              JUNIT(K)=10
          ENDIF
  100 CONTINUE
      WMOL(12)=1000.*WMOL(12)

!     MOLECULAR DENSITIES FOR CFC SPECIES:
      IF(JUNITX.LE.6)THEN
          DO 130 K=1,NMOLX
              PBNDLO=PMLOG(2,JUNITX)
              PBNDHI=PMLOG(LBNDM1,JUNITX)
              IF(PLOG.GE.PBNDLO)THEN

!                 INPUT PRESSURE IS ABOVE PRESSURE OF
!                 SECOND ATMOSPHERIC LEVEL:
                  WMOLX(K)=AMOLX(2,K)
                  WMOLX(K)=WMOLX(K)*(AMOLX(1,K)/WMOLX(K))               &
     &              **((PLOG-PBNDLO)/(PMLOG(1,JUNITX)-PBNDLO))
              ELSEIF(PLOG.LE.PBNDHI)THEN

!                 INPUT PRESSURE IS BELOW ABOVE PRESSURE OF
!                 SECOND-TO-LAST ATMOSPHERIC LEVEL:
                  WMOLX(K)=AMOLX(LBNDM1,K)
                  WMOLX(K)=WMOLX(K)*(AMOLX(LBOUND,K)/WMOLX(K))          &
     &              **((PLOG-PBNDHI)/(PMLOG(LBOUND,JUNITX)-PBNDHI))
              ELSE
                  IF(JUNITP.NE.JUNITX .AND. K.EQ.1)THEN

!                     P AND MOL FROM THE SAME MODEL ATMOSPHERE:
                      DO 110 I2=2,LBNDM1
                          IF(PLOG.GE.PMLOG(I2,JUNITX))GOTO120
  110                 CONTINUE
  120                 CONTINUE
                      I0=I2-2
                      I1=I2-1
                      I3=I2+1
                      P0=PMLOG(I0,JUNITX)
                      P1=PMLOG(I1,JUNITX)
                      P2=PMLOG(I2,JUNITX)
                      P3=PMLOG(I3,JUNITX)
                      LINIT=.TRUE.
                  ENDIF
                  Y0=AMOLX(I0,K)
                  Y1=AMOLX(I1,K)
                  Y2=AMOLX(I2,K)
                  Y3=AMOLX(I3,K)
                  WMOLX(K)=FLG4PT(PLOG,P0,P1,P2,P3,Y0,Y1,Y2,Y3,LINIT)
                  LINIT=.FALSE.
              ENDIF
  130     CONTINUE
      ENDIF
      RETURN
      END
      REAL FUNCTION HYSTAT(ZREF,TREF,PREF,T,P)

!     INTEGRATES HYDROSTATIC EQUATION ASSUMING TEMPERATURE
!     VARIES LINEARLY WITH LOG(PRESSURE):

!     ARGUMENTS:
!       ZREF     REFERENCE ALTITUDE [KM].
!       TREF     REFERENCE TEMPERATURE [K].
!       PREF     REFERENCE PRESSURE [MB].
!       T        TEMPERATURE AT UNKNOWN ALTITUDE [K].
!       P        PRESSURE AT UNKNOWN ALTITUDE [MB].
      REAL ZREF,TREF,PREF,T,P

!     PARAMETERS:
!       GRAV     GRAVITATIONAL CONSTANT OF EARTH [KM3/SEC2].
!       AIRMWT   MOLECULAR WEIGHT OF AIR [GM/MOL].
!       GASCON   GAS CONSTANT [GM KM2 / SEC2 MOL K].
!       HSTAT    CONSTANT USED IN HYDROSTATIC EQUATION [KM K].
      REAL GRAV,AIRMWT,GASCON,HSTAT
      PARAMETER(GRAV=398600.4418,AIRMWT=28.964,GASCON=.00831441,        &
     &  HSTAT=2*AIRMWT*GRAV/GASCON)

!     /CARD3/
!       H1      OBSERVER (SENSOR) ALTITUDE [KM].
!       H2      FINAL (TARGET) ALTITUDE [KM].
!       ANGLE   ZENITH ANGLE FROM H1 TO H2 [DEG].
!       RANGE   DISTANCE FROM H1 TO H2 [KM].
!       BETA    EARTH CENTER ANGLE BETWEEN H1 AND H2 [DEG].
!       REE     RADIUS OF THE EARTH [KM].
!       LENN    PATH LENGTH SWITCH (0=SHORT, 1=LONG).
      INTEGER LENN
      REAL H1,H2,ANGLE,RANGE,BETA,REE
      COMMON/CARD3/H1,H2,ANGLE,RANGE,BETA,REE,LENN

!     LOCAL VARIABLES:
!       RADIUS   EARTH CENTER DISTANCE [KM].
!       PRATLN   LOGARITHM OF THE PRESSURE RATIO.
      REAL RADIUS,PRATLN

!     INTEGRATE:
      RADIUS=REE+ZREF
      PRATLN=LOG(PREF/P)
      HYSTAT=ZREF+PRATLN*RADIUS/(HSTAT/(RADIUS*(TREF+T))-PRATLN)
      RETURN
      END
