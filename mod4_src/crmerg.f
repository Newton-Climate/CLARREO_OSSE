      SUBROUTINE CRMERG

!     THIS ROUTINE MERGES TOGETHER CLOUD/RAIN AND OLD ATMOSPHERIC
!     PROFILES.  WATER DROPLET DENSITIES [GM/M3], ICE PARTICLE
!     DENSITIES [GM/M3] AND RAIN RATES [MM/HR] ARE STORED IN
!     DENSTY(66,.), DENSTY(67,.) AND DENSTY(3,.), RESPECTIVELY.
!     LAYER BOUNDARIES ARE MERGED TOGETHER IF THEY DIFFER BY LESS
!     THAN "ZTOL" KM (HALF A METER).

!     PARAMETERS:
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:
      INCLUDE 'BASE.h'
      REAL P,T,WH,WCO2,WO,WN2O,WCO,WCH4,WO2
      COMMON/MDATA/P(LAYDIM),T(LAYDIM),WH(LAYDIM),WCO2(LAYDIM),         &
     &  WO(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),WCH4(LAYDIM),WO2(LAYDIM)
      REAL WNO,WSO2,WNO2,WNH3,WHNO3
      COMMON/MDATA1/WNO(LAYDIM),WSO2(LAYDIM),WNO2(LAYDIM),              &
     &  WNH3(LAYDIM),WHNO3(LAYDIM)
      REAL WMOLXT
      COMMON/MDATAX/WMOLXT(NMOLX,LAYDIM)
      INCLUDE 'IFIL.h'

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

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT
      INTEGER NBND
      REAL ZCLDRN(NZCLD),DRPWAT(NZCLD),PRTICE(NZCLD),RNPROF(NZCLD)
      COMMON/CLDRN/NBND,ZCLDRN,DRPWAT,PRTICE,RNPROF
      INTEGER NCRALT,NCRSPC
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      &
     &  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     DECLARE LOCAL VARIABLES
      INTEGER NEWBND,MBND,IBND,I,MLOLD
      LOGICAL LSET,WARNWD,WARNIP
      REAL FAC,TRATIO,CLDHUM

!     DECLARE FUNCTION NAMES
      REAL EXPINT

!     LIST DATA
!       ZTOL     TOLERANCE FOR MERGING TOGETHER LAYERS [KM]
!       TFREEZ   TEMPERATURE BELOW WHICH A WARNING IS GENERATED IF
!                LIQUID WATER DROPLET DENSITY IS POSITIVE [K]
!       TMELT    TEMPERATURE ABOVE WHICH A WARNING IS GENERATED IF
!                ICE PARTICLE DENSITY IS POSITIVE [K]
      REAL ZTOL,TFREEZ,TMELT
      DATA ZTOL/.0005/,TFREEZ/260./,TMELT/278./

!     CHECK THAT CLOUD/RAIN ALTITUDES ARE BOUNDED BY ZM ALTITUDES
      IF(ZCLDRN(1).LT.ZM(1) .OR. ZCLDRN(NBND).GT.ZM(ML))THEN
          WRITE(IPR,'(/A,2(F10.5,A),/14X,A,2(F10.5,A))')                &
     &      ' FATAL ERROR:  CLOUD/RAIN MODEL BOUNDING ALTITUDES (',     &
     &      ZCLDRN(1),' AND',ZCLDRN(NBND),' KM) ARE',                   &
     &      ' NOT BRACKETED BY ATMOSPHERE BOUNDING ALTITUDES (',        &
     &      ZM(1),' AND', ZM(ML),' KM).'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP
      ENDIF

!     CHECK THAT RAIN RATE AND CLOUD WATER DROPLET AND ICE PARTICLE
!     DENSITIES ARE INDEED ZERO AT THE TOP OF THE PROFILE.
      IF(DRPWAT(NBND).NE.0. .OR. PRTICE(NBND).NE.0.                     &
     &                      .OR. RNPROF(NBND).NE.0.)THEN
          WRITE(IPR,'(/2A,/14X,A)')' FATAL ERROR:  RAIN RATE AND',      &
     &      ' CLOUD WATER DROPLET AND ICE PARTICLE DENSITIES AND',      &
     &      ' ARE NOT ZERO AT THE TOP OF THE CLOUD/RAIN PROFILE.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP
      ENDIF

!     INITIALIZE WARNING MESSAGE VARIABLES
!     WARNWD   A LOGICAL VARIABLE THAT IS SET TO FALSE AFTER THE USER
!              HAS BEEN WARNED THAT WATER DROPLETS EXIST BELOW TFREEZ.
!     WARNIP   A LOGICAL VARIABLE THAT IS SET TO FALSE AFTER THE USER
!              HAS BEEN WARNED THAT ICE PARTICLES EXIST ABOVE TMELT.
      WARNWD=.TRUE.
      WARNIP=.TRUE.

!     DETERMINE THE RELATIVE HUMIDITY [%] WITHIN THE CLOUD.
!     THE RELATIVE HUMIDITY MUST BE POSITIVE (A DRY ATMOSPHERE
!     IS NOT ALLOWED), AND SUPER-SATURATION IS LIMITED TO 5%.
      IF(CHUMID.GT.0.)THEN
          CLDHUM=CHUMID
      ELSE
          CLDHUM=100.
      ENDIF

!     TEMPORARILY TRANSLATE BOUNDARY LAYER DATA TO MAXIMUM
!     LAYER INDICES TO MAKE SPACE FOR ADDITIONAL LAYERS.
      NEWBND=LAYDIM+1
      DO 20 MBND=ML,2,-1
          NEWBND=NEWBND-1
          ZM(NEWBND)=ZM(MBND)
          P(NEWBND)=P(MBND)
          T(NEWBND)=T(MBND)
          RELHUM(NEWBND)=RELHUM(MBND)
          WH(NEWBND)=WH(MBND)
          LRHSET(NEWBND)=LRHSET(MBND)
          WCO2(NEWBND)=WCO2(MBND)
          WO(NEWBND)=WO(MBND)
          WN2O(NEWBND)=WN2O(MBND)
          WCO(NEWBND)=WCO(MBND)
          WCH4(NEWBND)=WCH4(MBND)
          WO2(NEWBND)=WO2(MBND)
          WHNO3(NEWBND)=WHNO3(MBND)
          WNO(NEWBND)=WNO(MBND)
          WSO2(NEWBND)=WSO2(MBND)
          WNO2(NEWBND)=WNO2(MBND)
          WNH3(NEWBND)=WNH3(MBND)
          DO 10 I=1,NMOLX
              WMOLXT(I,NEWBND)=WMOLXT(I,MBND)
   10     CONTINUE
          DENSTY(7,NEWBND)=DENSTY(7,MBND)
          DENSTY(12,NEWBND)=DENSTY(12,MBND)
          DENSTY(13,NEWBND)=DENSTY(13,MBND)
          DENSTY(14,NEWBND)=DENSTY(14,MBND)
          DENSTY(15,NEWBND)=DENSTY(15,MBND)
          DENSTY(16,NEWBND)=DENSTY(16,MBND)
   20 CONTINUE

!     INITIALIZE RAIN RATE AND WATER DROPLET AND
!     ICE PARTICLE DENSITIES AT THE GROUND.
      ML=1
      DENSTY(66,1)=0.
      DENSTY(67,1)=0.
      DENSTY(3,1)=0.

!     LOOP OVER OLD ATMOSPHERIC LAYERS
      IBND=1
      DO 60 MBND=NEWBND,LAYDIM
          LSET=.FALSE.
   30     CONTINUE
          IF(ZCLDRN(IBND).LE.ZM(MBND))THEN

!             CLOUD LAYER BOUNDARY WITHIN OLD ATMOSPHERIC LAYER
              IF(ZCLDRN(IBND).LT.ZM(ML)+ZTOL)THEN

!                 MERGE INTO LOWER LAYER BOUNDARY:
                  DENSTY(66,ML)=DRPWAT(IBND)
                  DENSTY(67,ML)=PRTICE(IBND)
                  DENSTY(3,ML)=RNPROF(IBND)
                  IF(CLDHUM.LE.105. .AND. (DRPWAT(IBND).GT.0. .OR.      &
     &              PRTICE(IBND).GT.0. .OR. RNPROF(IBND).GT.0.))THEN
                      RELHUM(ML)=CLDHUM
                      TRATIO=273.15/T(ML)
                      WH(ML)=.01*CLDHUM*TRATIO*                         &
     &                  EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
                      LRHSET(ML)=.TRUE.
                      IF(WARNWD .AND. DRPWAT(IBND).GT.0.                &
     &                          .AND. T(ML).LT.TFREEZ)THEN
                          WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,2A))')   &
     &                      ' WARNING:  AT ALTITUDE',ZM(ML),            &
     &                      ' KM, THE LIQUID WATER DROPLET DENSITY',    &
     &                      ' IS POSITIVE [',DRPWAT(IBND),' GM/M3]',    &
     &                      ' EVEN THOUGH',' THE TEMPERATURE IS',       &
     &                      T(ML)-273.15,' DEGREES',' CELSIUS.'
                          WARNWD=.FALSE.
                      ELSEIF(WARNIP .AND. PRTICE(IBND).GT.0.            &
     &                              .AND. T(ML).GT.TMELT)THEN
                          WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,2A))')   &
     &                      ' WARNING:  AT ALTITUDE',ZM(ML),            &
     &                      ' KM, THE ICE PARITCLE DENSITY',            &
     &                      ' IS POSITIVE [',PRTICE(IBND),' GM/M3]',    &
     &                      ' EVEN THOUGH',' THE TEMPERATURE IS',       &
     &                      T(ML)-273.15,' DEGREES',' CELSIUS.'
                          WARNIP=.FALSE.
                      ENDIF
                  ENDIF
              ELSEIF(ZCLDRN(IBND).GT.ZM(MBND)-ZTOL)THEN

!                 MERGE INTO UPPER LAYER BOUNDARY
                  LSET=.TRUE.
                  DENSTY(66,MBND)=DRPWAT(IBND)
                  DENSTY(67,MBND)=PRTICE(IBND)
                  DENSTY(3,MBND)=RNPROF(IBND)
                  IF(CLDHUM.LE.105. .AND. (DRPWAT(IBND).GT.0. .OR.      &
     &              PRTICE(IBND).GT.0. .OR. RNPROF(IBND).GT.0.))THEN
                      RELHUM(MBND)=CLDHUM
                      TRATIO=273.15/T(MBND)
                      WH(MBND)=.01*CLDHUM*TRATIO*                       &
     &                  EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
                      LRHSET(ML)=.TRUE.
                  ENDIF
              ELSEIF(ML+1.LT.MBND)THEN

!                 ADD NEW ZM LAYER BOUNDARY
                  MLOLD=ML
                  ML=ML+1
                  ZM(ML)=ZCLDRN(IBND)
                  DENSTY(66,ML)=DRPWAT(IBND)
                  DENSTY(67,ML)=PRTICE(IBND)
                  DENSTY(3,ML)=RNPROF(IBND)
                  FAC=(ZM(ML)-ZM(MLOLD))/(ZM(MBND)-ZM(MLOLD))
                  P(ML)=EXPINT(P(MLOLD),P(MBND),FAC)
                  T(ML)=EXPINT(T(MLOLD),T(MBND),FAC)
                  TRATIO=273.15/T(ML)
                  IF(CLDHUM.LE.105. .AND. (DRPWAT(IBND).GT.0. .OR.      &
     &              PRTICE(IBND).GT.0. .OR. RNPROF(IBND).GT.0.))THEN
                      RELHUM(ML)=CLDHUM
                      WH(ML)=.01*CLDHUM*TRATIO*                         &
     &                  EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
                      LRHSET(ML)=.TRUE.
                      IF(WARNWD .AND. DRPWAT(IBND).GT.0.                &
     &                  .AND. T(ML).LT.TFREEZ)THEN
                          WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,2A))')   &
     &                      ' WARNING:  AT ALTITUDE',ZM(ML),            &
     &                      ' KM, THE LIQUID WATER DROPLET DENSITY',    &
     &                      ' IS POSITIVE [',DRPWAT(IBND),' GM/M3]',    &
     &                      ' EVEN THOUGH',' THE TEMPERATURE IS',       &
     &                      T(ML)-273.15,' DEGREES',' CELSIUS.'
                          WARNWD=.FALSE.
                      ELSEIF(WARNIP .AND. PRTICE(IBND).GT.0.            &
     &                              .AND. T(ML).GT.TMELT)THEN
                          WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,2A))')   &
     &                      ' WARNING:  AT ALTITUDE',ZM(ML),            &
     &                      ' KM, THE ICE PARITCLE DENSITY',            &
     &                      ' IS POSITIVE [',PRTICE(IBND),' GM/M3]',    &
     &                      ' EVEN THOUGH',' THE TEMPERATURE IS',       &
     &                      T(ML)-273.15,' DEGREES',' CELSIUS.'
                          WARNIP=.FALSE.
                      ENDIF
                  ELSE
                      RELHUM(ML)=EXPINT(RELHUM(MLOLD),RELHUM(MBND),FAC)
                      WH(ML)=.01*RELHUM(ML)*TRATIO*                     &
     &                  EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
                      LRHSET(ML)=.FALSE.
                  ENDIF
                  WCO2(ML)=EXPINT(WCO2(MLOLD),WCO2(MBND),FAC)
                  WO(ML)=EXPINT(WO(MLOLD),WO(MBND),FAC)
                  WN2O(ML)=EXPINT(WN2O(MLOLD),WN2O(MBND),FAC)
                  WCO(ML)=EXPINT(WCO(MLOLD),WCO(MBND),FAC)
                  WCH4(ML)=EXPINT(WCH4(MLOLD),WCH4(MBND),FAC)
                  WO2(ML)=EXPINT(WO2(MLOLD),WO2(MBND),FAC)
                  WHNO3(ML)=EXPINT(WHNO3(MLOLD),WHNO3(MBND),FAC)
                  WNO(ML)=EXPINT(WNO(MLOLD),WNO(MBND),FAC)
                  WSO2(ML)=EXPINT(WSO2(MLOLD),WSO2(MBND),FAC)
                  WNO2(ML)=EXPINT(WNO2(MLOLD),WNO2(MBND),FAC)
                  WNH3(ML)=EXPINT(WNH3(MLOLD),WNH3(MBND),FAC)
                  DO 40 I=1,NMOLX
                      WMOLXT(I,ML)=                                     &
     &                  EXPINT(WMOLXT(I,MLOLD),WMOLXT(I,MBND),FAC)
   40             CONTINUE
                  DENSTY(7,ML)=                                         &
     &              EXPINT(DENSTY(7,MLOLD),DENSTY(7,MBND),FAC)
                  DENSTY(12,ML)=                                        &
     &              EXPINT(DENSTY(12,MLOLD),DENSTY(12,MBND),FAC)
                  DENSTY(13,ML)=                                        &
     &              EXPINT(DENSTY(13,MLOLD),DENSTY(13,MBND),FAC)
                  DENSTY(14,ML)=                                        &
     &              EXPINT(DENSTY(14,MLOLD),DENSTY(14,MBND),FAC)
                  DENSTY(15,ML)=                                        &
     &              EXPINT(DENSTY(15,MLOLD),DENSTY(15,MBND),FAC)
                  DENSTY(16,ML)=                                        &
     &              EXPINT(DENSTY(16,MLOLD),DENSTY(16,MBND),FAC)
              ELSE

!                 NO MORE SPACE IN ARRAYS FOR AN ADDITIONAL LAYER
                  WRITE(IPR,'(/3A,/14X,A,2(A,I4))')' FATAL ERROR: ',    &
     &              ' FILE "PARAMS.h" PARAMETER "LAYDIM" MUST',         &
     &              ' BE INCREASED.',' IT SUFFICES TO INCREASE',        &
     &              ' LAYDIM FROM',LAYDIM,' TO',LAYDIM+NBND-IBND+1
                  IF(LJMASS)CALL WRTBUF(FATAL)
                  STOP
              ENDIF

!             EXIT LOOP IF ALL CLOUD BOUNDARIES HAVE BEEN INTEGRATED.
              IF(IBND.GE.NBND)GOTO70

!             INCREMENT CLOUD BOUNDARY INDEX AND START AGAIN
              IBND=IBND+1
              GOTO30
          ENDIF

!         LINEARLY INTERPOLATE WATER PARTICLE DENSITIES AND RAIN RATE.
          IF(.NOT.LSET)THEN
              FAC=(ZM(MBND)-ZM(ML))/(ZCLDRN(IBND)-ZM(ML))
              DENSTY(66,MBND)=DENSTY(66,ML)                             &
     &          +FAC*(DRPWAT(IBND)-DENSTY(66,ML))
              DENSTY(67,MBND)=DENSTY(67,ML)                             &
     &          +FAC*(PRTICE(IBND)-DENSTY(67,ML))
              DENSTY(3,MBND)=DENSTY(3,ML)                               &
     &          +FAC*(RNPROF(IBND)-DENSTY(3,ML))
              IF(CLDHUM.LE.105. .AND. (DENSTY(66,MBND).GT.0. .OR.       &
     &          DENSTY(67,MBND).GT.0. .OR. DENSTY( 3,MBND).GT.0.))THEN
                  RELHUM(MBND)=CLDHUM
                  TRATIO=273.15/T(MBND)
                  WH(MBND)=.01*CLDHUM*TRATIO*                           &
     &              EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
                  LRHSET(ML)=.TRUE.
              ENDIF
          ENDIF

!         TRANSLATE LAYER BOUNDARY DATA TO NEW BOUNDARY INDEX.
          ML=ML+1
          ZM(ML)=ZM(MBND)
          DENSTY(66,ML)=DENSTY(66,MBND)
          DENSTY(67,ML)=DENSTY(67,MBND)
          DENSTY(3,ML)=DENSTY(3,MBND)
          P(ML)=P(MBND)
          T(ML)=T(MBND)
          IF(WARNWD .AND. DENSTY(66,ML).GT.0. .AND. T(ML).LT.TFREEZ)THEN
              WRITE(IPR,'(//A,F9.4,A,3(/10X,A,F9.4,A))')                &
     &          ' WARNING:  AT ALTITUDE',ZM(ML),                        &
     &          ' KM, THE LIQUID WATER DROPLET DENSITY',                &
     &          ' IS POSITIVE [',DENSTY(66,ML),' GM/M3] EVEN THOUGH',   &
     &          ' THE TEMPERATURE IS',T(ML)-273.15,' DEGREES CELSIUS.'
              WARNWD=.FALSE.
          ENDIF
          IF(WARNIP .AND. DENSTY(67,ML).GT.0. .AND. T(ML).GT.TMELT)THEN
              WRITE(IPR,'(//2A,F9.4,A,3(/10X,A,F9.4,A))')' WARNING: ',  &
     &          ' AT ALTITUDE',ZM(ML),' KM, THE ICE PARITCLE DENSITY',  &
     &          ' IS POSITIVE [',DENSTY(67,ML),' GM/M3] EVEN THOUGH',   &
     &          ' THE TEMPERATURE IS',T(ML)-273.15,' DEGREES CELSIUS.'
              WARNIP=.FALSE.
          ENDIF
          RELHUM(ML)=RELHUM(MBND)
          WH(ML)=WH(MBND)
          LRHSET(ML)=LRHSET(MBND)
          WCO2(ML)=WCO2(MBND)
          WO(ML)=WO(MBND)
          WN2O(ML)=WN2O(MBND)
          WCO(ML)=WCO(MBND)
          WCH4(ML)=WCH4(MBND)
          WO2(ML)=WO2(MBND)
          WHNO3(ML)=WHNO3(MBND)
          WNO(ML)=WNO(MBND)
          WSO2(ML)=WSO2(MBND)
          WNO2(ML)=WNO2(MBND)
          WNH3(ML)=WNH3(MBND)
          DO 50 I=1,NMOLX
              WMOLXT(I,ML)=WMOLXT(I,MBND)
   50     CONTINUE
          DENSTY(7,ML)=DENSTY(7,MBND)
          DENSTY(12,ML)=DENSTY(12,MBND)
          DENSTY(13,ML)=DENSTY(13,MBND)
          DENSTY(14,ML)=DENSTY(14,MBND)
          DENSTY(15,ML)=DENSTY(15,MBND)
          DENSTY(16,ML)=DENSTY(16,MBND)
   60 CONTINUE
   70 CONTINUE

!     TRANSLATE REMAINING LAYER BOUNDARY DATA TO NEW BOUNDARY INDEX
!     SETTING CLOUD PARTICLE DENSITIES AND RAIN RATES TO ZERO.
      NEWBND=MBND
      DO 90 MBND=NEWBND,LAYDIM
          ML=ML+1
          ZM(ML)=ZM(MBND)
          DENSTY(66,ML)=0.
          DENSTY(67,ML)=0.
          DENSTY(3,ML)=0.
          P(ML)=P(MBND)
          T(ML)=T(MBND)
          RELHUM(ML)=RELHUM(MBND)
          WH(ML)=WH(MBND)
          LRHSET(ML)=LRHSET(MBND)
          WCO2(ML)=WCO2(MBND)
          WO(ML)=WO(MBND)
          WN2O(ML)=WN2O(MBND)
          WCO(ML)=WCO(MBND)
          WCH4(ML)=WCH4(MBND)
          WO2(ML)=WO2(MBND)
          WHNO3(ML)=WHNO3(MBND)
          WNO(ML)=WNO(MBND)
          WSO2(ML)=WSO2(MBND)
          WNO2(ML)=WNO2(MBND)
          WNH3(ML)=WNH3(MBND)
          DO 80 I=1,NMOLX
              WMOLXT(I,ML)=WMOLXT(I,MBND)
   80     CONTINUE
          DENSTY(7,ML)=DENSTY(7,MBND)
          DENSTY(12,ML)=DENSTY(12,MBND)
          DENSTY(13,ML)=DENSTY(13,MBND)
          DENSTY(14,ML)=DENSTY(14,MBND)
          DENSTY(15,ML)=DENSTY(15,MBND)
          DENSTY(16,ML)=DENSTY(16,MBND)
   90 CONTINUE
      RETURN
      END
