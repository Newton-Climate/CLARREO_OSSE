      SUBROUTINE NOVMRG(ALTNOV,RHNOV,PNOV,TNOV,DENNOV,NLNOV)

!     THIS ROUTINE MERGES TOGETHER NOVAM LAYERS WITH CLOUD/RAIN
!     AND OLD ATMOSPHERIC PROFILES.
!     LAYER BOUNDARIES ARE MERGED TOGETHER IF THEY DIFFER BY LESS
!     THAN "ZTOL" KM (HALF A METER).

!     THE NOVAM DENSITIES ARE DENSTY(68,*) THROUGH DENSITY(77,*)
!     THESE OVERLAP

!     LIST PARAMETERS
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:
      INCLUDE 'BASE.h'
      REAL P,T,WH,WCO2,WO,WN2O,WCO,WCH4,WO2
      COMMON/MDATA/P(LAYDIM),T(LAYDIM),WH(LAYDIM),WCO2(LAYDIM),         &
     &     WO(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),WCH4(LAYDIM),WO2(LAYDIM)
      REAL WNO,WSO2,WNO2,WNH3,WHNO3
      COMMON/MDATA1/WNO(LAYDIM),WSO2(LAYDIM),WNO2(LAYDIM),              &
     &     WNH3(LAYDIM),WHNO3(LAYDIM)
      REAL WMOLXT
      COMMON/MDATAX/WMOLXT(NMOLX,LAYDIM)
      REAL GNDALT
      COMMON/GRAUND/GNDALT

      INTEGER NNOV,NLNOV
      REAL ALTNOV(NLNOV),RHNOV(NLNOV),TNOV(NLNOV),PNOV(NLNOV),          &
     &     DENNOV(MNOV,MLNOV)
      INCLUDE 'IFIL.h'
      REAL ZM,PM,TM,RFNDX,DENSTY
      LOGICAL LRHSET
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    &
     &     RFNDX(LAYDIM),DENSTY(MEXT,LAYDIM),LRHSET(LAYDIM)

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     DECLARE LOCAL VARIABLES
      INTEGER NEWBND,IL,IMERGE,I,J,II,MLOLD,INOV,JNOV
      LOGICAL LSET
      REAL FAC,TRATIO

!     DECLARE FUNCTION NAMES
      REAL EXPINT
      REAL ALTOFF

!     LIST DATA
!     ZTOL     TOLERANCE FOR MERGING TOGETHER LAYERS [KM]
      REAL ZTOL
      DATA ZTOL/.0005/
      IF(GNDALT.NE.0)                                                   &
     &   WRITE(IPR,'(/2A)')' WARNING:  GROUND ALTITUDE',                &
     &     ' NOT AT SEA LEVEL, BUT NOVAM AEROSOLS ARE BEING USED.'
      ALTOFF=ALTNOV(1)-GNDALT
      ALTNOV(1)=GNDALT
      DO 5 IMERGE=2,NLNOV
          ALTNOV(IMERGE)=ALTNOV(IMERGE)-ALTOFF
    5 CONTINUE

!     CHECK THAT ALTNOV, NOVAM ALTITUDES, ARE BOUNDED BY ZM
      IF(ALTNOV(1).LT.ZM(1) .OR. ALTNOV(NLNOV).GT.ZM(ML))THEN
         WRITE(IPR,'(/A,2(F10.5,A),/14X,A,2(F10.5,A))')                 &
     &        ' FATAL ERROR:  NOVAM BOUNDING ALTITUDES (',              &
     &        ALTNOV(1),' AND',ALTNOV(NLNOV),' KM) ARE',                &
     &        ' NOT BRACKETED BY ATMOSPHERE BOUNDING ALTITUDES (',      &
     &        ZM(1),' AND', ZM(ML),' KM).'
         IF(LJMASS)CALL WRTBUF(FATAL)
         STOP
      ENDIF

!     COMPUTE NUMBER OF NOVAM AEROSOLS; IT'S IN THE COMMON BLOCK COMNOV
      NNOV=(NLNOV-2)/2

!     TEMPORARILY TRANSLATE BOUNDARY LAYER DATA TO MAXIMUM
!     LAYER INDICES TO MAKE SPACE FOR ADDITIONAL LAYERS.
!     KEEP Z(1) WHERE IT IS BECAUSE IT WILL ALWAYS BE THERE.
!     THAT IS, IT IS THE LOWEST LAYER.
      NEWBND=LAYDIM+1
      DO 20 IL=ML,2,-1
         NEWBND=NEWBND-1
         ZM(NEWBND)=ZM(IL)
         P(NEWBND)=P(IL)
         T(NEWBND)=T(IL)
         RELHUM(NEWBND)=RELHUM(IL)
         WH(NEWBND)=WH(IL)
         WCO2(NEWBND)=WCO2(IL)
         WO(NEWBND)=WO(IL)
         WN2O(NEWBND)=WN2O(IL)
         WCO(NEWBND)=WCO(IL)
         WCH4(NEWBND)=WCH4(IL)
         WO2(NEWBND)=WO2(IL)
         WHNO3(NEWBND)=WHNO3(IL)
         WNO(NEWBND)=WNO(IL)
         WSO2(NEWBND)=WSO2(IL)
         WNO2(NEWBND)=WNO2(IL)
         WNH3(NEWBND)=WNH3(IL)
         DO 10 I=1,NMOLX
            WMOLXT(I,NEWBND)=WMOLXT(I,IL)
 10      CONTINUE

!        INDICES 7, 12-15 = AEROSOLS, 16 = CIRRUS CLOUD
         DENSTY(7,NEWBND)=DENSTY(7,IL)
         DENSTY(12,NEWBND)=DENSTY(12,IL)
         DENSTY(13,NEWBND)=DENSTY(13,IL)
         DENSTY(14,NEWBND)=DENSTY(14,IL)
         DENSTY(15,NEWBND)=DENSTY(15,IL)
         DENSTY(16,NEWBND)=DENSTY(16,IL)

!        66 = WATER DROPLETS (G/M**3), 67 = ICE PARTICLES (G/M**3)
!        3 = RAIN RATE (MM/HR)
         DENSTY(66,NEWBND)=DENSTY(66,IL)
         DENSTY(67,NEWBND)=DENSTY(67,IL)
         DENSTY(3,NEWBND)=DENSTY(3,IL)

!        INITIALIZE NOVAM DENSTY (ACTUALLY EXTINCTION COEFF)
         DENSTY(68,NEWBND)=0.0
         DENSTY(69,NEWBND)=0.0
         DENSTY(70,NEWBND)=0.0
         DENSTY(71,NEWBND)=0.0
         DENSTY(72,NEWBND)=0.0
         DENSTY(73,NEWBND)=0.0
         DENSTY(74,NEWBND)=0.0
         DENSTY(75,NEWBND)=0.0
         DENSTY(76,NEWBND)=0.0
         DENSTY(77,NEWBND)=0.0
 20   CONTINUE

!     ONE MORE INITIALIZATION
      DENSTY(68,1)=0.0
      DENSTY(69,1)=0.0
      DENSTY(70,1)=0.0
      DENSTY(71,1)=0.0
      DENSTY(72,1)=0.0
      DENSTY(73,1)=0.0
      DENSTY(74,1)=0.0
      DENSTY(75,1)=0.0
      DENSTY(76,1)=0.0
      DENSTY(77,1)=0.0

!     LOOP OVER OLD LAYERS (NOW AT 1, NEWBND, NEWBND+1, ... LAYDIM)
!     TO SEE WHERE NOVAM LAYERS FIT IN
      IMERGE=1
      ML = 1
      DO 60 IL=NEWBND,LAYDIM
         LSET=.FALSE.
 30      CONTINUE
         IF(ALTNOV(IMERGE).LE.ZM(IL))THEN

!           NOVAM LAYER BOUNDARY WITHIN OLD ATMOSPHERIC LAYER
            IF(ALTNOV(IMERGE).LT.ZM(ML)+ZTOL)THEN

!              MERGE INTO LOWER LAYER BOUNDARY
               RELHUM(ML)=RHNOV(IMERGE)
!              DENSTY(68,ML)=EXTNOV(4,IMERGE)
               DO 160 INOV=1,NNOV
                  JNOV=INOV-1
                  DENSTY(68+JNOV,ML)=DENNOV(INOV,IMERGE)
 160           CONTINUE
               DENSTY(7,ML)=0.0
               T(ML)=TNOV(IMERGE)
               P(ML)=PNOV(IMERGE)
               TRATIO=273.15/T(ML)
!              THE FORMULA FOR WH IS FROM AERNSM (WH IS IN G/M**3)
               WH(ML)=.01*RELHUM(ML)*TRATIO*                            &
     &              EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
            ELSEIF(ALTNOV(IMERGE).GT.ZM(IL)-ZTOL)THEN

!              MERGE INTO UPPER LAYER BOUNDARY
               LSET=.TRUE.
               RELHUM(IL)=RHNOV(IMERGE)
!              DENSTY(68,IL)=EXTNOV(4,IMERGE)
               DO 170 INOV=1,NNOV
                  JNOV=INOV-1
                  DENSTY(68+JNOV,IL)=DENNOV(INOV,IMERGE)
 170           CONTINUE
               DENSTY(7,IL)=0.0
               T(IL)=TNOV(IMERGE)
               P(IL)=PNOV(IMERGE)
               TRATIO=273.15/T(IL)
               WH(IL)=.01*RELHUM(IL)*TRATIO*                            &
     &              EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
            ELSEIF(ML+1.LT.IL)THEN

!              ADD NEW ZM LAYER BOUNDARY
               MLOLD=ML
               ML=ML+1
               RELHUM(ML)=RHNOV(IMERGE)
               ZM(ML)=ALTNOV(IMERGE)
!              DENSTY(68,ML)=EXTNOV(4,IMERGE)
               DO 180 INOV=1,NNOV
                  JNOV=INOV-1
                  DENSTY(68+JNOV,ML)=DENNOV(INOV,IMERGE)
 180           CONTINUE
               DENSTY(7,ML)=0.0
               T(ML)=TNOV(IMERGE)
               P(ML)=PNOV(IMERGE)

               FAC=(ZM(ML)-ZM(MLOLD))/(ZM(IL)-ZM(MLOLD))
               TRATIO=273.15/T(ML)
               WH(ML)=.01*RELHUM(ML)*TRATIO*                            &
     &              EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)

               WCO2(ML)=EXPINT(WCO2(MLOLD),WCO2(IL),FAC)
               WO(ML)=EXPINT(WO(MLOLD),WO(IL),FAC)
               WN2O(ML)=EXPINT(WN2O(MLOLD),WN2O(IL),FAC)
               WCO(ML)=EXPINT(WCO(MLOLD),WCO(IL),FAC)
               WCH4(ML)=EXPINT(WCH4(MLOLD),WCH4(IL),FAC)
               WO2(ML)=EXPINT(WO2(MLOLD),WO2(IL),FAC)
               WHNO3(ML)=EXPINT(WHNO3(MLOLD),WHNO3(IL),FAC)
               WNO(ML)=EXPINT(WNO(MLOLD),WNO(IL),FAC)
               WSO2(ML)=EXPINT(WSO2(MLOLD),WSO2(IL),FAC)
               WNO2(ML)=EXPINT(WNO2(MLOLD),WNO2(IL),FAC)
               WNH3(ML)=EXPINT(WNH3(MLOLD),WNH3(IL),FAC)
               DO 40 I=1,NMOLX
                  WMOLXT(I,ML)=                                         &
     &                 EXPINT(WMOLXT(I,MLOLD),WMOLXT(I,IL),FAC)
 40            CONTINUE
               DENSTY(12,ML)=                                           &
     &              EXPINT(DENSTY(12,MLOLD),DENSTY(12,IL),FAC)
               DENSTY(13,ML)=                                           &
     &              EXPINT(DENSTY(13,MLOLD),DENSTY(13,IL),FAC)
               DENSTY(14,ML)=                                           &
     &              EXPINT(DENSTY(14,MLOLD),DENSTY(14,IL),FAC)
               DENSTY(15,ML)=                                           &
     &              EXPINT(DENSTY(15,MLOLD),DENSTY(15,IL),FAC)
               DENSTY(16,ML)=                                           &
     &              EXPINT(DENSTY(16,MLOLD),DENSTY(16,IL),FAC)
            ELSE

!              NO MORE SPACE IN ARRAYS FOR AN ADDITIONAL LAYER
               WRITE(IPR,'(/3A,/14X,A,2(A,I4))')' FATAL ERROR: ',       &
     &              ' FILE "PARAMS.h" PARAMETER "LAYDIM" MUST',         &
     &              ' BE INCREASED.',' IT SUFFICES TO INCREASE',        &
     &              ' LAYDIM FROM',LAYDIM,' TO',LAYDIM+NLNOV-IMERGE+1
               IF(LJMASS)CALL WRTBUF(FATAL)
               STOP
            ENDIF

!           EXIT LOOP IF ALL NOVAM BOUNDARIES HAVE BEEN INTEGRATED.
            IF(IMERGE.GE.NLNOV)GOTO70

!           INCREMENT NOVAM BOUNDARY INDEX AND START AGAIN
            IMERGE=IMERGE+1
            GOTO 30
         ENDIF

!        LINEARLY INTERPOLATE WATER PARTICLE DENSITIES
         IF(.NOT.LSET)THEN
            FAC=(ZM(IL)-ZM(ML))/(ALTNOV(IMERGE)-ZM(ML))
!           DENSTY(68,IL)
!           =DENSTY(68,ML)+FAC*(EXTNOV(4,IMERGE)-DENSTY(68,ML))
            DO 190 INOV=1,NNOV
               JNOV=INOV-1
               DENSTY(68+JNOV,IL)=DENNOV(INOV,IMERGE)
 190        CONTINUE
            DENSTY(7,IL)=0.0
            RELHUM(IL)=RELHUM(ML)                                       &
     &           +FAC*(RHNOV(IMERGE)-RELHUM(ML))
            T(IL)=T(ML)                                                 &
     &           +FAC*(TNOV(IMERGE)-T(ML))
            TRATIO=273.15/T(IL)
            WH(IL)=.01*RELHUM(IL)*TRATIO*                               &
     &           EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
         ENDIF

!        TRANSLATE LAYER BOUNDARY DATA TO NEW BOUNDARY INDEX.
         ML=ML+1
         ZM(ML)=ZM(IL)
         DENSTY(66,ML)=DENSTY(66,IL)
         DENSTY(67,ML)=DENSTY(67,IL)
         DENSTY(3,ML)=DENSTY(3,IL)
         P(ML)=P(IL)
         T(ML)=T(IL)
         RELHUM(ML)=RELHUM(IL)
         WH(ML)=WH(IL)
         WCO2(ML)=WCO2(IL)
         WO(ML)=WO(IL)
         WN2O(ML)=WN2O(IL)
         WCO(ML)=WCO(IL)
         WCH4(ML)=WCH4(IL)
         WO2(ML)=WO2(IL)
         WHNO3(ML)=WHNO3(IL)
         WNO(ML)=WNO(IL)
         WSO2(ML)=WSO2(IL)
         WNO2(ML)=WNO2(IL)
         WNH3(ML)=WNH3(IL)
         DO 50 I=1,NMOLX
            WMOLXT(I,ML)=WMOLXT(I,IL)
 50      CONTINUE

         DENSTY(68,ML)=DENSTY(68,IL)
         DENSTY(69,ML)=DENSTY(69,IL)
         DENSTY(70,ML)=DENSTY(70,IL)
         DENSTY(71,ML)=DENSTY(71,IL)
         DENSTY(72,ML)=DENSTY(72,IL)
         DENSTY(73,ML)=DENSTY(73,IL)
         DENSTY(74,ML)=DENSTY(74,IL)
         DENSTY(75,ML)=DENSTY(75,IL)
         DENSTY(76,ML)=DENSTY(76,IL)
         DENSTY(77,ML)=DENSTY(77,IL)

         DENSTY(7,ML)=0.0
         DENSTY(12,ML)=DENSTY(12,IL)
         DENSTY(13,ML)=DENSTY(13,IL)
         DENSTY(14,ML)=DENSTY(14,IL)
         DENSTY(15,ML)=DENSTY(15,IL)
         DENSTY(16,ML)=DENSTY(16,IL)
 60   CONTINUE
 70   CONTINUE

!     TRANSLATE REMAINING LAYER BOUNDARY DATA TO NEW BOUNDARY INDEX
      NEWBND=IL
      DO 90 IL=NEWBND,LAYDIM
         ML=ML+1
         ZM(ML)=ZM(IL)
         DENSTY(66,ML)=DENSTY(66,IL)
         DENSTY(67,ML)=DENSTY(67,IL)
         DENSTY(3,ML)=DENSTY(3,IL)
         P(ML)=P(IL)
         T(ML)=T(IL)
         RELHUM(ML)=RELHUM(IL)
         WH(ML)=WH(IL)
         WCO2(ML)=WCO2(IL)
         WO(ML)=WO(IL)
         WN2O(ML)=WN2O(IL)
         WCO(ML)=WCO(IL)
         WCH4(ML)=WCH4(IL)
         WO2(ML)=WO2(IL)
         WHNO3(ML)=WHNO3(IL)
         WNO(ML)=WNO(IL)
         WSO2(ML)=WSO2(IL)
         WNO2(ML)=WNO2(IL)
         WNH3(ML)=WNH3(IL)
         DO 80 I=1,NMOLX
            WMOLXT(I,ML)=WMOLXT(I,IL)
 80      CONTINUE
         DENSTY(68,ML)=DENSTY(68,IL)
         DENSTY(69,ML)=DENSTY(69,IL)
         DENSTY(70,ML)=DENSTY(70,IL)
         DENSTY(71,ML)=DENSTY(71,IL)
         DENSTY(72,ML)=DENSTY(72,IL)
         DENSTY(73,ML)=DENSTY(73,IL)
         DENSTY(74,ML)=DENSTY(74,IL)
         DENSTY(75,ML)=DENSTY(75,IL)
         DENSTY(76,ML)=DENSTY(76,IL)
         DENSTY(77,ML)=DENSTY(77,IL)
         DENSTY(7,ML)=DENSTY(7,IL)
         DENSTY(12,ML)=DENSTY(12,IL)
         DENSTY(13,ML)=DENSTY(13,IL)
         DENSTY(14,ML)=DENSTY(14,IL)
         DENSTY(15,ML)=DENSTY(15,IL)
         DENSTY(16,ML)=DENSTY(16,IL)
 90   CONTINUE

!     MAKE THE LOWEST NOVAM LAYER THE BOTTOM; THAT IS EQUAL TO ZM(1)
      DO 100 IL=1, ML
         IF (ABS(ALTNOV(1)-ZM(IL)).LE.ZTOL) THEN
            IF (IL.EQ.1)RETURN
            GO TO 200
         ENDIF
 100  CONTINUE
 200  CONTINUE
      DO 300 I=IL, ML
         J=I-IL+1
         ZM(J)=ZM(I)
         DENSTY(66,J)=DENSTY(66,I)
         DENSTY(67,J)=DENSTY(67,I)
         DENSTY(3,J)=DENSTY(3,I)
         P(J)=P(I)
         T(J)=T(I)
         RELHUM(J)=RELHUM(I)
         WH(J)=WH(I)
         WCO2(J)=WCO2(I)
         WO(J)=WO(I)
         WN2O(J)=WN2O(I)
         WCO(J)=WCO(I)
         WCH4(J)=WCH4(I)
         WO2(J)=WO2(I)
         WHNO3(J)=WHNO3(I)
         WNO(J)=WNO(I)
         WSO2(J)=WSO2(I)
         WNO2(J)=WNO2(I)
         WNH3(J)=WNH3(I)
         DO 400 II=1,NMOLX
            WMOLXT(II,J)=WMOLXT(II,I)
 400     CONTINUE
         DENSTY(68,J)=DENSTY(68,I)
         DENSTY(69,J)=DENSTY(69,I)
         DENSTY(70,J)=DENSTY(70,I)
         DENSTY(71,J)=DENSTY(71,I)
         DENSTY(72,J)=DENSTY(72,I)
         DENSTY(73,J)=DENSTY(73,I)
         DENSTY(74,J)=DENSTY(74,I)
         DENSTY(75,J)=DENSTY(75,I)
         DENSTY(76,J)=DENSTY(76,I)
         DENSTY(77,J)=DENSTY(77,I)
         DENSTY(7,J)=DENSTY(7,I)
         DENSTY(12,J)=DENSTY(12,I)
         DENSTY(13,J)=DENSTY(13,I)
         DENSTY(14,J)=DENSTY(14,I)
         DENSTY(15,J)=DENSTY(15,I)
         DENSTY(16,J)=DENSTY(16,I)
 300  CONTINUE
      ML=ML-IL+1
      RETURN
      END
