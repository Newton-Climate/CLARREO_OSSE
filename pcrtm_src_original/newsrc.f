      LOGICAL FUNCTION NEWSRC(ZMAX,H1SAV,H2SAV,OZENSV,RNGSAV,           &
     &  BETASV,LENNSV,IPARM,PARM1,PARM2,SRCLAT,SRCLON,TRUEAZ,BETAH2)

!     NEWSRC IS CALLED FOR SOLAR (SINGLE OR MULTIPLE) SCATTER
!     CALCULATIONS (IEMSCT).  IF THE MULTIPLE LINE-OF-SIGHT OPTION
!     IS USED, NEWSRC IS CALLED FOR THE FIRST LINE-OF-SIGHT ONLY.
!     NEWSRC RETURNS TRUE IF THE FOLLOWING TASKS ARE SUCCESSFUL:
!        (1)  RELATIVE SOLAR GEOMETRY PARAMETERS ARE REPLACED BY
!             ABSOLUTE SOLAR GEOMETRY PARAMETERS.
!        (2)  IF THE OBSERVER IS ABOVE THE TOP OF THE ATMOSPHERE,
!             I.E. H1ALT>ZMAX,  THIS ROUTINE DEFINES THE SOLAR/LUNAR
!             ANGLES FOR THE POINT AT WHICH THE LINE-OF-SIGHT PATH
!             ENTERS THE ATMOSPHERE.
!        (3)  IF IMULT EQUALS -1, INDICATING THAT MULTIPLE SCATTERING
!             CALCULATIONS ARE TO BE PERFORMED AT H2ALT INSTEAD
!             H1ALT, THIS ROUTINE RETURNS BETAH2, THE EARTH CENTER
!             ANGLE BETWEEN H1ALT AND H2ALT FOR USE IN ROUTINE H2SRC.
!        (4)  THIS ROUTINE DETERMINES THE SOLAR ANGLES AT H1ALT AND
!             RESETS IPARM IF SOLAR ANGLES ARE DEFINED AT H2ALT,
!             I.E., IPARM=10, 11 OR 12.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       ZMAX     MAXIMUM ATMOSPHERIC PROFILE ALTITUDE [KM].
!       H1SAV    ORIGINAL OBSERVER ALTITUDE [KM]
!       H2SAV    ORIGINAL FINAL ALTITUDE [KM]
!       OZENSV   PATH ZENITH ANGLE AT THE OBSERVER [DEG]
!       RNGSAV   PATH RANGE [KM]
!       BETASV   EARTH CENTER ANGLE SUBTENDED BY PATH [DEG]
!       LENNSV   FLAG EQUAL TO 1 FOR PATHS THROUGH A TANGENT HEIGHT
!       IPARM    SOLAR/LUNAR GEOMETRY SPECIFICATION FLAG
!       PARM1    SENSOR (H1SAV) LAT IF IPARM=0 OR =1 [DEG N OF EQ].
!                TARGET (H2SAV) LAT IF IPARM=10 OR =11 [DEG N OF EQ].
!                SUN REL AZ AT SENSOR (H1SAV) IF IPARM=2 [DEG E OF N].
!                SUN REL AZ AT TARGET (H1SAV) IF IPARM=12 [DEG E OF N].
!       PARM2    SENSOR (H1SAV) LON IF IPARM=0 OR =1 [DEG W OF GNWCH].
!                TARGET (H2SAV) LON IF IPARM=10 OR =11 [DEG W OF GNWCH].
!                SOLAR ZENITH AT SENSOR IF IPARM IS 2 [DEG].
!                SOLAR ZENITH AT TARGET IF IPARM IS 12 [DEG].
!       SRCLAT   SOURCE LATITUDE [DEG NORTH OF EQUATOR].
!       SRCLON   SOURCE LONGITUDE [DEG WEST OF GREENWICH].
!       TRUEAZ   H1SAV TO H2SAV AZM IF IPARM=0 OR =1 [DEG E OF N].
!                H2SAV TO H1SAV AZM IF IPARM=10 OR =11 [DEG E OF N].

!     OUTPUT ARGUMENTS:
!       BETAH2   EARTH CENTER ANGLE BETWEEN H1ALT AND H2SAV.
      DOUBLE PRECISION ZMAX,H1SAV,H2SAV,OZENSV,RNGSAV,BETASV,BETAH2,    &
     &  PARM1,PARM2,SRCLAT,SRCLON,TRUEAZ
      INTEGER LENNSV,IPARM

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /CARD1/
!       MODEL    MODEL ATMOSPHERE INDEX.
!       ITYPE    SLANT PATH TYPE.
!       IEMSCT   RADIATIVE TRANSFER MODE.
!                  0 FOR TRANSMITTANCE
!                  1 FOR THERMAL EMISSION ONLY
!                  2 FOR THERMAL EMISSION PLUS SOLAR SCATTER
!                  3 FOR TRANSMITTED SOLAR IRRADIANCE
!                  4 FOR SOLAR SCATTER ONLY
!       M1       MODEL ATMOSPHERE FOR PRESSURE & TEMPERATURE PROFILES.
!       M2       MODEL ATMOSPHERE FOR H2O PROFILE.
!       M3       MODEL ATMOSPHERE FOR O3 PROFILE.
!       I_RD2C   READ CARD 2C, 2C1, ... IF EQUAL 1; SKIP IF EQUAL TO 0.
!       NOPRNT   PRINT FLAG.
!       MODTRN   MODTRAN BAND MODEL FLAG.
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,I_RD2C,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,I_RD2C,NOPRNT,MODTRN

!     /CARD3/
!       H1ALT    OBSERVER (SENSOR) ALTITUDE [KM].
!       H2ALT    FINAL (TANGENT FOR LIMB PATH) ALTITUDE [KM].
!       OBSZEN   OBSERVER ZENITH ANGLE (H1ALT TO H2ALT) [DEG].
!       HRANGE   DISTANCE FROM H1ALT TO H2ALT [KM].
!       BETA     EARTH CENTER ANGLE BETWEEN H1ALT AND H2ALT [DEG].
!       REARTH   RADIUS OF THE EARTH [KM].
!       HMIN     PATH MINIMUM ALTITUDE [KM].
!       HMAX     PATH MAXIMUM ALTITUDE [KM].
!       CKRANG   MAXIMUM PATH RANGE FOR K-DISTRIBUTION OUTPUT
!                (=0. FOR TOTAL PATH ONLY; <0. FOR ALL RANGES).
!       BCKZEN   ZENITH ANGLE FOR BACKWARD (H2ALT TO H1ALT) PATH [DEG].
!       ANGLEM   LUNAR PHASE ANGLE [0 TO 180 DEG].
!       LENN     PATH LENGTH SWITCH (0=SHORT, 1=LONG).
!       IDAY     DAY OF YEAR [0-366, DEFAULT (0) IS DAY 91].
!       ISOURC   SOURCE FLAG [0=SUN AND 1=MOON].
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
      INTEGER LENN,IDAY,ISOURC,NLOS
      DOUBLE PRECISION H1ALT,H2ALT,OBSZEN,HRANGE,BETA,REARTH,           &
     &  HMIN,HMAX,CKRANG,BCKZEN,ANGLEM
      COMMON/CARD3/H1ALT(MLOS),H2ALT(MLOS),OBSZEN(MLOS),HRANGE(MLOS),   &
     &  BETA(MLOS),HMIN(MLOS),HMAX(MLOS),CKRANG(MLOS),BCKZEN(MLOS),     &
     &  REARTH,ANGLEM,LENN(MLOS),IDAY,ISOURC,NLOS

!     /CNTRL/
!       NSEG     NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
!       THERML   FLAG TO CALCULATE THERMAL SCATTER.
      INTEGER NSEG,ML,MLFLX,IMULT
      LOGICAL THERML
      COMMON/CNTRL/NSEG(0:MLOSP1),ML,MLFLX,IMULT,THERML

!     /SMALL2/
      LOGICAL LSMALL
      COMMON/SMALL2/LSMALL

!     /DCNSTN/
!       DRIGHT   SMALLEST DOUBLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
!       DPDEG    NUMBER OF DEGREES IN ONE RADIAN IN DOUBLE PRECISION
      DOUBLE PRECISION DRIGHT,DPDEG
      COMMON/DCNSTN/DRIGHT,DPDEG

!     /LATLON/
!       H1LAT    ORIGINAL OBSERVER LATITUDE [DEG NORTH OF EQUATOR].
!       H1LON    ORIGINAL OBSERVER LONGITUDE [DEG WEST OF GREENWICH].
!       H2LAT    ORIGINAL TARGET LATITUDE [DEG NORTH OF EQUATOR].
!       H2LON    ORIGINAL TARGET LONGITUDE [DEG WEST OF GREENWICH].
      DOUBLE PRECISION H1LAT,H1LON,H2LAT,H2LON
      COMMON/LATLON/H1LAT,H1LON,H2LAT,H2LON

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
      DOUBLE PRECISION ECA,STORE,H1ALTX,H2ALTX,HMINX,BEGZEN,ENDZEN,     &
     &  DPBETA,DPRANG,BEND,PSINEW,PRM1NW,PRM2NW
      INTEGER IERROR,NPRSAV,LENGTH

!     GEOINP DETERMINES (H1ALT,H2ALT,HMIN,OBSZEN,BCKZEN,LENN).
      NEWSRC=.TRUE.
      BETAH2=0.D0
      IF(ITYPE.EQ.4)THEN

!         USER-DEFINED REFRACTIVE PATH:
          ITYPE=2
      ELSE
          IERROR=0
          NPRSAV=NPR
          NPR=8
          CALL GEOINP(REARTH,ZMAX,H1ALT(1),H2ALT(1),OBSZEN(1),HRANGE(1),&
     &      BETA(1),ITYPE,LENN(1),HMIN(1),BCKZEN(1),IERROR,0)
          NPR=NPRSAV

!         RESET VALUES AND RETURN IF AN INCONSISTENCY WAS DETECTED.
          IF(IERROR.NE.0)THEN
              H1ALT(1)=H1SAV
              H2ALT(1)=H2SAV
              OBSZEN(1)=OZENSV
              HRANGE(1)=RNGSAV
              BETA(1)=BETASV
              LENN(1)=LENNSV
              NEWSRC=.FALSE.
              RETURN
          ENDIF
      ENDIF

!     CHECK FOR RELATIVE SOLAR GEOMETRY PARAMETERS.
      IF(IPARM.EQ.2 .OR. IPARM.EQ.12)THEN

!         RELATIVE SOLAR AZIMUTH AND SOLAR ZENITH ANGLES WERE INPUTS.
!         THE OBSERVER AND SUN ARE BOTH PLACED ON THE EQUATOR WITH THE
!         OBSERVER IS AT 0 DEG LONGITUDE AND THE SUN TO THE WEST.
!         STEP 1:  DETERMINE THE TRUE PATH AZIMUTH
!                  [TRUEAZ = TRUE SOLAR AZIMUTH (270 DEG EAST OF NORTH)
!                          - RELATIVE SOLAR AZIMUTH (PARM1)].
          TRUEAZ=270-PARM1

!         STEP 2:  DETERMINE SOLAR ZENITH EXITING THE ATMOSPHERE
!                  OR UPON HITTING THE HARD EARTH (180-ENDZEN).
          IF(IPARM.EQ.2)THEN

!             OBSERVER / SENSOR (H1) BASED ANGLES:
              H1ALTX=H1SAV
              H1LAT=0.D0
              H1LON=0.D0
          ELSE

!             FINAL ALTITUDE (H2) BASED ANGLES:
              H1ALTX=H2SAV
              H2LAT=0.D0
              H2LON=0.D0
          ENDIF
          BEGZEN=PARM2
          H2ALTX=ZMAX
          LENGTH=1
          IERROR=0
          CALL FINDMN(REARTH,H1ALTX,BEGZEN,H2ALTX,                      &
     &      LENGTH,0,HMINX,ENDZEN,IERROR,.FALSE.)

!         STEP 3:  SET SOLAR LONGITUDE TO SOLAR ZENITH PLUS BENDING.
          IF(IERROR.EQ.0)THEN

!             THE SOLAR PATH FROM H1SAV PASSES THROUGH THE ATMOSPHERE:
              LSMALL=.FALSE.
              IF(H1ALTX.GT.ZMAX)THEN

!                 REFERENCE ALTITUDE ABOVE TOP-OF-ATMOSPHERE:
                  H1ALTX=ZMAX
                  BEGZEN=ENDZEN
              ENDIF
              CALL RFPATH(REARTH,H1ALTX,H2ALTX,BEGZEN,LENGTH,           &
     &          HMINX,.FALSE.,ENDZEN,DPBETA,DPRANG,BEND)
              SRCLON=PARM2+BEND
          ELSE

!             SOLAR PATH IS ABOVE ATMOSPHERE OR INTERSECTS EARTH:
              SRCLON=PARM2
              IERROR=0
          ENDIF

!         STEP 4:  WRITE OUT CHANGE.
          IF(.NOT.LJMASS)THEN
              IF(IPARM.EQ.2)THEN
                  WRITE(IPR,'(/A,/(F13.5,A))')' REPLACING'              &
     &              //' IPARM=2 SOLAR/LUNAR GEOMETRY PARAMETERS',       &
     &              PARM1,'   RELATIVE SOLAR/LUNAR AZIMUTH AT H1ALT'    &
     &              //' (DEG EAST OF NORTH)',                           &
     &              PARM2,'   SOLAR/LUNAR ZENITH ANGLE AT H1ALT (DEG)'
                  WRITE(IPR,'(4(/2A),/(F13.5,2A))')'      WITH',        &
     &              ' IPARM=0 SOLAR/LUNAR GEOMETRY PARAMETERS',         &
     &              '      0.00000   OBSERVER (H1ALT) LATITUDE',        &
     &              ' (DEG NORTH OF EQUATOR)',                          &
     &              '      0.00000   OBSERVER (H1ALT) LONGITUDE',       &
     &              ' (DEG WEST OF GREENWICH)',                         &
     &              '      0.00000   SOURCE LATITUDE',                  &
     &              ' (DEG NORTH OF EQUATOR)',                          &
     &              SRCLON,'   SOURCE LONGITUDE',                       &
     &              ' (DEG WEST OF GREENWICH)',                         &
     &              TRUEAZ,'   TRUE PATH AZIMUTH FROM H1ALT TO H2ALT',  &
     &              ' (DEG EAST OF NORTH)'
              ELSE
                  WRITE(IPR,'(/A,/(F13.5,2A))')' REPLACING'//           &
     &              ' IPARM=12 SOLAR/LUNAR GEOMETRY PARAMETERS',        &
     &              PARM1,'   RELATIVE SOLAR/LUNAR AZIMUTH AT H2ALT',   &
     &              ' (DEG EAST OF NORTH)',                             &
     &              PARM2,'   SOLAR/LUNAR ZENITH ANGLE AT H2ALT (DEG)'
                  WRITE(IPR,'(4(/2A),/(F13.5,2A))')'      WITH',        &
     &              ' IPARM=10 SOLAR/LUNAR GEOMETRY PARAMETERS',        &
     &              '      0.00000   TARGET (H2ALT) LATITUDE',          &
     &              ' (DEG NORTH OF EQUATOR)',                          &
     &              '      0.00000   TARGET (H2ALT) LONGITUDE',         &
     &              ' (DEG WEST OF GREENWICH)',                         &
     &              '      0.00000   SOURCE LATITUDE',                  &
     &              ' (DEG NORTH OF EQUATOR)',                          &
     &              SRCLON,'   SOURCE LONGITUDE',                       &
     &              ' (DEG WEST OF GREENWICH)',                         &
     &              TRUEAZ,'   TRUE PATH AZIMUTH FROM H2ALT TO H1ALT',  &
     &              ' (DEG EAST OF NORTH)'
              ENDIF
          ENDIF

!         STEP 5:  SET ABSOLUTE SOLAR GEOMETRY VALUES
          IPARM=IPARM-2
          PARM1=0.D0
          PARM2=0.D0
          SRCLAT=0.D0
      ELSEIF(IPARM.LT.2)THEN

!         SAVE OBSERVER / SENSOR LATITUDE AND LONGITUDE:
          H1LAT=PARM1
          H1LON=PARM2
      ELSE

!         SAVE TARGET LATITUDE AND LONGITUDE:
          H2LAT=PARM1
          H2LON=PARM2
      ENDIF

!     UPDATE ON PARAMETER VALUES:
!       PARM1    SENSOR (H1SAV) LAT IF IPARM=(0,1,2) [DEG N OF EQ].
!                TARGET (H2SAV) LAT IF IPARM=(10,11,12) [DEG N OF EQ].
!       PARM2    SENSOR (H1SAV) IF IPARM=(0,1,2) LON [DEG W OF GNWCH].
!                TARGET (H2SAV) IF IPARM=(10,11,12) LON [DG W OF GNWCH].
!       TRUEAZ   H1SAV TO H2SAV AZM IF IPARM=(0,1,2) [DEG E OF N].
!                H2SAV TO H1SAV AZM IF IPARM=(10,11,12) [DEG E OF N].

!     IF H1ALT IS BELOW THE TOP-OF-ATMOSPHERE AND SOLAR GEOMETRY IS
!     DEFINED RELATIVE TO SENSOR (IPARM < 10) RETURN UNLESS IMULT=-1
!     (IF IMULT=-1, THE EARTH CENTER ANGLE TO H2ALT MUST BE DETERMINED):
      IF(H1SAV.LE.ZMAX .AND. IMULT.NE.-1 .AND. IPARM.LT.10)RETURN

!     CALCULATE THE EARTH CENTER ANGLE, ECA, BETWEEN THE ORIGINAL
!     OBSERVER ALTITUDE, H1SAV, AND THE NEW OBSERVER ALTITUDE, H1ALT.
!     ASSUME NO REFRACTION ABOVE, H1ALT, THE TOP OF THE ATMOSPHERE.
      STORE=SIN(OBSZEN(1)/DPDEG)/(REARTH+H1SAV)
      IF(H1SAV.GT.ZMAX)THEN
          ECA=180-OBSZEN(1)-DPDEG*ASIN((REARTH+ZMAX)*STORE)
      ELSE
          ECA=0.D0
      ENDIF

!     IF MULTIPLE SCATTERING IS TO BE CALCULATED AT THE LAT/LONG OF H2ALT
!     OR SOLAR GEOMETRY PARAMETERS ARE DEFINED AT THE LAT/LONG OF H2ALT,
!     DETERMINE THE EARTH CENTER ANGLE BETWEEN H1ALT AND H2SAV.
      IF(IMULT.EQ.-1 .OR. IPARM.GE.10)THEN
          IF(ITYPE.EQ.4)THEN
              BETAH2=BETASV
          ELSE
              LSMALL=.FALSE.
              IF(ITYPE.EQ.3 .AND. H2SAV.NE.0.D0)THEN

!                 H2ALT IS THE TANGENT POINT (CASE 3B).
                  H1ALTX=HMIN(1)
                  H2ALTX=H1ALT(1)
                  BEGZEN=DBLE(90.)
                  ENDZEN=OBSZEN(1)
                  LENGTH=0
              ELSE

!                 H2ALT IS THE FINAL ALTITUDE (NO MORE THAN ZMAX).
                  BEGZEN=OBSZEN(1)
                  H1ALTX=H1ALT(1)
                  H2ALTX=H2ALT(1)
                  ENDZEN=BCKZEN(1)
                  LENGTH=LENN(1)
              ENDIF
              CALL RFPATH(REARTH,H1ALTX,H2ALTX,BEGZEN,LENGTH,           &
     &          HMIN(1),.FALSE.,ENDZEN,DPBETA,DPRANG,BEND)
              IF(H2SAV.GT.ZMAX)THEN

!                 H2ALT WAS REDUCED TO THE TOP OF ATMOSPHERE.  ADD TO
!                 THE EARTH CENTER ANGLE THE H2ALT TO H2SAV CONTRIBUTION
                  BETAH2=DPBETA+180-ENDZEN-DPDEG*ASIN                   &
     &              ((REARTH+ZMAX)*SIN(ENDZEN/DPDEG)/(REARTH+H2SAV))
              ELSE
                  BETAH2=DPBETA
              ENDIF
          ENDIF

!         DETERMINE OBSERVER / SENSOR LATITUDE AND LONGITUDE:
          IF(IPARM.GE.10)THEN

!             CONVERT TARGET (H2ALT) BASED SOLAR GEOMETRY PARAMETERS
!             TO OBSERVER (H1ALT) BASED SOLAR GEOMETRY PARAMETERS.
              IF(.NOT.LJMASS)THEN
                  WRITE(IPR,'(/A,/(F13.5,2A))')' CONVERT FROM IPARM=10' &
     &              //' TARGET (H2ALT) SOLAR/LUNAR GEOMETRY PARAMETERS',&
     &              PARM1,'   TARGET (H2ALT) LATITUDE',                 &
     &              ' (DEG NORTH OF EQUATOR)',                          &
     &              PARM2,'   TARGET (H2ALT) LONGITUDE',                &
     &             ' (DEG WEST OF GREENWICH)'
                  WRITE(IPR,'(F13.5,2X,A)')TRUEAZ,' TRUE PATH'//        &
     &              ' AZIMUTH FROM H2ALT TO H1ALT (DEG EAST OF NORTH)'
              ENDIF
              CALL LOCATE(PARM1,PARM2,TRUEAZ,BETAH2,PRM1NW,PRM2NW)
              H1LAT=PRM1NW
              H1LON=PRM2NW
              CALL PSIECA(PRM1NW,PRM2NW,PARM1,PARM2,TRUEAZ,BETAH2)
              PARM1=PRM1NW
              PARM2=PRM2NW
              IF(.NOT.LJMASS)WRITE(IPR,'(/A,/(F13.5,A))')' TO IPARM=0'  &
     &          //' OBSERVER (H1ALT) SOLAR/LUNAR GEOMETRY PARAMETERS',  &
     &          PARM1,'   OBSERVER (H1ALT) LATITUDE'                    &
     &          //' (DEG NORTH OF EQUATOR)',                            &
     &          PARM2,'   OBSERVER (H1ALT) LONGITUDE'                   &
     &          //' (DEG WEST OF GREENWICH)',                           &
     &          TRUEAZ,'   TRUE PATH AZIMUTH FROM H1ALT'                &
     &          //' TO H2ALT (DEG EAST OF NORTH)'
          ELSE

!             DETERMINE H2LAT AND H2LON SINCE IMULT=-1 AND IPARM<10:
              CALL LOCATE(PARM1,PARM2,TRUEAZ,BETAH2+ECA,PRM1NW,PRM2NW)  &
     &
              H2LAT=PRM1NW
              H2LON=PRM2NW
          ENDIF

!     UPDATE ON PARAMETER VALUES:
!       PARM1    SENSOR (H1SAV) LAT IF IPARM=(0,1,2) [DEG N OF EQ].
!                SENSOR (H1ALT) LAT IF IPARM=(10,11,12) [DEG N OF EQ].
!       PARM2    SENSOR (H1SAV) IF IPARM=(0,1,2) LON [DEG W OF GNWCH].
!                SENSOR (H1ALT) IF IPARM=(10,11,12) LON [DG W OF GNWCH].
!       TRUEAZ   H1SAV TO H2SAV AZM IF IPARM=(0,1,2) [DEG E OF N].
!                H1ALT TO H2SAV AZM IF IPARM=(10,11,12) [DEG E OF N].
          IF(H1SAV.LE.ZMAX)THEN
              H1ALT(1)=H1SAV
              H2ALT(1)=H2SAV
              OBSZEN(1)=OZENSV
              HRANGE(1)=RNGSAV
              BETA(1)=BETASV
              LENN(1)=LENNSV
              IPARM=0
              RETURN
          ENDIF
      ENDIF

!     SLANT PATH BRANCHING:
      IF(ITYPE.EQ.2)THEN

!         IF PATH WAS DEFINED BY CASE 2B (H1ALT,OBSZEN,HRANGE)
!         DETERMINE THE REDUCED HRANGE.  IF PATH GEOMETRY WAS DEFINED BY
!         CASE 2C (H1ALT,H2ALT,HRANGE) OR CASE 2D (H1ALT,H2ALT,BETA),
!         CHANGE TO CASE 2A (H1ALT,H2ALT,OBSZEN).
          BETA(1)=0.D0
          HRANGE(1)=0.D0
          IF(BETASV.EQ.0.D0 .AND. RNGSAV.NE.0.D0 .AND. OZENSV.NE.0.D0)  &
     &      HRANGE(1)=RNGSAV-SIN(ECA/DPDEG)/STORE
      ENDIF

!     WRITE OUT CHANGES IN GEOMETRY:
      IF(.NOT.LJMASS)THEN
          WRITE(IPR,'(/A,//2A,/2A,/A,2F13.5)')' OBSERVER ALTITUDE'      &
     &    //' WAS LOWERED TO THE TOP OF THE ATMOSPHERE:',               &
     &    ' NAME                   DESCRIPTION                ',        &
     &    '      INPUT      REVISION',                                  &
     &    ' -----   ------------------------------------------',        &
     &    '   ----------   ----------',                                 &
     &    ' H1ALT   OBSERVER ALTITUDE (KM)                    ',        &
     &    H1SAV,H1ALT(1)
          IF(ITYPE.EQ.3)THEN
              WRITE(IPR,'(A,21X,2F13.5)')                               &
     &          ' H2ALT   TANGENT ALTITUDE (KM)',H2SAV,H2SAV
          ELSE
              WRITE(IPR,'(A,21X,2F13.5)')                               &
     &          ' H2ALT   FINAL ALTITUDE (KM)  ',H2SAV,H2SAV
          ENDIF
          WRITE(IPR,'(A,18X,3(2F13.5,/A,18X),2I13)')                    &
     &      ' ANGLE   OBSERVER ZENITH (DEG)   ',OZENSV,OBSZEN(1),       &
     &      ' HRANGE  SLANT HRANGE (KM)       ',RNGSAV,HRANGE(1),       &
     &      ' BETA    EARTH CENTER ANGLE (DEG)',BETASV,BETA(1),         &
     &      ' LENN    SHORT/LONG PATH FLAG    ',LENNSV,LENN(1)
      ENDIF

!     BRANCH BASED ON SOLAR/LUNAR GEOMETRY SPECIFICATION FLAG:
      IF(IPARM.LE.2)THEN

!         OBSERVER LATITUDE, LONGITUDE, AND PATH AZIMUTH WERE INPUTS.
!         SAVE ORIGINAL VALUES AND DETERMINE THEIR VALUES AT H1ALT.
          H1LAT=PARM1
          H1LON=PARM2
          CALL LOCATE(PARM1,PARM2,TRUEAZ,ECA,PRM1NW,PRM2NW)
          CALL PSIECA(PRM1NW,PRM2NW,PARM1,PARM2,PSINEW,ECA)
          PSINEW=PSINEW+180
          IF(PSINEW.GE.360.D0)PSINEW=PSINEW-360
          IF(.NOT.LJMASS)WRITE(IPR,'((A,2F13.5))')                      &
     &      ' PARM1   OBSERVER LATITUDE (DEG NORTH OF EQUATOR)  ',      &
     &      PARM1,PRM1NW,                                               &
     &      ' PARM2   OBSERVER LONGITUDE (DEG WEST OF GREENWICH)',      &
     &      PARM2,PRM2NW,                                               &
     &      ' TRUEAZ  PATH AZIMUTH (DEG EAST OF NORTH)          ',      &
     &      TRUEAZ,PSINEW
          PARM1=PRM1NW
          PARM2=PRM2NW
          TRUEAZ=PSINEW
      ELSE

!         DETERMINE THE TRUE SENSOR LATITUDE AND LONGITUDE:
          CALL LOCATE(PARM1,PARM2,180+TRUEAZ,ECA,PRM1NW,PRM2NW)
          H1LAT=PRM1NW
          H1LON=PRM2NW
      ENDIF

!     UPDATE ON PARAMETER VALUES:
!       PARM1    SENSOR (H1ALT) LAT [DEG N OF EQ].
!       PARM2    SENSOR (H1ALT) LON [DEG W OF GNWCH].
!       TRUEAZ   H1ALT TO H2SAV AZM [DEG E OF N].

!     SET H2ALT TO ITS ORIGINAL VALUE.
      H2ALT(1)=H2SAV

!     REPLACE PATH GEOMETRY INPUTS:
      H1SAV=H1ALT(1)
      OZENSV=OBSZEN(1)
      RNGSAV=HRANGE(1)
      BETASV=BETA(1)
      LENNSV=LENN(1)
      IPARM=0

!     RETURN TO GEODRV:
      RETURN
      END