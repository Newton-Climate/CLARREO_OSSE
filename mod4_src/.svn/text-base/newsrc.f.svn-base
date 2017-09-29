      SUBROUTINE NEWSRC(H1SAV,H2SAV,ANGSAV,RNGSAV,BETASV,LENNSV,        &
     &  IPARM,PARM1,PARM2,PARM3,PARM4,PSIPO,BETAH2,IERROR)

!     THIS ROUTINE PERFORMS THE FOLLOWING TASKS:
!          (1)  RELATIVE SOLAR GEOMETRY PARAMETERS ARE REPLACED BY
!               ABSOLUTE SOLAR GEOMETRY PARAMETERS.
!          (2)  IF THE OBSERVER IS ABOVE THE TOP OF THE ATMOSPHERE,
!               I.E. H1>ZMAX,  THIS ROUTINE DEFINES THE SOLAR/LUNAR
!               ANGLES FOR THE POINT AT WHICH THE LINE-OF-SIGHT PATH
!               ENTERS THE ATMOSPHERE.
!          (3)  IF IMULT EQUALS -1, INDICATING THAT MULTIPLE SCATTERING
!               CALCULATIONS ARE TO BE PERFORMED AT H2 INSTEAD H1, THIS
!               ROUTINE RETURNS BETAH2, THE EARTH CENTER ANGLE BETWEEN
!               H1 AND H2 FOR USE IN ROUTINE H2SRC.
!          (4)  IF SOLAR ANGLES ARE DEFINED AT H2 (IPARM=10, 11 OR 12),
!               THIS ROUTINE DETERMINES THE SOLAR ANGLES AT H1 AND
!               RESETS IPARM.

!     DECLARE INPUTS/OUTPUTS
!       H1SAV    OBSERVER ALTITUDE [KM]
!       H2SAV    FINAL ALTITUDE [KM]
!       ANGSAV   PATH ZENITH ANGLE AT THE OBSERVER [DEG]
!       RNGSAV   PATH RANGE [KM]
!       BETASV   EARTH CENTER ANGLE SUBTENDED BY PATH [DEG]
!       LENNSV   FLAG EQUAL TO 1 FOR PATHS THROUGH A TANGENT HEIGHT
!       IPARM    SOLAR/LUNAR GEOMETRY SPECIFICATION FLAG
!       PARM1    OBSERVER (H1) LATITUDE [DEG NORTH OF EQUATOR] IF IPARM
!                IS 0 OR 1; TARGET (H2) LATITUDE IF IPARM IS 10 OR 11.
!                RELATIVE SOLAR AZIMUTH [DEG EAST OF NORTH] AT OBSERVER
!                (H1) IF IPARM IS 2 AND AT TARGET (H2) IF IPARM IS 12.
!       PARM2    OBSERVER LONGITUDE [DEG WEST OF GREENWICH]  IF IPARM
!                IS 0 OR 1; TARGET LONGITUDE IF IPARM IS 10 OR 11.
!                SOLAR ZENITH ANGLE [DEG] AT OBSERVER IF IPARM IS 2 AND
!                AT TARGET IF IPARM IS 12.
!       PARM3    SOURCE LATITUDE [DEG NORTH OF EQUATOR].
!       PARM4    SOURCE LONGITUDE [DEG WEST OF GREENWICH].
!       PSIPO    PATH TRUE AZIMUTH ANGLE [DEG EAST OF NORTH]
!       BETAH2   EARTH CENTER ANGLE BETWEEN H1 AND H2
      REAL H1SAV,H2SAV,ANGSAV,RNGSAV,BETASV,PARM1,PARM2,PARM3,PARM4,    &
     &  PSIPO,BETAH2
      INTEGER LENNSV,IPARM

!     PARAMETERS:
      INCLUDE 'ERROR.h'

!     COMMONS:

!     /CARD1/
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,MODTRN

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

!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT
      REAL RE,ZMAX
      INTEGER IPATH
      COMMON/PARMTR/RE,ZMAX,IPATH
      INCLUDE 'IFIL.h'
      LOGICAL LSMALL
      COMMON/SMALL2/LSMALL

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     DECLARE LOCAL VARIABLES
      REAL HMIN,PHI,STORE,ECA,PSINEW,PRM1NW,PRM2NW,H2MX
      DOUBLE PRECISION DPH1,DPH2,DPANG,DPPHI,DPHMIN,DPBETA,DPRANG,DPBEND
      INTEGER IERROR,NPRSAV,IAMT,LENGTH

!     DEFINE DATA
      LOGICAL LOGLOS
      INTEGER ISLCT,ITER
      LOGICAL LWARN
      DATA LOGLOS/.FALSE./,ISLCT/0/,ITER/0/,LWARN/.FALSE./

!     CALL GEOINP TO DETERMINE THE SET (H1,H2,HMIN,ANGLE,PHI,LENN,H2MX).
      IERROR=0
      NPRSAV=NPR
      NPR=2
      CALL GEOINP(H1,H2,ANGLE,RANGE,BETA,ITYPE,                         &
     &  LENN,HMIN,PHI,IERROR,ISLCT,LOGLOS,H2MX)
      NPR=NPRSAV

!     RESET VALUES AND RETURN IF AN INCONSISTENCY WAS DETECTED.
      BETAH2=0.
      IF(IERROR.NE.0)THEN
          H1=H1SAV
          H2=H2SAV
          ANGLE=ANGSAV
          RANGE=RNGSAV
          BETA=BETASV
          LENN=LENNSV
          RETURN
      ENDIF

!     CHECK FOR RELATIVE SOLAR GEOMETRY PARAMETERS.
      IF(IPARM.EQ.2 .OR. IPARM.EQ.12)THEN

!         RELATIVE SOLAR AZIMUTH AND SOLAR ZENITH ANGLES WERE
!         INPUTS.  ASSUME THAT THE OBSERVER AND SUN ARE BOTH ON
!         THE EQUATOR, THAT THE OBSERVER IS AT 0 DEG LONGITUDE,
!         AND THAT THE SUN IS WEST OF THE OBSERVER.
!         STEP 1:  DETERMINE THE TRUE PATH AZIMUTH [PSIPO = TRUE SOLAR
!                  AZIMUTH (270 DG) - RELATIVE SOLAR AZIMUTH (PARM1)].
          PSIPO=270.-PARM1

!         STEP 2:  DETERMINE SOLAR ZENITH EXITING THE ATMOSPHERE
!                  OR UPON HITTING THE HARD EARTH (180-DPPHI).
          IF(IPARM.EQ.2)THEN
              DPH1=DBLE(H1SAV)
          ELSE
              DPH1=DBLE(H2MX)
          ENDIF
          DPANG=DBLE(PARM2)
          DPH2=DBLE(ZMAX)
          LENGTH=1
          IERROR=0
          CALL DPFNMN(DPH1,DPANG,DPH2,LENGTH,                           &
     &      ITER,DPHMIN,DPPHI,IERROR,LWARN)

!         STEP 3:  SET SOLAR LONGITUDE TO SOLAR ZENITH PLUS BENDING.
          IF(IERROR.EQ.0)THEN

!             THE SOLAR PATH FROM H1SAV PASSES THROUGH THE ATMOSPHERE.
              IAMT=2
              LSMALL=.FALSE.
              IF(DPH1.GT.DBLE(ZMAX))THEN
                  DPH1=DBLE(ZMAX)
                  DPANG=DPPHI
              ENDIF
              CALL DPRFPA(DPH1,DPH2,DPANG,DPPHI,LENGTH,                 &
     &          DPHMIN,IAMT,DPBETA,DPRANG,DPBEND)
              PARM4=PARM2+REAL(DPBEND)
          ELSE

!             THE SOLAR PATH FROM H1SAV IS ABOVE THE ATMOSPHERE.
              PARM4=PARM2
              IERROR=0
          ENDIF

!         STEP 4:  WRITE OUT CHANGE.
          IF(.NOT.LJMASS)THEN
              IF(IPARM.EQ.2)THEN
                  WRITE(IPR,'(/A,/(F13.5,2A))')' REPLACING'//           &
     &              ' IPARM=2 SOLAR/LUNAR GEOMETRY PARAMETERS',         &
     &              PARM1,'   RELATIVE SOLAR/LUNAR AZIMUTH AT H1',      &
     &              ' (DEG EAST OF NORTH)',                             &
     &              PARM2,'   SOLAR/LUNAR ZENITH ANGLE AT H1 (DEG)'
                  WRITE(IPR,'(4(/2A),/(F13.5,2A))')'      WITH',        &
     &              ' IPARM=0 SOLAR/LUNAR GEOMETRY PARAMETERS',         &
     &              '      0.00000   OBSERVER (H1) LATITUDE',           &
     &              ' (DEG NORTH OF EQUATOR)',                          &
     &              '      0.00000   OBSERVER (H1) LONGITUDE',          &
     &              ' (DEG WEST OF GREENWICH)',                         &
     &              '      0.00000   SOURCE LATITUDE',                  &
     &              ' (DEG NORTH OF EQUATOR)',                          &
     &              PARM4,'   SOURCE LONGITUDE',                        &
     &              ' (DEG WEST OF GREENWICH)',                         &
     &              PSIPO,'   TRUE PATH AZIMUTH FROM H1 TO H2',         &
     &              ' (DEG EAST OF NORTH)'
              ELSE
                  WRITE(IPR,'(/A,/(F13.5,2A))')' REPLACING'//           &
     &              ' IPARM=12 SOLAR/LUNAR GEOMETRY PARAMETERS',        &
     &              PARM1,'   RELATIVE SOLAR/LUNAR AZIMUTH AT H2',      &
     &              ' (DEG EAST OF NORTH)',                             &
     &              PARM2,'   SOLAR/LUNAR ZENITH ANGLE AT H2 (DEG)'
                  WRITE(IPR,'(4(/2A),/(F13.5,2A))')'      WITH',        &
     &              ' IPARM=10 SOLAR/LUNAR GEOMETRY PARAMETERS',        &
     &              '      0.00000   TARGET (H2) LATITUDE',             &
     &              ' (DEG NORTH OF EQUATOR)',                          &
     &              '      0.00000   TARGET (H2) LONGITUDE',            &
     &              ' (DEG WEST OF GREENWICH)',                         &
     &              '      0.00000   SOURCE LATITUDE',                  &
     &              ' (DEG NORTH OF EQUATOR)',                          &
     &              PARM4,'   SOURCE LONGITUDE',                        &
     &              ' (DEG WEST OF GREENWICH)',                         &
     &              PSIPO,'   TRUE PATH AZIMUTH FROM H2 TO H1',         &
     &              ' (DEG EAST OF NORTH)'
              ENDIF
          ENDIF

!         STEP 5:  SET ABSOLUTE SOLAR GEOMETRY VALUES
          IPARM=IPARM-2
          PARM1=0.
          PARM2=0.
          PARM3=0.
      ENDIF

!     RETURN IF H1 IS BELOW THE TOP OF THE ATMOSPHERE
!     UNLESS IMULT EQUALS -1 (IF IMULT EQUALS -1, THE
!     EARTH CENTER ANGLE TO H2 MUST BE DETERMINED).
      IF(H1.LE.ZMAX .AND. IMULT.NE.-1 .AND. IPARM.LE.2)RETURN

!     IF MULTIPLE SCATTERING IS TO BE CALCULATED AT H2
!     OR SOLAR GEOMETRY PARAMETER ARE DEFINED AT H2,
!     DETERMINE THE EARTH CENTER ANGLE BETWEEN H1 AND H2.
      IF(IMULT.EQ.-1 .OR. IPARM.GE.10)THEN
          IAMT=2
          LSMALL=.FALSE.
          DPHMIN=DBLE(HMIN)
          IF(ITYPE.EQ.3 .AND. H2SAV.NE.0.)THEN

!             H2 IS THE TANGENT POINT (CASE 3B).
              DPH1=DPHMIN
              DPH2=DBLE(H1)
              DPANG=DBLE(90.)
              DPPHI=DBLE(ANGLE)
              LENGTH=0
          ELSE

!             H2 IS THE FINAL ALTITUDE (NO MORE THAN ZMAX).
              DPH1=DBLE(H1)
              DPH2=DBLE(H2)
              DPANG=DBLE(ANGLE)
              DPPHI=DBLE(PHI)
              LENGTH=LENN
          ENDIF
          CALL DPRFPA(DPH1,DPH2,DPANG,DPPHI,LENGTH,                     &
     &      DPHMIN,IAMT,DPBETA,DPRANG,DPBEND)
          BETAH2=REAL(DPBETA)
          IF(H2MX.GT.H2)THEN

!             H2 WAS REDUCED TO THE TOP OF THE ATMOSPHERE.  ADD TO
!             THE EARTH CENTER ANGLE THE H2 TO H2MX CONTRIBUTION.
              STORE=SIN(PHI/DEG)/(RE+H2SAV)
              BETAH2=BETAH2+(180.-(DEG*ASIN((RE+ZMAX)*STORE)+PHI))
          ENDIF
          IF(IPARM.GE.10)THEN

!             CONVERT FROM TARGET (H2) BASED SOLAR GEOMETRY PARAMETERS
!             TO OBSERVER (H1) BASED SOLAR GEOMETRY PARAMETERS.
              IF(.NOT.LJMASS)THEN
                  WRITE(IPR,'(/A,/(F13.5,2A))')' CONVERT FROM IPARM=10' &
     &              //' TARGET (H2) SOLAR/LUNAR GEOMETRY PARAMETERS',   &
     &              PARM1,'   TARGET (H2) LATITUDE',                    &
     &              ' (DEG NORTH OF EQUATOR)',                          &
     &              PARM2,'   TARGET (H2) LONGITUDE',                   &
     &              ' (DEG WEST OF GREENWICH)'
                  WRITE(IPR,'(F13.5,2X,A)')PSIPO,' TRUE PATH'//         &
     &              ' AZIMUTH FROM H2 TO H1 (DEG EAST OF NORTH)'
              ENDIF
              CALL LOCATE(PARM1,PARM2,PSIPO,BETAH2,PRM1NW,PRM2NW)
              CALL PSIECA(PRM1NW,PRM2NW,PARM1,PARM2,PSIPO,BETAH2)
              PARM1=PRM1NW
              PARM2=PRM2NW
              IF(.NOT.LJMASS)WRITE(IPR,'(/A,/(F13.5,A))')' TO IPARM=0'  &
     &          //' OBSERVER (H1) SOLAR/LUNAR GEOMETRY PARAMETERS',     &
     &          PARM1,'   OBSERVER (H1) LATITUDE'//                     &
     &          ' (DEG NORTH OF EQUATOR)',                              &
     &          PARM2,'   OBSERVER (H1) LONGITUDE'//                    &
     &          ' (DEG WEST OF GREENWICH)',                             &
     &          PSIPO,'   TRUE PATH AZIMUTH FROM H1 TO H2'//            &
     &          ' (DEG EAST OF NORTH)'
          ENDIF
          IF(H1SAV.LE.ZMAX)THEN
              H1=H1SAV
              H2=H2SAV
              ANGLE=ANGSAV
              RANGE=RNGSAV
              BETA=BETASV
              LENN=LENNSV
              IPARM=0
              RETURN
          ENDIF
      ENDIF

!     SET H2 TO ITS ORIGINAL VALUE
      H2=H2SAV

!     CALCULATE THE EARTH CENTER ANGLE, ECA, BETWEEN THE ORIGINAL
!     OBSERVER ALTITUDE, H1SAV, AND THE NEW OBSERVER ALTITUDE, H1.
!     ASSUME NO REFRACTION ABOVE, H1, THE TOP OF THE ATMOSPHERE.
      STORE=SIN(ANGLE/DEG)/(RE+H1SAV)
      ECA=180.-(DEG*ASIN((RE+H1)*STORE)+ANGLE)
      IF(ECA.LT.0.)ECA=0.

!     SLANT PATH BRANCHING
      IF(ITYPE.EQ.2)THEN

!         IF PATH GEOMETRY WAS DEFINED BY CASE 2C (H1,H2,RANGE) OR
!         CASE 2D (H1,H2,BETA), CHANGE TO CASE 2A (H1,H2,ANGLE).
!         IF PATH WAS DEFINED BY CASE 2B (H1,ANGLE,RANGE) DETERMINE
!         THE NEW RANGE.
          BETA=0.
          RANGE=0.
          IF(BETASV.EQ.0. .AND. RNGSAV.NE.0. .AND. ANGSAV.NE.0.)        &
     &      RANGE=RNGSAV-SIN(ECA/DEG)/STORE
      ENDIF

!     WRITE OUT CHANGES IN GEOMETRY.
      IF(.NOT.LJMASS)THEN
          WRITE(IPR,'(/A,//2A,/2A,/A,2F13.5)')' OBSERVER ALTITUDE'      &
     &    //' WAS LOWERED TO THE TOP OF THE ATMOSPHERE:',               &
     &    ' NAME                   DESCRIPTION                ',        &
     &    '      INPUT      REVISION',                                  &
     &    ' -----   ------------------------------------------',        &
     &    '   ----------   ----------',                                 &
     &    ' H1      OBSERVER ALTITUDE (KM)                    ',        &
     &    H1SAV,H1
          IF(ITYPE.EQ.3)THEN
              WRITE(IPR,'(2A,2F13.5)')' H2      ',                      &
     &          'TANGENT ALTITUDE (KM)                     ',H2SAV,H2
          ELSE
              WRITE(IPR,'(2A,2F13.5)')' H2      ',                      &
     &          'FINAL ALTITUDE (KM)                       ',H2SAV,H2
          ENDIF
          WRITE(IPR,'(3(A,2F13.5,/),A,2I13)')                           &
     &      ' ANGLE   OBSERVER ZENITH (DEG)                     ',      &
     &      ANGSAV,ANGLE,                                               &
     &      ' RANGE   SLANT RANGE (KM)                          ',      &
     &      RNGSAV,RANGE,                                               &
     &      ' BETA    EARTH CENTER ANGLE (DEG)                  ',      &
     &      BETASV,BETA,                                                &
     &      ' LENN    SHORT/LONG PATH FLAG                      ',      &
     &      LENNSV,LENN
      ENDIF

!     BRANCH BASED ON SOLAR/LUNAR GEOMETRY SPECIFICATION FLAG
      IF(IPARM.LT.2)THEN

!         OBSERVER LATITUDE, LONGITUDE, AND PATH AZIMUTH WERE
!         INPUTS.  DETERMINE THEIR VALUES AT THE NEW H1.
          CALL LOCATE(PARM1,PARM2,PSIPO,ECA,PRM1NW,PRM2NW)
          CALL PSIECA(PRM1NW,PRM2NW,PARM1,PARM2,PSINEW,ECA)
          PSINEW=PSINEW+180.
          IF(PSINEW.GE.360.)PSINEW=PSINEW-360.
          IF(.NOT.LJMASS)WRITE(IPR,'((A,2F13.5))')                      &
     &      ' PARM1   OBSERVER LATITUDE (DEG NORTH OF EQUATOR)  ',      &
     &      PARM1,PRM1NW,                                               &
     &      ' PARM2   OBSERVER LONGITUDE (DEG WEST OF GREENWICH)',      &
     &      PARM2,PRM2NW,                                               &
     &      ' PSIPO   PATH AZIMUTH (DEG EAST OF NORTH)          ',      &
     &      PSIPO,PSINEW
          PARM1=PRM1NW
          PARM2=PRM2NW
          PSIPO=PSINEW
      ENDIF

!     REPLACE PATH GEOMETRY INPUTS
      H1SAV=H1
      H2SAV=H2
      ANGSAV=ANGLE
      RNGSAV=RANGE
      BETASV=BETA
      LENNSV=LENN
      IPARM=0

!     RETURN TO DRIVER
      RETURN
      END
