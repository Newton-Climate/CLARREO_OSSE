      SUBROUTINE GEO(IERROR,BENDNG,MSOFF,ICH1)
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     ROUTINE 'GEO' SERVES AS AN INTERFACE BETWEEN ROUTINE
!     'DRIVER' AND THE GEOMETRY SUBROUTINES INCLUDING 'GEOINP',
!     'REDUCE', 'FDBETA', 'EXPINT', 'DPEXNT', 'DPFNMN',
!     'DPFISH', 'DPSCHT', 'DPANDX', 'DPRARF', 'DPRFPA', 'DPFILL'
!     AND 'DPLAYR'.  THESE ROUTINES CALCULATE THE ABSORBER
!     AMOUNTS FOR A REFRACTED PATH THROUGH THE ATMOSPHERE.
      INTEGER IERROR,MSOFF,ICH1
      REAL BENDNG

!     PARAMETER MEXT DENOTES THE NUMBER OF MODTRAN "SPECIES".
!     THIS INCLUDES THE 12 ORIGINAL BAND MODEL PARAMETER MOLECULES
!     PLUS A HOST OF OTHER ABSORPTION AND/OR SCATTERING SOURCES.
      CHARACTER*8 CNAMEX
      COMMON/NAMEX/CNAMEX(NMOLX)
      REAL DNSTYX
      COMMON/MODELX/DNSTYX(NMOLX,LAYDIM)
      REAL DENPX,AMTPX
      COMMON/RFRPTX/DENPX(NMOLX,LAYDIM+1),AMTPX(NMOLX,LAYDIM+1)
      REAL WX
      COMMON/NONAME/WX(NMOLX)
      DOUBLE PRECISION DPH1,DPH2,DPANGL,DPPHI,DPHMIN,                   &
     &  DPBETA,DPBEND,DPRANG,SMMIN
      LOGICAL LSAVE,LNOGEO

!     SSI COMMENTS ON DOUBLE PRECISION VARIABLES:  /RFRPTH/ IS THE OLD
!     REFRACTED PATH COMMON BLOCK IN SINGLE PRECISION.  /DPRFRP/ IS
!     THE SAME COMMON BLOCK IN DOUBLE PRECISION AND IS NEW.  IN THIS
!     ROUTINE "DP" IS USED AS A PREFIX TO DENOTE THE DOUBLE PRECISION
!     VARIABLES OF DPRFRP.  THE FOLLOWING ARE THE EXCEPTIONS:
!       PPSUM    STANDS FOR THE OLD DOUBLE PRECISION SPPPSM
!       TPSUM    STANDS FOR THE OLD DOUBLE PRECISION SPTPSM
!       RHOPSM   STANDS FOR THE OLD DOUBLE PRECISION SPRHOP

!     SOME OTHER VARIABLES WERE DECLARED DOUBLE PRECISION; THEY ARE ALSO
!     IDENTIFIED BY THE "DP" PREFIX, BUT ARE NOT COMMON BLOCK VARIABLES.
      INCLUDE 'BASE.h'
      INCLUDE 'IFIL.h'

!     /CARD1/
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,MODTRN

!     /SURFWV/
!       LAMBER  LOGICAL FLAG, .TRUE. FOR LAMBERTIAN SURFACE.
!       TPTEMP  TARGET-PIXEL SURFACE TEMPERATURES [K].
!       TPHDIR  TARGET-PIXEL HEMISPHERE DIRECTIONAL REFLECTANCE AT
!               VIEWING ANGLE.
!       TPBRDF  TARGET-PIXEL BIDIRECTIONAL REFLECTANCE DISTRIBUTION
!               FUNCTION AT VIEWING AND SUN ANGLE.
!       AATEMP  AREA-AVERAGED GROUND SURFACE TEMPERATURES [K].
!       AASALB  AREA-AVERAGED GROUND SURFACE ALBEDO.
!       AADREF  AREA-AVERAGED GROUND SURFACE DIRECTIONAL REFLECTIVITY
!               AT THE SOLAR ZENITH ANGLE.
!       EMU     GROUND DIRECTIONAL EMISSIVITY AT VIEWING ANGLE.
!       BEM     GROUND DIRECTIONAL EMISSIVITY AT QUADRATURE ANGLE.
!       RMU     GROUND BRDF AZIMUTH COMPONENTS AT VIEWING ANGLE
!               AND AT SUN (=0) OR QUADRATURE (>0) ANGLE.
!       BDR     GROUND BRDF AZIMUTH COMPONENTS AT QUADRATURE ANGLE
!               AND AT SUN (=0) OR QUADRATURE (>0) ANGLE.
      LOGICAL LAMBER
      REAL TPTEMP,TPHDIR,TPBRDF,AATEMP,AASALB,AADREF,EMU,BEM,RMU,BDR
      COMMON/SURFWV/LAMBER,TPTEMP,TPHDIR,TPBRDF,AATEMP,AASALB,AADREF,   &
     &  EMU(MXUMU),BEM(MI),RMU(1:MXUMU,0:MI,0:MAZ),BDR(1:MI,0:MI,0:MAZ)

!     /ANGSRF/
!       CVWSRF  COSINE OF THEN VIEW ZENITH ANGLE FROM THE SURFACE [RAD].
!       CSNSRF  COSINE OF THE SOLAR (LUNAR) ZENITH AT SURFACE [RAD].
!       AZMSRF  RELATIVE AZIMUTH ANGLE (SUN - SENSOR AT SURFACE) [RAD].
!       UMU1    COSINE OF THE PATH NADIR ANGLE.
!               (AT H1 IF IMULT=1; AT OR "NEAR" H2 IF IMULT=-1)
!       UMU0    COSINE OF THE SOLAR ZENITH ANGLE.
!               (AT H1 IF IMULT=1; AT OR "NEAR" H2 IF IMULT=-1)
!       PHI1    RELATIVE AZIMUTH ANGLE (SUN - LOS PATH AT SENSOR) [DEG].
!               (AT H1 IF IMULT=1; AT OR "NEAR" H2 IF IMULT=-1)
!       CMU     COSINE OF THE NADIR ANGLES USED IN DISORT.
      REAL CVWSRF,CSNSRF,AZMSRF,UMU1,UMU0,PHI1,CMU
      COMMON/ANGSRF/CVWSRF,CSNSRF,AZMSRF,UMU1,UMU0,PHI1,CMU(MI)

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

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT

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
      REAL RE,ZMAX
      INTEGER IPATH
      COMMON/PARMTR/RE,ZMAX,IPATH

!     /RFRPTH/
      REAL SPZP,SPPP,SPTP,SPPPSM,SPTPSM,SPRHOP,SPDENP,SPAMTP,DRANGE
      COMMON/RFRPTH/SPZP(LAYDIM+1),SPPP(LAYDIM+1),SPTP(LAYDIM+1),       &
     &  SPPPSM(LAYDIM+1),SPTPSM(LAYDIM+1),SPRHOP(LAYDIM+1),             &
     &  SPDENP(MEXT,LAYDIM+1),SPAMTP(MEXT,LAYDIM+1),DRANGE(LAYDIM+1)
      DOUBLE PRECISION ZP,PP,TP,RFNDXP,PPSUM,TPSUM,RHOPSM,DENP,AMTP
      COMMON/DPRFRP/ZP(LAYDIM+1),PP(LAYDIM+1),TP(LAYDIM+1),             &
     &  RFNDXP(LAYDIM+1),PPSUM(LAYDIM+1),TPSUM(LAYDIM+1),               &
     &  RHOPSM(LAYDIM+1),DENP(MEXT,LAYDIM+1),AMTP(MEXT,LAYDIM+1)
      INCLUDE 'SOLS.h'

!     /PATH/
!       QTHETA  COSINE OF PATH ZENITH AT PATH BOUNDARIES.
!       AHT     ALTITUDES AT PATH BOUNDARIES [KM].
!       IHT     ALTITUDES AT PATH BOUNDARIES [M].
!       TPH     TEMPERATURE AT PATH BOUNDARIES [K].
!       IMAP    MAPPING FROM PATH SEGMENT MIDPOINT TO VERTICAL LAYER.
!       LOWAHT  INDEX OF VERTICAL LAYER BOUNDARY AT OR JUST BELOW AHT.
!       FACAHT  ALTITUDE INTERPOLATION FRACTION FOR AHT.
      INTEGER IHT,IMAP,LOWAHT
      REAL QTHETA,AHT,TPH,FACAHT
      COMMON/PATH/QTHETA(LAYTWO),AHT(LAYTWO),IHT(0:LAYTWO),             &
     &  TPH(LAYTWO),IMAP(LAYTWO),LOWAHT(LAYTWO),FACAHT(LAYTWO)
      DOUBLE PRECISION DHALFR,DPRNG2
      COMMON/SMALL1/DHALFR,DPRNG2
      LOGICAL LSMALL
      COMMON/SMALL2/LSMALL
      REAL SMALL
      COMMON/SMALL3/SMALL
      LOGICAL LPRINT
      COMMON/CPRINT/LPRINT
      REAL GNDALT
      COMMON/GRAUND/GNDALT

!     /AER/
!     THERE ARE "MAER=17" COMPONENTS:
!      1     AEROSOL 1 (NOMINALLY, BOUNDARY LAYER AEROSOL).
!      2     AEROSOL 2 (NOMINALLY, TROPOSPHERIC AEROSOL).
!      3     AEROSOL 3 (NOMINALLY, STRATOSPHERIC AEROSOL).
!      4     AEROSOL 4 (NOMINALLY, VOLCANIC AEROSOL).
!      5     CIRRUS CLOUD.
!      6     CLOUD 1 (NOMINALLY, WATER CLOUD).
!      7     CLOUD 2 (NOMINALLY, ICE CLOUD).
!     8-17   NOVAM (NAVY OCEANIC VERTICAL AEROSOL MODEL) AEROSOL LAYERS.

!       NAER     NUMBER OF ACTIVE AEROSOLS.
!       EXTV     SPECTRAL EXTINCTION (NORMALIZED TO 1 AT 550 NM).
!       ABSV     SPECTRAL ABSORPTION (1-ABSV/EXTV=SCATTERING ALBEDO).
!       ASYV     SPECTRAL LEGENDRE MOMENT (DIVIDED BY 2N+1).
      INTEGER NAER
      REAL EXTV,ABSV,ASYV
      COMMON/AER/NAER,EXTV(MAER),ABSV(MAER),ASYV(MXCMU,MAER)

!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /SMALL3/,/NAMEX/
      EXTERNAL DEVCBD,ATMCON,XMLATM

!     DECLARE FUNCTIONS
      INTEGER PSLCT
      REAL EXPINT

!     DECLARE LOCAL VARIABLES
      INTEGER JMAX,JMAXP1,I,J,K,L,ILO,IHI,IH1,IH2,JOFF,JOFFM1,ISLCT,    &
     &  LENNSV,IAMT,LDEL
      LOGICAL LOGLOS
      REAL PHI,FAC,WTEM,RANGSV,HMIN,PHISV,PRCNT,H2MX,H2SAV

!     DATA:
!       KMOL     A POINTER USED TO REORDER THE AMOUNTS WHEN PRINTING.
!       RAY550   RAYLEIGH SCATTERING COEFFICIENT (KM-1) AT 550 NM
!                (SEE ROUTINE C6DTA).
      INTEGER KMOL(17)
      REAL TLRNCE,PZERO,RAY550
      DATA KMOL/1,2,3,11,8,5,9,10,4,6,7,12,13,14,16,15,17/
      DATA TLRNCE/0.001/,PZERO/1013.25/,RAY550/.0121181/
      H2SAV=H2
      LSMALL=.FALSE.
      LPRINT=.TRUE.
      LNOGEO=.FALSE.

!     INITIALIZE CONSTANTS AND CLEAR CUMULATIVE VARIABLES.
      IF(.NOT.LJMASS .AND. GNDALT.GT.ZM(1) .AND. ISSGEO.EQ.0)           &
     &  WRITE(IPR,'(/A,2F10.6,I5)')                                     &
     &  ' GNDALT IS ABOVE FIRST PROFILE ALTITUDE:',GNDALT,ZM(1),MODEL
      IERROR=0

!     INITIALIZE CUMULATIVE VARIABLES
      DO 30 I=1,LAYDIM+1
          LJ(I)=0
          SPPPSM(I)=0.
          SPTPSM(I)=0.
          SPRHOP(I)=0.
          DRANGE(I)=0.
          PPSUM(I)=DBLE(0.)
          TPSUM(I)=DBLE(0.)
          RHOPSM(I)=DBLE(0.)
          DO 10 K=1,MEXT
              SPAMTP(K,I)=0.
              AMTP(K,I)=DBLE(0.)
   10     CONTINUE
          DO 20 K=1,NMOLX
              AMTPX(K,I)=0.
   20     CONTINUE
   30 CONTINUE
      IF(ITYPE.LT.2)THEN

!         HORIZONTAL PATH, MODEL EQ 1 TO 7:  INTERPOLATE PROFILE TO H1
          SPZP(1)=H1
          ZP(1)=DBLE(SPZP(1))
          IF(ML.EQ.1)THEN
              SPPP(1)=PM(1)
              PP(1)=DBLE(PM(1))
              SPTP(1)=TM(1)
              TPH(1)=SPTP(1)
              TP(1)=DBLE(TM(1))
              LJ(1)=1
              DO 40 K=1,MEXT
                  SPDENP(K,1)=DENSTY(K,1)
   40         CONTINUE
              DO 50 K=1,NMOLX
                  DENPX(K,1)=DNSTYX(K,1)
   50         CONTINUE
          ELSEIF(H1.LT.ZM(1))THEN
              WRITE(IPR,'(A,2F10.6)')                                   &
     &             ' ERROR HORIZ PATH H1 < ZM(1)',H1,ZM(1)
              IERROR=1
              RETURN
          ELSE
              ILO=1
              DO 60 IHI=2,ML
                  IF(H1.LT.ZM(IHI))GOTO70
   60         ILO=IHI
   70         FAC=(H1-ZM(ILO))/(ZM(IHI)-ZM(ILO))
              SPPP(1)=EXPINT(PM(ILO),PM(IHI),FAC)
              SPTP(1)=TM(ILO)+(TM(IHI)-TM(ILO))*FAC
              TPH(1)=SPTP(1)
              TP(1)=DBLE(TM(ILO)+(TM(IHI)-TM(ILO))*FAC)
              LJ(1)=ILO
              IF(FAC.GT..5)LJ(1)=IHI
              DO 80 K=1,MEXT
                  SPDENP(K,1)=EXPINT(DENSTY(K,ILO),DENSTY(K,IHI),FAC)
   80         CONTINUE
              DO 90 K=1,NMOLX
                  DENPX(K,1)=EXPINT(DNSTYX(K,ILO),DNSTYX(K,IHI),FAC)
   90         CONTINUE

!             USE LINEAR INTERPOLATION FOR CLOUDS
              SPDENP(16,1)=DENSTY(16,ILO)                               &
     &          +FAC*(DENSTY(16,IHI)-DENSTY(16,ILO))
              SPDENP(66,1)=DENSTY(66,ILO)                               &
     &          +FAC*(DENSTY(66,IHI)-DENSTY(66,ILO))
              SPDENP(67,1)=DENSTY(67,ILO)                               &
     &          +FAC*(DENSTY(67,IHI)-DENSTY(67,ILO))
          ENDIF

!         CALCULATE ABSORBER AMOUNTS FOR A HORIZONTAL PATH
          IF(.NOT.LJMASS) WRITE(IPR,'(/2(A,F11.3),A,I4)')               &
     &      ' HORIZONTAL PATH AT ALTITUDE =',H1,                        &
     &      ' KM WITH RANGE =',RANGE,' KM, MODEL =',MODEL
          IKMAX=1
          JMAX=1
          IF(MODEL.EQ.0)THEN
              SPTP(1)=TM(1)
              SPPP(1)=PM(1)
              TP(1)=DBLE(TM(1))
              PP(1)=DBLE(PM(1))
          ENDIF
          JOFF=MSOFF+1
          PATM(JOFF)=SPPP(1)/PZERO
          TBBY(JOFF)=SPTP(1)
          DRNG(JOFF)=RANGE
          DO 100 K=1,MEXT
              W(K)=SPDENP(K,1)*RANGE
              WPATH(JOFF,K)=W(K)
  100     CONTINUE
          DO 110 K=1,NMOLX
              WX(K)=DENPX(K,1)*RANGE
              WPATH(JOFF,MEXT+K)=WX(K)
  110     CONTINUE
          W(9)=W(5)*(296.-SPTP(1))/(296.-260.)
          WTEM=SPTP(1)-273.15
          W(59)=W(8)*.269*WTEM
          W(60)=W(59)*WTEM
          WPATH(JOFF,9)=W(9)
          WPATH(JOFF,59)=W(59)
          WPATH(JOFF,60)=W(60)
      ELSE

!         SLANT PATH SELECTED.  INTERPRET SLANT PATH PARAMETERS.

!         LOGLOS IS A LOGICAL VARIABLE THAT IS TRUE
!         ONLY THE OPTICAL PATH WITH ITYPE = 2.
          LOGLOS=.FALSE.
          IF(ITYPE.EQ.2 .AND. MSOFF.EQ.0)THEN
              LOGLOS=.TRUE.
              ISLCT=PSLCT(ANGLE,RANGE,BETA)

!             SPECIAL TREATMENT EXCEPT FOR CASE 2A (ISLCT=21)
              IF(ISLCT.GT.21)THEN

!                 IF RANGE IS SMALL, CONVERT TO CASE 2C (ISLCT=23)
                  CALL SMPREP(H1,H2,ANGLE,RANGE,BETA,ISLCT)
                  IF(RANGE.GT.0. .AND. RANGE.LE.SMALL)THEN
                      LSMALL=.TRUE.
                      RANGSV=RANGE
                      ISLCT= 23
                  ELSEIF(ISLCT.EQ.22)THEN

!                     CASE 2B:  (H1,ANGLE,RANGE)
!                     IF PATH TYPE IS CASE 2B, CHECK THAT THE RANGE
!                     USED IN THE CALCULATION EQUALS THE INPUT RANGE.
!                     DETERMINE H2 AND LENN USING ROUTINE NEWH2.
                      LENNSV=LENN
                      RANGSV=RANGE
                      DPH1=DBLE(H1)
                      DPANGL=DBLE(ANGLE)
                      DPRANG=DBLE(RANGE)
                      CALL NEWH2(DPH1,DPH2,DPANGL,                      &
     &                  DPRANG,DPBETA,LENN,DPHMIN,DPPHI)
                      IF(LENN.EQ.0)DPHMIN=MIN(DPH2,DPH1)
                      H2=REAL(DPH2)
                      HMIN=REAL(DPHMIN)
                      PHISV=REAL(DPPHI)
                      LPRINT=.FALSE.
                      IAMT=2
                      CALL DPRFPA(DPH1,DPH2,DPANGL,DPPHI,               &
     &                  LENN,DPHMIN,IAMT,DPBETA,DPRANG,DPBEND)
                      LPRINT=.TRUE.
                      PRCNT=100*ABS(REAL(DPRANG)-RANGSV)/RANGSV
                      IF(.NOT.LJMASS) THEN
                         WRITE(IPR,'((A))')'1',                         &
     &                     ' SOME INTERNAL DETAILS:'
                         WRITE(IPR,'(/(A))')                            &
     &                     ' LOS IS INTERNAL CASE 2B',                  &
     &                     ' (H1, ANGLE, RANGE).',                      &
     &                     ' USING H2 OBTAINED FROM SUBROUTINE NEWH2:'
                         WRITE(IPR,'(/(A,F12.5))')                      &
     &                     ' H1                        =',H1,           &
     &                     ' H2                        =',H2,           &
     &                     ' ANGLE                     =',ANGLE,        &
     &                     ' PHI                       =',PHISV,        &
     &                     ' BETA                      =',BETA,         &
     &                     ' HMIN (MINIMUM ALTITUDE)   =',HMIN,         &
     &                     ' RANGE (OUTPUT)            =',DPRANG,       &
     &                     ' RANGE (INPUT)             =',RANGSV,       &
     &                     ' PERCENT DIFFERENCE        =',PRCNT
                         WRITE(IPR,'(A,I12,/)')                         &
     &                     ' LENN                      =',LENN
                      ENDIF
                      IF(ANGLE.GT.90. .AND.                             &
     &                  ABS(H2-GNDALT).LT.TLRNCE)THEN
                          LNOGEO=.TRUE.
                          PHI=PHISV
                          ANGLE=REAL(DPANGL)
                          HMIN=REAL(DPHMIN)
                          BETA=REAL(DPBETA)
                          RANGE=REAL(DPRANG)
                          WRITE(IPR,'(/2A,/)')'***** WARNING ***** ',   &
     &                      ' PATH HITS THE EARTH.'
                      ELSEIF(ANGLE.LE.90. .AND.                         &
     &                  ABS(H2-ZMAX).LT.TLRNCE)THEN
                          LNOGEO=.TRUE.
                          PHI=PHISV
                          ANGLE=REAL(DPANGL)
                          HMIN=REAL(DPHMIN)
                          BETA=REAL(DPBETA)
                          RANGE=REAL(DPRANG)
                          WRITE(IPR,'(/2A,/)')'***** WARNING ***** ',   &
     &                      ' PATH HITS THE UPPERMOST LAYER BOUNDARY.'
                      ELSEIF(PRCNT.LE.1. .OR. H2.GE.ZMAX)THEN
                          IF(.NOT.LJMASS)                               &
     &                       WRITE(IPR,'(/(2A))')                       &
     &                         ' PERCENT DIFFERENCE IS THAN 1 OR THE',  &
     &                         ' PATH TERMINATES AT THE TOP OF THE',    &
     &                         ' ATMOSPHERE; THESE PATH PARAMETERS',    &
     &                         ' WILL BE USED WITHOUT CALLING GEOINP.'
                          LNOGEO=.TRUE.
                          PHI=PHISV
                          ANGLE=REAL(DPANGL)
                          HMIN=REAL(DPHMIN)
                          BETA=REAL(DPBETA)
                          RANGE=REAL(DPRANG)
                      ELSE
                          IF(.NOT.LJMASS)                               &
     &                       WRITE(IPR,'(/2A,/2A,/)')                   &
     &                         ' SINCE THE PERCENT',                    &
     &                         ' DIFFERENCE EXCEEDS 1,',' "EQUIVALENT"',&
     &                         ' INTERNAL CASE 2C WILL BE USED.'
                          LNOGEO=.FALSE.
                          RANGE=RANGSV
                          BETA=0.
                          HMIN=0.
                          PHI=0.
                          LENN=LENNSV
                          ANGLE=0.
                          ISLCT=PSLCT(ANGLE,RANGE,BETA)
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF
          LSAVE=LSMALL
          LSMALL=.FALSE.
          IF(.NOT.LNOGEO)CALL GEOINP(H1,H2,ANGLE,RANGE,                 &
     &      BETA,ITYPE,LENN,HMIN,PHI,IERROR,ISLCT,LOGLOS,H2MX)
          LSMALL=LSAVE

!         CHECK FOR IERROR
          IF(IERROR.NE.0)THEN
              IF(.NOT.LJMASS .AND. ISSGEO.EQ.0)                         &
     &          WRITE(IPR,'(2A)')'0GEO:  IERROR NE 0:',                 &
     &          ' END THIS CALCULATION AND SKIP TO THE NEXT CASE'
              RETURN
          ELSEIF(LSMALL)THEN
              RANGE=RANGSV
              CALL SMGEO(ANGLE,BETA,PHI,DHALFR,DPRNG2,BENDNG,LENN,SMMIN)
              HMIN=REAL(SMMIN)
          ENDIF

!         CALCULATE THE PATH THROUGH THE ATMOSPHERE
          IAMT=1
          DPH1=DBLE(H1)
          DPH2=DBLE(H2)
          DPANGL=DBLE(ANGLE)
          DPPHI=DBLE(PHI)
          DPHMIN=DBLE(HMIN)
          DPBETA=DBLE(BETA)
          DPRANG=DBLE(RANGE)

!         IF ANGLE IS NEAR 90., HMIN IS NEAR H1, AND RANGE EXCEEDS
!         SMALL, THEN YOU HAVE A HALF-TANGENT PATH I.E. H1, H2 AND
!         THE EARTH CENTER FORM A RIGHT TRIANGLE.  IN THIS CASE, DO
!         NOT LET LENN EQUAL 1 BECAUSE THAT MAY MAKE HMIN IDENTICALLY
!         EQUAL TO H1 PRODUCING PROBLEMS IN SUBROUTINE FILL.  THE
!         CHECK ON THE VARIABLE SMALL SIMPLY POINTS OUT THE FACT
!         SMALL PATHS ARE ALREADY DEALT WITH BY OTHER METHODS.
          IF(ABS(ANGLE-90.).LT..001 .AND. ABS(HMIN-H1).LT.TLRNCE        &
     &      .AND. RANGE.GT.SMALL)LENN=0
          CALL DPRFPA(DPH1,DPH2,DPANGL,DPPHI,                           &
     &      LENN,DPHMIN,IAMT,DPBETA,DPRANG,DPBEND)
          H1=REAL(DPH1)
          H2=REAL(DPH2)
          ANGLE=REAL(DPANGL)
          PHI=REAL(DPPHI)
          HMIN=REAL(DPHMIN)
          IF(HMIN.LT.GNDALT-TLRNCE .AND. (LOGLOS .OR. ITYPE.EQ.3))THEN
!             LOGLOS ONLY DEALS WITH ITYPE=2.
!             THEREFORE, THE CHECK FOR ITYPE=3.
              IF(.NOT.LJMASS)                                           &
     &          WRITE(IPR,'(2(/A,F12.6),/(A))')'GNDALT =',GNDALT,       &
     &          'HMIN =',HMIN,'HMIN IS LESS THAN GNDALT.',              &
     &          'THIS RUN ABORTED, NEXT RUN ATTEMPTED',' '
              IERROR=1
              RETURN
          ENDIF
          BETA=REAL(DPBETA)
          RANGE=REAL(DPRANG)
          BENDNG=REAL(DPBEND)

!         LOAD LAYER AMOUNTS IN SPAMTP INTO WPATH FROM H1 TO H2
          DO 120 I=1,IPATH
              IF(H1.EQ.SPZP(I))IH1=I
              IF(H2.EQ.SPZP(I))IH2=I
  120     CONTINUE
          JMAX=(IPATH-1)+LENN*(MIN0(IH1,IH2)-1)
          IF(JMAX.GT.LAYTWO) THEN
              WRITE(IPR,'(/A) ')                                        &
     &              'JMAX IS GREATER THAN PARAMETER LAYTWO'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP 'JMAX IS GREATER THAN PARAMETER LAYTWO'
          ENDIF
          IKMAX=JMAX

!         DETERMINE LJ(J), WHICH IS THE LAYER NUMBER L=LJ(J)
!         IN SPAMTP(K,L), STARTING FROM HMIN, WHICH CORRESPONDS
!         TO THE LAYER J IN WPATH(J+MSOFF,K), STARTING FROM
!         H1 INITIAL DIRECTION OF PATH IS DOWN.
          L=IH1
          LDEL=-1
          IF(LENN.EQ.0 .AND. H1.LE.H2)THEN

!             INITIAL DIRECTION OF PATH IS UP
              L=0
              LDEL=1
          ENDIF
          JTURN=0
          JMAXP1=JMAX+1
          DO 130 J=1,JMAXP1

!             TEST FOR REVERSING DIRECTION OF PATH FROM DOWN TO UP
              IF(L.EQ.1 .AND. LDEL.EQ.-1)THEN
                  JTURN=J
                  L=0
                  LDEL=1
              ENDIF
              L=L+LDEL
              LJ(J)=L
  130     CONTINUE

!         LOAD TBBY = THE DENSITY WEIGHTED MEAN TEMPERATURE
!         AND WPATH = THE INCREMENTAL LAYER AMOUNTS.
          DO 140 K=1,MEXT
              W(K)=0.
  140     CONTINUE
          DO 150 K=1,NMOLX
              WX(K)=0.
  150     CONTINUE
          DO 180 J=1,JMAX
              L=LJ(J)
              IF(L.GE.ML)L=ML
              JOFF=MSOFF+J
              TBBY(JOFF)=SPTPSM(L)/SPRHOP(L)
              PATM(JOFF)=SPPPSM(L)/(PZERO*SPRHOP(L))
              DRNG(JOFF)=DRANGE(L)
              SPAMTP(9,L)=SPAMTP(5,L)*(296.-TBBY(JOFF))/36.
              WTEM=TBBY(JOFF)-273.15
              SPAMTP(59,L)=.269*SPAMTP(8,L)*WTEM
              SPAMTP(60,L)=SPAMTP(59,L)*WTEM
              DO 160 K=1,MEXT
                  WPATH(JOFF,K)=SPAMTP(K,L)
                  W(K)=W(K)+SPAMTP(K,L)
  160         CONTINUE
              DO 170 K=1,NMOLX
                  WPATH(JOFF,MEXT+K)=AMTPX(K,L)
                  WX(K)=WX(K)+AMTPX(K,L)
  170         CONTINUE
  180     CONTINUE
          IKMAX=JMAX

!         INCLUDE BOUNDARY EMISSION IF CARD 1 CONTAINS A NON-ZERO
!         TPTEMP OR IF THE SLANT PATH INTERSECTS THE EARTH [TPTEMP
!         SET TO TEMPERATURE OF LOWEST BOUNDARY IN THE LATTER CASE].
          IF(TPTEMP.EQ.0. .AND. H2.LE.ZM(1))TPTEMP=TM(1)

!         PRINT LAYER PATH SEGMENT ABSORBER AMOUNTS
          IF(NPR.NE.2)THEN

!             LINEARLY INTERPOLATE TEMPERATURE AT H1:
              ILO=1
              DO 190 IHI=2,ML-1
                  IF(ZM(IHI).GE.H1)GOTO200
                  ILO=IHI
  190         CONTINUE
              IHI=ML
  200         CONTINUE
              TPH(1)=TM(ILO)                                            &
     &          +(TM(IHI)-TM(ILO))*(H1-ZM(ILO))/(ZM(IHI)-ZM(ILO))

!             LDEL=0 => GOING DOWN;    LDEL=1 => GOING UP.
              LDEL=1
              IF(LENN.EQ.1 .OR. H1.GT.H2)LDEL=0
              AHT(1)=H1
              IHT(0)=INT(1000*H1+.5)
              DO 210 J=1,JMAX
                  L=LJ(J)
                  IF(J.EQ.JTURN)LDEL=1
                  AHT(J+1)=SPZP(L+LDEL)
                  IHT(J)=INT(1000*SPZP(L+LDEL)+.5)
                  TPH(J+1)=SPTP(L+LDEL)
  210         CONTINUE

!             OUTPUTS AVAILABLE IN COMMON BLOCKS
              IF(.NOT.LJMASS .AND. NPR.NE.1)THEN
                  WRITE(IPR,'(////2A,//(3A))')'1  LAYER ABSORBER',      &
     &              ' AMOUNTS FOR THE PATH SEGMENT ENDING AT Z',        &
     &              '  J      Z        TBAR       HNO3         O3 UV ', &
     &              '    CNTMSLF1     CNTMSLF2      CNTMFRN ',          &
     &              '       O2         WAT DROP     ICE PART',          &
     &              '       (KM)        (K)     (ATM CM)     (ATM CM)', &
     &              '   (MOL CM-2)   (MOL CM-2)   (MOL CM-2)',          &
     &              '   (MOL CM-2)   (KM GM/M3)   (KM GM/M3)'
                  WRITE(IPR,'(/(I3,0P,F10.4,F9.2,1P,6E13.3,0P,2F13.7))')&
     &              (J,AHT(J+1),TBBY(MSOFF+J),                          &
     &              (WPATH(MSOFF+J,KMOL(K)),K=4,8),WPATH(MSOFF+J,58),   &
     &              (WPATH(MSOFF+J,K),K=66,67),J=1,JMAX)
                  WRITE(IPR,'(///A,T10,A,T19,A,T31,A,T47,A,T60,A,T73,   &
     &              A,T86,A,T99,A,T110,A,T123,A,/T9,A,T29,A,T109,A)')   &
     &              '1 J','Z','N2 CONT','MOL SCAT','AER 1','AER 2',     &
     &              'AER 3','AER 4','CIRRUS','WAT DROP','ICE PART',     &
     &              '(KM)','(550nm EXT)','(550 NM OPTICAL DEPTH)'

!                 FETCH CLOUD 550 NM EXTINCTION COEFFICIENTS:
                  CALL AEREXT(18181.82,6)
                  WRITE(IPR,'(/(I3,0P,F10.4,9F13.7))')(J,AHT(J+1),      &
     &              WPATH(MSOFF+J,4),RAY550*WPATH(MSOFF+J,6),           &
     &              (WPATH(MSOFF+J,KMOL(K)),K=11,15),WPATH(MSOFF+J,66)  &
     &              *EXTV(6),WPATH(MSOFF+J,67)*EXTV(7),J=1,JMAX)
              ENDIF
          ENDIF

!         PRINT PATH SUMMARY
          IF(ISSGEO.EQ.0)THEN
              IF(.NOT.LJMASS) THEN
                  IF(NPR.LT.1)THEN
                     WRITE(IPR,'(///3A)')                               &
     &                 '1    J    Z       H2O       O3        CO2 ',    &
     &                 '      CO        CH4       N2O       O2  ',      &
     &                 '      NH3       NO        NO2       SO2'
                     IF(MODTRN)THEN
                        WRITE(IPR,'(8X,A,44X,A,45X,A)')                 &
     &                     '(KM)    (          ','ATM CM',')'
                     ELSE
                        WRITE(IPR,'(8X,A,44X,A,45X,A)')                 &
     &                     '(KM)   (G/CM**2)  (','ATM CM',')'
                     ENDIF
                     WRITE(IPR,'(/(I4,0P,F10.4,1P,11E10.2))')           &
     &                       (J,AHT(J+1),WPATH(MSOFF+J,17),             &
     &                 WPATH(MSOFF+J,31),WPATH(MSOFF+J,36),             &
     &                 WPATH(MSOFF+J,44),WPATH(MSOFF+J,46),             &
     &                 WPATH(MSOFF+J,47),WPATH(MSOFF+J,50),             &
     &                 WPATH(MSOFF+J,52),WPATH(MSOFF+J,54),             &
     &                 WPATH(MSOFF+J,55),WPATH(MSOFF+J,56),J=1,JMAX)
                     WRITE(IPR,'(///A14,13(1X,A8:),/(14X,13(1X,A8)))')  &
     &                 '1  J      Z   ',(CNAMEX(K),K=1,NMOLX)
                     WRITE(IPR,'(8X,A,2(53X,A))')                       &
     &                 '(KM)    (','ATM CM',')'
                     DO 220 J=1,JMAX
                         JOFF=MSOFF+J
                         WRITE(IPR,'(I4,F10.4,1P,13E9.2:,               &
     &                     /(14X,1P,13E9.2:))')                         &
     &                     J,AHT(J+1),(WPATH(JOFF,MEXT+K),K=1,NMOLX)
 220                 CONTINUE
                  ENDIF
                  WRITE(IPR,'(//A,/8(/10X,A,F11.5,A),/10X,A,I11)')      &
     &              ' SUMMARY OF THE GEOMETRY CALCULATION',             &
     &              'H1      =',H1,   ' KM', 'H2      =',H2,    ' KM',  &
     &              'ANGLE   =',ANGLE,' DEG','RANGE   =',RANGE, ' KM',  &
     &              'BETA    =',BETA, ' DEG','PHI     =',PHI,   ' DEG', &
     &              'HMIN    =',HMIN, ' KM', 'BENDING =',BENDNG,' DEG', &
     &              'LENN    =',LENN
              ENDIF

!             SAVE ZENITH ANGLE FROM H2 TO H1 AND INITIALIZE THE COSINE
!             SOLAR ZENITH AND RELATIVE AZIMUTH ANGLES AT H2.
              CVWSRF=COS(PHI/DEG)
              CSNSRF=0.
              AZMSRF=PI
              IF(IMULT.EQ.1)THEN
                  UMU1=-COS(ANGLE/DEG)
              ELSE
                  UMU1=CVWSRF
              ENDIF

!             SINCE UMU1 IS USED IN THE PLANE-PARALLEL MULTIPLE
!             SCATTERING DISORT MODEL, IT CANNOT BE TOO CLOSE TO 0.
              IF(ABS(UMU1).LT..05 .OR.                                  &
     &          (ITYPE.EQ.3 .AND. H2SAV.GT.0.))THEN
                  IF(UMU1.GE.0.)THEN
                      UMU1=.05
                  ELSE
                      UMU1=-.05
                  ENDIF
              ENDIF
              UMU0=0.
          ENDIF
      ENDIF

!     CALCULATE THE AEROSOL WEIGHTED MEAN RELATIVE HUMIDITY (RH)
      IF(W(7).GT.0. .AND. ICH1.LE.7)THEN
          W(15)=100.-EXP(W(15)/W(7))
      ELSEIF(W(12).GT.0. .AND. ICH1.GT.7)THEN
          W(15)=100.-EXP(W(15)/W(12))
      ELSE
          W(15)=0.
      ENDIF

!     CONVERT MOLECULAR BAND WPATH TO CUMULATIVE
!     PATH AMOUNTS FOR LOWTRAN RUNS.
      IF(.NOT.MODTRN)THEN
          JOFFM1=MSOFF+1
          DO 240 J=2,JMAX
              JOFF=MSOFF+J
              DO 230 K=17,57
                  WPATH(JOFF,K)=WPATH(JOFFM1,K)+WPATH(JOFF,K)
  230         CONTINUE
  240     JOFFM1=JOFF
      ENDIF

!     PRINT TOTAL PATH AMOUNTS
      IF(H1.LT.ZM(1))THEN
          WRITE(IPR,'(/2A,F8.4,A,/14X,A,F8.4,A)')' FATAL ERROR: ',      &
     &         ' INITIAL ALTITUDE [H1 =',H1,' KM] IS BELOW THE',        &
     &         ' BOTTOM OF THE ATMOSPHERE [ZM(1) =',ZM(1),' KM].'
          IERROR=1
          RETURN
      ELSEIF(H2.LT.ZM(1) .AND. ITYPE.NE.1)THEN
          WRITE(IPR,'(/2A,F8.4,A,/14X,A,F8.4,A)')' FATAL ERROR: ',      &
     &         ' FINAL OR TANGENT ALTITUDE [H2 =',H2,' KM] IS BELOW',   &
     &         ' THE BOTTOM OF THE ATMOSPHERE [ZM(1) =',ZM(1),' KM].'
          IERROR=1
          RETURN
      ENDIF
      IF(ISSGEO.EQ.1)RETURN
      IF(.NOT.LJMASS) THEN
          IF(.NOT.MODTRN)THEN
             WRITE(IPR,'(//2A)')'   EQUIVALENT SEA',                    &
     &         ' LEVEL TOTAL ABSORBER AMOUNTS:'
          ELSEIF(MSOFF.EQ.0.)THEN
             WRITE(IPR,'(//2A)')'   TOTAL COLUMN ABSORBER',             &
     &         ' AMOUNTS FOR THE LINE-OF-SIGHT PATH:'
          ELSE
             WRITE(IPR,'(//2A)')'   TOTAL COLUMN ABSORBER',             &
     &         ' AMOUNTS FOR A VERTICAL PATH FROM GROUND TO SPACE:'
!lar         write(iplot,*)w(17)
          ENDIF
          WRITE(IPR,'(2(/15X,2A),/10X,1P,7E12.4,                        &
     &      /(//15X,2A:,/73X,A,/10X,0P,7F12.6,F12.2))')                 &
     &      '  HNO3      O3 UV      CNTMSLF1    CNTMSLF2    CNTMFRN  ', &
     &      '   N2 CONT     MOL SCAT',                                  &
     &      '(ATM CM)   (ATM CM)   (MOL CM-2)  (MOL CM-2)  (MOL CM-2)', &
     &      '             (550 NM EXT)',(W(KMOL(K)),K=4,9),RAY550*W(6), &
     &      ' AER 1       AER 2       AER 3      ',                     &
     &      ' AER 4      CIRRUS     WAT DROP    ICE PART   MEAN AER RH',&
     &      '(KM GM/M3)  (KM GM/M3)    (PRCNT)',                        &
     &      (W(KMOL(K)),K=11,15),W(66),W(67),W(KMOL(16)),' H2O     ',   &
     &      '    O3          CO2         CO          CH4         N2O'
            IF(MODTRN)THEN
               WRITE(IPR,'(13X,A,2(35X,A))')'(','ATM CM',')'
            ELSE
               WRITE(IPR,'(13X,A,2(22X,A))')                            &
     &           '(G/CM**2)     (','ATM CM',')'
            ENDIF
            WRITE(IPR,'((10X,1P,6E12.4,//2(/15X,A),2(22X,A)))')         &
     &        W(17),W(31),W(36),W(44),W(46),W(47),                      &
     &        ' O2          NH3         NO          NO2         SO2',   &
     &        '(','ATM CM',')',W(50),W(52),W(54),W(55),W(56)
          IHI=7
  250     CONTINUE
          IF(NMOLX.GT.IHI)THEN
              WRITE(IPR,'(//8X,7(4X,A),/13X,A,2(36X,A),/10X,1P,7E12.4)')&
     &          (CNAMEX(I),I=IHI-6,IHI),                                &
     &          '(','ATM CM',')',(WX(I),I=IHI-6,IHI)
              IHI=IHI+7
              GOTO250
          ENDIF
          WRITE(IPR,'(//8X,7(4X,A))')(CNAMEX(I),I=IHI-6,NMOLX)
          WRITE(IPR,'(13X,A,2(36X,A),/10X,1P,7E12.4)')                  &
     &      '(','ATM CM',')',(WX(I),I=IHI-6,NMOLX)
      ENDIF
      RETURN
      END
