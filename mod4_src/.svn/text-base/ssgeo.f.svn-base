      SUBROUTINE SSGEO(IERROR,IPH,PARM1,PARM2,PARM3,PARM4,              &
     &  PSIPO,G,MSOFF,ICH1,KNTRVL,SOLAR_ZENITH)

!     THIS ROUTINE CALLS THE LOWTRAN GEOMETRY ROUTINES REPEATEDLY
!     TO OBTAIN THE ABSORBER AMOUNTS FROM THE SCATTERING POINTS ON
!     THE OPTICAL PATH TO THE EXTRATERRESTRIAL SOURCE, AS IS NECESSARY
!     FOR THE LAYER BY LAYER SINGLE SCATTERING RADIANCE CALCULATION.

!     INPUT ARGUMENTS:
!       IERROR   ERROR FLAG (=0 FOR SUCCESSFUL CALLS)
!       IPH      PHASE FUNCTION FLAG (=0 FOR HENYEY-GREENSTEIN)
!                                    (=1 FOR USER-DEFINED)
!                                    (=2 FOR MIE-GENERATED)
!       PARM1    OBSERVER LATITUDE [DEG NORTH OF EQUATOR]
!       PARM2    OBSERVER LONGITUDE [DEG WEST OF GREENWICH]
!       PARM3    SOURCE LATITUDE [DEG NORTH OF EQUATOR]
!       PARM4    SOURCE LONGITUDE [DEG WEST OF GREENWICH]
!       PSIPO    PATH AZIMUTH ANGLE [DEG EAST OF NORTH]
!       G        HENYEY-GREENSTEIN ASYMMETRY FACTOR (IPH=0)
!       MSOFF    MULTIPLE SCATTERING LAYER OFFSET (EQUALS 0 FOR
!                OPTICAL PATH AND EQUALS LAYTWO FOR VERTICAL PATH)
!       ICH1     HAZE MODEL FLAG
!       KNTRVL   NUMBER OF CORRELATED-K SUB-INTERVALS
!                (=1 IF CORRELATED-K APPROACH IS NOT USED)
      INTEGER IERROR,IPH,MSOFF,ICH1,KNTRVL
      REAL PARM1,PARM2,PARM3,PARM4,PSIPO,G
      REAL SOLAR_ZENITH !DRF 

!     PARAMETERS:
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:
      REAL WX
      COMMON/NONAME/WX(NMOLX)
      REAL TBBYSX
      COMMON/SOLSX/TBBYSX(LAYTHR,NMOLX)
      INCLUDE 'BASE.h'
      INCLUDE 'SOLS.h'
      INCLUDE 'IFIL.h'

!     /CARD1/
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,MODTRN
      INTEGER NCRALT,NCRSPC
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      &
     &  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP

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
      REAL RE,ZMAX
      INTEGER IPATH
      COMMON/PARMTR/RE,ZMAX,IPATH

!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI.
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

!     /RFRPTH/
      REAL SPZP,SPPP,SPTP,SPPPSM,SPTPSM,SPRHOP,SPDENP,SPAMTP,DRANGE
      COMMON/RFRPTH/SPZP(LAYDIM+1),SPPP(LAYDIM+1),SPTP(LAYDIM+1),       &
     &  SPPPSM(LAYDIM+1),SPTPSM(LAYDIM+1),SPRHOP(LAYDIM+1),             &
     &  SPDENP(MEXT,LAYDIM+1),SPAMTP(MEXT,LAYDIM+1),DRANGE(LAYDIM+1)

!     /MSRD/
!       CSSCAT   COSINE OF THE SCATTERING ANGLE.
!                (AT H1 IF IMULT=1; AT OR "NEAR" H2 IF IMULT=-1)
!       SLEGEN   Nth LEGENDRE POLYNOMIAL EVALUATED AT THE COSINE OF THE
!                SCATTERING ANGLE TIMES (2N+1)/4pi (N=0 TO NSTR-1).
!       CSZEN0   LAYER BOUNDARY COSINE OF SOLAR/LUNAR ZENITH.
!       CSZEN    LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH.
!       CSZENX   AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING
!                (AWAY FROM EARTH) THE CURRENT LAYER.
!       BBGRND   THERMAL EMISSION (FLUX) AT THE GROUND [W CM-2 / CM-1].
!       BBNDRY   LAYER BOUNDARY THERMAL EMISSION (FLUX) [W CM-2 / CM-1].
!       TCONT    LAYER CONTINUUM OPTICAL DEPTH.
!       TAUT     LAYER TOTAL OPTICAL DEPTH.
!       GTSCAT   SUM OVER SCATTERING SOURCES OF SCATTERING OPTICAL DEPTH
!                AND PHASE FUNCTION LEGENDRE COEFFICIENT PRODUCTS.
!       COSBAR   LAYER EFFECTIVE SCATTERING ASYMMETRY FACTOR.
!       DEPRAT   FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN.
!       S0DEP    OPTICAL DEPTH FROM LAYER BOUNDARY TO SUN.
!       S0TRN    TRANSMITTED SOLAR IRRADIANCES [W CM-2 / CM-1]
!       UPF      LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].
!       DNF      LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].
!       UPFS     LAYER BOUNDARY UPWARD SOLAR FLUX [W CM-2 / CM-1].
!       DNFS     LAYER BOUNDARY DOWNWARD SOLAR FLUX [W CM-2 / CM-1].
      REAL CSSCAT,SLEGEN,CSZEN0,CSZEN,CSZENX,TCONT,TAUT,GTSCAT,COSBAR,  &
     &  BBGRND,BBNDRY,S0DEP,S0TRN,DEPRAT,UPF,DNF,UPFS,DNFS
      COMMON/MSRD/CSSCAT,SLEGEN(0:MAZ),                                 &
     &  CSZEN0(LAYDIM),CSZEN(LAYDIM),CSZENX(LAYDIM),TCONT(LAYDIM),      &
     &  TAUT(MXKSUB,LAYDIM),GTSCAT(0:MXCMU,1:LAYDIM),COSBAR(LAYDIM),    &
     &  BBGRND,BBNDRY(LAYDIM),S0DEP(MXKSUB,LAYTWO),S0TRN(MXKSUB,LAYTWO),&
     &  DEPRAT(MXKSUB,LAYDIM),UPF(MXKSUB,LAYDIM),DNF(MXKSUB,LAYDIM),    &
     &  UPFS(MXKSUB,LAYDIM),DNFS(MXKSUB,LAYDIM)

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

!     /DISRT/
!       DIS      LOGICAL FLAG, TRUE FOR DISORT MULTIPLE SCATTERING.
!       DISAZM   LOGICAL FLAG, TRUE FOR DISORT WITH AZIMUTH DEPENDENCE.
!       DISALB   LOGICAL FLAG, TRUE FOR DISORT SPHERICAL ALBEDO OPTION.
!       LDISCL   LOGICAL FLAG, TRUE FOR ISAACS SCALED TO DISORT.
!       NSTR     NUMBER OF DISCRETE ORDINATE STREAMS.
!       NAZ      NUMBER OF DISORT AZIMUTH COMPONENTS.
!       N2GAUS   ORDER OF DOUBLE-GAUSS QUADRATURES.
      LOGICAL DIS,DISAZM,DISALB,LDISCL
      INTEGER NSTR,NAZ,N2GAUS
      COMMON/DISRT/DIS,DISAZM,DISALB,LDISCL,NSTR,NAZ,N2GAUS

!     /CARD4/
!       IV1      LOWEST FREQUENCY OUTPUT [CM-1].
!       IV2      HIGHEST FREQUENCY OUTPUT [CM-1].
!       IDV      PRINTOUT FREQUENCY STEP SIZE [CM-1].
!       IFWHM    TRIANGULAR SLIT FULL-WIDTH-HALF-MAXIMUM [CM-1].
!       IVX      CURRENT COMPUTATION FREQUENCY [CM-1].
!       IVOFF    OFFSET BETWEEN COMPUTATION AND OUTPUT FREQUENCIES,
!                REQUIRED FOR SLIT FUNCTION [CM-1].
!       IWRITE   COMPUTATION FREQUENCY OF NEXT WRITE [CM-1].
!       NSPCDT   NUMBER OF OUTPUT SPECTRAL DATA POINTS.
!       NWGT     NUMBER OF SPECTRAL BINS CONTRIBUTING TO SLIT FUNCTION.
!       WGT      NORMALIZED WEIGHTS USED TO DEFINE THE SLIT FUNCTION.
      INTEGER IV1,IV2,IDV,IFWHM,IVX,IVOFF,IWRITE,NSPCDT,NWGT
      REAL WGT
      COMMON/CARD4/IV1,IV2,IDV,IFWHM,IVX,IVOFF,IWRITE,NSPCDT,           &
     &  NWGT,WGT(NBINS)

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

!     /CJM5/
!       AMOD3D   FLAG INDICATING OUTPUT DATABASE FILE TYPE:
      CHARACTER AMOD3D*1
      COMMON/CJM5/AMOD3D

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       IANGLS  LOOP INDEX FOR SCATTERING PHASE FUNCTION ANGLES.
!       PSIO    RELATIVE AZIMUTH ANGLE AT SENSOR [DEG].
!       AUMU1   ABSOLUTE VALUE OF COSINE OF THE PATH NADIR ANGLE.
!       UMU1L   COSINE OF THE PATH NADIR ANGLE AT SCATTERING POINT L.
!       AUMU1L  ABSOLUTE VALUE OF UMU1L.
!       ISTR    DISCRETE ORDINATE STREAM INDEX.
      CHARACTER MESSAG*48
      INTEGER IANGLS,MSOFFX,JTRNSV,IKMAX1,LENNSV,ITYPSV,J,MSOFFJ,       &
     &  K,IARBO,IARB,LM1,L,LJ0,JITER,MSOFFL,KP,IK,MSOFFK,LNMESS,ISTR
      REAL THETAO,PHIO,THETAS,PHIS,PSIO,DELO,H1SAV,H2SAV,ANGSAV,ANGMX,  &
     &  RNGSAV,BETASV,PSIPOS,BETAST,PSIST,ANGL0,RELH,THTST,ANGERR,      &
     &  BENDNG,ANGSV,ANGMN,PNUM,TNUM,DENOM,WPTH,AUMU1,UMU1L,AUMU1L

!     LOCAL ARRAYS:
!       MU       COSINE OF SCATTERING ANGLE.
!       WLH2O    SENSOR-LEVEL_BOUNDARY-SUN PATH H2O COLUMN [ATM CM].
!       WLCO2    SENSOR-LEVEL_BOUNDARY-SUN PATH CO2 COLUMN [ATM CM].
!       WLO3     SENSOR-LEVEL_BOUNDARY-SUN PATH O3 COLUMN [ATM CM].
!       WH2O     SENSOR-LEVEL_BOUNDARY PATH H2O COLUMN [ATM CM].
!       WCO2     SENSOR-LEVEL_BOUNDARY PATH CO2 COLUMN [ATM CM].
!       WO3      SENSOR-LEVEL_BOUNDARY PATH O3 COLUMN [ATM CM].
      INTEGER LJSAV(LAYTWO)
      REAL WPSAV(LAYTWO,MEXT),TBSAV(LAYTWO),WSAV(MEXT),ZPSAV(LAYDIM+1), &
     &  PSAV(LAYTWO),RHSAV(LAYDIM+1),WPSUM(MEXT),WPSAVX(LAYTHR,NMOLX),  &
     &  MU(1),WLH2O(LAYTWO),WLCO2(LAYTWO),WLO3(LAYTWO),                 &
     &  WH2O(LAYTWO),WCO2(LAYTWO),WO3(LAYTWO)

!     LOCAL FUNCTIONS:
      REAL PSI,DEL,COSSCT,HENGNS,DENLAY

!     INITIALIZATIONS:
      MSOFFX=0
      NPR=2
      ISSGEO=1
      !WRITE(*,*) 'INTO SSGEO,SOLAR_ZEN'
!     SPECIFY THE GEOMETRICAL CONFIGURATION
      THETAO=PARM1
      PHIO=PARM2
      THETAS=PARM3
      PHIS=PARM4
      IF(THETAO.LT.-89.5)THEN

!         OBSERVER IS AT OR NEAR THE SOUTH POLE.
!         REMAP TO THE EQUATOR.
          IF(.NOT.LJMASS) WRITE(IPR,'(/2A)')                            &
     &      '  THETAO < -89.5 (OBSERVER AT OR NEAR THE',                &
     &      ' SOUTH POLE).  PROBLEM HAS BEEN REMAPPED TO THE EQUATOR.'
          PSIPO=PSIPO-PHIS
          THETAO=0.
          PHIO=0.
          THETAS=0.
          PHIS=90.+THETAS
      ELSEIF(THETAO.GT.89.5)THEN

!         OBSERVER IS AT OR NEAR THE NORTH POLE.  REMAP TO EQUATOR.
          IF(.NOT.LJMASS) WRITE(IPR,'(/2A)')                            &
     &      '  THETAO > +89.5 (OBSERVER AT OR NEAR THE',                &
     &      ' NORTH POLE).  PROBLEM HAS BEEN REMAPPED TO THE EQUATOR.'
          PSIPO=PHIS-PSIPO
          THETAO=0.
          PHIO=0.
          THETAS=0.
          PHIS=90.-THETAS
      ENDIF

!     SAVE OPTICAL PATH PARAMETERS AND AMOUNTS
      JTRNSV=JTURN
      IKMAX1=IKMAX+1
      H1SAV=H1
      H2SAV=H2
      ANGSAV=ANGLE
      RNGSAV=RANGE
      BETASV=BETA
      BETA=0.
      LENNSV=LENN
      ITYPSV=ITYPE
      DO J=1,ML
          ZPSAV(J)=SPZP(J)
          RHSAV(J)=RELHUM(J)
      ENDDO
      DO J=1,IKMAX1
          MSOFFJ=MSOFF+J
          TBSAV(J)=TBBY(MSOFFJ)
          PSAV(J)=PATM(MSOFFJ)
          LJSAV(J)=LJ(J)
          IF(LJSAV(J).GT.ML)LJSAV(J)=ML
          DO K=1,MEXT
              WPSAV(J,K)=WPATH(MSOFFJ,K)
          ENDDO
          DO K=1,NMOLX
              WPSAVX(J,K)=WPATH(MSOFFJ,MEXT+K)
          ENDDO
      ENDDO
      DO K=1,MEXT
          WSAV(K)=W(K)
          WPSUM(K)=0.
      ENDDO

!     ESTABLISH PSIO AND DELO
      IARBO=0
      IF(ANGLE.LT..01 .OR. ANGLE.GT.179.99)IARBO=1
      IARB=IARBO
      BETAST=0.
      CALL PSIECA(THETAO,PHIO,THETAS,PHIS,PSIPOS,DELO)
      PSIO=PSIPOS-PSIPO
      IF(PSIO.GT.180.)PSIO=PSIO-360.
      IF(PSIO.LT.-180.)PSIO=PSIO+360.
      PSIST=PSIO
      ANGL0=DELO

!     LOOP OVER THE POINT TO SUN PATHS TO OBTAIN AMOUNTS
      IF(MSOFF.EQ.0)THEN
          IF(.NOT.LJMASS)WRITE(IPR,'(//(2A))')                          &
     &      '1 SINGLE SCATTER SOLAR PATH GEOMETRY TABLE',               &
     &      ' FOR THE LINE-OF-SIGHT PATH',                              &
     &      '  ----------------------------------------',               &
     &      '---------------------------'
      ELSE
          IF(.NOT.LJMASS)WRITE(IPR,'(//(3A))')                          &
     &      '1 SINGLE SCATTER SOLAR PATH GEOMETRY TABLE',               &
     &      ' FOR THE VERTICAL GROUND-TO-SPACE PATH',                   &
     &      ' (THIS PATH IS USED FOR MULTIPLE SCATTERING)',             &
     &      '  ----------------------------------------',               &
     &      '--------------------------------------',                   &
     &      '--------------------------------------------'
      ENDIF
      IF(.NOT.LJMASS)WRITE(IPR,'(/(2A))')' SCT    SCATTER  SUBTENDED',  &
     &  '      SOLAR       PATH   RELATIVE    SCATTER',                 &
     &                    ' PNT   ALTITUDE      ANGLE',                 &
     &  '     ZENITH     ZENITH    AZIMUTH      ANGLE',                 &
     &                    '           [KM]      [DEG]',                 &
     &  '      [DEG]      [DEG]      [DEG]      [DEG]'
      LM1=0
      DO L=1,IKMAX1
          IF(LENNSV.EQ.0 .AND. JTRNSV.EQ.0)THEN

!             SHORT PATH, UP
              IF(L.GE.2)BETAST=BETAST+ADBETA(LM1)
              H1=ZPSAV(L)
              RELH=RHSAV(L)
              THTST=ATHETA(L)
          ELSE

!             LONG PATH, OR SHORT PATH DOWN
              IF(L.GE.2)BETAST=BETAST+ADBETA(LJSAV(LM1))
              LJ0=LJSAV(L)
              IF(L.LT.JTRNSV)LJ0=LJ0+1
              THTST=ATHETA(LJ0)
              IF(L.LE.JTRNSV)THTST=180.-THTST
              H1=ZPSAV(LJ0)
              RELH=RHSAV(LJ0)
          ENDIF
          AH1(L)=H1
          ARH(L)=RELH
          IF(L.GE.2)THEN
              PSIST=PSI(PSIO,DELO,BETAST,IARB,IARBO)
              ANGL0=DEL(PSIO,DELO,BETAST,IARBO)
          ENDIF
          ANGERR=0.
          ANGMX=ANGL0
          BENDNG=0.
          ANGLE=ANGL0

!         RANGE=UNKNOWN
          ITYPE=3
          DO JITER=1,12
              LNMESS=1
              MESSAG=' '

!             SET H2 TO ZERO TO INSURE THAT ANGLE
!             IS USED TO DEFINE SOLAR PATH.
              H2=0.
              LENN=0
              IF(ANGLE.GT.90.)THEN
                  IF(ANGLE.LE.90.0001)THEN

!                     MULTIPLE SCATTERING CORRECTION FOR
!                     SCATTERING POINT TO SUN PATHS.
                      ANGLE=90.
                  ELSE
                      LENN=1
                      LNMESS=48
                      MESSAG(1:23)='THIS SOLAR PATH PASSES '
                      MESSAG(24:48)='THROUGH A TANGENT HEIGHT.'
                  ENDIF
              ENDIF
              IF(H1.GE.ZMAX .AND. LENN.EQ.0)THEN

!                 SCATTERING POINT IS AT OR ABOVE ZMAX
!                 AND LENN=0;  SET W(K)=0.0 AND CONTINUE.
                  DO K=1,MEXT
                      W(K)=0.
                  ENDDO
                  DO K=1,NMOLX
                      WX(K)=0.
                  ENDDO
                  GOTO 10
              ENDIF
              ANGSV=ANGLE
              CALL GEO(IERROR,BENDNG,MSOFFX,ICH1)
              ANGLE=ANGSV

!             IERROR=-5 IF SCATTERING POINT IS SHADED; SET W(36)=-5.
              IF(IERROR.EQ.-5)THEN
                  LNMESS=38
                  MESSAG='THIS SCATTERING POINT IS IN THE SHADE.'
                  W(36)=-5.
                  IERROR=0
                  GOTO 10
              ENDIF

!             SOLAR ZENITH ERROR
              ANGERR=ANGLE+BENDNG-ANGL0
              IF(ABS(ANGERR).LT..001)GOTO 10

!             CALCULATE ANGLE USING THE BISECTION METHOD.  ANGL0 IS
!             A MAXIMUM AND ANGL0 MINUS ITS BENDING IS A MINIMUM.
              IF(JITER.EQ.1)THEN
                  ANGMN=ANGL0-BENDNG
                  ANGLE=ANGMN
              ELSE
                  IF(ANGERR.GT.0)THEN
                      ANGMX=ANGLE
                  ELSE
                      ANGMN=ANGLE
                  ENDIF
                  ANGLE=.5*(ANGMX+ANGMN)
              ENDIF
          ENDDO
          LNMESS=40
          MESSAG='THE SOLAR ZENITH ANGLE DID NOT CONVERGE.'
          WRITE(IPR,'(2A,F12.5,A)')                                     &
     &      ' AFTER 12 ITERATIONS, THE SOLAR ZENITH EXITING',           &
     &      ' THE ATMOSPHERE IS STILL IN ERROR BY',ANGERR,'DEG'
   10     CONTINUE
          PHSCOS(L)=COSSCT(ANGLE,THTST,PSIST,IARB)
          PHSANG(L)=DEG*ACOS(PHSCOS(L))

!         BRANCH BASED ON SCATTERING PHASE FUNCTION TYPE:
          IF(IPH.EQ.0)THEN

!             LOAD HENYEY-GREENSTEIN PHASE FUNCTION ARRAY:
              PHASFN(L,2)=HENGNS(G,PHSCOS(L))
          ELSEIF(IPH.EQ.1)THEN

!             SET UP SCATTERING ANGLE INTERPOLATION:
              IF(PHSANG(L).LE.ANGF(1))THEN
                  IUPPHS(L)=2
                  PHSFAC(L)=0.
              ELSE
                  DO IANGLS=2,NANGLS
                      IF(PHSANG(L).LT.ANGF(IANGLS))THEN
                          IUPPHS(L)=IANGLS
                          PHSFAC(L)=(PHSANG(L)-ANGF(IANGLS-1))          &
     &                      /(ANGF(IANGLS)-ANGF(IANGLS-1))
                          GOTO 20
                      ENDIF
                  ENDDO
                  IUPPHS(L)=NANGLS
                  PHSFAC(L)=1.
              ENDIF
          ENDIF

!         LOAD WATER DROPLET HENYEY-GREENSTEIN PHASE FUNCTION.
   20     CONTINUE
          IF(ABS(ASYMWD).LT.1.)PHASFN(L,3)=HENGNS(ASYMWD,PHSCOS(L))

!         LOAD ICE PARTICLE HENYEY-GREENSTEIN PHASE FUNCTION.
          IF(ABS(ASYMIP).LT.1.)PHASFN(L,4)=HENGNS(ASYMIP,PHSCOS(L))

!         STORE SOLAR ZENITH DATA FOR MULTIPLE SCATTERING.
          IF(MSOFF.GT.0)CALL SOLZEN(L,IKMAX1,REE,DEG,ANGLE)

!         WRITE SOLAR SCATTER PATH DATA.
          IF(.NOT.LJMASS)WRITE(IPR,'(I4,6F11.5,3X,A)')                  &
     &      L,H1,BETAST,ANGLE,THTST,PSIST,PHSANG(L),MESSAG(1:LNMESS)
          SOLAR_ZENITH = ANGLE !DRF
          !WRITE(*,*) 'ANGLE = ',ANGLE

!         STORE SOLAR AND PATH ANGLES FOR DISORT CALCULATION:
          IF(MSOFF.EQ.0)THEN

!             LINE-OF-SIGHT PATH:
              IF(.NOT.DIS)THEN

!                 SET UMU0 FOR CALCULATION OF DIRECTIONAL REFLECTIVITY:
                  UMU0=CSZEN0(1)
              ELSEIF(IMULT.EQ.1 .AND. L.EQ.1)THEN

!                 ANGLES AT SENSOR (H1) ARE TO BE USED:
                  PHI1=ABS(PSIST)
                  UMU1=-COS(THTST/DEG)
                  UMU0=COS(ANGLE/DEG)
                  CSSCAT=PHSCOS(1)
                  MU(1)=PHSCOS(1)
                  CALL LEPOLY(1,0,MAZ,NSTR-1,MU,SLEGEN)
                  DO ISTR=0,NSTR-1
                      SLEGEN(ISTR)=SLEGEN(ISTR)*(ISTR+.5)/(2*pi)
                  ENDDO
              ELSEIF(IMULT.EQ.-1)THEN

!                 ANGLES AT H2 ARE TO BE USED:
                  IF(H2SAV.NE.0. .AND. ITYPE.EQ.3)THEN

!                     H2 IS A TANGENT POINT, BUT DISORT (BEING PLANE
!                     PARALLEL) WILL NOT ACCEPT A 90 DEG PATH_ZENITH.
!                     REQUIRE |COS(PATH_ZENITH)| TO EXCEED .05
!                     (PATH_ZENITH < 87.134 OR PATH_ZENITH > 92.866)
                      IF(L.EQ.1)THEN

!                         INITIALIZE THE ANGLES AT H1:
                          PHI1=ABS(PSIST)
                          UMU1=-COS(THTST/DEG)
                          UMU0=COS(ANGLE/DEG)
                          CSSCAT=PHSCOS(1)
                          MU(1)=PHSCOS(1)
                          CALL LEPOLY(1,0,MAZ,NSTR-1,MU,SLEGEN)
                          DO ISTR=0,NSTR-1
                              SLEGEN(ISTR)=SLEGEN(ISTR)*(ISTR+.5)/(2*pi)
                          ENDDO
                          AUMU1=ABS(UMU1)
                      ELSE
                          UMU1L=-COS(THTST/DEG)
                          AUMU1L=ABS(UMU1L)
                          IF(AUMU1L.GT..05 .AND. AUMU1L.LT.AUMU1)THEN
                              PHI1=ABS(PSIST)
                              UMU0=COS(ANGLE/DEG)
                              UMU1=UMU1L
                              CSSCAT=PHSCOS(L)
                              MU(1)=PHSCOS(L)
                              CALL LEPOLY(1,0,MAZ,NSTR-1,MU,SLEGEN)
                              DO ISTR=0,NSTR-1
                                  SLEGEN(ISTR)                          &
     &                              =SLEGEN(ISTR)*(ISTR+.5)/(2*pi)
                              ENDDO
                              AUMU1=ABS(UMU1)
                          ENDIF
                      ENDIF
                  ELSEIF(L.EQ.IKMAX1)THEN
                      PHI1=ABS(PSIST)
                      UMU1=-COS(THTST/DEG)
                      UMU0=COS(ANGLE/DEG)
                      CSSCAT=PHSCOS(IKMAX1)
                      MU(1)=PHSCOS(IKMAX1)
                      CALL LEPOLY(1,0,MAZ,NSTR-1,MU,SLEGEN)
                      DO ISTR=0,NSTR-1
                          SLEGEN(ISTR)=SLEGEN(ISTR)*(ISTR+.5)/(2*pi)
                      ENDDO
                  ENDIF
              ENDIF
          ENDIF

!         LOAD AMOUNTS FROM W(K) INTO WPATHS(L,K)
          MSOFFL=MSOFF+L
          IF(MSOFF.EQ.0 .AND. L.GT.1 .AND.                              &
     &      W(36).GE.0. .AND. KNTRVL.EQ.1)THEN

!             ADD OBSERVER TO SCATTERING POINT PATH AMOUNTS TO
!             W FOR ALL BUT THE MODTRAN BAND MODEL ABSORBERS.
              DO K=1,MEXT
                  WPSUM(K)=WPSUM(K)+WPSAV(LM1,K)
                  WPATHS(MSOFFL,K)=W(K)+WPSUM(K)
              ENDDO
              IF(MODTRN)THEN
                  WLH2O(L)=WPATHS(L,17)
                  WLCO2(L)=WPATHS(L,36)
                  WLO3(L)=WPATHS(L,31)
                  WH2O(L)=WPSUM(17)
                  WCO2(L)=WPSUM(36)
                  WO3(L)=WPSUM(31)
                  DO K=1,NMOL
                      KP=KPOINT(K)
                      WPATHS(MSOFFL,KP)=W(KP)
                  ENDDO
              ENDIF
          ELSE
              DO K=1,MEXT
                  WPATHS(MSOFFL,K)=W(K)
              ENDDO
          ENDIF
          DO K=1,NMOLX
              WPATHS(MSOFFL,MEXT+K)=WX(K)
          ENDDO

!         WHEN THE MODERATE RESOLUTION OPTION IS USED, A CURTIS-GODSON
!         AVERAGE PRESSURE, PATMS, AND TEMPERATURE, TBBYS, IS DEFINED
!         FOR THE SCATTERING POINT TO EXTRATERRESTRIAL SOURCE "LAYER"
          IF(MODTRN)THEN
              DO K=1,NMOLXT
                  KP=KPOINT(K)
                  TNUM=0.
                  DENOM=0.
                  IF(K.LE.NMOL)THEN
                      PNUM=0.
                      DO IK=1,IKMAX
                          MSOFFK=MSOFFX+IK
                          WPTH=WPATH(MSOFFK,KP)
                          IF(WPTH.GT.0.)THEN
                              PNUM=PNUM+PATM(MSOFFK)*WPTH
                              TNUM=TNUM+TBBY(MSOFFK)*WPTH
                              DENOM=DENOM+WPTH
                          ENDIF
                      ENDDO
                      IF(DENOM.GT.0.)THEN
                          PATMS(MSOFFL,K)=PNUM/DENOM
                          TBBYS(MSOFFL,K)=TNUM/DENOM
                      ELSE
                          PATMS(MSOFFL,K)=PATM(MSOFFX+1)
                          TBBYS(MSOFFL,K)=TBBY(MSOFFX+1)
                      ENDIF
                  ELSE
                      DO IK=1,IKMAX
                          MSOFFK=MSOFFX+IK
                          WPTH=WPATH(MSOFFK,KP)
                          IF(WPTH.GT.0.)THEN
                              TNUM=TNUM+TBBY(MSOFFK)*WPTH
                              DENOM=DENOM+WPTH
                          ENDIF
                      ENDDO
                      IF(DENOM.NE.0.)THEN
                          TBBYSX(MSOFFL,K-NMOL)=TNUM/DENOM
                      ELSE
                          TBBYSX(MSOFFL,K-NMOL)=TBBY(MSOFFX+1)
                      ENDIF
                  ENDIF
              ENDDO
          ENDIF
          LM1=L
      ENDDO
      IF(AMOD3D.EQ.'C')THEN
          WRITE(IDBOUT,'((I13,A))')                                     &
     &      INT(1.E7/(IV2+.5)+1.5),' MINIMUM WAVELENGTH [NM].',         &
     &      INT(1.E7/(IV1+.5)+0.5),' MAXIMUM WAVELENGTH [NM].'
          WRITE(IDBOUT,'(I13,A,/(F13.4,A))')                            &
     &      INT(1000*H1SAV+.5),' SENSOR ALTITUDE [M].',                 &
     &      THTST,' SENSOR ZENITH AT GROUND [DEG].',                    &
     &      ANGLE,' SOLAR  ZENITH AT GROUND [DEG].'
          WRITE(IDBOUT,'((I13,A))')                                     &
     &      IKMAX,' NUMBER OF ATMOSPHERE LAYERS.',                      &
     &      IKMAX1,' NUMBER OF PATH LEVELS.'
          WRITE(IDBOUT,'(/2A,/(I6,1X,0P,F12.7,1X,F9.3),                 &
     &      /(I6,1X,0P,F12.7,1X,F9.3,3(1X,1P,E15.7)))')                 &
     &      'ALT[M]   PRES[MBAR]   TEMP[K]',                            &
     &      '  H2O[ATM-CM/HM]  CO2[ATM-CM/HM]   O3[ATM-CM/HM]',         &
     &      INT(1000*ZM(ML)+.5),DBLE(PM(ML)),TM(ML),                    &
     &      (INT(1000*ZM(L)+.5),DBLE(PM(L)),TM(L),                      &
     &      .1*DENLAY(DENSTY(17,L),DENSTY(17,L+1)),                     &
     &      .1*DENLAY(DENSTY(36,L),DENSTY(36,L+1)),                     &
     &      .1*DENLAY(DENSTY(31,L),DENSTY(31,L+1)),L=IKMAX,1,-1)
          WRITE(IDBOUT,'(/2A,/(I6,6(1X,1P,E15.7)))')                    &
     &      'ALT[M]   L-H2O[ATM-CM]   L-CO2[ATM-CM]    L-O3[ATM-CM]',   &
     &            '     H2O[ATM-CM]     CO2[ATM-CM]      O3[ATM-CM]',   &
     &      INT(1000*AH1(1)+.5),WPATHS(1,17),WPATHS(1,36),WPATHS(1,31), &
     &      0.,0.,0.,(INT(1000*AH1(L)+.5),WLH2O(L),WLCO2(L),WLO3(L),    &
     &      WH2O(L),WCO2(L),WO3(L),L=2,IKMAX1)
      ENDIF
      ANGSUN=ANGLE
      CSNSRF=COS(ANGLE/DEG)
      IF(PSIST.GE.0.)THEN
          AZMSRF=PI-PSIST/DEG
      ELSE
          AZMSRF=-PSIST/DEG-PI
      ENDIF

!     RESTORE OPTICAL PATH AMOUNTS
      IKMAX=IKMAX1-1
      H1=H1SAV
      H2=H2SAV
      ANGLE=ANGSAV
      RANGE=RNGSAV
      BETA=BETASV
      LENN=LENNSV
      ITYPE=ITYPSV
      NPR=NOPRNT
      DO J=1,IKMAX1
          MSOFFJ=MSOFF+J
          TBBY(MSOFFJ)=TBSAV(J)
          PATM(MSOFFJ)=PSAV(J)
          LJ(J)=LJSAV(J)
          DO K=1,MEXT
              WPATH(MSOFFJ,K)=WPSAV(J,K)
          ENDDO
          DO K=1,NMOLX
              WPATH(MSOFFJ,MEXT+K)=WPSAVX(J,K)
          ENDDO
      ENDDO
      DO K=1,MEXT
          W(K)=WSAV(K)
      ENDDO

!     RETURN TO DRIVER OR H2SRC
      RETURN
      END
