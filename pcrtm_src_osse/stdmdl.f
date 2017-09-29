      SUBROUTINE STDMDL(ICH1,MARIC1,MARK,H2OAER,VRFRAC,LNOVAM)

!     THIS SUBROUTINE LOADS ONE OF THE 6 STANDARD ATMOSPHERIC
!     PROFILES INTO COMMON /MPROF/ AND CALCULATES THE
!     DENSITIES OF THE VARIOUS ABSORBING GASES AND AEROSOLS.
      IMPLICIT NONE

!     INPUT ARGUMENTS:
!       ICH1     NUMERIC LABEL FOR BOUNDARY LAYER AEROSOL MODEL.
!       MARIC1   MARINE AEROSOL REGION NUMBER (=0 IF NOT USED).
!       MARK     RELATIVE HUMIDITY LAYER INDEX FOR MARINE AEROSOL.
!       H2OAER   FLAG, TRUE IF AEROSOL DEFAULT PROFILE PROPERTIES
!                ARE UPDATED AFTER WATER COLUMN SCALING.
!       VRFRAC   SPECTRAL FREQUENCY USED FOR INDEX OF REFRACTION [CM-1].
!       LNOVAM   LOGICAL FLAG, .TRUE. IF NOVAM AEROSOLS ARE USED.
      INTEGER ICH1,MARIC1,MARK
      CHARACTER H2OAER*1
      REAL VRFRAC
      LOGICAL LNOVAM

!     PARAMETERS:
!       T0RAT    RATIO OF FASCODE STANDARD TEMPERATURE (296K)
!                TO MODTRAN STANDARD TEMPERATURE (273.15K)
!       V550NM   FREQUENCY EQUIVALENT OF 550 NM [CM-1].
!       H2O_P    CONVERSION FROM GM H2O / M3 TO PRESSURE [ATM] AT STP.
!       RAYDEN   RAYLEIGH SCATTERING CROSS SECTION DENOMINATOR FOR
!                "STANDARD" AIR (SEE USE BELOW).
      INCLUDE 'PARAMS.h'
      REAL T0RAT,V550NM,H2O_P,RAYDEN
      PARAMETER(T0RAT=296/273.15,H2O_P=STDVOL/(1E6*H2OMWT),             &
     &  V550NM=1E4/.55,RAYDEN=67.6*.781+60.5*.00934+3.5*.000018+        &
     &  158.1*.000360+151.2*.0000017+56.5*.209+13.4*.00000058+.0000052)

!     COMMONS:
      INCLUDE 'YPROP.h'
      INCLUDE 'BASE.h'
      INCLUDE 'IFIL.h'

!     /NAMEX/
!       CNAMEX   NAME OF CROSS-SECTION (X) SPECIES.
      CHARACTER CNAMEX*8
      COMMON/NAMEX/CNAMEX(NMOLX)

!     /MDATXY/
!       WMOLXT   CROSS-SECTION MOLECULE DENSITY PROFILE [PPMV].
!       WMOLYT   AUXILIARY (Y) MOLECULE DENSITY PROFILE [PPMV].
      REAL WMOLXT,WMOLYT
      COMMON/MDATXY/WMOLXT(NMOLX,LAYDIM),WMOLYT(MMOLY,LAYDIM)

!     /MPROF/
!       ZM       PROFILE LEVEL ALTITUDES [KM].
!       PM       PROFILE LEVEL PRESSURES [MBAR].
!       TM       PROFILE LEVEL TEMPERATURES [K].
!       RFNDX    PROFILE LEVEL REFRACTIVITIES.
!       LRHSET   FLAG, .TRUE. IF RELATIVE HUMIDITY IS NOT TO BE SCALED.
      DOUBLE PRECISION ZM
      REAL PM,TM,RFNDX
      LOGICAL LRHSET
      COMMON/MPROF/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    &
     &  RFNDX(LAYDIM),LRHSET(LAYDIM)

!     /DEN/
!       DENSTY   PROFILE LEVEL DENSITIES [ATM CM / KM FOR MOST SPECIES].
      REAL DENSTY
      COMMON/DEN/DENSTY(0:MEXTXY,1:LAYDIM)

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

!     /M_PTWO/
!       PPROF    PRESSURE PROFILE [MB].
!       TPROF    TEMPERATURE PROFILE [K].
!       WH2O     H2O VOLUME MIXING RATIO PROFILE [GM / M3].
!       WO3      O3 VOLUME MIXING RATIO PROFILE [PPMV].
      REAL PPROF,TPROF,WH2O,WO3
      COMMON/M_PTWO/PPROF(LAYDIM),TPROF(LAYDIM),WH2O(LAYDIM),WO3(LAYDIM)

!     /M_UMIX/
!       WCO2     CO2 VOLUME MIXING RATIO PROFILE [PPMV].
!       WN2O     N2O VOLUME MIXING RATIO PROFILE [PPMV].
!       WCO      CO VOLUME MIXING RATIO PROFILE [PPMV].
!       WCH4     CH4 VOLUME MIXING RATIO PROFILE [PPMV].
!       WO2      O2 VOLUME MIXING RATIO PROFILE [PPMV].
!       WNO      NO VOLUME MIXING RATIO PROFILE [PPMV].
!       WSO2     SO2 VOLUME MIXING RATIO PROFILE [PPMV].
!       WNO2     NO2 VOLUME MIXING RATIO PROFILE [PPMV].
!       WHNO3    HNO3 VOLUME MIXING RATIO PROFILE [PPMV].
      REAL WCO2,WN2O,WCO,WCH4,WO2,WNO,WSO2,WNO2,WNH3,WHNO3
      COMMON/M_UMIX/WCO2(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),              &
     &  WCH4(LAYDIM),WO2(LAYDIM),WNO(LAYDIM),WSO2(LAYDIM),              &
     &  WNO2(LAYDIM),WNH3(LAYDIM),WHNO3(LAYDIM)

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

!     /AER/ THERE ARE "MAER=17" PARTICULATE COMPONENTS:
!           1       AEROSOL 1 (NOMINALLY, BOUNDARY LAYER AEROSOL).
!           2       AEROSOL 2 (NOMINALLY, TROPOSPHERIC AEROSOL).
!           3       AEROSOL 3 (NOMINALLY, STRATOSPHERIC AEROSOL).
!           4       AEROSOL 4 (NOMINALLY, VOLCANIC AEROSOL).
!           5       CIRRUS CLOUD.
!           6       CLOUD 1 (NOMINALLY, WATER CLOUD).
!           7       CLOUD 2 (NOMINALLY, ICE CLOUD).
!           8-17    NOVAM (NAVY OCEANIC VERTICAL AEROSOL MODEL) LAYERS.
!       NAER     NUMBER OF ACTIVE AEROSOLS.
!       EXTV     SPECTRAL EXTINCTION (NORMALIZED TO 1 AT 550 NM).
!       ABSV     SPECTRAL ABSORPTION (1-ABSV/EXTV=SCATTERING ALBEDO).
!       ASYV     SPECTRAL LEGENDRE MOMENT (DIVIDED BY 2N+1).
!       FRAC5    5 CM-1 GRID SPECTRAL INTERPOLATION FRACTION.
!       ASYVLO   ASYMMETRY FACTOR FROM PREVIOUS SPECTRAL FREQUENCY.
      INTEGER NAER
      REAL EXTV,ABSV,ASYV,FRAC5,ASYVLO
      COMMON/AER/NAER,EXTV(MAER),ABSV(MAER),ASYV(MXCMU,MAER),           &
     &  FRAC5,ASYVLO(MAER)

!     /CARD2/
!       IHAZE    BOUNDARY LAYER AEROSOL MODEL NUMBER.
!       ISEASN   SEASON NUMBER (1=SPRING-SUMMER, 2=FALL-WINTER).
!       IVULCN   VOLCANIC AEROSOL MODEL NUMBER.
!       ICSTL    COASTAL AIRMASS MODEL NUMBER.
!       ICLD     CLOUD MODEL NUMBER.
!       IVSA     VERTICAL STRUCTURE ALGORITHM (0=OFF, 1=ON).
!       VIS      SURFACE VISIBILITY (GROUND METEOROLOGICAL RANGE) [KM].
!       WSS      CURRENT WIND SPEED (M/S).
!       WHH      24-HOUR WIND SPEED (M/S).
!       RAINRT   RAIN RATE (MM/HR)
!       LSAP     LOGICAL FLAG FOR SPECTRAL AEROSOL PROFILES INPUT.
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA
      REAL VIS,WSS,WHH,RAINRT
      LOGICAL LSAP
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 &
     &  VIS,WSS,WHH,RAINRT,LSAP

!     /S_PROF/
!       S_UMIX   SCALE FACTORS FOR UNIFORMLY MIXED MOLECULAR PROFILES.
!       S_XSEC   SCALE FACTORS FOR CROSS-SECTION MOLECULAR PROFILES.
!       S_TRAC   SCALE FACTORS FOR TRACE MOLECULAR PROFILES.
!       L_UMIX   LOGICAL, .TRUE. IF S_UMIX VALUES ARE TO BE READ IN.
!       L_XSEC   LOGICAL, .TRUE. IF S_XSEC VALUES ARE TO BE READ IN.
!       L_TRAC   LOGICAL, .TRUE. IF S_TRAC VALUES ARE TO BE READ IN.
      REAL S_UMIX(4:NMOL+1),S_XSEC(NMOLX),S_TRAC(16)
      LOGICAL L_UMIX,L_XSEC,L_TRAC
      COMMON/S_PROF/S_UMIX,S_XSEC,S_TRAC,L_UMIX,L_XSEC,L_TRAC

!     /MOL_FM/
!       UMX_FM   ATM LEVEL MOLE FRACTION RATIO, RELATIVE TO CO2.
!       AUX_FM   ATM LEVEL AUX MOLECULAR GASES RELATIVE MOLE FRACTION.
      REAL UMX_FM,AUX_FM
      COMMON/MOL_FM/UMX_FM(3:NMOLXT,1:LAYDIM),AUX_FM(MMOLY,LAYDIM)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /NAMEX/
      EXTERNAL DEVCBD,XMLATM

!     FUNCTIONS:
!       RFRACT  RETURNS THE REFRACTIVITY (1 - INDEX OF REFRACTION).
      REAL RFRACT

!     LOCAL VARIABLES:
!       LEV     LEVEL INDEX.
!       IPRES   LAST ATMOSPHERIC LEVEL WITH PRESSURE EXCEEDING 10 MBAR.
!       IPRES1  FIRST ATMOSPHERIC LEVEL WITH PRESSURE LESS THAN 10 MBAR.
!       WN2     N2 MOLE FRACTION.
!       WAR     AR MOLE FRACTION.
!       WNE     HE MOLE FRACTION.
!       H2OATM  WATER PARTIAL PRESSURE [ATM AT STP].
!       RHOACK  WATER PARTIAL PRESSURE [ATM CM / KM AT STP].
!       RHOH2O  WATER DENSITY [(MOLECULES/CM2) / KM) AT 296K].
!       WH100   WATER DENSITY AT 100% RELATIVE HUMIDITY [GM / M3].
!       AIRATM  STANDARD AIR DENSITY AT STP [ATM].
!       RHSCAL  RELATIVE HUMIDITY AFTER SCALING.
!       WAUXSM  SUM OF AUXILIARY MOLECULAR GASES DENSITIES (PPMV).
      INTEGER LEV,I,K,IPRES,IPRES1
      REAL WN2,WAR,WNE,PRATIO,TRATIO,H2OATM,H2OACK,RHOH2O,CONAIR,       &
     &  CON_O2,CONH2O,CON_O3,CONCO2,CON_CO,CONCH4,CONN2O,CONNH3,        &
     &  CON_NO,CONNO2,CONSO2,WH100,AIRATM,WAUXSM,RHSCAL(LAYDIM)

!     DATA:
      INTEGER MNAER
      DATA MNAER/6/
!       MNAER    MINIMUM AEROSOL INDEX USED IN CALL TO AEREXT.

!     SET NUMBER OF ACTIVE AEROSOLS:
      IF(LNOVAM)THEN

!         INCLUDE NOVAM AEROSOLS:
          NAER=MAER
      ELSE

!         NO NOVAM AEROSOLS:
          NAER=MAER-10
      ENDIF

!     N2 MOLE FRACTION:
      IF(L_UMIX)THEN
          WN2=.781*S_UMIX(NMOL+1)
          WAR=.00934*S_UMIX(NMOL+1)
          WNE=.000018*S_UMIX(NMOL+1)
      ELSE
          WN2=.781
          WAR=.00934
          WNE=.000018
      ENDIF

!     LOAD ATMOSPHERE PROFILE INTO /MPROF/
      IPRES=0
      DO LEV=1,ML
          PM(LEV)=PPROF(LEV)
          TM(LEV)=TPROF(LEV)
          PRATIO=PM(LEV)/PZERO
          TRATIO=TZERO/TM(LEV)
          AIRATM=PRATIO*TRATIO

!         MULTIPLY AIRATM BY 1E5 TO CONVERT TO ATM CM / KM:
          DENSTY(0,LEV)=1E5*PRATIO*TRATIO

!         SCALE RAYLEIGH MOLECULAR SCATTERING BASED ON MOLE FRACTIONS
!         (RELATIVE SCATTERING CROSS SECTIONS AT 589.3NM ARE TAKEN
!         FROM TABLE 4 OF SHARDANAND AND A.D. PRASAD RAO, "ABSOLUTE
!         RAYLEIGH SCATTERING CROSS SECTIONS OF GASES AND FREONS
!         OF STRATOSPHERIC INTEREST IN THE VISIBLE AND ULTRAVIOLET
!         REGIONS," NASA TECHNICAL NOTE D-8422, MARCH 1977):
          DENSTY(6,LEV)=AIRATM*(67.6*WN2+60.5*WAR+3.5*WNE               &
     &      +(158.1*WCO2(LEV)+151.2*WCH4(LEV)+56.5*WO2(LEV)             &
     &      +13.4*WMOLXT(14,LEV)+WMOLXT(15,LEV))/1E6)/RAYDEN

!         MULTIPLY AIR DENSITY AIRATM BY (1.E5 CM/KM) / (1.E6 PPMV):
          CONAIR=AIRATM/10

!         UV OZONE (CHANGED FROM G/M**3 TO PPMV).
          DENSTY(8,LEV)=CONAIR*WO3(LEV)

!         N2 CONTINUUM (USES 296K REFERENCE TEMPERATURE):
          DENSTY(4,LEV)=WN2*(T0RAT*AIRATM)**2

!         H2O concentration in ATM at STP (H2OATM), in ATM CM / KM at
!         STP (H2OACK) and in (MOLECULES/CM2) / KM (RHOH2O) at 296K:
          H2OATM=H2O_P*WH2O(LEV)
          H2OACK=1E5*H2OATM
          RHOH2O=T0RAT*LOSCHM*H2OACK

!         SELF (5) AND FOREIGN (10) BROADENED H2O DENSITY TERMS:
          DENSTY(5,LEV)=RHOH2O*H2OATM
          DENSTY(10,LEV)=RHOH2O*(AIRATM-H2OATM)

!         TEMPERATURE DEPENDENCY OF WATER SET IN GEO.F
          DENSTY(9,LEV)=0.

!         SATURATED WATER VAPOR DENSITY [GM / M3]
          WH100=TRATIO*EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
          RHSCAL(LEV)=100*WH2O(LEV)/WH100
          IF(H2OAER.EQ.'T')THEN
              RELHUM(LEV)=RHSCAL(LEV)

!             IF MARINE AEROSOL IS USED, UPDATE OPTICAL DATA.
              IF(MARK.EQ.LEV .AND. MARIC1.GT.0)CALL MARINE(VIS,MODEL,   &
     &          RELHUM(LEV),WSS,WHH,ICSTL,EXTC,ABSC,MARIC1)
          ENDIF

!         BOUNDARY LAYER AEROSOL (0 TO 2 KM)
!         LOG WEIGHTING OF RELATIVE HUMIDITY
          IF(RELHUM(LEV).GE.99.)THEN
              DENSTY(15,LEV)=0.
          ELSEIF(ICH1.LE.7)THEN
              DENSTY(15,LEV)=LOG(100-RELHUM(LEV))*DENSTY(7,LEV)
          ELSE
              DENSTY(15,LEV)=LOG(100-RELHUM(LEV))*DENSTY(12,LEV)
          ENDIF

!         HNO3 IN ATM * CM /KM
          DENSTY(11,LEV)=CONAIR*WHNO3(LEV)

!         O2 TEMPERATURE DEPENDENCE
          CON_O2=CONAIR*WO2(LEV)
          DENSTY(63,LEV)=CON_O2*PRATIO
          DENSTY(1,LEV)=DENSTY(63,LEV)*TM(LEV)
          DENSTY(2,LEV)=DENSTY(63,LEV)*(TM(LEV)-220)**2

!         O2 1.27 MICRON CONTINUUM ABSORPTION:
          DENSTY(82,LEV)=CON_O2*AIRATM*(.7*WO2(LEV)+300000)/446000

!         O2 1.06 MICRON CONTINUUM ABSORPTION:
          DENSTY(83,LEV)=(1.E-20*LOSCHM)*T0RAT*CON_O2**2/20900

!         O2 O2 VISIBLE SPECTRAL REGION ABSORPTION:
!                       L
!                     /         2
!         tau = alpha |    p(O2)  dl    with L in km.
!                     /
!                       0
!         ABSORPTION COEFFICIENT alpha [cm-1 atm-2] IS COMPUTED FROM
!         89.5CM PATH ABSORBANCE (Anu) DATA AT 296K & 55ATM:

!                                 Anu
!         alpha = --------------------------------
!                                                2
!                 89.5cm ( 55atm 273.15K / 296K )

!         PARTIAL PRESSURE p(O2) [atm] AT TEMPERATURE T AND PRESSURE P
!         IS COMPUTED AS:

!                 chi[PPMV]
!         p(O2) = --------- P[atm] 273.15K / T[K]
!                 1E6 PPMV

!         Including a 1E5cm/km conversion factor, results in
          DENSTY(84,LEV)=(WO2(LEV)/55*PRATIO*296/TM(LEV))**2/8.95E8

!         MICROWAVE TEMPERATURE DEPENDENT RAIN
          DENSTY(61,LEV)=0.
          DENSTY(62,LEV)=0.
          IF(DENSTY(3,LEV).GT.0.)THEN
              DENSTY(3,LEV)=DENSTY(3,LEV)**.63
              DENSTY(61,LEV)=DENSTY(3,LEV)*TPROF(LEV)
              DENSTY(62,LEV)=1.
          ENDIF

!         NO2 AND SO2 DENSITIES:
          DENSTY(64,LEV)=CONAIR*WNO2(LEV)
          DENSTY(65,LEV)=CONAIR*WSO2(LEV)

!         MODTRAN USES TRUE DENSITIES [AMAGATS-CM/KM].
!         LOWTRAN USES BAND-DEPENDENT SCALED DENSITIES.
          IF(MODTRN)THEN

!             STORE WATER PARTIAL PRESSURE AT STP IN ATM CM / KM:
              DENSTY(17,LEV)=H2OACK
              DENSTY(18,LEV)=0.
              DENSTY(19,LEV)=0.
              DENSTY(20,LEV)=0.
              DENSTY(21,LEV)=0.
              DENSTY(22,LEV)=0.
              DENSTY(23,LEV)=0.
              DENSTY(24,LEV)=0.
              DENSTY(25,LEV)=0.
              DENSTY(26,LEV)=0.
              DENSTY(27,LEV)=0.
              DENSTY(28,LEV)=0.
              DENSTY(29,LEV)=0.
              DENSTY(30,LEV)=0.
              DENSTY(31,LEV)=CONAIR*WO3(LEV)
              DENSTY(32,LEV)=0.
              DENSTY(33,LEV)=0.
              DENSTY(34,LEV)=0.
              DENSTY(35,LEV)=0.

!             CO2
              DENSTY(36,LEV)=CONAIR*WCO2(LEV)
              DENSTY(37,LEV)=0.
              DENSTY(38,LEV)=0.
              DENSTY(39,LEV)=0.
              DENSTY(40,LEV)=0.
              DENSTY(41,LEV)=0.
              DENSTY(42,LEV)=0.
              DENSTY(43,LEV)=0.
              DENSTY(44,LEV)=CONAIR*WCO(LEV)
              DENSTY(45,LEV)=0.
              DENSTY(46,LEV)=CONAIR*WCH4(LEV)
              DENSTY(47,LEV)=CONAIR*WN2O(LEV)
              DENSTY(48,LEV)=0.
              DENSTY(49,LEV)=0.
              DENSTY(50,LEV)=CON_O2
              DENSTY(51,LEV)=CON_O2*PRATIO**.9353*TRATIO**(.1936)
              DENSTY(52,LEV)=CONAIR*WNH3(LEV)
              DENSTY(53,LEV)=0.
              DENSTY(54,LEV)=CONAIR*WNO(LEV)
              DENSTY(55,LEV)=CONAIR*WNO2(LEV)
              DENSTY(56,LEV)=CONAIR*WSO2(LEV)
              DENSTY(57,LEV)=0.
          ELSE

!             LOWTRAN BAND MODEL SCALED DENSITIES:
!             --- FOR H2O
              CONH2O=.1*WH2O(LEV)
              DENSTY(17,LEV)=CONH2O*PRATIO**0.9810*TRATIO**( 0.3324)
              DENSTY(18,LEV)=CONH2O*PRATIO**1.1406*TRATIO**(-2.6343)
              DENSTY(19,LEV)=CONH2O*PRATIO**0.9834*TRATIO**(-2.5294)
              DENSTY(20,LEV)=CONH2O*PRATIO**1.0443*TRATIO**(-2.4359)
              DENSTY(21,LEV)=CONH2O*PRATIO**0.9681*TRATIO**(-1.9537)
              DENSTY(22,LEV)=CONH2O*PRATIO**0.9555*TRATIO**(-1.5378)
              DENSTY(23,LEV)=CONH2O*PRATIO**0.9362*TRATIO**(-1.6338)
              DENSTY(24,LEV)=CONH2O*PRATIO**0.9233*TRATIO**(-0.9398)
              DENSTY(25,LEV)=CONH2O*PRATIO**0.8658*TRATIO**(-0.1034)
              DENSTY(26,LEV)=CONH2O*PRATIO**0.8874*TRATIO**(-0.2576)
              DENSTY(27,LEV)=CONH2O*PRATIO**0.7982*TRATIO**( 0.0588)
              DENSTY(28,LEV)=CONH2O*PRATIO**0.8088*TRATIO**( 0.2816)
              DENSTY(29,LEV)=CONH2O*PRATIO**0.6642*TRATIO**( 0.2764)
              DENSTY(30,LEV)=CONH2O*PRATIO**0.6656*TRATIO**( 0.5061)

!             --- FOR O3
              CON_O3=CONAIR*WO3(LEV)
              DENSTY(31,LEV)=CON_O3*PRATIO**0.4200*TRATIO**( 1.3909)
              DENSTY(32,LEV)=CON_O3*PRATIO**0.4221*TRATIO**( 0.7678)
              DENSTY(33,LEV)=CON_O3*PRATIO**0.3739*TRATIO**( 0.1225)
              DENSTY(34,LEV)=CON_O3*PRATIO**0.1770*TRATIO**( 0.9827)
              DENSTY(35,LEV)=CON_O3*PRATIO**0.3921*TRATIO**( 0.1942)

!             --- FOR CO2
              CONCO2=CONAIR*WCO2(LEV)
              DENSTY(36,LEV)=CONCO2*PRATIO**0.6705*TRATIO**(-2.2560)
              DENSTY(37,LEV)=CONCO2*PRATIO**0.7038*TRATIO**(-5.0768)
              DENSTY(38,LEV)=CONCO2*PRATIO**0.7258*TRATIO**(-1.6740)
              DENSTY(39,LEV)=CONCO2*PRATIO**0.6982*TRATIO**(-1.8107)
              DENSTY(40,LEV)=CONCO2*PRATIO**0.8867*TRATIO**(-0.5327)
              DENSTY(41,LEV)=CONCO2*PRATIO**0.7883*TRATIO**(-1.3244)
              DENSTY(42,LEV)=CONCO2*PRATIO**0.6899*TRATIO**(-0.8152)
              DENSTY(43,LEV)=CONCO2*PRATIO**0.6035*TRATIO**( 0.6026)

!             --- FOR CO
              CON_CO=CONAIR*WCO(LEV)
              DENSTY(44,LEV)=CON_CO*PRATIO**0.7589*TRATIO**( 0.6911)
              DENSTY(45,LEV)=CON_CO*PRATIO**0.9267*TRATIO**( 0.1716)

!             --- FOR CH4
              CONCH4=CONAIR*WCH4(LEV)
              DENSTY(46,LEV)=CONCH4*PRATIO**0.7139*TRATIO**(-0.4185)

!             --- FOR N2O
              CONN2O=CONAIR*WN2O(LEV)
              DENSTY(47,LEV)=CONN2O*PRATIO**0.3783*TRATIO**( 0.9399)
              DENSTY(48,LEV)=CONN2O*PRATIO**0.7203*TRATIO**(-0.1836)
              DENSTY(49,LEV)=CONN2O*PRATIO**0.7764*TRATIO**( 1.1931)

!             --- FOR O2
              DENSTY(50,LEV)=CON_O2*PRATIO**1.1879*TRATIO**( 2.9738)
              DENSTY(51,LEV)=CON_O2*PRATIO**0.9353*TRATIO**( 0.1936)

!             --- FOR NH3
              CONNH3=CONAIR*WNH3(LEV)
              DENSTY(52,LEV)=CONNH3*PRATIO**0.8023*TRATIO**(-0.9111)
              DENSTY(53,LEV)=CONNH3*PRATIO**0.6968*TRATIO**( 0.3377)

!             --- FOR NO
              CON_NO=CONAIR*WNO(LEV)
              DENSTY(54,LEV)=CON_NO*PRATIO**0.5265*TRATIO**(-0.4702)

!             --- FOR NO2
              CONNO2=CONAIR*WNO2(LEV)
              DENSTY(55,LEV)=CONNO2*PRATIO**0.3956*TRATIO**(-0.0545)

!             --- FOR SO2
              CONSO2=CONAIR*WSO2(LEV)
              DENSTY(56,LEV)=CONSO2*PRATIO**0.2943*TRATIO**( 1.2316)
              DENSTY(57,LEV)=CONSO2*PRATIO**0.2135*TRATIO**( 0.0733)
          ENDIF

!         CROSS-SECTION SPECIES:
          DO K=1,NMOLX
              DENSTY(MEXT+K,LEV)=CONAIR*WMOLXT(K,LEV)
          ENDDO

!         H2-H2 and H2-He Dimers:
!         Multiply density by additional H2 density, and
!         divide (atm cm / km)^2 by 1E5 cm / km to get atm^2 cm / km.
          DENSTY(MEXT+15,LEV)                                           &
     &      =DENSTY(MEXT+15,LEV)*DENSTY(MEXT+14,LEV)/100000
          DENSTY(MEXT+14,LEV)=DENSTY(MEXT+14,LEV)**2/100000

!         MOLE FRACTIONS RELATIVE TO CO2:
          IF(WCO2(LEV).GT.0.)THEN
!             UMX_F(3,*) IS N2: INCLUDE FOREIGN BROADENED
!             AIR AND CONVERSION FROM KM-1 TO CM-1.
              UMX_FM( 3,LEV)=1.E-5*T0RAT*2*AIRATM*WN2/WCO2(LEV)
              UMX_FM( 4,LEV)=WN2O(LEV)/WCO2(LEV)
              UMX_FM( 5,LEV)=WCO(LEV)/WCO2(LEV)
              UMX_FM( 6,LEV)=WCH4(LEV)/WCO2(LEV)
              UMX_FM( 7,LEV)=WO2(LEV)/WCO2(LEV)
              UMX_FM( 8,LEV)=WNO(LEV)/WCO2(LEV)
              UMX_FM( 9,LEV)=WSO2(LEV)/WCO2(LEV)
              UMX_FM(10,LEV)=WNO2(LEV)/WCO2(LEV)
              UMX_FM(11,LEV)=WNH3(LEV)/WCO2(LEV)
              UMX_FM(12,LEV)=WHNO3(LEV)/WCO2(LEV)
              DO K=1,NMOLX
                  UMX_FM(NMOL+K,LEV)=WMOLXT(K,LEV)/WCO2(LEV)
              ENDDO
          ELSE
              DO K=3,NMOLXT
                  UMX_FM(K,LEV)=0.
              ENDDO
          ENDIF

!         Y-SPECIES:
          WAUXSM=0.
          DO K=1,NMOLY
              DENSTY(MEXTX+K,LEV)=CONAIR*WMOLYT(K,LEV)
              WAUXSM=WAUXSM+WMOLYT(K,LEV)
          ENDDO
          IF(WAUXSM.GT.0.)THEN

!             NORMALIZE THE AUXILIARY SPECIES MOLE FRACTIONS:
              DO K=1,NMOLY
                  AUX_FM(K,LEV)=WMOLYT(K,LEV)/WAUXSM
              ENDDO
          ELSEIF(LEV.GT.1)THEN

!             USE PREVIOUS LEVEL MOLE FRACTION RATIOS IF SUM IS ZERO:
              DO K=1,NMOLY
                  AUX_FM(K,LEV)=AUX_FM(K,LEV-1)
              ENDDO
          ELSE
              DO K=1,NMOLY
                  AUX_FM(K,LEV)=0.
              ENDDO
          ENDIF

!         HERZBERG OXYGEN CONTINUUM PRESSURE DEPENDENCE
!         CALCULATION, SHARDANAND 1977 AND YOSHINO ET AL 1988
          DENSTY(58,LEV)=(1.+.83*AIRATM)*CON_O2
          DENSTY(59,LEV)=0.
          DENSTY(60,LEV)=0.

!         UV-VIS NO2 AND SO2.
          DENSTY(64,LEV)=CONAIR*WNO2(LEV)
          DENSTY(65,LEV)=CONAIR*WSO2(LEV)

!         RFNDX    REFRACTIVITY (1-INDEX OF REFRACTION):
          RFNDX(LEV)=RFRACT(VRFRAC,PRATIO,TM(LEV),WH2O(LEV),WCO2(LEV))

!         PRESSURE PRINT FLAG:
          IF(PM(LEV).GT.9.9995)IPRES=LEV
      ENDDO
      IPRES1=IPRES+1
      IF(NPR.GE.1 .OR. LJMASS)RETURN
      WRITE(IPR,'(//A,//(3A))')' ATMOSPHERIC PROFILES',                 &
     &  '   I     Z        P       T        N2     ',                   &
     &  '  CNTMSLF   CNTMFRN MOL SCAT     N-1   ',                      &
     &  '  O3 (UV)   O2 (UV)   WAT DROP  ICE PART  RAIN RATE',          &
     &  '        (KM)     (MB)    (K)              ',                   &
     &  ' (  MOL/CM2 KM  )      (-)       (-)   ',                      &
     &  ' (  ATM CM/KM  )       (GM/M3)   (GM/M3)   (MM/HR)'
      IF(IPRES.GT.0)WRITE(IPR,'((I4,0P,F9.4,F10.3,F7.2,1X,1P,7E10.3,    &
     &  0P,3F10.3))')(I,ZM(I),PM(I),TM(I),DENSTY(4,I),DENSTY(5,I),      &
     &  DENSTY(10,I),DENSTY(6,I),RFNDX(I),DENSTY(8,I),DENSTY(58,I),     &
     &  (DENSTY(K,I),K=66,67),DENSTY(3,I)**1.5873,I=1,IPRES)
      IF(IPRES.LT.ML)WRITE(IPR,'((I4,0P,F9.4,F10.7,F7.2,1X,1P,7E10.3,   &
     &  0P,3F10.3))')(I,ZM(I),PM(I),TM(I),DENSTY(4,I),DENSTY(5,I),      &
     &  DENSTY(10,I),DENSTY(6,I),RFNDX(I),DENSTY(8,I),DENSTY(58,I),     &
     &  (DENSTY(K,I),K=66,67),DENSTY(3,I)**1.5873,I=IPRES1,ML)

!     CALCULATE 550 NM CLOUD EXTINCTION COEFFICIENTS (KM-1 M3/GM).
      CALL AEREXT(V550NM,MNAER)
      IF(H2OAER.EQ.'T')THEN
          WRITE(IPR,'(//A,//(3A))')' ATMOSPHERIC PROFILES',             &
     &      '   I     Z        P       T    ',                          &
     &      ' AEROSOL 1 AEROSOL 2 AEROSOL 3 AEROSOL 4 ',                &
     &      ' AER1*RH     RH (%)   CIRRUS   WAT DROP  ICE PART',        &
     &      '        (KM)     (MB)    (K)   ',                          &
     &      '(       550nm EXTINCTION [KM-1]        ) ',                &
     &      ' (AFTER H2O SCALING)   (-)     (550nm VIS [KM-1])'
          IF(IPRES.GT.0)WRITE(IPR,'((I4,0P,F9.4,F10.3,F7.2,1X,1P,       &
     &      6E10.3,0P,3F10.5))')(I,ZM(I),PM(I),TM(I),DENSTY(7,I),       &
     &      (DENSTY(K,I),K=12,15),RELHUM(I),DENSTY(16,I),               &
     &      EXTV(6)*DENSTY(66,I),EXTV(7)*DENSTY(67,I),I=1,IPRES)
          IF(IPRES.LT.ML)WRITE(IPR,'((I4,0P,F9.4,F10.7,F7.2,1X,1P,      &
     &      6E10.3,0P,3F10.5))')(I,ZM(I),PM(I),TM(I),DENSTY(7,I),       &
     &      (DENSTY(K,I),K=12,15),RELHUM(I),DENSTY(16,I),               &
     &      EXTV(6)*DENSTY(66,I),EXTV(7)*DENSTY(67,I),I=IPRES1,ML)
      ELSE
          WRITE(IPR,'(//A,//(3A))')' ATMOSPHERIC PROFILES',             &
     &      '   I     Z        P       T     AEROSOL 1 AEROSOL 2',      &
     &      ' AEROSOL 3 AEROSOL 4  AER1*RH     RH (%)    RH (%)',       &
     &      '   CIRRUS   WAT DROP  ICE PART',                           &
     &      '        (KM)     (MB)    (K)    (      550nm EXTINC',      &
     &      'TION [KM-1]      )  (BEFORE H2O SCALING)   (AFTER)',       &
     &      '    (-)     (550nm VIS [KM-1])'
          IF(IPRES.GT.0)WRITE(IPR,'((I4,0P,F9.4,F10.3,F7.2,1X,1P,       &
     &      7E10.3,0P,3F10.5))')(I,ZM(I),PM(I),TM(I),DENSTY(7,I),       &
     &      (DENSTY(K,I),K=12,15),RELHUM(I),RHSCAL(I),DENSTY(16,I),     &
     &      EXTV(6)*DENSTY(66,I),EXTV(7)*DENSTY(67,I),I=1,IPRES)
          IF(IPRES.LT.ML)WRITE(IPR,'((I4,0P,F9.4,F10.7,F7.2,1X,1P,      &
     &      7E10.3,0P,3F10.5))')(I,ZM(I),PM(I),TM(I),DENSTY(7,I),       &
     &      (DENSTY(K,I),K=12,15),RELHUM(I),RHSCAL(I),DENSTY(16,I),     &
     &      EXTV(6)*DENSTY(66,I),EXTV(7)*DENSTY(67,I),I=IPRES1,ML)
      ENDIF
      WRITE(IPR,'(//A,//3A)')                                           &
     &  ' ATMOSPHERIC PROFILES (AFTER COLUMN SCALING)',                 &
     &   '   I      Z        P       H2O      O3       ',               &
     &  'CO2      CO       CH4      N2O      O2       ',                &
     &  'NH3      NO       NO2      SO2      HNO3'
      IF(MODTRN)THEN
          WRITE(IPR,'(3A)')'         (KM)     (MB)  (         ',        &
     &      '                                    ATM CM / KM',          &
     &      '                                                )'
      ELSE
          WRITE(IPR,'(3A)')'         (KM)     (MB)  (         ',        &
     &      '                   LOWTRAN SCALED DENSITIES FOR',          &
     &      ' ONE OF THE BANDS                               )'
      ENDIF
      IF(IPRES.GT.0)WRITE(IPR,'((I4,0P,F9.4,F10.3,1X,1P,12E9.2))')(I,   &
     &  ZM(I),PM(I),DENSTY(17,I),DENSTY(31,I),DENSTY(36,I),DENSTY(44,I),&
     &  DENSTY(46,I),DENSTY(47,I),DENSTY(50,I),DENSTY(52,I),            &
     &  DENSTY(54,I),DENSTY(55,I),DENSTY(56,I),DENSTY(11,I),I=1,IPRES)
      IF(IPRES.LT.ML)WRITE(IPR,'((I4,0P,F9.4,F10.7,1X,1P,12E9.2))')(I,  &
     &  ZM(I),PM(I),DENSTY(17,I),DENSTY(31,I),DENSTY(36,I),DENSTY(44,I),&
     &  DENSTY(46,I),DENSTY(47,I),DENSTY(50,I),DENSTY(52,I),            &
     &  DENSTY(54,I),DENSTY(55,I),DENSTY(56,I),DENSTY(11,I),I=IPRES1,ML)

!M3D  MOD3D PROFILES:
!M3D  WRITE(69,'(2A)')'LAY Z[KM] P[ATM] T[K] H2O[GM/M3] CO2[GM/M3]',
!M3D 1  ' O3[GM/M3] AER[KM-1 @ 550nm] CLOUD[GM/M3] RAIN [(MM/HR)**.63]'
!M3D  WRITE(69,'((0P,F7.3,1P,E10.3,0P,F7.2,1P,4E11.4,0P,2F4.1))')
!M3D 1  (ZM(I),PM(I)/1013.25,TM(I),DENSTY(17,I)*18.015/2241.383,
!M3D 2  DENSTY(36,I)*44.010/2241.383,DENSTY(31,I)*47.998/2241.383,
!M3D 3  DENSTY(7,I)+DENSTY(12,I)+DENSTY(13,I)+DENSTY(14,I),
!M3D 4  DENSTY(66,I),DENSTY(3,I),I=1,ML)

      WRITE(IPR,'(//A,//A,13(1X,A8:),/(14X,13(1X,A8:)))')               &
     &  ' ATMOSPHERIC PROFILES','   I      Z   ',(CNAMEX(K),K=1,NMOLX)
      WRITE(IPR,'((9X,A,55X,A,50X,A))')'(KM)  (','ATM CM/KM',')',       &
     &                                 '      (  ATM^2 CM/KM  )'
      DO I=1,ML
          WRITE(IPR,'(I4,0P,F9.4,1X,1P,13E9.2:,/(14X,13E9.2:))')        &
     &      I,ZM(I),(DENSTY(K,I),K=MEXT+1,MEXTX)
      ENDDO
      RETURN
      END