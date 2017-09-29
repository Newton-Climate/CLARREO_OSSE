      SUBROUTINE STDMDL(ICH1,MARIC1,MARK,H2OAER,LNOVAM)

!     THIS SUBROUTINE LOADS ONE OF THE 6 STANDARD ATMOSPHERIC
!     PROFILES INTO COMMON /MODEL/ AND CALCULATES THE
!     DENSITIES OF THE VARIOUS ABSORBING GASES AND AEROSOLS.

!     INPUT ARGUMENTS:
!       MARIC1   MARINE AEROSOL REGION NUMBER (=0 IF NOT USED).
!       MARK     RELATIVE HUMIDITY LAYER INDEX FOR MARINE AEROSOL.
!       H2OAER   FLAG, TRUE IF AEROSOL DEFAULT PROFILE PROPERTIES
!                ARE UPDATED AFTER WATER COLUMN SCALING.
!       LNOVAM   LOGICAL FLAG, .TRUE. IF NOVAM AEROSOLS ARE USED.
      INTEGER ICH1,MARIC1,MARK
      CHARACTER*1 H2OAER
      LOGICAL LNOVAM

!     PARAMETERS:
!       T0RAT    RATIO OF FASCODE STANDARD TEMPERATURE (296K)
!                TO MODTRAN STANDARD TEMPERATURE (273.15K)
      REAL T0RAT
      PARAMETER(T0RAT=296./273.15)
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:
      CHARACTER*8 CNAMEX
      COMMON/NAMEX/CNAMEX(NMOLX)
      REAL WMOLXT
      COMMON/MDATAX/WMOLXT(NMOLX,LAYDIM)
      REAL DNSTYX
      COMMON/MODELX/DNSTYX(NMOLX,LAYDIM)
      INCLUDE 'BASE.h'

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
      INCLUDE 'IFIL.h'

!     /CARD1/
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,MODTRN

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

      REAL P,T,WH,WCO2,WO,WN2O,WCO,WCH4,WO2
      COMMON/MDATA/P(LAYDIM),T(LAYDIM),WH(LAYDIM),WCO2(LAYDIM),         &
     &  WO(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),WCH4(LAYDIM),WO2(LAYDIM)
      REAL WNO,WSO2,WNO2,WNH3,WHNO3
      COMMON/MDATA1/WNO(LAYDIM),WSO2(LAYDIM),WNO2(LAYDIM),              &
     &  WNH3(LAYDIM),WHNO3(LAYDIM)

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT

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
      REAL RE,ZMAX
      INTEGER IPATH
      COMMON/PARMTR/RE,ZMAX,IPATH

!     /CARD2/
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA
      REAL VIS,WSS,WHH,RAINRT
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 &
     &  VIS,WSS,WHH,RAINRT

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /NAMEX/
      EXTERNAL DEVCBD,XMLATM

!     FUNCTIONS:
!       RFRACT  RETURNS THE REFRACTIVITY (1 - INDEX OF REFRACTION).
      REAL RFRACT

!     LOCAL VARIABLES:
!       WH100   WATER DENSITY AT 100% RELATIVE HUMIDITY:
!       RHSCAL  RELATIVE HUMIDITY AFTER SCALING.
      INTEGER I,K
      REAL PRATIO,TRATIO,RHOAIR,RHOH2O,CONAIR,CONO2,                    &
     &  CONH2O,CONO3,CONCO2,CONCO,CONCH4,CONN2O,CONNH3,                 &
     &  CONNO,CONNO2,CONSO2,WH100,RHSCAL(LAYDIM)

!     DATA:
!       MNAER    MINIMUM AEROSOL INDEX USED IN CALL TO AEREXT.
!       V550NM   FREQUENCY EQUIVALENT OF 550 NM [CM-1].
!       XLOSCH   LOSCHMIDT'S NUMBER,MOLECULES CM-2 / KM.
!       CON      CONVERSION FROM GM H2O / M3 TO MOLECULES / (CM2 KM).
!       RV       GAS CONSTANT FOR WATER IN MB/(GM M-3 K).
      INTEGER MNAER
      REAL V550NM,PZERO,TZERO,XLOSCH,CON
      DATA MNAER/6/,V550NM/18181.818/,PZERO/1013.25/,TZERO/273.15/,     &
     &  XLOSCH/2.686754E24/,CON/3.3429E21/

!     /PARMTR/ VARIABLES.
      RE=REE
      ZMAX=ZM(ML)

!     SET NUMBER OF ACTIVE AEROSOLS:
      IF(LNOVAM)THEN

!         INCLUDE NOVAM AEROSOLS:
          NAER=MAER
      ELSE

!         NO NOVAM AEROSOLS:
          NAER=MAER-10
      ENDIF

!     LOAD ATMOSPHERE PROFILE INTO /MODEL/
      DO 20 I=1,ML
          PM(I)=P(I)
          TM(I)=T(I)
          PRATIO=PM(I)/PZERO
          TRATIO=TZERO/TM(I)
          RHOAIR=PRATIO*TRATIO

!         MOLECULAR SCATTERING
          DENSTY(6,I)=RHOAIR

!         MULTIPLY THE AIR DENSITY (RHOAIR) BY (1.E5 CM/KM)/(1.E6 PPMV):
          CONAIR=.1*RHOAIR

!         UV OZONE (CHANGED FROM G/M**3 TO PPMV).
          DENSTY(8,I)=CONAIR*WO(I)

!         N2 CONTINUUM (FASCODE APPROACH)
          DENSTY(4,I)=.781*(T0RAT*RHOAIR)**2

!         H2O CONTINUUM
          RHOH2O=CON*WH(I)
          DENSTY(5,I)=T0RAT*RHOH2O/XLOSCH*RHOH2O

!         H2O FOREIGN BROADENED
          DENSTY(10,I)=T0RAT*RHOH2O*(RHOAIR-RHOH2O/XLOSCH)

!         TEMPERATURE DEPENDENCY OF WATER SET IN GEO.F
          DENSTY(9,I)=0.

!         SATURATED WATER VAPOR DENSITY [GM / M3]
          WH100=TRATIO*EXP(18.9766-(14.9595+2.43882*TRATIO)*TRATIO)
          RHSCAL(I)=100.*WH(I)/WH100
          IF(H2OAER.EQ.'T')THEN
              RELHUM(I)=RHSCAL(I)

!             IF MARINE AEROSOL IS USED, UPDATE OPTICAL DATA.
              IF(MARK.EQ.I .AND. MARIC1.GT.0)CALL MARINE(VIS,MODEL,     &
     &          RELHUM(I),WSS,WHH,ICSTL,EXTC,ABSC,MARIC1)
          ENDIF

!         BOUNDARY LAYER AEROSOL (0 TO 2 KM)
!         LOG WEIGHTING OF RELATIVE HUMIDITY
          IF(RELHUM(I).GE.99.)THEN
              DENSTY(15,I)=0.
          ELSEIF(ICH1.LE.7)THEN
              DENSTY(15,I)=LOG(100.-RELHUM(I))*DENSTY(7,I)
          ELSE
              DENSTY(15,I)=LOG(100.-RELHUM(I))*DENSTY(12,I)
          ENDIF

!         HNO3 IN ATM * CM /KM
          DENSTY(11,I)=RHOAIR*WHNO3(I)*1.E-4

!         O2 TEMPERATURE DEPENDENCE
          CONO2=CONAIR*WO2(I)
          DENSTY(63,I)=CONO2*PRATIO
          DENSTY(1,I)=DENSTY(63,I)*TM(I)
          DENSTY(2,I)=DENSTY(63,I)*(TM(I)-220)**2

!         O2 1.27 MICRON CONTINUUM ABSORPTION:
          DENSTY(82,I)=CONO2*RHOAIR*(.7*WO2(I)+300000)/446000

!         O2 1.06 MICRON CONTINUUM ABSORPTION:
          DENSTY(83,I)=(1.E-20*XLOSCH)*T0RAT*CONO2**2/.209E10

!         MICROWAVE TEMPERATURE DEPENDENT RAIN
          DENSTY(61,I)=0.
          DENSTY(62,I)=0.
          IF(DENSTY(3,I).GT.0.)THEN
              DENSTY(3,I)=DENSTY(3,I)**.63
              DENSTY(61,I)=DENSTY(3,I)*T(I)
              DENSTY(62,I)=1.
          ENDIF

!         NO2 AND SO2 DENSITIES:
          DENSTY(64,I)=CONAIR*WNO2(I)
          DENSTY(65,I)=CONAIR*WSO2(I)

!         MODTRAN USES TRUE DENSITIES [AMAGATS-CM/KM].
!         LOWTRAN USES BAND-DEPENDENT SCALED DENSITIES.
          IF(MODTRN)THEN

!             FOR WATER, THE GM / M3 VALUE, WH(I), IS CONVERTED
!             TO AMAGAT-CM BY MULTIPLYING BY 124.42
              DENSTY(17,I)=WH(I)*124.42
              DENSTY(18,I)=0.
              DENSTY(19,I)=0.
              DENSTY(20,I)=0.
              DENSTY(21,I)=0.
              DENSTY(22,I)=0.
              DENSTY(23,I)=0.
              DENSTY(24,I)=0.
              DENSTY(25,I)=0.
              DENSTY(26,I)=0.
              DENSTY(27,I)=0.
              DENSTY(28,I)=0.
              DENSTY(29,I)=0.
              DENSTY(30,I)=0.
              DENSTY(31,I)=CONAIR*WO(I)
              DENSTY(32,I)=0.
              DENSTY(33,I)=0.
              DENSTY(34,I)=0.
              DENSTY(35,I)=0.

!             CO2
              DENSTY(36,I)=CONAIR*WCO2(I)
              DENSTY(37,I)=0.
              DENSTY(38,I)=0.
              DENSTY(39,I)=0.
              DENSTY(40,I)=0.
              DENSTY(41,I)=0.
              DENSTY(42,I)=0.
              DENSTY(43,I)=0.
              DENSTY(44,I)=CONAIR*WCO(I)
              DENSTY(45,I)=0.
              DENSTY(46,I)=CONAIR*WCH4(I)
              DENSTY(47,I)=CONAIR*WN2O(I)
              DENSTY(48,I)=0.
              DENSTY(49,I)=0.
              DENSTY(50,I)=CONO2
              DENSTY(51,I)=CONO2*PRATIO**.9353*TRATIO**(.1936)
              DENSTY(52,I)=CONAIR*WNH3(I)
              DENSTY(53,I)=0.
              DENSTY(54,I)=CONAIR*WNO(I)
              DENSTY(55,I)=CONAIR*WNO2(I)
              DENSTY(56,I)=CONAIR*WSO2(I)
              DENSTY(57,I)=0.
          ELSE

!             LOWTRAN BAND MODEL SCALED DENSITIES:
!             --- FOR H2O
              CONH2O=.1*WH(I)
              DENSTY(17,I)=CONH2O*PRATIO**0.9810*TRATIO**( 0.3324)
              DENSTY(18,I)=CONH2O*PRATIO**1.1406*TRATIO**(-2.6343)
              DENSTY(19,I)=CONH2O*PRATIO**0.9834*TRATIO**(-2.5294)
              DENSTY(20,I)=CONH2O*PRATIO**1.0443*TRATIO**(-2.4359)
              DENSTY(21,I)=CONH2O*PRATIO**0.9681*TRATIO**(-1.9537)
              DENSTY(22,I)=CONH2O*PRATIO**0.9555*TRATIO**(-1.5378)
              DENSTY(23,I)=CONH2O*PRATIO**0.9362*TRATIO**(-1.6338)
              DENSTY(24,I)=CONH2O*PRATIO**0.9233*TRATIO**(-0.9398)
              DENSTY(25,I)=CONH2O*PRATIO**0.8658*TRATIO**(-0.1034)
              DENSTY(26,I)=CONH2O*PRATIO**0.8874*TRATIO**(-0.2576)
              DENSTY(27,I)=CONH2O*PRATIO**0.7982*TRATIO**( 0.0588)
              DENSTY(28,I)=CONH2O*PRATIO**0.8088*TRATIO**( 0.2816)
              DENSTY(29,I)=CONH2O*PRATIO**0.6642*TRATIO**( 0.2764)
              DENSTY(30,I)=CONH2O*PRATIO**0.6656*TRATIO**( 0.5061)

!             --- FOR O3
              CONO3=CONAIR*WO(I)
              DENSTY(31,I)=CONO3 *PRATIO**0.4200*TRATIO**( 1.3909)
              DENSTY(32,I)=CONO3 *PRATIO**0.4221*TRATIO**( 0.7678)
              DENSTY(33,I)=CONO3 *PRATIO**0.3739*TRATIO**( 0.1225)
              DENSTY(34,I)=CONO3 *PRATIO**0.1770*TRATIO**( 0.9827)
              DENSTY(35,I)=CONO3 *PRATIO**0.3921*TRATIO**( 0.1942)

!             --- FOR CO2
              CONCO2=CONAIR*WCO2(I)
              DENSTY(36,I)=CONCO2*PRATIO**0.6705*TRATIO**(-2.2560)
              DENSTY(37,I)=CONCO2*PRATIO**0.7038*TRATIO**(-5.0768)
              DENSTY(38,I)=CONCO2*PRATIO**0.7258*TRATIO**(-1.6740)
              DENSTY(39,I)=CONCO2*PRATIO**0.6982*TRATIO**(-1.8107)
              DENSTY(40,I)=CONCO2*PRATIO**0.8867*TRATIO**(-0.5327)
              DENSTY(41,I)=CONCO2*PRATIO**0.7883*TRATIO**(-1.3244)
              DENSTY(42,I)=CONCO2*PRATIO**0.6899*TRATIO**(-0.8152)
              DENSTY(43,I)=CONCO2*PRATIO**0.6035*TRATIO**( 0.6026)

!             --- FOR CO
              CONCO=CONAIR*WCO (I)
              DENSTY(44,I)=CONCO *PRATIO**0.7589*TRATIO**( 0.6911)
              DENSTY(45,I)=CONCO *PRATIO**0.9267*TRATIO**( 0.1716)

!             --- FOR CH4
              CONCH4=CONAIR*WCH4(I)
              DENSTY(46,I)=CONCH4*PRATIO**0.7139*TRATIO**(-0.4185)

!             --- FOR N2O
              CONN2O=CONAIR*WN2O(I)
              DENSTY(47,I)=CONN2O*PRATIO**0.3783*TRATIO**( 0.9399)
              DENSTY(48,I)=CONN2O*PRATIO**0.7203*TRATIO**(-0.1836)
              DENSTY(49,I)=CONN2O*PRATIO**0.7764*TRATIO**( 1.1931)

!             --- FOR O2
              DENSTY(50,I)=CONO2 *PRATIO**1.1879*TRATIO**( 2.9738)
              DENSTY(51,I)=CONO2 *PRATIO**0.9353*TRATIO**( 0.1936)

!             --- FOR NH3
              CONNH3=CONAIR*WNH3(I)
              DENSTY(52,I)=CONNH3*PRATIO**0.8023*TRATIO**(-0.9111)
              DENSTY(53,I)=CONNH3*PRATIO**0.6968*TRATIO**( 0.3377)

!             --- FOR NO
              CONNO=CONAIR*WNO(I)
              DENSTY(54,I)=CONNO *PRATIO**0.5265*TRATIO**(-0.4702)

!             --- FOR NO2
              CONNO2=CONAIR*WNO2(I)
              DENSTY(55,I)=CONNO2*PRATIO**0.3956*TRATIO**(-0.0545)

!             --- FOR SO2
              CONSO2=CONAIR*WSO2(I)
              DENSTY(56,I)=CONSO2*PRATIO**0.2943*TRATIO**( 1.2316)
              DENSTY(57,I)=CONSO2*PRATIO**0.2135*TRATIO**( 0.0733)
          ENDIF

!         CROSS-SECTION SPECIES:
          DO 10 K=1,NMOLX
              DNSTYX(K,I)=CONAIR*WMOLXT(K,I)
   10     CONTINUE

!         HERZBERG OXYGEN CONTINUUM PRESSURE DEPENDENCE
!         CALCULATION, SHARDANAND 1977 AND YOSHINO ET AL 1988
          DENSTY(58,I)=(1.+.83*RHOAIR)*CONO2
          DENSTY(59,I)=0.
          DENSTY(60,I)=0.

!         UV-VIS NO2 AND SO2.
          DENSTY(64,I)=CONAIR*WNO2(I)
          DENSTY(65,I)=CONAIR*WSO2(I)

!         RFNDX    REFRACTIVITY (1-INDEX OF REFRACTION):
          RFNDX(I)=RFRACT(.5*(IV1+IV2),PRATIO,TM(I),WH(I),WCO2(I))
   20 CONTINUE
      IF(NPR.GE.1 .OR. LJMASS)RETURN
      WRITE(IPR,'(//A,//(3A))')'1 ATMOSPHERIC PROFILES',                &
     &  '   I     Z       P       T        N2     ',                    &
     &  '  CNTMSLF   CNTMFRN MOL SCAT     N-1   ',                      &
     &  '  O3 (UV)   O2 (UV)   WAT DROP  ICE PART  RAIN RATE',          &
     &  '        (KM)    (MB)    (K)              ',                    &
     &  ' (  MOL/CM2 KM  )      (-)       (-)   ',                      &
     &  ' (  ATM CM/KM  )       (GM/M3)   (GM/M3)   (MM/HR)'
      WRITE(IPR,'(/(I4,0P,F9.4,F9.3,F7.1,1X,1P,7E10.3,0P,3F10.3))')     &
     &  (I,ZM(I),PM(I),TM(I),DENSTY(4,I),DENSTY(5,I),DENSTY(10,I),      &
     &  DENSTY(6,I),RFNDX(I),DENSTY(8,I),DENSTY(58,I),                  &
     &  (DENSTY(K,I),K=66,67),DENSTY(3,I)**1.5873,I=1,ML)

!     CALCULATE 550 NM CLOUD EXTINCTION COEFFICIENTS (KM-1 M3/GM).
      CALL AEREXT(V550NM,MNAER)
      IF(H2OAER.EQ.'T')THEN
          WRITE(IPR,'(//A,//(3A))')'1 ATMOSPHERIC PROFILES',            &
     &      '   I     Z       P       T    ',                           &
     &      ' AEROSOL 1 AEROSOL 2 AEROSOL 3 AEROSOL 4 ',                &
     &      ' AER1*RH     RH (%)   CIRRUS   WAT DROP  ICE PART',        &
     &      '        (KM)    (MB)    (K)   ',                           &
     &      '(       550nm EXTINCTION [KM-1]        ) ',                &
     &      ' (AFTER H2O SCALING)   (-)     (550nm VIS [KM-1])'
          WRITE(IPR,'(/(I4,0P,F9.4,F9.3,F7.1,1X,1P,6E10.3,0P,3F10.5))') &
     &      (I,ZM(I),PM(I),TM(I),DENSTY(7,I),                           &
     &      (DENSTY(K,I),K=12,15),RELHUM(I),DENSTY(16,I),               &
     &      EXTV(6)*DENSTY(66,I),EXTV(7)*DENSTY(67,I),I=1,ML)
      ELSE
          WRITE(IPR,'(//A,//(3A))')'1 ATMOSPHERIC PROFILES',            &
     &      '   I     Z       P       T     AEROSOL 1 AEROSOL 2',       &
     &      ' AEROSOL 3 AEROSOL 4  AER1*RH     RH (%)    RH (%)',       &
     &      '   CIRRUS   WAT DROP  ICE PART',                           &
     &      '        (KM)    (MB)    (K)    (      550nm EXTINC',       &
     &      'TION [KM-1]      )  (BEFORE H2O SCALING)   (AFTER)',       &
     &      '    (-)     (550nm VIS [KM-1])'
          WRITE(IPR,'(/(I4,0P,F9.4,F9.3,F7.1,1X,1P,7E10.3,0P,3F10.5))') &
     &      (I,ZM(I),PM(I),TM(I),DENSTY(7,I),                           &
     &      (DENSTY(K,I),K=12,15),RELHUM(I),RHSCAL(I),DENSTY(16,I),     &
     &      EXTV(6)*DENSTY(66,I),EXTV(7)*DENSTY(67,I),I=1,ML)
      ENDIF
      WRITE(IPR,'(//A,//3A)')                                           &
     &  '1 ATMOSPHERIC PROFILES (AFTER COLUMN SCALING)',                &
     &   '   I      Z       P       H2O      O3       ',                &
     &  'CO2      CO       CH4      N2O      O2       ',                &
     &  'NH3      NO       NO2      SO2      HNO3'
      IF(MODTRN)THEN
          WRITE(IPR,'(3A)')'         (KM)    (MB)  (         ',         &
     &      '                                    ATM CM / KM',          &
     &      '                                                )'
      ELSE
          WRITE(IPR,'(3A)')'         (KM)    (MB)  (         ',         &
     &      '                   LOWTRAN SCALED DENSITIES FOR',          &
     &      ' ONE OF THE BANDS                               )'
      ENDIF
      WRITE(IPR,'((I4,0P,F9.4,F9.3,1X,1P,12E9.2))')(I,ZM(I),PM(I),      &
     &  DENSTY(17,I),DENSTY(31,I),DENSTY(36,I),DENSTY(44,I),            &
     &  DENSTY(46,I),DENSTY(47,I),DENSTY(50,I),DENSTY(52,I),            &
     &  DENSTY(54,I),DENSTY(55,I),DENSTY(56,I),DENSTY(11,I),I=1,ML)
!69   write(69,'((0P,F10.5,20X,1P,3E10.4,A,/8E10.4,/E10.4))')(ZM(I),    &
!69  &  DENSTY(17,I)*2.686754E14,                                       &
!69  &  DENSTY(36,I)*2.686754E14,                                       &
!69  &  DENSTY(31,I)*2.686754E14,                                       &
!69  &  '  BBBBBBBBBBBBBBBBBB',                                         &
!69  &  DENSTY(47,I)*2.686754E14,                                       &
!69  &  DENSTY(44,I)*2.686754E14,                                       &
!69  &  DENSTY(46,I)*2.686754E14,                                       &
!69  &  DENSTY(50,I)*2.686754E14,                                       &
!69  &  DENSTY(54,I)*2.686754E14,                                       &
!69  &  DENSTY(56,I)*2.686754E14,                                       &
!69  &  DENSTY(55,I)*2.686754E14,                                       &
!69  &  DENSTY(52,I)*2.686754E14,                                       &
!69  &  DENSTY(11,I)*2.686754E14,I=1,ML)

!     MOD3D PROFILES:
!     WRITE(69,'(2A)')'LAY Z[KM] P[ATM] T[K] H2O[GM/M3] CO2[GM/M3]',
!    1  ' O3[GM/M3] AER[KM-1 @ 550nm] CLOUD[GM/M3] RAIN [(MM/HR)**.63]'
!     WRITE(69,'((A,0P,F5.1,1P,E10.3,0P,F6.1,1P,4E9.2,0P,2F4.1))')
!    1  ('!',ZM(I),PM(I)/1013.25,TM(I),DENSTY(17,I)*18.015/2241.383,
!    2  DENSTY(36,I)*44.010/2241.383,DENSTY(31,I)*47.998/2241.383,
!    3  DENSTY(7,I)+DENSTY(12,I)+DENSTY(13,I)+DENSTY(14,I),
!    4  DENSTY(66,I),DENSTY(3,I),I=1,ML)
      WRITE(IPR,'(//A,//A,13(1X,A8:),/(14X,13(1X,A8:)))')               &
     &  '1 ATMOSPHERIC PROFILES','   I      Z   ',                      &
     &  (CNAMEX(K),K=1,NMOLX)
      WRITE(IPR,'(9X,A,50X,A,50X,A)')'(KM)   (','ATM CM/KM',')'
      DO 30 I=1,ML
          WRITE(IPR,'(I4,0P,F9.4,1X,1P,13E9.2:,/(14X,13E9.2:))')        &
     &      I,ZM(I),(DNSTYX(K,I),K=1,NMOLX)
   30 CONTINUE
      RETURN
      END
