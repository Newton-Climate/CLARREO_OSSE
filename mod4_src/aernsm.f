      SUBROUTINE AERNSM(JPRT,GNDALT,MARIC1,MARK,ICH,LMODEL,FULL_ATM,    &
     &                  A_FLAG,M_AEROSOL,O_AEROSOL,AC_WAVELEN,          &
     &                  NUM_LEVS,NUM_AC_WVL)

!     ROUTINE AERNSM DEFINES ALTITUDE, PRESSURE, TEMPERATURE,
!     MOLECULAR, AEROSOL, CLOUD, AND RAIN PROFILES.

!     ARGUMENTS:
      INTEGER JPRT,MARIC1,MARK,ICH(4)
      LOGICAL LMODEL
      REAL GNDALT
      INTEGER NUM_LEVS,NUM_AC_WVL
      REAL FULL_ATM(NUM_LEVS,28)      !DRF
      LOGICAL A_FLAG
      REAL M_AEROSOL(NUM_LEVS,4),O_AEROSOL(NUM_AC_WVL,12)
      REAL AC_WAVELEN(NUM_AC_WVL)

!     PARAMETERS:
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:
      INCLUDE 'BASE.h'
      REAL P,T,WH,WCO2,WO,WN2O,WCO,WCH4,WO2
      COMMON/MDATA/P(LAYDIM),T(LAYDIM),WH(LAYDIM),WCO2(LAYDIM),         &
     &  WO(LAYDIM),WN2O(LAYDIM),WCO(LAYDIM),WCH4(LAYDIM),WO2(LAYDIM)
      REAL WMOLXT
      COMMON/MDATAX/WMOLXT(NMOLX,LAYDIM)
      REAL WNO,WSO2,WNO2,WNH3,WHNO3
      COMMON/MDATA1/WNO(LAYDIM),WSO2(LAYDIM),WNO2(LAYDIM),              &
     &  WNH3(LAYDIM),WHNO3(LAYDIM)
      INCLUDE 'IFIL.h'

!     /CARD1/
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,MODTRN
      INTEGER M4,M5,M6,MDEF,IRD1,IRD2
      COMMON/CARD1A/M4,M5,M6,MDEF,IRD1,IRD2

!     /CARD1B/
      INTEGER JUNITP,JUNITT,JUNIT,JLOW
      REAL WMOL
      COMMON/CARD1B/JUNITP,JUNITT,JUNIT(13),WMOL(12),JLOW

!     /JM1B/
      CHARACTER*16 JCHAR
      COMMON/JM1B/JCHAR

      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA
      REAL VIS,WSS,WHH,RAINRT
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 &
     &  VIS,WSS,WHH,RAINRT
      INTEGER NCRALT,NCRSPC
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      &
     &  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      INTEGER IREG,IREGC
      REAL ALTB
      COMMON/CARD2D/IREG(4),ALTB(4),IREGC(4)

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
      REAL ZVSA,RHVSA,AHVSA
      INTEGER IHVSA
      COMMON/ZVSALY/ZVSA(10),RHVSA(10),AHVSA(10),IHVSA(10)
      REAL HMDLZ
      COMMON/MDLZ/HMDLZ(8)
      CHARACTER*20 HHAZE(16),HSEASN(2),HVULCN(8),HMET(2),HMODEL(8)
      CHARACTER*26 HTRRAD(4)
      COMMON/TITL/HHAZE,HSEASN,HVULCN,HMET,HMODEL,HTRRAD

!     /CRD1BX/
      INTEGER JUNITX
      REAL WMOLX
      COMMON/CRD1BX/JUNITX,WMOLX(NMOLX)

!     /JM2C3/
      REAL AHAZE(4),EQLWCZ,RRATZ
      INTEGER IHA1,ICLD1,IVUL1,ISEA1,ICHR
      COMMON/JM2C3/AHAZE,EQLWCZ,RRATZ,IHA1,ICLD1,                       &
     &   IVUL1,ISEA1,ICHR

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /TITL/
      EXTERNAL DEVCBD,TITLE

!     FUNCTIONS:
!       JOU      TRANSLATE UNIT SPECIFIER CHARACTER INTO INTEGER LABEL.
      INTEGER JOU
      REAL CLDPRF,AERPRF

!     LOCAL VARIABLES:
!       JDEF     VALUE OF JUNIT FOR MDEF SPECIES.
!       LCIRZ    FLAG, TRUE IF CIRRUS ALTITUDES NEED TO BE DEFINED.
      CHARACTER*20 HHOL,AHOL1,AHOL2
      INTEGER JDEF,ICLDS,ITYAER,IC1,NUMAER,IRD0,K,I
      LOGICAL LCIRZ,LDESRT
      REAL WH100,RH,VIS1,CLDD,CLD0,CLD1,CLD2,CLD3,FAC,HAZ1,HAZ2
      LOGICAL DRF_FLAG

!     LOCAL ARRAYS:
      INTEGER ITY1(LAYDIM),IH1(LAYDIM),IS1(LAYDIM),IVL1(LAYDIM)
      REAL CLDTOP(10),AHAST(LAYDIM),ZGN(LAYDIM)

!     DATA:
      CHARACTER*20 AHAHOL(13)
      DATA AHAHOL/             'CUMULUS             ',                  &
     &  'ALTOSTRATUS         ','STRATUS             ',                  &
     &  'STRATUS STRATO CUM  ','NIMBOSTRATUS        ',                  &
     &  'DRIZZLE 2.0 MM/HR   ','LT RAIN 5.0 MM/HR   ',                  &
     &  'MOD RAIN 12.5 MM/HR ','HEAVY RAIN 25 MM/HR ',                  &
     &  'EXTREME RAIN 75MM/HR','USER ATMOSPHERE     ',                  &
     &  'USER RAIN NO CLOUD  ','CIRRUS CLOUD        '/
      DATA CLDTOP/3.,3.,1.,2.,.66,1.,.66,.66,3.,3./

!     INITIALIZATIONS:

      DRF_FLAG = .TRUE.
      IC1=1
      NUMAER=7
      IRD0=1
      IF(LMODEL)THEN
          IRD0=0
      ELSEIF(IM.NE.1)THEN
          RETURN
      ELSEIF(IVSA.EQ.1)THEN

!         VERTICAL STRUCTURE ALGORITHM:
          IF(MODEL.EQ.0)THEN
              WRITE(IPR,'(/2A)')' ERROR: ',                             &
     &          ' MODEL equals 0 and army (VSA) model cannot mix.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP ' MODEL equals 0 and army (VSA) model cannot mix.'
          ENDIF
          IRD0=0
          IRD1=0
          IRD2=0
          ML=ML+10-JLOW
          IF(ML.GT.LAYDIM)THEN
              WRITE(IPR,'(/2A)')' WARNING:  ML exceeds parameter',      &
     &          ' LAYDIM and army (VSA) model top layer was truncated.'
              ML=LAYDIM
          ENDIF
          ZVSA(10)=ZVSA(9)+.01
          RHVSA(10)=0.
          AHVSA(10)=0.
          IHVSA(10)=0
      ENDIF

!     INITIALIZE COMMON /CARD2D/:
      IREGC(1)=0
      IREGC(2)=0
      IREGC(3)=0
      IREGC(4)=0
      !IREGC(5)=0
      !IREGC(6)=0
      !IREGC(7)=0
      !IREGC(8)=0
      !IREGC(9)=0
      !IREGC(10)=0
      !IREGC(11)=0
      !IREGC(12)=0
      !IREGC(13)=0
      !IREGC(14)=0
      ALTB(1)=0.
      ALTB(2)=0.
      ALTB(3)=0.
      ALTB(4)=0.
      !ALTB(5)=0.
      !ALTB(6)=0.
      !ALTB(7)=0.
      !ALTB(8)=0.
      !ALTB(9)=0.
      !ALTB(10)=0.
      !ALTB(11)=0.
      !ALTB(12)=0.
      !ALTB(13)=0.
      !ALTB(14)=0.
      LDESRT=.TRUE.
      IF(ICLD.EQ.18 .OR. ICLD.EQ.19)CALL CIRR18(ICLD,LCIRZ)

!     DEFAULTS:
      IF(IVULCN.LE.0)IVULCN=1
      IF(ISEASN.LE.0)ISEASN=1
      IF(.NOT.LJMASS)                                                   &
     &    WRITE(IPR,'(/A,I3)')' MODEL ATMOSPHERE NO.',MODEL
      IF(LMODEL)THEN
          CALL FLAYZ(ML,ICLD,GNDALT,IVSA)
      ELSE
          IF(MDEF.EQ.0)THEN
              JDEF=10
          ELSE
              JDEF=6
          ENDIF
          IF(.NOT.LJMASS)WRITE(IPR,'(/10X,A)')                          &
     &      ' MODEL 0, 7, or 8 USER INPUT DATA:'
      ENDIF

!     LOOP OVER LAYERS:
      IF(ML.GT.LAYDIM)THEN
          WRITE(IPR,'(/2A,I5,A,/19X,A,I5,A)')' Error in AERNSM: ',      &
     &      ' Number of atmospheric levels (ML =',ML,') exceeds',       &
     &      'parameter LAYDIM (=',LAYDIM,').  LAYDIM must be increased.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP ' Error:  Parameter LAYDIM must be increased.'
      ENDIF
      DO 50 K=1,ML
          LRHSET(K)=.FALSE.
          RH=0.
          WH(K)=0.
          WO(K)=0.
          IHA1=0
          ICLD1=0
          ISEA1=0
          IVUL1=0
          VIS1=0.
          AHAZE(1)=0.
          AHAZE(2)=0.
          AHAZE(3)=0.
          AHAZE(4)=0.
          !AHAZE(5)=0.
          !AHAZE(6)=0.
          !AHAZE(7)=0.
          !AHAZE(8)=0.
          !AHAZE(9)=0.
          !AHAZE(10)=0.
          !AHAZE(11)=0.
          !AHAZE(12)=0.
          !AHAZE(13)=0.
          !AHAZE(14)=0.
          EQLWCZ=0.
          RRATZ=0.
          ICHR=0
          DENSTY(16,K)=0.
          WCO2(K)=0.
          WCO(K)=0.
          WCH4(K)=0.
          WN2O(K)=0.
          WO2(K) =0.
          WNH3(K)=0.
          WNO (K)=0.
          WNO2(K)=0.
          WSO2(K)=0.
          WHNO3(K)= 0.
          WMOL(1)=0.
          WMOL(2)=0.
          WMOL(3)=0.
          WMOL(4)=0.
          WMOL(5)=0.
          WMOL(6)=0.
          WMOL(7)=0.
          WMOL(8)=0.
          WMOL(9)=0.
          WMOL(10)=0.
          WMOL(11)=0.
          WMOL(12)=0.
          JCHAR='                '
          IF(IRD0.EQ.1)THEN

!             READ IN USER-SPECIFIED ATMOSPHERE.
!             FOR MOLECULAR SPECIES, JCHAR IS DEFINED AS FOLLOWS:
!               JCHAR   JUNIT
!               -----   -----
!               " ",A     10    VOLUME MIXING RATIO (PPMV)
!                 B       11    NUMBER DENSITY (CM-3)
!                 C       12    MASS MIXING RATIO (GM(K)/KG(AIR))
!                 D       13    MASS DENSITY (GM M-3)
!                 E       14    PARTIAL PRESSURE (MB)
!                 F       15    DEW POINT TEMPERATURE (K) - H2O ONLY
!                 G       16    DEW POINT TEMPERATURE (C) - H2O ONLY
!                 H       17    RELATIVE HUMIDITY (%) - H2O ONLY
!                1-6     1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE

!             OTHER 'JCHAR' SPECIFICATIONS -
!               JCHAR   JUNIT
!               -----   -----
!               " ",A     10    PRESSURE IN (MB)
!                 B       11    PRESSURE IN (ATM)
!                 C       12    PRESSURE IN (TORR)
!                1-6     1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE

!               " ",A     10    AMBIENT TEMPERATURE (K)
!                 B       11    AMBIENT TEMPERATURE (C)
!                 C       12    AMBIENT TEMPERATURE (F)
!                1-6     1-6    DEFAULT TO SPECIFIED MODEL ATMOSPHERE
              IF(LJMASS)THEN
                  CALL INITCARD( 'CARD2C1' )
              ELSE
                 IF(DRF_FLAG) THEN
                    ZM(K) = FULL_ATM(K,1)
                    P(K) =  FULL_ATM(K,2)
                    T(K) =  FULL_ATM(K,3)
                    WMOL(1) = FULL_ATM(K,4)
                    WMOL(2) = FULL_ATM(K,5)
                    WMOL(3) = FULL_ATM(K,6)
                    JCHAR = 'AAAAAAAAAAAAAA A'
		    !WRITE(IPR,*) 'K = ',K
		    !WRITE(IPR,*) 'ML = ',ML
		    !WRITE(IPR,*) 'ZM(K) = ',ZM(K)
                    !WRITE(IPR,*) 'P(K) = ',P(K)
                    !WRITE(IPR,*) 'T(K) = ',T(K)
                 ELSE
                  READ(IRD,'(F10.0,5E10.3,A16)')                        &
     &              ZM(K),P(K),T(K),WMOL(1),WMOL(2),WMOL(3),JCHAR
                  ENDIF
                  WRITE(IPR,'(F10.5,1P,5E10.3,10X,A16)')                &
     &              ZM(K),P(K),T(K),WMOL(1),WMOL(2),WMOL(3),JCHAR
              ENDIF
              IF(IRD1.EQ.1)THEN
                  IF(LJMASS)THEN
                      CALL INITCARD( 'CARD2C2' )
                  ELSE
                     IF(DRF_FLAG) THEN
                        WMOL(4) = FULL_ATM(K,7)
                        WMOL(5) = FULL_ATM(K,8)
                        WMOL(6) = FULL_ATM(K,9)
                        WMOL(7) = FULL_ATM(K,10)
                        WMOL(8) = FULL_ATM(K,11)
                        WMOL(9) = FULL_ATM(K,12)
                        WMOL(10) = FULL_ATM(K,13)
                        WMOL(11) = FULL_ATM(K,14)
                        WMOL(12) = FULL_ATM(K,15)
                     ELSE
                      READ(IRD,'((8F10.0))')(WMOL(I),I=4,12)
                     ENDIF
                      WRITE(IPR,'((1P,8E10.3))')(WMOL(I),I=4,12)
                  ENDIF
                  IF(MDEF.EQ.2)THEN

!                     THE EXTRA SPECIES (I.E. THE WMOLX SPECIES) ARE
!                     READ IF MDEF=2 AND IRD1=1, 8 SPECIES PER LINE:
                      IF(LJMASS)THEN
                          CALL INITCARD( 'CARD2C2X' )
                      ELSE
                         IF(DRF_FLAG) THEN
                            WMOLX(1) = FULL_ATM(K,16)
                            WMOLX(2) = FULL_ATM(K,17)
                         ELSE
                          READ(IRD,'((8F10.0))')(WMOLX(I),I=1,NMOLX)
                         ENDIF 
                          WRITE(IPR,'((1P,8E10.3))')(WMOLX(I),I=1,NMOLX)
                      ENDIF
                  ENDIF
              ENDIF
              IF(IRD2.EQ.1)THEN

!                 CARD2C3:  ONLY ONE OF IHA1, ICLD1  OR IVUL1 IS
!                           ALLOWED. IF IHA1>0, OTHERS IGNORED; IF
!                           IHA1=0 AND ICLD1=0, USE ICLD1.  IF AHAZE
!                           AND EQLWCZ ARE BOTH ZERO, DEFAULT PROFILES
!                           ARE LOADED FROM IHA1, ICLD1, OR IVUL1.

!                   AHAZE    AEROSOL EXTINCTION AT 550 NM [KM-1].
!                   EQLWCZ   LIQUID WATER CONTENT AT ALTITUDE Z FOR
!                            AEROSOL, CLOUD OR FOG MODELS [PPMV].
!                   RRATZ    RAIN RATE AT ALTITUDE Z [MM/HR].
!                   IHA1     BOUNDARY LAYER AEROSOL MODEL USED FOR
!                            SPECTRAL EXTINCTION.
!                   IVUL1    STRATOSPHERIC AEROSOL MODEL USED FOR
!                            SPECTRAL EXTINCTION.
!                   ICLD1    CLOUD MODEL USED FOR SPECTRAL EXTINCTION.
!                   ISEA1    AEROSOL SEASON CONTROL FOR ALTITUDE Z.
!                   ICHR     AEROSOL PROFILE REGION SWITCH FOR IHA=7.
                  IF(LJMASS)THEN
                      CALL INITCARD( 'CARD2C3' )
                  ELSE
                     IF(.NOT.DRF_FLAG) THEN
                      READ(IRD,'(10X,3F10.0,5I5)')AHAZE(1),             &
     &                  EQLWCZ,RRATZ,IHA1,ICLD1,IVUL1,ISEA1,ICHR
                      ENDIF
                      WRITE(IPR,'(10X,3F10.3,5I5)')AHAZE(1),            &
     &                  EQLWCZ,RRATZ,IHA1,ICLD1,IVUL1,ISEA1,ICHR
                  ENDIF
                  AHAZE(2)=0.
                  AHAZE(3)=0.
                  AHAZE(4)=0.
                  !AHAZE(5)=0.
                  !AHAZE(6)=0.
                  !AHAZE(7)=0.
                  !AHAZE(8)=0.
                  !AHAZE(9)=0.
                  !AHAZE(10)=0.
                  !AHAZE(11)=0.
                  !AHAZE(12)=0.
                  !AHAZE(13)=0.
                  !AHAZE(14)=0.
              ELSEIF(IRD2.EQ.2)THEN
!                 READ 4 AEROSOL PROFILES:
                  IF(LJMASS)THEN
                      CALL INITCARD( 'CARD2C3' )
                  ELSE
                     IF(DRF_FLAG)THEN
                      IF (A_FLAG) THEN
                         AHAZE(1) = M_AEROSOL(K,1)
                         AHAZE(2) = M_AEROSOL(K,2)
                         AHAZE(3) = M_AEROSOL(K,3)
                         AHAZE(4) = M_AEROSOL(K,4)
                         !AHAZE(5) = M_AEROSOL(K,5)
                         !AHAZE(6) = M_AEROSOL(K,6)
                         !AHAZE(7) = M_AEROSOL(K,7)
                         !AHAZE(8) = M_AEROSOL(K,8)
                         !AHAZE(9) = M_AEROSOL(K,9)
                         !AHAZE(10) = M_AEROSOL(K,10)
                         !AHAZE(11) = M_AEROSOL(K,11)
                         !AHAZE(12) = M_AEROSOL(K,12)
                         !AHAZE(13) = M_AEROSOL(K,13)
                         !AHAZE(14) = M_AEROSOL(K,14)
                      ELSE
                         AHAZE(1) = 0
                         AHAZE(2) = 0
                         AHAZE(3) = 0
                         AHAZE(4) = 0
                         !AHAZE(5) = 0
                         !AHAZE(6) = 0
                         !AHAZE(7) = 0
                         !AHAZE(8) = 0
                         !AHAZE(9) = 0
                         !AHAZE(10) = 0
                         !AHAZE(11) = 0
                         !AHAZE(12) = 0
                         !AHAZE(13) = 0
                         !AHAZE(14) = 0
                      ENDIF
                      RRATZ = 0
                     ELSE
                      READ(IRD,'(10X,F10.3,10X,4F10.3)')                &
     &                  AHAZE(1),RRATZ,AHAZE(2),AHAZE(3),AHAZE(4)
                     ENDIF
                      WRITE(IPR,'(10X,F10.3,10X,4F10.3)')               &
     &                  AHAZE(1),RRATZ,AHAZE(2),AHAZE(3),AHAZE(4)
                  ENDIF
              ENDIF
          ENDIF
          IF((IHA1.EQ.0 .AND. ICLD1.NE.11) .OR. IHA1.NE.7)ICHR=0
          IF(MODEL.EQ.0)THEN

!             HORIZONTAL (CONSTANT PRESSURE) PATH.
              HMDLZ(1)=ZM(K)
              HMDLZ(2)=P(K)
              HMDLZ(3)=T(K)
              HMDLZ(4)=WMOL(1)
              HMDLZ(5)=WMOL(2)
              HMDLZ(6)=WMOL(3)
              HMDLZ(7)=AHAZE(1)
          ENDIF

!         JUNITP, JUNITT, JUNIT() AND JUNITX:
          IF(IRD0.EQ.0)THEN
              JUNITP=M1
              JUNITT=M1
              JUNIT(1)=M2
              JUNIT(2)=6
              JUNIT(3)=M3
              JUNIT(4)=M5
              JUNIT(5)=M6
              JUNIT(6)=M4
              JUNIT(7)=6
              JUNIT(8)=6
              JUNIT(9)=6
              JUNIT(10)=6
              JUNIT(11)=6
              JUNIT(12)=6
              JUNIT(13)=6
              JUNITX=6
          ELSE
              JUNITP=JOU(JCHAR(1:1),M1)
              JUNITT=JOU(JCHAR(2:2),M1)
              JUNIT(1)=JOU(JCHAR(3:3),M2)
              JUNIT(2)=JOU(JCHAR(4:4),JDEF)
              JUNIT(3)=JOU(JCHAR(5:5),M3)
              JUNIT(4)=JOU(JCHAR(6:6),M5)
              JUNIT(5)=JOU(JCHAR(7:7),M6)
              JUNIT(6)=JOU(JCHAR(8:8),M4)
              JUNIT(7)=JOU(JCHAR(9:9),JDEF)
              JUNIT(8)=JOU(JCHAR(10:10),JDEF)
              JUNIT(9)=JOU(JCHAR(11:11),JDEF)
              JUNIT(10)=JOU(JCHAR(12:12),JDEF)
              JUNIT(11)=JOU(JCHAR(13:13),JDEF)
              JUNIT(12)=JOU(JCHAR(14:14),JDEF)
              JUNIT(13)=JOU(JCHAR(15:15),JDEF)
              JUNITX=JOU(JCHAR(16:16),JDEF)
          ENDIF
          IF(IVSA.EQ.1 .AND. .NOT.LMODEL)THEN
              IF(MODEL.EQ.8)THEN
                  WRITE(IPR,'(/2A,/(18X,A))')' Error in AERNSM:  The',  &
     &              ' pressure-dependent radiosonde data (MODEL=8)',    &
     &              ' option cannot be used with the Army Vertical',    &
     &              ' Structure Algorithm (IVSA=1).'
                  IF(LJMASS)CALL WRTBUF(FATAL)
                  STOP 'Error: MODEL=8 cannot be coupled with IVSA=1.'
              ENDIF
              CALL VSANSM(K,AHAZE(1),IHA1,ZM(K))
          ELSE
              CALL CHECKP(P(K),JUNITP)
              CALL CHECKT(T(K),JUNITT)
              IF(MODEL.EQ.8)THEN
                  CALL DFLTP(K,GNDALT)
              ELSE
                  CALL DEFALT(ZM(K),P(K),T(K))
              ENDIF
              CALL CONVRT(P(K),T(K))
              WH(K)=WMOL(1)
              WCO2(K)=WMOL(2)
              WO(K)=WMOL(3)
              WN2O(K)=WMOL(4)
              WCO(K)=WMOL(5)
              WCH4(K)=WMOL(6)
              WO2(K)=WMOL(7)
              WNO(K)=WMOL(8)
              WSO2(K)=WMOL(9)
              WNO2(K)=WMOL(10)
              WNH3(K)=WMOL(11)
              WHNO3(K)=WMOL(12)
              DO 10 I=1,NMOLX
                  WMOLXT(I,K)=WMOLX(I)
   10         CONTINUE
          ENDIF

!         WITH ALTITUDE SET, RAIN RATE CAN BE DEFINED:
          IF(IRD2.EQ.0 .AND. ZM(K).LE.6.)RRATZ=RAINRT

!         SANITY CHECKS ON PRESSURE AND TEMPERATURE:
          IF(P(K).LE.0.)THEN
              WRITE(IPR,'(/2A,I3,A,1P,E12.4,A)')' Error in routine',    &
     &          ' AERNSM:  The pressure at level',K,' is',P(K),' mbars.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP ' Non-positive pressure encountered.'
          ELSEIF(K.GT.1)THEN
              IF(P(K).GE.P(K-1))                                        &
     &          WRITE(IPR,'(/2A,I3,/30X,2(A,1P,E12.4,A,I3))')           &
     &          ' Warning from routine AERNSM:  A pressure',            &
     &          ' inversion occurs between level',K-1,' (P=',           &
     &          P(K-1),' mbar) and level',K,' (P=',P(K),' mbar).'
          ENDIF
          IF(T(K).LE.0.)THEN
              WRITE(IPR,'(/2A,I3,A,F9.3,A)')' Error in routine',        &
     &          ' AERNSM:  The temperature at level',K,' is',T(K),' K.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP ' Non-positive temperature encountered.'
          ELSEIF(T(K).LT.100. .OR. T(K).GT.400.)THEN
              WRITE(IPR,'(/2A,I3,A,F9.3,A)')' Warning from routine',    &
     &          ' AERNSM:  The temperature at level',K,' is',T(K),' K.'
          ENDIF

!         CHECK GROUND ALTITUDE:
          IF(GNDALT.NE.ZM(1))THEN

!             FIX GROUND ALTITUDE TO THE BOTTOM OF ATMOSPHERIC PROFILE.
              IF(ABS(GNDALT-ZM(1)).GT..0001 .AND. .NOT.LJMASS)          &
     &          WRITE(IPR,'(2A,F8.4,A,/10X,A,F8.4,A)')' WARNING:  THE', &
     &          ' INPUT GROUND ALTITUDE (',GNDALT,'KM) IS BEING SET',   &
     &          ' TO THE BOTTOM OF THE ATMOSPHERIC PROFILE,',ZM(1),'KM.'
              GNDALT=ZM(1)
          ENDIF

!         PLACE ORIGINAL ALTITUDES IN ZGN
          IF(ZM(K).LT.6.)THEN
              ZGN(K)=6.*(ZM(K)-GNDALT)/(6.-GNDALT)
          ELSE
              ZGN(K)=ZM(K)
          ENDIF
          ICLDS=ICLD1
          IF(ICLD1.EQ.0)ICLD1=ICLD
          IF(ICLD1.GT.11)ICLD1=0
          IF(IHA1.NE.0)THEN
              IVUL1=0
              ICLD1=0
          ELSEIF(ICLD1.NE.0)THEN
              IVUL1=0
          ENDIF
          IF(AHAZE(1).EQ.0. .AND. EQLWCZ.EQ.0.)THEN
              IF(IVSA.EQ.1 .AND. ICLD1.EQ.0)THEN
                  IF(MODEL.LT.7)CALL LAYVSA(K,RH,AHAZE(1),IHA1)
              ELSE
                  CALL LAYCLD(K,EQLWCZ,RRATZ,ICLD1,GNDALT)
                  IF(RAINRT.GT.0. AND. ZM(K).LT.6.)RRATZ=RAINRT
                  IF(ICLD1.GE.1 .AND. ICLD1.LE.10)THEN
                      IF(ZM(K).GT.CLDTOP(ICLD1)+GNDALT)RRATZ=0.
                  ENDIF
              ENDIF
          ENDIF

!         DENSTY(16,K):
          IF(ICLDS.EQ.18 .OR. ICLDS.EQ.19)THEN
              IF(AHAZE(1).GT.0.)THEN
                  DENSTY(16,K)=AHAZE(1)
                  AHAZE(1)=0.
              ELSEIF(EQLWCZ.GT.0)THEN
                  IF(ICLDS.EQ.18)THEN
                      DENSTY(16,K)=EQLWCZ/5.811E-2
                  ELSE
                      DENSTY(16,K)=EQLWCZ/3.446E-3
                  ENDIF
                  EQLWCZ=0.
              ENDIF
          ELSEIF(ICLDS.EQ.0 .AND. (ICLD.EQ.18 .OR. ICLD.EQ.19))THEN
              IF(LCIRZ)THEN
                  CLDD=.1*CTHIK
                  CLD0=CALT-.5*CLDD
                  IF(CLD0.LE.GNDALT)CLD0=GNDALT
                  CLD1=CLD0+CLDD
                  CLD2=CLD0+CTHIK
                  CLD3=CLD1+CTHIK
                  LCIRZ=.FALSE.
              ENDIF
              IF(ZM(K).LE.CLD0 .OR. ZM(K).GE.CLD3)THEN
                  DENSTY(16,K)=0.
              ELSEIF(ZM(K).LT.CLD1)THEN
                  DENSTY(16,K)=CEXT*(ZM(K)-CLD0)/CLDD
              ELSEIF(ZM(K).LE.CLD2)THEN
                  DENSTY(16,K)=CEXT
              ELSE
                  DENSTY(16,K)=CEXT*(CLD3-ZM(K))/CLDD
              ENDIF
          ENDIF
          DENSTY(66,K)=EQLWCZ
          IF(ICLDS.EQ.0 .AND. EQLWCZ.EQ.0.)ICLD1=0
          DENSTY(3,K)=RRATZ
          IF(LMODEL .AND. (EQLWCZ.GT.0. .OR. RRATZ.GT.0.))RH=100.
          AHAST(K)=AHAZE(1)

!           IHA1    IHAZE FOR THIS LAYER
!           ISEA1   ISEASN FOR THIS LAYER
!           IVUL1   IVULCN FOR THE LAYER
          IF(ISEA1.EQ.0)ISEA1=ISEASN
          IF(IHA1.GT.0)THEN
              ITYAER=IHA1
          ELSE
              ITYAER=IHAZE
          ENDIF
          IF(IVUL1.GT.0)THEN
              IVULCN=IVUL1
          ELSE
              IVUL1=IVULCN
          ENDIF
          IF(K.GT.1)THEN
              IF(ICHR.NE.1)THEN
                  IF(ICLD1.EQ.IREGC(IC1))THEN
                      IF(IHA1.EQ.0 .AND. ICLD1.EQ.0)THEN
                          IF(ZGN(K).GT.30.00001)THEN
                              ITYAER=19
                          ELSEIF(ZGN(K).GT.10.00001)THEN
                              ITYAER=IVULCN+10
                          ELSEIF(ZGN(K).GT.2.00001)THEN
                              ITYAER=6
                          ENDIF
                          IF(ITYAER.EQ.ICH(IC1))GOTO20
                      ELSE
                          NUMAER=7
                          IF(IC1.GT.1)NUMAER=IC1+10
                          IF(IHA1.EQ.0 .OR. IHA1.EQ.ICH(IC1))GOTO20
                      ENDIF
                  ELSEIF(ICLD1.NE.0)THEN
                      IF(ICLD1.EQ.IREGC(1))THEN
                          NUMAER=7
                          ALTB(1)=ZM(K)
                          GOTO30
                      ELSEIF(IC1.GT.1 .AND. ICLD1.EQ.IREGC(2))THEN
                          NUMAER=12
                          ALTB(2)=ZM(K)
                          GOTO30
                      ELSEIF(IC1.GT.2 .AND. ICLD1.EQ.IREGC(3))THEN
                          NUMAER=13
                          ALTB(3)=ZM(K)
                          GOTO30
                      ENDIF
                  ELSE
                      IF(IHA1.EQ.0 .AND. ICLD1.EQ.0)THEN
                          IF(ZGN(K).GT.30.00001)THEN
                              ITYAER=19
                          ELSEIF(ZGN(K).GT.10.00001)THEN
                              ITYAER=IVULCN+10
                          ELSEIF(ZGN(K).GT.2.00001)THEN
                              ITYAER=6
                          ENDIF
                      ENDIF
                      IF(ITYAER.EQ.ICH(1))THEN
                          NUMAER=7
                          ALTB(1)=ZM(K)
                          GOTO30
                      ELSEIF(IC1.GT.1 .AND. ITYAER.EQ.ICH(2))THEN
                          NUMAER=12
                          ALTB(2)=ZM(K)
                          GOTO30
                      ELSEIF(IC1.GT.2 .AND. ITYAER.EQ.ICH(3))THEN
                          NUMAER=13
                          ALTB(3)=ZM(K)
                          GOTO30
                      ENDIF
                  ENDIF
              ENDIF
              IF(IC1.LT.4)THEN
                  IC1=IC1+1
                  NUMAER=IC1+10
              ELSE
                  IC1=4
                  NUMAER=14
                  ITYAER=ICH(IC1)
              ENDIF
          ENDIF
   20     CONTINUE
          ICH(IC1)=ITYAER
          IREGC(IC1)=ICLD1
          ALTB(IC1)=ZM(K)
   30     CONTINUE

!         SATURATED WATER VAPOR DENSITY [GM / M3]
          WH100=273.15/T(K)
          WH100=WH100*EXP(18.9766-(14.9595+2.43882*WH100)*WH100)
          IF(RH.GT.0.)WH(K)=.01*RH*WH100
          DENSTY(7,K)=0.
          DENSTY(12,K)=0.
          DENSTY(13,K)=0.
          DENSTY(14,K)=0.
          DENSTY(15,K)=0.
          RELHUM(K)=0.
          IF(WH(K).GT.0.)THEN
              RELHUM(K)=100.*WH(K)/WH100
              IF(RELHUM(K).GT.100.)THEN
                  IF(RELHUM(K).GT.100.05)WRITE(IPR,'(/A,2(F12.5,A))')   &
     &              ' WARNING:  The relative humidity at',ZM(K),        &
     &              'km was reset from',RELHUM(K),'% to 100%'
                  RELHUM(K)=100.
                  WH(K)=WH100
              ELSEIF(RELHUM(K).LT.0.)THEN
                  WRITE(IPR,'(/A,2(F12.5,A))')                          &
     &              ' WARNING:  The relative humidity at',ZM(K),        &
     &              'km was reset from',RELHUM(K),'% to 0%'
                  RELHUM(K)=0.
                  WH(K)=0.
              ENDIF
          ENDIF
          IF(VIS1.LE.0.)VIS1=VIS
          IF(IRD2.EQ.2)THEN
              DENSTY(7,K)=AHAZE(1)
              DENSTY(12,K)=AHAZE(2)
              DENSTY(13,K)=AHAZE(3)
              DENSTY(14,K)=AHAZE(4)
	      !WRITE(*,*) 'in here' !YES IT GOES HERE
          ELSEIF(AHAZE(1).GT.0.)THEN
              DENSTY(NUMAER,K)=AHAZE(1)
              IF(ITYAER.NE.3 .AND. ITYAER.NE.10)GOTO40
          ENDIF
          IF(ITYAER.EQ.3 .AND. MARIC1.EQ.0)THEN
              CALL MARINE(VIS1,MODEL,RELHUM(K),                         &
     &          WSS,WHH,ICSTL,EXTC,ABSC,IC1)
              IREG(IC1)=1
              VIS=VIS1
              MARIC1=IC1
              MARK=K
          ELSEIF(ITYAER.EQ.10 .AND. LDESRT)THEN
              CALL DESATT(WSS,VIS1)
              IREG(IC1)=1
              VIS=VIS1
              LDESRT=.FALSE.
          ENDIF
          IF(IRD2.NE.2 .AND. AHAZE(1).LE.0.)THEN
              IF(IHA1.LE.0)IHA1=IHAZE
              IF(EQLWCZ.GT.0.)THEN
                  DENSTY(NUMAER,K)                                      &
     &              =CLDPRF(EQLWCZ,RELHUM(K),ICLD1,IHA1,IC1,ICH(IC1))
              ELSE
                  IF(ZM(K).LT.6.)THEN
                      I=IFIX(ZGN(K)+1.E-6)+1
                      FAC=ZGN(K)-FLOAT(I-1)
                  ELSEIF(ZM(K).GE.70.)THEN
                      I=32
                      FAC=(ZM(K)-70.)/30.
                  ELSEIF(ZM(K).GE.50.)THEN
                      I=31
                      FAC=(ZM(K)-50.)/20.
                  ELSEIF(ZM(K).GE.25.)THEN
                      FAC=(ZM(K)-25.)/5.
                      I=IFIX(FAC)
                      FAC=FAC-I
                      I=I+26
                  ELSE
                      I=IFIX(ZM(K))
                      FAC=ZM(K)-FLOAT(I)
                      I=I+1

                      IF (I.LE.0)THEN
!                        THIS FIX IS FOR THE ALTITUDES -1.0 KM.  IN THIS
!                        CASE, I .LE.0 AND FAC SHOULD BE ZM(K)-0.0.
!                        SET I T0 1 AND FAC TO THE APPROPRIATE VALUE.
                         FAC=ZM(K)-0.0
                         I=1
                      ENDIF

                  ENDIF
                  HAZ1=AERPRF(I,VIS1,IHA1,ISEA1,IVUL1)
                  DENSTY(NUMAER,K)=0.
                  IF(HAZ1.GT.0.)THEN
                      I=I+1
                      HAZ2=AERPRF(I,VIS1,IHA1,ISEA1,IVUL1)
                      IF(HAZ2.GT.0.)                                    &
     &                  DENSTY(NUMAER,K)=HAZ1*(HAZ2/HAZ1)**FAC
                  ENDIF
              ENDIF
          ENDIF
   40     CONTINUE
          ITY1(K)=ITYAER
          IF(AHAZE(1).EQ.0.)THEN
              IH1(K)=IHA1
          ELSE
              IH1(K)=-99
          ENDIF
          IS1(K)=ISEA1
          IVL1(K)=IVUL1
   50 CONTINUE
      JPRT=1
      IF(LJMASS .OR. (LMODEL .AND. IVSA.EQ.0 .AND. ICLD.EQ.0            &
     &  .AND. RAINRT.EQ.0. .AND. GNDALT.EQ.0.))RETURN
      JPRT=0
      IF(IVSA.EQ.1)THEN
          HHOL='VSA DEFINED         '
      ELSEIF(ICLD.GE.18)THEN
          HHOL=AHAHOL(13)
      ELSEIF(ICLD.LE.0 .OR. ICLD.GT.12)THEN
          HHOL=AHAHOL(12)
      ELSE
          HHOL=AHAHOL(ICLD)
      ENDIF
      IF(ICLD.NE.0)THEN
          WRITE(IPR,'(/2A)')' CLOUD AND/OR RAIN TYPE CHOSEN IS ',HHOL
          WRITE(IPR,'(//(2A))')                                         &
     &      '1     Z         P        T     REL H    H2O   ',           &
     &      '  CLD AMT   RAIN RATE                      AEROSOL',       &
     &      '     (KM)      (MB)     (K)     (%)  (GM / M3)',           &
     &      ' (GM / M3)  (MM / HR) TYPE                 PROFILE',       &
     &      '                              [Before scaling]'
      ELSE
          WRITE(IPR,'(//(2A))')                                         &
     &      '      Z         P        T     REL H    H2O   ',           &
     &      '                       AEROSOL',                           &
     &      '     (KM)      (MB)     (K)     (%)  (GM / M3)',           &
     &      '  TYPE                 PROFILE',                           &
     &      '                              [Before scaling]'
      ENDIF
      DO 60 K=1,ML
          IF(ITY1(K).LE.0)THEN
              ITYAER=1
          ELSEIF(ITY1(K).EQ.18)THEN
              ITYAER=13
          ELSEIF(ITY1(K).GE.16 .AND. ITY1(K).LE.19)THEN
              ITYAER=11
          ELSE
              ITYAER=ITY1(K)
          ENDIF
          IHA1=IH1(K)
          ISEA1=IS1(K)
          IVUL1=IVL1(K)
          IF((IVSA.EQ.1 .AND. K.LE.9) .OR. DENSTY(66,K).GT.0. .OR.      &
     &      DENSTY(3,K).GT.0. .OR. IHAZE.EQ.0)THEN
              AHOL1=HHOL
          ELSE
              AHOL1=HHAZE(ITYAER)
          ENDIF
          IF(AHAST(K).EQ.0.)THEN
              AHOL2=AHOL1
          ELSEIF(DENSTY(66,K).GT.0. .OR. DENSTY(3,K).GT.0.)THEN
              AHOL2=HHOL
          ELSE
              AHOL2='USER DEFINED        '
          ENDIF
          IF(ICLD.NE.0)THEN
              IF(ZGN(K).GT.2.00001)THEN
                  WRITE(IPR,'(2F10.3,2F8.2,1P,3E10.3,1X,3A)')           &
     &              ZM(K),P(K),T(K),RELHUM(K),WH(K),DENSTY(66,K),       &
     &              DENSTY(3,K),AHOL1,AHOL2,HSEASN(ISEA1)(1:13)
              ELSE
                  WRITE(IPR,'(2F10.3,2F8.2,1P,3E10.3,1X,2A)')           &
     &              ZM(K),P(K),T(K),RELHUM(K),WH(K),DENSTY(66,K),       &
     &              DENSTY(3,K),AHOL1,AHOL2
              ENDIF
          ELSE
              IF(ZGN(K).GT.2.00001)THEN
                  WRITE(IPR,'(2F10.3,2F8.2,1P,E10.3,1X,3A)')ZM(K),P(K), &
     &              T(K),RELHUM(K),WH(K),AHOL1,AHOL2,HSEASN(ISEA1)(1:13)
              ELSE
                  WRITE(IPR,'(2F10.3,2F8.2,1P,E10.3,1X,2A)')ZM(K),P(K), &
     &              T(K),RELHUM(K),WH(K),AHOL1,AHOL2
              ENDIF
          ENDIF
   60 CONTINUE
      RETURN
      END
