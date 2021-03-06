      SUBROUTINE MCMOL(VCEN,IK,IKMAX,NMWAVE,                            &
     &  T_H2O,T_UMX,T_O3,T_AUX,TL_H2O,TL_UMX,TL_O3,TL_AUX)

!     MCMOL WRITES MOLECULAR TRANSMITTANCE DATA FOR MONTE-CARLO
!     SIMULATIONS.  EITHER BAND MODEL RESOLUTION OR 1 NANOMETER SPECTRAL
!     BINS.  IN EITHER CASE, INTEGRATIONS ARE OVER SPECTRAL FREQUENCY.
      IMPLICIT NONE

!     PARAMETERS:
!       C_FH2O   CONVERSION FACTOR FOR H2O FOREIGN CONTINUUM.
!       C2       SECOND RADIATION CONSTANT [CM K].
      INCLUDE 'PARAMS.h'
      REAL C_FH2O,C2,C2_260,C2_296
      PARAMETER(C_FH2O=296*LOSCHM,                                      &
     &  C2=1.438786,C2_260=C2/260,C2_296=C2/296)

!     INPUT ARGUMENTS:
!       VCEN     SPECTRAL BAND CENTER FREQUENCY [CM-1].
!       IK       PATH SEGMENT INDEX.
!       IKMAX    NUMBER OF PATH SEGMENTS.
!       NMWAVE   WAVELENGTH BIN OF PREVIOUS DATA [NM].
!       T_H2O    DIRECT PATH MOLECULAR H2O TRANSMITTANCE.
!       T_UMX    DIRECT PATH UNIFORMLY MIXED GASES TRANSMITTANCE.
!       T_O3     DIRECT PATH MOLECULAR O3 TRANSMITTANCE.
!       T_AUX    DIRECT PATH AUXILIARY SPECIES TRANSMITTANCE.
!       TL_H2O   L-SHAPED PATH MOLECULAR H2O TRANSMITTANCE.
!       TL_UMX   L-SHAPED PATH UNIFORMLY MIXED GASES TRANSMITTANCE.
!       TL_O3    L-SHAPED PATH MOLECULAR O3 TRANSMITTANCE.
!       TL_AUX   L-SHAPED PATH AUXILIARY SPECIES TRANSMITTANCE.

!     OUTPUT ARGUMENTS:
!       NMWAVE   WAVELENGTH BIN OF CURRENT DATA [NM].
      INTEGER IK,IKMAX,NMWAVE
      REAL VCEN,T_H2O,T_UMX,T_O3,T_AUX,TL_H2O,TL_UMX,TL_O3,TL_AUX

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BMHEAD.h'
      INCLUDE 'BMDAT.h'

!     /JM5/
!       IRPT     REPEAT INPUT FLAG (0=NONE, 1=ALL, 3=GEOM, 4=SPEC).
!       IFAC     CURRENT COLUMN SCALING FACTOR INDEX.
!       NFACMN   NUMBER OF COLUMN SCALING FACTOR LESS THAN 1.
!       NFACMX   NUMBER OF COLUMN SCALING FACTOR GREATER THAN 1.
!       FACMC    CURRENT COLUMN SCALING FACTOR.
!       SCALMN   MINIMUM COLUMN SCALING FACTOR.
!       SCALMX   MAXIMUM COLUMN SCALING FACTOR.
!       LBMRES   LOGICAL FLAG, .TRUE. FOR BAND MODEL AND .FALSE.
!                FOR 1 NM SPECTRAL RESOLUTION OUTPUT.
      INTEGER IRPT,IFAC,NFACMN,NFACMX
      REAL FACMC
      DOUBLE PRECISION SCALMN,SCALMX
      LOGICAL LBMRES
      COMMON/JM5/SCALMN,SCALMX,IRPT,IFAC,NFACMN,NFACMX,FACMC,LBMRES

!     /ACTIVE/
!       NACT     NUMBER OF ACTIVE MOLECULES FOR CURRENT FREQUENCY BIN.
!       NACTBM   NUMBER OF ACTIVE BAND MODEL MOLECULES FOR CURRENT FREQ.
!       NACTX    NUMBER OF ACTIVE X & Y CROSS-SECTION MOLECULES.
!       MACTBM   LIST OF ACTIVE BAND MODEL MOLECULES FOR CURRENT FREQ.
!       MACTX    LIST OF ACTIVE X & Y CROSS-SECTION MOLECULES.
      INTEGER NACT,NACTBM,NACTX,MACTBM,MACTX
      COMMON/ACTIVE/NACT,NACTBM,NACTX,MACTBM(MMOLYT),MACTX(NMOLX)

!     /MC_MOL/
!       JT_MC    BAND MODEL TEMP INTERPOLATION INDEX AT LOS LEVELS.
!       FT_MC    BAND MODEL TEMP INTERPOLATION FRACTION AT LOS LEVELS.
!       KT_MC    CROSS-SECTION TEMP INTERPOLATION INDEX AT LOS LEVELS.
!       GT_MC    X-SECTION TEMP INTERPOLATION FRACTION AT LOS LEVELS.
!       P2_MC    PRESSURE^2 INTERPOLATION FRACTION AT PATH BOUNDARIES.
      INTEGER JT_MC,KT_MC
      REAL FT_MC,GT_MC,P2_MC
      COMMON/MC_MOL/JT_MC(0:LAYTWO),FT_MC(0:LAYTWO),                    &
     &              KT_MC(0:LAYTWO),GT_MC(0:LAYTWO),P2_MC(0:LAYTWO)

!     /MOL_F/
!       UMX_F    LOS LEVEL MOLE FRACTION RATIO, RELATIVE TO CO2.
!       AUX_F    LOS LEVEL AUX MOLECULAR GASES RELATIVE MOLE FRACTION.
      REAL UMX_F,AUX_F
      COMMON/MOL_F/UMX_F(2:NMOLXT,0:LAYTWO),AUX_F(1:MMOLY,0:LAYTWO)

!     /PATH/
!       PTHCOS   COSINE OF PATH ZENITH AT PATH BOUNDARIES.
!       PTHZEN   PATH ZENITH AT PATH BOUNDARIES [DEG].
!       PTHECA   SENSOR TO PATH EARTH CENTER ANGLE [DEG].
!       PTHALT   ALTITUDES AT PATH BOUNDARIES [KM].
!       PTH_MS   ALTITUDES AT PATH BOUNDARIES FOR THE MS PATH.
!       PTHSEG   PATH SEGMENT LENGTH [KM].
!       PTHRNG   SENSOR TO PATH BOUNDARY RANGE [KM].
!       JMAX     NUMBER OF DISTINCT LOS PATH SEGMENT ENDPOINT ALTITUDES.
!       IKHMIN   PATH BOUNDARY INDEX OF PATH MINIMUM ALTITUDE.
!       IKHMAX   PATH BOUNDARY INDEX OF PATH MAXIMUM ALTITUDE.
!       IKOUT    NUMBER OF PATH BOUNDARIES K DATA IS OUTPUT.
!       NTKDIS   RECORD NUMBER FOR K-DISTRIBUTION TRANSMITTANCE FILE.
!       NRKDIS   RECORD NUMBER FOR K-DISTRIBUTION RADIANCE FILE.
!       MAPPTH   MAPPING FROM PATH SEGMENT MIDPOINT TO VERTICAL LAYER.
!       IPTHHT   ALTITUDES (HEIGHTS) AT PATH BOUNDARIES [M].
!       LOWALT   VERTICAL LAYER BOUNDARY INDEX AT OR JUST BELOW PTHALT.
!       FACALT   ALTITUDE INTERPOLATION FRACTION FOR PTHALT
!       PATH_T   TEMPERATURE AT PATH BOUNDARIES [K].
!       PATH_P   PRESSURE AT PATH BOUNDARIES [ATM].
!       PTHRH    RELATIVE HUMIDITY AT PATH BOUNDARIES [K].
!       LSSGEO   LOGICAL FLAG, .TRUE. FOR SOLAR PATHS.
!       LTANMX   LOGICAL FLAG, .TRUE. IF PATH HAS A TANGENT MAXIMUM.
      DOUBLE PRECISION PTHCOS,PTHZEN,PTHECA,PTHALT,PTH_MS,PTHSEG,PTHRNG
      INTEGER JMAX,IKHMIN,IKHMAX,IKOUT,NTKDIS,NRKDIS,MAPPTH,            &
     &  IPTHHT,LOWALT
      REAL FACALT,PATH_T,PATH_P,PTHRH
      LOGICAL LSSGEO,LTANMX
      COMMON/PATH/PTHCOS(0:LAYTWO),PTHZEN(0:LAYTWO),PTHECA(0:LAYTWO),   &
     &  PTHALT(0:LAYTWO,1:MLOS),PTH_MS(0:LAYDIM),PTHSEG(LAYTWO),        &
     &  PTHRNG(0:LAYTWO,1:MLOS),JMAX,IKHMIN(MLOS),IKHMAX(MLOS),         &
     &  IKOUT(MLOS),MAPPTH(LAYTWO,1:MLOS),IPTHHT(0:LAYTWO),NTKDIS,      &
     &  NRKDIS,LOWALT(0:LAYTWO,1:MLOS),FACALT(0:LAYTWO,1:MLOS),         &
     &  PATH_T(0:LAYTWO,1:MLOS),PATH_P(0:LAYTWO,1:MLOS),                &
     &  PTHRH(0:LAYTWO,1:MLOS),LSSGEO,LTANMX

!     FUNCTIONS:
!       N2DATA   N2 CONTINUUM DATA [KM-1 / (ATM-AIR ATM-N2) AT 296K].
!       HERTDA   O2 HERZBERG BAND [CM-1/ATM AT 273.15].
!       O2MATE   O2 1.27 MICRON CONTINUUM ABSORPTION.
!       O3UV     O3 UV ABSORPTION COEFFICIENTS.
!       O3HHT0   O3 HARTLEY-HUGGINS DATA
!       O3HHT1   O3 HARTLEY-HUGGINS DATA
!       O3HHT2   O3 HARTLEY-HUGGINS DATA
!       NO2XS    NO2 SPECTRAL ABSORPTION CROSS-SECTION [CM-1/ATM].
!       SO2XS    SO2 SPECTRAL ABSORPTION CROSS-SECTION [CM-1/ATM].
!       CH4HI    0.4NM CH4 ABSORPTION FROM 520.1 TO 995.0 NM [CM-1/ATM].
!       CH4LO    1 NM CH4 ABSORPTION FROM 400.4 TO 1050.0 NM [CM-1/ATM].
      REAL N2DATA,HERTDA,O2MATE,O3UV,O3HHT0,O3HHT1,O3HHT2,NO2XS,SO2XS,  &
     &  CH4HI,CH4LO

!     LOCAL VARIABLES:
!       NM       CURRENT WAVELENGTH [NM].
!       I        PATH SEGMENT INDEX.
!       IACT     LOOP INDEX FOR ACTIVE MOLECULES.
!       IMOL     MOLECULAR INDEX.
!       JT_UP    TEMPERATURE INTERPOLATION UPPER INDEX.
!       JT_LO    TEMPERATURE INTERPOLATION LOWER INDEX.
!       FT       TEMPERATURE INTERPOLATION FRACTION.
!       P2       PRESSURE-SQUARED INTERPOLATION FRACTION.
!       ABSM     MOLECULAR ABSORPTION COEFFICIENT [CM-1/ATM].
!       NORM     NORMALIZATION FOR LINE TAIL ABSORPTION COEFS [ATM].
!       ITLSUB   LOOP INDEX FOR LINE TAIL ABSORPTION COEFFICIENTS.
!       TAILUP   LINE TAIL ABS COEF OVER PRES @ UPPER TEMP [CM-1/ATM^2].
!       TAILLO   LINE TAIL ABS COEF OVER PRES @ LOWER TEMP [CM-1/ATM^2].
!       TAIL     SUM OF LINE TAIL ABSORPTION COEFFICIENTS [CM-1/ATM].
!       AH2O     H2O WEAK-LINE ABSORPTION COEFFICIENT [CM-1/ATM].
!       AH2O_S   SELF-BROADENED H2O ABSORPTION COEFFICIENT [CM-2/ATM^2].
!       AUMX     UMIX GASES WEAK-LINE ABSORPTION COEFFICIENT [CM-1/ATM].
!       AO3      O3 WEAK-LINE ABSORPTION COEFFICIENT [CM-1/ATM].
!       AAUX     AUX GASES WEAK-LINE ABSORPTION COEFFICIENT [CM-1/ATM].
!       T_IK     LEVEL TEMPERATURE [K].
!       P_IK     LEVEL PRESSURE [ATM].
!       PoverT   P_IK OVER T_IK [ATM/K].
!       AIRATM   STANDARD AIR DENSITY [ATM].
!       EXPON    EXPONENTIAL IN RADIATION TERM.
!       DT220K   LEVEL TEMPERATURE MINUS 220K [K].
!       DT273K   LEVEL TEMPERATURE MINUS 273.15K [K].
      INTEGER NM,I,IACT,IMOL,JT_UP,JT_LO,ITLSUB
      REAL FT,P2,ABSM,NORM,TAILUP,TAILLO,TAIL,AH2O,AH2O_S,AUMX,AO3,AAUX,&
     &  SH2OT0,SH2OT,T_IK,P_IK,POVERT,AIRATM,EXPON,DT220K,DT273K

!     SAVED VARIABLES:
!       IREC     RECORD NUMBER.
!       N2CONT   N2 CONTINUUM ABS [KM-1 / (ATM-AIR ATM-N2) AT 296K].
!       FH2O     FOREIGN BROADENED H2O CONTINUUM ABS [K CM-1 / ATM].
!       O2HERZ   HERZBERG O2 ABSORPTION [CM-1/ATM AT 273.15].
!       O2_127   O2 1.27um CONT ABS [CM-1/(ATM-AIR ATM-O2) AT 273.15].
!       O2SIG    O2 6.4um CONT ABS [CM-1/(ATM-AIR ATM-O2) AT 273.15K].
!       O2SIGA   O2 6.4um CONT ABS [CM-1/(K ATM-AIR ATM-O2) AT 273.15K].
!       O2SIGB   O2 6.4um CONT ABS [CM-1/(K^2 ATMAIR ATMO2) AT 273.15K].
!       CRSNO2   NO2 SPECTRAL ABSORPTION CROSS-SECTION [CM-1/ATM].
!       CRSSO2   SO2 SPECTRAL ABSORPTION CROSS-SECTION [CM-1/ATM].
!       A_CH4    CH4 VIS/NIR ABSORPTION COEFFICIENT [CM-1/ATM].
!       WT       SPECTRAL BIN NORMALIZATION.
!       DWTBEG   BEGINNING (SHORTER WAVELENGTH) SPECTRAL BIN WEIGHT.
!       DWTEND   ENDING (LONGER WAVELENGTH) SPECTRAL BIN WEIGHT.
!       TLH2O    SENSOR-SUN BENT PATH MOLECULAR H2O TRANSMITTANCE.
!       TH2O     DIRECT PATH MOLECULAR H2O TRANSMITTANCE.
!       TLH2OS   SAVED SENSOR-SUN BENT PATH MOLECULAR H2O TRANSMITTANCE.
!       TH2OS    SAVED DIRECT PATH MOLECULAR H2O TRANSMITTANCE.
!       TLUMX    SENSOR-SUN BENT PATH UNIFORMLY MIXED GASES TRANS.
!       TUMX     DIRECT PATH UNIFORMLY MIXED GASES TRANSMITTANCE.
!       TLUMXS   SAVED SENSOR-SUN BENT PATH UNIFORMLY MIXED GASES TRANS.
!       TUMXS    SAVED DIRECT PATH UNIFORMLY MIXED GASES TRANSMITTANCE.
!       TLO3     SENSOR-SUN BENT PATH MOLECULAR O3 TRANSMITTANCE.
!       TO3      DIRECT PATH MOLECULAR O3 TRANSMITTANCE.
!       TLO3S    SAVED SENSOR-SUN BENT PATH MOLECULAR O3 TRANSMITTANCE.
!       TO3S     SAVED DIRECT PATH MOLECULAR O3 TRANSMITTANCE.
!       TLAUX    SENSOR-SUN BENT PATH AUXILIARY SPECIES TRANSMITTANCE.
!       TAUX     DIRECT PATH AUXILIARY SPECIES TRANSMITTANCE.
!       TLAUXS   SAVED SENSOR-SUN BENT PATH AUXILIARY SPECIES TRANS.
!       TAUXS    SAVED DIRECT PATH AUXILIARY SPECIES TRANSMITTANCE.
!       A_H2O    WEAK-LINE H2O ABSORPTION COEFFICIENT [CM-1/ATM].
!       A_H2OS   SAVED WEAK-LINE H2O ABSORPTION COEFFICIENT [CM-1/ATM].
!       A_UMX    WEAK-LINE UMX ABSORPTION COEFFICIENT [CM-1/ATM].
!       A_UMXS   SAVED WEAK-LINE UMX ABSORPTION COEFFICIENT [CM-1/ATM].
!       A_O3     WEAK-LINE O3  ABSORPTION COEFFICIENT [CM-1/ATM].
!       A_O3 S   SAVED WEAK-LINE O3  ABSORPTION COEFFICIENT [CM-1/ATM].
!       A_AUX    WEAK-LINE AUX ABSORPTION COEFFICIENT [CM-1/ATM].
!       A_AUXS   SAVED WEAK-LINE AUX ABSORPTION COEFFICIENT [CM-1/ATM].
      INTEGER IREC
      REAL N2CONT,FH2O,O2HERZ,O2_127,O2SIG,O2SIGA,O2SIGB,               &
     &  O3T0,O3T1,O3T2,CRSNO2,CRSSO2,A_CH4,WT,DWTBEG,DWTEND,            &
     &  TLH2O(0:LAYTWO),TH2O(0:LAYTWO),TLH2OS(0:LAYTWO),TH2OS(0:LAYTWO),&
     &  TLUMX(0:LAYTWO),TUMX(0:LAYTWO),TLUMXS(0:LAYTWO),TUMXS(0:LAYTWO),&
     &  TLO3(0:LAYTWO) ,TO3(0:LAYTWO) ,TLO3S(0:LAYTWO) ,TO3S(0:LAYTWO) ,&
     &  TLAUX(0:LAYTWO),TAUX(0:LAYTWO),TLAUXS(0:LAYTWO),TAUXS(0:LAYTWO),&
     &  A_H2O(0:LAYTWO),A_H2OS(0:LAYTWO),                               &
     &  A_UMX(0:LAYTWO),A_UMXS(0:LAYTWO),                               &
     &  A_O3(0:LAYTWO),A_O3S(0:LAYTWO),A_AUX(0:LAYTWO),A_AUXS(0:LAYTWO)
      SAVE IREC,N2CONT,FH2O,O2HERZ,O2_127,O2SIG,O2SIGA,O2SIGB,          &
     &  O3T0,O3T1,O3T2,CRSNO2,CRSSO2,A_CH4,WT,DWTBEG,DWTEND,            &
     &  TH2O,TUMX,TO3,TAUX,TH2OS,TUMXS,TO3S,TAUXS

!     DETERMINE ABSORPTION COEFFICIENTS IF SCALE FACTOR IS ZERO:
      IF(IFAC.EQ.-NFACMN)THEN

!         INITIALIZE SPECTRAL DATA DURING FIRST PASS:
          IF(IK.EQ.0)THEN

!             N2 CONTINUUM:
              N2CONT=N2DATA(VCEN)

!             H2O FOREIGN CONTINUUM (USE LOSCHMIDT NUMBER TO CONVERT
!             FH2O FROM CM2/MOL TO CM-1/ATM AND MULTIPLY FH2O BY 296K;
!             WITHIN THE LEVEL LOOP, MULTIPLY FH2O BY P[ATM] OVER T[K]
!             AND BY THE RADIATION TERM):
              CALL CKD(VCEN,SH2OT0,SH2OT,FH2O)
              IF(SH2OT0.GT.0.)THEN
                  SH2OT=SH2OT/SH2OT0
              ELSE
                  SH2OT=0.
              ENDIF
              SH2OT0=C_FH2O*SH2OT0
              FH2O=C_FH2O*FH2O

!             O2 HERZBERG BAND:
              O2HERZ=HERTDA(VCEN)

!             O2 1.27 MICRON CONTINUUM ABSORPTION:
              O2_127=300000*O2MATE(VCEN)/446000

!             O2 6.4 MICRON CONTINUUM ABSORPTION
              CALL O2CONT(VCEN,O2SIG,O2SIGA,O2SIGB)
              O2SIGA=O2SIG*O2SIGA
              O2SIGB=O2SIG*O2SIGB

!             DIFFUSE OZONE:
              IF(VCEN.GT.40799.)THEN
                  O3T0=.269*O3UV(VCEN)
                  O3T1=0.
                  O3T2=0.
              ELSEIF(VCEN.GT.24569.)THEN
                  O3T0=.269*O3HHT0(VCEN)
                  O3T1=O3T0*O3HHT1(VCEN)
                  O3T2=O3T0*O3HHT2(VCEN)
              ELSE
                  CALL O3CHAP(VCEN,O3T0,O3T1,O3T2)
                  O3T0=.269*O3T0
                  O3T1=O3T0*O3T1
                  O3T2=O3T0*O3T2
              ENDIF

!             NO2 CROSS-SECTIONS (14000 CM-1 TO 50000 CM-1):
              CRSNO2=NO2XS(VCEN)

!             SO2 CROSS-SECTIONS (14000 CM-1 TO 50000 CM-1):
              CRSSO2=SO2XS(VCEN)

!             CH4 CROSS-SECTIONS (400.4 TO 1050.0 NM):
              IF(VCEN.GE.10050.25 .AND. VCEN.LE.19227.07)THEN
                  A_CH4=CH4HI(VCEN)
              ELSEIF(VCEN.GE.9523.81 .AND. VCEN.LE.24975.02)THEN
                  A_CH4=CH4LO(VCEN)
              ELSE
                  A_CH4=0.
              ENDIF
          ENDIF
          P_IK=PATH_P(IK,1)
          T_IK=PATH_T(IK,1)
          DT220K=T_IK-220
          DT273K=T_IK-TZERO
          POVERT=P_IK/T_IK
          AIRATM=TZERO*POVERT
          AH2O=FH2O*POVERT
          IF(T_IK.LE.260.)THEN
              AH2O_S=SH2OT0*SH2OT
              IF(VCEN.LT.2100.)THEN
                  EXPON=EXP(-C2_260*VCEN)
                  EXPON=(1-EXPON)/(1+EXPON)
                  AH2O=AH2O*EXPON
                  AH2O_S=AH2O_S*EXPON
              ENDIF
          ELSEIF(T_IK.LT.296.)THEN
              AH2O_S=SH2OT0*SH2OT**((296-T_IK)/36)
              IF(VCEN.LT.2100.)THEN
                  EXPON=EXP(-C2*VCEN/T_IK)
                  EXPON=(1-EXPON)/(1+EXPON)
                  AH2O=AH2O*EXPON
                  AH2O_S=AH2O_S*EXPON
              ENDIF
          ELSE
              AH2O_S=SH2OT0
              IF(VCEN.LT.2100.)THEN
                  EXPON=EXP(-C2_296*VCEN)
                  EXPON=(1-EXPON)/(1+EXPON)
                  AH2O=AH2O*EXPON
                  AH2O_S=AH2O_S*EXPON
              ENDIF
          ENDIF
          AUMX=UMX_F(3,IK)*N2CONT                                       &
     &        +UMX_F(7,IK)*((1+.83*AIRATM)*O2HERZ                       &
     &                      +AIRATM*O2_127                              &
     &                      +P_IK*(O2SIG+DT220K*(O2SIGA+DT220K*O2SIGB)))&
     &        +UMX_F(9,IK)*CRSSO2                                       &
     &        +UMX_F(10,IK)*CRSNO2                                      &
     &        +UMX_F(6,IK)*A_CH4
          AO3=O3T0+DT273K*(O3T1+DT273K*O3T2)
          AAUX=0.

!         BAND MODEL LOOP:
          JT_UP=JT_MC(IK)
          JT_LO=JT_UP-1
          FT=FT_MC(IK)
          P2=P2_MC(IK)
          NORM=PATH_P(IK,1)/(2*NTLSUB)
          DO IACT=1,NACTBM
              IMOL=ABS(MACTBM(IACT))

!             CHECK IF LINE CENTER CONTRIBUTES TO ABSORPTION:
              IF(MACTBM(IACT).LT.0)THEN

!                 INTERPOLATE LINE CENTER DATA OVER TEMPERATURE:
                  ABSM=SDCN(JT_UP,IMOL)
                  ABSM=ABSM+FT*(SDCN(JT_LO,IMOL)-ABSM)
              ELSE
                  ABSM=0.
              ENDIF

!             LINE TAILS:
              TAIL=0.
              DO ITLSUB=1,NTLSUB-1
                  TAILUP=SDTL(ITLSUB,JT_UP,1,IMOL)
                  TAILUP=TAILUP+P2*(SDTL(ITLSUB,JT_UP,2,IMOL)-TAILUP)
                  TAILLO=SDTL(ITLSUB,JT_LO,1,IMOL)
                  TAILLO=TAILLO+P2*(SDTL(ITLSUB,JT_LO,2,IMOL)-TAILLO)
                  TAIL=TAIL+(TAILUP+FT*(TAILLO-TAILUP))
              ENDDO
              TAIL=2*TAIL
              TAILUP=SDTL(0,JT_UP,1,IMOL)
              TAILUP=TAILUP+P2*(SDTL(0,JT_UP,2,IMOL)-TAILUP)
              TAILLO=SDTL(0,JT_LO,1,IMOL)
              TAILLO=TAILLO+P2*(SDTL(0,JT_LO,2,IMOL)-TAILLO)
              TAIL=TAIL+(TAILUP+FT*(TAILLO-TAILUP))
              TAILUP=SDTL(NTLSUB,JT_UP,1,IMOL)
              TAILUP=TAILUP+P2*(SDTL(NTLSUB,JT_UP,2,IMOL)-TAILUP)
              TAILLO=SDTL(NTLSUB,JT_LO,1,IMOL)
              TAILLO=TAILLO+P2*(SDTL(NTLSUB,JT_LO,2,IMOL)-TAILLO)
              TAIL=TAIL+(TAILUP+FT*(TAILLO-TAILUP))
              TAIL=NORM*TAIL
              ABSM=ABSM+TAIL
              IF(IMOL.EQ.1)THEN

!                 H2O:
                  AH2O=AH2O+ABSM
              ELSEIF(IMOL.EQ.3)THEN

!                 O3:
                  AO3=AO3+ABSM
              ELSEIF(IMOL.GT.NMOLXT)THEN

!                 AUXILIARY (Y) BAND MODEL MOLECULES:
                  AAUX=AAUX+AUX_F(IMOL-NMOLXT,IK)*ABSM
              ELSE

!                 UNIFORMLY MIXED GASES:
                  AUMX=AUMX+UMX_F(IMOL,IK)*ABSM
              ENDIF
          ENDDO

!         CROSS-SECTION LOOP:
          JT_UP=KT_MC(IK)
          JT_LO=JT_UP-1
          FT=GT_MC(IK)
          DO IACT=1,NACTX
              IMOL=MACTX(IACT)

!             INTERPOLATE OVER TEMPERATURE:
              ABSM=SDTL(0,JT_UP,1,IMOL)
              ABSM=ABSM+FT*(SDTL(0,JT_LO,1,IMOL)-ABSM)
              IF(IMOL.GT.NMOLXT)THEN

!                 AUXILIARY (Y) BAND MODEL MOLECULES:
                  AAUX=AAUX+AUX_F(IMOL-NMOLXT,IK)*ABSM
              ELSE

!                 UNIFORMLY MIXED GASES:
                  AUMX=AUMX+UMX_F(IMOL,IK)*ABSM
              ENDIF
          ENDDO
      ENDIF

!     BAND MODEL (CM-1) OR WAVELENGTH OUTPUT?
      IF(LBMRES)THEN

!         BAND MODEL RESOLUTION:
          IF(IFAC.EQ.-NFACMN)THEN
              A_H2O(IK)=AH2O
              A_UMX(IK)=AUMX
              A_O3(IK)=AO3
              A_AUX(IK)=AAUX
              IF(IK.LT.IKMAX)RETURN

!             WRITE ABSORPTION COEFFICIENTS INSTEAD OF TRANSMITTANCES:
              NMWAVE=NMWAVE+1
              WRITE(IDBOUT,REC=NMWAVE)(1.,I=1,4*(IKMAX+1)),             &
     &          (A_H2O(I),I=0,IKMAX),(A_UMX(I),I=0,IKMAX),              &
     &          (A_O3(I) ,I=0,IKMAX),(A_AUX(I),I=0,IKMAX)
              RETURN
          ENDIF

!         STORE TRANSMITTANCES:
          TLH2O(IK)=TL_H2O
          TLUMX(IK)=TL_UMX
          TLO3(IK)=TL_O3
          TLAUX(IK)=TL_AUX
          TH2O(IK)=T_H2O
          TUMX(IK)=T_UMX
          TO3(IK)=T_O3
          TAUX(IK)=T_AUX
          IF(IK.LT.IKMAX)RETURN

!         WRITE L-SHAPED AND DIRECT PATH TRANSMITTANCES:
          NMWAVE=NMWAVE+1
          WRITE(IDBOUT,REC=NMWAVE)                                      &
     &      (TLH2O(I),I=0,IKMAX),(TLUMX(I),I=0,IKMAX),                  &
     &      (TLO3(I) ,I=0,IKMAX),(TLAUX(I),I=0,IKMAX),                  &
     &      (TH2O(I) ,I=0,IKMAX),(TUMX(I) ,I=0,IKMAX),                  &
     &      (TO3(I)  ,I=0,IKMAX),(TAUX(I) ,I=0,IKMAX)

!         WRITE THE SCALED BY ONE TRANSMITTANCE DATA:
          IF(IFAC.NE.0)RETURN
          WRITE(IPR1,'(/F8.2,A,I4,11F9.6:/(16F9.6))')VCEN,' CM-1:'//    &
     &      ' H2O  (L-SHAPE) - RECORD #',NMWAVE,(TLH2O(I),I=0,IKMAX)
          WRITE(IPR1,'( F8.2,A,I4,11F9.6:/(16F9.6))')VCEN,' CM-1:'//    &
     &      ' UMIX (L-SHAPE) - RECORD #',NMWAVE,(TLUMX(I),I=0,IKMAX)
          WRITE(IPR1,'( F8.2,A,I4,11F9.6:/(16F9.6))')VCEN,' CM-1:'//    &
     &      ' O3   (L-SHAPE) - RECORD #',NMWAVE,(TLO3(I),I=0,IKMAX)
          WRITE(IPR1,'( F8.2,A,I4,11F9.6:/(16F9.6))')VCEN,' CM-1:'//    &
     &      ' AUX  (L-SHAPE) - RECORD #',NMWAVE,(TLAUX(I),I=0,IKMAX)
          WRITE(IPR1,'( F8.2,A,I4,11F9.6:/(16F9.6))')VCEN,' CM-1:'//    &
     &      ' H2O   (DIRECT) - RECORD #',NMWAVE,(TH2O(I),I=0,IKMAX)
          WRITE(IPR1,'( F8.2,A,I4,11F9.6:/(16F9.6))')VCEN,' CM-1:'//    &
     &      ' UMIX  (DIRECT) - RECORD #',NMWAVE,(TUMX(I),I=0,IKMAX)
          WRITE(IPR1,'( F8.2,A,I4,11F9.6:/(16F9.6))')VCEN,' CM-1:'//    &
     &      ' O3    (DIRECT) - RECORD #',NMWAVE,(TO3(I),I=0,IKMAX)
          WRITE(IPR1,'( F8.2,A,I4,11F9.6:/(16F9.6))')VCEN,' CM-1:'//    &
     &      ' AUX   (DIRECT) - RECORD #',NMWAVE,(TAUX(I),I=0,IKMAX)
          RETURN
      ENDIF
      NM=NINT(1.E7/(VCEN+HBNDWD))
      IF(NMWAVE.EQ.NM)THEN

!         BAND MODEL DATA COMPLETELY WITHIN NM BIN.  INCREMENT & RETURN:
          IF(IK.EQ.IKMAX)WT=WT+1
          IF(IFAC.EQ.-NFACMN)THEN
              A_H2O(IK)=A_H2O(IK)+AH2O
              A_UMX(IK)=A_UMX(IK)+AUMX
              A_O3(IK)=A_O3(IK)+AO3
              A_AUX(IK)=A_AUX(IK)+AAUX
              RETURN
          ENDIF
          TLH2O(IK)=TLH2O(IK)+TL_H2O
          TLUMX(IK)=TLUMX(IK)+TL_UMX
          TLO3(IK)=TLO3(IK)+TL_O3
          TLAUX(IK)=TLAUX(IK)+TL_AUX
          TH2O(IK)=TH2O(IK)+T_H2O
          TUMX(IK)=TUMX(IK)+T_UMX
          TO3(IK)=TO3(IK)+T_O3
          TAUX(IK)=TAUX(IK)+T_AUX
          RETURN

!     IS THIS THE INITIATION OF THE FIRST NM BIN (I.E., DOES NMWAVE=0)?
      ELSEIF(NMWAVE.EQ.0)THEN

!         RETURN IF THE NM BIN DOESN'T BEGIN WITHIN THIS BAND MODEL BIN.
          IF(NINT(1.E7/(VCEN-HBNDWD)).EQ.NM)RETURN

!         WHAT FRACTION OF THE BAND MODEL BIN OVERLAPS WITH THE NM BIN?
          IF(IK.EQ.0)THEN

!             FIRST LEVEL:
              DWTBEG=(VCEN-1.E7/(NM+.5))/BNDWID+.5
              WT=DWTBEG
          ELSEIF(IK.EQ.IKMAX)THEN

!             LAST LEVEL:
              IREC=0
              NMWAVE=NM
          ENDIF
      ELSEIF(NM+1.NE.NMWAVE)THEN

!         A NM BIN HAS BEEN COMPLETELY SKIPPED:
          STOP 'Error in MCMOL:  Band model data is too coarse.'
      ELSE

!         WHAT FRACTION OF THE BAND MODEL BIN OVERLAPS WITH OLD NM BIN?
          IF(IK.EQ.0)THEN

!             FIRST LEVEL:
              DWTBEG=(VCEN-1.E7/(NM+.5))/BNDWID+.5
              DWTEND=1-DWTBEG
              WT=WT+DWTEND
          ENDIF
          IF(IFAC.EQ.-NFACMN)THEN
              A_H2OS(IK)=(A_H2O(IK)+DWTEND*AH2O)/WT
              A_UMXS(IK)=(A_UMX(IK)+DWTEND*AUMX)/WT
              A_O3S(IK)=(A_O3(IK)+DWTEND*AO3)/WT
              A_AUXS(IK)=(A_AUX(IK)+DWTEND*AAUX)/WT
              IF(IK.EQ.IKMAX)THEN

!                 LAST LEVEL:
                  IREC=IREC+1
                  WRITE(IDBOUT,REC=IREC)(1.,I=1,4*(IKMAX+1)),           &
     &              (A_H2OS(I),I=0,IKMAX),(A_UMXS(I),I=0,IKMAX),        &
     &              (A_O3S(I) ,I=0,IKMAX),(A_AUXS(I),I=0,IKMAX)
                  WT=DWTBEG
                  NMWAVE=NM
              ENDIF
          ELSE
              TLH2OS(IK)=(TLH2O(IK)+DWTEND*TL_H2O)/WT
              TLUMXS(IK)=(TLUMX(IK)+DWTEND*TL_UMX)/WT
              TLO3S(IK)=(TLO3(IK)+DWTEND*TL_O3)/WT
              TLAUXS(IK)=(TLAUX(IK)+DWTEND*TL_AUX)/WT
              TH2OS(IK)=(TH2O(IK)+DWTEND*T_H2O)/WT
              TUMXS(IK)=(TUMX(IK)+DWTEND*T_UMX)/WT
              TO3S(IK)=(TO3(IK)+DWTEND*T_O3)/WT
              TAUXS(IK)=(TAUX(IK)+DWTEND*T_AUX)/WT
              IF(IK.EQ.IKMAX)THEN

!                 LAST LEVEL:
                  IREC=IREC+1
                  WRITE(IDBOUT,REC=IREC)                                &
     &              (TLH2OS(I),I=0,IKMAX),(TLUMXS(I),I=0,IKMAX),        &
     &              (TLO3S(I),I=0,IKMAX),(TLAUXS(I),I=0,IKMAX),         &
     &              (TH2OS(I),I=0,IKMAX),(TUMXS(I),I=0,IKMAX),          &
     &              (TO3S(I),I=0,IKMAX),(TAUXS(I),I=0,IKMAX)
                  IF(IFAC.EQ.0)THEN
                      IF(BINOUT)THEN
                          CALL BNWT4(102,IKMAX,NMWAVE,IREC,TLH2OS)
                          CALL BNWT4(103,IKMAX,NMWAVE,IREC,TLUMXS)
                          CALL BNWT4(104,IKMAX,NMWAVE,IREC,TLO3S)
                          CALL BNWT4(105,IKMAX,NMWAVE,IREC,TH2OS)
                          CALL BNWT4(106,IKMAX,NMWAVE,IREC,TUMXS)
                          CALL BNWT4(107,IKMAX,NMWAVE,IREC,TO3S)
                          CALL BNWT4(108,IKMAX,NMWAVE,IREC,TLAUXS)
                          CALL BNWT4(109,IKMAX,NMWAVE,IREC,TAUXS)
                      ELSE

!                         WRITE THE SCALED BY ONE TRANSMITTANCE DATA:
                          WRITE(IPR1,'(/I8,A,I4,11F9.6:/(16F9.6))')     &
     &                      NMWAVE,' NM -  H2O  (L-SHAPE) -  RECORD #', &
     &                      IREC,(TLH2OS(I),I=0,IKMAX)
                          WRITE(IPR1,'(I8,A,I4,11F9.6:/(16F9.6))')      &
     &                      NMWAVE,' NM -  UMIX (L-SHAPE) -  RECORD #', &
     &                      IREC,(TLUMXS(I),I=0,IKMAX)
                          WRITE(IPR1,'(I8,A,I4,11F9.6:/(16F9.6))')      &
     &                      NMWAVE,' NM -  O3   (L-SHAPE) -  RECORD #', &
     &                      IREC,(TLO3S(I),I=0,IKMAX)
                          WRITE(IPR1,'(I8,A,I4,11F9.6:/(16F9.6))')      &
     &                      NMWAVE,' NM -  AUX  (L-SHAPE) -  RECORD #', &
     &                      IREC,(TLAUXS(I),I=0,IKMAX)
                          WRITE(IPR1,'(I8,A,I4,11F9.6:/(16F9.6))')      &
     &                      NMWAVE,' NM -  H2O   (DIRECT) -  RECORD #', &
     &                      IREC,(TH2OS(I),I=0,IKMAX)
                          WRITE(IPR1,'(I8,A,I4,11F9.6:/(16F9.6))')      &
     &                      NMWAVE,' NM -  UMIX  (DIRECT) -  RECORD #', &
     &                      IREC,(TUMXS(I),I=0,IKMAX)
                          WRITE(IPR1,'(I8,A,I4,11F9.6:/(16F9.6))')      &
     &                      NMWAVE,' NM -  O3    (DIRECT) -  RECORD #', &
     &                      IREC,(TO3S(I),I=0,IKMAX)
                          WRITE(IPR1,'(I8,A,I4,11F9.6:/(16F9.6))')      &
     &                      NMWAVE,' NM -  AUX   (DIRECT) -  RECORD #', &
     &                      IREC,(TAUXS(I),I=0,IKMAX)
                      ENDIF
                  ENDIF
                  WT=DWTBEG
                  NMWAVE=NM
              ENDIF
          ENDIF
      ENDIF
      IF(IFAC.EQ.-NFACMN)THEN
          A_H2O(IK)=DWTBEG*AH2O
          A_UMX(IK)=DWTBEG*AUMX
          A_O3(IK)=DWTBEG*AO3
          A_AUX(IK)=DWTBEG*AAUX
          RETURN
      ENDIF
      TLH2O(IK)=DWTBEG*TL_H2O
      TLUMX(IK)=DWTBEG*TL_UMX
      TLO3(IK)=DWTBEG*TL_O3
      TLAUX(IK)=DWTBEG*TL_AUX
      TH2O(IK)=DWTBEG*T_H2O
      TUMX(IK)=DWTBEG*T_UMX
      TO3(IK)=DWTBEG*T_O3
      TAUX(IK)=DWTBEG*T_AUX
      RETURN
      END
