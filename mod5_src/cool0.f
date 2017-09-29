      SUBROUTINE COOL0(IRPT)
      IMPLICIT NONE

!     THIS ROUTINE DEFINES THE /COOLRT/ ARRAYS WHICH CONTAIN THE
!     FREQUENCY INDEPENDENT DATA FOR THE COOLING RATE CALCULATIONS.
!     COOLING RATES, THE CHANGE IN TEMPERATURE T WITH TIME t, EQUAL
!     THE GRAVITATIONAL CONSTANT g OVER THE SPECIFIC HEAT OF AIR
!     AT CONSTANT PRESSURE Cp TIMES THE CHANGE IN NET UPWARD FLUX F
!     WITH PRESSURE P:

!            dT     g   dF     g   1      dF
!            --  =  --  --  =  --  -  ----------
!            dt     Cp  dP     Cp  P  d(ln P/Po)

!     Po IS DEFINED AS 1 MILLIBAR.  THE DERIVATIVE OF THE
!     NET FLUX WITH LOG PRESSURE IS DEFINED AT EACH LAYER BOUNDARY BY
!     FITTING THE VALUES OF F AND ln P/Po FROM TWO LAYER BOUNDARIES
!     BELOW TO TWO LAYER BOUNDARIES ABOVE THE CURRENT LAYER BOUNDARY
!     TO A POLYNOMIAL (FOURTH DEGREE EXCEPT AT THE TOP AND BOTTOM
!     OF THE ATMOSPHERE), AND ANALYTICALLY COMPUTING THE DERIVATIVE.
!     IN THIS ROUTINE, THE PRESSURE DEPENDENT FITTING COEFFICIENTS,
!     P5COEF, ARE CALCULATED.  THESE TERMS, UNLIKE THE NET FLUXES,
!     ARE FREQUENCY INDEPENDENT AND CAN BE PRE-CALCULATED.

!     ARGUMENTS:
!       IRPT    REPEAT CALCULATION INDEX.
      INTEGER IRPT

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /COOLRT/
!       ISPCCR  UNIT NUMBER FOR SPECTRAL COOLING RATE FILE.
!       P3COEF  PRESSURE DEPENDENT COEFFICIENTS FOR CALCULATING COOLING
!               RATES WITH 3 POINT SPLINE FIT [(K/DAY) / (W/CM2)].
!       P5COEF  PRESSURE DEPENDENT COEFFICIENTS FOR CALCULATING COOLING
!               RATES WITH 5 POINT SPLINE FIT [(K/DAY) / (W/CM2)].
!       CLRTSM  LAYER BOUNDARY IN-BAND COOLING RATES [K/DAY].
!       CLRTS0  LAYER BOUNDARY IN-BAND COOLING RATES
!               EXCLUDING DIRECT SOLAR CONTRIBUTION [K/DAY].
!       NTFSM   LAYER BOUNDARY IN-BAND NET (THERMAL PLUS SCATTERED
!               SOLAR PLUS DIRECT SOLAR) UPWARD FLUX [W/CM2].
!       UPFSM   LAYER BOUNDARY IN-BAND UPWARD THERMAL (PLUS SCATTERED
!               SOLAR IF DISORT & NO AZIMUTH DEPENDENCE) FLUX [W/CM2].
!       DNFSM   LAYER BOUNDARY IN-BAND DOWNWARD THERMAL (PLUS SCATTERED
!               SOLAR IF DISORT & NO AZIMUTH DEPENDENCE) FLUX [W/CM2].
!       UPFSSM  LAYER BOUNDARY IN-BAND UPWARD SCATTERED SOLAR FLUX
!               (USED WITH DISORT ONLY IF AZIMUTH DEPENDENT) [W/CM2].
!       DNFSSM  LAYER BOUNDARY IN-BAND DOWNWARD SCATTERED SOLAR FLUX
!               (USED WITH DISORT ONLY IF AZIMUTH DEPENDENT) [W/CM2].
      DOUBLE PRECISION CLRTSM,CLRTS0,NTFSM,UPFSM,DNFSM,UPFSSM,DNFSSM
      REAL P3COEF,P5COEF
      INTEGER ISPCCR
      COMMON/COOLRT/CLRTSM(LAYDIM),CLRTS0(LAYDIM),NTFSM(LAYDIM),        &
     &  UPFSM(LAYDIM),DNFSM(LAYDIM),UPFSSM(LAYDIM),DNFSSM(LAYDIM),      &
     &  P3COEF(3,LAYDIM),P5COEF(5,LAYDIM),ISPCCR

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

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:
      CHARACTER*1 PFRMT

!     LOCAL VARIABLES:
!       Z1,Z2,Z3,Z4,Z5   LOG PRESSURE VALUES USED IN FIT.
!       LM2,L,LP2        BOTTOM, CURRENT & TOP LEVEL INDEX USED IN FIT.
!       COEF             g OVER P TIMES Cp [(K/day) / (W/cm2)].
!       LABEL            PRESSURE LABEL FOR HEADER TO COOLING RATE FILE.
!       I                INDEX FOR WRITING OUT PRESSURES.
!       ILO              CURRENT LOWER LIMIT FOR I.
!       IPOINT           POINTER FOR FORMAT STRING.
      INTEGER LM2,L,LP2,I,ILO,IPOINT
      REAL Z1,Z2,Z3,Z4,Z5,COEF
      CHARACTER*8 LABEL

!     DATA:
!       GCP     THE GRAVITATIONAL CONSTANT g OVER THE SPECIFIC HEAT OF
!               AIR AT CONSTANT PRESSURE Cp [84220 mb K cm2 day-1 W-1].
!       FRMT    FORMAT FOR OUTPUTING PRESSURES
      REAL GCP
      CHARACTER*103 FRMT
      DATA GCP/84220./,FRMT(1:49),FRMT(50:103)/                         &
     &  '(A8,1X,F10.9,1X,F10.9,1X,F10.9,1X,F10.9,1X,F10.9,',            &
     &  '1X,F10.9,1X,F10.9,1X,F10.9,1X,F10.9,1X,F10.9,1X,F10.9)'/

!     CHECK IF OLD ATMOSPHERE IS BEING RE-USED.
      IF(IRPT.EQ.3)THEN

!         REPEAT USE OF ATMOSPHERE.
!         RE-INITIALIZE BAND PASS COOLING RATE ARRAY AND RETURN.
          DO L=1,ML
              CLRTSM(L)=0.D0
              CLRTS0(L)=0.D0
              NTFSM(L)=0.D0
              UPFSM(L)=0.D0
              DNFSM(L)=0.D0
              UPFSSM(L)=0.D0
              DNFSSM(L)=0.D0
          ENDDO
          RETURN
      ENDIF

!     CHECK LAYER BOUNDARY NUMBER.
      IF(ML.LT.4)THEN
          WRITE(IPR,'(/2A)')' WARNING in COOL0:  THERE MUST BE AT',     &
     &      ' LEAST 4 LAYER BOUNDARIES FOR COOLING RATE CALCULATIONS.'
          RETURN
      ENDIF

!     ASSIGN UNIT NUMBER FOR SPECTRAL COOLING RATE FILE.
      ISPCCR=ICR
      ILO=1
      LABEL='  P(MB):'

!     DETERMINE FITTING COEFFICIENTS FOR THE FIRST LAYER.
      IPOINT=12
      FRMT(IPOINT:IPOINT)=PFRMT(PM(1))
      IPOINT=IPOINT+9
      FRMT(IPOINT:IPOINT)=PFRMT(PM(2))
      IPOINT=IPOINT+9
      FRMT(IPOINT:IPOINT)=PFRMT(PM(3))
      CLRTSM(1)=0.D0
      CLRTS0(1)=0.D0
      NTFSM(1)=0.D0
      UPFSM(1)=0.D0
      DNFSM(1)=0.D0
      UPFSSM(1)=0.D0
      DNFSSM(1)=0.D0
      CLRTSM(2)=0.D0
      CLRTS0(2)=0.D0
      NTFSM(2)=0.D0
      UPFSM(2)=0.D0
      DNFSM(2)=0.D0
      UPFSSM(2)=0.D0
      DNFSSM(2)=0.D0
      CLRTSM(3)=0.D0
      CLRTS0(3)=0.D0
      NTFSM(3)=0.D0
      UPFSM(3)=0.D0
      DNFSM(3)=0.D0
      UPFSSM(3)=0.D0
      DNFSSM(3)=0.D0
      Z1=LOG(PM(1))
      Z2=LOG(PM(2))
      Z3=LOG(PM(3))
      COEF=GCP/PM(1)
      P5COEF(1,1)=0.
      P5COEF(2,1)=0.
      P5COEF(3,1)=COEF/(Z1-Z2)+COEF/(Z1-Z3)
      P5COEF(4,1)=COEF*(Z1-Z3)/((Z2-Z1)*(Z2-Z3))
      P5COEF(5,1)=COEF*(Z1-Z2)/((Z3-Z1)*(Z3-Z2))

!     DETERMINE FITTING COEFFICIENTS FOR THE SECOND LAYER.
      IPOINT=IPOINT+9
      FRMT(IPOINT:IPOINT)=PFRMT(PM(4))
      CLRTSM(4)=0.D0
      CLRTS0(4)=0.D0
      NTFSM(4)=0.D0
      UPFSM(4)=0.D0
      DNFSM(4)=0.D0
      UPFSSM(4)=0.D0
      DNFSSM(4)=0.D0
      Z4=LOG(PM(4))
      COEF=GCP/PM(2)
      P3COEF(1,2)=COEF*(Z2-Z3)/((Z1-Z3)*(Z1-Z2))
      P3COEF(2,2)=COEF/(Z2-Z1)+COEF/(Z2-Z3)
      P3COEF(3,2)=COEF*(Z2-Z1)/((Z3-Z1)*(Z3-Z2))
      P5COEF(1,2)=0.
      P5COEF(2,2)=COEF*(Z2-Z3)*(Z2-Z4)/((Z1-Z2)*(Z1-Z3)*(Z1-Z4))
      P5COEF(3,2)=COEF/(Z2-Z1)+COEF/(Z2-Z3)+COEF/(Z2-Z4)
      P5COEF(4,2)=COEF*(Z2-Z1)*(Z2-Z4)/((Z3-Z1)*(Z3-Z2)*(Z3-Z4))
      P5COEF(5,2)=COEF*(Z2-Z1)*(Z2-Z3)/((Z4-Z1)*(Z4-Z2)*(Z4-Z3))

!     LOOP OVER INTERMEDIATE LAYERS
      LM2=1
      L=3
      DO LP2=5,ML

!         DETERMINE FITTING COEFFICIENTS.
          IPOINT=IPOINT+9
          FRMT(IPOINT:IPOINT)=PFRMT(PM(LP2))
          IF(IPOINT.GE.102)THEN

!             WRITE OUT PRESSURE HEADER FOR COOLING RATES FILE.
              IF(.NOT.LJMASS .AND. NPR.LE.-2)                           &
     &          WRITE(ISPCCR,FMT=FRMT)LABEL,(DBLE(PM(I)),I=ILO,LP2)
              ILO=LP2+1
              LABEL='      '
              IPOINT=3
          ENDIF
          CLRTSM(LP2)=0.D0
          CLRTS0(LP2)=0.D0
          NTFSM(LP2)=0.D0
          UPFSM(LP2)=0.D0
          DNFSM(LP2)=0.D0
          UPFSSM(LP2)=0.D0
          DNFSSM(LP2)=0.D0
          Z5=LOG(PM(LP2))
          COEF=GCP/PM(L)
          P3COEF(1,L)=COEF*(Z3-Z4)/((Z2-Z4)*(Z2-Z3))
          P3COEF(2,L)=COEF/(Z3-Z2)+COEF/(Z3-Z4)
          P3COEF(3,L)=COEF*(Z3-Z2)/((Z4-Z2)*(Z4-Z3))
          P5COEF(1,L)=COEF*(Z3-Z2)*(Z3-Z4)*(Z3-Z5)                      &
     &           /((Z1-Z3)*(Z1-Z2)*(Z1-Z4)*(Z1-Z5))
          P5COEF(2,L)=COEF*(Z3-Z1)*(Z3-Z4)*(Z3-Z5)                      &
     &           /((Z2-Z3)*(Z2-Z1)*(Z2-Z4)*(Z2-Z5))
          P5COEF(3,L)=COEF/(Z3-Z1)+COEF/(Z3-Z2)                         &
     &               +COEF/(Z3-Z4)+COEF/(Z3-Z5)
          P5COEF(4,L)=COEF*(Z3-Z1)*(Z3-Z2)*(Z3-Z5)                      &
     &           /((Z4-Z3)*(Z4-Z1)*(Z4-Z2)*(Z4-Z5))
          P5COEF(5,L)=COEF*(Z3-Z1)*(Z3-Z2)*(Z3-Z4)                      &
     &           /((Z5-Z3)*(Z5-Z1)*(Z5-Z2)*(Z5-Z4))
          LM2=LM2+1
          L=L+1
          Z1=Z2
          Z2=Z3
          Z3=Z4
          Z4=Z5
      ENDDO

!     DETERMINE FITTING COEFFICIENTS FOR THE SECOND TO LAST LAYER.
      COEF=GCP/PM(L)
      P3COEF(1,L)=COEF*(Z3-Z4)/((Z2-Z4)*(Z2-Z3))
      P3COEF(2,L)=COEF/(Z3-Z2)+COEF/(Z3-Z4)
      P3COEF(3,L)=COEF*(Z3-Z2)/((Z4-Z2)*(Z4-Z3))
      P5COEF(1,L)=COEF*(Z3-Z2)*(Z3-Z4)/((Z1-Z2)*(Z1-Z3)*(Z1-Z4))
      P5COEF(2,L)=COEF*(Z3-Z1)*(Z3-Z4)/((Z2-Z1)*(Z2-Z3)*(Z2-Z4))
      P5COEF(3,L)=COEF/(Z3-Z1)+COEF/(Z3-Z2)+COEF/(Z3-Z4)
      P5COEF(4,L)=COEF*(Z3-Z1)*(Z3-Z2)/((Z4-Z1)*(Z4-Z2)*(Z4-Z3))
      P5COEF(5,L)=0.

!     DETERMINE FITTING COEFFICIENTS FOR THE LAST LAYER.
      COEF=GCP/PM(ML)
      P5COEF(1,ML)=COEF*(Z4-Z3)/((Z2-Z3)*(Z2-Z4))
      P5COEF(2,ML)=COEF*(Z4-Z2)/((Z3-Z2)*(Z3-Z4))
      P5COEF(3,ML)=COEF/(Z4-Z2)+COEF/(Z4-Z3)
      P5COEF(4,ML)=0.
      P5COEF(5,ML)=0.

!     WRITE OUT PRESSURE HEADER FOR COOLING RATES FILE.
      IF(NPR.GT.-2)RETURN
      IF(.NOT.LJMASS) THEN
         IF(ILO.LE.ML)WRITE(ISPCCR,FMT=FRMT)LABEL,(DBLE(PM(I)),I=ILO,ML)
         WRITE(ISPCCR,'(/(A))')'    FREQ   SPECTRAL COOLING RATES',     &
     &                         '    CM-1   K DAY-1 / CM-1'
      ENDIF
      RETURN
      END
      CHARACTER*1 FUNCTION PFRMT(PMVAL)

!     FUNCTION USED TO DEFINE FORMAT FOR WRITING PM
      REAL PMVAL
      DOUBLE PRECISION DPMVAL
      DPMVAL=DBLE(PMVAL)
      IF(DPMVAL.GE.999.9999995D0)THEN
          PFRMT='5'
      ELSEIF(DPMVAL.GE.99.99999995D0)THEN
          PFRMT='6'
      ELSEIF(DPMVAL.GE.9.999999995D0)THEN
          PFRMT='7'
      ELSEIF(DPMVAL.GE..9999999995D0)THEN
          PFRMT='8'
      ELSE
          PFRMT='9'
      ENDIF
      RETURN
      END