      SUBROUTINE COOL(IV)

!     THIS ROUTINE CALCULATES AND WRITES OUT SPECTRAL COOLING RATES.
!     COOLING RATES, THE CHANGE IN TEMPERATURE T WITH TIME t, EQUAL
!     THE GRAVITATIONAL CONSTANT g OVER THE SPECIFIC HEAT OF AIR AT
!     CONSTANT PRESSURE Cp TIMES THE CHANGE IN NET UPWARD FLUX F
!     WITH PRESSURE P:

!            dT     g   dF     g   1      dF
!            --  =  --  --  =  --  -  ----------
!            dt     Cp  dP     Cp  P  d(ln P/Po)

!     Po IS STANDARD PRESSURE (Po = 1 ATM).  THE DERIVATIVE OF THE
!     NET FLUX WITH LOG PRESSURE IS DEFINED AT EACH LAYER BOUNDARY BY
!     FITTING THE VALUES OF F AND ln P/Po FROM TWO LAYER BOUNDARIES
!     BELOW TO TWO LAYER BOUNDARIES ABOVE THE CURRENT LAYER BOUNDARY
!     TO A POLYNOMIAL (FOURTH DEGREE EXCEPT AT THE TOP AND BOTTOM
!     OF THE ATMOSPHERE), AND ANALYTICALLY COMPUTING THE DERIVATIVE.

!     DECLARE ARGUMENTS
!       IV      SPECTRAL FREQUENCY [CM-1]
      INTEGER IV

!     INCLUDE PARAMETERS:
      INCLUDE 'PARAMS.h'

!     LIST COMMONS:
!       ISPCCR  UNIT NUMBER FOR SPECTRAL COOLING RATE FILE.
!       PCOEF   PRESSURE DEPENDENT COEFFICIENTS USED TO CALCULATE.
!               COOLING RATES FROM NET FLUXES [(K/DAY) / (W/CM2)].
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
      INTEGER ISPCCR
      REAL PCOEF,CLRTSM,CLRTS0,NTFSM,UPFSM,DNFSM,UPFSSM,DNFSSM
      COMMON/COOLRT/ISPCCR,PCOEF(5,LAYDIM),CLRTSM(LAYDIM),              &
     &  CLRTS0(LAYDIM),NTFSM(LAYDIM),UPFSM(LAYDIM),                     &
     &  DNFSM(LAYDIM),UPFSSM(LAYDIM),DNFSSM(LAYDIM)

!       SUBINT   SPECTRAL BIN "K" SUB-INTERVAL FRACTIONAL WIDTHS.
!       UPFLX    LAYER BOUNDARY UPWARD THERMAL SPECTRAL
!                FLUX (INCLUDES SCATTERED SOLAR IF DISORT
!                & NO AZIMUTH DEPENDENCE) [W CM-2 / CM-1].
!       DNFLX    LAYER BOUNDARY DOWNWARD THERMAL SPECTRAL
!                FLUX (INCLUDES SCATTERED SOLAR IF DISORT
!                & NO AZIMUTH DEPENDENCE) [W CM-2 / CM-1].
!       UPFLXS   LAYER BOUNDARY UPWARD SCATTERED SOLAR
!                SPECTRAL FLUX (USED WITH DISORT ONLY
!                IF AZIMUTH DEPENDENT) [W CM-2 / CM-1].
!       DNFLXS   LAYER BOUNDARY DOWNWARD SCATTERED SOLAR
!                SPECTRAL FLUX (USED WITH DISORT ONLY
!                IF AZIMUTH DEPENDENT) [W CM-2 / CM-1].
!       NTFLX    LAYER BOUNDARY NET (THERMAL PLUS SCATTERED SOLAR PLUS
!                DIRECT SOLAR) UPWARD SPECTRAL FLUX [W CM-2 / CM-1].
      REAL SUBINT,UPFLX,DNFLX,UPFLXS,DNFLXS,NTFLX
      COMMON/NETFLX/SUBINT(MXKSUB),UPFLX(LAYDIM),DNFLX(LAYDIM),         &
     &  UPFLXS(LAYDIM),DNFLXS(LAYDIM),NTFLX(LAYDIM)

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     DEFINE LOCAL VARIABLES:
!       LM2,LM1,L,LP1,LP2   LAYER BOUNDARY INDICES.
!       CLRATE              SPECTRAL COOLING RATES [K DAY-1 / CM-1].
!       NTFLX0              LAYER BOUNDARY NET UPWARD FLUX EXCLUDING
!                           DIRECT SOLAR [W CM-2 / CM-1]
      INTEGER LM2,LM1,L,LP1,LP2
      REAL CLRATE(LAYDIM),NTFLX0(LAYDIM)

!     GROUND SPECTRAL COOLING RATE:
      CLRATE(1)=PCOEF(3,1)*NTFLX(1)+PCOEF(4,1)*NTFLX(2)                 &
     &  +PCOEF(5,1)*NTFLX(3)
      CLRTSM(1)=CLRTSM(1)+CLRATE(1)
      UPFSM(1)=UPFSM(1)+UPFLX(1)
      DNFSM(1)=DNFSM(1)+DNFLX(1)
      UPFSSM(1)=UPFSSM(1)+UPFLXS(1)
      DNFSSM(1)=DNFSSM(1)+DNFLXS(1)
      NTFSM(1)=NTFSM(1)+NTFLX(1)

!     COOLING RATE WITH DIRECT SOLAR EXCLUDED:
      NTFLX0(1)=(UPFLX(1)-DNFLX(1))+(UPFLXS(1)-DNFLXS(1))
      NTFLX0(2)=(UPFLX(2)-DNFLX(2))+(UPFLXS(2)-DNFLXS(2))
      NTFLX0(3)=(UPFLX(3)-DNFLX(3))+(UPFLXS(3)-DNFLXS(3))
      CLRTS0(1)=CLRTS0(1)+(PCOEF(3,1)*NTFLX0(1)                         &
     &  +PCOEF(4,1)*NTFLX0(2)+PCOEF(5,1)*NTFLX0(3))

!     SECOND LAYER BOUNDARY SPECTRAL COOLING RATE:
      CLRATE(2)=PCOEF(2,2)*NTFLX(1)+PCOEF(3,2)*NTFLX(2)                 &
     &  +PCOEF(4,2)*NTFLX(3)+PCOEF(5,2)*NTFLX(4)
      CLRTSM(2)=CLRTSM(2)+CLRATE(2)
      UPFSM(2)=UPFSM(2)+UPFLX(2)
      DNFSM(2)=DNFSM(2)+DNFLX(2)
      UPFSSM(2)=UPFSSM(2)+UPFLXS(2)
      DNFSSM(2)=DNFSSM(2)+DNFLXS(2)
      NTFSM(2)=NTFSM(2)+NTFLX(2)

!     COOLING RATE WITH DIRECT SOLAR EXCLUDED:
      NTFLX0(4)=(UPFLX(4)-DNFLX(4))+(UPFLXS(4)-DNFLXS(4))
      CLRTS0(2)=CLRTS0(2)+(PCOEF(2,2)*NTFLX0(1)+PCOEF(3,2)*NTFLX0(2)    &
     &  +PCOEF(4,2)*NTFLX0(3)+PCOEF(5,2)*NTFLX0(4))

!     LOOP OVER INTERMEDIATE LAYERS:
      LM2=1
      LM1=2
      L=3
      LP1=4
      DO 10 LP2=5,ML

!         INTERMEDIATE LAYER BOUNDARY SPECTRAL COOLING RATES.
          CLRATE(L)=PCOEF(1,L)*NTFLX(LM2)+PCOEF(2,L)*NTFLX(LM1)         &
     &     +PCOEF(3,L)*NTFLX(L)+PCOEF(4,L)*NTFLX(LP1)                   &
     &     +PCOEF(5,L)*NTFLX(LP2)
          CLRTSM(L)=CLRTSM(L)+CLRATE(L)
          UPFSM(L)=UPFSM(L)+UPFLX(L)
          DNFSM(L)=DNFSM(L)+DNFLX(L)
          UPFSSM(L)=UPFSSM(L)+UPFLXS(L)
          DNFSSM(L)=DNFSSM(L)+DNFLXS(L)
          NTFSM(L)=NTFSM(L)+NTFLX(L)

!         COOLING RATE WITH DIRECT SOLAR EXCLUDED:
          NTFLX0(LP2)=(UPFLX(LP2)-DNFLX(LP2))+(UPFLXS(LP2)-DNFLXS(LP2))
          CLRTS0(L)=CLRTS0(L)+(PCOEF(1,L)*NTFLX0(LM2)                   &
     &      +PCOEF(2,L)*NTFLX0(LM1)+PCOEF(3,L)*NTFLX0(L)                &
     &      +PCOEF(4,L)*NTFLX0(LP1)+PCOEF(5,L)*NTFLX0(LP2))
          LM2=LM1
          LM1=L
          L=LP1
          LP1=LP2
   10 CONTINUE

!     SECOND TO LAST LAYER BOUNDARY SPECTRAL COOLING RATE.
      CLRATE(L)=PCOEF(1,L)*NTFLX(LM2)+PCOEF(2,L)*NTFLX(LM1)             &
     &  +PCOEF(3,L)*NTFLX(L)+PCOEF(4,L)*NTFLX(LP1)
      CLRTSM(L)=CLRTSM(L)+CLRATE(L)
      UPFSM(L)=UPFSM(L)+UPFLX(L)
      DNFSM(L)=DNFSM(L)+DNFLX(L)
      UPFSSM(L)=UPFSSM(L)+UPFLXS(L)
      DNFSSM(L)=DNFSSM(L)+DNFLXS(L)
      NTFSM(L)=NTFSM(L)+NTFLX(L)

!     COOLING RATE WITH DIRECT SOLAR EXCLUDED:
      CLRTS0(L)=CLRTS0(L)+(PCOEF(1,L)*NTFLX0(LM2)+PCOEF(2,L)*NTFLX0(LM1)&
     &  +PCOEF(3,L)*NTFLX0(L)+PCOEF(4,L)*NTFLX0(LP1))

!     TOP OF ATMOSPHERE SPECTRAL COOLING RATE.
      CLRATE(ML)=PCOEF(1,ML)*NTFLX(LM1)                                 &
     &  +PCOEF(2,ML)*NTFLX(L)+PCOEF(3,ML)*NTFLX(LP1)
      CLRTSM(ML)=CLRTSM(ML)+CLRATE(ML)
      UPFSM(ML)=UPFSM(ML)+UPFLX(ML)
      DNFSM(ML)=DNFSM(ML)+DNFLX(ML)
      UPFSSM(ML)=UPFSSM(ML)+UPFLXS(ML)
      DNFSSM(ML)=DNFSSM(ML)+DNFLXS(ML)
      NTFSM(ML)=NTFSM(ML)+NTFLX(ML)

!     COOLING RATE WITH DIRECT SOLAR EXCLUDED:
      CLRTS0(ML)=CLRTS0(ML)+(PCOEF(1,ML)*NTFLX0(LM1)                    &
     &  +PCOEF(2,ML)*NTFLX0(L)+PCOEF(3,ML)*NTFLX0(LP1))

!     WRITE OUT SPECTRAL COOLING RATE DATA.
      IF(NPR.LE.-2)WRITE(ISPCCR,'(I8,1P,11E11.3:,/(8X,1P,11E11.3))')    &
     &  IV,(-CLRATE(L),L=1,ML)
      RETURN
      END