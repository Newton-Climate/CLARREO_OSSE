      SUBROUTINE WTCOOL(IBNDWD,IEMSCT)

!     ROUTINE TO WRITE THE IN-BAND COOLING RATES AND VERTICAL
!     FLUXES.  VALUES BETWEEN LAYER BOUNDARIES ARE DETERMINED
!     BY FITTING A CUBIC POLYNOMIAL IN LOG PRESSURE TO UPPER
!     AND LOWER LAYER BOUNDARY VALUES AND THEIR FIRST DERIVATIVES.

!     DECLARE ROUTINE ARGUMENTS:
!       IBNDWD   SPECTRAL BIN WIDTH [CM-1].
!       IEMSCT   MODTRAN RADIATION TRANSPORT MODE FLAG
!                  0 FOR TRANSMITTANCE ONLY
!                  1 FOR THERMAL RADIANCE
!                  2 FOR THERMAL + SOLAR RADIANCE
!                  3 FOR TRANSMITTED SOLAR IRRADIANCE
      INTEGER IBNDWD,IEMSCT

!     PARAMETERS:
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'
      INTEGER NGRID
      PARAMETER(NGRID=3)

!     COMMONS:

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

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     DECLARE LOCAL VARIABLES:
      INTEGER I,L,LP1,LP2
      REAL CONVRT,PLOG1,PLOG2,PLOG3,PLOG21,PLOG31,PLOG32,DCOEF1,        &
     &  DCOEF2,DCOEF3,DZ1,DZ2,DCLRT1,DCLRT2,D0CRT1,D0CRT2,DUPF1,        &
     &  DUPF2,DDNF1,DDNF2,DUPFS1,DUPFS2,DDNFS1,DDNFS2,DNTF1,            &
     &  DNTF2,DENOM,RATIO,RATM1,COEF1,COEF2,RATIO2,P1LG21,P2LG21,       &
     &  TSTDIF,TSTMIN,TSTMAX,TST,COEF(2,NGRID),DCOEF(2,NGRID)
!9    INTEGER ICASE
!9    DATA ICASE/0/

!     THE IN-BAND FLUXES NEED TO BE MULTIPLIED BY THE
!     BAND WIDTH AND CONVERTED FROM W/CM2 TO W/M2.
      CONVRT=10000.*IBNDWD

!     WRITE HEADER AND GROUND DATA
      IF(.NOT.LJMASS) THEN
        IF(IEMSCT.EQ.1)THEN
          WRITE(IPR,'(/2A,2(I6,A),/3(/2A),/5X,A1,2F12.6,5F14.6)')       &
     &      ' INTEGRATED COOLING RATES AND VERTICAL',                   &
     &      ' FLUXES FROM',IV1,' TO',IV2,' CM-1:',                      &
     &      '  CALC                                  COOL',             &
     &      'ING RATE          TOTAL FLUX          TOTAL THERMAL',      &
     &      ' LAYER    ALTITUDE    PRESSURE        TOTAL ',             &
     &      ' NO DIRECT SUN     UP-DN-DIR       FLUX UP     FLUX DOWN', &
     &      '            (KM)          (MB)       (K/DAY)',             &
     &      '       (K/DAY)        (W/M2)        (W/M2)        (W/M2)', &
     &      '1',ZM(1),PM(1),-IBNDWD*CLRTSM(1),-IBNDWD*CLRTS0(1),        &
     &      CONVRT*NTFSM(1),CONVRT*UPFSM(1),CONVRT*DNFSM(1)
        ELSEIF(DIS .AND. .NOT.DISAZM)THEN
          WRITE(IPR,'(/2A,2(I6,A),/3(/2A),/5X,A1,2F12.6,5F14.6)')       &
     &      ' INTEGRATED COOLING RATES AND VERTICAL',                   &
     &      ' FLUXES FROM',IV1,' TO',IV2,' CM-1:',                      &
     &      '  CALC                                  COOL',             &
     &      'ING RATE          TOTAL FLUX   THERMAL + SCATTERED SOLAR', &
     &      ' LAYER    ALTITUDE    PRESSURE        TOTAL ',             &
     &      ' NO DIRECT SUN     UP-DN-DIR       FLUX UP     FLUX DOWN', &
     &      '            (KM)          (MB)       (K/DAY)',             &
     &      '       (K/DAY)        (W/M2)        (W/M2)        (W/M2)', &
     &      '1',ZM(1),PM(1),-IBNDWD*CLRTSM(1),-IBNDWD*CLRTS0(1),        &
     &      CONVRT*NTFSM(1),CONVRT*UPFSM(1),CONVRT*DNFSM(1)
        ELSE
          WRITE(IPR,'(/2A,2(I6,A),/3(/3A),/5X,A1,2F12.6,7F14.6)')       &
     &      ' INTEGRATED COOLING RATES AND VERTICAL',                   &
     &      ' FLUXES FROM',IV1,' TO',IV2,' CM-1:',                      &
     &      '  CALC                        ',                           &
     &      '          COOLING RATE          TOTAL FLUX',               &
     &      '          TOTAL THERMAL              SCATTERED SOLAR',     &
     &      ' LAYER    ALTITUDE    PRESSURE',                           &
     &      '        TOTAL  NO DIRECT SUN     UP-DN-DIR',               &
     &      '       FLUX UP     FLUX DOWN       FLUX UP     FLUX DOWN', &
     &      '            (KM)          (MB)',                           &
     &      '       (K/DAY)       (K/DAY)        (W/M2)',               &
     &      '        (W/M2)        (W/M2)        (W/M2)        (W/M2)', &
     &      '1',ZM(1),PM(1),-IBNDWD*CLRTSM(1),-IBNDWD*CLRTS0(1),        &
     &      CONVRT*NTFSM(1),CONVRT*UPFSM(1),CONVRT*DNFSM(1),            &
     &      CONVRT*UPFSSM(1),CONVRT*DNFSSM(1)
        ENDIF
      ENDIF

!     DETERMINE THE DERIVATIVES AT THE SECOND LAYER BOUNDARY.
      PLOG1=LOG(PM(1))
      PLOG2=LOG(PM(2))
      PLOG3=LOG(PM(3))
      PLOG21=PLOG2-PLOG1
      PLOG31=PLOG3-PLOG1
      PLOG32=PLOG3-PLOG2
      DCOEF1=-PLOG32/(PLOG21*PLOG31)
      DCOEF2=1/PLOG21-1/PLOG32
      DCOEF3=PLOG21/(PLOG31*PLOG32)
      DZ2=DCOEF1*ZM(1)+DCOEF2*ZM(2)+DCOEF3*ZM(3)
      DCLRT2=DCOEF1*CLRTSM(1)+DCOEF2*CLRTSM(2)+DCOEF3*CLRTSM(3)
      D0CRT2=DCOEF1*CLRTS0(1)+DCOEF2*CLRTS0(2)+DCOEF3*CLRTS0(3)
      DUPF2=DCOEF1*UPFSM(1)+DCOEF2*UPFSM(2)+DCOEF3*UPFSM(3)
      DDNF2=DCOEF1*DNFSM(1)+DCOEF2*DNFSM(2)+DCOEF3*DNFSM(3)
      DUPFS2=DCOEF1*UPFSSM(1)+DCOEF2*UPFSSM(2)+DCOEF3*UPFSSM(3)
      DDNFS2=DCOEF1*DNFSSM(1)+DCOEF2*DNFSSM(2)+DCOEF3*DNFSSM(3)
      DNTF2=DCOEF1*NTFSM(1)+DCOEF2*NTFSM(2)+DCOEF3*NTFSM(3)

!     DETERMINE VALUES WITHIN FIRST LAYER (BETWEEN PLOG1 AND PLOG2)
      DENOM=NGRID+1
      IF(CLRTSM(1).LT.CLRTSM(2))THEN
          TSTDIF=.25*(CLRTSM(2)-CLRTSM(1))
          TSTMIN=CLRTSM(1)-TSTDIF
          TSTMAX=CLRTSM(2)+TSTDIF
      ELSE
          TSTDIF=.25*(CLRTSM(1)-CLRTSM(2))
          TSTMIN=CLRTSM(2)-TSTDIF
          TSTMAX=CLRTSM(1)+TSTDIF
      ENDIF
      DO 10 I=1,NGRID
          RATIO=I/DENOM
          RATM1=RATIO-1
          COEF1=RATM1**2
          COEF2=RATIO*(2-RATIO)
          DCOEF2=RATIO*RATM1*PLOG21
          TST=COEF1*CLRTSM(1)+COEF2*CLRTSM(2)+DCOEF2*DCLRT2
          IF(TST.GT.TSTMIN .AND. TST.LT.TSTMAX)THEN
              IF(IEMSCT.EQ.1 .OR. (DIS .AND. .NOT.DISAZM))THEN
                  IF(.NOT.LJMASS) WRITE(IPR,'(6X,2F12.6,5F14.6)')       &
     &              (COEF1*ZM(1)+COEF2*ZM(2)+DCOEF2*DZ2),               &
     &              EXP(PLOG1+RATIO*PLOG21),-IBNDWD*TST,                &
     &              -IBNDWD*(COEF1*CLRTS0(1)                            &
     &                      +COEF2*CLRTS0(2)+DCOEF2*D0CRT2),            &
     &               CONVRT*(COEF1*NTFSM(1)                             &
     &                      +COEF2*NTFSM(2) +DCOEF2*DNTF2 ),            &
     &               CONVRT*(COEF1*UPFSM(1)                             &
     &                      +COEF2*UPFSM(2) +DCOEF2*DUPF2 ),            &
     &               CONVRT*(COEF1*DNFSM(1)                             &
     &                      +COEF2*DNFSM(2) +DCOEF2*DDNF2 )
              ELSE
                  IF(.NOT.LJMASS) WRITE(IPR,'(6X,2F12.6,7F14.6)')       &
     &              (COEF1*ZM(1)+COEF2*ZM(2)+DCOEF2*DZ2),               &
     &              EXP(PLOG1+RATIO*PLOG21),-IBNDWD*TST,                &
     &              -IBNDWD*(COEF1*CLRTS0(1)                            &
     &                      +COEF2*CLRTS0(2)+DCOEF2*D0CRT2),            &
     &               CONVRT*(COEF1*NTFSM(1)                             &
     &                      +COEF2*NTFSM(2) +DCOEF2*DNTF2 ),            &
     &               CONVRT*(COEF1*UPFSM(1)                             &
     &                      +COEF2*UPFSM(2) +DCOEF2*DUPF2 ),            &
     &               CONVRT*(COEF1*DNFSM(1)                             &
     &                      +COEF2*DNFSM(2) +DCOEF2*DDNF2 ),            &
     &               CONVRT*(COEF1*UPFSSM(1)                            &
     &                      +COEF2*UPFSSM(2)+DCOEF2*DUPFS2),            &
     &               CONVRT*(COEF1*DNFSSM(1)                            &
     &                      +COEF2*DNFSSM(2)+DCOEF2*DDNFS2)
               ENDIF
           ENDIF

!         DETERMINE COEFFICIENTS FOR INTERMEDIATE LAYERS
          RATIO2=RATIO**2
          COEF(1,I)=COEF1*(1+2*RATIO)
          COEF(2,I)=RATIO2*(3-2*RATIO)
          DCOEF(1,I)=COEF1*RATIO
          DCOEF(2,I)=RATIO2*RATM1
   10 CONTINUE

!     WRITE VALUES FOR TOP OF THE FIRST LAYER.
      IF(.NOT.LJMASS) THEN
        IF(IEMSCT.EQ.1 .OR. (DIS .AND. .NOT.DISAZM))THEN
          WRITE(IPR,'(5X,A1,2F12.6,5F14.6)')                            &
     &      '2',ZM(2),PM(2),-IBNDWD*CLRTSM(2),-IBNDWD*CLRTS0(2),        &
     &      CONVRT*NTFSM(2),CONVRT*UPFSM(2),CONVRT*DNFSM(2)
        ELSE
          WRITE(IPR,'(5X,A1,2F12.6,7F14.6)')                            &
     &      '2',ZM(2),PM(2),-IBNDWD*CLRTSM(2),-IBNDWD*CLRTS0(2),        &
     &      CONVRT*NTFSM(2),CONVRT*UPFSM(2),CONVRT*DNFSM(2),            &
     &      CONVRT*UPFSSM(2),CONVRT*DNFSSM(2)
        ENDIF
      ENDIF

!     LOOP OVER INTERMEDIATE LAYERS.
      L=2
      LP1=3
      DO 30 LP2=4,ML

!         DETERMINE UPPER BOUNDARY DERIVATIVES.
          DZ1=DZ2
          DCLRT1=DCLRT2
          D0CRT1=D0CRT2
          DUPF1=DUPF2
          DDNF1=DDNF2
          DUPFS1=DUPFS2
          DDNFS1=DDNFS2
          DNTF1=DNTF2
          PLOG1=PLOG2
          PLOG2=PLOG3
          PLOG3=LOG(PM(LP2))
          PLOG21=PLOG2-PLOG1
          PLOG31=PLOG3-PLOG1
          PLOG32=PLOG3-PLOG2
          DCOEF1=-PLOG32/(PLOG21*PLOG31)
          DCOEF2=1/PLOG21-1/PLOG32
          DCOEF3=PLOG21/(PLOG31*PLOG32)
          DZ2=DCOEF1*ZM(L)+DCOEF2*ZM(LP1)+DCOEF3*ZM(LP2)
          DCLRT2=DCOEF1*CLRTSM(L)+DCOEF2*CLRTSM(LP1)+DCOEF3*CLRTSM(LP2)
          D0CRT2=DCOEF1*CLRTS0(L)+DCOEF2*CLRTS0(LP1)+DCOEF3*CLRTS0(LP2)
          DUPF2=DCOEF1*UPFSM(L)+DCOEF2*UPFSM(LP1)+DCOEF3*UPFSM(LP2)
          DDNF2=DCOEF1*DNFSM(L)+DCOEF2*DNFSM(LP1)+DCOEF3*DNFSM(LP2)
          DUPFS2=DCOEF1*UPFSSM(L)+DCOEF2*UPFSSM(LP1)+DCOEF3*UPFSSM(LP2)
          DDNFS2=DCOEF1*DNFSSM(L)+DCOEF2*DNFSSM(LP1)+DCOEF3*DNFSSM(LP2)
          DNTF2=DCOEF1*NTFSM(L)+DCOEF2*NTFSM(LP1)+DCOEF3*NTFSM(LP2)
          IF(CLRTSM(L).LT.CLRTSM(LP1))THEN
              TSTDIF=.25*(CLRTSM(LP1)-CLRTSM(L  ))
              TSTMIN=CLRTSM(L  )-TSTDIF
              TSTMAX=CLRTSM(LP1)+TSTDIF
          ELSE
              TSTDIF=.25*(CLRTSM(L  )-CLRTSM(LP1))
              TSTMIN=CLRTSM(LP1)-TSTDIF
              TSTMAX=CLRTSM(L  )+TSTDIF
          ENDIF

!         WRITE VALUES WITHIN LAYER.
          DO 20 I=1,NGRID
              P1LG21=PLOG21*DCOEF(1,I)
              P2LG21=PLOG21*DCOEF(2,I)
              TST=COEF(1,I)*CLRTSM(L  )+P1LG21*DCLRT1+                  &
     &            COEF(2,I)*CLRTSM(LP1)+P2LG21*DCLRT2
              IF(TST.GT.TSTMIN .AND. TST.LT.TSTMAX)THEN
                IF(.NOT.LJMASS) THEN
                  IF(IEMSCT.EQ.1 .OR. (DIS .AND. .NOT.DISAZM))THEN
                      WRITE(IPR,'((6X,2F12.6,5F14.6))')                 &
     &                          (COEF(1,I)*ZM(L  )    +P1LG21*DZ1       &
     &                          +COEF(2,I)*ZM(LP1)    +P2LG21*DZ2),     &
     &                  EXP(PLOG1+I*PLOG21/DENOM),-IBNDWD*TST,          &
     &                  -IBNDWD*(COEF(1,I)*CLRTS0(L  )+P1LG21*D0CRT1    &
     &                          +COEF(2,I)*CLRTS0(LP1)+P2LG21*D0CRT2),  &
     &                   CONVRT*(COEF(1,I)* NTFSM(L  )+P1LG21*DNTF1     &
     &                          +COEF(2,I)* NTFSM(LP1)+P2LG21*DNTF2 ),  &
     &                   CONVRT*(COEF(1,I)* UPFSM(L  )+P1LG21*DUPF1     &
     &                          +COEF(2,I)* UPFSM(LP1)+P2LG21*DUPF2 ),  &
     &                   CONVRT*(COEF(1,I)* DNFSM(L  )+P1LG21*DDNF1     &
     &                          +COEF(2,I)* DNFSM(LP1)+P2LG21*DDNF2 )
                  ELSE
                      WRITE(IPR,'((6X,2F12.6,7F14.6))')                 &
     &                          (COEF(1,I)*ZM(L  )    +P1LG21*DZ1       &
     &                          +COEF(2,I)*ZM(LP1)    +P2LG21*DZ2),     &
     &                  EXP(PLOG1+I*PLOG21/DENOM),-IBNDWD*TST,          &
     &                  -IBNDWD*(COEF(1,I)*CLRTS0(L  )+P1LG21*D0CRT1    &
     &                          +COEF(2,I)*CLRTS0(LP1)+P2LG21*D0CRT2),  &
     &                   CONVRT*(COEF(1,I)* NTFSM(L  )+P1LG21*DNTF1     &
     &                          +COEF(2,I)* NTFSM(LP1)+P2LG21*DNTF2 ),  &
     &                   CONVRT*(COEF(1,I)* UPFSM(L  )+P1LG21*DUPF1     &
     &                          +COEF(2,I)* UPFSM(LP1)+P2LG21*DUPF2 ),  &
     &                   CONVRT*(COEF(1,I)* DNFSM(L  )+P1LG21*DDNF1     &
     &                          +COEF(2,I)* DNFSM(LP1)+P2LG21*DDNF2 ),  &
     &                   CONVRT*(COEF(1,I)*UPFSSM(L  )+P1LG21*DUPFS1    &
     &                          +COEF(2,I)*UPFSSM(LP1)+P2LG21*DUPFS2),  &
     &                   CONVRT*(COEF(1,I)*DNFSSM(L  )+P1LG21*DDNFS1    &
     &                          +COEF(2,I)*DNFSSM(LP1)+P2LG21*DDNFS2)
                  ENDIF
                ENDIF
              ENDIF
   20     CONTINUE

!         WRITE VALUES FOR TOP OF THE CURRENT LAYER.
          L=LP1
          LP1=LP2
          IF(.NOT.LJMASS) THEN
            IF(IEMSCT.EQ.1 .OR. (DIS .AND. .NOT.DISAZM))THEN
              WRITE(IPR,'(I6,2F12.6,7F14.6)')                           &
     &          L,ZM(L),PM(L),-IBNDWD*CLRTSM(L),-IBNDWD*CLRTS0(L),      &
     &          CONVRT*NTFSM(L),CONVRT*UPFSM(L),CONVRT*DNFSM(L)
            ELSE
              WRITE(IPR,'(I6,2F12.6,7F14.6)')                           &
     &          L,ZM(L),PM(L),-IBNDWD*CLRTSM(L),-IBNDWD*CLRTS0(L),      &
     &          CONVRT*NTFSM(L),CONVRT*UPFSM(L),CONVRT*DNFSM(L),        &
     &          CONVRT*UPFSSM(L),CONVRT*DNFSSM(L)
            ENDIF
          ENDIF
   30 CONTINUE

!     DETERMINE VALUES WITHIN LAST LAYER.
      IF(CLRTSM(L).LT.CLRTSM(ML))THEN
          TSTDIF=.25*(CLRTSM(ML)-CLRTSM(L ))
          TSTMIN=CLRTSM(L )-TSTDIF
          TSTMAX=CLRTSM(ML)+TSTDIF
      ELSE
          TSTDIF=.25*(CLRTSM(L )-CLRTSM(ML))
          TSTMIN=CLRTSM(ML)-TSTDIF
          TSTMAX=CLRTSM(L )+TSTDIF
      ENDIF
      DO 40 I=1,NGRID
          RATIO=I/DENOM
          COEF1=(1+RATIO)*(1-RATIO)
          DCOEF1=RATIO*(1-RATIO)*PLOG32
          COEF2=RATIO**2
          TST=COEF1*CLRTSM(L)+DCOEF1*DCLRT2+COEF2*CLRTSM(ML)
          IF(TST.GT.TSTMIN .AND. TST.LT.TSTMAX)THEN
            IF(.NOT.LJMASS) THEN
              IF(IEMSCT.EQ.1 .OR. (DIS .AND. .NOT.DISAZM))THEN
                WRITE(IPR,'(6X,2F12.6,5F14.6)')                         &
     &            (COEF1*ZM(L)+DCOEF1*DZ2+COEF2*ZM(ML)),                &
     &            EXP(PLOG2+RATIO*PLOG32),-IBNDWD*TST,                  &
     &            -IBNDWD*(COEF1*CLRTS0( L)+DCOEF1*D0CRT2               &
     &                    +COEF2*CLRTS0(ML)),                           &
     &             CONVRT*(COEF1*NTFSM( L) +DCOEF1*DNTF2                &
     &                    +COEF2*NTFSM(ML) ),                           &
     &             CONVRT*(COEF1*UPFSM( L) +DCOEF1*DUPF2                &
     &                    +COEF2*UPFSM(ML) ),                           &
     &             CONVRT*(COEF1*DNFSM( L) +DCOEF1*DDNF2                &
     &                    +COEF2*DNFSM(ML) )
              ELSE
                WRITE(IPR,'(6X,2F12.6,7F14.6)')                         &
     &            (COEF1*ZM(L)+DCOEF1*DZ2+COEF2*ZM(ML)),                &
     &            EXP(PLOG2+RATIO*PLOG32),-IBNDWD*TST,                  &
     &            -IBNDWD*(COEF1*CLRTS0( L)+DCOEF1*D0CRT2               &
     &                    +COEF2*CLRTS0(ML)),                           &
     &             CONVRT*(COEF1*NTFSM( L) +DCOEF1*DNTF2                &
     &                    +COEF2*NTFSM(ML) ),                           &
     &             CONVRT*(COEF1*UPFSM( L) +DCOEF1*DUPF2                &
     &                    +COEF2*UPFSM(ML) ),                           &
     &             CONVRT*(COEF1*DNFSM( L) +DCOEF1*DDNF2                &
     &                    +COEF2*DNFSM(ML) ),                           &
     &             CONVRT*(COEF1*UPFSSM( L)+DCOEF1*DUPFS2               &
     &                    +COEF2*UPFSSM(ML)),                           &
     &             CONVRT*(COEF1*DNFSSM( L)+DCOEF1*DDNFS2               &
     &                    +COEF2*DNFSSM(ML))
              ENDIF
            ENDIF
          ENDIF
   40 CONTINUE

!     WRITE VALUES FOR TOP OF ATMOSPHERE.
      IF(.NOT.LJMASS) THEN
        IF(IEMSCT.EQ.1 .OR. (DIS .AND. .NOT.DISAZM))THEN
          WRITE(IPR,'(I6,2F12.6,5F14.6)')                               &
     &      ML,ZM(ML),PM(ML),-IBNDWD*CLRTSM(ML),-IBNDWD*CLRTS0(ML),     &
     &      CONVRT*NTFSM(ML),CONVRT*UPFSM(ML),CONVRT*DNFSM(ML)
        ELSE
          WRITE(IPR,'(I6,2F12.6,7F14.6)')                               &
     &      ML,ZM(ML),PM(ML),-IBNDWD*CLRTSM(ML),-IBNDWD*CLRTS0(ML),     &
     &      CONVRT*NTFSM(ML),CONVRT*UPFSM(ML),CONVRT*DNFSM(ML),         &
     &      CONVRT*UPFSSM(ML),CONVRT*DNFSSM(ML)
        ENDIF
      ENDIF
      RETURN
      END
