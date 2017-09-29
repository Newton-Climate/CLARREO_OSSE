      SUBROUTINE BMOD(IK,IKMX,IPATH,IV,MSOFF)

!     THIS ROUTINE RETURNS THE TRANSMITTANCE AT A SPECTRAL
!     RESOLUTION OF 1 CM-1 FOR THE CFC'S AND THE "NMOL" SPECIES.
!           K = ( 1, 13)      H2O LINE (CENTERS, TAILS)
!             = ( 2, 14)      CO2 LINE (CENTERS, TAILS)
!             = ( 3, 15)       O3 LINE (CENTERS, TAILS)
!             = ( 4, 16)      N2O LINE (CENTERS, TAILS)
!             = ( 5, 17)       CO LINE (CENTERS, TAILS)
!             = ( 6, 18)      CH4 LINE (CENTERS, TAILS)
!             = ( 7, 19)       O2 LINE (CENTERS, TAILS)
!             = ( 8, 20)       NO LINE (CENTERS, TAILS)
!             = ( 9, 21)      SO2 LINE (CENTERS, TAILS)
!             = (10, 22)      NO2 LINE (CENTERS, TAILS)
!             = (11, 23)      NH3 LINE (CENTERS, TAILS)
!             = (12, 24)     HNO3 LINE (CENTERS, TAILS)

!     DECLARE INPUTS
      INTEGER IK,IKMX,IPATH,IV,MSOFF

!     INCLUDE PARAMETERS
      INCLUDE 'PARAMS.h'

!     LIST COMMONS
      INTEGER IBINX,IMOLX,IALFX
      REAL SDZX,ODZX
      COMMON/BMDCMX/SDZX(MXTEMP),ODZX(MXTEMP),IBINX,IMOLX,IALFX
      INCLUDE 'BASE.h'
      INCLUDE 'SOLS.h'

!     /CARD1/
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,MODTRN

!     /JM5/
!       IRPT     REPEAT INPUT FLAG (0=NONE, 1=ALL, 3=GEOM, 4=SPEC).
!       IFAC     CURRENT COLUMN SCALING FACTOR INDEX.
!       NFACMN   NUMBER OF COLUMN SCALING FACTOR LESS THAN 1.
!       NFACMX   NUMBER OF COLUMN SCALING FACTOR GREATER THAN 1.
!       FACMC    CURRENT COLUMN SCALING FACTOR.
!       SCALMN   MINIMUM COLUMN SCALING FACTOR.
!       SCALMX   MAXIMUM COLUMN SCALING FACTOR.
      INTEGER IRPT,IFAC,NFACMN,NFACMX
      REAL FACMC
      DOUBLE PRECISION SCALMN,SCALMX
      COMMON/JM5/SCALMN,SCALMX,IRPT,IFAC,NFACMN,NFACMX,FACMC
      INCLUDE 'BMHEAD.h'
      INCLUDE 'BMDAT.h'
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     DECLARE LOCAL FUNCTIONS
      REAL BMTRAN,BMTRN

!     DECLARE LOCAL VARIABLES
      INTEGER K,IT,I,LMAX,LS,KPNT,L,MSOFFL,J,JM1,ITLSUB,MDIVID
      REAL FREQ,V,F,STORE,ABSM,DINV,TAIL,WPTH,DEPTH,ODBAR,              &
     &  ADBAR,ACBAR,ACBAR2,TRANSM,TDEPTH(MTLSUB)

!     SAVED LOCAL VARIABLES AND DATA:
!     DOPFAC (UNITLESS) EQUALS SQRT(2 LN2 R T / M)/C WHERE T IS THE
!     STANDARD TEMPERATURE (273.15K) AND M IS MOLECULAR WEIGHT.
      INTEGER IPATH0
      LOGICAL MOLTRN
      REAL COLO3,SDSUM(NMOLXT),ODSUM(NMOLXT),DOPSUM(NMOLXT),            &
     &  COLSUM(NMOLXT),TAILSM(NMOLXT,MTLSUB),DOPFAC(NMOLXT)
      SAVE MOLTRN,COLO3,SDSUM,ODSUM,DOPSUM,COLSUM,TAILSM,DOPFAC
      DATA IPATH0/0/,(DOPFAC(K),K=1,NMOLXT)/                            &
     &  1.3945E-6,  0.8922E-6,  0.8543E-6,  0.8921E-6,  1.1183E-6,      &
     &  1.4777E-6,  1.0463E-6,  1.0805E-6,  0.7395E-6,  0.8726E-6,      &
     &  1.4342E-6,  0.7456E-6,  5.04981E-7, 5.38245E-7, 5.79088E-7,     &
     &  6.30907E-7, 6.36486E-7, 4.32375E-7, 4.52709E-7, 4.76212E-7,     &
     &  5.99528E-7, 6.65841E-7, 5.83391E-7, 4.7721E-7,  5.69488E-7/

!     IF IV IS BEYOND BAND MODEL TAPE RANGE, SET TRANSMITTANCES TO 1
      IF(IV.GT.MXFREQ)THEN
          DO 10 K=1,NMOL
              TX(KPOINT(K))=1.
   10     CONTINUE
          RETURN
      ENDIF

!     IK=0 IS THE INITIAL CALL FOR EACH WAVENUMBER AND IS MADE PRIOR
!     TO THE LOOP OVER LAYERS.  IK=-1 IS THE INITIAL CALL TO BMOD
!     AFTER THE MULTIPLE SCATTERING LOOP IS COMPLETE
      IF(IK.EQ.0)THEN

!         IS THERE MOLECULAR DATA FOR FREQUENCY IV?
          IF(IBIN.NE.IV)THEN
              MOLTRN=.FALSE.
              RETURN
          ENDIF
          MOLTRN=.TRUE.

!         FOR EACH SPECIES, SET THE 'NO DATA' INDICATOR
          DO 30 K=1,NMOLXT
              SD(1,K,0)=-99.9
              SD(1,K,1)=-99.9
   30     CONTINUE

!         LOAD MOLECULAR DATA FOR FREQUENCY IV
   40     CONTINUE
          CALL BMLOAD
          IF(IP.LT.LSTREC)THEN
              IP=IP+1
              READ(ITB,REC=IP)                                          &
     &          IBIN,IMOL,(SDZ(I),I=1,MXTEMP),IALF,(ODZ(I),I=1,MXTEMP)
              IF(IBIN.EQ.IV)GOTO40
          ENDIF

!         LOAD CFC DATA FOR FREQUENCY IV
   50     CONTINUE
          IF(IBINX.EQ.IV)THEN
              CALL BMXLD
              READ(ITBX,*,END=60)FREQ,IMOLX,(SDZX(IT),IT=1,NTEMP)
              IBINX=INT(FREQ+0.01)
              READ(ITBX,*)IALFX,(ODZX(IT),IT=1,NTEMP)
              GOTO50
          ENDIF
   60     CONTINUE

!         ZERO LAYER LOOP QUANTITIES; DEFINE STANDARD DOPPLER WIDTH
          V=FLOAT(IV)
          IF(IV.EQ.0)V=.25*IBNDWD
          DO 80 K=1,NMOLXT
              DOP0(K)=V*DOPFAC(K)
              SDSUM(K)=0.
              ODSUM(K)=0.
              COLSUM(K)=0.
              DOPSUM(K)=0.
              DO 70 ITLSUB=1,NTLSUB
                  TAILSM(K,ITLSUB)=0.
   70         CONTINUE
   80     CONTINUE
          COLO3=0.
          RETURN
      ELSEIF(IK.EQ.-1)THEN
          DO 100 K=1,NMOLXT
              SDSUM(K)=0.
              ODSUM(K)=0.
              COLSUM(K)=0.
              DOPSUM(K)=0.
              DO 90 ITLSUB=1,NTLSUB
                  TAILSM(K,ITLSUB)=0.
   90         CONTINUE
  100     CONTINUE
          COLO3=0.
          RETURN
      ENDIF
      IF(.NOT.MOLTRN)THEN
          DO 110 K=1,NMOL
              TX(KPOINT(K))=1.
  110     CONTINUE
          RETURN
      ENDIF

!     START CALCULATION OF MOLECULAR TRANSMITTANCE
      IF(IEMSCT.EQ.0 .OR. IEMSCT.EQ.3)THEN

!         LOOP OVER ALL LAYERS FOR TRANSMITTANCE CALCULATIONS
          LMAX=IKMX
      ELSEIF(IEMSCT.EQ.1)THEN

!         LOOP OVER SINGLE LAYER FOR RADIANCE CALCULATION WITHOUT SOLAR
          LMAX=IK
      ELSE
          GOTO(120,130,140),IPATH

!         SKIP LAYER LOOP COMPLETELY FOR FIRST SOLAR PATH WITH IPATH=1
  120     LMAX=0
          LS=1
          GOTO150

!         LOOP OVER SINGLE LAYER ONLY FOR "L PATH" WITH IPATH=2
  130     LMAX=IK
          IF(MSOFF.GT.0)LMAX=0
          LS=IK+1
          GOTO150

!         LOOP OVER SINGLE LAYER IF NOT PERFORMED WHEN IPATH EQUALED 2
  140     LMAX=0
          IF(MSOFF.GT.0 .OR. IPATH0.NE.2)LMAX=IK
  150     IPATH0=IPATH
      ENDIF

!     START SPECIES LOOP
      MDIVID=NTLSUB
      DO 280 K=1,NMOLXT
          KPNT=KPOINT(K)

!         FOR CROSS-SECTION SPECIES (K>NMOL) MDIVID IS 1.
          IF(K.EQ.NMOL+1)MDIVID=1

!         CHECK IF LINE CENTER CONTRIBUTES TO ABSORPTION
          IF(SD(1,K,0).GT.-99.)THEN

!             LOOP OVER LAYERS
              DO 170 L=IK,LMAX

!                 PATH AMOUNT
                  MSOFFL=MSOFF+L
                  WPTH=FACMC*WPATH(MSOFFL,KPNT)
                  IF(WPTH.GT.0.)THEN

!                     INTERPOLATE BAND MODEL PARAMETERS OVER TEMPERATURE
                      J=JJ(MSOFFL)
                      F=FF(MSOFFL)
                      JM1=J-1
                      ABSM=SD(J,K,0)
                      ABSM=ABSM+F*(SD(JM1,K,0)-ABSM)
                      DINV=OD(J,K)
                      DINV=DINV+F*(OD(JM1,K)-DINV)

!                     PERFORM CURTIS-GODSON SUMS
                      STORE=ABSM*WPTH
                      SDSUM(K)=SDSUM(K)+STORE
                      STORE=DINV*STORE
                      ODSUM(K)=ODSUM(K)+STORE
                      DOPSUM(K)=DOPSUM(K)+STORE*T5(MSOFFL)
                      STORE=STORE*PTM75(MSOFFL)
                      COLSUM(K)=COLSUM(K)+STORE
                      IF(K.EQ.3)COLO3=COLO3+STORE*DINV*PTM75(MSOFFL)
                      WPTH=WPTH*PATM(MSOFFL)
                      DO 160 ITLSUB=1,MDIVID
                          TAIL=SD(J,K,ITLSUB)
                          TAIL=TAIL+F*(SD(JM1,K,ITLSUB)-TAIL)
                          TAILSM(K,ITLSUB)=TAILSM(K,ITLSUB)+WPTH*TAIL
  160                 CONTINUE
                  ENDIF
  170         CONTINUE
              DEPTH=SDSUM(K)
              ODBAR=ODSUM(K)
              ADBAR=DOPSUM(K)
              ACBAR=COLSUM(K)
              ACBAR2=0.
              IF(K.EQ.3)ACBAR2=COLO3
              DO 180 ITLSUB=1,NTLSUB
                  TDEPTH(ITLSUB)=TAILSM(K,ITLSUB)
  180         CONTINUE

!             IF SOLAR PATH, CALCULATE ADDITIONAL LAYER
              IF(IEMSCT.EQ.2 .AND. IPATH.NE.3)THEN
                  MSOFFL=MSOFF+LS
                  WPTH=FACMC*WPATHS(MSOFFL,KPNT)
                  IF(WPTH.GT.0.)THEN
                      J=JJS(MSOFFL,K)
                      F=FFS(MSOFFL,K)
                      JM1=J-1
                      ABSM=SD(J,K,0)
                      ABSM=ABSM+F*(SD(JM1,K,0)-ABSM)
                      DINV=OD(J,K)
                      DINV=DINV+F*(OD(JM1,K)-DINV)
                      IF(MSOFF.EQ.0)THEN

!                         L-SHAPED PATH
                          STORE=ABSM*WPTH
                          DEPTH=DEPTH+STORE
                          STORE=DINV*STORE
                          ODBAR=ODBAR+STORE
                          ADBAR=ADBAR+STORE*T5S(MSOFFL,K)
                          STORE=STORE*PTM75S(MSOFFL,K)
                          ACBAR=ACBAR+STORE
                          IF(K.EQ.3)                                    &
     &                      ACBAR2=ACBAR2+STORE*DINV*PTM75S(MSOFFL,K)
                          WPTH=WPTH*PATMS(MSOFFL,K)
                          DO 190 ITLSUB=1,MDIVID
                              TAIL=SD(J,K,ITLSUB)
                              TAIL=TAIL+F*(SD(JM1,K,ITLSUB)-TAIL)
                              TDEPTH(ITLSUB)=TDEPTH(ITLSUB)+WPTH*TAIL
  190                     CONTINUE
                      ELSE

!                         SOLAR PATH ONLY
                          DEPTH=ABSM*WPTH
                          ODBAR=DINV*DEPTH
                          ADBAR=ODBAR*T5S(MSOFFL,K)
                          ACBAR=ODBAR*PTM75S(MSOFFL,K)
                          IF(K.EQ.3)ACBAR2=ACBAR*DINV*PTM75S(MSOFFL,K)
                          WPTH=WPTH*PATMS(MSOFFL,K)
                          DO 200 ITLSUB=1,MDIVID
                              TAIL=SD(J,K,ITLSUB)
                              TAIL=TAIL+F*(SD(JM1,K,ITLSUB)-TAIL)
                              TDEPTH(ITLSUB)=WPTH*TAIL
  200                     CONTINUE
                      ENDIF
                  ELSEIF(MSOFF.GT.0)THEN
                      TX(KPNT)=1.
                      GOTO280
                  ENDIF
              ENDIF

!             CHECK FOR WEAK LINE
              IF(DEPTH.LT.0.001)THEN
                  TRANSM=1.-DEPTH
              ELSE

!                 CALCULATE EQUIVALENT WIDTH TRANSMITTANCE
                  ODBAR=ODBAR/DEPTH
                  ADBAR=DOP0(K)*ADBAR/DEPTH
                  ACBAR=ALF0(K)*ACBAR/DEPTH
                  IF(ACBAR2.NE.0.)THEN
                      ACBAR2=ACBAR2/DEPTH*ALF0(K)**2
                      TRANSM=BMTRAN(DEPTH,ODBAR,ADBAR,ACBAR,ACBAR2)
                  ELSE
                      TRANSM=BMTRN(DEPTH,ODBAR,ADBAR,ACBAR)
                  ENDIF
              ENDIF

!             ADD LINE TAIL CONTRIBUTIONS.
              IF(TRANSM.GT.0.)THEN
                  TAIL=EXP(-TDEPTH(1))
                  DO 210 ITLSUB=2,MDIVID
                      TAIL=TAIL+EXP(-TDEPTH(ITLSUB))
  210             CONTINUE
                  TRANSM=TRANSM*TAIL/MDIVID
              ENDIF
              TX(KPNT)=TRANSM

!         CHECK IF LINE TAILS CONTRIBUTE TO ABSORPTION:
          ELSEIF(SD(1,K,1).GT.-99.)THEN

!             LOOP OVER LAYERS
              DO 230 L=IK,LMAX
                  MSOFFL=MSOFF+L

!                 PERFORM CURTIS-GODSON SUMS
                  WPTH=FACMC*WPATH(MSOFFL,KPNT)
                  IF(WPTH.GT.0.)THEN

!                     LINE TAILS ARE SCALED BY PRESSURE;
!                     EXTINCTION CROSS-SECTIONS ARE NOT.
                      IF(K.LE.NMOL)WPTH=WPTH*PATM(MSOFFL)

!                     INTERPOLATE BAND MODEL PARAMETERS OVER TEMPERATURE
                      J=JJ(MSOFFL)
                      F=FF(MSOFFL)
                      DO 220 ITLSUB=1,MDIVID
                          TAIL=SD(J,K,ITLSUB)
                          TAIL=TAIL+F*(SD(J-1,K,ITLSUB)-TAIL)
                          TAILSM(K,ITLSUB)=TAILSM(K,ITLSUB)+WPTH*TAIL
  220                 CONTINUE
                  ENDIF
  230         CONTINUE
              DO 240 ITLSUB=1,MDIVID
                  TDEPTH(ITLSUB)=TAILSM(K,ITLSUB)
  240         CONTINUE

!             IF SOLAR PATH, CALCULATE ADDITIONAL LAYER.
              IF(IEMSCT.EQ.2 .AND. IPATH.NE.3)THEN
                  MSOFFL=MSOFF+LS
                  WPTH=FACMC*WPATHS(MSOFFL,KPNT)
                  IF(WPTH.GT.0.)THEN

!                     LINE TAILS ARE SCALED BY PRESSURE;
!                     EXTINCTION CROSS-SECTIONS ARE NOT.
                      IF(K.LE.NMOL)WPTH=WPTH*PATMS(MSOFFL,K)
                      J=JJS(MSOFFL,K)
                      F=FFS(MSOFFL,K)
                      IF(MSOFF.EQ.0)THEN

!                         L-SHAPED PATH
                          DO 250 ITLSUB=1,MDIVID
                              TAIL=SD(J,K,ITLSUB)
                              TAIL=TAIL+F*(SD(J-1,K,ITLSUB)-TAIL)
                              TDEPTH(ITLSUB)=TDEPTH(ITLSUB)+WPTH*TAIL
  250                     CONTINUE
                      ELSE

!                         SOLAR PATH ONLY
                          DO 260 ITLSUB=1,MDIVID
                              TAIL=SD(J,K,ITLSUB)
                              TAIL=TAIL+F*(SD(J-1,K,ITLSUB)-TAIL)
                              TDEPTH(ITLSUB)=WPTH*TAIL
  260                     CONTINUE
                      ENDIF
                  ELSEIF(MSOFF.GT.0)THEN
                      TX(KPNT)=1.
                      GOTO280
                  ENDIF
              ENDIF

!             CALCULATE LINE TAIL TRANSMITTANCE
              TAIL=EXP(-TDEPTH(1))
              DO 270 ITLSUB=2,MDIVID
                  TAIL=TAIL+EXP(-TDEPTH(ITLSUB))
  270         CONTINUE
              TX(KPNT)=TAIL/MDIVID
              IF(TX(KPNT).GT.1.)TX(KPNT)=1.
          ELSE
              TX(KPNT)=1.
          ENDIF
  280 CONTINUE
      RETURN
      END
