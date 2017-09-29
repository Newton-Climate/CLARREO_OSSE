      SUBROUTINE SNCMS(NTAU,S0,NLOS,PMOM,SSALB,DTAUC,UTAU,UU,S0CMS)

!     SNCMS COMPUTES AZIMUTHALLY-DEPENDENT LAYER SOLAR SOURCE FUNCTIONS.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       MAXCOE  MAXIMUM NUMBER OF PHASE FUNCTION LEGENDRE COEFFICIENTS.
!       MXUMU   MAXIMUM NUMBER OF LINE-OF-SIGHT NADIR ANGLES.
!       NTAU    NUMBER OF LAYER BOUNDARIES.
!       S0      TOP-OF-ATMOSPHERE SOLAR/LUNAR SPECTRAL
!               IRRADIANCE [W CM-2 / CM-1].
!       NLOS    NUMBER OF LINE-OF-SIGHT PATHS.
!       PMOM    LEGENDRE POLYNOMIAL MOMENTS; PMOM(1,*) CONTAINS
!               THE SCATTERING PHASE FUNCTION ASYMMETRY FACTORS.
!       SSALB   SINGLE SCATTERING ALBEDOS ORDERED FROM SPACE-TO-GROUND.
!       DTAUC   LAYER VERTICAL OPTICAL DEPTHS FROM SPACE-TO-GROUND.
!       UTAU    CUMULATIVE VERTICAL OPTICAL DEPTHS FROM SPACE-TO-GROUND.
!       UU      LAYER BOUNDARY TO SPACE (GROUND) SPECTRAL RADIANCES,
!               COMPUTED FOR A PLANE-PARALLEL ATMOSPHERE AND ORDERED
!               FROM SPACE-TO-GROUND [W CM-2 SR-1 / CM-1].
      INTEGER NTAU,NLOS
      REAL S0,PMOM(0:MAXCOE,*),SSALB(*),DTAUC(*),UTAU(*),               &
     &  UU(MXUMU,LAYDIM,*)

!     OUTPUT ARGUMENTS:
!       S0CMS   AZIMUTHALLY-DEPENDENT LAYER SOLAR SOURCE
!               FUNCTIONS [W CM-2 SR-1 / CM-1].
      REAL S0CMS(MLOS,*)

!     COMMONS:

!     /DISRT/
!       UMU      MONOTONICALLY INCREASING LIST OF DISTINCT USER-PATH
!                COSINE POLAR ANGLES.
!       PHI      MONOTONICALLY INCREASING LIST OF DISTINCT RELATIVE
!                SOLAR AZIMUTH ANGLES [0 TO 180 DEG].
!       NSTR     NUMBER OF DISCRETE ORDINATE STREAMS.
!       NAZ      NUMBER OF DISORT AZIMUTH COMPONENTS.
!       N2GAUS   ORDER OF DOUBLE-GAUSS QUADRATURES.
!       NUMU     NUMBER OF DISTINCT USER LINE-OF-SIGHT POLAR ANGLES.
!       MAPUMU   MAPPING FROM LINE-OF-SIGHT INDEX TO UMU ARRAY ENTRY.
!       NPHI     NUMBER OF DISTINCT RELATIVE SOLAR AZIMUTH ANGLES.
!       MAPPHI   MAPPING FROM LINE-OF-SIGHT INDEX TO PHI ARRAY ENTRY.
!       DIS      LOGICAL FLAG, TRUE FOR DISORT MULTIPLE SCATTERING.
!       DISAZM   LOGICAL FLAG, TRUE FOR DISORT WITH AZIMUTH DEPENDENCE.
!       DISALB   LOGICAL FLAG, TRUE FOR DISORT SPHERICAL ALBEDO OPTION.
!       LDISCL   LOGICAL FLAG, TRUE FOR ISAACS SCALED TO DISORT.
      REAL UMU,PHI
      INTEGER NSTR,NAZ,N2GAUS,NUMU,MAPUMU,NPHI,MAPPHI
      LOGICAL DIS,DISAZM,DISALB,LDISCL
      COMMON/DISRT/UMU(MXUMU),PHI(MAXPHI),NSTR,NAZ,N2GAUS,NUMU,         &
     &  MAPUMU(MLOS),NPHI,MAPPHI(MLOS),DIS,DISAZM,DISALB,LDISCL
      SAVE /DISRT/

!     /ANGSRF/
!       CVWSRF  COSINE OF THE VIEW ZENITH ANGLE FROM THE SURFACE.
!       CSNSRF  COSINE OF THE SOLAR (LUNAR) ZENITH AT SURFACE.
!       AZMSRF  RELATIVE AZIMUTH ANGLE (SUN - SENSOR AT SURFACE) [RAD].
!       UMU1    COSINE OF THE PATH NADIR ANGLE.
!               (AT H1ALT IF IMULT=1; AT OR "NEAR" H2ALT IF IMULT=-1)
!       UMU0    COSINE OF THE SOLAR ZENITH ANGLE.
!               (AT H1ALT IF IMULT=1; AT OR "NEAR" H2ALT IF IMULT=-1)
!       PHI1    RELATIVE AZIMUTH ANGLE (SUN - LOS PATH AT SENSOR) [DEG].
!               (AT H1ALT IF IMULT=1; AT OR "NEAR" H2ALT IF IMULT=-1)
!       CMU     COSINE OF THE NADIR ANGLES USED IN DISORT.
      REAL CVWSRF,CSNSRF,AZMSRF,UMU1,UMU0,PHI1,CMU
      COMMON/ANGSRF/CVWSRF(MLOS),CSNSRF(MLOS),AZMSRF(MLOS),UMU1(MLOS),  &
     &  UMU0,PHI1(MLOS),CMU(MI)

!     /MSRD/
!       CSSCAT   COSINE OF THE SCATTERING ANGLE.
!                (AT H1ALT IF IMULT=1; AT OR "NEAR" H2ALT IF IMULT=-1)
!       SLEGEN   Nth LEGENDRE POLYNOMIAL EVALUATED AT THE COSINE OF THE
!                SCATTERING ANGLE TIMES (2N+1)/4pi (N=0 TO NSTR-1).
!       CSZEN0   LAYER BOUNDARY COSINE OF SOLAR/LUNAR ZENITH.
!       CSZEN    LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH.
!       CSZENX   AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING
!                (AWAY FROM EARTH) THE CURRENT LAYER.
!       BBGRND   THERMAL EMISSION (FLUX) AT THE GROUND [W CM-2 / CM-1].
!       BBNDRY   LAYER BOUNDARY THERMAL EMISSION (FLUX) [W CM-2 / CM-1].
!       TCONT    LAYER CONTINUUM OPTICAL DEPTH.
!       TAUT     LAYER TOTAL VERTICAL EXTINCTION OPTICAL DEPTH.
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
!       CO_LIN   TRUE IF LOS AND SOLAR PATHS ARE NEARLY IDENTICAL.
!bug  REAL CSSCAT,SLEGEN,CSZEN0,CSZEN,CSZENX,TCONT,TAUT,GTSCAT,COSBAR,  &
!bug &  BBGRND,BBNDRY,S0DEP,S0TRN,DEPRAT,UPF,DNF,UPFS,DNFS
!bug  LOGICAL CO_LIN
!bug  COMMON/MSRD/CSSCAT(MLOS),SLEGEN(0:MAZ,MLOS),CSZEN0(LAYDIM),       &
!bug &  CSZEN(LAYDIM),CSZENX(LAYDIM),TCONT(LAYDIM),TAUT(MXKSUB,LAYDIM), &
!bug &  GTSCAT(0:MXCMU,1:LAYDIM),COSBAR(LAYDIM),BBGRND,BBNDRY(LAYDIM),  &
!bug &  S0DEP(MXKSUB,LAYTWO),S0TRN(MXKSUB,LAYTWO),DEPRAT(MXKSUB,LAYDIM),&
!bug &  UPF(MXKSUB,LAYDIM),DNF(MXKSUB,LAYDIM),UPFS(MXKSUB,LAYDIM),      &
!bug &  DNFS(MXKSUB,LAYDIM),CO_LIN(MLOS)

!     FUNCTIONS:
!       HENGNS  HENYEY-GREENSTEIN SCATTERING PHASE FUNCTION [SR-1].
!       PFEXP   LEGENDRE EXPANSION SCATTERING PHASE FUNCTION [SR-1].
!bug  REAL HENGNS,PFEXP
      REAL PFEXP

!     LOCAL VARIABLES:
!       ILOS     LOOP INDEX FOR LINE-FO-SIGHT PATHS.
!       N        LAYER INDEX INCREASING SPACE-TO-GROUND.
!       NP1      N PLUS 1.
!       IUMU     INDEX OF PATH NADIR COSINE.
!       IPHI     INDEX OF PATH RELATIVE AZIMUTH.
!       FACTOR   OPTICAL DEPTH SCALING FACTOR REQUIRED FOR LAYER
!                SINGLE SCATTER SOLAR RADIANCE CALCULATION.
!       OD       LAYER OPTICAL DEPTH ALONG THE LINE-OF-SIGHT.
!       ARG      PRODUCT OF FACTOR AND OD.
!       TRANSM   LAYER TRANSMITTANCE ALONG THE LINE-OF-SIGHT.
!       EMIS     LAYER EMISSIVITY ALONG THE LINE-OF-SIGHT.
!       DELISS   LAYER CONTRIBUTION TO THE SINGLE SCATTER
!                SOLAR/LUNAR SPECTRAL RADIANCE, ORDERED
!                FROM SPACE-TO-GROUND [W CM-2 SR-1 / CM-1].
!       DELI     LAYER CONTRIBUTION TO THE SPECTRAL RADIANCE, ORDERED
!                FROM SPACE-TO-GROUND [W CM-2 SR-1 / CM-1].
      INTEGER ILOS,N,NP1,IUMU,IPHI
      REAL FACTOR,OD,ARG,TRANSM,EMIS,DELISS,DELI
      
!*************VINCENT ROSS ADDED FOR SOURCE FUNCTION IMPROVEMENT***************
#ifdef DISORT_BOUND_SRC
!       ODCUM   CUMULATIVE OPTICAL DEPTH ALONG DISORT SLANT PATH
!       TAUCAL  OPTICAL DEPTH FROM SEGMENT START WHERE
!               SOURCE IS CALCULATED
!       ODCAL   CUMULATIVE PATH OPTICAL DEPTH FROM SPACE 
!               WHERE SOURCE IS CALCULATED
!       ODOLD   PREVIOUS LAYER (N-1) ODCAL
!       SNEW    NEW SOURCE FUNCTION (TEMPORARY SWAP VARIABLE)
!       SOLD    OLD SOURCE FUNCTION (TEMPORARY SWAP VARIABLE)
      DOUBLE PRECISION ODCUM,TAUCAL,ODCAL,ODOLD,SNEW,SOLD
      ODCUM=0
#endif
!*************END VINCENT ROSS ADDITION****************************************

!     LOOP OVER DISTINCT LINE-OF-SIGHT NADIR COSINES.
      DO ILOS=1,NLOS

!         BRANCH BASED ON LOOK DIRECTION:
          FACTOR=1+UMU1(ILOS)/UMU0
          N=1
          IUMU=MAPUMU(ILOS)
          IPHI=MAPPHI(ILOS)
          IF(UMU1(ILOS).GT.0.)THEN

!             LOOKING DOWN; LOOP OVER LAYERS:
              DO NP1=2,NTAU

!                 N IS THE LAYER BOUNDARY AT THE SENSOR; NP1 IS THE
!                 NEXT LAYER BOUNDARY ALONG THE LINE-OF-SIGHT.
                  OD=DTAUC(N)/UMU1(ILOS)
                  IF(OD.LT.1.E-12)THEN

!                     LAYER OPTICAL DEPTH IS SMALL; SET MULTIPLE
!                     SCATTERING SOLAR SOURCE TERM TO ZERO.
                      S0CMS(ILOS,N)=0.
                  ELSE

!                     LAYER TRANSMITTANCE AND EMISSIVITY:
                      IF(OD.LT..001)THEN
                          EMIS=OD*(1-OD/2)
                          TRANSM=1-EMIS
                      ELSE
                          TRANSM=EXP(-OD)
                          EMIS=1-TRANSM
                      ENDIF

!                     SINGLE SCATTER SOLAR RADIANCE:
                      ARG=FACTOR*OD
                      DELISS=-UTAU(N)/UMU0
                      IF(ARG.LT..001)THEN
                          DELISS=OD*(1-ARG/2)*EXP(DELISS)
                      ELSE
                          DELISS=(EXP(DELISS)-EXP(DELISS-ARG))/FACTOR
                      ENDIF
                      DELISS=S0*SSALB(N)*PFEXP(ILOS,PMOM(0,N))*DELISS
!bug  DELISS=S0*SSALB(N)*HENGNS(PMOM(1,N),CSSCAT(ILOS))*DELISS

!                     LAYER MULTIPLE SCATTERING SOLAR SOURCE:
                      DELI=UU(IUMU,N,IPHI)-UU(IUMU,NP1,IPHI)*TRANSM
!     write(*,'(4i3,a,1p,2e12.5)')ilos,iumu,n,iphi,
!    &  ' uu',UU(IUMU,N,IPHI),UU(IUMU,NP1,IPHI)*TRANSM
                      IF(DELI.LT.1.001*DELISS)THEN

!                         SCATTERED SOLAR IS DOMINATED BY SINGLE SCATTER
!                         SO JUST SET THE MULTIPLE SCATTER PART TO ZERO.
                          S0CMS(ILOS,N)=0.
                      ELSE

!                         DIVIDE BY LAYER EMISSIVITY TO OBTAIN SOURCE.
                          S0CMS(ILOS,N)=(DELI-DELISS)/EMIS
                      ENDIF
                  ENDIF
!*************VINCENT ROSS CHANGED FOR SOURCE FUNCTION IMPROVEMENT*****************
#ifdef DISORT_BOUND_SRC
              
                  SNEW = S0CMS(ILOS,N)

!                 APPROXIMATE OPTICAL 'CENTER OF MASS' WHERE SNEW WAS CALCULATED
                  TAUCAL = (-TRANSM*(OD+1)+1)/(-TRANSM+1)
              
                  ODCAL = ODCUM+TAUCAL
              
!                 ON SECOND PASS, EXTRAPOLATE TO FIRST BOUNDARY
                  IF(N.EQ.2)THEN
                      S0CMS(ILOS,N-1) = SOLD + 
     &                (SNEW-SOLD)/(ODCAL-ODOLD)*(-ODOLD)
                  ENDIF

!                 INTERPOLATE TO CURRENT BOUNDARY
                  IF(N.GT.1)THEN
                      S0CMS(ILOS,N) = SOLD + 
     &                (SNEW-SOLD)/(ODCAL-ODOLD)*(ODCUM-ODOLD)
                  ENDIF

                  ODCUM=ODCUM+OD
              
!                 ON LAST PASS, EXTRAPOLATE TO LAST BOUNDARY
                  IF(NP1.EQ.NTAU) THEN
                      S0CMS(ILOS,NP1) = SOLD + 
     &                (SNEW-SOLD)/(ODCAL-ODOLD)*(ODCUM-ODOLD)
                  ENDIF
              
                  ODOLD=ODCAL
                  SOLD=SNEW
#endif
!*************END VINCENT ROSS*****************************************************
                  N=NP1

!             END LAYER LOOP:
              ENDDO
          ELSEIF(UMU1(ILOS).LT.0.)THEN

!             LOOKING UP; LOOP OVER LAYERS:
              DO NP1=2,NTAU

!                 NP1 IS THE LAYER BOUNDARY AT THE SENSOR; N IS THE
!                 NEXT LAYER BOUNDARY ALONG THE LINE-OF-SIGHT.
                  OD=-DTAUC(N)/UMU1(ILOS)
                  IF(OD.LT.1.E-12)THEN

!                     LAYER OPTICAL DEPTH IS SMALL; SET MULTIPLE
!                     SCATTERING SOLAR SOURCE TERM TO ZERO.
                      S0CMS(ILOS,N)=0.
                  ELSE

!                     LAYER TRANSMITTANCE AND EMISSIVITY:
                      IF(OD.LT..001)THEN
                          EMIS=OD*(1-OD/2)
                          TRANSM=1-EMIS
                      ELSE
                          TRANSM=EXP(-OD)
                          EMIS=1-TRANSM
                      ENDIF

!                     SINGLE SCATTER SOLAR RADIANCE:
                      ARG=FACTOR*OD
                      DELISS=-UTAU(NP1)/UMU0
                      IF(ABS(ARG).LT..001)THEN
                          DELISS=OD*(1-ARG/2)*EXP(DELISS)
                      ELSE
                          DELISS=(EXP(DELISS)-EXP(DELISS-ARG))/FACTOR
                      ENDIF
                      DELISS=S0*SSALB(N)*PFEXP(ILOS,PMOM(0,N))*DELISS
!bug  DELISS=S0*SSALB(N)*HENGNS(PMOM(1,N),CSSCAT(ILOS))*DELISS

!                     MULTIPLE SCATTERING SOURCE:
                      DELI=UU(IUMU,NP1,IPHI)-UU(IUMU,N,IPHI)*TRANSM
                      IF(DELI.LT.1.001*DELISS)THEN

!                         SCATTERED SOLAR IS DOMINATED BY SINGLE SCATTER
!                         SO JUST SET THE MULTIPLE SCATTER PART TO ZERO.
                          S0CMS(ILOS,N)=0.
                      ELSE

!                         DIVIDE BY LAYER EMISSIVITY TO OBTAIN SOURCE.
                          S0CMS(ILOS,N)=(DELI-DELISS)/EMIS
                      ENDIF
                  ENDIF
!*************VINCENT ROSS CHANGED FOR SOURCE FUNCTION IMPROVEMENT*****************
#ifdef DISORT_BOUND_SRC
          
                  SNEW = S0CMS(ILOS,N)                                  
 
!                 APPROXIMATE OPTICAL 'CENTER OF MASS' WHERE SNEW WAS CALCULATED
                  TAUCAL = (-TRANSM*(OD+1)+1)/(-TRANSM+1)

                  ODCAL = ODCUM+OD-TAUCAL
          
!                 ON SECOND PASS, EXTRAPOLATE TO FIRST BOUNDARY
                  IF(N.EQ.2)THEN
                      S0CMS(ILOS,N-1) = SOLD + 
     &                (SNEW-SOLD)/(ODCAL-ODOLD)*(-ODOLD)
                  ENDIF
                  
!                 INTERPOLATE TO CURRENT BOUNDARY
                  IF(N.GT.1)THEN
                      S0CMS(ILOS,N) = SOLD + 
     &                (SNEW-SOLD)/(ODCAL-ODOLD)*(ODCUM-ODOLD)
                  ENDIF
                
                  ODCUM=ODCUM+OD
                  
!                 ON LAST PASS, EXTRAPOLATE TO LAST BOUNDARY
                  IF(NP1.EQ.NTAU) THEN
                      S0CMS(ILOS,NP1) = SOLD + 
     &                (SNEW-SOLD)/(ODCAL-ODOLD)*(ODCUM-ODOLD)
                  ENDIF

                  ODOLD=ODCAL
                  SOLD=SNEW
#endif
!*************END VINCENT ROSS*****************************************************
                  N=NP1

!             END LAYER LOOP:
              ENDDO
          ENDIF

!     END LINE-OF-SIGHT LOOP:
      ENDDO

!     RETURN TO MSRAD:
      RETURN
      END
