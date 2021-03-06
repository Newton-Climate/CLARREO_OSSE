      SUBROUTINE SNCMS(MAXCMU,MAXUMU,MAXPHI,NLYR,NTAU,S0,               &
     &  UMU0,CSSCAT,UMU1,PMOM,SSALB,DTAUC,UTAU,UU,S0CMS)

!     SNCMS COMPUTES AZIMUTHALLY-DEPENDENT LAYER SOLAR SOURCE FUNCTIONS.

!     INPUT ARGUMENTS:
!       MAXCMU  MAXIMUM NUMBER OF DISCRETE ORDINATE STREAMS.
!       MAXUMU  MAXIMUM NUMBER OF LINE-OF-SIGHT NADIR ANGLES.
!       MAXPHI  MAXIMUM NUMBER OF LINE-OF-SIGHT AZIMUTH ANGLES.
!       NLYR    NUMBER OF LAYERS.
!       NTAU    NUMBER OF LAYER BOUNDARIES.
!       S0      TOP-OF-ATMOSPHERE SOLAR/LUNAR SPECTRAL
!               IRRADIANCE [W CM-2 / CM-1].
!       UMU0    COSINE OF THE SOLAR/LUNAR ZENITH ANGLE (> 0.).
!       CSSCAT  COSINE OF THE SCATTERING ANGLE.
!       UMU1    COSINE OF THE LINE-OF-SIGHT NADIR ANGLE.
!       PMOM    LEGENDRE POLYNOMIAL MOMENTS; PMOM(1,*) CONTAINS
!               THE SCATTERING PHASE FUNCTION ASYMMETRY FACTORS.
!       SSALB   SINGLE SCATTERING ALBEDOS ORDERED FROM SPACE-TO-GROUND.
!       DTAUC   LAYER VERTICAL OPTICAL DEPTHS FROM SPACE-TO-GROUND.
!       UTAU    CUMULATIVE VERTICAL OPTICAL DEPTHS FROM SPACE-TO-GROUND.
!       UU      LAYER BOUNDARY TO SPACE (GROUND) SPECTRAL RADIANCES,
!               COMPUTED FOR A PLANE-PARALLEL ATMOSPHERE AND ORDERED
!               FROM SPACE-TO-GROUND [W CM-2 SR-1 / CM-1].
      INTEGER MAXCMU,MAXUMU,MAXPHI,NLYR,NTAU
      REAL S0,UMU0,CSSCAT,UMU1,PMOM(0:MAXCMU,1:NLYR),SSALB(NLYR),       &
     &  DTAUC(NLYR),UTAU(NTAU),UU(MAXUMU,NTAU,MAXPHI)

!     OUTPUT ARGUMENTS:
!       S0CMS   AZIMUTHALLY-DEPENDENT LAYER SOLAR SOURCE
!               FUNCTIONS [W CM-2 SR-1 / CM-1].
      REAL S0CMS(MAXUMU,NTAU)

!     FUNCTIONS:
!       HENGNS  HENYEY-GREENSTEIN SCATTERING PHASE FUNCTION [SR-1].
!       PFEXP   LEGENDRE EXPANSION SCATTERING PHASE FUNCTION [SR-1].
      REAL HENGNS,PFEXP

!     LOCAL VARIABLES:
!       N       LAYER INDEX INCREASING SPACE-TO-GROUND.
!       NP1     N PLUS 1.
!       FACTOR  OPTICAL DEPTH SCALING FACTOR REQUIRED FOR LAYER
!               SINGLE SCATTER SOLAR RADIANCE CALCULATION.
!       OD      LAYER OPTICAL DEPTH ALONG THE LINE-OF-SIGHT.
!       ARG     PRODUCT OF FACTOR AND OD.
!       TRANS   LAYER TRANSMITTANCE ALONG THE LINE-OF-SIGHT.
!       EMIS    LAYER EMISSIVITY ALONG THE LINE-OF-SIGHT.
!       DELISS  LAYER CONTRIBUTION TO THE SINGLE SCATTER
!               SOLAR/LUNAR SPECTRAL RADIANCE, ORDERED
!               FROM SPACE-TO-GROUND [W CM-2 SR-1 / CM-1].
!       DELI    LAYER CONTRIBUTION TO THE SPECTRAL RADIANCE, ORDERED
!               FROM SPACE-TO-GROUND [W CM-2 SR-1 / CM-1].
      INTEGER N,NP1
      REAL FACTOR,OD,ARG,TRANS,EMIS,DELISS,DELI

!     BRANCH BASED ON LOOK DIRECTION:
      FACTOR=1.+UMU1/UMU0
      N=1
      IF(UMU1.GT.0.)THEN

!         LOOKING DOWN; LOOP OVER LAYERS:
          DO NP1=2,NTAU

!             N IS THE LAYER BOUNDARY AT THE SENSOR; NP1 IS THE
!             NEXT LAYER BOUNDARY ALONG THE LINE-OF-SIGHT.
              OD=DTAUC(N)/UMU1
              IF(OD.LT.1.E-6)THEN

!                 LAYER OPTICAL DEPTH IS SMALL; SET MULTIPLE
!                 SCATTERING SOLAR SOURCE TERM TO ZERO.
                  S0CMS(1,N)=0.
              ELSE

!                 LAYER TRANSMITTANCE AND EMISSIVITY:
                  IF(OD.LT..001)THEN
                      EMIS=OD*(1.-.5*OD)
                      TRANS=1.-EMIS
                  ELSE
                      TRANS=EXP(-OD)
                      EMIS=1.-TRANS
                  ENDIF

!                 SINGLE SCATTER SOLAR RADIANCE:
                  ARG=FACTOR*OD
                  DELISS=-UTAU(N)/UMU0
                  IF(ARG.LT..001)THEN
                      DELISS=OD*(1-ARG/2)*EXP(DELISS)
                  ELSE
                      DELISS=(EXP(DELISS)-EXP(DELISS-ARG))/FACTOR
                  ENDIF
                  DELISS=S0*SSALB(N)*PFEXP(PMOM(0,N))*DELISS
!bug &                              *HENGNS(PMOM(1,N),CSSCAT)*DELISS

!                 LAYER MULTIPLE SCATTERING SOLAR SOURCE:
                  DELI=UU(1,N,1)-UU(1,NP1,1)*TRANS
                  IF(DELI.LT.1.001*DELISS)THEN

!                     THE SCATTERED SOLAR IS DOMINATED BY SINGLE SCATTER
!                     SO JUST SET THE MULTIPLE SCATTER PART TO ZERO.
                      S0CMS(1,N)=0.
                  ELSE

!                     DIVIDE BY LAYER EMISSIVITY TO OBTAIN THE SOURCE.
                      S0CMS(1,N)=(DELI-DELISS)/EMIS
                  ENDIF
              ENDIF
              N=NP1

!         END LAYER LOOP:
          ENDDO
      ELSEIF(UMU1.LT.0.)THEN

!         LOOKING UP; LOOP OVER LAYERS:
          DO NP1=2,NTAU

!             NP1 IS THE LAYER BOUNDARY AT THE SENSOR; N IS THE
!             NEXT LAYER BOUNDARY ALONG THE LINE-OF-SIGHT.
              OD=-DTAUC(N)/UMU1
              IF(OD.LT.1.E-12)THEN

!                 LAYER OPTICAL DEPTH IS SMALL; SET MULTIPLE
!                 SCATTERING SOLAR SOURCE TERM TO ZERO.
                  S0CMS(1,N)=0.
              ELSE

!                 LAYER TRANSMITTANCE AND EMISSIVITY:
                  IF(OD.LT..001)THEN
                      EMIS=OD*(1-OD/2)
                      TRANS=1.-EMIS
                  ELSE
                      TRANS=EXP(-OD)
                      EMIS=1.-TRANS
                  ENDIF

!                 SINGLE SCATTER SOLAR RADIANCE:
                  ARG=FACTOR*OD
                  DELISS=-UTAU(NP1)/UMU0
                  IF(ABS(ARG).LT..001)THEN
                      DELISS=OD*(1.-.5*ARG)*EXP(DELISS)
                  ELSE
                      DELISS=(EXP(DELISS)-EXP(DELISS-ARG))/FACTOR
                  ENDIF
                  DELISS=S0*SSALB(N)*PFEXP(PMOM(0,N))*DELISS
!bug &                              *HENGNS(PMOM(1,N),CSSCAT)*DELISS

!                 MULTIPLE SCATTERING SOURCE:
                  DELI=UU(1,NP1,1)-UU(1,N,1)*TRANS
                  IF(DELI.LT.1.001*DELISS)THEN

!                     THE SCATTERED SOLAR IS DOMINATED BY SINGLE SCATTER
!                     SO JUST SET THE MULTIPLE SCATTER PART TO ZERO.
                      S0CMS(1,N)=0.
                  ELSE

!                     DIVIDE BY LAYER EMISSIVITY TO OBTAIN THE SOURCE.
                      S0CMS(1,N)=(DELI-DELISS)/EMIS
                  ENDIF
              ENDIF
              N=NP1

!         END LAYER LOOP:
          ENDDO
      ENDIF
      RETURN
      END
