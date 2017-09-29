      SUBROUTINE SHADE(IPH,ILOS,IK,IPATH,JNTRVL,VCEN,                   &
     &  PFMOL0,PFMOL2,TSNOBS,TSNREF,SUMSSS,RADCUM)

!     ROUTINE SHADE SETS THE SOLAR TRANSMISSION TERMS TO ZERO FOR
!     SHADED SCATTERING POINTS
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       IPH      SCATTERING PHASE FUNCTION SWITCH.
!                =2 FOR LOWTRAN MIE FUNCTIONS
!       ILOS     PATH INDEX.
!       IK       LAYER INDEX.
!       IPATH    PATH TYPE SWITCH.
!                =1 FOR DIRECT SUN TO OBSERVER PATH ONLY
!                =2 FOR SUN TO SCATTERING POINT PATH
!                =3 FOR OBSERVER TO SCATTERING POINT LINE-OF-SIGHT PATH
!       JNTRVL   NUMBER OF K-DISTRIBUTION INTERVALS FOR CURRENT FREQ.
!       VCEN     SPECTRAL BIN CENTRAL FREQUENCY [CM-1].
!       PFMOL0   MOLECULAR PHASE FUNCTION ISOTROPIC SCATTERING TERM.
!       PFMOL2   MOLECULAR PHASE FUNCTION ANISOTROPIC SCATTERING TERM.
      INTEGER IPH,ILOS,IK,IPATH,JNTRVL
      REAL VCEN,PFMOL0,PFMOL2

!     OUTPUT ARGUMENTS:
!       TSNOBS   SOLAR/LUNAR IRRADIANCE AT THE OBSERVER [W CM-2 / CM-1].
!       TSNREF   SOLAR/LUNAR IRRADIANCE ALONG L-SHAPED PATH.
!       SUMSSS   SINGLE SCATTER SOLAR RADIANCE SUM [W SR-1 CM-2 / CM-1].
!       RADCUM   CUMULATIVE PATH RADIANCE K-DATA [W SR-1 CM-2 / CM-1].
      REAL TSNOBS,TSNREF,SUMSSS,RADCUM(MXKSUB)

!     COMMONS:
      INCLUDE 'BASE.h'

!     /CORKDT/
!       WTKSUB   SPECTRAL BIN SUB-INTERVAL FRACTIONAL WIDTHS.
!       WTKSAV   SAVED SPECTRAL BIN SUB-INTERVAL FRACTIONAL WIDTHS.
!       DEPLAY   INCREMENTAL OPTICAL DEPTHS.
!       TRNLAY   INCREMENTAL TRANSMITTANCES.
!       TRNCUM   CUMULATIVE TRANSMITTANCES.
!       K2TAIL   POINTER FROM K BIN TO LINE TAIL SUB-BIN
!                (=0 IF MULTIPLE LINE TAIL SUB-BINS CONTRIBUTE).
!       CONTWT   WEIGHTS FOR PARTITIONING LINE TAILS INTO K'S
!                (ONLY USED IF K2TAIL IS 0).
      REAL WTKSUB(MXKSUB),WTKSAV(NTLSUB),DEPLAY(MXKSUB),                &
     &  TRNLAY(MXKSUB),TRNCUM(MXKSUB),CONTWT(NTLSUB,MXKSUB)
      INTEGER K2TAIL(MXKSUB)
      COMMON/CORKDT/K2TAIL,WTKSUB,WTKSAV,DEPLAY,TRNLAY,TRNCUM,CONTWT
      SAVE /CORKDT/

!     /RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

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
      REAL CSSCAT,SLEGEN,CSZEN0,CSZEN,CSZENX,TCONT,TAUT,GTSCAT,COSBAR,  &
     &  BBGRND,BBNDRY,S0DEP,S0TRN,DEPRAT,UPF,DNF,UPFS,DNFS
      LOGICAL CO_LIN
      COMMON/MSRD/CSSCAT(MLOS),SLEGEN(0:MAZ,MLOS),CSZEN0(LAYDIM),       &
     &  CSZEN(LAYDIM),CSZENX(LAYDIM),TCONT(LAYDIM),TAUT(MXKSUB,LAYDIM), &
     &  GTSCAT(0:MXCMU,1:LAYDIM),COSBAR(LAYDIM),BBGRND,BBNDRY(LAYDIM),  &
     &  S0DEP(MXKSUB,LAYTWO),S0TRN(MXKSUB,LAYTWO),DEPRAT(MXKSUB,LAYDIM),&
     &  UPF(MXKSUB,LAYDIM),DNF(MXKSUB,LAYDIM),UPFS(MXKSUB,LAYDIM),      &
     &  DNFS(MXKSUB,LAYDIM),CO_LIN(MLOS)

!     LOCAL VARIABLES:
!       INTRVL   LOOP INDEX FOR K-DISTRIBUTIONS.
!       SSAP     SCATTERING OPTICAL DEPTH OF SPECTRAL AEROSOL PROFILE.
      INTEGER INTRVL
      REAL S0,TMOL,BIGDEP,SSAP
      S0=0.
      IF(JNTRVL.EQ.1)THEN
          TMOL=0.
          CALL SSRAD(IPH,ILOS,IK,IPATH,VCEN,S0,PFMOL0,PFMOL2,           &
     &      TMOL,SSAP,TX(1),TSNOBS,TSNREF,SUMSSS)
      ELSE

!         DIVIDE BIGNUM BY THE LAYER NUMBER TO ARTIFICIALLY REQUIRE
!         THAT OPTICAL DEPTH DECREASES WITH INCREASING LAYER NUMBER.
          BIGDEP=BIGNUM/IK
          DO INTRVL=1,JNTRVL
              TRNLAY(INTRVL)=0.
              DEPLAY(INTRVL)=BIGDEP
          ENDDO
          IF(ILOS.EQ.MLOSP1)THEN

!             FOR MULTIPLE SCATTERING, STORE THE DECREASE IN WEAK-LINE
!             OPTICAL DEPTH TO THE SUN ACROSS THE CURRENT LAYER.
              IF(IPATH.EQ.1)THEN
                  DO INTRVL=1,JNTRVL
                      DEPRAT(INTRVL,IK)=BIGDEP
                  ENDDO
              ELSE
                  DO INTRVL=1,JNTRVL
                      DEPRAT(INTRVL,IK+1)=BIGDEP
                      DEPRAT(INTRVL,IK)=1-BIGDEP/DEPRAT(INTRVL,IK)
                  ENDDO
              ENDIF
          ENDIF
          CALL SSCORK(.FALSE.,IPH,ILOS,IK,IPATH,JNTRVL,VCEN,S0,         &
     &      PFMOL0,PFMOL2,SSAP,TX(1),TSNOBS,SUMSSS,RADCUM)
      ENDIF
      RETURN
      END
