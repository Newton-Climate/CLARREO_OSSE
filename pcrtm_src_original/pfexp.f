      REAL FUNCTION PFEXP(ILOS,PMOM)

!     CALCULATES SCATTERING PHASE FUNCTION FROM LEGENDRE EXPANSION.

!     INPUTS:
!       ILOS     LOOP INDEX FOR LINE-OF-SIGHT PATHS.
!       PMOM     PHASE FUNCTION LEGENDRE EXPANSION COEFFICIENTS.
      INTEGER ILOS
      REAL PMOM(0:*)

!     OUTPUTS:
!       PFLEG    SCATTERING PHASE FUNCTION [STER-1]

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

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

!     FUNCTIONS:
!bug  REAL HENGNS

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       ISTR     STREAM NUMBER INDEX.
!       PFLEG    LEGENDRE EXPANSION WITH NO HIGHER ORDER TERMS [SR-1].
!       DELTAM   DELTA-M LEGENDRE EXPANSION [SR-1].
!bug    HGSUM    PERTURBED HENYEY-GREENSTEIN LEGENDRE EXPANSION [SR-1].
!bug    GN       ASYMMETRY FACTOR RAISED TO THE Nth POWER [SR-1].
!bug  REAL HGSUM,GN
      REAL PFLEG,DELTAM
      INTEGER ISTR

!     DIRECT LEGENDRE EXPANSION:
      PFLEG=SLEGEN(0,ILOS)+PMOM(1)*SLEGEN(1,ILOS)

!     LEGENDRE EXPANSION AS PERTURBED HENYEY-GREENSTEIN:
!bug  HGSUM=HENGNS(PMOM(1),CSSCAT)
!bug  GN=PMOM(1)

!     DELTA-M LEGENDRE EXPANSION:
      DELTAM=PFLEG-PMOM(NSTR)*(SLEGEN(0,ILOS)+SLEGEN(1,ILOS))

      DO ISTR=2,NSTR-1
          PFLEG=PFLEG+PMOM(ISTR)*SLEGEN(ISTR,ILOS)
!bug      GN=GN*PMOM(1)
!bug      HGSUM=HGSUM+(PMOM(ISTR)-GN)*SLEGEN(ISTR,ILOS)
          DELTAM=DELTAM+(PMOM(ISTR)-PMOM(NSTR))*SLEGEN(ISTR,ILOS)
      ENDDO
      PFEXP=DELTAM
      RETURN
      END
