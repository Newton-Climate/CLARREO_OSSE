
!     MODTRAN INCLUDE PARAMETERS FILE.

!     LAYER PARAMETERS:
!       LAYDIM   MAXIMUM NUMBER OF LAYER BOUNDARY.
!       LAYTWO   TWICE LAYDIM.
!       LAYTHR   THREE TIMES LAYDIM.
!       NLAYX    CROSS-SECTION SPECIES PROFILE LAYER BOUNDARIES NUMBER.
      INTEGER LAYDIM,LAYTWO,LAYTHR,NLAYX
      PARAMETER(LAYDIM=502,LAYTWO=2*LAYDIM,LAYTHR=3*LAYDIM,NLAYX=50)

!     DISCRETE ORDINATE MULTIPLE SCATTERING PARAMETERS (DISORT NAMES):
!       MXCMU   MAXIMUM NUMBER OF DISORT STREAMS.
!       MAZ     MAXIMUM NUMBER OF DISORT AZIMUTH COMPONENTS.
!       MI      MAXIMUM ORDER OF DOUBLE-GAUSS QUADRATURE.
!       MXUMU   MAXIMUM NUMBER OF DISORT USER ANGLES.
      INTEGER MXCMU,MAZ,MI,MXUMU
      PARAMETER(MXCMU=16,MAZ=MXCMU-1,MI=MXCMU/2,MXUMU=2)

!     SPECIES NUMBER PARAMETERS:
!       NMOL     NUMBER OF MODTRAN MOLECULAR BAND MODEL PARAMETER
!                SPECIES IN THE BAND MODEL FILE.
!       NMOLX    NUMBER OF TEMPERATURE-DEPENDENT CROSS-SECTION SPECIES
!                (TERMED X-SPECIES, CODING IS DONE ASSUMING THAT FULL
!                 BAND MODEL PARAMETERS MAY BE AVAILABLE IN THE FUTURE).
!       NMOLXT   NMOL PLUS NMOLX.
!       MEXT     NUMBER OF EXTINCTION SOURCES EXCLUDING TEMPERATURE
!                DEPENDENT CROSS-SECTION SPECIES.
!       MEXTX    NUMBER OF MODTRAN SPECIES (EXTINCTION SOURCES).
!       NLAY5    NUMBER OF LOW SPECTRAL RESOLUTION (5 CM-1) SPECIES.
      INTEGER NMOL,NMOLX,NMOLXT,MEXT,MEXTX,NLAY5
      PARAMETER(NMOL=12,NMOLX=13,NMOLXT=NMOL+NMOLX,                     &
     &  MEXT=83,MEXTX=MEXT+NMOLX,NLAY5=33)

!     PARAMETERS FOR THE BAND MODEL EXTERNAL FILES:
!       MXTEMP   MAXIMUM NUMBER OF BAND MODEL PARAMETER TEMPERATURES.
!       MTLSUB   MAXIMUM NUMBER OF LINE TAIL PARAMETERS PER BAND.
      INTEGER MXTEMP,MTLSUB
      PARAMETER(MXTEMP=6,MTLSUB=4)

!     SPECTRAL SLIT FUNCTION PARAMETERS:
!       NBINS    MAXIMUM NUMBER OF SAVED SPECTRAL DATA POINTS.
      INTEGER NBINS
      PARAMETER(NBINS=99)

!     CLOUD/RAIN/AEROSOL PARAMETERS:
!       NZCLD    NUMBER OF CLOUD/RAIN LAYER BOUNDARIES.
!       MAER     NUMBER OF AEROSOL/CLOUD PROFILES.
!       NWAVLN   NUMBER OF WAVELENGTHS IN BUILT-IN SPECTRAL DATA.
!       MXWVLN   MAXIMUM NUMBER OF USER-DEFINED SPECTRAL DATA POINTS.
      INTEGER NZCLD,MAER,NWAVLN,MXWVLN
      PARAMETER(NZCLD=26,MAER=17,NWAVLN=788,MXWVLN=788)

!     CORRELATED-K APPROACH PARAMETERS:
!       MXKSUB   DIMENSION OF K-DISTRIBUTION SUB-INTERVAL ARRAY.
!       MXGAML   DIMENSION OF LORENTZ HALF-WIDTH ARRAY.
!       MXGAMD   DIMENSION OF DOPPLER HALF-WIDTH ARRAY.
!       MXNUML   DIMENSION OF EFFECTIVE NUMBER OF LINES ARRAY.
      INTEGER MXKSUB,MXGAMD,MXGAML,MXNUML
      PARAMETER(MXKSUB=33,MXGAMD=23,MXGAML=36,MXNUML=14)

!     NOVAM AEROSOLS
!       MNOV     MAXIMUM NUMBER OF ESSENTIAL NOVAM LAYER BOUNDARIES
!       MLNOV    MAXIMUM NUMBER OF ACTUAL NOVAM BOUNDARIES
!                (A PAIR OF VERY CLOSELY SPACED BOUNDARIES OF THE LATTER
!                MAKE UP A LAYER OF THE FORMER)
      INTEGER MNOV,MLNOV
      PARAMETER(MNOV=11,MLNOV=2*MNOV)

!     USER-DEFINED PHASE FUNCTION F(MAERF,MANGLS,MWLF):
!       MANGLS   NUMBER OF ANGLES
!       MAERF    NUMBER OF AEROSOLS
!       MWLF     NUMBER OF WAVELENGTHS
      INTEGER MANGLS,MWLF,MAERF
      PARAMETER(MANGLS=51,MWLF=15,MAERF=4)

!     SURFACE REFLECTANCE PARAMETERS:
!       MWVSRF   MAXIMUM NUMBER OF SPECTRAL DATA INPUT POINTS.
      INTEGER MWVSRF
      PARAMETER(MWVSRF=3000)

!     FILE NAME DATA:
!       NAMLEN   MAXIMUM FILE NAME LENGTH.
      INTEGER NAMLEN
      PARAMETER(NAMLEN=256)