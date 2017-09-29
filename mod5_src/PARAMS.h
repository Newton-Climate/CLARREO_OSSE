
!     MODTRAN INCLUDE PARAMETERS FILE.
!     (NOTE: BEFORE COMPILING, SET MMOLY SAME AS MYMOL IN IFIL.h)

!     LAYER PARAMETERS:
!       LAYDIM   MAXIMUM NUMBER OF LAYER BOUNDARY.
!       LAYDM1   LAYDIM+1
!       LAYTWO   TWICE LAYDIM.
!       NLEVEL   NUMBER OF LAYER BOUNDARIES IN MODEL ATMOSPHERES.
!       NLEVM1   NUMBER OF LAYERS IN MODEL ATMOSPHERE PROFILES - 1.
      INTEGER LAYDIM,LAYDM1,LAYTWO,NLEVEL,NLEVM1
!MCSc PARAMETER(LAYDIM=502,LAYDM1=LAYDIM+1,LAYTWO=2*LAYDIM,             &
      PARAMETER(LAYDIM=152,LAYDM1=LAYDIM+1,LAYTWO=2*LAYDIM,             &
     &  NLEVEL=50,NLEVM1=NLEVEL-1)

!     DISCRETE ORDINATE MULTIPLE SCATTERING PARAMETERS (DISORT NAMES):
!       MXCMU    MAXIMUM NUMBER OF DISORT STREAMS.
!       MAZ      MAXIMUM NUMBER OF DISORT AZIMUTH COMPONENTS.
!       MI       MAXIMUM ORDER OF DOUBLE-GAUSS QUADRATURE.
!       MAXCOE   MAXIMUM NUMBER OF PHASE FUNCTION LEGENDRE COEFFICIENTS.
!       MXUMU    MAXIMUM NUMBER OF DISORT USER NADIR ANGLES.
!       MAXPHI   MAXIMUM NUMBER OF LINE-OF-SIGHT AZIMUTH ANGLES.
!       M_SCAL   NUMBER OF ISAACS_TO_DISORT SCALING SPECTRAL POINTS.
!       MAXMSG   MAXIMUM NUMBER OF MESSAGES WRITTEN TO OUTPUT.
!       MLOS     MAXIMUM NUMBER OF LINE-OF-SIGHT PATHS.
!       MLOSP1   MLOS PLUS ONE.
      INTEGER MXCMU,MAZ,MI,MXUMU,MAXPHI,MAXCOE,M_SCAL,MAXMSG,MLOS,MLOSP1
      PARAMETER(MXCMU=32,MAZ=MXCMU-1,MI=MXCMU/2,MAXCOE=MXCMU,MXUMU=8,   &
     &  MAXPHI=12,MLOS=MXUMU*8,MLOSP1=MLOS+1,M_SCAL=12,MAXMSG=50)

!     SURFACE REFLECTANCE PARAMETERS:
!       MWVSRF   MAXIMUM NUMBER OF SPECTRAL DATA INPUT POINTS.
!                [ARRAYS GAZMOM(1:MI,0:MI,0:MAZ,1:MWVSRF) AND
!                VAZMOM(1:MXUMU,0:MI,0:MAZ,1:MWVSRF) CAN EXCEED
!                AVAILABLE MEMORY IF MXCMU IS TOO LARGE.  SINCE
!                DYNAMIC MEMORY ALLOCATION IS NOT BEING USED,
!                MWVSRF IS REDUCED WHEN MXCMU IS INCREASED.]
!       MBRDF    MAXIMUM NUMBER OF PARAMETES IN BRDF REPRESENTATION.
      INTEGER MWVSRF,MBRDF
!****************** VINCENT ROSS CHANGED FOR SEA BRDF SUPPORT ***********
!Preprocessor directives not supported in a .h file...
      PARAMETER(MWVSRF=((((98304/MXCMU)*32)/MXCMU)*32)/MXCMU,MBRDF=7)
!***************************** END VINCENT ROSS *************************

!     SPECIES NUMBER PARAMETERS:
!       NMOL     NUMBER OF MODTRAN BAND MODEL MOLECULES.
!       NMOLX    NUMBER OF MODTRAN TEMPERATURE-DEPENDENT CROSS-SECTION
!                (X) MOLECULES (H2 AND HE HAVE BEEN ADDED).
!       NMOLXT   NMOL PLUS NMOLX.
!       MMOLY    MAXIMUM NUMBER OF USER-SUPPLIED MOLECULES WITH
!                CROSS-SECTIONS OR BAND MODEL DATA (TERMED Y-SPECIES).
!       MMOLYT   NMOL+MMOLY (MAXIMUM NUMBER OF BAND MODEL MOLECULES).
!       MMOLT    NMOL+NMOLX+MMOLY (MAXIMUM TOTAL NUMBER OF MOLECULES).
!       MEXT     NUMBER OF EXTINCTION SOURCES EXCLUDING TEMPERATURE
!                DEPENDENT CROSS-SECTION SPECIES.
!       MEXTX    NUMBER OF MODTRAN SPECIES (EXTINCTION SOURCES)
!                EXCLUDING Y-SPECIES.
!       NSEG5    NUMBER OF LOW SPECTRAL RESOLUTION (5 CM-1) SPECIES.
      INTEGER NMOL,NMOLX,NMOLXT,MMOLY,MMOLYT,MMOLT,MEXT,MEXTX,          &
     &  MEXTXY,NSEG5
      PARAMETER(NMOL=12,NMOLX=15,NMOLXT=NMOL+NMOLX,MMOLY=50,            &
     &  MMOLYT=NMOL+MMOLY,MMOLT=NMOLXT+MMOLY,MEXT=84,                   &
     &  MEXTX=MEXT+NMOLX,MEXTXY=MEXTX+MMOLY,NSEG5=36)

!     PARAMETERS FOR THE BAND MODEL EXTERNAL FILES:
!       MPRES    MAXIMUM NUMBER OF BAND MODEL PARAMETER PRESSURES;
!       MTEMP    MAXIMUM NUMBER OF BAND MODEL PARAMETER TEMPERATURES.
!       NTEMPX   NUMBER OF MOLECULAR CROSS-SECTION TEMPERATURES.
!       NTLSUB   NUMBER OF LINE TAIL SUBINTERVAL SPECTRAL POINTS LESS 1.
!       RTLSUB   RECIPROCAL OF NTLSUB.
      INTEGER MPRES,MTEMP,NTEMPX,NTLSUB
      REAL RTLSUB
      PARAMETER(MTEMP=12,NTEMPX=6,MPRES=2,NTLSUB=4,RTLSUB=1./NTLSUB)

!     SPECTRAL SLIT FUNCTION PARAMETERS:
!       MWGT     MAXIMUM NUMBER OF SAVED SPECTRAL DATA POINTS.
      INTEGER MWGT
      PARAMETER(MWGT=201)

!     CLOUD/RAIN/AEROSOL PARAMETERS:
!       NZCLD    NUMBER OF CLOUD/RAIN LAYER BOUNDARIES.
!       MAER     NUMBER OF AEROSOL/CLOUD PROFILES.
!       NWAVLN   NUMBER OF WAVELENGTHS IN BUILT-IN SPECTRAL DATA.
!       MXWVLN   MAXIMUM NUMBER OF USER-DEFINED SPECTRAL DATA POINTS.
!       MSSALB   MAXIMUM NUMBER OF AEROSOL ALBEDO SPECTRAL GRID POINTS.
!       ZTOL     TOLERANCE FOR MERGING TOGETHER ALTITUDE LEVELS [KM].
      INTEGER NZCLD,MAER,NWAVLN,MXWVLN,MSSALB
      DOUBLE PRECISION ZTOL
      PARAMETER(NZCLD=21,MAER=21,NWAVLN=788,MXWVLN=788,MSSALB=20,       &
     &  ZTOL=.0005D0)

!     CORRELATED-K APPROACH PARAMETERS:
!       MXKSUB   DIMENSION OF K-DISTRIBUTION SUB-INTERVAL ARRAY.
!       MXGAML   DIMENSION OF LORENTZ HALF-WIDTH ARRAY.
!       MXGAMD   DIMENSION OF DOPPLER HALF-WIDTH ARRAY.
!       MXNUML   DIMENSION OF EFFECTIVE NUMBER OF LINES ARRAY.
      INTEGER MXKSUB,MXGAMD,MXGAML,MXNUML
      PARAMETER(MXKSUB=33,MXGAMD=25,MXGAML=37,MXNUML=15)

!     NOVAM AEROSOLS
!       MNOV     MAXIMUM NUMBER OF ESSENTIAL NOVAM LAYER BOUNDARIES
!       MLNOV    MAXIMUM NUMBER OF ACTUAL NOVAM BOUNDARIES
!                (A PAIR OF VERY CLOSELY SPACED BOUNDARIES OF THE LATTER
!                MAKE UP A LAYER OF THE FORMER)
      INTEGER MNOV,MLNOV
      PARAMETER(MNOV=11,MLNOV=2*MNOV+1)

!     USER-DEFINED PHASE FUNCTION F(MAERF,MANGLS,MWLF):
!       MANGLS   NUMBER OF ANGLES
!       MAERF    NUMBER OF AEROSOLS
!       MWLF     NUMBER OF WAVELENGTHS
      INTEGER MANGLS,MWLF,MAERF
      PARAMETER(MANGLS=94,MWLF=40,MAERF=4)

!     FILE NAME DATA:
!       NAMLEN   MAXIMUM FILE NAME LENGTH.
!       LENSUN   LENGTH OF DEFAULT SOLAR IRRADIANCE DATA FILE.
      INTEGER NAMLEN,LENSUN
      PARAMETER(NAMLEN=256,LENSUN=19)

!     JMASS LOGICAL
!       LJMASS   LOGICAL SWITCH FOR JMASS IF TRUE
      LOGICAL LJMASS
      INTEGER FATAL,WARNING,SUCCES
      PARAMETER(LJMASS=.FALSE.,FATAL=2,WARNING=1,SUCCES=0)

!     CONSTANTS:
!       PZERO    STANDARD PRESSURE [MB].
!       P0TORR   STANDARD PRESSURE [TORR].
!       TZERO    STANDARD TEMPERATURE [K].
!       H2OMWT   WATER MOLECULAR WEIGHT [GM/MOLE].
!       AIRMWT   EFFECTIVE WEIGHT OF AIR [GM/MOLE].
!       AVOGAD   AVOGADRO CONSTANT [MOLECULES/MOLE].
!       STDVOL   STANDARD VOLUME OF IDEAL GAS [CM3 ATM / MOLE].
!       LOSCHM   LOSCHMIDT NUMBER [ATM-1 MOLECULES / CM3].
      REAL PZERO,P0TORR,TZERO,H2OMWT,AVOGAD,STDVOL,LOSCHM
      DOUBLE PRECISION AIRMWT
      PARAMETER(PZERO=1013.25,P0TORR=760.,TZERO=273.15,H2OMWT=18.01528, &
     &  AIRMWT=28.964D0,AVOGAD=6.022045E+23,STDVOL=22413.83,            &
     &  LOSCHM=AVOGAD/STDVOL)

!     NUMERICAL CONSTANTS:
!       IONE     INTEGER ONE
      INTEGER IONE
      PARAMETER(IONE=1)

!     CONVOLUTION PARAMETERS:
!       MSPCDT   MAXIMUM NUMBER OF SPECTRAL POINTS IN *.tp7 (TAPE7).
!       MDDGRD   MAXIMUM NUMBER OF SPECTRAL POINTS IN *.7sc (TAPE7.SCN).
      INTEGER MSPCDT,MDDGRD
      PARAMETER(MSPCDT=273001,MDDGRD=100001)

!     CHANNEL RESPONSE PARAMETERS:
!       MXCHAN   MAXIMUM NUMBER OF CHANNELS READ FROM FILTER
!                FUNCTION FILE.
!       MNBIN    MINIMUM SPECTRAL BIN USED IN CHANNEL INTEGRATIONS
!                (CORRESPONDS TO FREQUENCY MNBIN*BNDWID cm-1).
!       MXBIN    MAXIMUM SPECTRAL BIN USED IN CHANNEL INTEGRATIONS
!                (CORRESPONDS TO FREQUENCY MXBIN*BNDWID cm-1).
!       MXNCHN   MAXIMUM NUMBER OF CHANNELS FOR ANY SPECTRAL BIN.
!       MOUT     MAXIMUM NUMBER OF CHANNEL OUTPUT VARIABLES.
      INTEGER MXCHAN,MNBIN,MXBIN,MXNCHN,MOUT
!sav  PARAMETER(MXCHAN=400,MNBIN=0/15,MXBIN=205,MXNCHN=10)
!sav  PARAMETER(MXCHAN=400,MNBIN=3900/15,MXBIN=28005,MXNCHN=10)
!sav  PARAMETER(MXCHAN=400,MNBIN=765/15,MXBIN=28005,MXNCHN=10)
!sav  PARAMETER(MXCHAN=400,MNBIN=0,MXBIN=28005,MXNCHN=10,               &
!sav &  MOUT=NMOLX+MMOLY+21)
!     MXNCHN increased from 10 to 60 to handle modis399_2176nm.flt:
!     PARAMETER(MXCHAN=400,MNBIN=0,MXBIN=28005,MXNCHN=60,               &
!    &  MOUT=NMOLX+MMOLY+21)
!     Increased MXCHAN from 400 to 2380 and increased MXNCHN
!     from 60 to 92 to handle airs.flt:
      PARAMETER(MXCHAN=2380,MNBIN=0,MXBIN=28005,MXNCHN=92,              &
     &  MOUT=NMOLX+MMOLY+21)

!     CLOUD ARRAY PARAMETERS:
!       NMDLS    NUMBER OF WATER CLOUD/RAIN MODELS.
!       NCLDS    NUMBER OF MODEL WATER CLOUDS.
!       NCIRS    NUMBER OF MODEL CIRRUS CLOUDS.
      INTEGER NMDLS,NCLDS,NCIRS
      PARAMETER(NMDLS=10,NCLDS=5,NCIRS=2)

!     SPECTRAL AEROSOL ARRAY PARAMETERS:
!       MWVSAP   MAXIMUM NUMBER OF AEROSOL SPECTRAL GRID POINTS.
!       MLGSAP   MAXIMUM AEROSOL PHASE FUNCTION LEGENDRE MOMENT.
!       MANSAP   MAXIMUM NUMBER OF AEROSOL PHASE FUNCTION ANGLES.
      INTEGER MWVSAP,MLGSAP,MANSAP
      PARAMETER(MWVSAP=40,MLGSAP=64,MANSAP=181)

!     SOLAR IRRADIANCE PARAMETERS:
!       MSUN15   NUMBER OF SOLAR IRRADIANCES TABULATED ON 15 CM-1 GRID.
!       MSUN05   NUMBER OF SOLAR IRRADIANCES TABULATED ON 5 CM-1 GRID.
!       MSUN01   NUMBER OF SOLAR IRRADIANCES TABULATED ON 1 CM-1 GRID.
!       MSUNP1   NUMBER OF SOLAR IRRADIANCES TABULATED ON 0.1 CM-1 GRID.
      INTEGER MSUN15,MSUN05,MSUN01,MSUNP1
      PARAMETER(MSUN15=3333,MSUN05=10000,MSUN01=50000,MSUNP1=500000)

!     BINARY OUTPUT:
!       BUFFSZ   BUFFER SIZE.
!       BINSPR   FILE SEPARATOR.
      INTEGER BUFFSZ,BINSPR
      PARAMETER(BUFFSZ=2046,BINSPR=-9999)

!     ASCII OUTPUT PARAMETER:
!       TXTSPR   SEPARATOR
      CHARACTER TXTSPR*7
      PARAMETER(TXTSPR=' -9999.')

!     SHORT RANGE PARAMETER:
!       RSMALL   SHORT RANGE CUTOFF [KM].
      DOUBLE PRECISION RSMALL
      PARAMETER(RSMALL=2.D0)
