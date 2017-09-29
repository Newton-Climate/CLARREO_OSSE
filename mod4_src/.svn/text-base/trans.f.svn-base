      SUBROUTINE TRANS(IPH,ISOURC,IDAY,ANGLEM,                          &
     &  GROUND,IMSMX,KNTRVL,LNFILT,SUNFUN,APPREF,WV_HRES,               &
     &  RADIANCE_HRES,NUM_HRES,SOLAR_FLUX,NUM_LRES,DIFFUSE_FLUX,LNFLRT, &
     &  FLRT,NUM_LEVS,BB_UPDIFFUSE,BB_DNDIFFUSE,BB_DNDIRECT,            &
     &  FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)

!     ROUTINE TRANS CALCULATES TRANSMITTANCE AND RADIANCE VALUES
!     BETWEEN IV1 AND IV2 FOR A GIVEN ATMOSPHERIC SLANT PATH.

!     ARGUMENTS:
!       LNFILT   LENGTH OF FILTER FUNCTION FILE NAME
!                (>0 IF FILTER FUNCTION FILE IS TO BE USED)
!       SUNFUN   TOP-OF-ATMOSPHERE SOLAR IRRADIANCE
!                FUNCTION [W CM-2 / CM-1].
      LOGICAL GROUND,APPREF
      INTEGER IPH,ISOURC,IDAY,IMSMX,KNTRVL,LNFILT
      REAL ANGLEM,SUNFUN
      EXTERNAL SUNFUN
      INTEGER NUM_HRES !DRF
      INTEGER NUM_LRES !DRF
      INTEGER NUM_LEVS !DRF 
      REAL WV_HRES(NUM_HRES), RADIANCE_HRES(NUM_HRES) !DRF
      REAL SOLAR_FLUX(NUM_LRES),DIFFUSE_FLUX(NUM_LRES) !DRF
      INTEGER LNFLRT
      CHARACTER FLRT*(*)
      REAL BB_UPDIFFUSE(NUM_LEVS)
      REAL BB_DNDIFFUSE(NUM_LEVS) 
      REAL BB_DNDIRECT(NUM_LEVS)
      REAL FLX_UPDIFFUSE(NUM_LEVS,NUM_LRES)
      REAL FLX_DNDIFFUSE(NUM_LEVS,NUM_LRES) 
      REAL FLX_DNDIRECT(NUM_LEVS,NUM_LRES)

!     LIST PARAMETERS
      INCLUDE 'PARAMS.h'
      INCLUDE 'JMASS.h'
      INCLUDE 'ERROR.h'
      INTEGER IPRINT,MAXV,MDISCL
      PARAMETER(IPRINT=50,MDISCL=12,MAXV=50000)

!     INCREASE FIRST DIMENSION IN SLIT TO MEXT TO INSURE THAT ALL
!     TRANSMITTANCES (TX ARRAY) ARE PASSED THROUGH THE SLIT FUNCTION.
      REAL SLIT(MEXT,NBINS)

!     PARAMETER NMOLX IS THE NUMBER OF CROSS-SECTION SPECIES.
!     PARAMETER MEXT DENOTES THE NUMBER OF MODTRAN "SPECIES".
!     THIS INCLUDES THE 12 ORIGINAL BAND MODEL PARAMETER MOLECULES
!     PLUS A HOST OF OTHER ABSORPTION AND/OR SCATTERING SOURCES.

!     COMMONS:

!     COMMON /CFLAGS/
!       YFLAG    Y COORDINATE FLAG FOR PLOT.DAT FILE
!                  = "T" FOR TRANSMITTANCE
!                  = "R" FOR RADIANCE (IRRADIANCE FOR IEMSCT=3)
!                  = "N" FOR NO PLOT.DAT OUTPUT
!       XFLAG    X COORDINATE FLAG FOR PLOT.DAT FILE
!                  = "W" FOR FREQUENCY IN WAVENUMBERS (CM-1) AND
!                        RADIANCE IN W SR-1 CM-2 / CM-1
!                  = "M" FOR WAVELENGTH IN MICRONS AND
!                        RADIANCE IN W SR-1 CM-2 / MICRON
!                  = "N" FOR WAVELENGTH IN NANOMETERS AND
!                        RADIANCE IN MICRO-WATTS SR-1 CM-2 / NANOMETER
!       DLIMIT   DELIMITER CHARACTER STRING BETWEEN MODTRAN RUNS
!       FLAGS    SCANNING FUNCTION FLAGS.
      CHARACTER YFLAG*1,XFLAG*1,DLIMIT*8,FLAGS*7
      COMMON/CFLAGS/YFLAG,XFLAG,DLIMIT,FLAGS
      CHARACTER*8 CNAMEX
      COMMON/NAMEX/CNAMEX(NMOLX)
      SAVE /NAMEX/
!NAMX EXTERNAL XMLATM
      LOGICAL IVTEST,LOOP0
      INCLUDE 'IFIL.h'
      INCLUDE 'BASE.h'

!     /CARD1/
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,MODTRN

!     /SRFACE/
!       NWVSRF  NUMBER OF WAVELENGTH GRID POINTS.
!       WVSURF  WAVELENGTH GRID FOR TARGET-PIXEL (1) AND
!               AREA-AVERAGE (2) SURFACES [MICRONS].
!       SALB    AREA-AVERAGE GROUND SURFACE ALBEDO SPECTRAL ARRAY.
!       GDIRRF  AREA-AVERAGE GROUND SURFACE DIRECTIONAL REFLECTIVITY
!               AT SOLAR ANGLE.
!       HDIR    TARGET-PIXEL HEMISPHERE DIRECTIONAL REFLECTANCES
!               AT VIEWING ANGLE.
!       BRDF    TARGET-PIXEL BIDIRECTIONAL REFLECTANCE DISTRIBUTION
!               FUNCTION AT VIEWING, SUN AND RELATIVE AZIMUTH ANGLES.
!       GDIREM  AREA-AVERAGED GROUND DIRECTIONAL EMISSIVITIES
!               AT VIEWING AND GAUSSIAN QUADRATURE ANGLES.
!       GNDMOM  AREA-AVERAGED GROUND BRDF AZIMUTH FOURIER MOMENTS.
      INTEGER NWVSRF
      REAL WVSURF,SALB,GDIRRF,HDIR,BRDF,GDIREM,GNDMOM
      COMMON/SRFACE/NWVSRF(2),WVSURF(2,MWVSRF),SALB(MWVSRF),            &
     &  GDIRRF(MWVSRF),HDIR(MWVSRF),BRDF(MWVSRF),                       &
     &  GDIREM(0:MI,1:MWVSRF),GNDMOM(0:MI,0:MI,0:MAZ,1:MWVSRF)

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

!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT
      INCLUDE 'SOLS.h'

!     /SAVEMS/
!       LUSEMS  LOGICAL, TRUE IF MULTIPLE SCATTERING DATA IS REUSED.
!       LSAVMS  LOGICAL, TRUE IF MULTIPLE SCATTERING DATA IS  SAVED.
      LOGICAL LUSEMS,LSAVMS
      COMMON/SAVEMS/LUSEMS,LSAVMS
      INCLUDE 'BMHEAD.h'

!     /M3DSPC/
!       LASCII   ASCII FILE OUTPUT FLAG FOR MOD3D FILES.
!       NATM     COUNTER FOR NUMBER OF DIFFERENT MOLECULAR ATMOSPHERES.
!       IVXMIN   MINIMUM COMPUTATION SPECTRAL FREQUENCY [CM-1].
!       IVXMAX   MAXIMUM COMPUTATION SPECTRAL FREQUENCY [CM-1].
      LOGICAL LASCII
      INTEGER NATM,IVXMIN,IVXMAX
      COMMON/M3DSPC/LASCII,NATM,IVXMIN,IVXMAX

!     /SURFWV/
!       LAMBER  LOGICAL FLAG, .TRUE. FOR LAMBERTIAN SURFACE.
!       TPTEMP  TARGET-PIXEL SURFACE TEMPERATURES [K].
!       TPHDIR  TARGET-PIXEL HEMISPHERE DIRECTIONAL REFLECTANCE AT
!               VIEWING ANGLE.
!       TPBRDF  TARGET-PIXEL BIDIRECTIONAL REFLECTANCE DISTRIBUTION
!               FUNCTION AT VIEWING AND SUN ANGLE.
!       AATEMP  AREA-AVERAGED GROUND SURFACE TEMPERATURES [K].
!       AASALB  AREA-AVERAGED GROUND SURFACE ALBEDO.
!       AADREF  AREA-AVERAGED GROUND SURFACE DIRECTIONAL REFLECTIVITY
!               AT THE SOLAR ZENITH ANGLE.
!       EMU     GROUND DIRECTIONAL EMISSIVITY AT VIEWING ANGLE.
!       BEM     GROUND DIRECTIONAL EMISSIVITY AT QUADRATURE ANGLE.
!       RMU     GROUND BRDF AZIMUTH COMPONENTS AT VIEWING ANGLE
!               AND AT SUN (=0) OR QUADRATURE (>0) ANGLE.
!       BDR     GROUND BRDF AZIMUTH COMPONENTS AT QUADRATURE ANGLE
!               AND AT SUN (=0) OR QUADRATURE (>0) ANGLE.
      LOGICAL LAMBER
      REAL TPTEMP,TPHDIR,TPBRDF,AATEMP,AASALB,AADREF,EMU,BEM,RMU,BDR
      COMMON/SURFWV/LAMBER,TPTEMP,TPHDIR,TPBRDF,AATEMP,AASALB,AADREF,   &
     &  EMU(MXUMU),BEM(MI),RMU(1:MXUMU,0:MI,0:MAZ),BDR(1:MI,0:MI,0:MAZ)

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

!     /DISCAL/
      REAL VDISCL(MDISCL),RDISCL(2,MDISCL)
      INTEGER NDISCL
      COMMON/DISCAL/NDISCL,VDISCL,RDISCL

!     /CJM5/
!       AMOD3D   FLAG INDICATING OUTPUT DATABASE FILE TYPE:
      CHARACTER*1 AMOD3D
      COMMON/CJM5/AMOD3D

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:
      REAL BBFN,SOURCE

!     LOCAL VARIABLES:
!       TRANSM  FLAG [TRUE FOR TRANSMITTANCE ONLY CALCULATIONS].
!       LSURF   SPECTRAL FLAG [TRUE IF REFLECTANCES ARE CHANGING].
!       SINIT   FLAG, TRUE FOR INITIAL CALL TO SOURCE FROM TRANS.
!       IWIDTH  FWHM IN BAND WIDTHS (=IFWHM/IBNDWD).
!       WNORM   SLIT FUNCTION NORMALIZATION FACTOR (= 1/IWIDTH**2).
!       DEMIS   DIRECTIONAL EMISSIVITY CONVOLVED WITH SLIT FUNCTION.
!       disalbsv saves DISORT logical disalb for possible repeat run.
!       disazmsv saves DISORT logical disazm for possible repeat run.
!       dissv    saves DISORT logical dis for possible repeat run.
!       n2gaussv saves DISORT variable n2gaus for possible repeat run.
!       nstrsv   saves DISORT variable nstr for possible repeat run.
!       nazsv    saves DISORT variable naz for possible repeat run.
      LOGICAL TRANSM,LSURF,SINIT,DISALBSV,DISAZMSV,DISSV
      INTEGER I,J,IKMX,IDV5,IV,IVXBIN,IWIDTH,                           &
     &  NWGTM1,K,IP1,ICOUNT,IVRMIN,IVRMAX,JLO,NAZSV,NSTRSV,N2GAUSSV
      REAL YPLTMX,RADMIN,RADMAX,SUMA,RADSUM,SSOL,BBG,SUMSSS,SUMSSR,     &
     &  RFLSS,RFSURF,S0,TS0,FACTOR,SUMTMS,SUMMS,UNIF,TRACE,TRANSX,      &
     &  RADCUM,V,CONVRT,WNORM,STORE,ALAM,ALTX9,SUMT,TSNOBS,TSNREF,      &
     &  FDNSRT,FDNTRT,DEMIS,SCALE,FRAC
      INTEGER ICOUNT2

      IF(LDISCL)THEN
         DISAZMSV=DISAZM
         DISALBSV=DISALB
         NAZSV=NAZ
         N2GAUSSV=N2GAUS
         NSTRSV=NSTR
         DISSV=DIS
         CALL TNOWRT(IPH,ISOURC,IDAY,ANGLEM,                            &
     &     GROUND,IMSMX,KNTRVL,SUNFUN,APPREF)
      ENDIF

!     INITIALIZE PLOT.DAT FILE MAXIMUM.
      YPLTMX=0.

!     INCREASE DIMENSION IN SLIT FUNCTION INITIALIZATION
      DO 10 I=1,MEXT
          DO 10 J=1,NBINS
   10 SLIT(I,J)=0.

!     INITIALIZE RADIANCE MINIMUM AND MAXIMUM PARAMETERS
      RADMIN=BIGNUM
      RADMAX=0.

!     STORE THE NUMBER OF PATH LAYERS IN IKMX
      IKMX=IKMAX

!     INITIALIZE INTEGRATED ABSORPTION, RADIANCE, SOLAR IRRADIANCE SUMS.
      SUMA=0.
      RADSUM=0.
      SSOL=0.

!     INITIALIZE TRANSMITTANCE/RADIANCE/IRRADIANCE TERMS:
      UNIF=0.
      TRACE=0.
      BBG=0.
      SUMSSS=0.
      SUMSSR=0.
      RFLSS=0.
      RFSURF=0.
      S0=0.
      TS0=0.
      TSNOBS=0.
      TSNREF=0.
      SUMTMS=0.

!DRF  INITIALIZE ICOUNT 2
      ICOUNT2=0
!     INITIALIZE INTEGRATION WEIGHTING FACTOR
      FACTOR=.5

!     INITIALIZE ICOUNT, USED TO DETERMINE WHEN HEADER MUST BE PRINTED
      NSPCDT=0
      IF(LJMASS)THEN

!         NEVER PRINT HEADER FOR JMASS:
          ICOUNT=IPRINT+1
      ELSE

!         PRINT HEADER AT FIRST FREQUENCY:
          ICOUNT=IPRINT
      ENDIF

!     DEFINE SPECTRAL INCREMENT VARIABLES:
!       IDV5     INCREMENT FOR RETRIEVING 5 CM-1 RESOLUTION DATA [CM-1].
!       IV       INITIAL FREQUENCY USED IN CALCULATIONS [CM-1].
      IDV5=5*((IBNDWD+4)/5)
      IWIDTH=IFWHM/IBNDWD
      IVOFF=(IWIDTH-1)*IBNDWD
      IV=5*((IV1-IVOFF)/5)
      IF(IV.LT.0)IV=0

!     CALL ROUTINE "BMDATA" TO INITIALIZE BAND MODEL PARAMETERS.
      IF(MODTRN .AND. IV1.LE.MXFREQ)CALL BMDATA(IV,IKMX,IMSMX)

!     OFFSET FREQUENCY VARIABLES IVX AND IV SO THAT
!     WHEN THEY ARE INCREMENTED THEY BEGIN IN SINC.
      IVXMIN=IV
      IVX=IV-IBNDWD
      IVXBIN=IV/IBNDWD-1
      IV=IV-IDV5
      IVXMAX=IV2+IVOFF
      IF(IVXMAX.GT.MAXV)IVXMAX=MAXV
      IVXMAX=(IVXMAX/IBNDWD)*IBNDWD
      IWRITE=IV1+IVOFF

!     PERFORM TRIANGULAR SLIT INITIALIZATION.  TRANSMITTANCES AT A
!     GIVEN FREQUENCY CONTRIBUTE TO 2*IWIDTH-1 TRIANGULAR SLITS.
!     THESE CONTRIBUTIONS ARE STORED IN ARRAY SLIT.
      NWGT=2*IWIDTH
      WNORM=FLOAT(IWIDTH*IWIDTH)
      DO 20 I=1,IWIDTH
          WGT(I)=I/WNORM
   20 WGT(NWGT-I)=WGT(I)
      NWGT=NWGT-1
      NWGTM1=NWGT-1

!     INITIALIZE SINIT (= 0 FOR INITIAL CALL TO ROUTINE SOURCE)
      SINIT=.TRUE.

!     INITIALIZE TRANSMITTANCE AND SURFACE REFLECTANCE FLAGS:
      IF(IEMSCT.EQ.1 .OR. IEMSCT.EQ.2)THEN
          TRANSM=.FALSE.
          LSURF=.TRUE.
      ELSE
          TRANSM=.TRUE.
      ENDIF

!     PRINT HEADERS IF NOT JMASS:
      IF(.NOT.LJMASS)THEN
          IF(IEMSCT.EQ.0)THEN
              WRITE(IPU,'(100A)')'    FREQ COMBIN    H2O   UMIX',       &
     &          '     O3  TRACE     N2    H2O MOLEC AER+CLD',           &
     &          '  HNO3 AER+CLD    -LOG    CO2     CO    CH4    N2O',   &
     &          '     O2    NH3     NO    NO2    SO2',                  &
     &          (CNAMEX(K)(2:8),K=1,NMOLX)
              WRITE(IPU,'(100A)')'    CM-1  TRANS  TRANS  TRANS',       &
     &          '  TRANS  TRANS   CONT   CONT   SCAT  TRANS',           &
     &          '  TRANS abTRNS  COMBIN',('  TRANS',K=1,NMOLX+9)
              WRITE(IPUSCR,'(100A)')'    FREQ COMBIN    H2O   UMIX',    &
     &          '     O3  TRACE     N2    H2O MOLEC AER+CLD',           &
     &          '  HNO3 AER+CLD    -LOG    CO2     CO    CH4    N2O',   &
     &          '     O2    NH3     NO    NO2    SO2',                  &
     &          (CNAMEX(K)(2:8),K=1,NMOLX)
              WRITE(IPUSCR,'(100A)')'    CM-1  TRANS  TRANS  TRANS',    &
     &          '  TRANS  TRANS   CONT   CONT   SCAT  TRANS',           &
     &          '  TRANS abTRNS  COMBIN',('  TRANS',K=1,NMOLX+9)
          ELSEIF(IEMSCT.EQ.3)THEN
              WRITE(IPU,   '(A)')'    FREQ   TRANS     SOL TR  SOLAR'
              WRITE(IPUSCR,'(A)')'    FREQ   TRANS     SOL TR  SOLAR'
          ELSE
              WRITE(IPU,   '(3A)')'    FREQ  TOT TRANS  PTH THRML  ',   &
     &          'THRML SCT  SURF EMIS   SOL SCAT  SING SCAT  GRND RFLT',&
     &          '  DRCT RFLT  TOTAL RAD  REF SOL  SOL@OBS   DEPTH'
              WRITE(IPUSCR,'(3A)')'    FREQ  TOT TRANS  PTH THRML  ',   &
     &          'THRML SCT  SURF EMIS   SOL SCAT  SING SCAT  GRND RFLT',&
     &          '  DRCT RFLT  TOTAL RAD  REF SOL  SOL@OBS   DEPTH'
          ENDIF
          IF(NOPRNT.LE.-1)THEN
              IF(DIS .AND. .NOT.DISAZM)THEN
                  WRITE(IPR1,'((3A))')                                  &
     &              ' FREQ     ALT     TOTAL    DELTA   UP-DN-DIR',     &
     &              '  DIFFUS_UP  DIFFUS_DN  THRML_SRC  THRML_SUM',     &
     &              '   NOT_USED   NOT_USED  SOLAR_SRC  SOLAR_SUM',     &
     &              '(CM-1)   (KM)     TRANS    TRANS (          ',     &
     &              '                                 W CM-2 / CM',     &
     &              '-1                                         )'
              ELSEIF(IMULT.NE.0)THEN
                  WRITE(IPR1,'((3A))')                                  &
     &              ' FREQ     ALT     TOTAL    DELTA   UP-DN-DIR',     &
     &              '   THRML_UP   THRML_DN  THRML_SRC  THRML_SUM',     &
     &              '   SOLAR_UP   SOLAR_DN  SOLAR_SRC  SOLAR_SUM',     &
     &              '(CM-1)   (KM)     TRANS    TRANS (          ',     &
     &              '                                 W CM-2 / CM',     &
     &              '-1                                         )'
              ELSEIF(IEMSCT.GT.0)THEN
                  WRITE(IPR1,'((1X,2A))')                               &
     &              '             ALTITUDES                 B(V,T)  ',  &
     &              '           TRANSMISSION            RADIANCE',      &
     &              '  FREQ  BEGINNING  ENDING  INT   LAYER     BOUN',  &
     &              'DARY   TO BEGIN   IN LAYER    LAYER      TOTAL',   &
     &              '(CM-1)    (KM)     (KM)         (W SR-1 CM-2 / ',  &
     &              'CM-1)                       (W SR-1 CM-2 / CM-1)'
              ENDIF
          ENDIF
      ENDIF

!     INITIALIZE THE SPECTRAL FLUX TABLE:
      IF(IMULT.NE.0 .AND. MODTRN)CALL FLXTBL

!     INITIALIZE LAYER LOOP VARIABLES
      LOOP0=.TRUE.
      CALL LOOP(LOOP0,IV,IVX,IDV5,IKMX,SUMTMS,SUMMS,TRANSM,             &
     &  IPH,SUMSSS,IVTEST,UNIF,TRACE,TRANSX,RADCUM,S0,                  &
     &  GROUND,TSNOBS,TSNREF,FDNTRT,FDNSRT,KNTRVL,NUM_LRES,SOLAR_FLUX,  &
     &  DIFFUSE_FLUX,LNFLRT,FLRT,NUM_LEVS,                              &
     &  FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
      LOOP0=.FALSE.

!         DETERMINE LAYER LOOP MAXIMUM
          IF(TRANSM)THEN

!             FOR TRANSMISSION CALCULATIONS, SKIP OVER LAYER LOOP IN TRA
              IKMAX=1
          ELSEIF(IMULT.NE.0 .AND. .NOT.LUSEMS)THEN

!             FOR MULTIPLE SCATTERING SET IKMAX TO IMSMX
              IKMAX=IMSMX
          ELSE

!             IF NOT MULTIPLE SCATTERING, RESET IKMAX TO ORIGINAL VALUE
              IKMAX=IKMX
          ENDIF

!     END INITIALIZATION, BEGIN OF FREQUENCY LOOP

!     "IVX" IS THE FREQUENCY AT WHICH TRANSMITTANCE WILL BE CALCULATED.
!     DURING THE FIRST PASS, "IVX" AND "IV" MUST BE EQUAL.
   30 CONTINUE
          IVX=IVX+IBNDWD
          IVXBIN=IVXBIN+1
          IF(IV.LT.IVX)THEN
              IV=IV+IDV5
              IVTEST=.TRUE.
          ELSE
              IVTEST=.FALSE.
          ENDIF

!         EXTRA-TERRESTRIAL SOLAR IRRADIANCE:
          IF(IEMSCT.GE.2)THEN
              IF(APPREF .AND. ANGSUN.LT.90.)THEN
                  S0=PI/COS(ANGSUN/DEG)
              ELSE
                  S0=SOURCE(SINIT,IVX,ISOURC,IDAY,ANGLEM,SUNFUN)
              ENDIF
          ENDIF

!         PRINT HEADER EVERY 50 (IPRINT) LINES (NEVER IF JMASS):
          IF(ICOUNT.EQ.IPRINT)THEN

!             REINITIALIZE COUNTER AND PRINT HEADER
              ICOUNT=0
              IF(IEMSCT.EQ.0)THEN
                  WRITE(IPR,'(A,/(3A))')'1',                            &
     &              '  FREQ    WAVLEN    TOTAL     H2O     CO2+  ',     &
     &              '   OZONE    TRACE  N2 CONT  H2O CONT MOL SCAT',    &
     &              '  AER-HYD  HNO3    AER-HYD  INTEGRATED',           &
     &              '  1/CM   MICRONS    TRANS    TRANS    TRANS ',     &
     &              '   TRANS    TRANS   TRANS    TRANS    TRANS  ',    &
     &              '   TRANS   TRANS   TAU-ABS  ABSORPTION'
!NAMX             WRITE(IPR,111)'ALL MINOR SPECIES',
!NAMX1                 (CNAMEX(KX),KX=1,11),
!NAMX2                 (CNAMEX(KX),KX=12,MIN(22,NMOLX)),
!NAMX3                 (' ',KX=MIN(22,NMOLX)+1,22),
!NAMX4                 ('TRANS',KX=1,12)
!111              FORMAT(8X,A,11(1X,A),/,25X,11(1X,A),
!NAMX1                 /,20X,A,8(4X,A),2(3X,A),4X,A/)
              ELSEIF(IEMSCT.EQ.1)THEN
                  WRITE(IPR,'(53X,A,3(/3A))')                           &
     &              'RADIANCE(WATTS/CM2-STER-XXX)',                     &
     &              '  FREQ    WAVLEN  DIREC      PATH THERMAL   ',     &
     &              'SCAT PART    SURFACE EMISSION   SURFACE REFLECTED',&
     &              '     TOTAL RADIANCE   INTEGRAL    TOTAL',          &
     &              '  1/CM   MICRONS   EMIS    (CM-1)   (MICRN) ',     &
     &              '   (CM-1)    (CM-1)   (MICRN)    (CM-1)   (MICRN)',&
     &              '    (CM-1)   (MICRN)              TRANS'
              ELSEIF(IEMSCT.EQ.3)THEN
                  WRITE(IPR,'(/A,22X,A,//(2A))')                        &
     &              '1','IRRADIANCE (WATTS/CM2-XXXX)',                  &
     &              '   FREQ   WAVLEN      TRANSMITTED     ',           &
     &              '      SOLAR           INTEGRATED         TOTAL',   &
     &              '  (CM-1) (MICRN)   (CM-1)    (MICRN)  ',           &
     &              ' (CM-1)    (MICRN)   TRANS.    SOLAR     TRANS'
              ELSEIF(IMULT.EQ.0)THEN
                  WRITE(IPR,'(A,27X,2A,/3(/3A))')                       &
     &              '1','RADIANCE (WATTS/CM2-STER-XXX)',                &
     &              '   [No STER in TOA SOLAR IRRADIANCE]',             &
     &              '  FREQ   WAVLEN    PATH_THERMAL    SURFACE_',      &
     &              'EMISSION  SOL_IRR    SINGLE_SCATTER  GROUND_',     &
     &              'REFLECTED   TOTAL_RADIANCE  INTEGRAL    TOTAL',    &
     &              ' (CM-1) (MICRN)   (CM-1)  (MICRN)   (CM-1) ',      &
     &              ' (MICRN)   (CM-1)   (CM-1)  (MICRN)   (CM-1)',     &
     &              '  (MICRN)   (CM-1)  (MICRN)             TRANS'
              ELSEIF(DIS)THEN
                  WRITE(IPR,'(A,45X,A,/4(/3A))')                        &
     &              '1','RADIANCE(WATTS/CM2-STER-XXX)',                 &
     &              '0 FREQ   WAVLEN   PATH THERMAL    SURFACE   ',     &
     &              '   PATH SCATTERED SOLAR   GROUND REFLECTED ',      &
     &              'RADIANCE   TOTAL RADIANCE   INTEGRAL    TOTAL',    &
     &              '                                  EMISSION  ',     &
     &              '   TOTAL RAD      SINGLE        TOTAL      ',      &
     &              '  DIRECT                                TRANS',    &
     &              ' (CM-1) (MICRN)   (CM-1)  (MICRN)   (CM-1)  ',     &
     &              ' (CM-1)  (MICRN)   (CM-1)   (CM-1)  (MICRN)',      &
     &              '   (CM-1)   (CM-1)  (MICRN)'
              ELSE
                  WRITE(IPR,'(A,52X,A,/4(/3A))')                        &
     &              '1','RADIANCE(WATTS/CM2-STER-XXX)',                 &
     &              '0 FREQ   WAVLEN  DIREC          PATH THERMAL     ',&
     &              ' SURFACE   PATH SCAT SOLAR   GROUND REFLECTED',    &
     &              '     TOTAL RADIANCE   INTEGRAL   TOTAL',           &
     &              '                  EMIS       TOTAL      SCATTERED',&
     &              ' EMISSION    TOTAL   SINGLE    TOTAL   DIRECT',    &
     &              '                                 TRANS',           &
     &              ' (CM-1) (MICRN)          (CM-1)  (MICRN)   (CM-1)',&
     &              '   (CM-1)   (CM-1)   (CM-1)   (CM-1)   (CM-1)',    &
     &              '    (CM-1)   (MICRN)'
              ENDIF
          ENDIF

!         DETERMINE LAYER LOOP MAXIMUM
          IF(TRANSM)THEN

!             LAYER LOOP IS SKIPPED FOR TRANSMISSION CALCULATIONS.
              IKMAX=1
          ELSE
              IF(IMULT.NE.0 .AND. .NOT.LUSEMS)THEN

!                 FOR MULTIPLE SCATTERING SET IKMAX TO IMSMX.
                  IKMAX=IMSMX
              ELSE

!                 IF NOT MULTIPLE SCATTERING, RESET IKMAX TO IKMX.
                  IKMAX=IKMX
              ENDIF

!             SURFACE SPECTRAL REFLECTANCE/EMISSION INTEGRALS:
              V=IVX
              IF(IVX.EQ.0)V=.5*IBNDWD
              IF(LSURF)CALL GTSURF(10000./V,LSURF)
          ENDIF

!         INITIALIZE TRANSMISSION ARRAY
          TX(1)=1.
          TX(2)=1.
          TX(3)=1.
          DO 40 K=4,MEXT
   40     TX(K)=0.

!         CALL LAYER LOOP ROUTINE
          CALL LOOP(LOOP0,IV,IVX,IDV5,IKMX,SUMTMS,SUMMS,TRANSM,         &
     &      IPH,SUMSSS,IVTEST,UNIF,TRACE,TRANSX,RADCUM,S0,              &
     &      GROUND,TSNOBS,TSNREF,FDNTRT,FDNSRT,KNTRVL,NUM_LRES,         &
     &      SOLAR_FLUX,DIFFUSE_FLUX,LNFLRT,FLRT,NUM_LEVS,               &
     &      FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
          IF(IEMSCT.EQ.0)THEN
              IF(LNFILT.GT.0)CALL FILTSM(IVXBIN,1.-TX(9))
          ELSEIF(IEMSCT.EQ.3)THEN

!             TRANSMITTED SOLAR IRRADIANCE
              TS0=TX(9)*S0
              IF(LNFILT.GT.0)CALL FILTSM(IVXBIN,TS0)
          ELSE

!             THERMAL BOUNDARY EMISSION:
              IF(TPTEMP.GT.0.)BBG=BBFN(TPTEMP,V)*TX(9)*(1.-TPHDIR)

!             SURFACE REFLECTED THERMAL SCATTERED RADIANCE:
              RFSURF=0.
              IF(IMULT.NE.0 .AND. GROUND)RFSURF=FDNTRT
              IF(IEMSCT.EQ.2)THEN

!                 SINGLE + MULTIPLE SOLAR SCATTERED RADIANCE
                  IF(LDISCL)THEN

!                     SCALE ISAACS BY DISORT-ISAACS RATIO
                      IF(VDISCL(1).LE.V .AND. V.LE.VDISCL(NDISCL))THEN
                          CALL HUNT(VDISCL,NDISCL,V,JLO)
                          FRAC=                                         &
     &                      (V-VDISCL(JLO))/(VDISCL(JLO+1)-VDISCL(JLO))
                          SCALE=RDISCL(1,JLO)                           &
     &                      +FRAC*(RDISCL(1,JLO+1)-RDISCL(1,JLO))
                          SUMMS=SUMMS*SCALE
                          SUMTMS=SUMTMS*SCALE
                          SCALE=RDISCL(2,JLO)                           &
     &                      +FRAC*(RDISCL(2,JLO+1)-RDISCL(2,JLO))
                          FDNSRT=FDNSRT*SCALE
                          FDNTRT=FDNTRT*SCALE
                      ENDIF
                  ENDIF

                  SUMSSR=SUMSSS+SUMMS

!                 SURFACE REFLECTED SOLAR SCATTER RADIANCES
                  RFLSS=0.
                  IF(GROUND)THEN

!                     DIRECT TERM
                      IF(TSNREF.GT.0.)                                  &
     &                  RFLSS=TSNREF*TPBRDF*COS(ANGSUN/DEG)/PI

!                     SOLAR + THERMAL SURFACE REFLECTED RADIANCE
                      RFSURF=RFSURF+RFLSS
                      IF(IMULT.NE.0)RFSURF=RFSURF+FDNSRT
                  ENDIF
              ENDIF
              IF(LNFILT.GT.0)                                           &
     &          CALL FILTSM(IVXBIN,RADCUM+BBG+RFSURF+SUMSSR)
          ENDIF

!         TRANSMITTANCES, IRRADIANCES [W CM-2 / CM-1], AND RADIANCES
!         [W SR-1 CM-2 / CM-1] ARE TEMPORARILY STORED IN "TX" SO THAT
!         THEIR CONVOLUTION OVER THE TRIANGULAR SLIT CAN BE CALCULATED.
          TX(1)=BBG
          TX(2)=UNIF
          TX(3)=TRACE
          TX(8)=RADCUM
          TX(12)=SUMSSR
          TX(13)=SUMSSS
          TX(14)=TSNREF
          TX(15)=RFSURF
          TX(18)=TSNOBS
          TX(19)=RFLSS
          TX(20)=TS0
          TX(21)=S0
          TX(22)=SUMTMS
          TX(23)=1.-TPHDIR
          DO 60 K=1,MEXT
              IP1=NWGT
              DO 50 I=NWGTM1,1,-1
                  SLIT(K,IP1)=SLIT(K,I)+WGT(IP1)*TX(K)
                  IP1=I
   50         CONTINUE
              SLIT(K,1)=WGT(1)*TX(K)
              TX(K)=SLIT(K,NWGT)
   60     CONTINUE

!         CHECK IF VALUES ARE TO BE PRINTED
          IF(IVX.LT.IWRITE)GOTO30
          IWRITE=IWRITE+IDV
          IF(IWRITE.GT.IVXMAX)FACTOR=.5
          ICOUNT=ICOUNT+1
	  ICOUNT2 = ICOUNT2+1 !DRF
          NSPCDT=NSPCDT+1

!         RENORMALIZE IF TRIANGULAR SLIT EXTENDS TO NEGATIVE FREQUENCIES
          IF(IVX.LT.NWGTM1)THEN
              STORE=1.-.5*(NWGTM1-IVX)*(NWGTM1-IVX+1)/WNORM
              DO 70 K=1,MEXT
   70         TX(K)=TX(K)/STORE
          ENDIF
          BBG=TX(1)
          UNIF=TX(2)
          TRACE=TX(3)
          RADCUM=TX(8)
          SUMSSR=TX(12)
          SUMSSS=TX(13)
          TSNREF=TX(14)
          RFSURF=TX(15)
          TSNOBS=TX(18)
          RFLSS=TX(19)
          TS0=TX(20)
          S0=TX(21)
          SUMTMS=TX(22)
          DEMIS=TX(23)
          V=FLOAT(IVX-IVOFF)
          ALAM=10000./(V+.000001)
          SUMA=SUMA+FACTOR*IDV*(1.0-TX(9))

!         ALTX9 IS NOW OUTPUT USING AN F FORMAT WHEN IEMSCT = 1 OR 2,
!         SO THE MAXIMUM IS REDUCED TO 999.999 (STILL ABSURDLY LARGE)
          ALTX9=999.999
          IF(TX(9).GT.0.)ALTX9=-LOG(TX(9))
          GOTO(80,90,90,100),IEMSCT+1

!         TRANSMITTANCE ONLY
   80     CONTINUE
          TX(7)=TX(7)*TX(16)
          IF(.NOT.LJMASS)THEN

!             WRITE TO OUTPUT FILES:
              WRITE(IPR,'(F6.0,F10.3,11F9.4,F12.3)')                    &
     &          V,ALAM,TX(9),TX(17),UNIF,TX(31),TRACE,                  &
     &          TX(4),TX(5),TX(6),TX(7),TX(11),TX(10),SUMA
              WRITE(IPU,'(F8.2,11(1X,F6.4),1X,F7.3,99(1X,F6.4))')V,     &
     &          TX(9),TX(17),UNIF,TX(31),TRACE,TX(4),TX(5),TX(6),TX(7), &
     &          TX(11),TX(10),ALTX9,TX(36),TX(44),TX(46),TX(47),TX(50), &
     &          TX(52),TX(54),TX(64),TX(65),(TX(IEXT),IEXT=MEXT+1,MEXTX)
              WRITE(IPUSCR,'(F8.2,11(2X,F6.4),1X,F7.3,99(2X,F6.4))')V,  &
     &          TX(9),TX(17),UNIF,TX(31),TRACE,TX(4),TX(5),TX(6),TX(7), &
     &          TX(11),TX(10),ALTX9,TX(36),TX(44),TX(46),TX(47),TX(50), &
     &          TX(52),TX(54),TX(64),TX(65),(TX(IEXT),IEXT=MEXT+1,MEXTX)
          ELSEIF(AMOD3D.NE.'T')THEN

!             FILL JMASS COMMON BLOCKS:
              OUTPUTCOUNT = OUTPUTCOUNT + 1
              IF (OUTPUTCOUNT .GT. ARRAYSIZE)THEN
                  ERRORCODE = FATAL
                  WRITE(ERRORBUF,*)                                     &
     &              'MODTRAN:  Parameter ARRAYSIZE exceeded.'
                  RETURN
              ENDIF
              FREQS(OUTPUTCOUNT) = V
              WAVELENGTHS(OUTPUTCOUNT) = ALAM
              TOTTRANS(OUTPUTCOUNT) = TX(9)
              H2OTRANS(OUTPUTCOUNT) = TX(17)
              CO2PTRANS(OUTPUTCOUNT) = UNIF
              OZONETRANS(OUTPUTCOUNT) = TX(31)
              TRACETRANS(OUTPUTCOUNT) = TRACE
              N2CONTTRANS(OUTPUTCOUNT) = TX(4)
              H2OCONTTRANS(OUTPUTCOUNT) = TX(5)
              MOLSCATTRANS(OUTPUTCOUNT) = TX(6)
              TOTAERTRANS(OUTPUTCOUNT) = TX(7)
              HNO3TRANS(OUTPUTCOUNT) = TX(11)
              AERABSORP(OUTPUTCOUNT) = TX(10)
              CO2TRANS(OUTPUTCOUNT) = TX(36)
              COTRANS(OUTPUTCOUNT) = TX(44)
              CH4TRANS(OUTPUTCOUNT) = TX(46)
              N2OTRANS(OUTPUTCOUNT) = TX(47)
              O2TRANS(OUTPUTCOUNT) = TX(50)
              NH3TRANS(OUTPUTCOUNT) = TX(52)
              NOTRANS(OUTPUTCOUNT) = TX(54)
              NO2TRANS(OUTPUTCOUNT) = TX(55)
              SO2TRANS(OUTPUTCOUNT) = TX(56)
              LOGTOTTRANS(OUTPUTCOUNT) = ALTX9
              INTABS(OUTPUTCOUNT) = SUMA
          ENDIF
          GOTO110

!         RADIANCE PATHS
   90     CONTINUE

!         CONVRT IS THE CONVERSION FROM (W SR-1 CM-2 / CM-1)
!         TO (W SR-1 CM-2 / MICRON).
          CONVRT=1.E-4*V**2
          IF(V.EQ.0.)CONVRT=1.E-4*(.25*IBNDWD)**2
          IF(IEMSCT.EQ.1)THEN

!             SUMT IS THE TOTAL SPECTRAL RADIANCE, I.E. SUMT EQUALS
!             THE SUM OF THE DIRECT + MULTIPLY SCATTERED THERMAL PATH
!             RADIANCE (RADCUM), THE SURFACE EMISSION (BBG), AND THE
!             REFLECTED SURFACE TERM (RFSURF).  SUMTMS IS THE MULTIPLE
!             SCATTERING CONTRIBUTION TO RADCUM AND IS ONLY SEPARATED
!             OUT FOR ISAACS (NOT DISORT) CALCULATIONS.  EACH OF THESE
!             TERMS HAS UNITS (W SR-1 CM-2 / CM-1).  RADSUM IS THE
!             SPECTRALLY INTEGRATED TOTAL RADIANCE (W SR-1 CM-2).
              SUMT=RADCUM+BBG+RFSURF
              RADSUM=RADSUM+IDV*FACTOR*SUMT

!             THERMAL RADIANCE DATA:
              IF(.NOT.LJMASS)THEN

!                 WRITE TO OUTPUT FILES:
                  IF(DIS)THEN
                      WRITE(IPR,'(F6.0,F10.3,                           &
     &                  F7.3,1P,2E10.2,10X,7E10.2,0P,F9.5)')            &
     &                  V,ALAM,DEMIS,RADCUM,CONVRT*RADCUM,              &
     &                  BBG,CONVRT*BBG,RFSURF,CONVRT*RFSURF,            &
     &                  SUMT,CONVRT*SUMT,RADSUM,TX(9)
                      WRITE(IPU,'(0P,F8.2,1P,E11.4,2(E11.4,11X),        &
     &                  2(11X,E11.4),18X,0P,F8.3)')                     &
     &                  V,TX(9),RADCUM,BBG,RFSURF,SUMT,ALTX9
                      WRITE(IPUSCR,'(0P,F8.2,1P,E11.4,2(E11.4,11X),     &
     &                  2(11X,E11.4),18X,0P,F8.3)')                     &
     &                  V,TX(9),RADCUM,BBG,RFSURF,SUMT,ALTX9
                  ELSE
                      WRITE(IPR,'(F6.0,F10.3,F7.3,1P,10E10.2,0P,F9.5)') &
     &                  V,ALAM,DEMIS,RADCUM,CONVRT*RADCUM,SUMTMS,       &
     &                  BBG,CONVRT*BBG,RFSURF,CONVRT*RFSURF,            &
     &                  SUMT,CONVRT*SUMT,RADSUM,TX(9)
                      WRITE(IPU,'(0P,F8.2,1P,4E11.4,11X,                &
     &                  2(11X,E11.4),18X,0P,F8.3)')                     &
     &                  V,TX(9),RADCUM,SUMTMS,BBG,RFSURF,SUMT,ALTX9
                      WRITE(IPUSCR,'(0P,F8.2,1P,4E11.4,11X,             &
     &                  2(11X,E11.4),18X,0P,F8.3)')                     &
     &                  V,TX(9),RADCUM,SUMTMS,BBG,RFSURF,SUMT,ALTX9
                  ENDIF
              ELSEIF(AMOD3D.NE.'T')THEN

!                 FILL JMASS COMMON BLOCKS:
                  OUTPUTCOUNT = OUTPUTCOUNT + 1
                  IF (OUTPUTCOUNT .GT. ARRAYSIZE)THEN
                      ERRORCODE = FATAL
                      WRITE(ERRORBUF,*)                                 &
     &                  'MODTRAN:  Parameter ARRAYSIZE exceeded.'
                      RETURN
                  ENDIF
                  FREQS(OUTPUTCOUNT) = V
                  WAVELENGTHS(OUTPUTCOUNT) = ALAM
                  ATMOSRADCM(OUTPUTCOUNT) = RADCUM
                  ATMOSRADMIC(OUTPUTCOUNT) = CONVRT*RADCUM
                  GRRADCM(OUTPUTCOUNT) = RFSURF
                  GRRADMIC(OUTPUTCOUNT) = CONVRT*RFSURF
                  TRADCM(OUTPUTCOUNT) = SUMT
                  TRADMIC(OUTPUTCOUNT) = CONVRT*SUMT
                  INTRAD(OUTPUTCOUNT) = RADSUM
                  TOTTRANS(OUTPUTCOUNT) = TX(9)
                  LOGTOTTRANS(OUTPUTCOUNT) = ALTX9
                  SRADCM(OUTPUTCOUNT) = SUMSSR
                  SRADMIC(OUTPUTCOUNT) = CONVRT*SUMSSR
                  SRADSSCM(OUTPUTCOUNT) = SUMSSS
                  GRRADDIRCM(OUTPUTCOUNT) = RFLSS
                  GRRADDIRMIC(OUTPUTCOUNT) = CONVRT*RFLSS
              ENDIF
          ELSE

!             SUMT IS THE TOTAL SPECTRAL RADIANCE, I.E. SUMT
!             EQUALS THE SUM OF THE DIRECT + MULTIPLY SCATTERED
!             THERMAL PATH RADIANCE (RADCUM), THE SURFACE EMISSION
!             (BBG), THE SINGLE + MULTIPLE SOLAR SCATTERED RADIANCE
!             (SUMSSR) AND THE REFLECTED SURFACE TERM (RFSURF).  EACH
!             OF THESE TERMS HAS UNITS (W SR-1 CM-2 / CM-1).  RADSUM IS
!             THE SPECTRALLY INTEGRATED TOTAL RADIANCE (W SR-1 CM-2).
              SUMT=RADCUM+BBG+SUMSSR+RFSURF
              RADSUM=RADSUM+IDV*FACTOR*SUMT
              RADCUM = ABS(RADCUM)  !DRF

!             SOLAR + THERMAL RADIANCE DATA:
              IF(.NOT.LJMASS)THEN

!                 WRITE TO OUTPUT FILES:
                  IF(IMULT.EQ.0)THEN
                      WRITE(IPR,'(F7.0,F8.3,1P,12E9.2,0P,F9.5)')        &
     &                  V,ALAM,RADCUM,CONVRT*RADCUM,BBG,CONVRT*BBG,     &
     &                  S0,SUMSSR,CONVRT*SUMSSR,RFSURF,CONVRT*RFSURF,   &
     &                  SUMT,CONVRT*SUMT,RADSUM,TX(9)
                      WRITE(IPU,'(0P,F8.2,F11.8,1P,8E11.4,2E9.2,        &
     &                  0P,F8.3)')V,TX(9),RADCUM,SUMTMS,BBG,SUMSSR,     &
     &                  SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,TSNOBS,ALTX9
                      WRITE(IPUSCR,'(0P,F8.2,F11.8,1P,8E11.4,2E9.2,     &
     &                  0P,F8.3)')V,TX(9),RADCUM,SUMTMS,BBG,SUMSSR,     &
     &                  SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,TSNOBS,ALTX9
                  ELSEIF(DIS)THEN
                      WRITE(IPR,'(F7.0,F8.3,1P,11E9.2,E10.2,0P,F8.5)')  &
     &                  V,ALAM,RADCUM,CONVRT*RADCUM,BBG,SUMSSR,         &
     &                  CONVRT*SUMSSR,SUMSSS,RFSURF,CONVRT*RFSURF,      &
     &                  RFLSS,SUMT,CONVRT*SUMT,RADSUM,TX(9)
                      WRITE(IPU,'(0P,F8.2,F11.8,1P,E11.4,11X,6E11.4,    &
     &                  2E9.2,0P,F8.3)')V,TX(9),RADCUM,BBG,SUMSSR,      &
     &                  SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,TSNOBS,ALTX9
                      WRITE(IPUSCR,'(0P,F8.2,F11.8,1P,E11.4,11X,6E11.3, &
     &                  2E9.2,0P,F8.3)')V,TX(9),RADCUM,BBG,SUMSSR,      &
     &                  SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,TSNOBS,ALTX9
                  ELSE
                      WRITE(IPR,                                        &
     &                  '(F7.0,F8.3,F7.3,1P,8E9.2,3E10.3,0P,F8.5)')     &
     &                  V,ALAM,DEMIS,RADCUM,CONVRT*RADCUM,              &
     &                  SUMTMS,BBG,SUMSSR,SUMSSS,RFSURF,RFLSS,          &
     &                  SUMT,CONVRT*SUMT,RADSUM,TX(9)
                      WRITE(IPU,'(0P,F8.2,F11.8,1P,8E11.4,2E9.2,        &
     &                  0P,F8.3)')V,TX(9),RADCUM,SUMTMS,BBG,SUMSSR,     &
     &                  SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,TSNOBS,ALTX9
                      WRITE(IPUSCR,'(0P,F8.2,F11.8,1P,8E11.4,2E9.2,     &
     &                  0P,F8.3)')V,TX(9),RADCUM,SUMTMS,BBG,SUMSSR,     &
     &                  SUMSSS,RFSURF,RFLSS,SUMT,TSNREF,TSNOBS,ALTX9
                  ENDIF
              ELSEIF(AMOD3D.NE.'T')THEN

!                 FILL JMASS COMMON BLOCKS:
                  OUTPUTCOUNT = OUTPUTCOUNT + 1
                  IF (OUTPUTCOUNT .GT. ARRAYSIZE)THEN
                      ERRORCODE = FATAL
                      WRITE(ERRORBUF,*)                                 &
     &                  'MODTRAN:  Parameter ARRAYSIZE exceeded.'
                      RETURN
                  ENDIF
                  FREQS(OUTPUTCOUNT) = V
                  WAVELENGTHS(OUTPUTCOUNT) = ALAM
                  ATMOSRADCM(OUTPUTCOUNT) = RADCUM
                  ATMOSRADMIC(OUTPUTCOUNT) = CONVRT*RADCUM
                  GRRADCM(OUTPUTCOUNT) = RFSURF
                  GRRADMIC(OUTPUTCOUNT) = CONVRT*RFSURF
                  TRADCM(OUTPUTCOUNT) = SUMT
                  TRADMIC(OUTPUTCOUNT) = CONVRT*SUMT
                  INTRAD(OUTPUTCOUNT) = RADSUM
                  TOTTRANS(OUTPUTCOUNT) = TX(9)
                  LOGTOTTRANS(OUTPUTCOUNT) = ALTX9
                  SRADCM(OUTPUTCOUNT) = SUMSSR
                  SRADMIC(OUTPUTCOUNT) = CONVRT*SUMSSR
                  SRADSSCM(OUTPUTCOUNT) = SUMSSS
                  GRRADDIRMIC(OUTPUTCOUNT) = RFLSS
                  SOLARIRRADCM(OUTPUTCOUNT) = S0
                  SOLARIRRADMIC(OUTPUTCOUNT) = CONVRT*S0
              ENDIF
          ENDIF
          GOTO110

!         DIRECTLY TRANSMITTED SOLAR IRRADIANCE [WATTS/(CM2 MICROMETER)]
  100     CONTINUE

!         CONVRT IS THE CONVERSION FROM W/CM2-(CM-1) TO W/CM2-UM
          CONVRT=1.E-4*V**2
          IF(V.EQ.0.)CONVRT=1.E-4*(.25*IBNDWD)**2

!         RADSUM IS THE INTEGRATED TRANSMITTED SOLAR IRRADIANCE AND SSOL
!         IS THE INTEGRATED EXTRA-TERRESTRIAL SOLAR IRRADIANCE (W/CM2).
          RADSUM=RADSUM+TS0*IDV*FACTOR
          SSOL=SSOL+S0*IDV*FACTOR

!         SOLAR IRRADIANCE DATA:
          IF(.NOT.LJMASS)THEN

!             WRITE TO OUTPUT FILES:
              WRITE(IPR,'(F8.0,F8.3,1P,6E10.2,0P,F9.4)')                &
     &          V,ALAM,TS0,CONVRT*TS0,S0,CONVRT*S0,RADSUM,SSOL,TX(9)
              WRITE(IPU,'(F8.2,F8.4,1P,2E9.2,T96,E10.3)')               &
     &          V,TX(9),TS0,S0,ALTX9
              WRITE(IPUSCR,'(F8.2,F8.4,1P,2E9.2,T96,E10.3)')            &
     &          V,TX(9),TS0,S0,ALTX9
              SUMT=TS0
          ELSEIF(AMOD3D.NE.'T')THEN

!             FILL JMASS COMMON BLOCKS:
              OUTPUTCOUNT = OUTPUTCOUNT + 1
              IF (OUTPUTCOUNT .GT. ARRAYSIZE)THEN
                  ERRORCODE = FATAL
                  WRITE(ERRORBUF,*)                                     &
     &              'MODTRAN:  Parameter ARRAYSIZE exceeded.'
                  RETURN
              ENDIF
              FREQS(OUTPUTCOUNT) = V
              WAVELENGTHS(OUTPUTCOUNT) = ALAM
              TRANSIRRADCM(OUTPUTCOUNT) = TS0
              TRANSIRRADMIC(OUTPUTCOUNT) = CONVRT*TS0
              SOLARIRRADCM(OUTPUTCOUNT) = S0
              SOLARIRRADMIC(OUTPUTCOUNT) = CONVRT*S0
!????????     INTTIRRAD(OUTPUTCOUNT) = STSOL
              INTRAD(OUTPUTCOUNT) = RADSUM
              INTSIRRAD(OUTPUTCOUNT) = SSOL
              TOTTRANS(OUTPUTCOUNT) = TX(9)
              LOGTOTTRANS(OUTPUTCOUNT) = ALTX9
          ENDIF
  110     CONTINUE

!         WRITE OUT PLOT.DAT FILE:
          IF(YFLAG.EQ.'T')THEN
              IF(TX(9).GT.YPLTMX)YPLTMX=TX(9)
              IF(XFLAG.EQ.'N')THEN
                  WRITE(IPLOT,'(0P,F15.3,F15.8)')1000*ALAM,TX(9)
!		  IPLOT(ICOUNT2,1) = 1000.*ALAM
!		  IPLOT(ICOUNT2,2) = TX(9)
              ELSEIF(XFLAG.EQ.'M')THEN
                  WRITE(IPLOT,'(0P,F15.6,F15.8)')ALAM,TX(9)
!		  IPLOT(ICOUNT2,1) = ALAM
!		  IPLOT(ICOUNT2,2) = TX(9)
              ELSE
                  WRITE(IPLOT,'(0P,F15.0,F15.8)')V,TX(9)
!		  IPLOT(ICOUNT2,1) = V
!		  IPLOT(ICOUNT2,2) = TX(9)
              ENDIF
          ELSEIF(YFLAG.EQ.'R')THEN
              IF(XFLAG.EQ.'N')THEN
                  CONVRT=1000.*CONVRT
                  IF(CONVRT*SUMT.GT.YPLTMX)YPLTMX=CONVRT*SUMT
                  WRITE(IPLOT,'(0P,F15.3,1P,E15.5)')                    &
     &              1000.*ALAM,CONVRT*SUMT
                  !WRITE(*,*) 'ICOUNT2 = ',ICOUNT2
                  !WRITE(*,*) 'WVL = ',1000.*ALAM
                  !WRITE(*,*) 'WV_HRES(1) = ',WV_HRES
                  WV_HRES(ICOUNT2) = 1000.*ALAM
                  RADIANCE_HRES(ICOUNT2) = REAL(CONVRT*SUMT)
                  !WRITE(*,*) 'RAD = ',RADIANCE_HRES(ICOUNT2)
!                  RADIANCE_HRES(1) = 1000.*ALAM
!DRF		  IPLOT(ICOUNT2,1) = 1000.*ALAM  
!DRF		  IPLOT(ICOUNT2,2) = CONVRT*SUMT
              ELSEIF(XFLAG.EQ.'M')THEN
                  IF(CONVRT*SUMT.GT.YPLTMX)YPLTMX=CONVRT*SUMT
                  WRITE(IPLOT,'(0P,F15.6,1P,E15.5)')ALAM,CONVRT*SUMT
!		  IPLOT(ICOUNT2,1) = ALAM
!		  IPLOT(ICOUNT2,2) = CONVRT*SUMT
              ELSE
                  IF(SUMT.GT.YPLTMX)YPLTMX=SUMT
                  WRITE(IPLOT,'(0P,F15.0,1P,2E15.5)')V,SUMT
!		  IPLOT(ICOUNT2,1) = V
!		  IPLOT(ICOUNT2,2) = SUMT
              ENDIF
          ENDIF
          IF(IEMSCT.NE.0)THEN
              IF(SUMT.GE.RADMAX)THEN
                  IVRMAX=INT(V+.5)
                  RADMAX=SUMT
              ENDIF
              IF(SUMT.LE.RADMIN)THEN
                  IVRMIN=INT(V+.5)
                  RADMIN=SUMT
              ENDIF
          ENDIF
          FACTOR=1.
      IF(IWRITE.LE.IVXMAX)GOTO30

!     END OF FREQUENCY LOOP; PRINT PLOT.DAT FILE DELIMITER.
      IF(YFLAG.EQ.'T')THEN
          IF (DLIMIT.EQ.'        ')THEN
                  WRITE(IPR,*) 'DLIMIT = blank'
          ELSEIF(XFLAG.EQ.'N')THEN
              WRITE(IPR,'(2A,1P,E15.5,A)')DLIMIT,                       &
     &          '   THE MAXIMUM TRANSMITTANCE IS',                      &
     &          YPLTMX,' (WAVELENGTHS IN NANOMETERS)'
          ELSEIF(XFLAG.EQ.'M')THEN
              WRITE(IPR,'(2A,1P,E15.5,A)')DLIMIT,                       &
     &          '   THE MAXIMUM TRANSMITTANCE IS',                      &
     &          YPLTMX,' (WAVELENGTHS IN MICRONS)'
          ELSE
              WRITE(IPR,'(2A,1P,E15.5,A)')DLIMIT,                       &
     &          '   THE MAXIMUM TRANSMITTANCE IS',                      &
     &          YPLTMX,' (FREQUENCIES IN CM-1)'
          ENDIF
      ELSEIF(YFLAG.EQ.'R')THEN
          IF(IEMSCT.EQ.3)THEN
              IF (DLIMIT.EQ.'        ')THEN
                  WRITE(IPR,*) 'DLIMIT = blank'
              ELSEIF(XFLAG.EQ.'N')THEN
                  WRITE(IPR,'(2A,1P,E15.5,A)')DLIMIT,                   &
     &              '   THE MAXIMUM TRANSMITTED SOLAR IRRADIANCE IS',   &
     &              YPLTMX,' MICRO-WATTS CM-2 / NANOMETER'
              ELSEIF(XFLAG.EQ.'M')THEN
                  WRITE(IPR,'(2A,1P,E15.5,A)')DLIMIT,                   &
     &              '   THE MAXIMUM TRANSMITTED SOLAR IRRADIANCE IS',   &
     &              YPLTMX,' W CM-2 / MICRON'
              ELSE
                  WRITE(IPR,'(2A,1P,E15.5,A)')DLIMIT,                   &
     &              '   THE MAXIMUM TRANSMITTED SOLAR IRRADIANCE IS',   &
     &              YPLTMX,' W CM-2 / CM-1'
              ENDIF
          ELSE
              IF (DLIMIT.EQ.'        ')THEN
                  WRITE(IPR,*) 'DLIMIT = blank'
              ELSEIF(XFLAG.EQ.'N')THEN
                  !DRF, usually write to IPLOT instead
                  WRITE(IPR,'(2A,1P,E15.5,A)')DLIMIT,                   &
     &              '   THE MAXIMUM SPECTRAL RADIANCE IS',              &
     &              YPLTMX,' MICRO-WATTS SR-1 CM-2 / NANOMETER'
              ELSEIF(XFLAG.EQ.'M')THEN
                  WRITE(IPR,'(2A,1P,E15.5,A)')DLIMIT,                   &
     &              '   THE MAXIMUM SPECTRAL RADIANCE IS',              &
     &              YPLTMX,' W SR-1 CM-2 / MICRON'
              ELSE
                  WRITE(IPR,'(2A,1P,E15.5,A)')DLIMIT,                   &
     &              '   THE MAXIMUM SPECTRAL RADIANCE IS',              &
     &              YPLTMX,' W SR-1 CM-2 / CM-1'
              ENDIF
          ENDIF
      ENDIF
      IF(.NOT.LJMASS .AND. IMULT.NE.0)THEN

!         INTEGRATED COOLING RATES.
          WRITE(IPR,'(/A)')'1 MULTIPLE SCATTERING CALCULATION RESULTS:'
          IF(ML.GE.4)CALL WTCOOL(IBNDWD,IEMSCT)

!         WRITE OUT IN-BAND LAYER BOUNDARY FLUX VALUES:
          !IF(MODTRN)CALL FLXSUM
          IF(MODTRN)CALL FLXSUM(NUM_LEVS,BB_UPDIFFUSE,BB_DNDIFFUSE,     &
     &                          BB_DNDIRECT) 
          !WRITE(*,*) 'BLAH FLXSUM'
      ENDIF

      IF (LDISCL) THEN
!        Reclaim saved DISORT values in case there is a repeat run.
         DISAZM=DISAZMSV
         DISALB=DISALBSV
         NAZ=NAZSV
         N2GAUS=N2GAUSSV
         NSTR=NSTRSV
         DIS=DISSV
      ENDIF

!     BAND PASS DATA:
      IF(LJMASS)THEN

!         FILL JMASS COMMON BLOCKS:
          MINIMUMRAD = RADMIN
          MAXIMUMRAD = RADMAX
          MINRADVALUE = IVRMIN
          MAXRADVALUE = IVRMAX
          BOUNDEMIS = 1.-AASALB
      ELSE

!         WRITE TO OUTPUT FILES:
          WRITE(IPR,'(/A,2(I6,A),F12.4,A,/A,F7.4)')                     &
     &      ' INTEGRATED ABSORPTION FROM',IV1,' TO',INT(V+.5),' CM-1 =',&
     &      SUMA,' CM-1',' AVERAGE TRANSMITTANCE =',1.-SUMA/(V-IV1)
          IF(IEMSCT.EQ.3)THEN
              WRITE(IPR,'(/A,1P,E15.6,2(A,I6),A,/(A,E15.6,A,I6,A))')    &
     &          ' INTEGRATED IRRADIANCE =',RADSUM,                      &
     &          ' WATTS CM-2 (FROM',IV1,' TO',IV2,' CM-1)',             &
     &          ' MINIMUM IRRADIANCE    =',RADMIN,                      &
     &          ' WATTS CM-2 AT',IVRMIN,' CM-1',                        &
     &          ' MAXIMUM IRRADIANCE    =',RADMAX,                      &
     &          ' WATTS CM-2 AT',IVRMAX,' CM-1'
          ELSEIF(IEMSCT.NE.0)THEN
              WRITE(IPR,'(/A,1P,E15.6,A,2(A,I6),A/(A,E15.6,A,I6,A))')   &
     &          ' INTEGRATED TOTAL RADIANCE =',RADSUM,                  &
     &          ' WATTS CM-2 STER-1',' (FROM',IV1,' TO',IV2,' CM-1 )',  &
     &          ' MINIMUM SPECTRAL RADIANCE =',RADMIN,                  &
     &          ' WATTS CM-2 STER-1 / CM-1  AT',IVRMIN,' CM-1',         &
     &          ' MAXIMUM SPECTRAL RADIANCE =',RADMAX,                  &
     &          ' WATTS CM-2 STER-1 / CM-1  AT',IVRMAX,' CM-1'
              WRITE(IPR,'(/(A,F11.3))')                                 &
     &          ' TARGET-PIXEL (H2) SURFACE TEMPERATURE [K] =',TPTEMP,  &
     &          ' AREA-AVERAGED GROUND TEMPERATURE [K]      =',AATEMP
              IF(NWVSRF(1).EQ.1)WRITE(IPR,'(A,F11.3)')                  &
     &          ' TARGET-PIXEL (H2) DIRECTIONAL EMISSIVITY  =',1.-TPHDIR
              IF(NWVSRF(2).EQ.1)WRITE(IPR,'(A,F11.3)')                  &
     &          ' AREA-AVERAGED GROUND EMISSIVITY           =',1.-AASALB
          ENDIF
          IF(LNFILT.GT.0)CALL WRTFLT(IEMSCT,IVXMIN,IVX)
      ENDIF
      RETURN
      END

      SUBROUTINE TNOWRT(IPH,ISOURC,IDAY,ANGLEM,                         &
     &  GROUND,IMSMX,KNTRVL,SUNFUN,APPREF)

!     ROUTINE TRANS CALCULATES TRANSMITTANCE AND RADIANCE VALUES
!     BETWEEN IV1 AND IV2 FOR A GIVEN ATMOSPHERIC SLANT PATH.

!     ARGUMENTS:
!       SUNFUN   TOP-OF-ATMOSPHERE SOLAR IRRADIANCE
!                FUNCTION [W CM-2 / CM-1].
      LOGICAL GROUND,APPREF
      INTEGER IPH,ISOURC,IDAY,IMSMX,KNTRVL
      REAL ANGLEM,SUNFUN
      EXTERNAL SUNFUN

!     PARAMETERS:
      INCLUDE 'PARAMS.h'
      INTEGER MDISCL
      PARAMETER(MDISCL=12)

!     COMMONS:
      INCLUDE 'BMHEAD.h'
      INCLUDE 'SOLS.h'
      INCLUDE 'BASE.h'

!     /CARD1/
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,MODTRN

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

!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT

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

!     /DISCAL/
      REAL VDISCL(MDISCL),RDISCL(2,MDISCL)
      INTEGER NDISCL,IDISCL
      COMMON/DISCAL/NDISCL,VDISCL,RDISCL
      INTEGER JLO,JHI,ISWITCH

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:
      REAL SOURCE

!     LOCAL VARIABLES:
!       TRANSM  FLAG [TRUE FOR TRANSMITTANCE ONLY CALCULATIONS].
!       LSURF   SPECTRAL FLAG [TRUE IF REFLECTANCES ARE CHANGING].
!       SINIT   FLAG, TRUE FOR INITIAL CALL TO SOURCE FROM TRANS.
!       IWIDTH  FWHM IN BAND WIDTHS (=IFWHM/IBNDWD).
!       WNORM   SLIT FUNCTION NORMALIZATION FACTOR (= 1/IWIDTH**2).
!       SUMTMS   SCATTERED THERMAL EMISSION [W SR-1 CM-2 / CM-1].
!       TSNREF   SENSOR-FINAL_ALTITUDE-SUN TRANSMITTANCE CONVOLVED WITH
!                THE TOP-OF-ATMOSPHERE SOLAR IRRADIANCE [W CM-2 / CM-1].
!       TSNOBS   TRANSMITTED SOLAR IRRADIANCE AT SENSOR [W CM-2 / CM-1].
      LOGICAL TRANSM,LSURF,SINIT,DAZMSV,IVTEST,LOOP0
      INTEGER IV,I,IKMX,IDV5,K,NSTRSV
      REAL V,SUMSSS,S0,SUMTMS,UNIF,TRACE,TRANSX,RADCUM,TSNOBS,TSNREF,   &
     &  SUMMSD,FDNSRD,FDNTRD,SUMMS,FDNSRT,FDNTRT

!     DATA:
      REAL VSCALE(MDISCL)
      SAVE VSCALE
      DATA VSCALE/4665.0,  5955.0,  8040.0,  9840.0, 11625.0, 15585.0,  &
     &  18000.0, 22020.0, 29400.0, 30975.0, 32550.0, 49995.0/

!     BEFORE SCALING, LDISCL IS SET TO .FALSE.
      LDISCL=.FALSE.

!     ALL VALUES OF VSCALE ARE MULTIPLES OF 15.
      CALL HUNT(VSCALE,MDISCL,IV1+0.01,JLO)
      CALL HUNT(VSCALE,MDISCL,FLOAT(IV2),JHI)
      IF(JLO.EQ.0)JLO=1
      IF(JHI.EQ.0)JHI=1
      NDISCL=JHI-JLO+1
      DO I=JLO,JHI
         VDISCL(I-JLO+1)=VSCALE(I)
      ENDDO
      IF(JHI+1.LE.MDISCL)THEN
         NDISCL=NDISCL+1
         VDISCL(NDISCL)=VSCALE(JHI+1)
      ENDIF

!     IDV5 INCREMENT FOR RETRIEVING 5 CM-1 RESOLUTION DATA [CM-1].
      IDV5=5*INT((IBNDWD+4.999)/5)

      NSTRSV=NSTR
      DAZMSV=DISAZM

      DO IDISCL=1, NDISCL
         DIS=.TRUE.
         DISAZM=DAZMSV
         NSTR=NSTRSV
         NAZ=ABS(NSTR)-1
         N2GAUS=ABS(NSTR)/2
         DO ISWITCH=1, 2

!           STORE THE NUMBER OF PATH LAYERS IN IKMX
            IKMX=IKMAX

!           INITIALIZE TRANSMITTANCE/RADIANCE/IRRADIANCE TERMS:
            UNIF=0.
            TRACE=0.
            SUMSSS=0.
            S0=0.
            TSNOBS=0.
            TSNREF=0.
            SUMTMS=0.

            IVX=NINT(VDISCL(IDISCL))
            IV=IVX
            IVTEST=.TRUE.

!           CALL ROUTINE "BMDATA" TO INITIALIZE BAND MODEL PARAMETERS.
            IF(MODTRN .AND. IVX.LE.MXFREQ)CALL BMDATA(IV,IKMX,IMSMX)

!           INITIALIZE TRANSMITTANCE AND SURFACE REFLECTANCE FLAGS:
            TRANSM=.FALSE.
            LSURF=.TRUE.

!           INITIALIZE LAYER LOOP VARIABLES

            LOOP0=.TRUE.
      CALL LOOP(LOOP0,IV,IVX,IDV5,IKMX,SUMTMS,SUMMS,TRANSM,             &
     &  IPH,SUMSSS,IVTEST,UNIF,TRACE,TRANSX,RADCUM,S0,                  &
     &  GROUND,TSNOBS,TSNREF,FDNTRT,FDNSRT,KNTRVL,NUM_LRES,SOLAR_FLUX,  &
     &  DIFFUSE_FLUX,LNFLRT,FLRT,NUM_LEVS,                              &
     &  FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
            LOOP0=.FALSE.

!           END INITIALIZATION, BEGIN OF FREQUENCY LOOP
            SINIT=.TRUE.
            IF(IEMSCT.GE.2)THEN
               IF(APPREF .AND. ANGSUN.LT.90.)THEN
                  S0=PI/COS(ANGSUN/DEG)
               ELSE
                  S0=SOURCE(SINIT,IVX,ISOURC,IDAY,ANGLEM,SUNFUN)
               ENDIF
            ENDIF

!           FOR MULTIPLE SCATTERING SET IKMAX TO IMSMX.
            IKMAX=IMSMX

!           SURFACE SPECTRAL REFLECTANCE/EMISSION INTEGRALS:
            V=IVX
            IF(IVX.EQ.0)V=.5*IBNDWD
            IF(LSURF)CALL GTSURF(10000./V,LSURF)

!           INITIALIZE TRANSMISSION ARRAY
            TX(1)=1.
            TX(2)=1.
            TX(3)=1.
            DO K=4,MEXT
                TX(K)=0.
            ENDDO

!           CALL LAYER LOOP ROUTINE
            IF(ISWITCH.EQ.1)THEN

!              DISORT:
          CALL LOOP(LOOP0,IV,IVX,IDV5,IKMX,SUMTMS,SUMMSD,TRANSM,        &
     &      IPH,SUMSSS,IVTEST,UNIF,TRACE,TRANSX,RADCUM,S0,              &
     &      GROUND,TSNOBS,TSNREF,FDNTRD,FDNSRD,KNTRVL,NUM_LRES,         &
     &      SOLAR_FLUX,DIFFUSE_FLUX,LNFLRT,FLRT,NUM_LEVS,               &
     &      FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)

!              SET UP FOR ISAACS (ISWITCH=2):
               DIS=.FALSE.
               DISAZM=.FALSE.
               NSTR=1
               NAZ=-1
               N2GAUS=0
            ELSE

!              ISAACS:
          CALL LOOP(LOOP0,IV,IVX,IDV5,IKMX,SUMTMS,SUMMS,TRANSM,         &
     &      IPH,SUMSSS,IVTEST,UNIF,TRACE,TRANSX,RADCUM,S0,              &
     &      GROUND,TSNOBS,TSNREF,FDNTRT,FDNSRT,KNTRVL,NUM_LRES,         &
     &      SOLAR_FLUX,DIFFUSE_FLUX,LNFLRT,FLRT,NUM_LEVS,               &
     &      FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
               SUMMS=SUMMS+SUMTMS
            ENDIF
         ENDDO

!        NEED THESE DISORT AND ISAAC RATIOS:
         IF(SUMMS.GT.1.D-30*DBLE(SUMMSD))THEN
             RDISCL(1,IDISCL)=SUMMSD/SUMMS
         ELSEIF(IDISCL.GT.1)THEN
             RDISCL(1,IDISCL)=RDISCL(1,IDISCL-1)
         ELSE
             RDISCL(1,IDISCL)=1.
         ENDIF
         IF(FDNSRT+FDNTRT.GT.1.D-30*DBLE(FDNSRD+FDNTRD))THEN
             RDISCL(2,IDISCL)=(FDNSRD+FDNTRD)/(FDNSRT+FDNTRT)
         ELSEIF(IDISCL.GT.1)THEN
             RDISCL(2,IDISCL)=RDISCL(2,IDISCL-1)
         ELSE
             RDISCL(2,IDISCL)=1.
         ENDIF
      ENDDO

!     RESET LDISCL TO .TRUE.
      LDISCL=.TRUE.
      RETURN
      END
