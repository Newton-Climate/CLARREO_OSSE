      SUBROUTINE FLXTBL

!     FLXTBL PERFORMS INITIALIZATIONS AND
!     WRITES HEADER FOR THE FLUX TABLE.

!     LIST PARAMETERS:
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /WTFLX/
!       NFLUX   SPECTRAL BIN COUNTER FOR FLUX TABLE.
!       UPDIFF  BOUNDARY UPWARD DIFFUSE SPECTRAL FLUX [W CM-2 / CM-1].
!       DNDIFF  BOUNDARY DOWNWARD DIFFUSE SPECTRAL FLUX [W CM-2 / CM-1].
!       DNDRCT  BOUNDARY DIRECT SOLAR SPECTRAL FLUX [W CM-2 / CM-1].
!       NTERMS  NUMBER OF TERMS IN FLUX SPECTRAL SUM.
!       SMUPDF  LAYER BOUNDARY UPWARD DIFFUSE IN-BAND FLUX [W CM-2].
!       SMDNDF  LAYER BOUNDARY DOWNWARD DIFFUSE IN-BAND FLUX [W CM-2].
!       SMDNDR  LAYER BOUNDARY DIRECT SOLAR IN-BAND FLUX [W CM-2].
      INTEGER NFLUX,NTERMS
      REAL UPDIFF,DNDIFF,DNDRCT,SMUPDF,SMDNDF,SMDNDR
      COMMON/WTFLX/NFLUX,UPDIFF(0:NBINS,1:LAYDIM),                      &
     &  DNDIFF(0:NBINS,1:LAYDIM),DNDRCT(0:NBINS,1:LAYDIM),              &
     &  NTERMS,SMUPDF(LAYDIM),SMDNDF(LAYDIM),SMDNDR(LAYDIM)

!     /WTFLXC/
!       FRMT    FORMAT USED IN FLUX TABLE.
      CHARACTER*50 FRMT
      COMMON/WTFLXC/FRMT

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT

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

!     COMMON/SCAN/
!       V1       LOWER BOUND ON SPECTRAL RANGE [CHUNIT DEFINES UNIT].
!       V2       UPPER BOUND ON SPECTRAL RANGE [CHUNIT DEFINES UNIT].
!       DV       SPECTRAL STEP SIZE FOR OUTPUT [CHUNIT DEFINES UNIT].
!       FWHM     FULL-WIDTH-AT-HALF-MAXIMUM [CHUNIT DEFINES UNIT].
!       FWHMSQ   TRIANGULAR SLIT NORMALIZATION FACTOR
!                (EQUALS FWHM SQUARED) [CHUNIT DEFINES UNIT].
!       VOUT     CURRENT SPECTRAL OUTPUT [CHUNIT DEFINES UNIT].
      REAL V1,V2,DV,FWHM,FWHMSQ,VOUT
      COMMON/SCAN/V1,V2,DV,FWHM,FWHMSQ,VOUT

!     COMMON/CSCAN/
!       CHUNIT   UNIT FLAG ('W'=WAVENUMBERS;'M'=MICRONS;'N'=NANOMETERS).
!       RELABS   SPECTRAL RESOLUTION FLAG('A'=ABSOLUTE;'R'=RELATIVE[%]).
!       LNFEED   LINE FEED FLAG FOR .FLX FILE ('T' FOR 80 CHARACTER
!                  LINES, 'F' FOR LONG LINES, ' ' FOR NO .FLX FILE).
      CHARACTER*1 CHUNIT,RELABS,LNFEED
      COMMON/CSCAN/CHUNIT,RELABS,LNFEED

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
!       IK       LAYER INDEX.
!       IBIN     SPECTRAL INDEX FOR SLIT FUNCTION.
!       UNIT     UNIT FLAG ('W'=WAVENUMBERS;'M'=MICRONS;'N'=NANOMETERS).
      INTEGER IK,IBIN
      CHARACTER*1 UNIT

!     BRANCH BASED ON SLIT FUNCTION TYPE:
      !WRITE(*,*) 'into flxtbl'
      NTERMS=0
      IF(FLAGS(1:4).EQ.'    ')THEN

!         INITIALIZE SPECTRAL BIN COUNTER:
          NFLUX=1

!         INITIALIZE DOWNWARD DIRECT SOLAR FLUX ARRAY:
          DO 20 IK=1,ML
              DO 10 IBIN=1,NWGT
                  DNDRCT(IBIN,IK)=0.
   10         CONTINUE
              SMUPDF(IK)=0.
              SMDNDF(IK)=0.
              SMDNDR(IK)=0.
   20     CONTINUE

!         RETURN IF FLUX FILE IS NOT TO BE WRITTEN:
          IF(LNFEED.EQ.' ')RETURN

!         WRITE TABLE HEADER:
          IF(XFLAG.EQ.'M')THEN
              IF(.NOT.LJMASS) WRITE(IFLUX,'(///(A))')                   &
     &          ' SPECTRAL VERTICAL FLUX TABLE (W CM-2 / MICRON)',      &
     &          ' ----------------------------------------------'
              UNIT='M'
          ELSEIF(XFLAG.EQ.'N')THEN
              IF(.NOT.LJMASS) WRITE(IFLUX,'(///(A))')                   &
     &          ' SPECTRAL VERTICAL FLUX TABLE (W CM-2 / NM)',          &
     &          ' ------------------------------------------'
              UNIT='N'
          ELSE
              IF(.NOT.LJMASS) WRITE(IFLUX,'(///(A))')                   &
     &          ' SPECTRAL VERTICAL FLUX TABLE (W CM-2 / CM-1)',        &
     &          ' --------------------------------------------'
              UNIT='W'
          ENDIF
          IF(.NOT.LJMASS)THEN
              WRITE(IFLUX,'(/2(I3,A))')                                 &
     &          MLFLX,' LEVELS (',3*MLFLX+1,' TABLE COLUMNS)'
              IF(IFWHM.EQ.1)THEN
                  WRITE(IFLUX,'(/A17,A,/)')' RECTANGULAR SLIT',         &
     &              ' FULL-WIDTH-AT-HALF-MAXIMUM:      1.00 CM-1.'
              ELSE
                  WRITE(IFLUX,'(/A17,A,F10.2,A,/)')' TRIANGULAR  SLIT', &
     &              ' FULL-WIDTH-AT-HALF-MAXIMUM:',FLOAT(IFWHM),' CM-1.'
              ENDIF
          ENDIF
      ELSE

!         CURRENT FLUX DATA IS STORED IN FIRST ELEMENT:
          NFLUX=0

!         INITIALIZE ENTIRE FLUX ARRAYS:
          DO 40 IK=1,ML
              DO 30 IBIN=0,NBINS
                  UPDIFF(IBIN,IK)=0.
                  DNDIFF(IBIN,IK)=0.
                  DNDRCT(IBIN,IK)=0.
   30         CONTINUE
              SMUPDF(IK)=0.
              SMDNDF(IK)=0.
              SMDNDR(IK)=0.
   40     CONTINUE

!         RETURN IF FLUX FILE IS NOT TO BE WRITTEN:
          IF(LNFEED.EQ.' ')RETURN

!         WRITE TABLE HEADER:
          IF(CHUNIT.EQ.'M')THEN
              IF(.NOT.LJMASS) WRITE(IFLUX,'(///(A))')                   &
     &          ' SPECTRAL VERTICAL FLUX TABLE (W CM-2 / MICRON)',      &
     &          ' ----------------------------------------------'
              UNIT='M'
              VOUT=V2
          ELSEIF(CHUNIT.EQ.'N')THEN
              IF(.NOT.LJMASS) WRITE(IFLUX,'(///(A))')                   &
     &          ' SPECTRAL VERTICAL FLUX TABLE (W CM-2 / NM)',          &
     &          ' ------------------------------------------'
              UNIT='N'
              VOUT=V2
          ELSE
              IF(.NOT.LJMASS) WRITE(IFLUX,'(///(A))')                   &
     &          ' SPECTRAL VERTICAL FLUX TABLE (W CM-2 / CM-1)',        &
     &          ' --------------------------------------------'
              UNIT='W'
              VOUT=V1
          ENDIF
          IF(.NOT.LJMASS) WRITE(IFLUX,'(/2(I3,A))')                     &
     &      MLFLX,' LEVELS (',3*MLFLX+1,' TABLE COLUMNS)'
          IF(RELABS.EQ.'R')THEN
              IF(.NOT.LJMASS)                                           &
     &          WRITE(IFLUX,'(/A17,A,F10.5,A,/)')' TRIANGULAR  SLIT',   &
     &          ' FULL-WIDTH-AT-HALF-MAXIMUM:',FWHM,' PERCENT.'
          ELSE
              FWHMSQ=FWHM**2
              IF(.NOT.LJMASS)THEN
                  IF(CHUNIT.EQ.'M')THEN
                      WRITE(IFLUX,'(/A17,A,F10.6,A,/)')                 &
     &                  ' TRIANGULAR  SLIT',                            &
     &                  ' FULL-WIDTH-AT-HALF-MAXIMUM:',FWHM,' MICRONS.'
                  ELSEIF(CHUNIT.EQ.'N')THEN
                      WRITE(IFLUX,'(/A17,A,F10.3,A,/)')                 &
     &                  ' TRIANGULAR  SLIT',                            &
     &                  ' FULL-WIDTH-AT-HALF-MAXIMUM:',FWHM,            &
     &                  ' NANOMETERS.'
                  ELSE
                      WRITE(IFLUX,'(/A17,A,F10.2,A,/)')                 &
     &                  ' TRIANGULAR  SLIT',                            &
     &                  ' FULL-WIDTH-AT-HALF-MAXIMUM:',FWHM,' CM-1.'
                  ENDIF
              ENDIF
          ENDIF
      ENDIF
      IF(.NOT.LJMASS)THEN
        IF(LNFEED.EQ.'T')THEN
          WRITE(IFLUX,'(/A11,F17.5,A3,F33.5,A3:,/(F28.5,A3,F33.5,A3))') &
     &      ' ALTITUDES:',(ZM(IK),' KM',IK=1,MLFLX)
          WRITE(IFLUX,'(A8,2A36)')'        ',                           &
     &      ('  ----------------------------------',IK=1,2)
          IF(UNIT.EQ.'M')THEN
              WRITE(IFLUX,'(A8,2A36)')'  WAVLEN',                       &
     &          ('      UPWARD    DOWNWARD      DIRECT',IK=1,2)
              WRITE(IFLUX,'(A8,2A36)')'(MICRON)',                       &
     &          ('     DIFFUSE     DIFFUSE       SOLAR',IK=1,2)
          ELSEIF(UNIT.EQ.'N')THEN
              WRITE(IFLUX,'(A8,2A36)')'  WAVLEN',                       &
     &          ('      UPWARD    DOWNWARD      DIRECT',IK=1,2)
              WRITE(IFLUX,'(A8,2A36)')'   (NM) ',                       &
     &          ('     DIFFUSE     DIFFUSE       SOLAR',IK=1,2)
          ELSE
              WRITE(IFLUX,'(A8,2A36)')'   FREQ ',                       &
     &          ('      UPWARD    DOWNWARD      DIRECT',IK=1,2)
              WRITE(IFLUX,'(A8,2A36)')'  (CM-1)',                       &
     &          ('     DIFFUSE     DIFFUSE       SOLAR',IK=1,2)
          ENDIF
          WRITE(IFLUX,'(A8,2A36)')' -------',                           &
     &      (' ----------- ----------- -----------',IK=1,2)
          FRMT='(F8.2,1P,  2(3(1X,E11.5)):,/(8X,  2(3(1X,E11.5))))'
        ELSE
          WRITE(IFLUX,'(/A11,F17.5,A3,124(F33.5,A3):,                   &
     &      /(125(F28.5,A3:,5X)))')                                     &
     &      ' ALTITUDES:',(ZM(IK),' KM',IK=1,MLFLX)
          WRITE(IFLUX,'((8X,125A36))')                                  &
     &      ('  ----------------------------------',IK=1,MLFLX)
          IF(UNIT.EQ.'M')THEN
              WRITE(IFLUX,'(A8,125A36:,/(8X,125A36))')'  WAVLEN',       &
     &          ('      UPWARD    DOWNWARD      DIRECT',IK=1,MLFLX)
              WRITE(IFLUX,'(A8,125A36:,/(8X,125A36))')'(MICRON)',       &
     &          ('     DIFFUSE     DIFFUSE       SOLAR',IK=1,MLFLX)
          ELSEIF(UNIT.EQ.'N')THEN
              WRITE(IFLUX,'(A8,125A36:,/(8X,125A36))')'  WAVLEN',       &
     &          ('      UPWARD    DOWNWARD      DIRECT',IK=1,MLFLX)
              WRITE(IFLUX,'(A8,125A36:,/(8X,125A36))')'   (NM) ',       &
     &          ('     DIFFUSE     DIFFUSE       SOLAR',IK=1,MLFLX)
          ELSE
              WRITE(IFLUX,'(A8,125A36:,/(8X,125A36))')'   FREQ ',       &
     &          ('      UPWARD    DOWNWARD      DIRECT',IK=1,MLFLX)
              WRITE(IFLUX,'(A8,125A36:,/(8X,125A36))')'  (CM-1)',       &
     &          ('     DIFFUSE     DIFFUSE       SOLAR',IK=1,MLFLX)
          ENDIF
          WRITE(IFLUX,'(A8,125A36:,/(8X,125A36))')' -------',           &
     &      ('  ----------  ----------  ----------',IK=1,MLFLX)
          FRMT='(F8.2,1P,125(3(1X,E11.5)):,/(8X,125(3(1X,E11.5))))'
        ENDIF
      ENDIF

!     WRITE FLUX TABLE DATA FORMAT:
      IF(UNIT.EQ.'M')THEN
          IF(IV1.GT.100)THEN
             FRMT(5:5)='5'
          ELSEIF(IV1.GT.10)THEN
             FRMT(5:5)='4'
          ELSE
             FRMT(5:5)='3'
          ENDIF
      ELSEIF(UNIT.EQ.'N')THEN
          IF(IV1.GT.100)THEN
             FRMT(5:5)='2'
          ELSEIF(IV1.GT.10)THEN
             FRMT(5:5)='1'
          ELSE
             FRMT(5:5)='0'
          ENDIF
      ENDIF

!     RETURN TO TRANS.
      RETURN
      END
