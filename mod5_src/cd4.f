      SUBROUTINE CD4(ICH1,MARIC1,MARK,NLOS,LNFILT,I_SCAN,FWSCAN,DEGALL, &
     &  PLTOUT,SPFLUX,NMFLRT,LNFLRT,FLRT,                               &
     &         TP5_FLAG,TP5_FILE,TP5_LINENUMBER)

!     PROCESS CARD4 INPUTS:
      IMPLICIT NONE

!     PARAMETERS:
!       R8LN2    RECIPROCAL OF 8 LN2.
      REAL R8LN2
      PARAMETER(R8LN2=.1803369)
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       ICH1     NUMERIC LABEL FOR BOUNDARY LAYER AEROSOL MODEL.
!       MARIC1   MARINE AEROSOL REGION NUMBER (=0 IF NOT USED).
!       MARK     RELATIVE HUMIDITY LAYER INDEX FOR MARINE AEROSOL.
!       NLOS     NUMBER OF LINE-OF-SIGHT PATHS.
!       LNFILT   LENGTH OF FILTER RESPONSE FUNCTION FILE NAME.
!       PLTOUT   NAME OF .plt PLOT OUTPUT FILE.
!       SPFLUX   NAME OF .flx SPECTRALLY FILTERED FLUX FILE.
!       NMFLRT   CASE NUMBER.
!       LNFLRT   LENGTH OF FILE ROOT NAME.
!       FLRT     FILE ROOT NAME.
!DRF    TP5_FLAG LOGICAL FLAG IF READING IN TP5 FILE
!DRF    TP5_FILE STRING CONTAINING TP5 FILE INFORMATION
!DRF    TP5_LINENUMBER (I/O) CURRENT LINENUMBER FOR TP5 READING
      INTEGER TP5_LINENUMBER
      LOGICAL TP5_FLAG
      CHARACTER TP5_FILE(1000)*110
      INTEGER ICH1,MARIC1,MARK,NLOS,LNFILT,NMFLRT,LNFLRT
      CHARACTER PLTOUT*(NAMLEN),SPFLUX*(NAMLEN),FLRT*(NAMLEN-4)

!     OUTPUT ARGUMENTS:
!       I_SCAN   SCANNING FUNCTION NUMERICAL LABEL.
!       FWSCAN   DOMAIN OF SCANNING FUNCTION OVER ITS FWHM.
!       DEGALL   LOGICAL FLAG, .TRUE. TO SCAN (DEGRADE) ALL TAPE7 DATA.
      INTEGER I_SCAN
      REAL FWSCAN
      LOGICAL DEGALL

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BMHEAD.h'

!     /CARD1/
!       MODEL    MODEL ATMOSPHERE INDEX.
!       ITYPE    SLANT PATH TYPE.
!       IEMSCT   RADIATIVE TRANSFER MODE.
!                  0 FOR TRANSMITTANCE
!                  1 FOR THERMAL EMISSION ONLY
!                  2 FOR THERMAL EMISSION PLUS SOLAR SCATTER
!                  3 FOR TRANSMITTED SOLAR IRRADIANCE
!                  4 FOR SOLAR SCATTER ONLY
!       M1       MODEL ATMOSPHERE FOR PRESSURE & TEMPERATURE PROFILES.
!       M2       MODEL ATMOSPHERE FOR H2O PROFILE.
!       M3       MODEL ATMOSPHERE FOR O3 PROFILE.
!       I_RD2C   READ CARD 2C, 2C1, ... IF EQUAL 1; SKIP IF EQUAL TO 0.
!       NOPRNT   PRINT FLAG.
!       MODTRN   MODTRAN BAND MODEL FLAG.
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,I_RD2C,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,I_RD2C,NOPRNT,MODTRN

!     /CARD4/
!       IV1      LOWEST SPECTRAL FREQUENCY OUTPUT [CM-1].
!       IV2      HIGHEST SPECTRAL FREQUENCY OUTPUT [CM-1].
!       IDV      PRINTOUT SPECTRAL FREQUENCY STEP SIZE [CM-1].
!       IFWHM    TRIANGULAR SLIT FULL-WIDTH-HALF-MAXIMUM [CM-1].
!       VBAND    CURRENT COMPUTATION BAND FREQUENCY [CM-1].
!                (EQUALS BAND CENTER FOR 1, 5 & 15 CM-1 BAND MODELS;
!                EQUALS THE MINIMUM BAND VALUE FOR 0.1 CM-1 BAND MODEL)
!       IBINPT   BIN NUMBER OF CURRENT SPECTRAL POINT.
!                (CENTER FREQUENCY = IBINPT * BNDWID + OSHIFT).
!       IBINLO   BIN NUMBER OF (PADDED) SPECTRAL RANGE LOWER BOUND.
!       IBINHI   BIN NUMBER OF (PADDED) SPECTRAL RANGE UPPER BOUND.
!       IBINMN   BIN NUMBER OF MINIMUM COMPUTATION SPECTRAL POINT.
!       IBINMX   BIN NUMBER OF MAXIMUM COMPUTATION SPECTRAL POINT.
!       IBINDL   BIN NUMBER INCREMENT FOR SPECTRAL PRINTOUT.
!       IBINRS   BIN NUMBER INCREMENT EQUAL TO SPECTRAL RESOLUTION.
!       IBINOS   BIN NUMBER OFFSET BETWEEN CURRENT & OUTPUT SPC POINTS.
!       IBINWR   BIN NUMBER OF NEXT SPECTRAL DATA WRITE.
!       MBINPT   BIN NUMBER MAXIMUM FOR CURRENT BAND MODEL RESOLUTION.
!       IDBIN5   SPECTRAL BIN NUMBER STEP SIZE FOR 5 CM-1 GRID.
!       ISTEP5   INCREMENT FOR RETRIEVING 5 CM-1 RESOLUTION DATA [CM-1].
!       NSPCDT   NUMBER OF OUTPUT SPECTRAL DATA POINTS.
      DOUBLE PRECISION IDV
      REAL IV1,IV2,IFWHM,VBAND
      INTEGER IBINPT,IBINLO,IBINHI,IBINMN,IBINMX,IBINDL,                &
     &  IBINRS,IBINOS,IBINWR,MBINPT,IDBIN5,ISTEP5,NSPCDT
      COMMON/CARD4/IDV,IV1,IV2,IFWHM,VBAND,IBINPT,IBINLO,IBINHI,IBINMN, &
     &  IBINMX,IBINDL,IBINRS,IBINOS,IBINWR,MBINPT,IDBIN5,ISTEP5,NSPCDT

!     /CFLAGS/
!       YFLAG    Y COORDINATE FLAG FOR plot.dat FILE
!                  = "T" FOR TRANSMITTANCE
!                  = "R" FOR RADIANCE (IRRADIANCE FOR IEMSCT=3)
!                  = "N" FOR NO plot.dat OUTPUT
!       XFLAG    X COORDINATE FLAG FOR plot.dat FILE
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

!     /CNTRL/
!       NSEG     NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
!       THERML   FLAG TO CALCULATE THERMAL SCATTER.
      INTEGER NSEG,ML,MLFLX,IMULT
      LOGICAL THERML
      COMMON/CNTRL/NSEG(0:MLOSP1),ML,MLFLX,IMULT,THERML

!     /CARD2/
!       IHAZE    BOUNDARY LAYER AEROSOL MODEL NUMBER.
!       ISEASN   SEASON NUMBER (1=SPRING-SUMMER, 2=FALL-WINTER).
!       IVULCN   VOLCANIC AEROSOL MODEL NUMBER.
!       ICSTL    COASTAL AIRMASS MODEL NUMBER.
!       ICLD     CLOUD MODEL NUMBER.
!       IVSA     VERTICAL STRUCTURE ALGORITHM (0=OFF, 1=ON).
!       VIS      SURFACE VISIBILITY (GROUND METEOROLOGICAL RANGE) [KM].
!       WSS      CURRENT WIND SPEED (M/S).
!       WHH      24-HOUR WIND SPEED (M/S).
!       RAINRT   RAIN RATE (MM/HR).
!       LSAP     LOGICAL FLAG FOR SPECTRAL AEROSOL PROFILES INPUT.
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA
      REAL VIS,WSS,WHH,RAINRT
      LOGICAL LSAP
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 &
     &  VIS,WSS,WHH,RAINRT,LSAP

!     /SCAN/
!       V1       LOWER BOUND ON SPECTRAL RANGE [CHUNIT DEFINES UNIT].
!       V2       UPPER BOUND ON SPECTRAL RANGE [CHUNIT DEFINES UNIT].
!       DV       SPECTRAL STEP SIZE FOR OUTPUT [CHUNIT DEFINES UNIT].
!       FWHM     FULL-WIDTH-AT-HALF-MAXIMUM [CHUNIT DEFINES UNIT].
!       FWHMSQ   TRIANGULAR SLIT NORMALIZATION FACTOR
!                (EQUALS FWHM SQUARED) [CHUNIT DEFINES UNIT].
!       VOUT     CURRENT SPECTRAL OUTPUT [CHUNIT DEFINES UNIT].
!       M2_FAC   FACTOR IN 2ND MOMENT DIFFERENCE INTEGRAL [CHUNIT UNIT].
      DOUBLE PRECISION V1,V2,DV,VOUT,FWHMSQ
      REAL FWHM,M2_FAC
      COMMON/SCANFN/V1,V2,DV,VOUT,FWHMSQ,FWHM,M2_FAC

!     /JM1A1/
!       RESCHR   BAND MODEL RESOLUTION CHARACTER STRING.
!       DISSTR   CHARACTER STRING USED TO READ IN DISORT LOGICALS.
!       H2OSTR   VERTICAL WATER COLUMN CHARACTER STRING (IF THE
!                FIRST NON-BLANK CHARACTER IS A "G", THE WATER
!                COLUMN IN GM/CM2 FOLLOWS "G"; IF THE FIRST
!                NON-BLANK CHARACTER IS AN "A", THE WATER COLUMN
!                IN ATM-CM AT 273.15K FOLLOWS "A"; OTHERWISE THE
!                STRING CONTAINS A WATER COLUMN SCALING FACTOR).
!       O3STR    VERTICAL OZONE COLUMN CHARACTER STRING (IF THE
!                FIRST NON-BLANK CHARACTER IS A "G", THE OZONE
!                COLUMN IN GM/CM2 FOLLOWS "G"; IF THE FIRST
!                NON-BLANK CHARACTER IS AN "A", THE OZONE COLUMN
!                IN ATM-CM AT 273.15K FOLLOWS "A"; OTHERWISE THE
!                STRING CONTAINS AN OZONE COLUMN SCALING FACTOR).
!       USRSUN   USER-SPECIFIED TOP-OF-ATMOSPHERE SOLAR IRRADIANCE FILE.
!       BMROOT   PREFIX OF MOLECULAR BAND MODEL PARAMETERS FILE.
!       FILTNM   NAME OF FILTER RESPONSE FUNCTION FILE.
!       H2OAER   FLAG, TRUE IF DEFAULT AEROSOL PROPERTIES ARE REVISED
!                BASED ON WATER COLUMN SCALING.
!       DATDIR   NAME OF THE MODTRAN DATA DIRECTORY.
      CHARACTER RESCHR*2,DISSTR*3,H2OSTR*10,O3STR*10,USRSUN*(NAMLEN),   &
     &  FILTNM*(NAMLEN),BMROOT*(NAMLEN),H2OAER*1,DATDIR*(NAMLEN-LENSUN)
      COMMON/JM1A1/RESCHR,DISSTR,H2OSTR,O3STR,USRSUN,                   &
     &  FILTNM,BMROOT,H2OAER,DATDIR

!     /COMNOV/
!       LNOVAM   LOGICAL FLAG, .TRUE. IF NOVAM AEROSOLS ARE USED.
      LOGICAL LNOVAM
      REAL EXTNOV(MNOV,MXWVLN),ABSNOV(MNOV,MXWVLN),                     &
     &  ASMNOV(MNOV,MXWVLN),WLNOV(MXWVLN)
      INTEGER NNOV,NWLNOV
      COMMON/COMNOV/LNOVAM,EXTNOV,ABSNOV,ASMNOV,WLNOV,NNOV,NWLNOV

!     /CSCAN/
!       CHUNIT   UNIT FLAG ('W'=WAVENUMBERS;'M'=MICRONS;'N'=NANOMETERS).
!       RELABS   SPECTRAL RESOLUTION FLAG('A'=ABSOLUTE;'R'=RELATIVE[%]).
!       LNFEED   LINE FEED FLAG FOR .FLX FILE ('T' FOR 80 CHARACTER
!                  LINES, 'F' FOR LONG LINES, ' ' FOR NO .FLX FILE).
      CHARACTER CHUNIT*1,RELABS*1,LNFEED*1
      COMMON/CSCAN/CHUNIT,RELABS,LNFEED

!     /SPCCTL/
!       DODRGD   .TRUE. FOR NON-DEFAULT SPECTRAL ANALYSIS.
      LOGICAL DODGRD
      COMMON/SPCCTL/DODGRD

!     /VRANGE/
!       VBNDMN   COMPUTATIONAL BANDPASS MINIMUM FREQUENCY [CM-1].
!       VBNDMX   COMPUTATIONAL BANDPASS MAXIMUM FREQUENCY [CM-1].
      REAL VBNDMN,VBNDMX
      COMMON/VRANGE/VBNDMN,VBNDMX

!     FUNCTIONS:
!       FILTER   .TRUE. IF FILTER FUNCTION IS SUCCESSFULLY PROCESSED.
      LOGICAL FILTER

!     LOCAL VARIABLES:
!       LOPEN    LOGICAL FLAG, .TRUE. IF FILE IS OPEN.
!       ALAMMX   MAXIMUM WAVELENGTH (MICRONS).
!       VRFRAC   SPECTRAL FREQUENCY USED FOR INDEX OF REFRACTION [CM-1].
      LOGICAL LOPEN
      REAL ALAMMX,VRFRAC

!     DATA:
      CHARACTER FLIST(7)*11
      SAVE FLIST
      DATA FLIST/     'TRIANGULAR ','RECTANGULAR','GAUSSIAN   ',        &
     &  'SINC       ','SINC**2    ','HAMMING    ','USER       '/

!     READ CARD4 INPUTS:
      IF(LJMASS)THEN

!         AUXILIARY OUTPUT FILES ARE NOT GENERATED FOR JMASS:
          CALL INITCD('CARD4')
          YFLAG=' '
          XFLAG=' '
          DLIMIT='        '
          FLAGS='       '
          MLFLX=ML-1
          VRFRAC=0.
      ELSE
          IF (TP5_FLAG) THEN
          READ(IRD,'(4F10.0,2A1,A8,A7,I3,F10.0)')                       &
     &      V1,V2,DV,FWHM,YFLAG,XFLAG,DLIMIT,FLAGS,MLFLX,VRFRAC
          ELSE
          READ(TP5_FILE(TP5_LINENUMBER),'(4F10.0,2A1,A8,A7,I3,F10.0)')  &
     &      V1,V2,DV,FWHM,YFLAG,XFLAG,DLIMIT,FLAGS,MLFLX,VRFRAC
            TP5_LINENUMBER=TP5_LINENUMBER+1
          ENDIF
          CALL UPCASE(FLAGS)
          IF(MLFLX.LE.0 .OR. MLFLX.GE.ML)MLFLX=ML-1
      ENDIF

!     FLAG NON-DEFAULT SPECTRAL SCENARIO.
      DODGRD=(FLAGS(1:4).NE.'    ')

!     INTERPRET FLAGS(2:2) INSTRUMENT RESPONSE (SCAN) FUNCTION INPUT.
      CALL UPCASE(FLAGS(2:2))
      IF(FLAGS(2:2).EQ.'T' .OR. FLAGS(2:2).EQ.'1' .OR.                  &
     &                          FLAGS(2:2).EQ.' ')THEN

!         TRIANGULAR:
          I_SCAN=1
          FWSCAN=1.
          IF(FLAGS(3:3).EQ.'R')THEN
              M2_FAC=FWHM**2/60000
          ELSE
              M2_FAC=FWHM**2/6
          ENDIF
      ELSEIF(FLAGS(2:2).EQ.'R' .OR. FLAGS(2:2).EQ.'2')THEN

!         RECTANGULAR:
          I_SCAN=2
          FWSCAN=.5
          IF(FLAGS(3:3).EQ.'R')THEN
              M2_FAC=FWHM**2/120000
          ELSE
              M2_FAC=FWHM**2/12
          ENDIF
      ELSEIF(FLAGS(2:2).EQ.'G' .OR. FLAGS(2:2).EQ.'3')THEN

!         GAUSSIAN:
          I_SCAN=3
          FWSCAN=5.
          IF(FLAGS(3:3).EQ.'R')THEN
              M2_FAC=R8LN2*FWHM**2/10000
          ELSE
              M2_FAC=R8LN2*FWHM**2
          ENDIF
      ELSEIF(FLAGS(2:2).EQ.'S' .OR. FLAGS(2:2).EQ.'4')THEN

!         SINC:
          I_SCAN=4
          FWSCAN=10.
          M2_FAC=0.
      ELSEIF(FLAGS(2:2).EQ.'C' .OR. FLAGS(2:2).EQ.'5')THEN

!         SINC-SQUARED:
          I_SCAN=5
          FWSCAN=10.
          M2_FAC=0.
      ELSEIF(FLAGS(2:2).EQ.'H' .OR. FLAGS(2:2).EQ.'6')THEN

!         HAMMING:
          I_SCAN=6
          FWSCAN=10.
          IF(FLAGS(3:3).EQ.'R')THEN
              M2_FAC=2.500957E-4*FWHM**2
          ELSE
              M2_FAC=2.500957*FWHM**2
          ENDIF
      ELSEIF(FLAGS(2:2).EQ.'U' .OR. FLAGS(2:2).EQ.'7')THEN

!         USER-DEFINED:
          I_SCAN=7
          FWSCAN=10.
          M2_FAC=0.
      ELSE

!         TRIANGULAR:
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'Error in CD4:  Bad input for I_SCAN in CARD 4 input.'
      ENDIF

!     TEST THE VALIDITY OF SPECTRAL INPUTS, DETERMINE THE PADDING
!     NECESSARY FOR V1 & V2, AND SET IBINLO, IBINHI, IBINDL & IBINRS.
      CALL CHKRES(RESCHR,FWSCAN,DEGALL)

      IV1=IBINLO*BNDWID
      IV2=IBINHI*BNDWID
      IDV=IBINDL*DBLE(BNDWID)
      IFWHM=IBINRS*BNDWID
      IF(MODTRN .AND. IBINRS.GT.50)THEN
          WRITE(IPR,'((A,F12.4,A))')' Warning from CD4:  The'//         &
     &      ' spectral resolution (IFWHM =',IFWHM,' CM-1) cannot',      &
     &      ' exceed the MODTRAN band model resolution (BNDWID =',      &
     &      BNDWID,' CM-1) by more than a factor of 50.  IFWHM',        &
     &      ' is being reduced to',50*BNDWID,' CM-1.'
          IFWHM=50*BNDWID
      ENDIF
      IF(IBINHI.LT.IBINLO+IBINDL)THEN
          WRITE(IPR,'(/A)')' IV2 WAS LESS THAN IV1 + IDV AND RESET.'
          IBINHI=IBINLO+IBINDL
          IV2=IBINHI*BNDWID
      ENDIF
      IF(.NOT.LJMASS)THEN
          WRITE(IPR,'(/A,4F10.1,2A)')' CARD 4  *****',                  &
     &      BNDWID*IBINLO,BNDWID*IBINHI,BNDWID*IBINDL,                  &
     &      IFWHM,'   SCANNING FUNCTION:',FLIST(I_SCAN)
          IF(IBINLO.GT.0)THEN
              ALAMMX=10000/(IBINLO*BNDWID)
          ELSE
              ALAMMX=999999.99
          ENDIF
          WRITE(IPR,'(/A,2(/13X,A,F10.1,A,F12.4,A),/(11X,A,F10.1,A))')  &
     &      ' FREQUENCY RANGE','IV1 =',BNDWID*IBINLO,' CM-1  (',        &
     &      ALAMMX,' MICRONS)','IV2 =',BNDWID*IBINHI,' CM-1  (',        &
     &      10000/(BNDWID*IBINHI),' MICRONS)','  IDV =',                &
     &      BNDWID*IBINDL,' CM-1','IFWHM =',IFWHM,' CM-1'
      ENDIF

!     CHECK MARINE AEROSOL:
      IF(IHAZE.EQ.3 .AND. IBINLO*BNDWID.LT.250.05)THEN
          WRITE(IPR,'(//A,/19X,A)')' Warning from MC4:  NAVY HAZE'//    &
     &      ' MODEL CANNOT BE USED BELOW 250 CM-1.',' PROGRAM WILL'     &
     &      //' SWITCH TO LOWTRAN 5 MARITIME HAZE MODEL (IHAZE=4).'
          IHAZE=4
      ENDIF

      IF(VRFRAC.LE.0.)THEN
          VRFRAC=(IBINLO+IBINHI)*HBNDWD
      ELSEIF(CHUNIT.EQ.'M')THEN
          VRFRAC=10000/VRFRAC
      ELSEIF(CHUNIT.EQ.'N')THEN
          VRFRAC=1.E7/VRFRAC
      ENDIF

!     CHECK PLOT.DAT FILE FLAGS
      IF(YFLAG.EQ.'R' .OR. YFLAG.EQ.'r')THEN

!         WRITE SPECTRAL RADIANCES (TRANSMITTANCES IF IEMSCT=0
!         OR TRANSMITTED SOLAR IRRADIANCES IF IEMSCT=3)
          YFLAG='R'
          IF(IEMSCT.EQ.0)YFLAG='T'
      ELSEIF(YFLAG.EQ.'T' .OR. YFLAG.EQ.'t')THEN

!         WRITE TRANSMITTANCES
          YFLAG='T'
      ELSE

!         DO NOT WRITE TO THE PLOT.DAT FILE.
          YFLAG='N'
      ENDIF
      IF(XFLAG.EQ.'N' .OR. XFLAG.EQ.'n')THEN
          XFLAG='N'
      ELSEIF(XFLAG.EQ.'M' .OR. XFLAG.EQ.'m')THEN
          XFLAG='M'
      ELSE
          XFLAG='W'
      ENDIF
      IF(YFLAG.NE.'N')THEN

!         bin out
          IF(BINOUT)THEN
              INQUIRE(JPLOT,OPENED=LOPEN)
              IF(.NOT.LOPEN)CALL OPNFL(JPLOT,0,PLTOUT,'UNKNOWN',        &
     &          'UNFORMATTED','CD4')
          ELSE
              INQUIRE(IPLOT,OPENED=LOPEN)
              IF(.NOT.LOPEN)                                            &
     &          CALL OPNFL(IPLOT,0,PLTOUT,'UNKNOWN','FORMATTED','CD4')
          ENDIF
      ENDIF
      IF(FLAGS(7:7).EQ.'T' .OR. FLAGS(7:7).EQ.'F')THEN
          IF(BINOUT)THEN

!             for binary output put all auxiliary files in tape8
              INQUIRE(JPR1,OPENED=LOPEN)
              IF(.NOT.LOPEN)THEN
                  IF(LNFLRT.LE.0)THEN
                      CALL OPNFL(JPR1,0,'tape8b',                       &
     &                  'UNKNOWN','UNFORMATTED','CD4')
                  ELSE
                      CALL OPNFL(JPR1,0,FLRT(1:LNFLRT)//'b.tp8',        &
     &                  'UNKNOWN','UNFORMATTED','CD4')
                  ENDIF
              ENDIF
          ELSE
              INQUIRE(IFLUX,OPENED=LOPEN)
              IF(.NOT.LOPEN)                                            &
     &          CALL OPNFL(IFLUX,0,SPFLUX,'UNKNOWN','FORMATTED','CD4')
          ENDIF
      ENDIF

!     LOAD ATMOSPHERIC PROFILE INTO /MPROF/
      CALL STDMDL(ICH1,MARIC1,MARK,H2OAER,VRFRAC,LNOVAM)

!     FILTER RESPONSE FUNCTION:
      IF(LNFILT.GT.0)THEN
          IF(FILTER(FILTNM,LNFILT,NMFLRT,NLOS))WRITE(IPR,'(/(A))')      &
     &      ' THE FILTER RESPONSE FUNCTION FILE ',                      &
     &      FILTNM(1:LNFILT),' WAS SUCCESSFULLY PROCESSED.'
      ENDIF

!     DEFINE SPECTRAL INCREMENT VARIABLES:
      IF(IBINRS.EQ.1)THEN
          IBINMN=IBINLO
      ELSE
          IBINMN=MAX(0,IBINLO-IBINRS)
      ENDIF

!     REQUIRE SPECTRAL BIN TIMES BNDWID TO BE A MULTIPLE OF 5 CM-1:
      IF(RESCHR.EQ.'p1')THEN

!         START VCEN FOR THE 'p1' BAND MODEL AT XXX4.95 OR XXX9.95:
          ISTEP5=5
          IDBIN5=50
          IBINMN=50*(IBINMN/50)-1
      ELSEIF(RESCHR.EQ.'01')THEN
          ISTEP5=5
          IDBIN5=5
          IBINMN=5*(IBINMN/5)
      ELSE
          ISTEP5=NINT(BNDWID)
          IDBIN5=1
      ENDIF

!     SLIT AND INTEGRATION WEIGHTING FACTOR INITIALIZATIONS:
      IF(IBINRS.GT.1)THEN

!         TRIANGULAR SLIT:
          IF(RESCHR.EQ.'p1')THEN
              IBINOS=IBINRS-1
          ELSE
              IBINOS=IBINRS
          ENDIF
          IBINMX=MIN(IBINHI+IBINOS,MBINPT)
      ELSE

!         SQUARE SLIT:
          IBINOS=0
          IF(RESCHR.EQ.'p1')THEN
              IBINMX=MIN(IBINHI+IBINOS-1,MBINPT)
          ELSE
              IBINMX=MIN(IBINHI+IBINOS,MBINPT)
          ENDIF
      ENDIF

!     FREQUENCY RANGE - MINIMUM AND MAXIMUM BANDPASS FREQUENCY:
      VBNDMN=IBINMN*BNDWID+OSHIFT
      IF(VBNDMN.LE.0.)VBNDMN=OSHIFT
      VBNDMN=VBNDMN-BNDWID/2
      VBNDMX=(IBINMX+.5)*BNDWID+OSHIFT

!     RETURN TO DRIVER:
      RETURN
      END
