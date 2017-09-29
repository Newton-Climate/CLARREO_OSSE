      SUBROUTINE CHKRES(FUNC,FWFUNC,IFUNC,DEGALL)

!     ARGUMENTS:
!       IFUNC    SCANNING FUNCTION NUMERICAL LABEL.
!       FWFUNC   NUMBER OF FULL-WIDTH AT HALF-MAXIMUM (FWHM)
!                USED TO INTEGRATE SCANNING FUNCTION
      INTEGER IFUNC
      REAL FWFUNC
      LOGICAL DEGALL
      CHARACTER FUNC*11

!     PARAMETERS:
!       THE FOLLOWING PARAMETER ALSO APPEARS IN SUBROUTINE DGRD
!       CHANGE IT THERE IF YOU CHANGE IT HERE
      INTEGER MDDGRD
      PARAMETER (MDDGRD=50000)
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BMHEAD.h'

!     /CARD1/
!       MODEL    ATMOSPHERIC PROFILE MODEL NUMBER.
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

!     COMMON /CFLAGS/
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

!     COMMON/CSCAN/
!       CHUNIT   UNIT FLAG ('W'=WAVENUMBERS;'M'=MICRONS;'N'=NANOMETERS).
!       RELABS   SPECTRAL RESOLUTION FLAG('A'=ABSOLUTE;'R'=RELATIVE[%]).
!       LNFEED   LINE FEED FLAG FOR .FLX FILE ('T' FOR 80 CHARACTER
!                  LINES, 'F' FOR LONG LINES, ' ' FOR NO .FLX FILE).
      CHARACTER CHUNIT*1,RELABS*1,LNFEED*1
      COMMON/CSCAN/CHUNIT,RELABS,LNFEED

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:
      INTEGER ISCAN

!     LOCAL VARIABLES:
!       WIDTH    MINIMUM SPECTRAL RESOLUTION [CM-1].
!       PAD      FILTER SPECTRAL PADDING [CHUNIT DEFINES UNIT].
!       BIN      BAND MODEL SPECTRAL RESOLUTION [CM-1].
      INTEGER NRD
      REAL WIDTH,PAD,BIN,WNMIN,WNMAX,DUMB
      LOGICAL LABS

!     DATA:
      INTEGER MAXF
      CHARACTER FLIST(7)*11
      REAL FWLST(7)
      DATA MAXF/50000/
      DATA FLIST(1)/'TRIANGULAR '/,FWLST(1)/ 1.0/,                      &
     &     FLIST(2)/'RECTANGULAR'/,FWLST(2)/ 0.5/,                      &
     &     FLIST(3)/'GAUSSIAN   '/,FWLST(3)/ 5.0/,                      &
     &     FLIST(4)/'SINC       '/,FWLST(4)/10.0/,                      &
     &     FLIST(5)/'SINC**2    '/,FWLST(5)/10.0/,                      &
     &     FLIST(6)/'HAMMING    '/,FWLST(6)/10.0/,                      &
     &     FLIST(7)/'USER       '/,FWLST(7)/10.0/
      SAVE MAXF,FLIST,FWLST

!     USER INPUT INSTRUMENT BANDPASS FREQUENCY (INPUTS):
!     V1 IS NUMERICALLY MINIMUM IN CM-1, MICRONS OR NM
!     V2 IS NUMERICALLY MAXIMUM

!     UNITS OF SCANNING FUNCTION, V1 AND V2 (INPUTS):
!     CHUNIT='W' (CM-1)
!     CHUNIT='M' (MICRON)
!     CHUNIT='N' (NM)

!     RELATIVE OR ABSOLUTE FWHM (INPUTS):
!     RELABS='A' OR BLANK; FWHM IS IN CM-1, MICRON OR NM
!     RELABS='R'; FWHM IS IN PERCENT

!     MODTRAN WILL BE RUN AT 1 CM-1 INTERVALS
!     IF NOT, MODTRAN WILL BE RUN AT THE GREATEST RES ALLOWED.
!     MUST CHECK IF THROUGHOUT THE SPECTRUM,
!     DESIRED FULL-WIDTH-AT-HALF-MAXIMUM IS .GE. IBNDWD (IN CM-1)

!     SUBROUTINE RETURNS (OUTPUTS):
!     STOPS IF CHECKS FAIL
!     FUNC (SCANNING FUNCTION)
!     IV1 (MINIMUM FREQUENCY FOR MODTRAN CALCULATION IN CM-1)
!     IV2 (MAXIMUM FREQUENCY FOR MODTRAN CALCULATION IN CM-1)
!     IDV
!     V1 AND V2 MAY BE REDEFINED IF REQUIRED

!     IF FLAGS(4:4)='A', DEGRADE EVERYTHING IN TAPE7:
      DEGALL=.FALSE.
      IF(FLAGS(4:4).EQ.'A')DEGALL=.TRUE.

!     CHECK RELATIVE/ABSOLUTE FLAG:
      RELABS=FLAGS(3:3)
      IF(RELABS.EQ.' ')THEN

!         SET RELABS TO DEFAULT VALUE:
          RELABS='A'
      ELSEIF(RELABS.NE.'A' .AND. RELABS.NE.'R')THEN

!         RELABS MUST BE 'R' OR 'A'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'BAD INPUT FOR RELABS IN CARD 4'
      ENDIF

!     CHECK RESPONSE FUNCTION FLAG:
      IFUNC=ISCAN(FLAGS(2:2))
      IF(IFUNC.EQ.-99)THEN
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'BAD INPUT FOR IFUNC IN CARD 4'
      ENDIF

!     CHECK UNIT FLAG:
      CHUNIT=FLAGS(1:1)
      IF(CHUNIT.EQ.' ')CHUNIT='W'
      IF(CHUNIT.NE.'W' .AND. CHUNIT.NE.'M' .AND. CHUNIT.NE.'N')THEN
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'BAD INPUT FOR IFUNC IN CARD 4'
      ENDIF

!     CARRIAGE RETURN FLAG (FOR .flx FILE).
      IF(FLAGS(7:7).EQ.'T' .OR. FLAGS(7:7).EQ.'F')THEN
          LNFEED=FLAGS(7:7)
      ELSE
          LNFEED=' '
      ENDIF

!     ANALYZE VALIDITY AND MEANING OF NUMBERS IN CARD4.
!     CONVERT TO EQUIVALENT INTEGERS IN CM-1 PRIOR TO THE MODS.
!     THE PURPOSE OF THESE MODS ARE ENUMERATED BELOW:

!     1.  ALLOW USE OF AN INSTRUMENT RESPONSE (SCANNING) FUNCTION.
!         THIS FUNCTION CAN BE (IFUNC VALUE IN PARENTHESES)
!         TRIANGLE (1 OR " "), SQUARE OR RECTANGULAR (2), GAUSSIAN (3),
!         SINC (4), SINC**2 (5), HAMMING (6) AND USER-DEFINED (7)

!     2.  CHUNIT DENOTES THE UNITS OF FREQUENCY AND RADIANCE
!         CHUNIT= 'W', V1, V2, ETC IN CM-1, RADIANCE IN W/CM^2/SR/CM-1
!         CHUNIT= 'M', V1, ETC IN MICRON, RADIANCE IN W/CM^2/SR/MICRON
!         CHUNIT= 'N', V1, ETC IN NM, RADIANCE IN MICROWATT/CM^2/SR/NM

!     3.  RELABS = 'A' OR " ", FULL-WIDTH-AT-HALF-MAX IS ABSOLUTE
!         RELABS='R', FWHM IS IN PERCENT

!     SWAP V2 AND V1 IF OUT OF ORDER:
      IF(V2.LT.V1)THEN
         DUMB=V1
         V1=V2
         V2=DUMB
      ENDIF

!     WNMIN (MINIMUM IN CM-1 OF THE INPUT (TAPE5) BANDPASS)
!     WNMAX (MAXIMUM IN CM-1 OF THE INPUT (TAPE5) BANDPASS)
      IF(CHUNIT.EQ.'W')THEN
         WNMIN=V1
         WNMAX=V2
      ELSEIF(CHUNIT.EQ.'M')THEN
         WNMIN=1.E4/V2
         WNMAX=1.E4/V1
      ELSEIF(CHUNIT.EQ.'N')THEN
         WNMIN=1.E7/V2
         WNMAX=1.E7/V1
      ENDIF
      IF(WNMIN.GT.MXFREQ .OR. .NOT.MODTRN)THEN

!        WE HAVE LOWTRAN
         WRITE(*,*)'CHKRES:  LOWTRAN (NOT MODTRAN) BAND MODEL'
         BIN=5.0
      ELSE
         BIN=IBNDWD
      ENDIF

      FUNC=FLIST(IFUNC)
      FWFUNC=FWLST(IFUNC)

      LABS=RELABS.EQ.CHAR(65) .OR. RELABS.EQ.CHAR(97)
!     LABS IS TRUE IF RELABS IS BLANK OR 'A'

      IF(LABS .AND. CHUNIT.EQ.'W')THEN

!        IN CM-1; ABSOLUTE RESOLUTION
         WIDTH=FWHM

         IF(.NOT.MODTRN .AND. FLAGS(1:4).EQ.'    ')WIDTH=BIN
!        IF(LOWTRAN .AND. FLAGS(1:4).EQ.'    ')THE PROGRAM WILL BE HERE
!        BUT THIS IS CONVENTIONAL LOWTRAN AND FWHM SHOULD BE IGNORED
!        SO TO MAKE SURE THAT WIDTH IS NOT .LT. BIN, SET WIDTH TO BIN
!        BIN IS ALREADY 5.0 AT THIS PT

         PAD=FWFUNC*FWHM
         IF(WIDTH.LT.BIN)THEN
            WRITE(IPR,'(/A,/17X,A,F12.5,A)')'Error in CHKRES: '//       &
     &        ' Specified spectral resolution (FWHM on CARD4) is too',  &
     &        ' fine.  It must be greater or equal to',BIN,' CM-1.'
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP 'Error in CHKRES:  Spectral resolution is too fine.'
         ENDIF

!        IF(FLAGS(1:4).EQ.0) JUST SKIP THIS ROUTINE AT THIS POINT
         IF(FLAGS(1:4).EQ.'    ')GOTO 100

!        THE USER WANTS TABLES OF NUMBERS FROM V1 TO V2
!        NEED TO PAD WNMIN AND WNMAX
!        MODTRAN WILL RUN FROM IV1 AND IV2
         IV1=INT(WNMIN-PAD-BIN)
         IV2=INT(WNMAX+PAD+BIN)
!        +/-BIN IS FOR EXTRA PADDING IN CASE SOMETHING GOES WRONG

!        FOR TRADITIONAL AND OTHER REASONS,
!        IV1 AND IV2 SHOULD BE MADE MULTIPLES OF 5
         IV1=5*INT(IV1/5)
         IF(MOD(INT(IV2+.01),5).NE.0) IV2=5*(INT(IV2/5)+1)

         IF(IV1.LT.0)THEN
!           MUST REDEFINE V1
            IV1=0
            V1=PAD+BIN
         ENDIF
         IF(IV2.GT.MAXF)THEN

!           MUST REDEFINE V2
!           FIRST CHOOSE IV2 TO BE MULTIPLE OF 5
            IV2=MAXF-MOD(MAXF,5)
            V2=IV2-PAD-BIN
         ENDIF

      ELSEIF(.NOT.LABS .AND. CHUNIT.EQ.'W')THEN

!        WNMIN AND WNMAX ARE IN CM-1; FWHM IS IN PERCENT
!        ABSOLUTE RESOLUTION WILL BE THE FINEST AT WNMIN
         WIDTH=WNMIN*(FWHM/100)
         IF(WIDTH.LT.BIN)THEN
            WRITE(IPR,'(/A,/17X,A,F12.5,A)')'Error in CHKRES: '//       &
     &        ' Specified spectral resolution (FWHM on CARD4) is too',  &
     &        ' fine.  It must be greater or equal to',100*BIN/WNMIN,'%'
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP 'Error in CHKRES:  Spectral resolution is too fine.'
         ENDIF
         PAD=WNMIN*(FWHM/100.)*FWFUNC
         IF(PAD.GT.WNMIN)THEN
!           NO HOPE BECAUSE (FWHM/100)*FWFUNC IS ALWAYS .GE. 1
!           SO PAD WILL ALWAYS EXCEED WNMIN NO MATTER WHAT
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP 'CHKRES: BAD CHOICE OF FREQUENCIES/RESOLUTIONS'
         ELSEIF(PAD.GE.WNMIN-BIN)THEN
            IV1=0
!           COULD NOT HAVE THE EXTRA SAFE PADDING BY AMOUNT EQUAL
!           TO BIN.  SHOULD WORK FOR MOST CASES; SO HOPE FOR BEST
         ELSE
            IV1=INT(WNMIN-PAD-BIN)
!           -BIN IS FOR EXTRA PADDING IN CASE SOMETHING GOES WRONG
         ENDIF
         PAD=WNMAX*(FWHM/100.)*FWFUNC
         IV2=INT(WNMAX+PAD+BIN)
!        +BIN IS FOR EXTRA PADDING IN CASE SOMETHING GOES WRONG

!        FOR TRADITIONAL AND OTHER REASONS,
!        IV1 & IV2 SHOULD BE MULTIPLES OF 5
         IV1=5*INT(IV1/5)
         IF(MOD(INT(IV2+.01),5).NE.0)IV2=5*(INT(IV2/5)+1)

         IF(IV2.GT.MAXF)THEN

!           REDEFINE V2
!           FIRST CHOOSE IV2 TO BE MULTIPLE OF 5
            IV2=MAXF-MOD(MAXF,5)
            V2=(IV2-BIN)/(1+(FWHM/100.)*FWFUNC)
         ENDIF

      ELSEIF(LABS .AND. CHUNIT.EQ.'M')THEN

!        FWHM IS A CONSTANT IN MICRONS.  ITS VALUE IN CM-1 IS SPECTRALLY
!        DEPENDENT: FWHM[uM] / WAVLEN[uM] = FWHM[cm-1] / FREQ[cm-1].
!        THE FINEST RESOLUTION IN CM-1 OCCURS AT FREQ[cm-1] = WNMIN.
         PAD=FWHM*FWFUNC

!        CHECK FOR SUFFICIENT SPECTRAL RESOLUTION:
         IF(FWHM.LT.10000*BIN/WNMIN**2)THEN
            WRITE(IPR,'(/A,/17X,A,F12.5,A)')'Error in CHKRES: '//       &
     &        ' Specified spectral resolution (FWHM on CARD4) is too',  &
     &        ' fine.  It must be greater or equal to',                 &
     &        10000*BIN/WNMIN**2,' MICRONS.'
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP 'Error in CHKRES:  Spectral resolution is too fine.'
         ENDIF

!        THE USER WANTS TABLES OF NUMBERS FROM V1 TO V2 [MICRONS].
!        PAD WNMIN AND WNMAX. MODTRAN WILL RUN FROM IV1 AND IV2
         IV1=INT(10000/(V2+PAD)-.05-BIN)
         IF(IV1.LT.0)THEN
            IV1=0
            V2=10000/BIN-PAD
         ENDIF

!        LOWER IV1 MAKING IT A MULTIPLE OF 5 (REQUIRED
!        FOR INTERPOLATION OF COARSE RESOLUTION DATA):
         IV1=5*INT(IV1/5)

         IF(MAXF*(V1-PAD).LT.10000.)THEN

!           SET IV2 TO A BIG NUMBER
            IV2=2*MAXF
         ELSE
            IV2=INT(10000/(V1-PAD)+.15+BIN)
         ENDIF

!        MAKE IV2 A MULTIPLE OF 5
         IF(MOD(INT(IV2+.01),5).NE.0)IV2=5*(INT(IV2/5)+1)
         IF(IV2.LT.0.)WRITE(*,'(/(A))')' NOT ENOUGH ROOM FOR PADDING',  &
     &     ' BAD FREQUENCY/RESOLUTION INPUTS',                          &
     &     ' LOWER WAVELENGTH WILL BE INCREASED'

         IF(IV2.GT.MAXF .OR. IV2.LE.0.)THEN

!           FIRST MAKE IV2 A MULTIPLE OF 5
            IV2=MAXF-MOD(MAXF,5)
            V1=10000/(IV2-BIN)+PAD
            IF(V1.GE.V2)THEN
               IF(LJMASS)CALL WRTBUF(FATAL)
               WRITE(*,'(/A,F10.2,A,F10.2)')' V1 =',V1,' V2 =',V2
               WRITE(IPR,'(/A,F10.2,A,F10.2)')' V1 =',V1,' V2 =',V2
               STOP 'CHKRES: V1 WAS ADJUSTED, MAKING V1 .GE. V2'
            ENDIF
         ENDIF

      ELSEIF(LABS .AND. CHUNIT.EQ.'N')THEN
         
!        FWHM IS ABSOLUTE AND IN NM
!        RESOLUTION IN CM-1 WILL BE DIFFERENT IN WNMIN AND WNMAX
!        RES_WL/WL=RES_WN/WN IMPLIES RES_WN=RES_WL*(WN/WL)
!        WL (WAVELENGTH) IS IN NM AND WN (WAVENUMBER) IS IN CM-1
!        THEREFORE, RES_WN=RES_WL*(WN^2/1.E7)
!        THIS MEANS THAT FINEST RESOLUTION IS AT WNMIN

         WIDTH=FWHM*WNMIN**2/1.E7
         PAD=FWHM*FWFUNC
         !PAD=15
         IF(WIDTH.LT.BIN)THEN
            WRITE(IPR,'(/A,/17X,A,F12.5,A)')'Error in CHKRES: '//       &
     &        ' Specified spectral resolution (FWHM on CARD4) is too',  &
     &        ' fine.  It must be greater or equal to',                 &
     &        1.E7*BIN/WNMIN**2,' NM.'
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP 'Error in CHKRES:  Spectral resolution is too fine.'
         ENDIF
!        THE USER WANTS TABLES OF NUMBERS FROM V1 TO V2
!        NEED TO PAD WNMIN AND WNMAX
!        MODTRAN WILL RUN FROM IV1 AND IV2
         IV1=INT(1.E7/(V2+PAD)-.05-BIN)

!        MAKE IV1 A MULTIPLE OF 5
         IV1=5*INT(IV1/5)

         IF(PAD.EQ.V1)THEN
            IV2=MAXF*2
!           IF(PAD.EQ.V1) IV2=A BIG NUMBER
         ELSE
            IV2=INT(1.E7/(V1-PAD)+.15+BIN)
         ENDIF
		!WRITE(*,*) 'PAD = ',PAD
		!WRITE(*,*) 'IV1 = ',IV1
		!WRITE(*,*) 'IV2 = ',IV2
!        +/-BIN IS FOR EXTRA PADDING IN CASE SOMETHING GOES WRONG

!        MAKE IV2 A MULTIPLE OF 5
         IF(MOD(INT(IV2+.01),5).NE.0)IV2=5*(INT(IV2/5)+1)

         IF(IV2.LT.0.)WRITE(*,'(/(A))')' NOT ENOUGH ROOM FOR PADDING',  &
     &     ' BAD FREQUENCY/RESOLUTION INPUTS',                          &
     &     ' LOWER WAVELENGTH WILL BE INCREASED'

         IF(IV1.LT.0)THEN
            IV1=0
            V2=1.E7/BIN-PAD
         ENDIF
         IF(IV2.GT.MAXF .OR. IV2.LT.0.)THEN

!           MAKE IV2 MULTIPLE OF 5
            IV2=MAXF-MOD(MAXF,5)
            V1=1.E7/(IV2-BIN)+PAD
            IF(V1.GT.V2)THEN
               IF(LJMASS) CALL WRTBUF(FATAL)
               WRITE(IPR,'(/A,F10.2,A,F10.2)')' V1 =',V1,' V2 =',V2
               STOP 'CHKRES: V1 WAS ADJUSTED, MAKING V1 .GE. V2'
            ENDIF
!DRF		WRITE(*,*) 'IV2 = ',IV2
         ENDIF

      ELSEIF(.NOT.LABS .AND. (CHUNIT.EQ.'M' .OR. CHUNIT.EQ.'N'))THEN

!        FWHM IS IN PERCENT MICRON OR NM
!        RES_WL/WL=PERCENT/100=RES_WN/WN
!        RES_WN=(PERCENT/100)*WN
!        THIS MEANS THAT FINEST RESOLUTION IS AT WNMIN
         WIDTH=(FWHM/100)*WNMIN
         IF(WIDTH.LT.BIN)THEN
            WRITE(IPR,'(/A,/17X,A,F12.5,A)')'Error in CHKRES: '//       &
     &        ' Specified spectral resolution (FWHM on CARD4) is too',  &
     &        ' fine.  It must be greater or equal to',100*BIN/WNMIN,'%'
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP 'Error in CHKRES:  Spectral resolution is too fine.'
         ENDIF
         PAD=V2*(FWHM/100.)*FWFUNC
         IF(CHUNIT.EQ.'M')THEN
            IV1=INT(1.E4/(V2+PAD)-.05-BIN)
            IF(IV1.LT.0)THEN
               IV1=0
!              V2+PAD=1.E4/BIN
!              V2*[1+(FWHM/100.)*FWFUNC]=1.E4/BIN
               V2=(1.E4/BIN)/(1+(FWHM/100.)*FWFUNC)
            ENDIF
         ELSEIF(CHUNIT.EQ.'N')THEN
            IV1=INT(1.E7/(V2+PAD)-.05-BIN)
            IF(IV1.LT.0)THEN
               IV1=0
!              V2+PAD=1.E7/BIN
!              V2*[1+(FWHM/100.)*FWFUNC]=1.E7/BIN
               V2=(1.E7/BIN)/(1+(FWHM/100.)*FWFUNC)
            ENDIF
         ENDIF

!        MAKE IV1 A MULTIPLE OF 5
         IV1=5*INT(IV1/5)

         PAD=V1*(FWHM/100.)*FWFUNC
         IF(PAD.GE.V1)THEN

!           (FWHM/100)*FWFUNC EXCEEDS OR EQUALS 1:
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP 'CHKRES: BAD CHOICE OF FREQUENCIES/RESOLUTIONS'
         ELSEIF(CHUNIT.EQ.'M')THEN
            IV2=INT(1.E4/(V1-PAD)+.15+BIN)

!           MAKE IV2 A MULTIPLE OF 5
            IF(MOD(INT(IV2+.01),5).NE.0)IV2=5*(INT(IV2/5)+1)

            IF(IV2.GT.MAXF)THEN

!              MUST REDEFINE V2
!              FIRST CHOOSE IV2 TO BE MULTIPLE OF 5
               IV2=MAXF-MOD(MAXF,5)
!              V1-PAD=V1*[1-(FWHM/100)*FWFUNC]=1.E4/(IV2-BIN)
               V1=(1.E4/(IV2-BIN))/(1-(FWHM/100)*FWFUNC)
               IF(V1.GE.V2)THEN
                  IF(LJMASS)CALL WRTBUF(FATAL)
                  WRITE(IPR,'(/A,F10.2,A,F10.2)')' V1 =',V1,' V2 =',V2
                  STOP 'CHKRES: V1 WAS ADJUSTED, MAKING V1 .GE. V2'
               ENDIF
            ENDIF
         ELSEIF(CHUNIT.EQ.'N')THEN
            IV2=INT(1.E7/(V1-PAD)+.15+BIN)
!DRF		WRITE(*,*) 'IV2 = ',IV2
!           MAKE IV2 A MULTIPLE OF 5
            IF(MOD(INT(IV2+.01),5).NE.0) IV2=5*(INT(IV2/5)+1)

            IF(IV2.GT.MAXF)THEN

!              MUST REDEFINE V2
!              FIRST CHOOSE IV2 TO BE MULTIPLE OF 5
               IV2=MAXF-MOD(MAXF,5)
!              V1-PAD=V1*[1-(FWHM/100)*FWFUNC]=1.E7/(IV2-BIN)
               V1=(1.E7/(IV2-BIN))/(1-(FWHM/100)*FWFUNC)
               IF(V1.GE.V2)THEN
                  IF(LJMASS)CALL WRTBUF(FATAL)
                  WRITE(IPR,'(/A,F10.2,A,F10.2)')' V1 =',V1,' V2 =',V2
                  STOP 'CHKRES: V1 WAS ADJUSTED, MAKING V1 .GE. V2'
               ENDIF
            ENDIF
!           +/-BIN IS FOR EXTRA PADDING IN CASE SOMETHING GOES WRONG
         ENDIF
      ENDIF

!     CHECK ON ARRAY SIZES

      NRD=INT((V2-V1)/DV)+1

      IF(NRD.EQ.1)THEN
          WRITE(IPR,'(/A)')'Error in CHKRES: Need more spectral points.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'Error in CHKRES: Need more spectral points.'
      ELSEIF(NRD.GT.MDDGRD)THEN
          WRITE(IPR,'(/A,2(I8,A),/17X,A)')'Error in CHKRES: '//         &
     &      ' Number of spectral points required for convolution (',NRD,&
     &      ') exceeds the maximum (MDDGRD =',MDDGRD,'.  Decrease',     &
     &      ' the frequency interval or increase parameter MDDGRD.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'Error in CHKRES:  Spectral resolution is too fine.'
      ENDIF

 100  CONTINUE
      IF(FLAGS(1:4).EQ.'    ')THEN
         IFWHM=INT(FWHM+0.001)
         IDV=INT(DV+0.001)
         IV1=INT(WNMIN+0.001)
         IV2=INT(WNMAX+0.001)
      ELSE
         IFWHM=INT(BIN)
         IDV=INT(BIN)
      ENDIF
!DRF	WRITE(*,*) 'IV2 = ',IV2
      END

      INTEGER FUNCTION ISCAN(C)
      CHARACTER C*1, UPCASE*1

!     INPUT IS A CHARACTER WHICH CAN BE
!     THE UPPER CASE ALPHABET, BLANK OR AN INTEGER.
!     C CAN BE:
!     'T', 1 OR BLANK (TRIANGULAR)
!     'R', 2 (RECTANGULAR)
!     'G', 3 (GAUSSIAN)
!     'S', 4 (SINC)
!     'C', 5 (SINC**2)
!     'H', 6 (HAMMING)
!     'U', 7 (USER)

!     DEPENDING ON C, ISCAN WILL RETURN THE INTEGER VALUE.
!     IF C IS NONE OF THE ABOVE ISCAN WILL RETURN -99 SIGNIFYING ERROR

      C=UPCASE(C)
      IF(C.EQ.'T' .OR. C.EQ.'1' .OR. C.EQ.' ')THEN
         ISCAN=1
      ELSEIF(C.EQ.'R' .OR. C.EQ.'2')THEN
         ISCAN=2
      ELSEIF(C.EQ.'G' .OR. C.EQ.'3')THEN
         ISCAN=3
      ELSEIF(C.EQ.'S' .OR. C.EQ.'4')THEN
         ISCAN=4
      ELSEIF(C.EQ.'C' .OR. C.EQ.'5')THEN
         ISCAN=5
      ELSEIF(C.EQ.'H' .OR. C.EQ.'6')THEN
         ISCAN=6
      ELSEIF(C.EQ.'U' .OR. C.EQ.'7')THEN
         ISCAN=7
      ELSE
         ISCAN=-99
      ENDIF
      END

      SUBROUTINE REUSE(IV1,IV2,IDV,LSKIP,IEMSCT,IPUSCR)

!     TO REUSE TAPE7.SCR, WE MUST READ OLD TAPE7.SCR.
!     IV1 AND IV2 INPUTS TO THIS ROUTINE ARE FOR CURRENT CARD 4.
!     BUT THESE MAY NOT MATCH THE TAPE7.SCR OF THE EARLIER RUN.
!     SO RESET IV1 AND IV2 TO TAPE7.SCR VALUES.
!     TAPE7.SCR CONTAINS VALUES FROM IV1 TO IV2 CM-1 IN INCREMENTS OF 1
!     UNLESS IV1 AND IV2 ARE KNOWN TAPE7.SCR CANNOT BE READ PROPERLY.

!     IDV INPUT MUST MATCH WHAT IS IN TAPE7.SCR.
!     IF NOT, RUN MODTRAN AS THE OLD TAPE7.SCR IS NO GOOD.

      CHARACTER HEADER*132
      INTEGER IV1,IV2,IDV,IEMSCT,ITMP1,ITMP2,ITMP3,IPUSCR,ICOUNT
      REAL TMP1,TMP2
      LOGICAL LSKIP

!     FIND WHERE HEADER ENDS IN TAPE7.SCR (UNIT=IPUSCR)
 1    READ(IPUSCR,'(A)')HEADER
      IF(HEADER(3:6).NE.'FREQ')GOTO 1
      IF(IEMSCT.EQ.0)READ(IPUSCR,*)

      READ(IPUSCR,*)TMP1
      ITMP1=INT(TMP1+0.001)
      ICOUNT=1
      IF(IEMSCT.EQ.0)THEN
 10      READ(IPUSCR,*,ERR=40,END=40)TMP2
         ICOUNT=ICOUNT+1
         GOTO 10
      ELSEIF(IEMSCT.EQ.1)THEN
 20      READ(IPUSCR,*,ERR=40,END=40)TMP2
         ICOUNT=ICOUNT+1
         GOTO 20
      ELSEIF(IEMSCT.EQ.2)THEN
 30      READ(IPUSCR,*,ERR=40,END=40)TMP2
         ICOUNT=ICOUNT+1
         GOTO 30
      ELSEIF(IEMSCT.EQ.3)THEN
 35      READ(IPUSCR,*,ERR=40,END=40)TMP2
         ICOUNT=ICOUNT+1
         GOTO 35
      ENDIF
 40   CONTINUE
      REWIND(IPUSCR)
      ITMP2=INT(TMP2+0.001)
      ITMP3=(ITMP2-ITMP1)/(ICOUNT-1)
      IF(IDV.NE.ITMP3)THEN
         LSKIP=.FALSE.
         RETURN
      ELSE
         LSKIP=.TRUE.
         IV1=ITMP1
         IV2=ITMP2
      ENDIF
      END
