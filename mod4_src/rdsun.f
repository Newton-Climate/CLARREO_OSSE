      SUBROUTINE RDSUN(IFWHM,SUNFL1,SUNFL2,SOLCON)

!     THIS ROUTINE READS 'SUNFL1', AND POSSIBLY, 'SUNFL2' FILES.
!     THESE SOLAR IRRADIANCE FILES ARE MERGED.  FINALLY THIS
!     ROUTINE PASSES A TRIANGULAR SLIT WITH IFWHM=|IFWHM|.
!     THE RESULTING DATA IS STORED IN ARRAY SUN OF /SOLAR1/.

!     PARAMETERS:
      INCLUDE 'ERROR.h'
      INTEGER MAXSUN
      PARAMETER(MAXSUN=50000)

!     ARGUMENTS:
!       IFWHM    FULL-WIDTH-AT-HALF-MAXIMUM USED IN TRIANGULAR SLIT
!                FUNCTION [IF NEGATIVE, ABS(IFWHM) IS USED IN THE
!                SLIT FUNCTION AND THE DATA IS OUTPUT TO A FILE].
!       SUNFL1   INPUT DATA FILE NAME.
!       SUNFL2   INPUT DATA FILE NAME.
!       SOLCON   SOLAR CONSTANT (WATTS/M^2) IF > 0.0; SCALE FACTOR IF <
      INTEGER IFWHM
      REAL SOLCON,SCALE
      CHARACTER SUNFL1*(*),SUNFL2*(*)

!     COMMONS:
      INCLUDE 'IFIL.h'
      REAL SUN
      COMMON/SOLAR1/SUN(0:MAXSUN)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /SOLAR1/
      EXTERNAL DEVCBD,SUNBD

!     FUNCTIONS:
!       NUNIT    RETURNS AN UNUSED FILE UNIT NUMBER.
      INTEGER NUNIT

!     LOCAL VARIABLES:
      INTEGER NSUN,IVLO,IVHI,IV,IRES,IRES2,IRESM1,IVOFF,I,J,IVOLD,LINE
      LOGICAL LEXIST,LOPEN
      CHARACTER*11 SUNDAT
      REAL SUNIV,SUM

!     DATA:
!       A        COEFFICIENT FOR LOW FREQUENCY POWER LAW APPROXIMATION.
!       B        EXPONENT FOR LOW FREQUENCY POWER LAW APPROXIMATION.
!       LBLKDT   LOGICAL FLAG, TRUE TO CREATE SUNBD.F FILE.
      REAL A,B
      LOGICAL LBLKDT
      DATA A,B/3.50187E-13,1.93281/,LBLKDT/.FALSE./

!     OPEN INPUT DATA FILE.

      INQUIRE(FILE=SUNFL1,EXIST=LEXIST,OPENED=LOPEN,NUMBER=NSUN)
      IF(.NOT.LEXIST)THEN
          WRITE(IPR,'(/3A)')                                            &
     &      ' ERROR in RDSUN:  File ',SUNFL1,' was not found.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP ' ERROR in RDSUN:  Solar irradiance file was not found.'
      ENDIF
      IF(LOPEN)THEN
          REWIND(NSUN)
      ELSE
          NSUN=NUNIT()
          CALL OPNFL(NSUN,SUNFL1,'OLD')
      ENDIF

!     READ FIRST LINE OF SOLAR DATA
!       IVLO     INITIAL FREQUENCY [CM-1]
!       SUN      SOLAR IRRADIANCE AT IVLO [W CM-2 / CM-1]
      READ(NSUN,*,ERR=100)
      READ(NSUN,*,ERR=100)
      READ(NSUN,*,ERR=100)IVLO,SUNIV
      IF(IVLO.LT.1 .OR. IVLO.GT.MAXSUN)GOTO100
      IVHI=IVLO
      SUN(IVHI)=SUNIV

!     READ REMAINING SOLAR DATA (1 CM-1 INCREMENTS).
   10 CONTINUE
      READ(NSUN,*,END=20,ERR=100)IV,SUNIV
      IF(IV.NE.IVHI+1)GOTO100
      IVHI=IV
      SUN(IVHI)=SUNIV
      IF(IVHI.LT.MAXSUN)GOTO10
   20 CONTINUE
      CLOSE(NSUN)

!     IF REQUESTED, READ ADDITIONAL "SUNFL2" AND MERGE WITH SUN ARRAY:
      IF(SUNFL2.NE.'1' .AND. SUNFL2.NE.'2' .AND.                        &
     &   SUNFL2.NE.'3' .AND. SUNFL2.NE.'4')CALL RDSUN2(SUNFL2,IVHI)

!     SCALE OR RENORMALIZE IF SOLCON (SOLAR CONSTANT) IS NOT ZERO
      IF (SOLCON.LT.0.0) THEN

!        TREAT ABSOLUTE VALUE AS A SCALE FACTOR.
!        NORMALLY IN THIS CASE THE ABSOLUTE VALUE SHOULD BE CLOSE TO UNI
         SCALE=ABS(SOLCON)
      ELSEIF (SOLCON.GT.0.0) THEN

!        HERE SOLCON IS THE SOLAR CONSTANT IN WATTS/M*2.
!        SO IF SOLCON IS POSITIVE, RENORMALIZE SOLAR DATA TO THIS NUMBER
!        FIND THE INTEGRAL OF THE SUN ARRAY BY SUMMING 1-CM-1 INCREMENT

!        FOR REFERENCE, THESE ARE THE SOLAR "CONSTANTS" FOR THE 4 DATAFI
!        1368.0 W/M**2 FOR NEWKUR.DAT
!        1359.8 FOR CHKUR.DAT
!        1362.1 W/M**2 FOR CEBCHKUR.DAT
!        1376.2 FOR THKUR.DAT

         SUM=0.0
         DO 25 I=IVLO,IVHI
            SUM=SUN(I)+SUM
 25      CONTINUE

!        MULTIPLY BY 1.0E4 TO CONVERT FROM W/CM^2 TO WATTS/M^2
         SUM=SUM*1.0E4
         SCALE=SOLCON/SUM
      ELSEIF (SOLCON.EQ.0.) THEN

!        DO NOTHING
         SCALE=1.0
      ENDIF

!     PASS TRIANGULAR SLIT FUNCTION OVER DATA.
      IRES=ABS(IFWHM)
      IF(IRES.GT.1)THEN
         IRES2=IRES*IRES
         IRESM1=IRES-1
         DO 40 IV=IVLO+IRESM1,IVHI-IRESM1
            SUM=IRES*SUN(IV)
            DO 30 IVOFF=1,IRESM1
               SUM=SUM+(IRES-IVOFF)*(SUN(IV+IVOFF)+SUN(IV-IVOFF))
 30         CONTINUE

!           STORE DATA OFFSET BY IRESM1 CM-1 TO
!           AVOID OVERWRITING REQUIRED DATA.
            SUN(IV-IRESM1)=SUM/IRES2
 40      CONTINUE
      ELSE
         IRESM1=0
      ENDIF

!     MOVE DATA TO PROPER LOCATION (FREQUENCY).
!     ALSO, MULTIPLY BY SCALE HERE.
      DO 50 IV=IVHI-IRESM1,IVLO+IRESM1,-1
          SUN(IV)=SUN(IV-IRESM1)*SCALE
   50 CONTINUE

!     DEFINE LOW FREQUENCY DATA [W CM-2 / CM-1] USING
!     POWER LAW APPROXIMATION.
      SUN(0)=0.
      DO 60 IV=1,IVLO+IRESM1-1
          SUN(IV)=A*FLOAT(IV)**B
   60 CONTINUE

!     INSERT A CONSTANT VALUE AT HIGH FREQUENCIES.
      DO 70 IV=IVHI-IRESM1+1,MAXSUN
          SUN(IV)=SUN(IVHI-IRESM1)
   70 CONTINUE

!     THE FOLLOWING CODING WAS USED TO CREATE THE BLOCK DATA SUNBD.
!     THE LOGICAL FLAG LBLKDT IS HARD-WIRED TO FALSE, BUT THE CODING IS
!     LEFT HERE TO AID IN INCORPORATION OF NEW SOLAR IRRADIANCE DATA.
      IF(LBLKDT)THEN

!         OPEN SUNBD.F FILE.
          CALL OPNFL(NSUN,'sunbd.f','UNKNOWN')

!         WRITE HEADER.
          WRITE(NSUN,'((A))')                                           &
     &      '      BLOCK DATA SUNBD',                                   &
     &      ' ',                                                        &
     &      '!     SOLAR IRRADIANCES [W CM-2 / CM-1] TABULATED EACH',   &
     &      '!     WAVENUMBER FROM 0 TO 50,000 CM-1 AT 5 CM-1 SPECTRAL',&
     &      '!     RESOLUTION (IFWHM TRIANGULAR SLIT).  THIS DATA IS',  &
     &      '!     DERIVED FROM THE 1 CM-1 TABULATION (FILE SUN2) OF',  &
     &      '!     KURUCZ HIGH SPECTRAL RESOLUTION CALCULATIONS.',      &
     &      '      INTEGER MAXSUN',                                     &
     &      '      PARAMETER(MAXSUN=50000)',                            &
     &      '      REAL SUN',                                           &
     &      '      COMMON/SOLAR1/SUN(0:MAXSUN)',                        &
     &      '      INTEGER I'

!         WRITE 0 TO 100 CM-1 DATA.
          WRITE(NSUN,'(A,2(/A),25X,1PE10.3,A)')                         &
     &      'C','C         0 TO   100 CM-1 DATA.',                      &
     &      '      DATA (SUN(I),I=     0,   100)/',SUN(0),','
          WRITE(NSUN,'((9(I6,6(1PE10.3,A),/),A,6(1PE10.3,A)))')         &
     &      (J,      (SUN(I),',',I=6*J-5,6*J),J=1,9),                   &
     &      '     &',(SUN(I),',',I=55   ,60 ),                          &
     &      (J-10,   (SUN(I),',',I=6*J-5,6*J),J=11,16),                 &
     &       J-10,   (SUN(I),',',I=97,99),SUN(100),'/'

!         WRITE 101 TO 100*INT(MAXSUN/100) CM-1 DATA.
          IVOLD=100
          DO 80 IV=200,MAXSUN,100
              WRITE(NSUN,'(A,2(/A,2(I6,A)))')                           &
     &          'C','C    ',IVOLD+1,' TO',IV,' CM-1 DATA.',             &
     &          '      DATA (SUN(I),I=',IVOLD+1,',',IV,')/'
              WRITE(NSUN,'((9(I6,6(1PE10.3,A),/),A,6(1PE10.3,A)))')     &
     &          (J      ,(SUN(IVOLD+I),',',I=6*J-5,6*J),J=1,9),         &
     &          '     &',(SUN(IVOLD+I),',',I=55,60),                    &
     &          (J-10   ,(SUN(IVOLD+I),',',I=6*J-5,6*J),J=11,16),       &
     &           J-10   ,(SUN(IVOLD+I),',',I=97,99),SUN(IV),'/'
              IVOLD=IV
   80     CONTINUE

!         WRITE REMAINING CM-1 DATA.
          IF(IVOLD.LT.MAXSUN)THEN
              WRITE(NSUN,'(A,2(/A,2(I6,A)))')                           &
     &          'C','C    ',IVOLD+1,' TO',MAXSUN,' CM-1 DATA.',         &
     &          '      DATA (SUN(I),I=',IVOLD+1,',',MAXSUN,')/'
              LINE=1
   90         CONTINUE
              IF(IVOLD+6.LT.MAXSUN)THEN
                  IF(LINE.LT.10)THEN
                      WRITE(NSUN,'(I6,6(1PE10.3,A))')LINE,              &
     &                  (SUN(IVOLD+I),',',I=1,6)
                      LINE=LINE+1
                  ELSE
                      WRITE(NSUN,'(A,6(1PE10.3,A))')'     &',           &
     &                  (SUN(IVOLD+I),',',I=1,6)
                      LINE=1
                  ENDIF
                  IVOLD=IVOLD+6
                  GOTO90
              ENDIF
              IF(LINE.LT.10)THEN
                  WRITE(NSUN,'(I6,6(1PE10.3,A))')LINE,                  &
     &              (SUN(I),',',I=IVOLD+1,MAXSUN-1),SUN(MAXSUN),'/'
              ELSE
                  WRITE(NSUN,'(A,6(1PE10.3,A))')'     &',               &
     &              (SUN(I),',',I=IVOLD+1,MAXSUN-1),SUN(MAXSUN),'/'
              ENDIF
          ENDIF
          WRITE(NSUN,'(A)')'      END'
          CLOSE(NSUN)
      ENDIF

!     RETURN UNLESS IFWHM IS NEGATIVE
      IF(IFWHM.GE.0)RETURN

!     OPEN OUTPUT DATA FILE.
      IF(SUNFL2.EQ.'1' .OR. SUNFL2.EQ.'2' .OR.                          &
     &   SUNFL2.EQ.'3' .OR. SUNFL2.EQ.'4')THEN
          SUNDAT='sn_0000.dat'
          SUNDAT(2:2)=SUNFL2
      ELSE
          SUNDAT='s6_0000.dat'
      ENDIF
      IF(IRES.LT.10)THEN
          WRITE(SUNDAT(7:7),'(I1)')IRES
      ELSEIF(IRES.LT.100)THEN
          WRITE(SUNDAT(6:7),'(I2)')IRES
      ELSEIF(IRES.LT.1000)THEN
          WRITE(SUNDAT(5:7),'(I3)')IRES
      ELSE
          WRITE(SUNDAT(4:7),'(I4)')IRES
      ENDIF
      CALL OPNFL(NSUN,SUNDAT,'UNKNOWN')

!     WRITE DATA TO SUNDAT IN ORIGINAL UNITS [W CM-2 / CM-1].
      WRITE(NSUN,'(A,/A,/(I7,1PE13.3))')                                &
     &  '  FREQ   SOLAR IRRADIANCE',                                    &
     &  ' (CM-1)  (W CM-2 / CM-1)',(IV,SUN(IV),IV=1,MAXSUN)
      CLOSE(NSUN)

!     RETURN TO DRIVER.
      RETURN

!     WRITE OUT ERROR MESSAGE AND STOP.
  100 CONTINUE
      WRITE(IPR,'(/3A)')                                                &
     &  ' ERROR in RDSUN:  Problem reading/using file ',SUNFL1,' data.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP                                                          &
     &    'ERROR in RDSUN:  Problem reading/using solar irradiance data'
      END

      SUBROUTINE RDSUN2(USRFIL,IVHI)

!     PARAMETERS:
      INCLUDE 'ERROR.h'
      CHARACTER*(*) USRFIL

!     THIS PROGRAM READS USER-SELECTED SOLAR IRRADIANCE DATA

!     THE DATA FILE MUST HAVE THIS FORM (SEE EXAMPLE BELOW):

!     HEADER INFORMATION ON 1ST LINE FOLLOWED BY 1 DATA PAIR PER LINE.
!     HEADER HAS TWO INTEGERS EACH OF WHICH IS 1, 2 OR 3.
!     THE FIRST DESIGNATES UNITS FOR X (WAVELENGTH OR FREQUENCY):
!     1 STANDS FOR X IN CM-1
!     2 STANDS FOR X IN NM
!     3 STANDS FOR X IN MICRONS

!     THE SECOND IS FOR UNITS OF Y (IRRADIANCE):
!     1 IS FOR IRRADIANCE IN (W/CM^2)/CM-1
!     2 FOR IRRADIANCE IN (PHOTONS/SEC)/CM^2/NM
!     3 FOR IRRADIANCE IN MILLIWATTS/M^2/NM WATTS/M^2/MICRON

!     E.G., FOR FREQUENCY IN CM-1 AND IRRADIANCE IN W CM-2 / CM-1,
!     THE USER-CHOSEN FILE SHOULD LOOK LIKE THIS:
!     1   1
!     51    7.453E-10
!     52    7.711E-10
!     53    7.974E-10
!     54    8.243E-10
!     55    8.516E-10
!     ...
!     49982    2.800E-09
!     49983    2.603E-09

      INTEGER MAXSUN,MDATA,IVHI,IUSRFL,NUNIT
      PARAMETER (MAXSUN=50000)
      PARAMETER (MDATA=100000)
      REAL SUN,TEMPX(MDATA),TEMPY(MDATA)
      COMMON/SOLAR1/SUN(0:MAXSUN)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /SOLAR1/
      EXTERNAL SUNBD

!     DECLARE FUNCTIONS

!     H IS IN JOULES.SEC, C IS IN CM PER S.

      REAL H,C

      INTEGER IUNITX,IUNITY,I,J,MAXIM,MINIM,NDATA
      REAL XTEMP,YTEMP,VAL,ENERGY

      DATA H,C/6.6262E-34,2.997925E10/

!     OPEN FILE WITH DATA BUT FIRST GET RID OF SPACES IN FILE NAME
      IUSRFL=NUNIT()
      CALL OPNFL(IUSRFL,USRFIL,'OLD')
      I=0
      READ(IUSRFL,*)IUNITX,IUNITY
   10 READ(IUSRFL,*,END=20,ERR=90)XTEMP,YTEMP
      I=I+1
      TEMPX(I)=XTEMP
      TEMPY(I)=YTEMP
      GO TO 10
   20 CONTINUE
      CLOSE(UNIT=IUSRFL)
      NDATA=I

!     CONVERT UNITS TO THOSE FLAGGED BY 1 (CM-1 AND W/CM+2/CM-1)

      DO 50 I =1, NDATA
         IF (IUNITX.EQ.1) THEN
            GO TO 30
         ELSEIF (IUNITX.EQ.2) THEN
            TEMPX(I)=1.0E7/TEMPX(I)
         ELSEIF (IUNITX.EQ.3) THEN
            TEMPX(I)=10000.0/TEMPX(I)
         ELSE
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP 'BAD UNIT FLAG IN USER-SUPPLIED SOLAR DATA FILE'
         ENDIF

   30    CONTINUE
         IF (IUNITY.EQ.1) THEN
!           IF (IUNITX.EQ.1) SKIP LOOP ALTOGETHER
            IF (IUNITX.EQ.1) GO TO 60
            GO TO 40
         ELSEIF (IUNITY.EQ.2) THEN

!           ENERGY IN JOULES.  NOTE THAT TEMPX(I) IS IN CM-1.
!           TEMPY*ENERGY IS IN WATTS/CM^2
!           MULTIPLY BY 1.0E7/TEMPX(I)**2 TO CONVERT FROM 1/NM TO 1/CM-1

            ENERGY=H*C*TEMPX(I)
            TEMPY(I)=TEMPY(I)*ENERGY*(1.0E7/TEMPX(I)**2)
         ELSEIF (IUNITY.EQ.3) THEN

!           MULTIPLY BY 1.0E-7 TO COVERT FROM MILLIWATTS/M^2 TO W/CM^2
!           MULTIPLY BY 1.0E7/TEMPX(I)**2 TO COVERT FROM 1/NM TO 1/CM-1
!           NOTE THAT TEMPX(I) IS IN CM-1

            TEMPY(I)=TEMPY(I)/(TEMPX(I)*TEMPX(I))
         ELSE
            IF(LJMASS)CALL WRTBUF(FATAL)
            STOP 'BAD UNITS FLAG FOR IRRADIANCE IN SOLAR DATA FILE'
         ENDIF
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE

!     ARRANGE IN INCREASING ORDER IN X WHICH IS NOW IN CM-1

      IF (TEMPX(1).GT.TEMPX(NDATA)) THEN
!        REVERSE
         DO 70 I=1, NDATA/2
            XTEMP=TEMPX(I)
            YTEMP=TEMPY(I)
            J=NDATA-I+1
            TEMPX(I)=TEMPX(J)
            TEMPY(I)=TEMPY(J)
            TEMPX(J)=XTEMP
            TEMPY(J)=YTEMP
   70    CONTINUE
      ENDIF

!     INTERPOLATE FOR INTEGRAL WAVENUMBERS.
!     THE FOLLOWING LOOP WILL ALSO MERGE W/ DEFAULT SOLAR FILE
!     BY OVERWRITING.

      MINIM=INT(TEMPX(1))+1
      MAXIM=INT(TEMPX(NDATA))
      IF (MAXIM.GT.IVHI) IVHI=MAXIM
      IF (MAXIM.GT.MAXSUN) IVHI=MAXSUN
      DO 80 I=MINIM, MIN(MAXIM,MAXSUN), 1
         CALL LINTRP(TEMPX,TEMPY,NDATA,REAL(I),VAL)
         SUN(I)=VAL
   80 CONTINUE
      RETURN
   90 CONTINUE
      IF(LJMASS)CALL WRTBUF(FATAL)
      STOP 'ERROR IN USER-SUPPLIED SOLAR FILE'
      END

      SUBROUTINE LINTRP(XA,YA,N,X,VAL)
      INTEGER N,J1,JHI
      REAL XA(N),YA(N),VAL,X,Y1,Y2,X1,X2,T

!     INPUTS ARE ARRAYS:  XA(N), YA(N) - YA IS THE FUNCTION VALUE AT XA.
!     FIND WHICH INDICES OF XA SANDWICH X:  XA(J1) .LE. X .LE. XA(J1+1).
!     NOTE THAT X IS THE VALUE ALONG THE XA LINE - IT IS NOT THE INDEX.
!     FIND VAL BY INTERPOLATING LINEARLY.

      CALL HUNT(XA,N,X,J1)
      JHI=J1+1
      Y1=YA(J1)
      Y2=YA(JHI)
      X1=XA(J1)
      X2=XA(JHI)

!     LINEAR INTERPOLATION BETWEEN (X1,Y1) AND (X2,Y2); FIND VALUE @ X.

      IF (Y1.EQ.Y2) THEN
         VAL=Y1
         RETURN
      ENDIF
      T=(X-X1)/(X2-X1)
      VAL=Y1+T*(Y2-Y1)
      END

      SUBROUTINE HUNT(XA,N,X,JLO)

!     MODIFIED FROM NUMERICAL RECIPES' ROUTINE HUNT.
!     LIKE LOCATE (FROM NR) BUT FASTER IF JLO IS A GOOD GUESS.

!     XA(1), XA(2), XA(3), ..., XA(N) ARE THE MONOTONIC GRID PTS.
!     RETURNS JLO=J, IF X IS IN BETWEEN XA(J) & XA(J+1).
!     IF X=XA(J), RETURNS JLO=J-1.
!     IF OUT OF RANGE, RETURNS 0 OR N.

      INTEGER N,JLO,JHI,INC,JM
      REAL XA(N),X
      LOGICAL ASCEND

      IF (X.EQ.XA(1)) THEN
         JLO=1
         RETURN
      ENDIF
      IF (X.EQ.XA(N)) THEN
         JLO=N-1
         RETURN
      ENDIF

      ASCEND=XA(N).GT.XA(1)
      IF(JLO.LE.0.OR.JLO.GT.N)THEN
         JLO=0
         JHI=N+1
         GO TO 3
      ENDIF
      INC=1
      IF(X.GE.XA(JLO).EQV.ASCEND)THEN
 1       JHI=JLO+INC
         IF(JHI.GT.N)THEN
            JHI=N+1
         ELSE IF(X.GE.XA(JHI).EQV.ASCEND)THEN
            JLO=JHI
            INC=INC+INC
            GO TO 1
         ENDIF
      ELSE
         JHI=JLO
 2       JLO=JHI-INC
         IF(JLO.LT.1)THEN
            JLO=0
         ELSE IF(X.LT.XA(JLO).EQV.ASCEND)THEN
            JHI=JLO
            INC=INC+INC
            GO TO 2
         ENDIF
      ENDIF

 3    IF(JHI-JLO.EQ.1) RETURN

      JM=(JHI+JLO)/2
      IF(X.GT.XA(JM).EQV.ASCEND)THEN
         JLO=JM
      ELSE
         JHI=JM
      ENDIF
      GO TO 3
      END
