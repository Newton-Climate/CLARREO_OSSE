      SUBROUTINE RDCORK(CKNAME,NSPEED,NTLSUB,KNTRVL)

!     ROUTINE TO INITIALIZE CORRELATED-K DISTRIBUTION DATA.

!     DECLARE ROUTINE ARGUMENT INPUTS.
!       CKNAME   NAME OF CORRELATED-K DISTRIBUTIONS FILE.
!       NSPEED   COMPUTATIONAL SPEED FLAG.
!       NTLSUB   NUMBER OF LINE TAIL PARAMETERS PER BAND.
      CHARACTER CKNAME*(*)
      INTEGER NSPEED,NTLSUB

!     DECLARE ROUTINE ARGUMENT OUTPUTS.
!       KNTRVL   NUMBER OF K-INTERVAL TO BE USED.
      INTEGER KNTRVL

!     PARAMETERS:
!       MXKSUB   DIMENSION OF K-DISTRIBUTION SUB-INTERVAL ARRAY.
!       MXGAML   DIMENSION OF LORENTZ HALF-WIDTH ARRAY.
!       MXGAMD   DIMENSION OF DOPPLER HALF-WIDTH ARRAY.
!       MXNUML   DIMENSION OF EFFECTIVE NUMBER OF LINES ARRAY.
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:

!     COMMON /CORKDT/
!       WTKSUB   SPECTRAL BIN SUB-INTERVAL FRACTIONAL WIDTHS.
!       DEPLAY   INCREMENTAL OPTICAL DEPTHS.
!       TRNLAY   INCREMENTAL TRANSMITTANCES.
!       TRNCUM   CUMULATIVE TRANSMITTANCES.
!       K2TAIL   POINTER FROM K BIN TO LINE TAIL SUB-BIN
!                (=0 IF MULTIPLE LINE TAIL SUB-BINS CONTRIBUTE).
!       CONTWT   WEIGHTS FOR PARTITIONING LINE TAILS INTO K'S
!                (ONLY USED IF K2TAIL IS 0).
      REAL WTKSUB,DEPLAY,TRNLAY,TRNCUM,CONTWT
      INTEGER K2TAIL
      COMMON/CORKDT/WTKSUB(MXKSUB),DEPLAY(MXKSUB),TRNLAY(MXKSUB),       &
     &  TRNCUM(MXKSUB),K2TAIL(MXKSUB),CONTWT(MTLSUB,MXKSUB)
      SAVE /CORKDT/

!     COMMON /CORKTB/
!       GAMLIN   LORENTZ HALF-WIDTH BOUNDARY VALUES [CM-1].
!       GAMDIN   DOPPLER HALF-WIDTH BOUNDARY VALUES [CM-1].
!       RLININ   EFFECTIVE NUMBER OF LINES BOUNDARY VALUES.
!       VAL      ABSORPTION COEFFICIENT TABLE [ATM-1 CM-1].
      INTEGER NGAML,NGAMD,NNUML
      REAL GAMLIN,GAMDIN,RLININ,VAL
      COMMON/CORKTB/NGAML,NGAMD,NNUML,GAMLIN(MXGAML),GAMDIN(MXGAMD),    &
     &  RLININ(MXNUML),VAL(1:MXGAML,1:MXGAMD,1:MXNUML,0:MXKSUB)
      SAVE /CORKTB/
      INCLUDE 'IFIL.h'

!     BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
      INTEGER IFILE,IKSUB,IKSBP1,IGAML,IGAMD,INUML,NKSUB,INTRVL,        &
     &  ITLSUB,NDIVM1,ILAST,INEXT
      LOGICAL LTEST
      REAL WEIGHT,TAILWT(MTLSUB)
!old  REAL SUMWTK,CUMWTK(MXKSUB)

!     DATA:
!       LSKIP    LIST OF K'S TO BE SKIPPED IF NSPEED IS 1 OR 2.
      LOGICAL LSKIP(MXKSUB,2)
      DATA (LSKIP(IKSUB,1),IKSUB=1,MXKSUB)     /.FALSE.,.FALSE.,.FALSE.,&
     &  .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,&
     &  .FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,&
     &  .FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,        &
     &  .TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.FALSE./
      DATA (LSKIP(IKSUB,2),IKSUB=1,MXKSUB)     /.TRUE. ,.FALSE.,.TRUE. ,&
     &  .FALSE.,.TRUE. ,.FALSE.,.FALSE.,.TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,&
     &  .FALSE.,.TRUE. ,.TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,.FALSE.,.TRUE. ,&
     &  .TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,.FALSE.,.TRUE. ,.TRUE. ,        &
     &  .TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.TRUE. ,.FALSE./

!     FIND AVAILABLE UNIT.
      DO IFILE=11,99
          INQUIRE(IFILE,OPENED=LTEST)
          IF(.NOT.LTEST)GOTO10
      ENDDO
      WRITE(IPR,'(/A)')                                                 &
     &  ' ERROR in routine RDCORK:  No unit numbers available?'
      IF(LJMASS)CALL WRTBUF(FATAL)
      STOP ' ERROR in routine RDCORK:  No unit numbers available?'

!     TEST FOR EXISTENCE OF FILE.
   10 CONTINUE
      INQUIRE(FILE=CKNAME,EXIST=LTEST)
      IF(.NOT.LTEST)THEN
          WRITE(IPR,'(/2A,/(11X,A))')                                   &
     &      ' WARNING:  The Correlated-k option was turned',            &
     &      ' off because file ',CKNAME,'was not found.'
          RETURN
      ENDIF

!     OPEN CORRELATED-K DATA FILE.
      OPEN(IFILE,FILE=CKNAME,STATUS='OLD',FORM='UNFORMATTED')

!     READ DATA FROM BINARY FILE.
      READ(IFILE)NKSUB,NGAML,NGAMD,NNUML

!     CHECK ARRAY DIMENSIONS.
      IF(NKSUB.GT.MXKSUB .OR. NGAML.GT.MXGAML .OR.                      &
     &  NGAMD.GT.MXGAMD .OR. NNUML.GT.MXNUML)THEN
          WRITE(IPR,'(/2A,/2(10X,2A),/(10X,2(A,I3)))')' WARNING: ',     &
     &      ' The Correlated-k option was turned off.',                 &
     &      ' The Correlated-k distributions file ',CKNAME,             &
     &      ' requires minimally the following assignments',            &
     &      ' in file PARAMS.h.',                                       &
     &      '   MXKSUB =',NKSUB,'        MXGAML =',NGAML,               &
     &      '   MXGAMD =',NGAMD,'        MXNUML =',NNUML
          RETURN
      ENDIF
      READ(IFILE)                                                       &
     &  (WTKSUB(IKSUB),IKSUB=1,NKSUB),(GAMLIN(IGAML),IGAML=1,NGAML),    &
     &  (GAMDIN(IGAMD),IGAMD=1,NGAMD),(RLININ(INUML),INUML=1,NNUML)
      DO 30 IKSUB=0,NKSUB
          READ(IFILE)(((VAL(IGAML,IGAMD,INUML,IKSUB),                   &
     &      IGAML=1,NGAML),IGAMD=1,NGAMD),INUML=1,NNUML)
   30 CONTINUE

!     CLOSE CORRELATED-K DATA FILE.
      CLOSE(IFILE)
      IF(.NOT.LJMASS) WRITE(IPR,'(2A)')                                 &
     &  ' SUCCESSFULLY READ CORRELATED-K DISTRIBUTIONS FILE:  ',CKNAME

!     TAILOR SELECTION OF K'S.
      IF(NSPEED.LE.0 .OR. NSPEED.GT.2)THEN

!         FOR HIGHEST RESOLUTION, USE ALL THE K'S.
          KNTRVL=NKSUB
      ELSE

!         THE FIRST AND LAST K'S CANNOT BE SKIPPED.
          KNTRVL=0
          IKSUB=1
          LSKIP(NKSUB,NSPEED)=.FALSE.
          DO 70 IKSBP1=2,NKSUB+1
              IF(LSKIP(IKSUB,NSPEED))THEN

!                 SKIP THE "IKSUB" K'S BY ADDING WEIGHT TO "IKSUB+1".
                  WTKSUB(IKSBP1)=WTKSUB(IKSBP1)+WTKSUB(IKSUB)
              ELSE

!                 INCLUDE THE "IKSUB" K'S.
                  KNTRVL=KNTRVL+1
                  IF(KNTRVL.NE.IKSUB)THEN

!                     MOVE WEIGHTS AND K'S TO NEW LOCATION.
                      WTKSUB(KNTRVL)=WTKSUB(IKSUB)
                      DO 60 INUML=1,NNUML
                          DO 50 IGAMD=1,NGAMD
                              DO 40 IGAML=1,NGAML
                                  VAL(IGAML,IGAMD,INUML,KNTRVL)         &
     &                              =VAL(IGAML,IGAMD,INUML,IKSUB)
   40                         CONTINUE
   50                     CONTINUE
   60                 CONTINUE
                  ENDIF
              ENDIF
   70     IKSUB=IKSBP1
      ENDIF

!     INITIALIZE TAIL WEIGHTS:
      DO 80 ITLSUB=1,NTLSUB
          TAILWT(ITLSUB)=1./NTLSUB
   80 CONTINUE

!     FIRST DETERMINE WTKSUB FITTING IN A SINGLE INTERVAL:
      NDIVM1=NTLSUB-1
      ILAST=NTLSUB
      DO 110 INTRVL=1,KNTRVL

!         INITIALIZE K2TAIL AND CONTWT ARRAYS:
          K2TAIL(INTRVL)=0
          DO 90 ITLSUB=1,NTLSUB
              CONTWT(ITLSUB,INTRVL)=0.
   90     CONTINUE

!         CHECK IF WTKSUB FITS IN A SINGLE INTERVAL:
          DO 100 ITLSUB=ILAST,ILAST+NDIVM1
              INEXT=MOD(ITLSUB,NTLSUB)+1
              IF(WTKSUB(INTRVL).LE.TAILWT(INEXT))THEN

!                 K SUB-INTERVAL INTRVL WILL USE LINE TAIL INEXT:
                  TAILWT(INEXT)=TAILWT(INEXT)-WTKSUB(INTRVL)
                  K2TAIL(INTRVL)=INEXT
                  CONTWT(INEXT,INTRVL)=1.
                  ILAST=INEXT
                  GOTO110
              ENDIF
  100     CONTINUE
  110 CONTINUE

!     LOOP OVER K SUB-INTERVAL REQUIRING PARTITIONING OF WTKSUB:
!old  SUMWTK=0.
      DO 130 INTRVL=1,KNTRVL
!old      SUMWTK=SUMWTK+WTKSUB(INTRVL)
!old      CUMWTK(INTRVL)=SUMWTK
          IF(K2TAIL(INTRVL).GT.0)GOTO130
          WEIGHT=WTKSUB(INTRVL)
          DO 120 ITLSUB=1,NTLSUB
              IF(WEIGHT.GT.TAILWT(ITLSUB))THEN

!                 USE REMAINDER OF TAILWT(ITLSUB):
                  CONTWT(ITLSUB,INTRVL)=TAILWT(ITLSUB)/WTKSUB(INTRVL)
                  WEIGHT=WEIGHT-TAILWT(ITLSUB)
                  TAILWT(ITLSUB)=0.
              ELSE

!                 REMAINDER OF WEIGHT FITS IN ITLSUB
                  CONTWT(ITLSUB,INTRVL)=WEIGHT/WTKSUB(INTRVL)
                  TAILWT(ITLSUB)=TAILWT(ITLSUB)-WEIGHT
                  GOTO130
              ENDIF
  120     CONTINUE
  130 CONTINUE

!     WRITE WTKSUB INTO HEADER OF UNIT=IDBOUT FILE:
!old  WRITE(IDBOUT,'(4A,/90X,I2,33(1X,F13.7))')'    P[ATM]',
!old 1  '   AEROSOL_SCT   AEROSOL_EXT   AER_G     CLOUD_SCT',
!old 2  '     CLOUD_EXT   CLD_G      RAIN_SCT      RAIN_EXT',
!old 3  '  RAIN_G  RAYLEIGH_SCT  MOLECULAR_ABSORPTION...',
!old 4  KNTRVL,(CUMWTK(INTRVL),INTRVL=1,KNTRVL)

!     RETURN TO DRIVER.
      RETURN
      END
