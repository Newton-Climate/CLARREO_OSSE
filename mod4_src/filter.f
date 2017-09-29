      LOGICAL FUNCTION FILTER(FILTNM,LNFILT)

!     THIS LOGICAL FUNCTION READS THE NAME OF A FILTER
!     FUNCTION FILE, OPENS THE FILE, AND PROCESSES THE FILTER
!     FOR EACH BAND TO ENABLE RAPID PROCESSING.  IF AN ERROR
!     OCCURS, FILTER IS RETURNED WITH A VALUE OF .FALSE.

!     DECLARE ARGUMENTS:
!       FILTNM   FILTER FUNCTION FILE NAME.
!       LNFILT   LENGTH OF FILTER FUNCTION FILE NAME.
      CHARACTER*(*) FILTNM
      INTEGER LNFILT

!     PARAMETERS:
      INCLUDE 'PARAMS.h'
      INCLUDE 'CHANLS.h'
      INCLUDE 'ERROR.h'
      INCLUDE 'IFIL.h'
      INCLUDE 'BMHEAD.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:
      INTEGER LENSTR,NUNIT
      LOGICAL WTCHAN

!     LOCAL VARIABLES:
!       IUNIT    UNIT NUMBER FOR A TEMPORARY FILE.
!       CONVRT   WAVELENGTH TO CM-1 CONVERSION FACTOR OVER THE
!                CALCULATION BIN WIDTH [NM OR MICRONS].
      CHARACTER*80 INLINE
      INTEGER IOTEST,IFILTR,LNTEST,LSTBIN,NSPEC,IBIN,NBNDWD,IUNIT
      LOGICAL LOPEN,LEXIST
      REAL SPECMN,WGTMN,SPECMX,WGTMX,CONVRT

!     DATA:
!       DATTST   CHARACTER STRING CONTAINING POSSIBLE INITIAL NON-BLANK
!                CHARACTERS FOR FREQUENCY/WAVELENGTH AND WEIGHT LINES.
      CHARACTER DATTST*11,FLTNMS*(NAMLEN)
      LOGICAL NOT0LO,NOT0HI
      SAVE NOT0LO,NOT0HI,NBNDWD,FLTNMS
      DATA DATTST/'.0123456789'/,NOT0LO,NOT0HI/2*.TRUE./,NBNDWD/0/,     &
     &  FLTNMS/' '/

!     INITIALIZE FILTER TO .FALSE.:
      FILTER=.FALSE.

!     CHECK IF FILTER FILE HAS ALREADY BEEN PROCESSED
!     WITH CURRENT SPECTRAL FREQUENCY STEP SIZE.
      IF(.NOT.LJMASS .AND. FILTNM.EQ.FLTNMS .AND. IBNDWD.EQ.NBNDWD)THEN
          WRITE(IPR,'(2(/2A),I4,A)')' FILTER FUNCTION FILE  ',          &
     &      FILTNM(1:LNFILT),' HAS ALREADY BEEN PROCESSED',             &
     &      ' WITH FREQUENCY STEP SIZE',IBNDWD,' CM-1.'
          RETURN
      ENDIF
      FLTNMS=' '

!     CHECK FILTER FILE STATUS:
      IF(LNFILT.LE.0)THEN
          WRITE(IPR,'(/A)')                                             &
     &      ' WARNING:  The input filter file name is BLANK.'
          RETURN
      ENDIF
      INQUIRE(FILE=FILTNM(1:LNFILT),EXIST=LEXIST,OPENED=LOPEN)
      IF(.NOT.LEXIST)THEN
          WRITE(IPR,'(/A,/(11X,A))')' WARNING:  The filter file',       &
     &      FILTNM(1:LNFILT),'does not exist.'
          LNFILT=0
          RETURN
      ELSEIF(LOPEN)THEN
          WRITE(IPR,'(/A,/(11X,A))')' WARNING:  The filter file',       &
     &      FILTNM(1:LNFILT),'is already opened.'
          LNFILT=0
          RETURN
      ENDIF

!     OPEN FILTER FUNCTION:
      IFILTR=NUNIT()
      OPEN(IFILTR,FILE=FILTNM(1:LNFILT),STATUS='OLD')
      IF(.NOT.LJMASS) WRITE(IPR,'(/2A)')                                &
     &   ' OPENED FILTER FILE:  ',FILTNM(1:LNFILT)

!     READ UNIT FLAG:
      READ(IFILTR,'(A80)',IOSTAT=IOTEST)INLINE
      IF(IOTEST.GT.0)THEN
          WRITE(IPR,'(/A)')                                             &
     &      ' WARNING:  Error reading FILTER FILE UNIT FLAG.'
          CLOSE(IFILTR,STATUS='KEEP')
          LNFILT=0
          RETURN
      ELSEIF(IOTEST.LT.0)THEN
          WRITE(IPR,'(/2A)')' WARNING:  End-Of-File',                   &
     &      ' occurred during read of FILTER FILE UNIT FLAG.'
          CLOSE(IFILTR,STATUS='KEEP')
          LNFILT=0
          RETURN
      ENDIF
      LNTEST=LENSTR(INLINE)
      IF(INLINE(1:1).EQ.'M' .OR. INLINE(1:1).EQ.'m')THEN

!         MICRONS:
          UNTFLG='M'
          CONVRT=1.E4/IBNDWD
      ELSEIF(INLINE(1:1).EQ.'N' .OR. INLINE(1:1).EQ.'n')THEN

!         NANOMETERS:
          UNTFLG='N'
          CONVRT=1.E7/IBNDWD
      ELSE

!         WAVENUMBERS:
          IF(INLINE(1:1).NE.'W' .AND. INLINE(1:1).NE.'w')               &
     &      WRITE(IPR,'(/(A))')                                         &
     &        ' WARNING:  Filter file unit flag does not begin with',   &
     &        '           "W", but WAVENUMBERS is being assumed.'
          UNTFLG='W'
      ENDIF

!     INITIALIZE NUMCHN ARRAY TO ZERO.
      DO 10 IBIN=MNBIN,MXBIN
          NUMCHN(IBIN)=0
   10 CONTINUE

!     READ TITLE OF FIRST CHANNEL AND THE FIRST TWO FREQUENCY/WAVELENGTH
!     AND WEIGHT PAIRS, AND THEN PROCESS THOSE WEIGHTS.
      NCHAN=1
      READ(IFILTR,'(A80)',IOSTAT=IOTEST)NMCHAN(1)
      IF(IOTEST.GT.0)THEN
          WRITE(IPR,'(/A)')' WARNING:  Error reading CHANNEL   1 NAME'
          CALL RMCHAN(IFILTR)
          RETURN
      ELSEIF(IOTEST.LT.0)THEN
          WRITE(IPR,'(/2A)')' WARNING:  End-Of-File',                   &
     &      ' occurred during read of CHANNEL   1 NAME'
          CALL RMCHAN(IFILTR)
          RETURN
      ENDIF
      LNCHAN(1)=LENSTR(NMCHAN(1))
      READ(IFILTR,*,IOSTAT=IOTEST)SPECMN,WGTMN
      IF(IOTEST.GT.0)THEN
          WRITE(IPR,'(/2A)')' WARNING:  Error reading data',            &
     &      ' pair   1 for channel   1 in filter function file.'
          CALL RMCHAN(IFILTR)
          RETURN
      ELSEIF(IOTEST.LT.0)THEN
          WRITE(IPR,'(/2A)')' WARNING:  End-Of-File - data',            &
     &      ' pair   1 for channel   1 in filter function file.'
          CALL RMCHAN(IFILTR)
          RETURN
      ELSEIF(WGTMN.NE.0. .AND. NOT0LO)THEN
          WRITE(IPR,'(/2A,/(11X,2A))')' WARNING: ',                     &
     &      ' First filter weight for channel   1,',                    &
     &      NMCHAN(1)(1:LNCHAN(1)),',','is non-zero.'
          NOT0LO=.FALSE.
      ENDIF
      READ(IFILTR,*,IOSTAT=IOTEST)SPECMX,WGTMX
      IF(IOTEST.GT.0)THEN
          WRITE(IPR,'(/2A)')' WARNING:  Error reading data',            &
     &      ' pair   2 for channel   1 in filter function file.'
          CALL RMCHAN(IFILTR)
          RETURN
      ELSEIF(IOTEST.LT.0)THEN
          WRITE(IPR,'(/2A)')' WARNING:  End-Of-File - data',            &
     &      ' pair   2 for channel   1 in filter function file.'
          CALL RMCHAN(IFILTR)
          RETURN
      ENDIF
      NSPEC=2
      LSTBIN=-1
      WDFREQ(1)=0.
      WDWAVE(1)=0.
      IF(.NOT.WTCHAN(LSTBIN,SPECMN,WGTMN,SPECMX,WGTMX))THEN
          WRITE(IPR,'(/(A))')                                           &
     &      ' WARNING:  Routine WTCHAN error for pairs   1 and   2',    &
     &      '           of channel   1 in filter function file.'
          CALL RMCHAN(IFILTR)
          RETURN
      ENDIF
      VLCHAN(1)=0.
      SPECLO(1)=SPECMN
      IF(UNTFLG.EQ.'W')THEN
          NFRQLO(1)=IBNDWD*INT(SPECMN/IBNDWD+.5)
      ELSE
          NFRQHI(1)=IBNDWD*INT(CONVRT/SPECMN+.5)
      ENDIF

!     READ AND PARSE NEXT LINE OF FILTER FILE.  THE FIRST NON-BLANK
!     CHARACTER OF A LINE CONTAINING A FREQUENCY/WAVELENGTH AND WEIGHT
!     MUST BE A "-", A ".", OR A NUMERIC CHARACTER (0, 1, ..., 9).
   20 CONTINUE
          READ(IFILTR,'(A80)',IOSTAT=IOTEST)INLINE
          IF(IOTEST.GT.0)THEN

!             ERROR READING FILTER RESPONSE FUNCTION FILE.
              WRITE(IPR,'(/2A,/(11X,A,I4,A))')' WARNING:  Error',       &
     &          ' reading filter function file;  the last line',        &
     &          'successfully read was from CHANNEL',NCHAN,             &
     &          ' entitled',NMCHAN(NCHAN)(1:LNCHAN(NCHAN))
              CALL RMCHAN(IFILTR)
              RETURN
          ELSEIF(IOTEST.LT.0)THEN

!             END-OF-FILE - NORMAL TERMINATION.
              IF(WGTMX.NE.0. .AND. NOT0HI)THEN
                  WRITE(IPR,'(/2A,I4,A,/(11X,2A))')' WARNING: ',        &
     &              ' Last filter weight for channel',NCHAN,',',        &
     &              NMCHAN(NCHAN)(1:LNCHAN(NCHAN)),',','is non-zero.'
                  NOT0HI=.FALSE.
              ENDIF
              SPECHI(NCHAN)=SPECMX
              IF(UNTFLG.EQ.'W')THEN
                  NFRQHI(NCHAN)=IBNDWD*INT(SPECMX/IBNDWD+.5)
              ELSE
                  NFRQLO(NCHAN)=IBNDWD*INT(CONVRT/SPECMX+.5)
              ENDIF
              FILTER=.TRUE.
              FLTNMS=FILTNM
              NBNDWD=IBNDWD
              IF(.NOT.LJMASS)WRITE(IPR,'(/I4,(2A))')NCHAN,              &
     &           ' channels were read',                                 &
     &           ' from the filter function file,',FILTNM(1:LNFILT)
              CLOSE(IFILTR,STATUS='KEEP')
              RETURN
          ENDIF
          LNTEST=LENSTR(INLINE)
          IF(INDEX(DATTST,INLINE(1:1)).GT.0)THEN

!             INLINE CONTAINS FREQUENCY/WAVELENGTH AND WEIGHT.
              SPECMN=SPECMX
              WGTMN=WGTMX
              NSPEC=NSPEC+1

!             SINCE READING FROM A CHARACTER STRING USING A FREE FORMAT,
!             I.E.  "READ(INLINE,FMT='(*)',IOSTAT=IOTEST)SPECMX,WGTMX"
!             IS NOT ASCII STANDARD AND NOT PERMITTED FOR SOME FORTRAN
!             SOME COMPILERS, INLINE IS WRITTEN TO A FILE AND THEN READ.
! JMASS: CREATE COMMON BLOCK FOR SPECMX,WGTMX??
              IUNIT=NUNIT()
              OPEN(IUNIT,FILE='ONELINE',                                &
     &          FORM='FORMATTED',STATUS='UNKNOWN')
              WRITE(IUNIT,'(A80)')INLINE
              REWIND(IUNIT)
              READ(IUNIT,*,IOSTAT=IOTEST)SPECMX,WGTMX
              CLOSE(IUNIT,STATUS='DELETE')
              IF(IOTEST.GT.0)THEN
                  WRITE(IPR,'(/3(A,I4))')                               &
     &              ' WARNING:  Error reading data pair',NSPEC,         &
     &              ' for channel',NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ELSEIF(IOTEST.LT.0)THEN
                  WRITE(IPR,'(/3(A,I4))')                               &
     &              ' WARNING:  End-Of-File - data pair',NSPEC,         &
     &              ' for channel',NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ENDIF
          ELSEIF(NCHAN.LT.MXCHAN)THEN

!             SAVE PREVIOUS CHANNEL MAXIMUM SPECTRAL POINT,
!             INCREMENT CHANNEL, ASSIGN TITLE, AND READ AND CHECK
!             THE FIRST TWO FREQUENCY/WAVELENGTH AND WEIGHT PAIRS.
              IF(WGTMX.NE.0. .AND. NOT0HI)THEN
                  WRITE(IPR,'(/2A,I4,A,/(11X,2A))')' WARNING: ',        &
     &              ' Last filter weight for channel',NCHAN,',',        &
     &              NMCHAN(NCHAN)(1:LNCHAN(NCHAN)),',','is non-zero.'
                  NOT0HI=.FALSE.
              ENDIF
              SPECHI(NCHAN)=SPECMX
              IF(UNTFLG.EQ.'W')THEN
                  NFRQHI(NCHAN)=IBNDWD*INT(SPECMX/IBNDWD+.5)
              ELSE
                  NFRQLO(NCHAN)=IBNDWD*INT(CONVRT/SPECMX+.5)
              ENDIF
              NCHAN=NCHAN+1
              NMCHAN(NCHAN)=INLINE
              LNCHAN(NCHAN)=LNTEST
              READ(IFILTR,*,IOSTAT=IOTEST)SPECMN,WGTMN
              IF(IOTEST.GT.0)THEN
                  WRITE(IPR,'(/2A,I4,A)')' WARNING: ',                  &
     &              ' Error reading data pair   1 for channel',         &
     &              NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ELSEIF(IOTEST.LT.0)THEN
                  WRITE(IPR,'(/2A,I4,A)')' WARNING: ',                  &
     &              ' End-Of-File - data pair   1 for channel',         &
     &              NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ELSEIF(WGTMN.NE.0. .AND. NOT0LO)THEN
                  WRITE(IPR,'(/2A,I4,A,/(11X,2A))')' WARNING: ',        &
     &              ' First filter weight for channel',NCHAN,',',       &
     &              NMCHAN(NCHAN)(1:LNCHAN(NCHAN)),',','is non-zero.'
                  NOT0LO=.FALSE.
              ENDIF
              READ(IFILTR,*,IOSTAT=IOTEST)SPECMX,WGTMX
              IF(IOTEST.GT.0)THEN
                  WRITE(IPR,'(/2A,I4,A)')' WARNING: ',                  &
     &              ' Error reading data pair   2 for channel',         &
     &              NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ELSEIF(IOTEST.LT.0)THEN
                  WRITE(IPR,'(/2A,I4,A)')' WARNING: ',                  &
     &              ' End-Of-File - data pair   2 for channel',         &
     &              NCHAN,' in filter function file.'
                  CALL RMCHAN(IFILTR)
                  RETURN
              ENDIF
              NSPEC=2
              LSTBIN=-1
              SPECLO(NCHAN)=SPECMN
              IF(UNTFLG.EQ.'W')THEN
                  NFRQLO(NCHAN)=IBNDWD*INT(SPECMN/IBNDWD+.5)
              ELSE
                  NFRQHI(NCHAN)=IBNDWD*INT(CONVRT/SPECMN+.5)
              ENDIF
              WDFREQ(NCHAN)=0.
              WDWAVE(NCHAN)=0.
              VLCHAN(NCHAN)=0.
          ELSE
              WRITE(IPR,'(/A,I4,A,/(10X,A))')                           &
     &          ' WARNING:  Only the first',NCHAN,                      &
     &          ' channels in the filter function',                     &
     &          ' file are being used; increase parameter MXCHAN',      &
     &          ' to accommodate all the channels in the file.'
              RETURN
          ENDIF

!         PROCESS DATA AND THEN RETURN TO 20 TO READ NEXT LINE.
          IF(.NOT.WTCHAN(LSTBIN,SPECMN,WGTMN,SPECMX,WGTMX))THEN
              WRITE(IPR,'(/A,2(A,I4),/10X,A,I4,A)')' WARNING: ',        &
     &          ' Routine WTCHAN error for pairs',NSPEC-1,' and',NSPEC, &
     &          ' of channel',NCHAN,' in filter function file.'
              CALL RMCHAN(IFILTR)
              RETURN
          ENDIF
          GOTO20
      END
