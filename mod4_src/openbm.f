      SUBROUTINE OPENBM(BMNAME,DATDIR,LNFLRT,FLRT)

!     OPENBM OPENS THE BAND MODEL DATA FILES, AND READS HEADER DATA.

!     PARAMETERS:
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     ARGUMENTS:
!       BMNAME   NAME OF BAND MODEL PARAMETER DATA FILE.
!       DATDIR   NAME OF MODTRAN DATA DIRECTORY.
      CHARACTER BMNAME*(NAMLEN),DATDIR*(*)
      INTEGER LNFLRT
      CHARACTER FLRT*(*)

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BMHEAD.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:
!       IRECLN   RETURNS THE RECORD LENGTH FOR A DIRECT ACCESS FILES
!                CONTAINING A KNOWN NUMBER OF VARIABLES IN EACH RECORD.
      INTEGER IRECLN,LENSTR

!     LOCAL VARIABLES:
!       LOPEN    FILE OPENED LOGICAL.
!       LEXIST   FILE EXISTS LOGICAL.
!       LRECLN   RECORD LENGTH OF DIRECT ACCESS FILE.
!       ITEMP    TEMPERATURE GRID INDEX
!       MNFREQ   MINIMUM FREQUENCY OF BAND MODEL DATA.
!       FULNAM   FULL PATH FILE NAME.
      LOGICAL LOPEN,LEXIST
      INTEGER LRECLN,ITEMP,MNFREQ
      CHARACTER*(NAMLEN) FULNAM

!     DATA:
!       EDGEMN   MINIMUM ACCEPTABLE VALUE FOR DEDGE [CM-1].
      REAL EDGEMN
      DATA EDGEMN/.2/

!     CLOSE FILE IF UNIT ALREADY CONNECTED:
      INQUIRE(ITB,OPENED=LOPEN)
      IF(LOPEN)CLOSE(ITB,STATUS='KEEP')

!     CHECK THAT FILE EXISTS:
      IF(LENSTR(BMNAME).EQ.0)BMNAME=DATDIR//'B2001_01.BIN'
      INQUIRE(FILE=BMNAME,EXIST=LEXIST)
      IF(.NOT.LEXIST)THEN
          WRITE(IPR,'(/3A,/18X,A)')' ERROR in OPENBM:  The molecular',  &
     &      ' band model data file ',BMNAME,' does not exist.'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'Error:  Molecular band model data file not found.'
      ENDIF

!     OPEN FILE:
      IF(.NOT.LJMASS)WRITE(IPR,'(/A,15X,A)')                            &
     &  ' MOLECULAR BAND MODEL DATA FILE:',BMNAME
      LRECLN=IRECLN(15,LNFLRT,FLRT)
      OPEN(ITB,FILE=BMNAME,STATUS='OLD',ACCESS='DIRECT',                &
     &  FORM='UNFORMATTED',RECL=LRECLN)

!     READ HEADER:
      READ(ITB,REC=1)NTEMP,(TBAND(ITEMP),ITEMP=1,NTEMP),                &
     &  IBNDWD,MNFREQ,MXFREQ,LSTREC,NTLSUB,EDGENR
      IF(NTLSUB.NE.1 .AND. NTLSUB.NE.4)THEN

!         OLD FORMAT:
          READ(ITB,REC=1)NTEMP,(TBAND(ITEMP),ITEMP=1,NTEMP),            &
     &      IBNDWD,MNFREQ,MXFREQ,LSTREC,EDGENR
          NTLSUB=1
      ENDIF
      IF(EDGENR.LT.EDGEMN)THEN
          EDGENR=EDGEMN
      ELSEIF(EDGENR.GT..5*IBNDWD)THEN
          EDGENR=.5*IBNDWD
      ENDIF
      EDGEFR=IBNDWD-EDGENR

!     OPEN THE FORMATTED BAND MODEL FILE, UNIT ITBX.
      INQUIRE(ITBX,OPENED=LOPEN)
      IF(LOPEN)CLOSE(ITBX,STATUS='KEEP')
      IF(IBNDWD.EQ.15)THEN

!         CHECK THAT FILE EXISTS:
          FULNAM=DATDIR//'CFC99_15.ASC'
          INQUIRE(FILE=FULNAM,EXIST=LEXIST)
          IF(.NOT.LEXIST)THEN
              WRITE(IPR,'(/2A,/18X,3A)')' ERROR in OPENBM: ',           &
     &          ' The CFC molecular band model data file',              &
     &          ' "',DATDIR,'CFC99_15.ASC" does not exist.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP                                                      &
     &          'Error:  CFC molecular band model data file not found.'
          ENDIF
          IF(.NOT.LJMASS)WRITE(IPR,'(/A,15X,A)')                        &
     &      ' CFC BAND MODEL DATA FILE:',FULNAM
          OPEN(ITBX,FILE=FULNAM,STATUS='OLD',FORM='FORMATTED')
      ELSEIF(IBNDWD.EQ.5)THEN

!         CHECK THAT FILE EXISTS:
          FULNAM=DATDIR//'CFC99_05.ASC'
          INQUIRE(FILE=FULNAM,EXIST=LEXIST)
          IF(.NOT.LEXIST)THEN
              WRITE(IPR,'(/2A,/18X,3A)')' ERROR in OPENBM: ',           &
     &          ' The CFC molecular band model data file',              &
     &          ' "',DATDIR,'CFC99_05.ASC" does not exist.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP                                                      &
     &          'Error:  CFC molecular band model data file not found.'
          ENDIF
          IF(.NOT.LJMASS)WRITE(IPR,'(/A,15X,A)')                        &
     &      ' CFC BAND MODEL DATA FILE:',FULNAM
          OPEN(ITBX,FILE=FULNAM,STATUS='OLD',FORM='FORMATTED')
      ELSE

!         CHECK THAT FILE EXISTS:
          FULNAM=DATDIR//'CFC99_01.ASC'
          INQUIRE(FILE=FULNAM,EXIST=LEXIST)
          IF(.NOT.LEXIST)THEN
              WRITE(IPR,'(/2A,/18X,3A)')' ERROR in OPENBM: ',           &
     &          ' The CFC molecular band model data file',              &
     &          ' "',DATDIR,'CFC99_01.ASC" does not exist.'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP                                                      &
     &          'Error:  CFC molecular band model data file not found.'
          ENDIF
          IF(.NOT.LJMASS)WRITE(IPR,'(/A,15X,A)')                        &
     &      ' CFC BAND MODEL DATA FILE:',FULNAM
          OPEN(ITBX,FILE=FULNAM,STATUS='OLD',FORM='FORMATTED')
      ENDIF

!     WRITE WATER CONTINUUM INFORMATION:
      IF(.NOT.LJMASS)WRITE(IPR,'(/A)')' Version 2.4 of the'//           &
     &  ' Clough-Kneizys Water Continuum Data from LBLRTM (24mar2000).'

!     RETURN TO DRIVER:
      RETURN
      END
