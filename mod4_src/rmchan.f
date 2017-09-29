      SUBROUTINE RMCHAN(IFILTR)

!     THIS ROUTINE REMOVES THE LAST CHANNEL ENCOUNTERED
!     DUE TO AN ERROR IN ITS PROCESSING.

!     ARGUMENTS:
!       IFILTR   UNIT NUMBER OF FILTER FUNCTION FILE.
      INTEGER IFILTR

!     PARAMETERS:
      INCLUDE 'CHANLS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
      INTEGER IBIN

!     APPEND WARNING:
      WRITE(IPR,'(10X,A,I4,A)')' Channel',NCHAN,                        &
     &  ' is being removed and no additional channels will be used.'

!     REMOVE NCHAN FROM COMMON/CHANEL/
      DO 10 IBIN=MNBIN,MXBIN
          IF(NUMCHN(IBIN).GT.0)THEN
              IF(LSTCHN(NUMCHN(IBIN),IBIN).EQ.NCHAN)                    &
     &          NUMCHN(IBIN)=NUMCHN(IBIN)-1
          ENDIF
   10 CONTINUE
      NCHAN=NCHAN-1
      CLOSE(IFILTR)
      RETURN
      END
