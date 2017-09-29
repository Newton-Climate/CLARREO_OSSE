      INTEGER FUNCTION NUNIT()

!     PARAMETERS:
      INCLUDE 'ERROR.h'

!     THIS INTEGER FUNCTION DETERMINES AN UNUSED
!     FILE UNIT NUMBER BETWEEN 10 AND 99.

!     COMMONS:
      INCLUDE 'IFIL.h'
      EXTERNAL DEVCBD

!     LOCAL VARIABLES
      LOGICAL LOPEN

!     FIND UNUSED UNIT NUMBER AND RETURN.
      DO 10 NUNIT=99,10,-1
          INQUIRE(UNIT=NUNIT,OPENED=LOPEN)
          IF(.NOT.LOPEN)RETURN
   10 CONTINUE

!     NO FILE UNIT NUMBER AVAILABLE
      WRITE(IPR ,'(/A)')                                                &
     &  ' Error in routine NUNIT:  No file unit numbers available?'
      IF(LJMASS)CALL WRTBUF(FATAL)
      STOP ' Error in routine NUNIT:  No file unit numbers available?'
      END
