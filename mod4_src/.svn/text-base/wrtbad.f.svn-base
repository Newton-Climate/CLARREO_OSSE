      LOGICAL FUNCTION  WRTBAD(VARNAM)

!                   Write names of erroneous variables and return 'TRUE'

!      INPUT :   VARNAM=name of erroneous variable to be written
!                         (character of any length)

      CHARACTER*(*)  VARNAM
      INTEGER        MAXMSG, NUMMSG
      SAVE           NUMMSG, MAXMSG
      DATA           NUMMSG / 0 /,  MAXMSG / 50 /

      WRTBAD=.TRUE.
      NUMMSG=NUMMSG + 1
      WRITE (*, '(3A)')  ' ****  INPUT VARIABLE  ', VARNAM,             &
     &                     '  IN ERROR  ****'
      IF (NUMMSG.EQ.MAXMSG)                                             &
     &   CALL  ERRMSG ('TOO MANY INPUT ERRORS.  ABORTING...$', .TRUE.)

      RETURN
      END
