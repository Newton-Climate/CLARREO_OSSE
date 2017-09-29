      LOGICAL FUNCTION  WRTDIM(DIMNAM, MINVAL)

!                     Write name of too-small symbolic dimension and the
!                        value it should be increased to;  return 'TRUE'

!      INPUT :  DIMNAM=name of symbolic dimension which is too small
!                        (character of any length)
!               MINVAL=value to which that dimension should be
!                        increased (at least)

      CHARACTER*(*)  DIMNAM
      INTEGER        MINVAL

      WRITE (*, '(3A,I7)')  ' ****  SYMBOLIC DIMENSION  ', DIMNAM,      &
     &                     '  SHOULD BE INCREASED TO AT LEAST ', MINVAL
      WRTDIM=.TRUE.

      RETURN
      END
