      INTEGER FUNCTION LENSTR(STRING)

!     TRIMS LEADING BLANKS FROM STRING AND RETURNS FINAL NON-BLANK
!     CHARACTER LOCATION.  FOR BLANK STRINGS, ZERO IS RETURNED.

!     DECLARE ARGUMENTS:
      CHARACTER STRING*(*)

!     LOCAL VARIABLES
      INTEGER IEND,IBEG

!     DETERMINE TRAILING NON-BLANK CHARACTER:
      DO 10 IEND=LEN(STRING),1,-1
          IF(STRING(IEND:IEND).NE.' ')GOTO20
   10 CONTINUE

!     BLANK STRING:
      LENSTR=0
      RETURN

!     CHECK IF ALREADY LEFT-JUSTIFIED:
   20 CONTINUE
      IF(STRING(1:1).NE.' ')THEN
          LENSTR=IEND
          RETURN
      ENDIF

!     DETERMINE LEADING NON-BLANK CHARACTER:
      DO 30 IBEG=2,IEND-1
          IF(STRING(IBEG:IBEG).NE.' ')THEN

!             LEFT-JUSTIFY STRING, SET LENGTH AND RETURN
              STRING=STRING(IBEG:IEND)
              LENSTR=IEND-IBEG+1
              RETURN
          ENDIF
   30 CONTINUE

!     SINGLE CHARACTER STRING:
      STRING=STRING(IEND:IEND)
      LENSTR=1
      RETURN
      END