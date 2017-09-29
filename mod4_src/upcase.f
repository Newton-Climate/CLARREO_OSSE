      CHARACTER*(*) FUNCTION UPCASE(STRING)
!     Created by:  Dr. William M. Cornette
!     Created on:  Tue Jul 28 14:49:15 1992
!     Revised on:  Mon Aug  2 11:06:27 1993

!**** This FUNCTION converts STRING from lower case to upper case.

! STRING - CHARACTER*(*) Variable - Input string

!**** INTRINSIC and EXTERNAL Declarations

      INTRINSIC        LEN,INDEX

!**** Argument Declaration

      CHARACTER*(*)    STRING

!**** Local Variable Declarations

      INTEGER          I,LOC
      CHARACTER*26     UPPER,LOWER

      DATA             UPPER /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA             LOWER /'abcdefghijklmnopqrstuvwxyz'/

      DO 101 I = 1,LEN(STRING)
         UPCASE(I:I)=STRING(I:I)

         LOC=INDEX(LOWER,STRING(I:I))
         IF (LOC.GT.0) UPCASE(I:I)=UPPER(LOC:LOC)
 101  CONTINUE

      RETURN
      END
