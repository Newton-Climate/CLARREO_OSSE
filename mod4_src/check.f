      SUBROUTINE CHECKP(P,JUNITP)

!     PARAMETERS:
      INCLUDE 'ERROR.h'

!     INPUT ARGUMENTS:
!       P        PRESSURE.
!       JUNITP   UNIT FLAG FOR P.

!     OUTPUT ARGUMENTS:
!       P        PRESSURE [MBARS].
      REAL P
      INTEGER JUNITP
      IF(JUNITP.EQ.11)THEN

!         CONVERT FROM ATM TO MBARS:
          P=1013.25*P
      ELSEIF(JUNITP.EQ.12)THEN

!         CONVERT FROM TORR TO MBARS
          P=(1013.25/760)*P
      ELSEIF(JUNITP.GT.12)THEN
           IF(LJMASS)CALL WRTBUF(FATAL)
           STOP 'CHECKP:  Cannot interpret unit specifier for pressure.'
      ENDIF
      RETURN
      END
      SUBROUTINE CHECKT(T,JUNITT)

!     PARAMETERS:
      INCLUDE 'ERROR.h'

!     INPUT ARGUMENTS:
!       T        TEMPERATURE.
!       JUNITT   UNIT FLAG FOR TEMPERATURE.

!     OUTPUT ARGUMENTS:
!       T        TEMPERATURE [K].
      REAL T
      INTEGER JUNITT
      IF(JUNITT.EQ.11)THEN

!         CONVERT FROM C TO T:
          T=T+273.15
      ELSEIF(JUNITT.GT.11)THEN
           IF(LJMASS)CALL WRTBUF(FATAL)
           STOP                                                         &
     &       'CHECKT:  Cannot interpret unit specifier for temperature.'
      ENDIF
      RETURN
      END
