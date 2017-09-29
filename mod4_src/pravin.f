      SUBROUTINE  PRAVIN(UMU, NUMU, MAXUMU, UTAU, NTAU, U0U)
!        PRINT AZIMUTHALLY AVERAGED INTENSITIES AT USER ANGLES

!     PARAMETERS:
      INCLUDE 'ERROR.h'

      REAL     UMU(*), UTAU(*), U0U(MAXUMU,*)

!     NO NEED TO PRINT IF JMASS
      IF(LJMASS) RETURN

      WRITE(*, '(//,A)')                                                &
     &         ' *********  AZIMUTHALLY AVERAGED INTENSITIES '          &
     &       // '(USER POLAR ANGLES)  *********'
      LENFMT=8
      NPASS=1 + NUMU / LENFMT
      IF (MOD(NUMU,LENFMT) .EQ. 0)  NPASS=NPASS - 1
      DO 10  NP=1, NPASS
         IUMIN=1 + LENFMT * (NP-1)
         IUMAX=MIN0(LENFMT*NP, NUMU)
         WRITE(*,101)  (UMU(IU), IU=IUMIN, IUMAX)
         DO 10  LU=1, NTAU
            WRITE(*,102) UTAU(LU), (U0U(IU,LU), IU=IUMIN,IUMAX)
 10   CONTINUE

      RETURN

101   FORMAT(/, 3X,'OPTICAL   POLAR ANGLE COSINES',                     &
     &        /, 3X,'  DEPTH', 8F14.5)
102   FORMAT(0P,F10.4, 1P,8E14.4)
      END
