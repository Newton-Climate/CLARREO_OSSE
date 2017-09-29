      SUBROUTINE  PRTINT(UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI,          &
     &                    MAXULV, MAXUMU)

!         PRINTS THE INTENSITY AT USER POLAR AND AZIMUTHAL ANGLES

!     ALL ARGUMENTS ARE DISORT INPUT OR OUTPUT VARIABLES

!+---------------------------------------------------------------------+

!     PARAMETERS:
      INCLUDE 'ERROR.h'

      REAL   PHI(*), UMU(*), UTAU(*), UU(MAXUMU, MAXULV, *)

!     NO NEED TO PRINT IF JMASS
      IF(LJMASS) RETURN

      WRITE (*, '(//,A)')                                               &
     &         ' *********  I N T E N S I T I E S  *********'
      LENFMT=10
      NPASS=1 + NPHI / LENFMT
      IF (MOD(NPHI,LENFMT) .EQ. 0)  NPASS=NPASS - 1
      DO 10  LU=1, NTAU
         DO 10  NP=1, NPASS
            JMIN=1 + LENFMT * (NP-1)
            JMAX=MIN0(LENFMT*NP, NPHI)
            WRITE(*,101)  (PHI(J), J=JMIN, JMAX)
            DO 10  IU=1, NUMU
               IF(IU.EQ.1)  WRITE(*,102)  UTAU(LU), UMU(IU),            &
     &           (UU(IU,LU,J), J=JMIN, JMAX)
               IF(IU.GT.1)  WRITE(*,103)  UMU(IU),                      &
     &           (UU(IU,LU,J), J=JMIN, JMAX)
10    CONTINUE

      RETURN

101   FORMAT(/, 3X,'          POLAR   AZIMUTH ANGLES (DEGREES)',        &
     &        /, 3X,'OPTICAL   ANGLE',                                  &
     &        /, 3X,' DEPTH   COSINE', 10F11.2)
102   FORMAT(F10.4, F8.4, 1P,10E11.3)
103   FORMAT(10X,   F8.4, 1P,10E11.3)
      END
