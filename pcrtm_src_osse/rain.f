      SUBROUTINE RAIN(VBAND5,WPATH,EXT,ABT,SCT,ASYMR)
      IMPLICIT NONE

!     INPUT ARGUMENTS:
      REAL VBAND5,WPATH(0:*)

!     OUTPUT ARGUMENTS:
      REAL EXT,ABT,SCT,ASYMR

!     COMMONS:

!     /RCNSTN/
!       PI       THE CONSTANT PI.
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

!     /CARD2/
!       IHAZE    BOUNDARY LAYER AEROSOL MODEL NUMBER.
!       ISEASN   SEASON NUMBER (1=SPRING-SUMMER, 2=FALL-WINTER).
!       IVULCN   VOLCANIC AEROSOL MODEL NUMBER.
!       ICSTL    COASTAL AIRMASS MODEL NUMBER.
!       ICLD     CLOUD MODEL NUMBER.
!       IVSA     VERTICAL STRUCTURE ALGORITHM (0=OFF, 1=ON).
!       VIS      SURFACE VISIBILITY (GROUND METEOROLOGICAL RANGE) [KM].
!       WSS      CURRENT WIND SPEED (M/S).
!       WHH      24-HOUR WIND SPEED (M/S).
!       RAINRT   RAIN RATE (MM/HR)
!       LSAP     LOGICAL FLAG FOR SPECTRAL AEROSOL PROFILES INPUT.
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA
      REAL VIS,WSS,WHH,RAINRT
      LOGICAL LSAP
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 &
     &  VIS,WSS,WHH,RAINRT,LSAP

!     LOCAL VARIABLES:
      INTEGER PHASE
      REAL TMPRN,VTEMP,RFD,RAINAV,CSSA

!     FUNCTIONS:
      REAL TNRAIN
      TMPRN=WPATH(61)/WPATH(3)
      VTEMP=1.438786*VBAND5

      IF(VTEMP.GE.TMPRN*BIGEXP)THEN
          RFD=VBAND5
      ELSE
          RFD=EXP(-VTEMP/TMPRN)
          RFD=VBAND5*(1-RFD)/(1+RFD)
      ENDIF
      RAINAV=(WPATH(3)/WPATH(62))**(1/.63)
      IF(VBAND5.LT.250.)THEN
          IF(ICLD.LE.10)THEN
              PHASE=1
          ELSE
              PHASE=2
          ENDIF

!         CALL SCATTERING ROUTINE TO OBTAIN ASYMMETRY
!         FACTOR AND RATIO OF ABSORPTION TO EXTINCTION DUE
!         TO RAIN WITHIN RANGE OF 19 TO 231 GHZ.
!         EXTRAPOLATE ABOVE AND BELOW THAT FREQ RANGE
          CALL RNSCAT(VBAND5,RAINAV,TMPRN,PHASE,1,CSSA,ASYMR)
      ELSE
          CSSA=.5
          ASYMR=.85
      ENDIF
      EXT=TNRAIN(RAINAV,VBAND5,TMPRN,RFD)*WPATH(62)
      ABT=EXT*CSSA
      SCT=EXT*(1-CSSA)
      RETURN
      END