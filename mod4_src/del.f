      FUNCTION   DEL(PSIO,DELO,BETA,IARBO)

!     FUNCTION DEL RETURNS THE VALUE OF THE SUN'S ZENITH ANGLE
!     AT ANY POINT ALONG THE OPTICAL PATH BASED UPON STRAIGHT
!     LINE GEOMETRY (NO REFRACTION). THIS ANGLE IS USED TO SPECIFY
!     THE SCATTERING POINT TO SUN PATHS. THE BENDING DUE TO REFRACTION
!     ALONG THIS PATH IS DETERMINED BY THE GEO ROUTINES. IF THE BENDING
!     IS GREATER THAN ONE DEGREE THE ZENITH ANGLE IS CORRECTED ACCORDING
!     AND THE PATH CALCULATION IS REPEATED.

!     COMMONS:
!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT
      IF(IARBO.EQ.0) GO TO 10
!     SPECIAL CASES IF PSIO IS ARBITRARY
      IF(IARBO.EQ.1) DEL=DELO
      IF(IARBO.EQ.2) DEL=BETA
      IF(IARBO.EQ.3) DEL=0.0
      RETURN
10    CONTINUE
      PSIOR=PSIO/DEG
      DELOR=DELO/DEG
      BETAR=BETA/DEG
!     GENERAL CASE
      X=COS(DELOR)*COS(BETAR)+SIN(DELOR)*SIN(BETAR)*COS(PSIOR)
      DEL=DEG*ACOS(X)
      RETURN
      END