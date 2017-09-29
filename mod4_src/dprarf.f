      DOUBLE PRECISION FUNCTION DPRARF(H,SH,GAMMA)

!     DOUBLE PRECISION (DP) VERSION OF THE ROUTINE PREVIOUSLY CALLED
!     RADREF.  COMPUTES THE RADIUS OF CURVATURE OF THE REFRACTED RAY FOR
!     HORIZONTAL PATH:  DPRARF = DPANDX/ D(DPANDX)/D(RADIUS)
!     DPANDX = DOUBLE PRECISION (DP) VERSION OF ANDEX
      DOUBLE PRECISION H,SH,GAMMA,HSH

!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT
      IF(SH.EQ.0.)GOTO20
      HSH=H/SH
      IF(HSH .GT. BIGEXP) GO TO 20
      DPRARF=SH*(DBLE(1.)+EXP(HSH )/GAMMA)
      RETURN
   20 DPRARF=DBLE(BIGNUM)
      RETURN
      END