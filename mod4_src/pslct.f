      INTEGER FUNCTION PSLCT(ANGLE,RANGE,BETA)

!     THIS ROUTINE RETURNS THE AN INTEGER VALUE INDICATING THE TYPE
!     OF SLANT (ITYPE=2) PATH.  THE FOLLOWING VALUES ARE RETURNED:
!          PSLCT = 21 FOR CASE 2A (H1,H2,ANGLE)
!          PSLCT = 22 FOR CASE 2B (H1,ANGLE,RANGE)
!          PSLCT = 23 FOR CASE 2C (H1,H2,RANGE)
!          PSLCT = 24 FOR CASE 2D (H1,H2,BETA)

!     H1     H2     ANGLE     RANGE      BETA        CASE  PSLCT
!--------------------------------------------     ----------------
!     X      X      X                                 2A     21

!     X             X         X                       2B     22

!     X      X                X                       2C     23

!     X      X                           X            2D     24

!     DECLARE INPUTS
      REAL ANGLE,RANGE,BETA
      IF(BETA.GT.0.)THEN

!         BETA > 0 IMPLIES CASE 2D (H1,H2,BETA)
          PSLCT=24
      ELSEIF(RANGE.GT.0.)THEN

!         RANGE > 0 AND BETA = 0 IMPLIES CASE 2B OR 2C.
          IF(ANGLE.EQ.0.)THEN

!             RANGE > 0, BETA = 0 AND ANGLE = 0 IMPLIES CASE 2C
              PSLCT=23
          ELSE

!             RANGE > 0, BETA = 0 AND ANGLE NON-ZERO IMPLIES CASE 2B
              PSLCT=22
          ENDIF
      ELSE

!         RANGE =0 AND BETA = 0 IMPLIES CASE 2A
          PSLCT=21
      ENDIF
      RETURN
      END
