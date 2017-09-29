      INTEGER FUNCTION  ISAMAX( N, SX, INCX )

!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR OF INTEREST
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

! --OUTPUT-- ISAMAX   FIRST I, I = 1 TO N, TO MAXIMIZE
!                         ABS(SX(1+(I-1)*INCX))

      REAL SX(*), SMAX, XMAG

      IF( N.LE.0 ) THEN
         ISAMAX = 0
      ELSE IF( N.EQ.1 ) THEN
         ISAMAX = 1
      ELSE
         SMAX = 0.0
         II = 1
         DO 20  I = 1, 1+(N-1)*INCX, INCX
            XMAG = ABS(SX(I))
            IF( SMAX.LT.XMAG ) THEN
               SMAX = XMAG
               ISAMAX = II
            ENDIF
            II = II + 1
   20    CONTINUE
      ENDIF

      RETURN
      END
