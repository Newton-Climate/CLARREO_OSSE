      SUBROUTINE DPEXNT(X,X1,X2,A)

!     DOUBLE PRECISION VERSION OF THE ROUTINE EXPINT,
!     EXPONENTIAL INTERPOLATION
      DOUBLE PRECISION X,X1,X2,A
      IF(X1.LT.0.)X1=DBLE(0.)
      IF(X2.LT.0.)X2=DBLE(0.)
      IF(X1.EQ.0.0 .OR. X2.EQ.0.0)  GO TO 100
      X = X1*(X2/X1)**A
      RETURN
  100 X = X1+(X2-X1)*A
      RETURN
      END
