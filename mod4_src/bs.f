      SUBROUTINE BS(I,A,B,N,S)
!**********************************************************************
      DIMENSION B(9)

!             THIS SUBROUTINE DOES THE BINARY SEARCH FOR THE INDEX I
!             SUCH THAT A IS IN BETWEEN B(I) AND B(I+1)
!             AND CALCULATES THE INTERPOLATION PARAMETER S
!             SUCH THAT A=S*B(I+1)+(1.-S)*B(I)

      I=1
      J=N
   10 M=(I+J)/2
      IF(A.LE.B(M)) THEN
      J=M
      ELSE
      I=M
      END IF
      IF(J.GT.I+1) GO TO 10
      S=(A-B(I))/(B(I+1)-B(I))
      RETURN
      END
