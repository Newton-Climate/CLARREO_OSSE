      DOUBLE PRECISION FUNCTION FACFNC(B,Z)

!     THIS ROUTINE EVALUATES THE FUNCTION

!         SQRT(PI)        2   2
!         --------  EXP( Z / B )  ERFC( Z / B )
!           2 B

!     THE VARIABLES "B" AND "Z" ARE ASSUMED TO BE GREATER THAN ZERO.
      IMPLICIT NONE
      DOUBLE PRECISION B,T,Z,DENOM,CUTOFF,P,A1,A2,A3,A4,A5
      DATA P,A1,A2,A3,A4,A5/.3275911D0,.225836846D0,-.252128668D0,      &
     &  1.259695129D0,-1.287822453D0,.940646070D0/,CUTOFF/1.883D0/

!     THE PARAMETER "CUTOFF" WAS SET TO 1.883 BECAUSE THE SELECTED
!     PADE-APPROXIMATE AND THE RATIONAL APPROXIMATION [ABRAMOWITZ AND
!     STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, EQ 7.1.26] CROSS
!     AT THAT VALUE.  ABRAMOWITZ AND STEGUN INDICATE THAT THE MAXIMUM
!     ABSOLUTE ERROR IN THE RATIONAL APPROXIMATION FOR THE ERROR
!     FUNCTION IS 1.5E-7.  ASSUMING B > .01, THE MAXIMUM POSSIBLE
!     ABSOLUTE ERROR OF THE RATIONAL APPROXIMATION IN THIS ROUTINE
!     IS 4.6E-4.  AT THIS MAXIMUM (B=.01, Z=.01883), THE VALUE OF THE
!     FUNCTION IS 23.8 CORRESPONDING TO A RELATIVE ERROR OF 1.9E-5.

      IF(Z.GT.CUTOFF*B)THEN
          T=(B/Z)**2
          FACFNC=(8+T*(140+T*(690+T*(975+T*192))))/                     &
     &      (Z*(16+T*(288+T*(1512+T*(2520+945*T)))))
      ELSE
          DENOM=B+P*Z
          T=B/DENOM
          FACFNC=(A1+T*(A2+T*(A3+T*(A4+T*A5))))/DENOM
      ENDIF
      RETURN
      END
