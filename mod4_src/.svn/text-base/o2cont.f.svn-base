      SUBROUTINE O2CONT(V,SIGMA,ALPHA,BETA)

!     THIS ROUTINE IS DRIVEN BY FREQUENCY, RETURNING ONLY THE
!     O2 COEFFICIENTS, INDEPENDENT OF TEMPERATURE.

!  *******************************************************************
!  *  THESE COMMENTS APPLY TO THE COLUME ARRAYS FOR:                 *
!  *       PBAR*UBAR(O2)                                             *
!  *       PBAR*UBAR(O2)*DT                                          *
!  *   AND PBAR*UBAR(O2)*DT*DT    WHERE:  DT=TBAR-220.               *
!  *  THAT HAVE BEEN COMPILED IN OTHER PARTS OF THE LOWTRAN CODE     *
!  *                                                                 *
!  *  LOWTRAN7 COMPATIBLE:                                           *
!  *  O2 CONTINUUM SUBROUTINE FOR 1395-1760CM-1                      *
!  *  MODIFIED BY G.P. ANDERSON, APRIL '88                           *
!  *                                                                 *
!  *  THE EXPONENTIAL TEMPERATURE EMPLOYED IN THE FASCOD2 ALGORITHM  *
!  *  (SEE BELOW) IS NOT READILY SUITABLE FOR LOWTRAN.  THEREFORE    *
!  *  THE EXPONENTIALS HAVE BEEN LINEARLY EXPANDED, KEEPING ONLY THE *
!  *  LINEAR AND QUADRATIC TERMS:                                    *
!  *                                                                 *
!  *  EXP(A*DT)=1.+ A*DT + (A*DT)**2/2. + ....                       *
!  *                                                                 *
!  *     EXP(B*DT*DT)=1.+ B*DT*DT + (B*DT*DT)**2/2. + ....           *
!  *                                                                 *
!  *  THE PRODUCT OF THE TWO TERMS IS:                               *
!  *                                                                 *
!  *     (1. + A*DT + (A*A/2. + B)*DT*DT )                           *
!  *                                                                 *
!  *  THIS EXPANSION ONLY WORKS WELL FOR SMALL VALUES OF X IN EXP(X) *
!  *                                                                 *
!  *  SINCE DT = T-220., THE APPROXIMATION IS VERY GOOD UNTIL        *
!  *  T.GT.260. OR DT.GT.40.   AT T=280, THE MAXIMUM ERRORS ARE STILL*
!  *  LESS THAN 10% BUT AT T=300, THOSE ERRORS ARE AS LARGE AS 20%   *
!  *******************************************************************

!     THE FOLLOWING COMMENTS ARE EXCERPTED DIRECTLY FROM FASCOD2

!      THIS SUBROUTINE CONTAINS THE ROGERS AND WALSHAW
!      EQUIVALENT COEFFICIENTS DERIVED FROM THE THEORETICAL
!      VALUES SUPPLIED BY ROLAND DRAYSON. THESE VALUES USE
!      THE SAME DATA AS TIMOFEYEV AND AGREE WITH TIMOFEYEV'S RESULTS.
!      THE DATA ARE IN THE FORM OF STRENGTHS(O2SO) AND TWO
!      COEFFICIENTS (O2A & O2B),  WHICH ARE USED TO CORRECT FOR
!      TEMPERATURE. THE DEPENDENCY ON PRESSURE SQUARED
!      IS CONTAINED IN THE P*WO2 PART OF THE CONSTANT.
!      NOTE THAT SINCE THE COEFFICIENTS ARE FOR AIR, THE
!      THE STRENGTHS ARE DIVIDED BY THE O2 MIXING RATIO FOR
!      DRY AIR OF 0.20946 (THIS IS ASSUMED CONSTANT).
!      ORIGINAL FORMULATION OF THE COEFFICIENTS WAS BY LARRY GORDLEY.
!      THIS VERSION WRITTEN BY EARL THOMPSON, JULY 1984.

      INTEGER NPTO2
      REAL O2S0,O2A,O2B,V1O2,DVO2
      COMMON/O2C/O2S0(74),O2A(74),O2B(74),V1O2,DVO2,NPTO2

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /O2C/
      EXTERNAL BO2C
      SIGMA =0.
      ALPHA =0.
      BETA  =0.
      IF(V .LT. 1395) GO TO 30
      IF(V .GT. 1760) GO TO 30

      CALL O2INT(V,V1O2,DVO2,NPTO2,C,O2S0,A,O2A,B,O2B)

!     OLD 'FASCOD2' TEMPERATURE DEPENDENCE USING BLOCK DATA ARRAYS

!     C(J)=O2S0(I)* EXP(O2A(I)*TD+O2B(I)*TD*TD) /(0.20946*VJ)

!     NEW COEFFICIENT DEFINITIONS FOR LOWTRAN FORMULATION

      ALPHA= A
      BETA=A**2/2.+B
      SIGMA=C/0.20946

!     NEW 'LOWTRAN7' TEMPERATURE DEPENDENCE

!     THIS WOULD BE THE CODING FOR THE LOWTRAN7 FORMULATION, BUT
!       BECAUSE THE T-DEPENDENCE IS INCLUDED IN THE AMOUNTS, ONLY
!       THE COEFFICIENTS (SIGMA, ALPHA & BETA) ARE BEING RETURNED

!     C(J)=SIGMA*(1.+ALPHA*TD+BETA*TD*TD)

!     THE COEFFICIENTS FOR O2 HAVE BEEN MULTIPLIED BY A FACTOR
!     OF 0.78  [RINSLAND ET AL, 1989: JGR 94; 16,303 - 16,322.].
      O2FAC = 0.78
      SIGMA=O2FAC*SIGMA
30    RETURN
      END
