      SUBROUTINE  UPBEAM(ARRAY, CC, CMU, DELM0, FBEAM, GL, IPVT,MAZIM,  &
     &                    MXCMU, NN, NSTR, PI, UMU0, WK, YLM0, YLMC, ZJ,&
     &                    ZZ)

!         FINDS THE INCIDENT-BEAM PARTICULAR SOLUTION  OF SS(18)

!     ROUTINES CALLED:  SGECO, SGESL

!   I N P U T    V A R I A B L E S:

!       CC     :  CAPITAL-C-SUB-IJ IN EQ. SS(5)
!       CMU    :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       DELM0  :  KRONECKER DELTA, DELTA-SUB-M0
!       GL     :  DELTA-M SCALED LEGENDRE COEFFICIENTS OF PHASE FUNCTION
!                    (INCLUDING FACTORS 2L+1 AND SINGLE-SCATTER ALBEDO)
!       MAZIM  :  ORDER OF AZIMUTHAL COMPONENT
!       YLM0   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
!                 AT THE BEAM ANGLE
!       YLMC   :  NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
!                 AT THE QUADRATURE ANGLES
!       (REMAINDER ARE 'DISORT' INPUT VARIABLES)

!   O U T P U T    V A R I A B L E S:

!       ZJ     :  RIGHT-HAND SIDE VECTOR CAPITAL-X-SUB-ZERO IN SS(19);
!                 ALSO THE SOLUTION VECTOR CAPITAL-Z-SUB-ZERO
!                 AFTER SOLVING THAT SYSTEM
!       ZZ     :  PERMANENT STORAGE FOR -ZJ-, BUT RE-ORDERED

!   I N T E R N A L    V A R I A B L E S:

!       ARRAY  :  COEFFICIENT MATRIX IN LEFT-HAND SIDE OF EQ. SS(19)
!       IPVT   :  INTEGER VECTOR OF PIVOT INDICES REQUIRED BY *LINPACK*
!       WK     :  SCRATCH ARRAY REQUIRED BY *LINPACK*
!+---------------------------------------------------------------------+

      INTEGER  IPVT(*)
      REAL     ARRAY(MXCMU,*), CC(MXCMU,*), CMU(*), GL(0:*),            &
     &         WK(*), YLM0(0:*), YLMC(0:MXCMU,*), ZJ(*), ZZ(*)
      INTEGER DUMMY1,DUMMY2 !DRF for totalview

      DO 40  IQ=1, NSTR

         DO 10  JQ=1, NSTR
            ARRAY(IQ,JQ)=- CC(IQ,JQ)
10       CONTINUE
         ARRAY(IQ,IQ)=1. + CMU(IQ) / UMU0 + ARRAY(IQ,IQ)

         SUM=0.
         DO 20  K=MAZIM, NSTR-1
            SUM=SUM + GL(K) * YLMC(K,IQ) * YLM0(K)
20       CONTINUE
         ZJ(IQ)=(2. - DELM0) * FBEAM * SUM / (4.0*PI)
40    CONTINUE
!                  ** FIND L-U (LOWER/UPPER TRIANGULAR) DECOMPOSITION
!                  ** OF -ARRAY- AND SEE IF IT IS NEARLY SINGULAR
!                  ** (NOTE:  -ARRAY- IS DESTROYED)
      RCOND=0.0
      CALL  SGECO(ARRAY, MXCMU, NSTR, IPVT, RCOND, WK)
      IF (1.0+RCOND .EQ. 1.0) THEN
                               CALL  ERRMSG                             &
     &   ('UPBEAM--SGECO SAYS MATRIX NEAR SINGULAR',.FALSE.)
         DUMMY1=0
         DUMMY2=0
         DUMMY1 = DUMMY1+DUMMY2
      ENDIF
!                ** SOLVE LINEAR SYSTEM WITH COEFF MATRIX -ARRAY-
!                ** (ASSUMED ALREADY L-U DECOMPOSED) AND R.H. SIDE(S)
!                ** -ZJ-;  RETURN SOLUTION(S) IN -ZJ-
      JOB=0
      CALL  SGESL(ARRAY, MXCMU, NSTR, IPVT, ZJ, JOB)

      DO 50  IQ=1, NN
         ZZ(IQ+NN)  =ZJ(IQ)
         ZZ(NN+1-IQ)=ZJ(IQ+NN)
50    CONTINUE

      RETURN
      END
