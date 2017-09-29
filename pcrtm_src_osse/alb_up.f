      REAL FUNCTION ALB_UP(GU,KK,LAYTOP,LL,NLYR,                        &
     &  NN,NSTR,IU,TAUCPR,MU,TAUTOP)

!     COMPUTE INTENSITY COMPONENTS AT USER OUTPUT ANGLES FOR AZIMUTHAL
!     EXPANSION TERMS IN EQ. SD(2).  ALL THE EXPONENTIAL FACTORS
!     (EXP1, EXPN,... ETC.) COME FROM THE SUBSTITUTION OF CONSTANTS OF
!     INTEGRATION IN EQ. SC(12) INTO EQS. S1(8-9).  THEY ALL HAVE
!     NEGATIVE ARGUMENTS SO THERE SHOULD NEVER BE OVERFLOW PROBLEMS.
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       GU       EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES
!                (I.E., G IN EQ. SC(1)).
!       KK       EIGENVALUES OF COEFFICIENT MATRIX IN EQ. SS(7).
!       LAYTOP   LAYER INDEX FOR TOP OF SUB-MEDIUM.
!       LL       CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED BY
!                SOLVING SCALED VERSION OF EQ. SC(5); EXPONENTIAL TERM
!                OF EQ. SC(12) NOT INCLUDED.
!       NLYR     NUMBER OF LAYERS IN MEDIUM.
!       NN       ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2).
!       NSTR     NUMBER OF DISCRETE ORDINATE STREAMS.
!       IU       INDEX FOR USER ANGLE COSINE.
!       TAUCPR   CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED).
!       MU       USER ANGLE COSINE.
!       TAUTOP   OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
!                COORDINATES;  EQUAL TO -UTAU- IF NO DELTA-M.
      INTEGER LAYTOP,NN,NLYR,NSTR,IU
      REAL GU(0:MXUMU,MXCMU,*),KK(MXCMU,*),LL(MXCMU,*),TAUCPR(0:*),MU,  &
     &  TAUTOP

!     INTERNAL ARGUMENTS/VARIABLES:
!       LC       LOOP INDEX FOR COMPUTATIONAL LAYERS.
!       LCM1     LC MINUS ONE.
!       DTAU     NADIR OPTICAL DEPTH OF CURRENT COMPUTATIONAL LAYER.
!       DTAU2    NEGATIVE OF NADIR OPTICAL DEPTH FROM USER TO TOP OF A
!                LOWER ALTITUDE LAYER.
!       WK       SCRATCH SCALAR FOR SAVING 'EXP' EVALUATIONS.
      INTEGER LC,LCM1,IQ,JQ
      REAL DTAU,EXP1,EXP2,DTAU2,WK,DTAU1,KK_IQ,MU_KK

!     FOR DOWNWARD INTENSITY, INTEGRATE FROM TOP TO 'LAYTOP-1' IN EQ.
!     S1(8); FOR UPWARD, INTEGRATE FROM BOTTOM TO 'LAYTOP+1' IN S1(9):
      ALB_UP=0.
      LC=NLYR
      DTAU2=TAUTOP-TAUCPR(LC)
      EXP2=EXP(DTAU2/MU)
      DO LCM1=NLYR-1,LAYTOP,-1
          DTAU=TAUCPR(LC)-TAUCPR(LCM1)
          EXP1=EXP2
          DTAU2=TAUTOP-TAUCPR(LCM1)
          EXP2=EXP(DTAU2/MU)
          JQ=NSTR
          DO IQ=1,NN

!             KK_IQ NEGATIVE:
              KK_IQ=KK(IQ,LC)
              WK=EXP(KK_IQ*DTAU)
              MU_KK=MU*KK_IQ
              IF(ABS(1+MU_KK).LT..0001)THEN

!                 L'Hospital limit:
                  ALB_UP=ALB_UP+GU(IU,IQ,LC)*LL(IQ,LC)*EXP1*DTAU/MU
              ELSE
                  ALB_UP=ALB_UP                                         &
     &              +GU(IU,IQ,LC)*LL(IQ,LC)*(EXP2*WK-EXP1)/(1+MU_KK)
              ENDIF

!             MU IS POSITIVE, SO DENOMINATOR IS NOT NEAR ZERO.
              ALB_UP=ALB_UP+                                            &
     &          GU(IU,JQ,LC)*LL(JQ,LC)*(EXP2-EXP1*WK)/(1-MU_KK)
              JQ=JQ-1
          ENDDO
          LC=LCM1
      ENDDO

!     CALCULATE CONTRIBUTION FROM USER OUTPUT
!     LEVEL TO NEXT COMPUTATIONAL LEVEL:
      IF(ABS(DTAU2).LT.1.E-6)RETURN
      DTAU1=TAUTOP-TAUCPR(LAYTOP-1)
      JQ=NSTR
      DO IQ=1,NN

!         KK_IQ NEGATIVE:
          KK_IQ=KK(IQ,LAYTOP)
          WK=EXP(-KK_IQ*DTAU2)
          MU_KK=MU*KK_IQ
          IF(ABS(1+MU_KK).LT..0001)THEN
              ALB_UP=ALB_UP-GU(IU,IQ,LAYTOP)*LL(IQ,LAYTOP)*EXP2*DTAU2/MU
          ELSE
              ALB_UP=ALB_UP                                             &
     &          +GU(IU,IQ,LAYTOP)*LL(IQ,LAYTOP)*(WK-EXP2)/(1+MU_KK)
          ENDIF

!         MU IS POSITIVE, SO DENOMINATOR IS NOT NEAR ZERO:
          ALB_UP=ALB_UP+GU(IU,JQ,LAYTOP)*LL(JQ,LAYTOP)                  &
     &      *EXP(KK_IQ*DTAU1)*(1-WK*EXP2)/(1-MU_KK)
          JQ=JQ-1
      ENDDO
      RETURN
      END