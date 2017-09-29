      SUBROUTINE SPALTR(CMU,CWT,GC,KK,LL,MXCMU,NLYR,NN,NSTR,TAUCPR,     &
     &  SFLUP,SFLDN)

!       CALCULATES SPHERICAL ALBEDO AND TRANSMISSIVITY FOR THE ENTIRE
!       MEDIUM FROM THE M=0 INTENSITY COMPONENTS
!       (THIS IS A VERY SPECIALIZED VERSION OF 'FLUXES')

!    I N P U T    V A R I A B L E S:

!       CMU     :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       CWT     :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       KK      :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
!       GC      :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
!       LL      :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED
!                  BY SOLVING SCALED VERSION OF EQ. SC(5);
!                  EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED
!       NN      :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
!       (REMAINDER ARE 'DISORT' INPUT VARIABLES)

!    O U T P U T   V A R I A B L E S:

!       SFLUP   :  UP-FLUX AT TOP (EQUIVALENT TO SPHERICAL ALBEDO DUE TO
!                  RECIPROCITY).  FOR ILLUMINATION FROM BELOW IT GIVES
!                  SPHERICAL TRANSMISSIVITY
!       SFLDN   :  DOWN-FLUX AT BOTTOM (FOR SINGLE LAYER
!                  EQUIVALENT TO SPHERICAL TRANSMISSIVITY
!                  DUE TO RECIPROCITY)

!    I N T E R N A L   V A R I A B L E S:

!       UP_INT   :  INTENSITY OF M=0 CASE, IN EQ. SC(1)
!       DN_INT   :  INTENSITY OF M=0 CASE, IN EQ. SC(1)
!+----------------------------------------------------------------------
      REAL CMU(*),CWT(*),GC(MXCMU,MXCMU,*),KK(MXCMU,*),LL(MXCMU,*),     &
     &  TAUCPR(0:*)

      DTAUPR=TAUCPR(NLYR)-TAUCPR(NLYR-1)
      LQ=NSTR
      SFLDN=0.
      SFLUP=0.
      DO IQ=1,NN
          KQ=NSTR
          DN_INT=0.
          UP_INT=0.
          DO JQ=1,NN
              DN_INT=DN_INT+GC(IQ,JQ,NLYR)*LL(JQ,NLYR)                  &
     &          +GC(IQ,KQ,NLYR)*LL(KQ,NLYR)*EXP(-KK(KQ,NLYR)*DTAUPR)
              UP_INT=UP_INT+GC(LQ,KQ,1)*LL(KQ,1)                        &
     &          +GC(LQ,JQ,1)*LL(JQ,1)*EXP(KK(JQ,1)*TAUCPR(1))
              KQ=KQ-1
          ENDDO
          WGT=CWT(LQ-NN)*CMU(LQ-NN)
          SFLDN=SFLDN+WGT*DN_INT
          SFLUP=SFLUP+WGT*UP_INT
          LQ=LQ-1
      ENDDO
      SFLUP=2*SFLUP
      SFLDN=2*SFLDN
      RETURN
      END