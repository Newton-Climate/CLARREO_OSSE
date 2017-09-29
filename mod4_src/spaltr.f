      SUBROUTINE  SPALTR(CMU, CWT, GC, KK, LL, MXCMU, NLYR,             &
     &                    NN, NSTR, TAUCPR, SFLUP, SFLDN)

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

!       ZINT    :  INTENSITY OF M=0 CASE, IN EQ. SC(1)
!+----------------------------------------------------------------------

      REAL  CMU(*), CWT(*), GC(MXCMU,MXCMU,*), KK(MXCMU,*),             &
     &      LL(MXCMU,*), TAUCPR(0:*)

      SFLUP=0.0
      DO 20  IQ=NN+1, NSTR
         ZINT =0.0
         DO 10   JQ=1, NN
            ZINT=ZINT + GC(IQ,JQ,1) * LL(JQ,1) *                        &
     &                    EXP(KK(JQ,1) * TAUCPR(1))
10       CONTINUE
         DO 11  JQ=NN+1, NSTR
            ZINT=ZINT + GC(IQ,JQ,1) * LL(JQ,1)
11       CONTINUE

         SFLUP=SFLUP + CWT(IQ-NN) * CMU(IQ-NN) * ZINT
20    CONTINUE

      SFLDN =0.0
      DO 40  IQ=1, NN
         ZINT  =0.0
         DO 30  JQ=1, NN
             ZINT=ZINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR)
30       CONTINUE
         DO 31  JQ=NN+1, NSTR
             ZINT=ZINT + GC(IQ,JQ,NLYR) * LL(JQ,NLYR) *                 &
     &              EXP(- KK(JQ,NLYR)*(TAUCPR(NLYR) - TAUCPR(NLYR-1)))
31       CONTINUE

         SFLDN=SFLDN + CWT(NN+1-IQ) * CMU(NN+1-IQ) * ZINT
40    CONTINUE

      SFLUP=2.0 * SFLUP
      SFLDN=2.0 * SFLDN

      RETURN
      END
