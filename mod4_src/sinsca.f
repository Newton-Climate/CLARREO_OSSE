      REAL FUNCTION  SINSCA(EPSIL, LAYRU, NLYR, PHASE, OMEGA, TAU, UMU, &
     &                       UMU0, UTAU)

!                       Singly scattered intensity of EQS. STW (65b,d,e)

!                I N P U T   V A R I A B L E S

!        EPSIL    10 times machine precision
!        LAYRU    index of UTAU in multi-layered system
!        NLYR     number of sublayers
!        PHASE    phase functions of sublayers
!        OMEGA    single scattering albedos of sublayers
!        TAU      optical thicknesses of sublayers
!        UMU      cosine of emergent angle
!        UMU0     cosine of incident zenith angle
!        UTAU     user defined optical depth for output intensity

      REAL        PHASE(*), OMEGA(*), TAU(0:*)
!                                                         Initialization
      SINSCA=0.
      EXP0  =EXP(-UTAU/UMU0)
!                                      EQ. STW (65e): direct transmitted

      IF (ABS(UMU+UMU0).LE.EPSIL)  THEN
         DO 10  LYR=1, LAYRU-1
10       SINSCA=SINSCA + OMEGA(LYR)*PHASE(LYR)*(TAU(LYR)-TAU(LYR-1))
         SINSCA=EXP0 / UMU0 * (SINSCA + OMEGA(LAYRU)*PHASE(LAYRU)       &
     &                                   * (UTAU-TAU(LAYRU-1)))
         RETURN
      END IF
!                                               EQ. STW (65b): reflected
      IF (UMU.GT.0.)  THEN
         DO 20  LYR=LAYRU, NLYR
            EXP1=EXP(- ((TAU(LYR)-UTAU)/UMU + TAU(LYR)/UMU0))
            SINSCA=SINSCA + OMEGA(LYR) * PHASE(LYR) * (EXP0-EXP1)
            EXP0=EXP1
20       CONTINUE
      ELSE
!                                     EQ. STW (65d): diffuse transmitted
         DO 30  LYR=LAYRU, 1, -1
            EXP1=EXP(- ((TAU(LYR-1)-UTAU)/UMU + TAU(LYR-1)/UMU0))
            SINSCA=SINSCA + OMEGA(LYR) * PHASE(LYR) * (EXP0-EXP1)
            EXP0=EXP1
30       CONTINUE
      END IF

      SINSCA=SINSCA / (1. + UMU/UMU0)

      RETURN
      END
