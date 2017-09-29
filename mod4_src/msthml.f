      SUBROUTINE MSTHML(CMU,CWT,GC,GU,KK,LAYRU,LL,LYRCUT,MAXUMU,MXCMU,  &
     &                  MXUMU,NCUT,NN,NSTR,NTAU,NUMU,OPRIM,PI,TAUCPR,   &
     &                  UTAUPR,XR0,XR1,Z0UMS,Z1UMS,ZPLK0,ZPLK1,         &
     &                  FDNTRT,T0CMS)

      INTEGER IQ,IU,JQ,LU,LYU,MAXUMU,MXCMU,MXUMU,NCUT,NN,NSTR,NTAU,NUMU
      REAL PI,ZINT
!    I N P U T     V A R I A B L E S:

!       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       GC       :  EIGENVECTORS AT POLAR QUADRATURE ANGLES, SC(1)
!       GU       :  EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES
!                   (I.E., -G- IN EQ. SC(1))
!       KK       :  EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)
!       LAYRU    :  LAYER NUMBER OF USER LEVEL -UTAU-
!       LL       :  CONSTANTS OF INTEGRATION IN EQ. SC(1), OBTAINED
!                   BY SOLVING SCALED VERSION OF EQ. SC(5);
!                   EXPONENTIAL TERM OF EQ. SC(12) NOT INCLUDED
!       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
!       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
!       NCUT     :  NUMBER OF COMPUTATIONAL LAYER WHERE ABSORPTION
!                     OPTICAL DEPTH EXCEEDS -ABSCUT-
!       OPRIM    :  SINGLE SCATTERING ALBEDO
!       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
!       UTAUPR   :  OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
!                     COORDINATES;  EQUAL TO  -UTAU- IF NO DELTA-M
!       XR0      :  EXPANSION OF THERMAL SOURCE FUNCTION
!       XR1      :  EXPANSION OF THERMAL SOURCE FUNCTION EQS.SS(14-16)
!       ZPLK0    :  THERMAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)
!       ZPLK1    :  THERMAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)

!   I N T E R N A L       V A R I A B L E S:

!       ZINT     :  INTENSITY OF M = 0 CASE, IN EQ. SC(1)

!   O U T P U T    V A R I A B L E S:

!       T0CMS    :  MULTIPLE SCATTERING THERMAL SOURCE FUNCTION
!      FDNTRT    :  DOWNWARD DIFFUSE THERMAL FLUX AT SURFACE

      LOGICAL LYRCUT
      REAL T0CMS(MAXUMU,*)
      INTEGER LAYRU(*)
      REAL T0C,CMU(*),CWT(*),GC(MXCMU,MXCMU,*),GU(MXUMU,MXCMU,*),       &
     &       KK(MXCMU,*),LL(MXCMU,*),OPRIM(*),TAUCPR(0:*),UTAUPR(*),    &
     &       XR0(*),XR1(*),Z0UMS(MXUMU,*),Z1UMS(MXUMU,*),FDNTRT,        &
     &       ZPLK0(MXCMU,*),ZPLK1(MXCMU,*)
!                                               ** LOOP OVER USER LEVELS
      DO 140 LU = 1,NTAU
         LYU = LAYRU(LU)
         IF(.NOT.(LYRCUT.AND.LYU.GT.NCUT))THEN

!                                     ** NO RADIATION REACHES THIS LEVEL

            T0C = (1.-OPRIM(LYU))*(XR0(LYU)+XR1(LYU)*UTAUPR(LU))
!                                ** LOOP OVER USER ANGLES
            DO 130 IU = 1,NUMU
               ZINT = 0.0
               DO 110 JQ = 1,NN
                  ZINT = ZINT+GU(IU,JQ,LYU)*LL(JQ,LYU)                  &
     &                   *EXP(-KK(JQ,LYU)*(UTAUPR(LU)-TAUCPR(LYU)))
  110          CONTINUE
               DO 120 JQ = NN+1,NSTR
                  ZINT = ZINT+GU(IU,JQ,LYU)*LL(JQ,LYU)                  &
     &                   *EXP(-KK(JQ,LYU)*(UTAUPR(LU)-TAUCPR(LYU-1)))
  120          CONTINUE

!                   **  MS SOURCE FUNCTIONS CALCULATED STW(30) M.S. TERM

               T0CMS(IU,LU) = T0C+ZINT+Z0UMS(IU,LYU)+Z1UMS(IU,LYU)      &
     &                        *UTAUPR(LU)
  130       CONTINUE
         ENDIF
  140 CONTINUE
!                                  **  LAYER AVERAGE OF T0CMS AS MODTRAN
      DO 160 IU = 1,NUMU
         DO 150 LU = 1,NTAU-1

            T0CMS(IU,LU) = (T0CMS(IU,LU)+T0CMS(IU,LU+1))/2.
  150    CONTINUE
  160 CONTINUE
      LU = NTAU

      LYU = LAYRU(LU)

      FDNTRT = 0.0

      DO 190 IQ = 1,NN

         ZINT = 0.0
         DO 170 JQ = 1,NN
            ZINT = ZINT+GC(IQ,JQ,LYU)*LL(JQ,LYU)                        &
     &             *EXP(-KK(JQ,LYU)*(UTAUPR(LU)-TAUCPR(LYU)))
  170    CONTINUE
         DO 180 JQ = NN+1,NSTR
            ZINT = ZINT+GC(IQ,JQ,LYU)*LL(JQ,LYU)                        &
     &             *EXP(-KK(JQ,LYU)*(UTAUPR(LU)-TAUCPR(LYU-1)))
  180    CONTINUE
         FDNTRT = FDNTRT+CWT(NN+1-IQ)*CMU(NN+1-IQ)                      &
     &            *(ZINT+ZPLK0(IQ,LYU)+ZPLK1(IQ,LYU)*UTAUPR(LU))
  190 CONTINUE

      FDNTRT = 2.0*PI*FDNTRT
      IF(FDNTRT.LT.0.)FDNTRT = 0.
      RETURN
      END
