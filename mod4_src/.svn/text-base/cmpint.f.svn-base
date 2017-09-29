      SUBROUTINE  CMPINT(FBEAM, GC, KK, LAYRU, LL, LYRCUT, MAZIM,       &
     &                    MXCMU, MXUMU, NCUT, NN, NSTR, PLANK, NTAU,    &
     &                    TAUCPR, UMU0, UTAUPR, ZZ, ZPLK0, ZPLK1,       &
     &                    UUM)

!          Calculates the Fourier intensity components at the quadrature
!              angles for azimuthal expansion terms (MAZIM) in EQ. SD(2)

!                  I N P U T    V A R I A B L E S:

!       KK      :  eigenvalues of coeff. matrix in EQ. SS(7)
!       GC      :  eigenvectors at polar quadrature angles in EQ. SC(1)
!       LL      :  constants of integration in EQ. SC(1), obtained
!                  by solving scaled version of EQ. SC(5);
!                  exponential term of EQ. SC(12) not included
!       LYRCUT  :  logical flag for truncation of computational layer
!       MAZIM   :  order of azimuthal component
!       NCUT    :  number of computational layer where absorption
!                  optical depth exceeds -ABSCUT-
!       NN      :  order of double-Gauss quadrature (NSTR/2)
!       TAUCPR  :  cumulative optical depth (delta-M-scaled)
!       UTAUPR  :  optical depths of user output levels in delta-M
!                  coordinates;  equal to -UTAU- if no delta-M
!       ZZ      :  beam source vectors in EQ. SS(19)
!       ZPLK0   :  thermal source vectors -Z0-, by solving EQ. SS(16)
!       ZPLK1   :  thermal source vectors -Z1-, by solving EQ. SS(16)
!       (Remainder are 'DisORT' input variables)

!                  O U T P U T   V A R I A B L E S:

!       UUM     :  Fourier components of the intensity in EQ. SD(12)
!                  (at polar quadrature angles)

!                  I N T E R N A L   V A R I A B L E S:

!       FACT    :  EXP(- UTAUPR / UMU0)
!       ZINT    :  intensity of M=0 case, in EQ. SC(1)
!+----------------------------------------------------------------------

         LOGICAL  LYRCUT, PLANK
         INTEGER  LAYRU(*)
         REAL     UUM(MXUMU,*)
         REAL     GC(MXCMU,MXCMU,*), KK(MXCMU,*), LL(MXCMU,*),          &
     &          TAUCPR(0:*), UTAUPR(*), ZZ(MXCMU,*),                    &
     &          ZPLK0(MXCMU,*), ZPLK1(MXCMU,*)

!                                                  Loop over user levels
         DO 40 LU=1, NTAU

            LYU=LAYRU(LU)
            IF (LYRCUT .AND. LYU.GT.NCUT)  GO TO 40

            DO 30 IQ=1, NSTR
               ZINT=0.0
               DO 10 JQ=1, NN
                 ZINT=ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *               &
     &                    EXP(- KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU)))
10           CONTINUE
               DO 20 JQ=NN+1, NSTR
                  ZINT=ZINT + GC(IQ,JQ,LYU) * LL(JQ,LYU) *              &
     &                  EXP(- KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU-1)))
20           CONTINUE

               UUM(IQ,LU)=ZINT
               IF (FBEAM.GT.0.0)                                        &
     &            UUM(IQ,LU)=ZINT +                                     &
     &                         ZZ(IQ,LYU) * EXP(- UTAUPR(LU) / UMU0)
               IF (PLANK .AND. MAZIM.EQ.0)                              &
     &            UUM(IQ,LU)=UUM(IQ,LU) + ZPLK0(IQ,LYU) +               &
     &                         ZPLK1(IQ,LYU) * UTAUPR(LU)
30        CONTINUE

40    CONTINUE

        RETURN
        END
