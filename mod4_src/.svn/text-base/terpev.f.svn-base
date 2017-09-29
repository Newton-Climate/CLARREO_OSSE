      SUBROUTINE  TERPEV(CWT, EVECC, GL, GU, MAZIM, MXCMU, MXUMU,       &
     &                    NN, NSTR, NUMU, WK, YLMC, YLMU)

!         INTERPOLATE EIGENVECTORS TO USER ANGLES; EQ SD(8)

      REAL  CWT(*), EVECC(MXCMU,*), GL(0:*), GU( MXUMU,*), WK(*),       &
     &      YLMC( 0:MXCMU,*), YLMU( 0:MXCMU,*)

      DO 50  IQ=1, NSTR

         DO 20  L=MAZIM, NSTR-1
!                                       ** INNER SUM IN SD(8) TIMES ALL
!                                   ** FACTORS IN OUTER SUM BUT PLM(MU)
            SUM=0.0
            DO 10  JQ=1, NSTR
               SUM=SUM + CWT(JQ) * YLMC(L,JQ) * EVECC(JQ,IQ)
10          CONTINUE
            WK(L+1)=0.5 * GL(L) * SUM
20       CONTINUE
!                                    ** FINISH OUTER SUM IN SD(8)
!                                    ** AND STORE EIGENVECTORS
         DO 40  IU=1, NUMU
            SUM=0.
            DO 30  L=MAZIM, NSTR-1
               SUM=SUM + WK(L+1) * YLMU(L,IU)
30          CONTINUE
            IF (IQ.LE.NN)  GU(IU, IQ+NN    )=SUM
            IF (IQ.GT.NN)  GU(IU, NSTR+1-IQ)=SUM
40       CONTINUE

50    CONTINUE

      RETURN
      END
