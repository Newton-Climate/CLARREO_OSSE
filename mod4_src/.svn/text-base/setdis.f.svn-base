      SUBROUTINE  SETDIS(CMU, CWT, DELTAM, DTAUC, EXPBEA, FBEAM, FLYR,  &
     &                    GL, IBCND, LAYRU, LYRCUT, MAXCOE,             &
     &                    MAXUMU, MXCMU, NCUT, NLYR, NTAU, NN, NSTR,    &
     &                    PLANK, NUMU, ONLYFL, OPRIM, PMOM, SSALB, TAUC,&
     &                    TAUCPR, UTAU, UTAUPR, UMU, UMU0, USRTAU,      &
     &                    USRANG)

!                            Perform miscellaneous setting-up operations

!       Routines called:  ERRMSG, QGAUSN

!       INPUT :  all are DISORT input variables (see DOC file)

!       OUTPUT:  NTAU,UTAU   if USRTAU=FALSE
!                NUMU,UMU    if USRANG=FALSE
!                CMU,CWT     computational polar angles and
!                            corresponding quadrature weights
!                EXPBEA      transmission of direct beam
!                FLYR        truncated fraction in delta-M method
!                GL          phase function Legendre coefficients multi-
!                            plied by (2L+1) and single-scatter albedo
!                HLPR        Legendre moments of surface bidirectional
!                            reflectivity, times 2K+1
!                LAYRU       Computational layer in which UTAU falls
!                LYRCUT      flag as to whether radiation will be zeroed
!                            below layer NCUT
!                NCUT        computational layer where absorption
!                            optical depth first exceeds  ABSCUT
!                NN          NSTR / 2
!                OPRIM       delta-M-scaled single-scatter albedo
!                TAUCPR      delta-M-scaled optical depth
!                UTAUPR      delta-M-scaled version of  UTAU

      LOGICAL  DELTAM, LYRCUT, PLANK, ONLYFL, USRTAU, USRANG
      INTEGER  LAYRU(*)
      REAL     CMU(*), CWT(*), DTAUC(*), EXPBEA(0:*),                   &
     &         FLYR(*), GL(0:MXCMU,*), OPRIM(*),                        &
     &         PMOM(0:MAXCOE,*), SSALB(*), TAUC(0:*),                   &
     &         TAUCPR(0:*), UTAU(*), UTAUPR(*), UMU(*)
      DATA     ABSCUT / 10. /

!                    Set output levels at computational layer boundaries

      IF (.NOT.USRTAU)  THEN
         NTAU=NLYR + 1
         DO 10 LC=0, NTAU-1
10         UTAU(LC+1)=TAUC(LC)
      END IF
!                             Apply delta-M scaling and move description
!                             of computational layers to local variables
      EXPBEA(0)=1.0
      TAUCPR(0)=0.0
      ABSTAU=0.0
      DO 40 LC=1, NLYR
         PMOM(0,LC)=1.0
         IF (ABSTAU.LT.ABSCUT)  NCUT=LC
         ABSTAU=ABSTAU + (1. - SSALB(LC)) * DTAUC(LC)

         IF (.NOT.DELTAM)  THEN
            OPRIM(LC)=SSALB(LC)
            TAUCPR(LC)=TAUC(LC)
            DO 20 K=0, NSTR-1
20            GL(K,LC)=(2*K+1) * OPRIM(LC) * PMOM(K,LC)
            F=0.0
         ELSE
!                                                 delta-M transformation
            F=PMOM(NSTR,LC)
            OPRIM(LC)=SSALB(LC) * (1. - F) / (1. - F * SSALB(LC))
            TAUCPR(LC)=TAUCPR(LC-1) + (1. - F*SSALB(LC)) * DTAUC(LC)
            DO 30 K=0, NSTR-1
30            GL(K,LC)=(2*K+1) * OPRIM(LC) * (PMOM(K,LC)-F) / (1.-F)
         END IF

         FLYR(LC)=F
         EXPBEA(LC)=0.0
         IF (FBEAM.GT.0.0)  EXPBEA(LC)=EXP(- TAUCPR(LC) / UMU0)
40    CONTINUE
!                If no thermal emission, cut off medium below absorption
!               optical depth=ABSCUT (note that delta-M transformation
!                 leaves absorption optical depth invariant).  Not wortH
!                            the trouble for one-layer problems, though.
      LYRCUT=.FALSE.
      IF (ABSTAU.GE.ABSCUT .AND. .NOT.PLANK .AND. IBCND.NE.1            &
     &     .AND. NLYR.GT.1)  LYRCUT =.TRUE.
      IF (.NOT.LYRCUT)  NCUT=NLYR

!                     Set arrays defining location of user output levels
!                               within delta-M-scaled computational mesh
      DO 70 LU=1, NTAU
         DO 50 LC=1, NLYR
            IF (UTAU(LU).GE.TAUC(LC-1) .AND. UTAU(LU).LE.TAUC(LC))      &
     &           GO TO 60
50       CONTINUE
         LC=NLYR

60       UTAUPR(LU)=UTAU(LU)
         IF(DELTAM) UTAUPR(LU)=TAUCPR(LC-1) + (1.-SSALB(LC)*FLYR(LC))   &
     &                                        * (UTAU(LU) - TAUC(LC-1))
         LAYRU(LU)=LC
70    CONTINUE
!                        Calculate computational polar angle cosines and
!                             associated quadrature weights for Gaussian
!                              quadrature on the interval (0,1) (upward)
      NN=NSTR / 2
      CALL  QGAUSN(NN, CMU, CWT)
!                                      Downward (neg) angles and weights
      DO 80 IQ=1, NN
         CMU(IQ+NN)=- CMU(IQ)
         CWT(IQ+NN)=  CWT(IQ)
80    CONTINUE
!                             Compare beam angle to computational angles

      IF (FBEAM.GT.0.0)  THEN
         DO 90 IQ=1, NN
            !WRITE(*,*) "CMU = ",CMU(IQ)  !DRF
            !WRITE(*,*) "UMU0 = ",UMU0 !DRF
            !WRITE(*,*) "HERE" !DRF
            IF (ABS(UMU0-CMU(IQ))/UMU0 .LT. 1.E-4)  CALL ERRMSG         &
     &         ('SETDIS--beam angle=computational angle; change NSTR',  &
     &            .TRUE.)
90       CONTINUE
      END IF
!                  Set output polar angles to computational polar angles

      IF (.NOT.USRANG .OR. (ONLYFL .AND. MAXUMU.GE.NSTR))  THEN
!m
!m       ALWAYS SKIP THIS BRANCH FOR MODTRAN:
!m       (This use to happen automatically since USRANG is .TRUE. and
!m       since MAXUMU equaled 1 and NSTR was at least 2; to accommodate
!m       IBCND=1 calculations, MAXUMU has been increased to 2).
!m       NUMU=NSTR
!m       DO 100 IU=1, NN
!m100       UMU(IU)=- CMU(NN+1-IU)
!m       DO 110 IU=NN+1, NSTR
!m110       UMU(IU)=CMU(IU-NN)
      END IF
!                             Shift positive user angle cosines to upper
!                         locations and put negatives in lower locations

      IF (USRANG .AND. IBCND.EQ.1)  THEN
         DO 120 IU=1, NUMU
120         UMU(IU+NUMU)=UMU(IU)
         DO 130 IU=1, NUMU
130         UMU(IU)=- UMU(2*NUMU+1-IU)
         NUMU=2*NUMU
      END IF

      RETURN
      END
