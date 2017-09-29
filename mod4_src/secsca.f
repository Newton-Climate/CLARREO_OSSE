      REAL FUNCTION  SECSCA(CS, FLYR, LAYRU, MAXCOE, NCOEF, NSTR,       &
     &                       PMOM, SSALB, TAUC, UMU, UMU0, UTAU)

!                          Secondary scattered intensity of EQ. STW (71)

!                I N P U T   V A R I A B L E S

!        CS      cosine of scattering angle
!        FLYR    truncated fraction in Delta-M method
!        LAYRU   index of UTAU in multi-layered system
!        MAXCOE  maximum number of phase function moment coefficients
!        NCOEF   number of phase function Legendre coefficients supplied
!        NSTR    number of polar quadrature angles
!        PMOM    phase function Legendre coefficients (K, LC)
!                K=0 to NCOEF, LC=1 to NLYR, with PMOM(0,LC)=1
!        SSALB   single scattering albedo of computational layers
!        TAUC    cumulative optical depth at computational layers
!        UMU     cosine of emergent angle
!        UMU0    cosine of incident zenith angle
!        UTAU    user defined optical depth for output intensity

      REAL       FLYR(*), PMOM(0:MAXCOE,*), SSALB(*), TAUC(0:*)

!              Weighting optical properties of single scattering albedo,
!          truncated fraction and incident zenith angle of EQ. STW (71b)

      TBAR=UTAU - TAUC(LAYRU-1)
      WBAR=SSALB(LAYRU) * TBAR
      FBAR=FLYR(LAYRU) * WBAR
      DO 10  LYR=1, LAYRU-1
         DTAU=TAUC(LYR) - TAUC(LYR-1)
         TBAR=TBAR + DTAU
         WBAR=WBAR + SSALB(LYR) * DTAU
         FBAR=FBAR + SSALB(LYR) * DTAU * FLYR(LYR)
10    CONTINUE

      FBAR=FBAR / WBAR
      WBAR=WBAR / TBAR
      UMU0P= UMU0 / (1.-FBAR*WBAR)
!                                               Weighting phase function
      PSPIKE=1.
      PHAT=1.
      PL1=1.
      PL2=0.
!                                        Legendre poly. recurrence: L<2N
      DO 20  K=1, NSTR-1
         PL =((2*K-1) * CS * PL1 - (K-1) * PL2) / K
         PL2=PL1
         PL1=PL
         PSPIKE=PSPIKE + (2.*PHAT - PHAT**2) * (2*K+1) * PL
20    CONTINUE
!                                      Legendre poly. recurrence: L>2N-1
      DO 40  K=NSTR, NCOEF
         PL =((2*K-1) * CS * PL1 - (K-1) * PL2) / K
         PL2=PL1
         PL1=PL
         DTAU=UTAU - TAUC(LAYRU-1)
         PHAT=PMOM(K,LAYRU) * SSALB(LAYRU) * DTAU
         DO 30  LYR=1, LAYRU-1
            DTAU=TAUC(LYR) - TAUC(LYR-1)
            PHAT=PHAT + PMOM(K,LYR) * SSALB(LYR) * DTAU
30       CONTINUE
         PHAT=PHAT / (FBAR*WBAR*TBAR)
         PSPIKE=PSPIKE + (2.*PHAT - PHAT**2) * (2*K+1) * PL
40    CONTINUE
!                                                  UHAT of EQ. STW (71a)

      SECSCA=(FBAR*WBAR)**2 / (1.-FBAR*WBAR) * PSPIKE *                 &
     &           XIFUNC(-UMU, UMU0P, UMU0P, UTAU)

      RETURN
      END
