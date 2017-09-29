      SUBROUTINE  INTCOR(EPSIL, FBEAM, FLYR, LAYRU, LYRCUT, MAXCOE,     &
     &                    MAXULV, MAXUMU, NCOEF, NCUT, NPHI, NSTR,      &
     &                    NTAU, NUMU, OPRIM, PHASA, PHASE, PHASM,       &
     &                    PHIRAD, PI, PMOM, SSALB, TAUC, TAUCPR,        &
     &                    UMU, UMU0, UTAU, UTAUPR, UU)

!             Correct intensity field by using Nakajima-Tanaka algorithm
!         (1988).  For more details, see Section 3.6 of STW NASA report.

!                I N P U T   V A R I A B L E S

!        EPSIL   10 times machine precision
!        FBEAM   incident beam radiation at top
!        FLYR    truncated fraction in delta-M method.
!        LAYRU   index of UTAU in multi-layered system
!        LYRCUT  logical flag for truncation of computational layer
!        NCOEF   number of phase function Legendre coefficients supplied
!        NCUT    total number of computational layers considered
!        NPHI    number of user azimuthal angles
!        NSTR    number of polar quadrature angles
!        NTAU    number of user-defined optical depths
!        NUMU    number of user polar angles
!        OPRIM   delta-M-scaled single-scatter albedo
!        PHIRAD  azimuthal angles
!        PMOM    phase function Legendre coefficients (K, LC)
!                    K=0 to NCOEF, LC=1 to NLYR with PMOM(0,LC)=1
!        SSALB   single scattering albedo at computational layers
!        TAUC    optical thickness at computational levels
!        TAUCPR  delta-M-scaled optical thickness
!        UMU     cosine of emergent angle
!        UMU0    cosine of incident zenith angle
!        UTAU    user defined optical depths
!        UTAUPR  delta-M-scaled version of UTAU

!                O U T P U T   V A R I A B L E S

!        UU      corrected intensity field; UU(IU,LU,J)
!                          IU=1,NUMU; LU=1,NTAU; J=1,NPHI

!                I N T E R N A L   V A R I A B L E S

!        CS      cosine of scattering angle
!        IMTHD   method flag:  1: MS, 2: TMS, 3, IMS
!        PHASA   actual phase function
!        PHASM   delta-M-scaled phase function
!        PHASE   scratch array, phase function multiplied by
!                various factors

!   Routines called:  SECSCA, SINSCA

      LOGICAL  LYRCUT
      INTEGER  LAYRU(*)
      REAL     FLYR(*), OPRIM(*), PHASA(*), PHASE(*),                   &
     &         PHASM(*), PHIRAD(*), PMOM(0:MAXCOE,*),                   &
     &         SSALB(*), TAUC(0:*), TAUCPR(0:*), UMU(*),                &
     &         UTAU(*), UTAUPR(*), UU(MAXUMU,MAXULV,*)

!                                   Loop over polar and azimuthal angles

      DO 200  IU=1, NUMU
!                          Method of correction (STW sec 3.6 for detail)
         IMTHD=0
         IF (UMU(IU).GT.0.)  IMTHD=2
         IF (UMU(IU).LT.0.)  THEN
            IMTHD=3
!                                              set aureole +/- 10 degree
            IF (ABS(UMU(IU)+UMU0) .LE. COS(10.))  THEN
               IF (TAUC(NCUT).GT.3. .AND. TAUC(NCUT).LT.8.)  IMTHD=1
               IF (TAUC(NCUT).GE.8.)  IMTHD=0
            END IF
         END IF
         IF (IMTHD.EQ.0)  GO TO 200

         DO 100  JP=1, NPHI
!                                                  COS(SCATTERING ANGLE)

            CS=-UMU0*UMU(IU) + SQRT((1.-UMU0**2)*(1.-UMU(IU)**2)) *     &
     &                           COS(PHIRAD(JP))

!                                              Initialize phase function
            DO 10  LC=1, NCUT
               PHASA(LC)=1.
               PHASM(LC)=1.
10          CONTINUE
!                                   Initialize Legendre poly. recurrence
            PL1=1.
            PL2=0.
            DO 40  K=1, NCOEF
!                                              Legendre poly. recurrence

               PL =((2*K-1) * CS * PL1 - (K-1) * PL2) / K
               PL2=PL1
               PL1=PL
!                                        Calculate actual phase function
               DO 20  LC=1, NCUT
20             PHASA(LC)=PHASA(LC) + (2*K+1) * PL * PMOM(K,LC)

!                                       Delta-M truncated phase function
               IF(K.LT.NSTR)  THEN
                  DO 30  LC=1, NCUT
30                PHASM(LC)=PHASM(LC) + (2*K+1) * PL *                  &
     &                        (PMOM(K,LC)-FLYR(LC)) / (1.-FLYR(LC))
               END IF
40          CONTINUE
!                                                        ** MS method **
!                                               US & UST of EQ. STW (67)
            IF (IMTHD.EQ.1)  THEN
               DO 50  LC=1, NCUT
50             PHASE(LC)=(1.-FLYR(LC)) * PHASM(LC)
               DO 60  LU=1, NTAU
                  IF (LYRCUT .AND. LAYRU(LU).GT.NCUT)  GO TO 60
                  US =FBEAM/(4.*PI) * SINSCA(EPSIL, LAYRU(LU), NCUT,    &
     &                  PHASA, SSALB, TAUC, UMU(IU), UMU0, UTAU(LU))
                  UTS=FBEAM/(4.*PI) * SINSCA(EPSIL, LAYRU(LU), NCUT,    &
     &                  PHASE, SSALB, TAUC, UMU(IU), UMU0, UTAU(LU))

                  UU(IU,LU,JP)=UU(IU,LU,JP) + US - UTS
60             CONTINUE

               GO TO 100
            END IF
!                                                       ** TMS method **
!                                             UPTS & UPS of EQ. STW (68)
            DO 70  LC=1, NCUT
70          PHASE(LC)=PHASA(LC) / (1.-FLYR(LC)*SSALB(LC))
            DO 80  LU=1, NTAU
               IF (LYRCUT .AND. LAYRU(LU).GT.NCUT)  GO TO 80
               UPTS=FBEAM/(4.*PI) * SINSCA(EPSIL, LAYRU(LU), NCUT,      &
     &                PHASE, SSALB, TAUCPR, UMU(IU), UMU0, UTAUPR(LU))
               UPS =FBEAM/(4.*PI) * SINSCA(EPSIL, LAYRU(LU), NCUT,      &
     &                PHASM, OPRIM, TAUCPR, UMU(IU), UMU0, UTAUPR(LU))

               UU(IU,LU,JP)=UU(IU,LU,JP) + UPTS - UPS
80          CONTINUE
!                                                       ** IMS method **
!                     Correction of secondary scattering below top level

            IF (IMTHD.EQ.3)  THEN
               LTAU=1
               IF (UTAU(1).LE.EPSIL)  LTAU=2
               DO 90 LU=LTAU, NTAU
                  IF (LYRCUT .AND. LAYRU(LU).GT.NCUT)  GO TO 90
                  UHAT=FBEAM/(4.*PI) * SECSCA(CS, FLYR, LAYRU(LU),      &
     &                   MAXCOE, NCOEF, NSTR, PMOM, SSALB, TAUC,        &
     &                   UMU(IU), UMU0, UTAU(LU))

                  UU(IU,LU,JP)=UU(IU,LU,JP) - UHAT
90             CONTINUE
            END IF
!                                                   End of angular loops
100      CONTINUE
200   CONTINUE

      RETURN
      END
