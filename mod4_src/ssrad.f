      SUBROUTINE SSRAD(IPH,IK,MSOFF,IPATH,V,S0,                         &
     &  PFMOL0,PFMOL2,TMOL,TSNOBS,TSNREF,SUMSSS)

!     SSRAD PERFORMS THE SINGLE SCATTERING SOLAR/LUNAR
!     RADIANCE LAYER SUM AND STORES TRANSMITTED SOLAR/LUNAR
!     IRRADIANCES FOR MULTIPLE SCATTERING CALCULATIONS.

!     INPUT ARGUMENTS:
!       IPH     SCATTERING PHASE FUNCTION SWITCH.
!               =2 FOR LOWTRAN MIE FUNCTIONS
!       IK      LAYER INDEX.
!       MSOFF   LAYER INDEX OFFSET (>0 FOR MULTIPLE SCATTERING PATH).
!       IPATH   PATH TYPE SWITCH.
!               =1 FOR DIRECT SUN TO OBSERVER PATH ONLY
!               =2 FOR SUN TO SCATTERING POINT (MSOFF>0)
!                  FOR SUN TO SCATTERING POINT TO OBSERVER (MSOFF=0)
!               =3 FOR OBSERVER TO SCATTERING POINT LINE-OF-SIGHT PATH
!       V       SPECTRAL FREQUENCY [CM-1].
!       S0      EXTRA-TERRESTRIAL SOURCE FUNCTION [W/(CM2-MICRON)].
!       PFMOL0  MOLECULAR PHASE FUNCTION ISOTROPIC SCATTERING TERM.
!       PFMOL2  MOLECULAR PHASE FUNCTION ANISOTROPIC SCATTERING TERM.
!       TMOL    SCATTERING POINT TO SUN (MOON) MOLECULAR TRANSMITTANCE.
!               [THE TOTAL TRANSMITTANCE IS TMOL*EXP(-TX(14))]
      INTEGER IPH,IK,MSOFF,IPATH
      REAL V,S0,PFMOL0,PFMOL2,TMOL

!     OUTPUT ARGUMENTS:
!       TSNOBS   SOLAR/LUNAR IRRADIANCE AT THE OBSERVER.
!       TSNREF   SOLAR/LUNAR IRRADIANCE ALONG L-SHAPED PATH.
!       SUMSSS   SOLAR/LUNAR SINGLE SCATTERING RADIANCE SUM.
      REAL TSNOBS,TSNREF,SUMSSS
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BASE.h'
      INCLUDE 'SOLS.h'
      INTEGER NCRALT,NCRSPC
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      &
     &  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP

!     /AER/
!     THERE ARE "MAER=17" COMPONENTS:
!      1     AEROSOL 1 (NOMINALLY, BOUNDARY LAYER AEROSOL).
!      2     AEROSOL 2 (NOMINALLY, TROPOSPHERIC AEROSOL).
!      3     AEROSOL 3 (NOMINALLY, STRATOSPHERIC AEROSOL).
!      4     AEROSOL 4 (NOMINALLY, VOLCANIC AEROSOL).
!      5     CIRRUS CLOUD.
!      6     CLOUD 1 (NOMINALLY, WATER CLOUD).
!      7     CLOUD 2 (NOMINALLY, ICE CLOUD).
!     8-17   NOVAM (NAVY OCEANIC VERTICAL AEROSOL MODEL) AEROSOL LAYERS.

!       NAER     NUMBER OF ACTIVE AEROSOLS.
!       EXTV     SPECTRAL EXTINCTION (NORMALIZED TO 1 AT 550 NM).
!       ABSV     SPECTRAL ABSORPTION (1-ABSV/EXTV=SCATTERING ALBEDO).
!       ASYV     SPECTRAL LEGENDRE MOMENT (DIVIDED BY 2N+1).
      INTEGER NAER
      REAL EXTV,ABSV,ASYV
      COMMON/AER/NAER,EXTV(MAER),ABSV(MAER),ASYV(MXCMU,MAER)

!     /MSRD/
!       CSSCAT   COSINE OF THE SCATTERING ANGLE.
!                (AT H1 IF IMULT=1; AT OR "NEAR" H2 IF IMULT=-1)
!       SLEGEN   Nth LEGENDRE POLYNOMIAL EVALUATED AT THE COSINE OF THE
!                SCATTERING ANGLE TIMES (2N+1)/4pi (N=0 TO NSTR-1).
!       CSZEN0   LAYER BOUNDARY COSINE OF SOLAR/LUNAR ZENITH.
!       CSZEN    LAYER AVERAGE COSINE OF SOLAR/LUNAR ZENITH.
!       CSZENX   AVERAGE SOLAR/LUNAR COSINE ZENITH EXITING
!                (AWAY FROM EARTH) THE CURRENT LAYER.
!       BBGRND   THERMAL EMISSION (FLUX) AT THE GROUND [W CM-2 / CM-1].
!       BBNDRY   LAYER BOUNDARY THERMAL EMISSION (FLUX) [W CM-2 / CM-1].
!       TCONT    LAYER CONTINUUM OPTICAL DEPTH.
!       TAUT     LAYER TOTAL OPTICAL DEPTH.
!       GTSCAT   SUM OVER SCATTERING SOURCES OF SCATTERING OPTICAL DEPTH
!                AND PHASE FUNCTION LEGENDRE COEFFICIENT PRODUCTS.
!       COSBAR   LAYER EFFECTIVE SCATTERING ASYMMETRY FACTOR.
!       DEPRAT   FRACTIONAL DECREASE IN WEAK-LINE OPTICAL DEPTH TO SUN.
!       S0DEP    OPTICAL DEPTH FROM LAYER BOUNDARY TO SUN.
!       S0TRN    TRANSMITTED SOLAR IRRADIANCES [W CM-2 / CM-1]
!       UPF      LAYER BOUNDARY UPWARD THERMAL FLUX [W CM-2 / CM-1].
!       DNF      LAYER BOUNDARY DOWNWARD THERMAL FLUX [W CM-2 / CM-1].
!       UPFS     LAYER BOUNDARY UPWARD SOLAR FLUX [W CM-2 / CM-1].
!       DNFS     LAYER BOUNDARY DOWNWARD SOLAR FLUX [W CM-2 / CM-1].
      REAL CSSCAT,SLEGEN,CSZEN0,CSZEN,CSZENX,TCONT,TAUT,GTSCAT,COSBAR,  &
     &  BBGRND,BBNDRY,S0DEP,S0TRN,DEPRAT,UPF,DNF,UPFS,DNFS
      COMMON/MSRD/CSSCAT,SLEGEN(0:MAZ),                                 &
     &  CSZEN0(LAYDIM),CSZEN(LAYDIM),CSZENX(LAYDIM),TCONT(LAYDIM),      &
     &  TAUT(MXKSUB,LAYDIM),GTSCAT(0:MXCMU,1:LAYDIM),COSBAR(LAYDIM),    &
     &  BBGRND,BBNDRY(LAYDIM),S0DEP(MXKSUB,LAYTWO),S0TRN(MXKSUB,LAYTWO),&
     &  DEPRAT(MXKSUB,LAYDIM),UPF(MXKSUB,LAYDIM),DNF(MXKSUB,LAYDIM),    &
     &  UPFS(MXKSUB,LAYDIM),DNFS(MXKSUB,LAYDIM)

!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI.
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

!     /USRDTA/
      INTEGER NAERF,NWLF
      REAL WLF,F
      COMMON/USRDTA/NAERF,NWLF,WLF(MWLF),F(MAERF,MANGLS,MWLF)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:
      REAL PHASEF,HENGNS,EXPINT

!     LOCAL VARIABLES:
!       PANGLO   PHASE ANGLE LOWER INTERPOLATION INDEX FOR ANGF ARRAY.
!       PANGUP   PHASE ANGLE UPPER INTERPOLATION INDEX FOR ANGF ARRAY.
      INTEGER IKP1,PANGLO,PANGUP
      REAL STRANL,DELDEP,PNOV

!     SAVED PHASE FUNCTION VARIABLES:
      REAL PMOLB,PMOLF,PAERB(4),PAERF(4),                               &
     &  PCIRB,PCIRF,PCLDB,PCLDF,PICEB,PICEF,PRAINB,PRAINF
      SAVE PMOLB,PMOLF,PAERB,PAERF,                                     &
     &  PCIRB,PCIRF,PCLDB,PCLDF,PICEB,PICEF,PRAINB,PRAINF

!     SAVE WAVELENGTH INTERPOLATION VARIABLES USED WITH IPH=1:
      INTEGER IWLFM1,IWLF
      REAL WAVFAC
      SAVE IWLFM1,IWLF,WAVFAC

!     DATA:
      LOGICAL LWARN
      SAVE LWARN
      DATA LWARN/.TRUE./

!     BRANCH BASED ON PATH TYPE
      GOTO(10,40,50),IPATH

!     IPATH=1 (DIRECT SUN TO OBSERVER PATH).
   10 CONTINUE

!     STORE THE TRANSMITTED SOLAR/LUNAR IRRADIANCE TO THE
!     OBSERVER AND THEN RETURN IF MULTIPLE SCATTERING PATH.
      IF(TMOL.GT.0.)THEN
          S0DEP(1,1)=TX(14)-LOG(TMOL)
          S0TRN(1,1)=S0*TMOL*EXP(-TX(14))
      ELSE
          S0DEP(1,1)=BIGNUM
          S0TRN(1,1)=0.
      ENDIF
      IF(MSOFF.GT.0)RETURN
      TSNOBS=S0TRN(1,1)

!     INITIALIZE THE SOLAR/LUNAR SINGLE SCATTERING RADIANCE LAYER SUM.
!     STORE THE SUN/MOON TO OBSERVER OPTICAL DEPTH.
      SUMSSS=0.

!     STORE RAYLEIGH (PMOLB), FOUR AEROSOL (PAERB), CIRRUS (PCIRB),
!     CLOUD WATER DROPLET (PCLDB) AND CLOUD ICE PARTICLE (PICEB)
!     PHASE FUNCTIONS AT THE OBSERVER.
      PMOLB=PFMOL0+PFMOL2*PHSCOS(1)**2
      IF(IPH.EQ.0)THEN

!        HENYEY-GREENSTEIN
         PAERB(1)=PHASFN(1,2)
         PAERB(2)=PAERB(1)
         PAERB(3)=PAERB(1)
         PAERB(4)=PAERB(1)
         PCIRB=PAERB(1)
      ELSEIF(IPH.EQ.1)THEN

!        /USRDTA/ USER-SUPPLIED PHASE FUNCTION; DETERMINE ANGLE INDICES:
         WAVFAC=1.E4/(V+1.E-5)
         IWLFM1=1
         IF(WAVFAC.LE.WLF(1))THEN
             IWLF=1
             WAVFAC=0.
         ELSE
             DO 20 IWLF=2,NWLF
                 IF(WAVFAC.LT.WLF(IWLF))THEN
                     WAVFAC=(WLF(IWLF)-WAVFAC)/(WLF(IWLF)-WLF(IWLFM1))
                     GOTO30
                 ENDIF
                 IWLFM1=IWLF
   20        CONTINUE
             IWLF=NWLF
             WAVFAC=0.
   30        CONTINUE
         ENDIF
         PANGUP=IUPPHS(1)
         PANGLO=PANGUP-1
         PAERB(1)=EXPINT(F(1,PANGLO,IWLF),F(1,PANGUP,IWLF),PHSFAC(1))
         PAERB(2)=EXPINT(F(2,PANGLO,IWLF),F(2,PANGUP,IWLF),PHSFAC(1))
         PAERB(3)=EXPINT(F(3,PANGLO,IWLF),F(3,PANGUP,IWLF),PHSFAC(1))
         PAERB(4)=EXPINT(F(4,PANGLO,IWLF),F(4,PANGUP,IWLF),PHSFAC(1))
         IF(IWLFM1.NE.IWLF)THEN
             PAERB(1)=PAERB(1)+WAVFAC*(EXPINT(F(1,PANGLO,IWLFM1),       &
     &         F(1,PANGUP,IWLFM1),PHSFAC(1))-PAERB(1))
             PAERB(2)=PAERB(2)+WAVFAC*(EXPINT(F(2,PANGLO,IWLFM1),       &
     &         F(2,PANGUP,IWLFM1),PHSFAC(1))-PAERB(2))
             PAERB(3)=PAERB(3)+WAVFAC*(EXPINT(F(3,PANGLO,IWLFM1),       &
     &         F(3,PANGUP,IWLFM1),PHSFAC(1))-PAERB(3))
             PAERB(4)=PAERB(4)+WAVFAC*(EXPINT(F(4,PANGLO,IWLFM1),       &
     &         F(4,PANGUP,IWLFM1),PHSFAC(1))-PAERB(4))
         ENDIF
         PCIRB=HENGNS(ASYV(1,5),PHSCOS(1))
      ELSEIF(IPH.EQ.2)THEN

!        MIE DATA-BASE
         PAERB(1)=PHASEF(1,V,AH1(1),PHSANG(1),ARH(1))
         PAERB(2)=PHASEF(2,V,AH1(1),PHSANG(1),ARH(1))
         PAERB(3)=PHASEF(3,V,AH1(1),PHSANG(1),ARH(1))
         PAERB(4)=PHASEF(4,V,AH1(1),PHSANG(1),ARH(1))
         PCIRB=HENGNS(ASYV(1,5),PHSCOS(1))
      ENDIF
      IF(NCRSPC.EQ.1)THEN

!         ANGULAR INTERPOLATE TABULATED SCATTERING PHASE FUNCTIONS:
          CALL GETPF(PHSCOS(1),PCLDB,PICEB)
      ELSE
          IF(ABS(ASYMWD).GE.1.)THEN

!             SPECTRALLY VARYING WATER DROPLET SCATTERING PHASE FUNCTION
              PCLDB=HENGNS(ASYV(1,6),PHSCOS(1))
          ELSE

!             SPECTRALLY INDEPENDENT WATER DROPLET PHASE FUNCTION
              PCLDB=PHASFN(1,3)
          ENDIF
          IF(ABS(ASYMIP).GE.1.)THEN

!             SPECTRALLY VARYING ICE PARTICLE SCATTERING PHASE FUNCTION
              PICEB=HENGNS(ASYV(1,7),PHSCOS(1))
          ELSE

!             SPECTRALLY INDEPENDENT ICE PARTICLE PHASE FUNCTION
              PICEB=PHASFN(1,4)
          ENDIF
      ENDIF
      PRAINB=HENGNS(TX(79),PHSCOS(1))
      RETURN

!     IPATH=2 (SUN TO SCATTERING POINT TO OBSERVER PATH).
   40 CONTINUE

!     STORE THE IPATH=2 TRANSMITTED SOLAR/LUNAR IRRADIANCE
!     RETURN IF MULTIPLE SCATTERING PATH.  CALCULATE THE
!     CHANGE IN AND STORE THE CURRENT IPATH=2 OPTICAL DEPTH.
      IKP1=IK+1
      IF(TMOL.GT.0.)THEN
          S0DEP(1,IKP1)=TX(14)-LOG(TMOL)
          S0TRN(1,IKP1)=S0*TMOL*EXP(-TX(14))
      ELSE

!         DIVIDE BIGNUM BY BOUNDARY NUMBER SO THAT DEPTH
!         TO THE SUN DECREASES WITH ALTITUDE IN MULTIPLE
!         SCATTERING FLUX ADDING ROUTINES BMFLUX AND FLXADD.
          S0DEP(1,IKP1)=BIGNUM/IKP1
          S0TRN(1,IKP1)=0.
      ENDIF
      IF(MSOFF.GT.0)THEN

!         STORE FRACTIONAL DECREASE IN OPTICAL DEPTH ACROSS THE LAYER.
          DEPRAT(1,IK)=0.
          IF(S0DEP(1,IK).LT.S0DEP(1,IKP1))THEN
              IF(LWARN .AND. S0DEP(1,IK).LT..99*S0DEP(1,IKP1))THEN
                  WRITE(IPR,'(/2A,2(/10X,2A,I3,A,1P,E10.3),/(10X,2A))') &
     &              ' WARNING:  WEAK-LINE OPTICAL DEPTH TO',            &
     &              ' THE SUN IS INCREASING WITH ALTITUDE.',            &
     &              ' THE DEPTH TO THE SUN FROM THE BOTTOM',            &
     &              ' OF LAYER',IK,' IS ',S0DEP(1,IK),                  &
     &              ' THE DEPTH TO THE SUN FROM THE TOP   ',            &
     &              ' OF LAYER',IK,' WAS',S0DEP(1,IKP1),                &
     &              ' THIS ANOMALY CAN OCCUR BECAUSE OF',               &
     &              ' CURVED-EARTH EFFECTS OR BECAUSE',                 &
     &              ' OF THE CURTIS-GODSON APPROXIMATION. ',            &
     &              ' THE DEPTH FROM LAYER TOP',                        &
     &              ' WAS DECREASED TO MATCH THE',                      &
     &              ' DEPTH FROM LAYER BOTTOM.',                        &
     &              ' ***  THIS WARNING WILL NOT BE REPEATED  ***'
                  LWARN=.FALSE.
              ENDIF
              S0DEP(1,IKP1)=S0DEP(1,IK)
              S0TRN(1,IKP1)=S0TRN(1,IK)
          ELSEIF(S0DEP(1,IK).GT.0.)THEN
              DEPRAT(1,IK)=1.-S0DEP(1,IKP1)/S0DEP(1,IK)
          ENDIF
          RETURN
      ENDIF
      TSNREF=S0TRN(1,IKP1)

!     STORE RAYLEIGH (PMOLB), FOUR AEROSOL (PAERB), CIRRUS (PCIRB),
!     CLOUD WATER DROPLET (PCLDB) AND CLOUD ICE PARTICLE (PICEB)
!     PHASE FUNCTIONS AT THE CURRENT (FRONT-SIDE) AND PREVIOUS
!     (BACK-SIDE) SCATTERING POINT.
      PMOLF=PMOLB
      PAERF(1)=PAERB(1)
      PAERF(2)=PAERB(2)
      PAERF(3)=PAERB(3)
      PAERF(4)=PAERB(4)
      PCIRF=PCIRB
      PCLDF=PCLDB
      PICEF=PICEB
      PRAINF=PRAINB
      PMOLB=PFMOL0+PFMOL2*PHSCOS(IKP1)**2
      IF(IPH.EQ.0)THEN
         PAERB(1)=PHASFN(IKP1,2)
         PAERB(2)=PAERB(1)
         PAERB(3)=PAERB(1)
         PAERB(4)=PAERB(1)
         PCIRB=PAERB(1)
      ELSEIF(IPH.EQ.1)THEN

!        USER-SUPPLIED F-FUNCTION (SEE COMMON BLOCK USRDTA)
         PANGUP=IUPPHS(IKP1)
         PANGLO=PANGUP-1
         PAERB(1)=EXPINT(F(1,PANGLO,IWLF),F(1,PANGUP,IWLF),PHSFAC(IKP1))
         PAERB(2)=EXPINT(F(2,PANGLO,IWLF),F(2,PANGUP,IWLF),PHSFAC(IKP1))
         PAERB(3)=EXPINT(F(3,PANGLO,IWLF),F(3,PANGUP,IWLF),PHSFAC(IKP1))
         PAERB(4)=EXPINT(F(4,PANGLO,IWLF),F(4,PANGUP,IWLF),PHSFAC(IKP1))
         IF(IWLFM1.NE.IWLF)THEN
             PAERB(1)=PAERB(1)+WAVFAC*(EXPINT(F(1,PANGLO,IWLFM1),       &
     &         F(1,PANGUP,IWLFM1),PHSFAC(IKP1))-PAERB(1))
             PAERB(2)=PAERB(2)+WAVFAC*(EXPINT(F(2,PANGLO,IWLFM1),       &
     &         F(2,PANGUP,IWLFM1),PHSFAC(IKP1))-PAERB(2))
             PAERB(3)=PAERB(3)+WAVFAC*(EXPINT(F(3,PANGLO,IWLFM1),       &
     &         F(3,PANGUP,IWLFM1),PHSFAC(IKP1))-PAERB(3))
             PAERB(4)=PAERB(4)+WAVFAC*(EXPINT(F(4,PANGLO,IWLFM1),       &
     &         F(4,PANGUP,IWLFM1),PHSFAC(IKP1))-PAERB(4))
         ENDIF
         PCIRB=HENGNS(ASYV(1,5),PHSCOS(IKP1))
      ELSEIF(IPH.EQ.2)THEN
         PAERB(1)=PHASEF(1,V,AH1(IKP1),PHSANG(IKP1),ARH(IKP1))
         PAERB(2)=PHASEF(2,V,AH1(IKP1),PHSANG(IKP1),ARH(IKP1))
         PAERB(3)=PHASEF(3,V,AH1(IKP1),PHSANG(IKP1),ARH(IKP1))
         PAERB(4)=PHASEF(4,V,AH1(IKP1),PHSANG(IKP1),ARH(IKP1))
         PCIRB=HENGNS(ASYV(1,5),PHSCOS(IKP1))
      ENDIF
      IF(NCRSPC.EQ.1)THEN

!         ANGULAR INTERPOLATE TABULATED SCATTERING PHASE FUNCTIONS:
          CALL GETPF(PHSCOS(IKP1),PCLDB,PICEB)
      ELSE
          IF(ABS(ASYMWD).GE.1.)THEN

!             SPECTRALLY VARYING WATER DROPLET SCATTERING PHASE FUNCTION
              PCLDB=HENGNS(ASYV(1,6),PHSCOS(IKP1))
          ELSE

!             SPECTRALLY INDEPENDENT WATER DROPLET PHASE FUNCTION
              PCLDB=PHASFN(IKP1,3)
          ENDIF
          IF(ABS(ASYMIP).GE.1.)THEN

!             SPECTRALLY VARYING ICE PARTICLE SCATTERING PHASE FUNCTION
              PICEB=HENGNS(ASYV(1,7),PHSCOS(IKP1))
          ELSE

!             SPECTRALLY INDEPENDENT ICE PARTICLE PHASE FUNCTION
              PICEB=PHASFN(IKP1,4)
          ENDIF
      ENDIF
      PRAINB=HENGNS(TX(79),PHSCOS(IKP1))
      RETURN

!     IPATH=3 (OBSERVER TO SCATTERING POINT).
   50 CONTINUE

!     NO PREPARATION REQUIRED FOR MULTIPLE SCATTERING PATH.
      IF(MSOFF.GT.0)RETURN

!     FOR NOVAM AEROSOLS, THE ASYMMETRY FACTOR IS A FUNCTION OF
!     BOTH WAVELENGTH AND ALTITUDE.  CALCULATE THE SCATTERING
!     PHASE FUNCTION FOR THE LAYER AVERAGE SCATTERING ANGLE COSINE.
      IKP1=IK+1
      IF(TX(72).LE.0.)THEN
          PNOV=0.
      ELSE
          PNOV=HENGNS(TX(73),.5*(PHSCOS(IK)+PHSCOS(IKP1)))
      ENDIF

!                           1
!                           /                                  Z
!     STRANL = S0TRN(1,IK)  |  [ S0TRN(1,IKP1) / S0TRN(1,IK) ]    DZ
!                           /
!                           0
      DELDEP=S0DEP(1,IK)-S0DEP(1,IKP1)
      IF(ABS(DELDEP).LT..001)THEN
          STRANL=.5*(S0TRN(1,IKP1)+S0TRN(1,IK))
      ELSE
          STRANL=(S0TRN(1,IKP1)-S0TRN(1,IK))/DELDEP
      ENDIF

!     SUM THE SOLAR/LUNAR SINGLE SCATTER RADIANCE CONTRIBUTIONS.
      SUMSSS=SUMSSS+STRANL*(TX(72)*PNOV+                                &
     &  .5*(TX(68)*(PAERB(1)+PAERF(1))+TX(69)*(PAERB(2)+PAERF(2))       &
     &     +TX(70)*(PAERB(3)+PAERF(3))+TX(71)*(PAERB(4)+PAERF(4))       &
     &     +TX(15)*(PMOLB   +PMOLF   )+TX(74)*(PCIRB   +PCIRF   )       &
     &     +TX(66)*(PCLDB   +PCLDF   )+TX(67)*(PICEB   +PICEF   )       &
     &     +TX(77)*(PRAINB  +PRAINF  )))
      RETURN
      END