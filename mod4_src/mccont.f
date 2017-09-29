      SUBROUTINE MCCONT(IVX,IK,IKMAX,NMWAVE)

!     MCCONT WRITES OUT CONTINUUM DATA FOR MONTE-CARLO SIMULATIONS.
!     (ASSUMES NM BINS AND 1 CM-1 BAND MODEL)

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     INPUT ARGUMENTS:
!       IVX      SPECTRAL FREQUENCY [CM-1].
!       IK       PATH SEGMENT INDEX.
!       IKMAX    NUMBER OF PATH SEGMENTS.
!       NMWAVE   WAVELENGTH BIN OF PREVIOUS DATA [NM].

!     OUTPUT ARGUMENTS:
!       NMWAVE   WAVELENGTH BIN OF CURRENT DATA [NM].
      INTEGER IVX,IK,IKMAX,NMWAVE

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /BASE/
!         TX ARRAY
!           1  AEROSOL SCATTERING DEPTH WEIGHTED ASYMMETRY PARAMETER.
!           2  INCREMENTAL AEROSOL SCATTERING OPTICAL DEPTH.
!           3  TOTAL O2 CONTINUUM TRANSMITTANCE.
!           4  N2 CONTINUUM TRANSMITTANCE.
!           5  TOTAL H2O CONTINUUM TRANSMITTANCE.
!           6  RAYLEIGH MOLECULAR SCATTERED TRANSMITTANCE.
!           7  AEROSOL EXTINCTION.
!           8  TOTAL OZONE CONTINUUM TRANSMITTANCE.
!           9  PRODUCT OF ALL CONTINUUM TRANSMITTANCES EXCEPT O2 & HNO3.
!          10  AEROSOL ABSORPTION.
!          11  HNO3 TRANSMITTANCE.
!          12  MOLECULAR CONTINUUM OPTICAL DEPTH.
!          13  INCREMENTAL AEROSOL EXTINCTION OPTICAL DEPTH.
!          14  TOTAL CONTINUUM OPTICAL DEPTH.
!          15  LAYER RAYLEIGH MOLECULAR SCATTERING OPTICAL DEPTH.
!          16  CIRRUS CLOUD TRANSMITTANCE (ICLD=20 ONLY).
!          64  UV/VIS NO2 TRANSMITTANCE.
!          65  UV/VIS SO2 TRANSMITTANCE.
!          66  INCREMENTAL WATER DROPLET SCATTERING OPTICAL DEPTH.
!          67  INCREMENTAL ICE PARTICLE SCATTERING OPTICAL DEPTH.
!          74  INCREMENTAL STD/SUB-VIS CIRRUS SCATTERING OPTICAL DEPTH.
!          75  INCREMENTAL CLOUD (WATER+ICE) EXTINCTION OPTICAL DEPTH.
!          76  CLOUD ASYMMETRY PARAMETER WEIGHTED BY SCATTERING DEPTH.
!          77  INCREMENTAL RAIN EXTINCTION OPTICAL DEPTH
!          78  INCREMENTAL RAIN SCATTERING OPTICAL DEPTH
!          79  INCREMENTAL RAIN ASYMMETRY PARAMETER.
!          --  MOLECULAR LINE CENTER TRANSMITTANCE  --
!          17=H2O  36=CO2  31=O3   47=N2O  44=CO   46=CH4
!          50=O2   54=NO   56=SO2  55=NO2  52=NH3  11=HNO3
      INCLUDE 'BASE.h'

!     /PATH/
!       QTHETA  COSINE OF PATH ZENITH AT PATH BOUNDARIES.
!       AHT     ALTITUDES AT PATH BOUNDARIES [KM].
!       IHT     ALTITUDES AT PATH BOUNDARIES [M].
!       TPH     TEMPERATURE AT PATH BOUNDARIES [K].
!       IMAP    MAPPING FROM PATH SEGMENT MIDPOINT TO VERTICAL LAYER.
!       LOWAHT  INDEX OF VERTICAL LAYER BOUNDARY AT OR JUST BELOW AHT.
!       FACAHT  ALTITUDE INTERPOLATION FRACTION FOR AHT.
      INTEGER IHT,IMAP,LOWAHT
      REAL QTHETA,AHT,TPH,FACAHT
      COMMON/PATH/QTHETA(LAYTWO),AHT(LAYTWO),IHT(0:LAYTWO),             &
     &  TPH(LAYTWO),IMAP(LAYTWO),LOWAHT(LAYTWO),FACAHT(LAYTWO)

!     LOCAL VARIABLES:
!       NM       CURRENT WAVELENGTH [NM].
!       DWT      SPECTRAL BIN WEIGHT.
      INTEGER NM,I,J
      REAL DWT

!     SAVED VARIABLES:
!       WT       SPECTRAL BIN NORMALIZATION.
!       SMOL     RAYLEIGH SCATTERING PATH SEGMENT OPTICAL DEPTH.
!       SAER     AEROSOL SCATTERING PATH SEGMENT OPTICAL DEPTH.
!       XAER     AEROSOL EXTINCTION PATH SEGMENT OPTICAL DEPTH.
!       GAER     AEROSOL SCATTERING ASYMMETRY FACTOR.
!       SCLD     CLOUD SCATTERING PATH SEGMENT OPTICAL DEPTH.
!       XCLD     CLOUD EXTINCTION PATH SEGMENT OPTICAL DEPTH.
!       GCLD     CLOUD SCATTERING ASYMMETRY FACTOR.
!       SRAIN    RAIN SCATTERING PATH SEGMENT OPTICAL DEPTH.
!       XRAIN    RAIN EXTINCTION PATH SEGMENT OPTICAL DEPTH.
!       GRAIN    RAIN SCATTERING ASYMMETRY FACTOR.
      INTEGER IREC
      REAL WT,SMOL(LAYTWO),SAER(LAYTWO),XAER(LAYTWO),GAER(LAYTWO),      &
     &  SCLD(LAYTWO),XCLD(LAYTWO),GCLD(LAYTWO),                         &
     &  SRAIN(LAYTWO),XRAIN(LAYTWO),GRAIN(LAYTWO),STORE(0:9,1:LAYTWO)
      SAVE WT,SMOL,SAER,XAER,GAER,SCLD,XCLD,GCLD,SRAIN,XRAIN,GRAIN,     &
     &  STORE,IREC
      NM=INT(1.E7/(IVX+.5)+.5)
      IF(NMWAVE.NE.NM)THEN

!         NEW SPECTRAL BIN:
          DWT=1.E7/(NM+.5)-(IVX-.5)
          IF(NMWAVE.NE.0)THEN

!             WRITE PREVIOUS WAVELENGTH BIN DATA:
              SMOL(IK)=SMOL(IK)+DWT*TX(15)
              SAER(IK)=SAER(IK)+DWT*TX(2)
              XAER(IK)=XAER(IK)+DWT*TX(13)
              GAER(IK)=GAER(IK)+DWT*TX(1)
              SCLD(IK)=SCLD(IK)+DWT*(TX(66)+TX(67)+TX(74))
              XCLD(IK)=XCLD(IK)+DWT*TX(75)
              GCLD(IK)=GCLD(IK)+DWT*TX(76)
              SRAIN(IK)=SRAIN(IK)+DWT*TX(78)
              XRAIN(IK)=XRAIN(IK)+DWT*TX(77)
              GRAIN(IK)=GRAIN(IK)+DWT*TX(79)
              IF(IK.EQ.1)WT=WT+DWT
              STORE(0,IK)=SMOL(IK)/WT
              STORE(1,IK)=SAER(IK)/WT
              STORE(2,IK)=XAER(IK)/WT
              IF(SAER(IK).GT.0.)THEN
                  STORE(3,IK)=GAER(IK)/SAER(IK)
              ELSE
                  STORE(3,IK)=0.
              ENDIF
              STORE(4,IK)=SCLD(IK)/WT
              STORE(5,IK)=XCLD(IK)/WT
              IF(SCLD(IK).GT.0.)THEN
                  STORE(6,IK)=GCLD(IK)/SCLD(IK)
              ELSE
                  STORE(6,IK)=0.
              ENDIF
              STORE(7,IK)=SRAIN(IK)/WT
              STORE(8,IK)=XRAIN(IK)/WT
              STORE(9,IK)=GRAIN(IK)/WT
              IF(IREC.EQ.0)THEN
                  IF(IK.EQ.1)THEN
                      WRITE(IPR1,'(4A)')'ALT(M)      MOL_SCAT',         &
     &                  '      AER_SCAT       AER_EXT   AER_G',         &
     &                  '      CLD_SCAT       CLD_EXT   CLD_G',         &
     &                  '     RAIN_SCAT      RAIN_EXT  RAIN_G'
                      WRITE(IPR1,'(I6,5X,I9,A,/I6,1X,1PE13.7,           &
     &                  3(2(1X,1PE13.7),1X,0PF7.5))')                   &
     &                  IHT(0),NMWAVE,' NM',IHT(1),(STORE(I,1),I=0,9)
                  ELSE
                      WRITE(IPR1,'(I6,1X,1PE13.7,3(2(1X,1PE13.7),       &
     &                  1X,0PF7.5))')IHT(IK),(STORE(I,IK),I=0,9)
                  ENDIF
              ENDIF
              IF(IK.EQ.IKMAX)THEN
                  IREC=IREC+1
                  WRITE(IDBOUT,REC=IREC)((STORE(I,J),I=0,9),J=1,IKMAX)
              ENDIF
          ELSE
              IREC=0
          ENDIF
          DWT=1.-DWT
          SMOL(IK)=DWT*TX(15)
          SAER(IK)=DWT*TX(2)
          XAER(IK)=DWT*TX(13)
          GAER(IK)=DWT*TX(1)
          SCLD(IK)=DWT*(TX(66)+TX(67)+TX(74))
          XCLD(IK)=DWT*TX(75)
          GCLD(IK)=DWT*TX(76)
          SRAIN(IK)=DWT*TX(78)
          XRAIN(IK)=DWT*TX(77)
          GRAIN(IK)=DWT*TX(79)
          IF(IK.EQ.IKMAX)THEN
              WT=DWT
              NMWAVE=NM
          ENDIF
      ELSE
          SMOL(IK)=SMOL(IK)+TX(15)
          SAER(IK)=SAER(IK)+TX(2)
          XAER(IK)=XAER(IK)+TX(13)
          GAER(IK)=GAER(IK)+TX(1)
          SCLD(IK)=SCLD(IK)+TX(66)+TX(67)+TX(74)
          XCLD(IK)=XCLD(IK)+TX(75)
          GCLD(IK)=GCLD(IK)+TX(76)
          SRAIN(IK)=SRAIN(IK)+TX(78)
          XRAIN(IK)=XRAIN(IK)+TX(77)
          GRAIN(IK)=GRAIN(IK)+TX(79)
          IF(IK.EQ.IKMAX)WT=WT+1.
      ENDIF
      RETURN
      END