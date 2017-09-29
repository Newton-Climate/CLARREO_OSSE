      SUBROUTINE M3DWRT(SMOL,I3DGEN,I3DCNT,I3DDAT,IHAZE,ICLD,M1,        &
     &  IBINMN,IBINMX,BNDWID,JNTRVL,SUBINT,IK,IKMAX,NATM,PRES,DENH2O,   &
     &  DENCO2,DENO3,SAER,XAER,GAER,SCLD,XCLD,GCLD,SRAIN,XRAIN,         &
     &  GRAIN,XMOL,LASCII,LFIRST,M3DGEN,M3DCNT,M3DDAT,M3DMOL,LABEL)

!     M3DWRT WRITES OUT DATA FILES FOR THE MOD3D PROGRAM.
      IMPLICIT NONE

!     ARGUMENTS:
!       I3DGEN   GENERAL INFORMATION FILE UNIT NUMBER FOR MOD3D.
!       I3DCNT   CONTINUA DATA FILE UNIT NUMBER FOR MOD3D.
!       I3DDAT   MOLECULAR EXTINCTION DATA FILE UNIT NUMBER FOR MOD3D.
!       IHAZE    INDEX FOR BOUNDARY LAYER AEROSOL TYPE.
!       ICLD     INDEX FOR CLOUD TYPE.
!       M1       INDEX FOR TEMPERATURE PROFILE.
!       IBINMN   BIN NUMBER OF MINIMUM COMPUTATION SPECTRAL POINT.
!       IBINMX   BIN NUMBER OF MAXIMUM COMPUTATION SPECTRAL POINT.
!       BNDWID   SPECTRAL WIDTH OF BAND MODEL INTERVAL [CM-1].
!       JNTRVL   NUMBER OF K-DISTRIBUTION INTERVALS FOR CURRENT FREQ.
!       SUBINT   SPECTRAL BIN "K" SUB-INTERVAL FRACTIONAL WIDTHS.
!       IK       CURRENT PATH SEGMENT.
!       IKMAX    NUMBER OF PATH SEGMENTS.
!       NATM     COUNTER FOR NUMBER OF DIFFERENT MOLECULAR ATMOSPHERES.
!       PRES     AVERAGE LAYER PRESSURE [ATM].
!       DENH2O   H2O AVERAGE LAYER DENSITY PROFILE [GM/M3].
!       DENCO2   CO2 AVERAGE LAYER DENSITY PROFILE [GM/M3].
!       DENO3    O3 AVERAGE LAYER DENSITY PROFILE [GM/M3].
!       SAER     AEROSOL SEGMENT SPECTRAL SCATTERING OPTICAL DEPTH
!                NORMALIZED BY THE VALUE AT 550 NM.
!       XAER     AEROSOL SEGMENT SPECTRAL EXTINCTION OPTICAL DEPTH
!                NORMALIZED BY THE VALUE AT 550 NM.
!       GAER     AEROSOL SEGMENT SPECTRAL ASYMMETRY FACTOR.
!       SCLD     CLOUD SEGMENT SPECTRAL SCATTERING
!                COEFFICIENT [KM-1 M3/GM].
!       XCLD     CLOUD SEGMENT SPECTRAL EXTINCTION
!                COEFFICIENT [KM-1 M3/GM].
!       GCLD     CLOUD SEGMENT SPECTRAL ASYMMETRY FACTOR.
!       SRAIN    RAIN SEGMENT SPECTRAL SCATTERING
!                COEFFICIENT [KM-1 (HR/MM)**0.63].
!       XRAIN    RAIN SEGMENT SPECTRAL EXTINCTION
!                COEFFICIENT [KM-1 (HR/MM)**0.63].
!       GRAIN    RAIN SEGMENT SPECTRAL ASYMMETRY FACTOR.
!       SMOL     RAYLEIGH OPTICAL DEPTH TIMES TEMPERATURE AND OVER THE
!                PRODUCT OF PRESSURE AND PATH LENGTH [KM-1 K/ATM].
!       XMOL     K SUB-INTERVAL MOLECULAR EXTINCTION OPTICAL DEPTH
!                TIMES TEMPERATURE AND OVER THE PRODUCT OF PRESSURE
!                AND PATH LENGTH [KM-1 K/ATM].
!       LASCII   ASCII FILE OUTPUT FLAG FOR MOD3D FILES.
!       LFIRST   LOGICAL INDICATING FIRST ENTRY WITH NEW ATMOSPHERE.
!       M3DGEN   NAME OF GENERAL INFORMATION FILE FOR MOD3D.
!       M3DCNT   NAME OF CONTINUA DATA FILE FOR MOD3D.
!       M3DDAT   NAME OF MOLECULAR EXTINCTION FILE SORTED BY ATMOSPHERE.
!       M3DMOL   NAME OF MOLECULAR EXTINCTION DATA FILE FOR MOD3D.
      INTEGER I3DGEN,I3DCNT,I3DDAT,IHAZE,ICLD,M1,IBINMN,IBINMX,JNTRVL,  &
     &  IK,IKMAX,NATM
      REAL BNDWID,SUBINT(*),PRES(*),DENH2O(*),DENCO2(*),                &
     &  DENO3(*),SAER(*),XAER(*),GAER(*),SCLD(*),XCLD(*),GCLD(*),       &
     &  SRAIN(*),XRAIN(*),GRAIN(*),SMOL(*),XMOL(*)
      LOGICAL LASCII,LFIRST
      CHARACTER M3DGEN*(*),M3DCNT*(*),M3DDAT*(*),M3DMOL*(*),LABEL*6
!     FUNCTIONS:
!       IRECLN   RETURNS THE RECORD LENGTH FOR A DIRECT ACCESS FILES
!                CONTAINING A KNOWN NUMBER OF VARIABLES IN EACH RECORD.
      INTEGER IRECLN

!     LOCAL VARIABLES:
!       LRECLN   RECORD LENGTH OF A DIRECT ACCESS FILE.
!       INTRVL   ABSORPTION COEFFICIENT (K) SUB-INTERVAL INDEX.
!       JK       LAYER INDEX.
      INTEGER LRECLN,INTRVL,JK

!     SAVED VARIABLES:
!       H2OSAV   SAVED H2O DENSITY IN FIRST ATMOSPHERIC LAYER [GM/M3].
!       CO2SAV   SAVED CO2 DENSITY IN FIRST ATMOSPHERIC LAYER [GM/M3].
!       O3SAV    SAVED O3 DENSITY IN FIRST ATMOSPHERIC LAYER [GM/M3].
!       NRECC    CONTINUA BINARY FILE RECORD NUMBER.
!       NRECM    MOLECULAR BINARY FILE RECORD NUMBER.
      INTEGER NRECC,NRECM
      REAL H2OSAV,CO2SAV,O3SAV
      SAVE NRECC,NRECM,H2OSAV,CO2SAV,O3SAV

!     DATA:
!       LOPEN    OPEN FILE FLAG, SET TO FALSE ONCE FILES ARE OPENED.
!       FRMT     FORMAT OF MOLECULAR EXTINCTION DATA.
!       CHAZE    BOUNDARY LAYER AEROSOL TYPE.
!       CCLD     CLOUD TYPE.
!       CTEMP    TEMPERATURE PROFILE LABEL.
      LOGICAL LOPEN
      CHARACTER FRMT*18,CHAZE(0:10)*16,CCLD(0:19)*21,CTEMP(0:6)*61
      SAVE LOPEN,FRMT,CHAZE,CCLD,CTEMP
      DATA LOPEN/.TRUE./,FRMT/'( 17(1X,1P,E12.6))'/
      DATA CHAZE/                                                       &
     &  'NONE            ','RURAL           ','RURAL           ',       &
     &  'NAVY MARITIME   ','LOWTRAN MARITIME','URBAN           ',       &
     &  'TROPOSPHERIC    ','USER-DEFINED    ','ADVECTION FOG   ',       &
     &  'RADIATIVE FOG   ','DESERT          '/
      DATA CCLD/'NO CLOUDS            ','CUMULUS CLOUD        ',        &
     &          'ALTOSTRATUS CLOUD    ','STRATUS CLOUD        ',        &
     &          'STRATUS/STRATOCUMULUS','NIMBOSTRATUS CLOUD   ',        &
     &          'STRATUS CLOUD        ','NIMBOSTRATUS CLOUD   ',        &
     &          'NIMBOSTRATUS CLOUD   ','CUMULUS CLOUD        ',        &
     &          'CUMULUS CLOUD        ','USER-DEFINED CLOUD   ',        &
     &          'UNDEFINED            ','UNDEFINED            ',        &
     &          'UNDEFINED            ','UNDEFINED            ',        &
     &          'UNDEFINED            ','UNDEFINED            ',        &
     &          'STANDARD CIRRUS      ','SUB-VISUAL CIRRUS    '/
      DATA CTEMP/                                                       &
     &  'USER-SPECIFIED TEMPERATURE PROFILE                           ',&
     &  'TROPICAL ATMOSPHERE (15dg NORTH LATITUDE) TEMPERATURE PROFILE',&
     &  'MID-LATITUDE SUMMER (45dg NORTH LATITUDE) TEMPERATURE PROFILE',&
     &  'MID-LATITUDE WINTER (45dg NORTH LATITUDE) TEMPERATURE PROFILE',&
     &  'SUB-ARCTIC SUMMER (45dg NORTH LATITUDE) TEMPERATURE PROFILE  ',&
     &  'SUB-ARCTIC WINTER (45dg NORTH LATITUDE) TEMPERATURE PROFILE  ',&
     &  '1976 US STANDARD ATMOSPHERE TEMPERATURE PROFILE              '/

!     FIRST CALL CHECK:
      IF(LOPEN)THEN

!         SET REENTRY AND ASCII/BINARY FLAGS:
          LOPEN=.FALSE.
          LASCII=.FALSE.
!ASCII    LASCII=.true.

!         OPEN MOD3D GENERAL INFORMATION DATA FILE:
          CALL OPNFL(I3DGEN,0,M3DGEN,'UNKNOWN','FORMATTED','M3DWRT')

!         FILE NAMES:
          WRITE(I3DGEN,'((A))')                                         &
     &      'MOD3D General, Continua, and Molecular data file names:',  &
     &      M3DGEN,M3DCNT,M3DMOL

!         AEROSOL, CLOUD, AND SPECTRAL DATA:
          WRITE(I3DGEN,'(3(/I8,3X,A),3(/F8.2,3X,A),/I8,3X,A)')          &
     &      IHAZE,CHAZE(IHAZE),ICLD,CCLD(ICLD),M1,CTEMP(M1),            &
     &      IBINMN*BNDWID,'MINIMUM SPECTRAL FREQUENCY [CM-1]',          &
     &      IBINMX*BNDWID,'MAXIMUM SPECTRAL FREQUENCY [CM-1]',          &
     &      BNDWID,'SPECTRAL STEP SIZE (BAND MODEL WIDTH) [CM-1]',      &
     &      JNTRVL-1,                                                   &
     &      'NUMBER OF K (ABSORPTION COEFFICIENT) SUB-INTERVALS MINUS 1'

!         SPECTRAL BIN "K" SUB-INTERVAL FRACTIONAL WIDTHS.
          WRITE(I3DGEN,'(/A,/(6F10.6))')                                &
     &      'SPECTRAL BIN "K" SUB-INTERVAL FRACTIONAL WIDTHS',          &
     &      (SUBINT(INTRVL),INTRVL=1,JNTRVL)

!         NUMBER OF ATMOSPHERIC LAYERS:
          WRITE(I3DGEN,'(/I12,3X,A)')                                   &
     &      IKMAX-1,'NUMBER OF ATMOSPHERIC LAYERS MINUS ONE'

!         BASELINE (MAXIMUM) PROFILES:
          WRITE(I3DGEN,'(2(/A),/(1P,4E12.5))')                          &
     &      '    PRESSURE       H2O         CO2         O3',            &
     &      '      [ATM]      [GM/M3]     [GM/M3]     [GM/M3]',         &
     &      (PRES(JK),DENH2O(JK),DENCO2(JK),DENO3(JK),JK=1,IKMAX)

!         OPEN MOD3D CONTINUA AND MOLECULAR ABSORPTION DATA FILES:
          IF(LASCII)THEN

!             ASCII DATA FILES:
              CALL OPNFL(I3DCNT,0,M3DCNT,'UNKNOWN','FORMATTED','M3DWRT')
              WRITE(FRMT(2:4),'(I3)')JNTRVL
              CALL OPNFL(I3DDAT,0,M3DDAT,'UNKNOWN','FORMATTED','M3DWRT')
          ELSE

!             BINARY DATA FILES:
              NRECC=0
              NRECM=0
              LRECLN=IRECLN(10*IKMAX,LABEL)
              CALL OPNFL(I3DCNT,LRECLN,M3DCNT,'UNKNOWN','UNFORMATTED',  &
     &          'M3DWRT')
              LRECLN=IRECLN(JNTRVL,LABEL)
              CALL OPNFL(I3DDAT,LRECLN,M3DDAT,'UNKNOWN','UNFORMATTED',  &
     &           'M3DWRT')
          ENDIF

!         SAVE BASELINE MOLECULAR DENSITIES AT SURFACE LAYER:
          H2OSAV=DENH2O(1)
          CO2SAV=DENCO2(1)
          O3SAV=DENO3(1)

!         WRITE OUT MOLECULAR PROFILE SCALE FACTORS:
          NATM=1
          WRITE(I3DGEN,'(/(A40))')                                      &
     &      ' ATM   H2O_SCALE   CO2_SCALE    O3_SCALE',                 &
     &      '   1    1.000000    1.000000    1.000000'
      ELSEIF(LFIRST)THEN

!         CONTINUA DATA FILE CAN BE CLOSED FOR SUBSEQUENT ATMOSPHERES.
          IF(NATM.EQ.1)CLOSE(UNIT=I3DCNT,STATUS='KEEP')

!         WRITE OUT MOLECULAR PROFILE SCALE FACTORS:
          NATM=NATM+1
          WRITE(I3DGEN,'(I4,3F12.6)')                                   &
     &      NATM,DENH2O(1)/H2OSAV,DENCO2(1)/CO2SAV,DENO3(1)/O3SAV
      ENDIF

!     ASCII OR BINARY?
      IF(LASCII)THEN

!         WRITE CONTINUA DATA:
          IF(NATM.EQ.1)WRITE(I3DCNT,'(10(1X,1P,E12.6))')                &
     &      SAER(IK),XAER(IK),GAER(IK),SCLD(IK),XCLD(IK),GCLD(IK),      &
     &      SRAIN(IK),XRAIN(IK),GRAIN(IK),SMOL(IK)

!         WRITE MOLECULAR EXTINCTION DATA:
          WRITE(I3DDAT,FRMT)(XMOL(INTRVL),INTRVL=1,JNTRVL)
      ELSE

!         WRITE CONTINUA DATA:
          IF(NATM.EQ.1 .AND. IK.EQ.IKMAX)THEN
              NRECC=NRECC+1
              WRITE(I3DCNT,REC=NRECC)                                   &
     &          (SAER(JK),XAER(JK),GAER(JK),SCLD(JK),XCLD(JK),GCLD(JK), &
     &          SRAIN(JK),XRAIN(JK),GRAIN(JK),SMOL(JK),JK=1,IKMAX)
          ENDIF

!         WRITE MOLECULAR EXTINCTION DATA:
          NRECM=NRECM+1
          WRITE(I3DDAT,REC=NRECM)(XMOL(INTRVL),INTRVL=1,JNTRVL)
      ENDIF
      RETURN
      END
