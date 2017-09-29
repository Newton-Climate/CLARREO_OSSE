      SUBROUTINE CD2A(LMODEL,TP5_FLAG,TP5_FILE,TP5_LINENUMBER)

!     PROCESS CARD2A INPUTS:
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
!       LMODEL   FLAG, .TRUE. IF MODEL ATMOSPHERE IS USED.
!DRF    TP5_FLAG LOGICAL FLAG IF READING IN TP5 FILE
!DRF    TP5_FILE STRING CONTAINING TP5 FILE INFORMATION
!DRF    TP5_LINENUMBER (I/O) CURRENT LINENUMBER FOR TP5 READING
      INTEGER TP5_LINENUMBER
      LOGICAL TP5_FLAG
      CHARACTER TP5_FILE(1000)*110
      LOGICAL LMODEL

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /CARD2A/
      INTEGER NCRALT,NCRSPC
      DOUBLE PRECISION CTHIK,CALT
      REAL CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,CWAVLN,CCOLWD,        &
     &  CCOLIP,CHUMID,ASYMWD,ASYMIP

!     /CARD2/
!       IHAZE    BOUNDARY LAYER AEROSOL MODEL NUMBER.
!       ISEASN   SEASON NUMBER (1=SPRING-SUMMER, 2=FALL-WINTER).
!       IVULCN   VOLCANIC AEROSOL MODEL NUMBER.
!       ICSTL    COASTAL AIRMASS MODEL NUMBER.
!       ICLD     CLOUD MODEL NUMBER.
!       IVSA     VERTICAL STRUCTURE ALGORITHM (0=OFF, 1=ON).
!       VIS      SURFACE VISIBILITY (GROUND METEOROLOGICAL RANGE) [KM].
!       WSS      CURRENT WIND SPEED (M/S).
!       WHH      24-HOUR WIND SPEED (M/S).
!       RAINRT   RAIN RATE (MM/HR)
!       LSAP     LOGICAL FLAG FOR SPECTRAL AEROSOL PROFILES INPUT.
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA
      REAL VIS,WSS,WHH,RAINRT
      LOGICAL LSAP
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 &
     &  VIS,WSS,WHH,RAINRT,LSAP

!     CARD2A INPUTS:
      IF(ICLD.EQ.18 .OR. ICLD.EQ.19)THEN

!         CARD 2A MODEL CIRRUS
          IF(LJMASS)THEN
              CALL INITCD('CARD2A')
              IF(CTHIK.LT.0.D0)CTHIK=0.D0
              IF(CALT.LT.0.)CALT=0.D0
              IF(CEXT.LT.0.)CEXT=0.
          ELSE
              IF(.NOT.LMODEL)WRITE(IPR,'(/A,/(19X,A))')                 &
     &          'Warning from CD2A:  When the cirrus option (ICLD=18 '//&
     &          'or ICLD=19) is combined with a user-specified profile',&
     &          ' (MODEL=7 or MODEL=8), the requested cirrus'//         &
     &          ' vertical optical depth (= CTHIK * CALT)',             &
     &          ' generally will not be obtained; this can be'//        &
     &          ' corrected by entering the cirrus cloud',              &
     &          ' profile explicitly using ICLD=1 and NCRALT>0'//       &
     &          ' (CARDS 2A and 2E1).'
              IF (TP5_FLAG) THEN
                READ(IRD,'(3F8.0)')CTHIK,CALT,CEXT
              ELSE
                READ(TP5_FILE(TP5_LINENUMBER),'(3F8.0)')CTHIK,CALT,CEXT
                TP5_LINENUMBER=TP5_LINENUMBER+1
              ENDIF
              IF(CTHIK.LT.0.D0)CTHIK=0.D0
              IF(CALT.LT.0.D0)CALT=0.D0
              IF(CEXT.LT.0.)CEXT=0.
              WRITE(IPR,'(/A,3F8.3)')' CARD 2A *****',CTHIK,CALT,CEXT
          ENDIF
          NCRALT=-99
          NCRSPC=-99
          CWAVLN=-99.
          CCOLWD=-99.
          CCOLIP=-99.
          CHUMID=-99.
          ASYMWD=-99.
          ASYMIP=-99.
      ELSEIF(ICLD.GE.1 .AND. ICLD.LE.10)THEN

!         CARD 2A MODEL CLOUDS
          IF(LJMASS)THEN
              CALL INITCD('CARD2A')
          ELSE
              IF (TP5_FLAG) THEN
                READ(IRD,'(3F8.0,2I4,6F8.0)')CTHIK,CALT,CEXT,NCRALT,    &
     &          NCRSPC,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
              ELSE
                READ(TP5_FILE(TP5_LINENUMBER),'(3F8.0,2I4,6F8.0)')      &
     &          CTHIK,CALT,CEXT,NCRALT,                                 &
     &          NCRSPC,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
                TP5_LINENUMBER=TP5_LINENUMBER+1
              ENDIF
              WRITE(IPR,'(/A,3F8.3,2I4,6F8.3)')                         &
     &          ' CARD 2A *****',CTHIK,CALT,CEXT,NCRALT,NCRSPC,         &
     &          CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
          ENDIF
      ELSE
          ICLD=0
          CTHIK=-99.D0
          CALT=-99.D0
          CEXT=-99.
          NCRALT=-99
          NCRSPC=-99
          CWAVLN=-99.
          CCOLWD=-99.
          CCOLIP=-99.
          CHUMID=-99.
          ASYMWD=-99.
          ASYMIP=-99.
      ENDIF

!     RETURN TO DRIVER:
      RETURN
      END
