      SUBROUTINE CD2E(TP5_FLAG,TP5_FILE,TP5_LINENUMBER)

!     PROCESS CARD2E INPUTS:
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!       ARGUMENTS
!DRF    TP5_FLAG LOGICAL FLAG IF READING IN TP5 FILE
!DRF    TP5_FILE STRING CONTAINING TP5 FILE INFORMATION
!DRF    TP5_LINENUMBER (I/O) CURRENT LINENUMBER FOR TP5 READING
      INTEGER TP5_LINENUMBER
      LOGICAL TP5_FLAG
      CHARACTER TP5_FILE(1000)*110

!     COMMONS:

!     /USSPC/
      INTEGER NARSPC
      REAL VARSPC
      LOGICAL LARUSS
      COMMON/USSPC/NARSPC(4),VARSPC(4,NWAVLN),LARUSS

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

!     READ CARD2E:
      IF(LARUSS)THEN
          CALL ARUEXA(TP5_FLAG,TP5_FILE,TP5_LINENUMBER)
      ELSE
          NARSPC(1)=0
          NARSPC(2)=0
          NARSPC(3)=0
          NARSPC(4)=0
          IF(IHAZE.EQ.7 .OR. ICLD.EQ.11)THEN
            CALL RDEXA(TP5_FLAG,TP5_FILE,TP5_LINENUMBER)
          ENDIF
      ENDIF

!     RETURN TO DRIVER:
      RETURN
      END
