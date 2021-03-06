      REAL FUNCTION GETVIS(AER,GNDALT,SUMMER,VOLCAN)
      IMPLICIT NONE

!     GETVIS RETURNS SURFACE METEOROLOGICAL RANGE (VISIBILITY) FOR THE
!     INPUT AEROSOL MODEL AND VERTICAL COLUMN 550 NM OPTICAL DEPTH.

!     PARAMETERS:
!       VIS1     FIRST  VISIBILITY GRID POINT [ 5KM].
!       VIS2     SECOND VISIBILITY GRID POINT [10KM].
!       VIS3     THIRD  VISIBILITY GRID POINT [23KM].
!       OFF0S    UNDER 6KM SUMMER 550NM DEPTH OFFSET FOR      VIS<VIS1.
!       OFF0W    UNDER 6KM WINTER 550NM DEPTH OFFSET FOR      VIS<VIS1.
!       OFF1S    UNDER 6KM SUMMER 550NM DEPTH OFFSET FOR VIS1<VIS<VIS2.
!       OFF1W    UNDER 6KM WINTER 550NM DEPTH OFFSET FOR VIS1<VIS<VIS2.
!       OFF2S    UNDER 6KM SUMMER 550NM DEPTH OFFSET FOR VIS2<VIS<VIS3.
!       OFF2W    UNDER 6KM WINTER 550NM DEPTH OFFSET FOR VIS2<VIS<VIS3.
!       OFF3S    UNDER 6KM SUMMER 550NM DEPTH OFFSET FOR VIS3<VIS.
!       OFF3W    UNDER 6KM WINTER 550NM DEPTH OFFSET FOR VIS3<VIS.
!       COEF0    UNDER 6KM 550NM DEPTH COEF FOR      VIS<VIS1 [KM].
!       COEF1    UNDER 6KM 550NM DEPTH COEF FOR VIS1<VIS<VIS2 [KM].
!       COEF2    UNDER 6KM 550NM DEPTH COEF FOR VIS2<VIS<VIS3 [KM].
!       COEF3S   UNDER 6KM SUMMER 550NM DEPTH COEF FOR VIS3<VIS [KM].
!       COEF3W   UNDER 6KM WINTER 550NM DEPTH COEF FOR VIS3<VIS [KM].
!       BCKGND   INDEX FOR BACKGROUND STRATOSPHERIC VOLCANIC AEROSOLS.
!       MODER    INDEX FOR MODERATE STRATOSPHERIC VOLCANIC AEROSOLS.
!       HIGH     INDEX FOR HIGH STRATOSPHERIC VOLCANIC AEROSOLS.
!       XTREME   INDEX FOR EXTREME STRATOSPHERIC VOLCANIC AEROSOLS.
!       BCKSUM   ABOVE 6KM SUMMER BACKGROUND 550NM VERTICAL DEPTH.
!       MODSUM   ABOVE 6KM SUMMER MODERATE   550NM VERTICAL DEPTH.
!       HISUM    ABOVE 6KM SUMMER HIGH       550NM VERTICAL DEPTH.
!       XTRSUM   ABOVE 6KM SUMMER EXTREME    550NM VERTICAL DEPTH.
!       BCKWIN   ABOVE 6KM WINTER BACKGROUND 550NM VERTICAL DEPTH.
!       MODWIN   ABOVE 6KM WINTER MODERATE   550NM VERTICAL DEPTH.
!       HIWIN    ABOVE 6KM WINTER HIGH       550NM VERTICAL DEPTH.
!       XTRWIN   ABOVE 6KM WINTER EXTREME    550NM VERTICAL DEPTH.
      INTEGER BCKGND,MODER,HIGH,XTREME
      REAL VIS1,VIS2,VIS3,OFF0S,OFF0W,OFF1S,OFF1W,                      &
     &  OFF2S,OFF2W,OFF3S,OFF3W,COEF0,COEF1,COEF2,COEF3S,COEF3W,        &
     &  BCKSUM,MODSUM,HISUM,XTRSUM,BCKWIN,MODWIN,HIWIN,XTRWIN,          &
     &  VIS4,U6SUM4,U6WIN4,U6SUM3,U6WIN3,U6SUM2,U6WIN2,U6SUM1,U6WIN1
      PARAMETER(  OFF0S=.190806, OFF0W=.170150, COEF0=4.781519,         &
     &  VIS1= 5., OFF1S=.147658, OFF1W=.127002, COEF1=4.997261,         &
     &  VIS2=10., OFF2S=.035200, OFF2W=.014544, COEF2=6.121839,         &
     &  VIS3=23., OFF3S=-.003179,               COEF3S=7.004562,        &
     &            OFF3W=-.010791,               COEF3W=6.704554,        &
     &  BCKGND=1,MODER=2,HIGH=3,XTREME=8,                               &
     &  BCKSUM=.023372,MODSUM=.044187,HISUM=.111495,XTRSUM=.280859,     &
     &  BCKWIN=.014820,MODWIN=.034590,HIWIN=.082906,XTRWIN=.266033,     &
     &  VIS4=331.3,U6SUM4=OFF3S+COEF3S/VIS4,U6WIN4=OFF3W+COEF3W/VIS4,   &
     &             U6SUM3=OFF3S+COEF3S/VIS3,U6WIN3=OFF3W+COEF3W/VIS3,   &
     &             U6SUM2=OFF2S+COEF2 /VIS2,U6WIN2=OFF2W+COEF2 /VIS2,   &
     &             U6SUM1=OFF1S+COEF1 /VIS1,U6WIN1=OFF1W+COEF1 /VIS1)

!     INPUT ARGUMENTS:
!       AER      550 NM VERTICAL COLUMN AEROSOL OPTICAL DEPTH.
!       GNDALT   GROUND ALTITUDE (<6.) [KM].
!       SUMMER   TRUE FOR SPRING/SUMMER.
!       VOLCAN   VOLCANIC AEROSOL MODEL NUMBER.
      REAL AER,GNDALT
      LOGICAL SUMMER
      INTEGER VOLCAN

!     LOCAL VARIABLES:
!       UNDER6   AEROSOL 550NM VERTICAL COLUMN OPTICAL DEPTH, 0 TO 6 KM.
      REAL UNDER6

!     DATA:
!       VOLEXT   VOLCANIC AEROSOL EXTINCTION MODEL:
      INTEGER VOLEXT(0:8)
      DATA VOLEXT/BCKGND,BCKGND,MODER,HIGH,HIGH,MODER,MODER,HIGH,XTREME/

!     BRANCH BASED ON SEASON FOR HIGHER ALTITUDE AEROSOLS:
      IF(SUMMER)THEN

!         BRANCH BASED ON VOLCANIC AEROSOL CONTENT:
          IF(VOLEXT(VOLCAN).EQ.BCKGND)THEN
              UNDER6=6*(AER-BCKSUM)/(6-GNDALT)
          ELSEIF(VOLEXT(VOLCAN).EQ.MODER)THEN
              UNDER6=6*(AER-MODSUM)/(6-GNDALT)
          ELSEIF(VOLEXT(VOLCAN).EQ.HIGH)THEN
              UNDER6=6*(AER-HISUM)/(6-GNDALT)
          ELSEIF(VOLEXT(VOLCAN).EQ.XTREME)THEN
              UNDER6=6*(AER-XTRSUM)/(6-GNDALT)
          ELSE
              STOP 'BAD VOLCAN VALUE'
          ENDIF
          IF(UNDER6.GT.U6SUM2)THEN
              IF(UNDER6.GT.U6SUM1)THEN
                  GETVIS=COEF0/(UNDER6-OFF0S)
              ELSE
                  GETVIS=COEF1/(UNDER6-OFF1S)
              ENDIF
          ELSEIF(UNDER6.GT.U6SUM3)THEN
              GETVIS=COEF2/(UNDER6-OFF2S)
          ELSE
              IF(UNDER6.LT.U6SUM4)STOP                                  &
     &          'INPUT AEROSOL + RAYLEIGH OPTICAL DEPTH TOO LOW.'
              GETVIS=COEF3S/(UNDER6-OFF3S)
          ENDIF
      ELSE
          IF(VOLEXT(VOLCAN).EQ.BCKGND)THEN
              UNDER6=6*(AER-BCKWIN)/(6-GNDALT)
          ELSEIF(VOLEXT(VOLCAN).EQ.MODER)THEN
              UNDER6=6*(AER-MODWIN)/(6-GNDALT)
          ELSEIF(VOLEXT(VOLCAN).EQ.HIGH)THEN
              UNDER6=6*(AER-HIWIN)/(6-GNDALT)
          ELSEIF(VOLEXT(VOLCAN).EQ.XTREME)THEN
              UNDER6=6*(AER-XTRWIN)/(6-GNDALT)
          ELSE
              STOP 'BAD VOLCAN VALUE'
          ENDIF
          IF(UNDER6.GT.U6WIN2)THEN
              IF(UNDER6.GT.U6WIN1)THEN
                  GETVIS=COEF0/(UNDER6-OFF0W)
              ELSE
                  GETVIS=COEF1/(UNDER6-OFF1W)
              ENDIF
          ELSEIF(UNDER6.GT.U6WIN3)THEN
              GETVIS=COEF2/(UNDER6-OFF2W)
          ELSE
              IF(UNDER6.LT.U6WIN4)STOP                                  &
     &          'INPUT AEROSOL + RAYLEIGH OPTICAL DEPTH TOO LOW.'
              GETVIS=COEF3W/(UNDER6-OFF3W)
          ENDIF
      ENDIF
      RETURN
      END
