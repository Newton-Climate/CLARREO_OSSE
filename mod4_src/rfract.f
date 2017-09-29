      REAL FUNCTION RFRACT(FREQ,PRATIO,TEMP,WH2O,WCO2)

!     RFRACT RETURNS THE REFRACTIVITY, DEFINED AS ONE MINUS THE REAL
!     PART OF THE INDEX OF REFRACTION.  DERIVED FROM BODHAINE, ET AL.,
!     "Note on Rayleigh Optical Depth Calculations," Submitted to
!     JTECH April 26, 1999. (bbodhaine@cmdl.noaa.gov).

!     INPUT ARGUMENTS:
!       FREQ     SPECTRAL FREQUENCY [CM-1].
!       PRATIO   ATMOSPHERIC PRESSURE [ATM].
!       TEMP     ATMOSPHERIC TEMPERATURE [K].
!       WH2O     ATMOSPHERIC WATER VAPOR DENSITY [GM / M3].
!       WCO2     ATMOSPHERIC CARBON DIOXIDE DENSITY [PPMV].
      REAL FREQ,PRATIO,TEMP,WH2O,WCO2

!     PARAMETERS:
!       CCO2     COEFFICIENT FOR CO2 DENSITY CORRECTION.
!       AH2O     FREQUENCY INDEPENDENT H2O DENSITY CORRECTION.
!       BH2O     FREQUENCY DEPENDENT H2O DENSITY CORRECTION.
!       DRY0     DRY AIR REFRACTIVITY EXPANSION PARAMETER.
!       DRY1     DRY AIR REFRACTIVITY EXPANSION PARAMETER.
!       DRY2     DRY AIR REFRACTIVITY EXPANSION PARAMETER.
!       DRY3     DRY AIR REFRACTIVITY EXPANSION PARAMETER.
!       DRY4     DRY AIR REFRACTIVITY EXPANSION PARAMETER.
      REAL CCO2,AH2O,BH2O,DRY0,DRY1,DRY2,DRY3,DRY4
      PARAMETER(CCO2=5.40335E-07,AH2O=1.9809E-10,BH2O=3.1759E-21,       &
     &  DRY0=.0232226,DRY1=7.147815E+08,DRY2=1.32274E+10,               &
     &  DRY3=5.029045E+06,DRY4=3.932957E+09)

!     LOCAL VARIABLES
!       FREQ2    SPECTRAL FREQUENCY SQUARED [CM-2].
      REAL FREQ2

!     REFRACTIVITY:
      FREQ2=FREQ*FREQ
      RFRACT=PRATIO*(DRY0+DRY1/(DRY2-FREQ2)+DRY3/(DRY4-FREQ2))/TEMP     &
     &  *(1.+CCO2*WCO2)-(AH2O-BH2O*FREQ2)*WH2O*TEMP
      RETURN
      END