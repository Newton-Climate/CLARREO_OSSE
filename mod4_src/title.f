      BLOCK DATA TITLE

!     TITLE INFORMATION
      CHARACTER*20 HHAZE(16),HSEASN(2),HVULCN(8),HMET(2),HMODEL(8)
      CHARACTER*26 HTRRAD(4)
      COMMON/TITL/HHAZE,HSEASN,HVULCN,HMET,HMODEL,HTRRAD
      REAL VSB
      COMMON/VSBD/VSB(10)
      SAVE /TITL/,/VSBD/
      DATA VSB/23.,5.,0.,23.,5.,50.,23.,0.2,0.5,0./
      DATA HHAZE/              'RURAL               ',                  &
     &  'RURAL               ','NAVY MARITIME       ',                  &
     &  'MARITIME            ','URBAN               ',                  &
     &  'TROPOSPHERIC        ','USER DEFINED        ',                  &
     &  'FOG1 (ADVECTTION)   ','FOG2(RADIATI0N)     ',                  &
     &  'DESERT AEROSOL      ','BACKGROUND STRATO   ',                  &
     &  'AGED VOLCANIC       ','FRESH VOLCANIC      ',                  &
     &  'AGED VOLCANIC       ','FRESH VOLCANIC      ',                  &
     &  'METEORIC DUST       '/
      DATA HSEASN/                                                      &
     &  'SPRING-SUMMER       ','FALL-WINTER         ' /
      DATA HVULCN/                                                      &
     &  'BACKGROUND STRATO   ','MODERATE VOLCANIC   ',                  &
     &  'HIGH     VOLCANIC   ','HIGH     VOLCANIC   ',                  &
     &  'MODERATE VOLCANIC   ','MODERATE VOLCANIC   ',                  &
     &  'HIGH     VOLCANIC   ','EXTREME  VOLCANIC   '/
      DATA HMET/                                                        &
     &  'NORMAL              ','TRANSITION          '/
      DATA HMODEL /                                                     &
     &  'TROPICAL MODEL      ','MIDLATITUDE SUMMER  ',                  &
     &  'MIDLATITUDE WINTER  ','SUBARCTIC   SUMMER  ',                  &
     &  'SUBARCTIC   WINTER  ','1976 U S STANDARD   ',                  &
     &  '                    ','MODEL =0HORIZONTAL  '/
      DATA HTRRAD/                                                      &
     & ' TRANSMITTANCE            ',' RADIANCE                 ',       &
     & ' RADIANCE+SOLAR SCATTERING',' TRANSMITTED SOLAR IRRAD. '/
      END
