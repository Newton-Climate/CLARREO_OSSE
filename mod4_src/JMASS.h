! JMASS.h header file for jmass
!+++JMASS+++Created common blocks to hold the output data+++++++++++++++
!     ARRAY SIZE USED BY JMASS COMMON BLOCK
!       ARRAYSIZE   SIZE OF AN ARRAY USED BY THE COMMON BLOCK IN JMASS
      INTEGER ARRAYSIZE
!true
!true FOR LJMASS .TRUE.
!true PARAMETER (ARRAYSIZE = 3000)
!fals
!fals FOR LJMASS .FALSE.
      PARAMETER (ARRAYSIZE = 1)
      REAL      FREQS(ARRAYSIZE), WAVELENGTHS(ARRAYSIZE),               &
     &          TOTTRANS(ARRAYSIZE), H2OTRANS(ARRAYSIZE),               &
     &          CO2PTRANS(ARRAYSIZE), OZONETRANS(ARRAYSIZE),            &
     &          TRACETRANS(ARRAYSIZE), N2CONTTRANS(ARRAYSIZE),          &
     &          H2OCONTTRANS(ARRAYSIZE), MOLSCATTRANS(ARRAYSIZE),       &
     &          TOTAERTRANS(ARRAYSIZE), HNO3TRANS(ARRAYSIZE),           &
     &          AERABSORP(ARRAYSIZE), INTABS(ARRAYSIZE),                &
     &          CO2TRANS(ARRAYSIZE), COTRANS(ARRAYSIZE),                &
     &          CH4TRANS(ARRAYSIZE), N2OTRANS(ARRAYSIZE),               &
     &          O2TRANS(ARRAYSIZE), NH3TRANS(ARRAYSIZE),                &
     &          NOTRANS(ARRAYSIZE), NO2TRANS(ARRAYSIZE),                &
     &          SO2TRANS(ARRAYSIZE), LOGTOTTRANS(ARRAYSIZE)
      COMMON /OUTTRANS/ FREQS, WAVELENGTHS,                             &
     &                  TOTTRANS, H2OTRANS, CO2PTRANS,                  &
     &                  OZONETRANS, TRACETRANS, N2CONTTRANS,            &
     &                  H2OCONTTRANS, MOLSCATTRANS, TOTAERTRANS,        &
     &                  HNO3TRANS, AERABSORP, INTABS,                   &
     &                  CO2TRANS, COTRANS, CH4TRANS,                    &
     &                  N2OTRANS, O2TRANS, NH3TRANS,                    &
     &                  NOTRANS, NO2TRANS, SO2TRANS,                    &
     &                  LOGTOTTRANS
      REAL      ATMOSRADCM(ARRAYSIZE),ATMOSRADMIC(ARRAYSIZE),           &
     &          SRADCM(ARRAYSIZE), SRADMIC(ARRAYSIZE),                  &
     &          SRADSSCM(ARRAYSIZE), GRRADCM(ARRAYSIZE),                &
     &          GRRADMIC(ARRAYSIZE), GRRADDIRCM(ARRAYSIZE),             &
     &          TRADCM(ARRAYSIZE), TRADMIC(ARRAYSIZE),                  &
     &          INTRAD(ARRAYSIZE), GRRADDIRMIC(ARRAYSIZE)
      COMMON /OUTRADS/ ATMOSRADCM, ATMOSRADMIC, SRADCM,                 &
     &                 SRADMIC, SRADSSCM, GRRADCM, GRRADMIC,            &
     &                 GRRADDIRCM, TRADCM, TRADMIC, INTRAD,             &
     &                 GRRADDIRMIC

      REAL      TRANSIRRADCM(ARRAYSIZE), TRANSIRRADMIC(ARRAYSIZE),      &
     &          SOLARIRRADCM(ARRAYSIZE), SOLARIRRADMIC(ARRAYSIZE),      &
     &          INTTIRRAD(ARRAYSIZE), INTSIRRAD(ARRAYSIZE)
      COMMON /OUTIRRADS/ TRANSIRRADCM, TRANSIRRADMIC,                   &
     &                   SOLARIRRADCM, SOLARIRRADMIC,                   &
     &                   INTTIRRAD, INTSIRRAD
      REAL      AVGTRANS, MINIMUMRAD, MAXIMUMRAD, MINRADVALUE,          &
     &          MAXRADVALUE, BOUNDEMIS, TRABBIK, TRADTAU, TRATLNEW

      COMMON /OUTMISC/ AVGTRANS, MINIMUMRAD, MAXIMUMRAD, MINRADVALUE,   &
     &                 MAXRADVALUE, BOUNDEMIS, TRABBIK,                 &
     &                 TRADTAU, TRATLNEW
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
