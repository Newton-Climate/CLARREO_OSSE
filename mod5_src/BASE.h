
!     /BASE/
!       MPOINT   ARRAY MAPPING MOLECULAR SPECIES (MOLSPC) TO THEIR
!                LOCATION IN THE COLUMN AMOUNT ARRAYS (E.G. WTOTAL).
!       WAVLEN   AEROSOL OPTICAL PROPERTIES WAVELENGTH GRID [MICRONS].
!       TX       SPECIES TRANSMITTANCES AND OPTICAL DEPTH.
      INTEGER MPOINT
      REAL AWCCON,WAVLEN,EXTC,ABSC,ASYM,TX,RELHUM
      COMMON/BASE/MPOINT(MMOLT),AWCCON(MAER),WAVLEN(MXWVLN),            &
     &  EXTC(MAER,MXWVLN),ABSC(MAER,MXWVLN),ASYM(MAER,MXWVLN),          &
     &  TX(MEXTXY),RELHUM(LAYDIM)