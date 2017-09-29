
!     /SOLS/
!       PHSANG   SOLAR SCATTERING PHASE ANGLE ALONG LINE-OF-SIGHT [DEG].
!       PHSCOS   COSINE OF SOLAR SCATTERING PHASE ANGLE.
!       IUPPHS   PHASE ANGLE UPPER INTERPOLATION INDEX FOR ANGF ARRAY.
!       PHSFAC   PHASE ANGLE INTERPOLATION FRACTION FOR ANGF ARRAY.
!       NANGLS   NUMBER OF ANGLES IN USER-DEFINED PHASE FUNCTIONS.
!       ANGF     GRID OF ANGLES FOR USER-DEFINED PHASE FUNCTIONS [DEG].
!       PHASFN   SPECTRALLY INDEPENDENT SCAT. PHASE FUNCTIONS [SR-1].
!       PATMS    AVERAGED PRESSURE FOR BAND MODEL SPECIES
!                ALONG SOLAR PATH [ATM].
!       WPATHS   SCATTERING POINT TO SUN MOLECULAR COLUMN AMOUNTS
!                FOR BAND MODEL SPECIES [ATM CM]
      INTEGER JTURN,LJ,IUPPHS,NANGLS
      REAL ATHETA,ADBETA,PHSFAC,PHSANG,PHSCOS,ANGF,PHASFN,AH1,ARH,      &
     &  ANGSUN,TBBYS,PATMS,WPATHS
      COMMON/SOLS/JTURN,LJ(LAYTWO+1),ATHETA(LAYDIM+1),ADBETA(LAYDIM+1), &
     &  PHSANG(LAYTWO),PHSCOS(LAYTWO),IUPPHS(LAYTWO),PHSFAC(LAYTWO),    &
     &  NANGLS,ANGF(MANGLS),PHASFN(LAYTWO,4),AH1(LAYTWO),ARH(LAYTWO),   &
     &  ANGSUN,TBBYS(LAYTHR,12),PATMS(LAYTHR,12),WPATHS(LAYTHR,MEXTX)
