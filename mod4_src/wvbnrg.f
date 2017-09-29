      BLOCK DATA WVBNRG
!>    BLOCK DATA
!     WAVENUMBER-LOW AND WAVENUMBER-HIGH SPECIFY A BAND REGION
!     FOR A MOLECULAR ABSORBER.
!     THE UNIT FOR WAVENUMBER IS 1/CM.
!     -999 IS AN INDICATOR TO INDICATE THE END OF ABSORPTION BANDS
!     FOR ANY SPECIFIC ABSORBER.
      COMMON /WNLOHI/                                                   &
     &   IWLH2O(15),IWLO3 ( 6),IWLCO2(11),IWLCO ( 4),IWLCH4( 5),        &
     &   IWLN2O(12),IWLO2 ( 7),IWLNH3( 3),IWLNO ( 2),IWLNO2( 4),        &
     &   IWLSO2( 5),                                                    &
     &   IWHH2O(15),IWHO3 ( 6),IWHCO2(11),IWHCO ( 4),IWHCH4( 5),        &
     &   IWHN2O(12),IWHO2 ( 7),IWHNH3( 3),IWHNO ( 2),IWHNO2( 4),        &
     &   IWHSO2( 5)
      SAVE /WNLOHI/

      DATA IWLH2O/   0,    350,   1005,   1645,   2535,   3425,   4315, &
     &    6155,   8005,   9620,  11545,  13075,  14865,  16340,   -999/
      DATA IWHH2O/ 345,   1000,   1640,   2530,   3420,   4310,   6150, &
     &    8000,   9615,  11540,  13070,  14860,  16045,  17860,   -999/

      DATA IWLO3 /   0,    515,   1630,   2670,   2850,   -999/
      DATA IWHO3 / 200,   1275,   2295,   2845,   3260,   -999/

      DATA IWLCO2/ 425,    840,   1805,   3070,   3760,   4530,   5905, &
     &    7395,   8030,   9340,   -999/
      DATA IWHCO2/ 835,   1440,   2855,   3755,   4065,   5380,   7025, &
     &    7785,   8335,   9670,   -999/

      DATA IWLCO /   0,   1940,   4040,   -999/
      DATA IWHCO / 175,   2285,   4370,   -999/

      DATA IWLCH4/1065,   2345,   4110,   5865,   -999/
      DATA IWHCH4/1775,   3230,   4690,   6135,   -999/

      DATA IWLN2O/   0,    490,    865,   1065,   1545,   2090,   2705, &
     &    3245,   4260,   4540,   4910,   -999/
      DATA IWHN2O/ 120,    775,    995,   1385,   2040,   2655,   2865, &
     &    3925,   4470,   4785,   5165,   -999/

      DATA IWLO2 /   0,   7650,   9235,  12850,  14300,  15695,   -999/
      DATA IWHO2 / 265,   8080,   9490,  13220,  14600,  15955,   -999/

      DATA IWLNH3/   0,    390,   -999/
      DATA IWHNH3/ 385,   2150,   -999/

      DATA IWLNO /1700,   -999/
      DATA IWHNO /2005,   -999/

      DATA IWLNO2/ 580,   1515,   2800,   -999/
      DATA IWHNO2/ 925,   1695,   2970,   -999/

      DATA IWLSO2/   0,    400,    950,   2415,   -999/
      DATA IWHSO2/ 185,    650,   1460,   2580,   -999/
      END
