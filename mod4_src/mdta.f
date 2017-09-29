      BLOCK DATA MDTA

!     CLOUD AND RAIN DATA
      INCLUDE 'PARAMS.h'

!     /CLDRR/
!       PCLD     CLOUD PROFILE PRESSURE [MB].
      REAL PCLD,ZCLD,CLD,CLDICE,RR
      COMMON/CLDRR/PCLD,ZCLD(1:NZCLD,0:1),CLD(1:NZCLD,0:5),             &
     &  CLDICE(1:NZCLD,0:1),RR(1:NZCLD,0:5)
      SAVE /CLDRR/

!     PROFILE ALTITUDES [KM]
      DATA ZCLD/NZCLD*0.,                                               &
     &    .00,  .16,  .33,  .66, 1.00, 1.50, 2.00,2.4,2.7,3.00,3.5,     &
     &   4.00, 4.50, 5.00, 5.50, 6.00, 6.5, 7.0, 7.5,8.0,8.5,9.0,9.5,   &
     &   10.0,10.5,11.0/

!     WATER DROPLET DENSITIES [GM/M3]
      DATA CLD/NZCLD*0.,                                                &
     &    .00,  .00,  .00,  .20,  .35, 1.00, 1.00,1.0, .3, .15,.0,15*.0,&
     &    .00,  .00,  .00,  .00,  .00,  .00,  .00, .3, .4, .30,.0,15*.0,&
     &    .00,  .00,  .15,  .30,  .15,  .00,  .00, .0, .0, .00,.0,15*.0,&
     &    .00,  .00,  .00,  .10,  .15,  .15,  .10, .0, .0, .00,.0,15*.0,&
     &    .00,  .30,  .65,  .40,  .00,  .00,  .00, .0, .0, .00,.0,15*.0/

!     ICE PARTICLE DENSITIES [GM/M3]
      DATA CLDICE/NZCLD*0.,NZCLD*0./

!     RAIN RATES [MM/HR]
      DATA RR/NZCLD*0.,                                                 &
     &   2.00, 1.78, 1.43, 1.22,  .86,  .22,  .00, .0, .0, .00,.0,15*.0,&
     &   5.00, 4.00, 3.40, 2.60,  .80,  .20,  .00, .0, .0, .00,.0,15*.0,&
     &  12.50,10.50, 8.00, 6.00, 2.50,  .80,  .20, .0, .0, .00,.0,15*.0,&
     &  25.00,21.50,17.50,12.00, 7.50, 4.20, 2.50,1.0, .7, .20,.0,15*.0,&
     &  75.00,70.00,65.00,60.00,45.00,20.00,12.50,7.0,3.5,1.00,.2,15*.0/
      END
