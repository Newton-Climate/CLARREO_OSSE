      BLOCK DATA SF296

!               06/28/82
!               UNITS OF (CM**3/MOL) * 1.E-20

!     COMMONS:

!     /S296DT/
!       V1S296   INITIAL SPECTRAL FREQUENCY [CM-1].
!       DVS296   SPECTRAL FREQUENCY STEP SIZE [CM-1].
!       N_S296   NUMBER OF SPECTRAL POINTS.
!       S296     SELF BROADENED H2O CONTINUUM DATA [1E-20*CM**3/MOL].
      INTEGER N_S296
      REAL V1S296,DVS296,S296
      COMMON/S296DT/V1S296,DVS296,N_S296,S296(2003)
      SAVE /S296DT/

!     LOCAL VARIABLES:
      INTEGER I

!     DATA:
      DATA V1S296,DVS296,N_S296/-20.,10.,2003/
      DATA (S296(I),I=1,2)/ 1.1109E-01 ,1.0573E-01/
      DATA (S296(I),I=3,52)/                                            &
     & 1.0162E-01, 1.0573E-01, 1.1109E-01, 1.2574E-01, 1.3499E-01,      &
     & 1.4327E-01, 1.5065E-01, 1.5164E-01, 1.5022E-01, 1.3677E-01,      &
     & 1.3115E-01, 1.2253E-01, 1.1271E-01, 1.0070E-01, 8.7495E-02,      &
     & 8.0118E-02, 6.9940E-02, 6.2034E-02, 5.6051E-02, 4.7663E-02,      &
     & 4.2450E-02, 3.6690E-02, 3.3441E-02, 3.0711E-02, 2.5205E-02,      &
     & 2.2113E-02, 1.8880E-02, 1.6653E-02, 1.4626E-02, 1.2065E-02,      &
     & 1.0709E-02, 9.1783E-03, 7.7274E-03, 6.7302E-03, 5.6164E-03,      &
     & 4.9089E-03, 4.1497E-03, 3.5823E-03, 3.1124E-03, 2.6414E-03,      &
     & 2.3167E-03, 2.0156E-03, 1.7829E-03, 1.5666E-03, 1.3928E-03,      &
     & 1.2338E-03, 1.0932E-03, 9.7939E-04, 8.8241E-04, 7.9173E-04/
      DATA (S296(I),I=53,102)/                                          &
     & 7.1296E-04, 6.4179E-04, 5.8031E-04, 5.2647E-04, 4.7762E-04,      &
     & 4.3349E-04, 3.9355E-04, 3.5887E-04, 3.2723E-04, 2.9919E-04,      &
     & 2.7363E-04, 2.5013E-04, 2.2876E-04, 2.0924E-04, 1.9193E-04,      &
     & 1.7618E-04, 1.6188E-04, 1.4891E-04, 1.3717E-04, 1.2647E-04,      &
     & 1.1671E-04, 1.0786E-04, 9.9785E-05, 9.2350E-05, 8.5539E-05,      &
     & 7.9377E-05, 7.3781E-05, 6.8677E-05, 6.3993E-05, 5.9705E-05,      &
     & 5.5788E-05, 5.2196E-05, 4.8899E-05, 4.5865E-05, 4.3079E-05,      &
     & 4.0526E-05, 3.8182E-05, 3.6025E-05, 3.4038E-05, 3.2203E-05,      &
     & 3.0511E-05, 2.8949E-05, 2.7505E-05, 2.6170E-05, 2.4933E-05,      &
     & 2.3786E-05, 2.2722E-05, 2.1736E-05, 2.0819E-05, 1.9968E-05/
      DATA (S296(I),I=103,152)/                                         &
     & 1.9178E-05, 1.8442E-05, 1.7760E-05, 1.7127E-05, 1.6541E-05,      &
     & 1.5997E-05, 1.5495E-05, 1.5034E-05, 1.4614E-05, 1.4230E-05,      &
     & 1.3883E-05, 1.3578E-05, 1.3304E-05, 1.3069E-05, 1.2876E-05,      &
     & 1.2732E-05, 1.2626E-05, 1.2556E-05, 1.2544E-05, 1.2604E-05,      &
     & 1.2719E-05, 1.2883E-05, 1.3164E-05, 1.3581E-05, 1.4187E-05,      &
     & 1.4866E-05, 1.5669E-05, 1.6717E-05, 1.8148E-05, 2.0268E-05,      &
     & 2.2456E-05, 2.5582E-05, 2.9183E-05, 3.3612E-05, 3.9996E-05,      &
     & 4.6829E-05, 5.5055E-05, 6.5897E-05, 7.5360E-05, 8.7213E-05,      &
     & 1.0046E-04, 1.1496E-04, 1.2943E-04, 1.5049E-04, 1.6973E-04,      &
     & 1.8711E-04, 2.0286E-04, 2.2823E-04, 2.6780E-04, 2.8766E-04/
      DATA (S296(I),I=153,202)/                                         &
     & 3.1164E-04, 3.3640E-04, 3.6884E-04, 3.9159E-04, 3.8712E-04,      &
     & 3.7433E-04, 3.4503E-04, 3.1003E-04, 2.8027E-04, 2.5253E-04,      &
     & 2.3408E-04, 2.2836E-04, 2.4442E-04, 2.7521E-04, 2.9048E-04,      &
     & 3.0489E-04, 3.2646E-04, 3.3880E-04, 3.3492E-04, 3.0987E-04,      &
     & 2.9482E-04, 2.8711E-04, 2.6068E-04, 2.2683E-04, 1.9996E-04,      &
     & 1.7788E-04, 1.6101E-04, 1.3911E-04, 1.2013E-04, 1.0544E-04,      &
     & 9.4224E-05, 8.1256E-05, 7.3667E-05, 6.2233E-05, 5.5906E-05,      &
     & 5.1619E-05, 4.5140E-05, 4.0273E-05, 3.3268E-05, 3.0258E-05,      &
     & 2.6440E-05, 2.3103E-05, 2.0749E-05, 1.8258E-05, 1.6459E-05,      &
     & 1.4097E-05, 1.2052E-05, 1.0759E-05, 9.1400E-06, 8.1432E-06/
      DATA (S296(I),I=203,252)/                                         &
     & 7.1460E-06, 6.4006E-06, 5.6995E-06, 4.9372E-06, 4.4455E-06,      &
     & 3.9033E-06, 3.4740E-06, 3.1269E-06, 2.8059E-06, 2.5558E-06,      &
     & 2.2919E-06, 2.0846E-06, 1.8983E-06, 1.7329E-06, 1.5929E-06,      &
     & 1.4631E-06, 1.3513E-06, 1.2461E-06, 1.1519E-06, 1.0682E-06,      &
     & 9.9256E-07, 9.2505E-07, 8.6367E-07, 8.0857E-07, 7.5674E-07,      &
     & 7.0934E-07, 6.6580E-07, 6.2580E-07, 5.8853E-07, 5.5333E-07,      &
     & 5.2143E-07, 4.9169E-07, 4.6431E-07, 4.3898E-07, 4.1564E-07,      &
     & 3.9405E-07, 3.7403E-07, 3.5544E-07, 3.3819E-07, 3.2212E-07,      &
     & 3.0714E-07, 2.9313E-07, 2.8003E-07, 2.6777E-07, 2.5628E-07,      &
     & 2.4551E-07, 2.3540E-07, 2.2591E-07, 2.1701E-07, 2.0866E-07/
      DATA (S296(I),I=253,302)/                                         &
     & 2.0082E-07, 1.9349E-07, 1.8665E-07, 1.8027E-07, 1.7439E-07,      &
     & 1.6894E-07, 1.6400E-07, 1.5953E-07, 1.5557E-07, 1.5195E-07,      &
     & 1.4888E-07, 1.4603E-07, 1.4337E-07, 1.4093E-07, 1.3828E-07,      &
     & 1.3569E-07, 1.3270E-07, 1.2984E-07, 1.2714E-07, 1.2541E-07,      &
     & 1.2399E-07, 1.2102E-07, 1.1878E-07, 1.1728E-07, 1.1644E-07,      &
     & 1.1491E-07, 1.1305E-07, 1.1235E-07, 1.1228E-07, 1.1224E-07,      &
     & 1.1191E-07, 1.1151E-07, 1.1098E-07, 1.1068E-07, 1.1109E-07,      &
     & 1.1213E-07, 1.1431E-07, 1.1826E-07, 1.2322E-07, 1.3025E-07,      &
     & 1.4066E-07, 1.5657E-07, 1.7214E-07, 1.9449E-07, 2.2662E-07,      &
     & 2.6953E-07, 3.1723E-07, 3.7028E-07, 4.4482E-07, 5.3852E-07/
      DATA (S296(I),I=303,352)/                                         &
     & 6.2639E-07, 7.2175E-07, 7.7626E-07, 8.7248E-07, 9.6759E-07,      &
     & 1.0102E-06, 1.0620E-06, 1.1201E-06, 1.2107E-06, 1.2998E-06,      &
     & 1.3130E-06, 1.2856E-06, 1.2350E-06, 1.1489E-06, 1.0819E-06,      &
     & 1.0120E-06, 9.4795E-07, 9.2858E-07, 9.8060E-07, 1.0999E-06,      &
     & 1.1967E-06, 1.2672E-06, 1.3418E-06, 1.3864E-06, 1.4330E-06,      &
     & 1.4592E-06, 1.4598E-06, 1.4774E-06, 1.4726E-06, 1.4820E-06,      &
     & 1.5050E-06, 1.4984E-06, 1.5181E-06, 1.5888E-06, 1.6850E-06,      &
     & 1.7690E-06, 1.9277E-06, 2.1107E-06, 2.3068E-06, 2.5347E-06,      &
     & 2.8069E-06, 3.1345E-06, 3.5822E-06, 3.9051E-06, 4.3422E-06,      &
     & 4.8704E-06, 5.5351E-06, 6.3454E-06, 7.2690E-06, 8.2974E-06/
      DATA (S296(I),I=353,402)/                                         &
     & 9.7609E-06, 1.1237E-05, 1.3187E-05, 1.5548E-05, 1.8784E-05,      &
     & 2.1694E-05, 2.5487E-05, 3.0092E-05, 3.5385E-05, 4.2764E-05,      &
     & 4.9313E-05, 5.5800E-05, 6.2968E-05, 7.1060E-05, 7.7699E-05,      &
     & 8.7216E-05, 8.9335E-05, 9.2151E-05, 9.2779E-05, 9.4643E-05,      &
     & 9.7978E-05, 1.0008E-04, 1.0702E-04, 1.1026E-04, 1.0828E-04,      &
     & 1.0550E-04, 1.0432E-04, 1.0428E-04, 9.8980E-05, 9.4992E-05,      &
     & 9.5159E-05, 1.0058E-04, 1.0738E-04, 1.1550E-04, 1.1229E-04,      &
     & 1.0596E-04, 1.0062E-04, 9.1742E-05, 8.4492E-05, 6.8099E-05,      &
     & 5.6295E-05, 4.6502E-05, 3.8071E-05, 3.0721E-05, 2.3297E-05,      &
     & 1.8688E-05, 1.4830E-05, 1.2049E-05, 9.6754E-06, 7.9192E-06/
      DATA (S296(I),I=403,452)/                                         &
     & 6.6673E-06, 5.6468E-06, 4.8904E-06, 4.2289E-06, 3.6880E-06,      &
     & 3.2396E-06, 2.8525E-06, 2.5363E-06, 2.2431E-06, 1.9949E-06,      &
     & 1.7931E-06, 1.6164E-06, 1.4431E-06, 1.2997E-06, 1.1559E-06,      &
     & 1.0404E-06, 9.4300E-07, 8.4597E-07, 7.6133E-07, 6.8623E-07,      &
     & 6.2137E-07, 5.6345E-07, 5.1076E-07, 4.6246E-07, 4.1906E-07,      &
     & 3.8063E-07, 3.4610E-07, 3.1554E-07, 2.8795E-07, 2.6252E-07,      &
     & 2.3967E-07, 2.1901E-07, 2.0052E-07, 1.8384E-07, 1.6847E-07,      &
     & 1.5459E-07, 1.4204E-07, 1.3068E-07, 1.2036E-07, 1.1095E-07,      &
     & 1.0237E-07, 9.4592E-08, 8.7530E-08, 8.1121E-08, 7.5282E-08,      &
     & 6.9985E-08, 6.5189E-08, 6.0874E-08, 5.6989E-08, 5.3530E-08/
      DATA (S296(I),I=453,502)/                                         &
     & 5.0418E-08, 4.7745E-08, 4.5367E-08, 4.3253E-08, 4.1309E-08,      &
     & 3.9695E-08, 3.8094E-08, 3.6482E-08, 3.4897E-08, 3.3500E-08,      &
     & 3.2302E-08, 3.0854E-08, 2.9698E-08, 2.8567E-08, 2.7600E-08,      &
     & 2.6746E-08, 2.5982E-08, 2.5510E-08, 2.5121E-08, 2.4922E-08,      &
     & 2.4909E-08, 2.5013E-08, 2.5216E-08, 2.5589E-08, 2.6049E-08,      &
     & 2.6451E-08, 2.6978E-08, 2.7687E-08, 2.8600E-08, 2.9643E-08,      &
     & 3.0701E-08, 3.2058E-08, 3.3695E-08, 3.5558E-08, 3.7634E-08,      &
     & 3.9875E-08, 4.2458E-08, 4.5480E-08, 4.8858E-08, 5.2599E-08,      &
     & 5.7030E-08, 6.2067E-08, 6.7911E-08, 7.4579E-08, 8.1902E-08,      &
     & 8.9978E-08, 9.9870E-08, 1.1102E-07, 1.2343E-07, 1.3732E-07/
      DATA (S296(I),I=503,552)/                                         &
     & 1.5394E-07, 1.7318E-07, 1.9383E-07, 2.1819E-07, 2.4666E-07,      &
     & 2.8109E-07, 3.2236E-07, 3.7760E-07, 4.4417E-07, 5.2422E-07,      &
     & 6.1941E-07, 7.4897E-07, 9.2041E-07, 1.1574E-06, 1.4126E-06,      &
     & 1.7197E-06, 2.1399E-06, 2.6266E-06, 3.3424E-06, 3.8418E-06,      &
     & 4.5140E-06, 5.0653E-06, 5.8485E-06, 6.5856E-06, 6.8937E-06,      &
     & 6.9121E-06, 6.9005E-06, 6.9861E-06, 6.8200E-06, 6.6089E-06,      &
     & 6.5809E-06, 7.3496E-06, 8.0311E-06, 8.3186E-06, 8.4260E-06,      &
     & 9.0644E-06, 9.4965E-06, 9.4909E-06, 9.0160E-06, 9.1494E-06,      &
     & 9.3629E-06, 9.5944E-06, 9.5459E-06, 8.9919E-06, 8.6040E-06,      &
     & 7.8613E-06, 7.1567E-06, 6.2677E-06, 5.1899E-06, 4.4188E-06/
      DATA (S296(I),I=553,602)/                                         &
     & 3.7167E-06, 3.0636E-06, 2.5573E-06, 2.0317E-06, 1.6371E-06,      &
     & 1.3257E-06, 1.0928E-06, 8.9986E-07, 7.4653E-07, 6.1111E-07,      &
     & 5.1395E-07, 4.3500E-07, 3.7584E-07, 3.2633E-07, 2.8413E-07,      &
     & 2.4723E-07, 2.1709E-07, 1.9294E-07, 1.7258E-07, 1.5492E-07,      &
     & 1.3820E-07, 1.2389E-07, 1.1189E-07, 1.0046E-07, 9.0832E-08,      &
     & 8.2764E-08, 7.4191E-08, 6.7085E-08, 6.0708E-08, 5.4963E-08,      &
     & 4.9851E-08, 4.5044E-08, 4.0916E-08, 3.7220E-08, 3.3678E-08,      &
     & 3.0663E-08, 2.7979E-08, 2.5495E-08, 2.3286E-08, 2.1233E-08,      &
     & 1.9409E-08, 1.7770E-08, 1.6260E-08, 1.4885E-08, 1.3674E-08,      &
     & 1.2543E-08, 1.1551E-08, 1.0655E-08, 9.8585E-09, 9.1398E-09/
      DATA (S296(I),I=603,652)/                                         &
     & 8.4806E-09, 7.8899E-09, 7.3547E-09, 6.8670E-09, 6.4131E-09,      &
     & 5.9930E-09, 5.6096E-09, 5.2592E-09, 4.9352E-09, 4.6354E-09,      &
     & 4.3722E-09, 4.1250E-09, 3.9081E-09, 3.7118E-09, 3.5372E-09,      &
     & 3.3862E-09, 3.2499E-09, 3.1324E-09, 3.0313E-09, 2.9438E-09,      &
     & 2.8686E-09, 2.8050E-09, 2.7545E-09, 2.7149E-09, 2.6907E-09,      &
     & 2.6724E-09, 2.6649E-09, 2.6642E-09, 2.6725E-09, 2.6871E-09,      &
     & 2.7056E-09, 2.7357E-09, 2.7781E-09, 2.8358E-09, 2.9067E-09,      &
     & 2.9952E-09, 3.1020E-09, 3.2253E-09, 3.3647E-09, 3.5232E-09,      &
     & 3.7037E-09, 3.9076E-09, 4.1385E-09, 4.3927E-09, 4.6861E-09,      &
     & 5.0238E-09, 5.4027E-09, 5.8303E-09, 6.3208E-09, 6.8878E-09/
      DATA (S296(I),I=653,702)/                                         &
     & 7.5419E-09, 8.3130E-09, 9.1952E-09, 1.0228E-08, 1.1386E-08,      &
     & 1.2792E-08, 1.4521E-08, 1.6437E-08, 1.8674E-08, 2.1160E-08,      &
     & 2.4506E-08, 2.8113E-08, 3.2636E-08, 3.7355E-08, 4.2234E-08,      &
     & 4.9282E-08, 5.7358E-08, 6.6743E-08, 7.8821E-08, 9.4264E-08,      &
     & 1.1542E-07, 1.3684E-07, 1.6337E-07, 2.0056E-07, 2.3252E-07,      &
     & 2.6127E-07, 2.9211E-07, 3.3804E-07, 3.7397E-07, 3.8205E-07,      &
     & 3.8810E-07, 3.9499E-07, 3.9508E-07, 3.7652E-07, 3.5859E-07,      &
     & 3.6198E-07, 3.7871E-07, 4.0925E-07, 4.2717E-07, 4.8241E-07,      &
     & 5.2008E-07, 5.6530E-07, 5.9531E-07, 6.1994E-07, 6.5080E-07,      &
     & 6.6355E-07, 6.9193E-07, 6.9930E-07, 7.3058E-07, 7.4678E-07/
      DATA (S296(I),I=703,752)/                                         &
     & 7.9193E-07, 8.3627E-07, 9.1267E-07, 1.0021E-06, 1.1218E-06,      &
     & 1.2899E-06, 1.4447E-06, 1.7268E-06, 2.0025E-06, 2.3139E-06,      &
     & 2.5599E-06, 2.8920E-06, 3.3059E-06, 3.5425E-06, 3.9522E-06,      &
     & 4.0551E-06, 4.2818E-06, 4.2892E-06, 4.4210E-06, 4.5614E-06,      &
     & 4.6739E-06, 4.9482E-06, 5.1118E-06, 5.0986E-06, 4.9417E-06,      &
     & 4.9022E-06, 4.8449E-06, 4.8694E-06, 4.8111E-06, 4.9378E-06,      &
     & 5.3231E-06, 5.7362E-06, 6.2350E-06, 6.0951E-06, 5.7281E-06,      &
     & 5.4585E-06, 4.9032E-06, 4.3009E-06, 3.4776E-06, 2.8108E-06,      &
     & 2.2993E-06, 1.7999E-06, 1.3870E-06, 1.0750E-06, 8.5191E-07,      &
     & 6.7951E-07, 5.5336E-07, 4.6439E-07, 4.0243E-07, 3.5368E-07/
      DATA (S296(I),I=753,802)/                                         &
     & 3.1427E-07, 2.7775E-07, 2.4486E-07, 2.1788E-07, 1.9249E-07,      &
     & 1.7162E-07, 1.5115E-07, 1.3478E-07, 1.2236E-07, 1.1139E-07,      &
     & 1.0092E-07, 9.0795E-08, 8.2214E-08, 7.4691E-08, 6.7486E-08,      &
     & 6.0414E-08, 5.4584E-08, 4.8754E-08, 4.3501E-08, 3.8767E-08,      &
     & 3.4363E-08, 3.0703E-08, 2.7562E-08, 2.4831E-08, 2.2241E-08,      &
     & 1.9939E-08, 1.8049E-08, 1.6368E-08, 1.4863E-08, 1.3460E-08,      &
     & 1.2212E-08, 1.1155E-08, 1.0185E-08, 9.3417E-09, 8.5671E-09,      &
     & 7.8292E-09, 7.1749E-09, 6.5856E-09, 6.0588E-09, 5.5835E-09,      &
     & 5.1350E-09, 4.7395E-09, 4.3771E-09, 4.0476E-09, 3.7560E-09,      &
     & 3.4861E-09, 3.2427E-09, 3.0240E-09, 2.8278E-09, 2.6531E-09/
      DATA (S296(I),I=803,852)/                                         &
     & 2.4937E-09, 2.3511E-09, 2.2245E-09, 2.1133E-09, 2.0159E-09,      &
     & 1.9330E-09, 1.8669E-09, 1.8152E-09, 1.7852E-09, 1.7752E-09,      &
     & 1.7823E-09, 1.8194E-09, 1.8866E-09, 1.9759E-09, 2.0736E-09,      &
     & 2.2083E-09, 2.3587E-09, 2.4984E-09, 2.6333E-09, 2.8160E-09,      &
     & 3.0759E-09, 3.3720E-09, 3.6457E-09, 4.0668E-09, 4.4541E-09,      &
     & 4.7976E-09, 5.0908E-09, 5.4811E-09, 6.1394E-09, 6.3669E-09,      &
     & 6.5714E-09, 6.8384E-09, 7.1918E-09, 7.3741E-09, 7.2079E-09,      &
     & 7.2172E-09, 7.2572E-09, 7.3912E-09, 7.6188E-09, 8.3291E-09,      &
     & 8.7885E-09, 9.2412E-09, 1.0021E-08, 1.0752E-08, 1.1546E-08,      &
     & 1.1607E-08, 1.1949E-08, 1.2346E-08, 1.2516E-08, 1.2826E-08/
      DATA (S296(I),I=853,902)/                                         &
     & 1.3053E-08, 1.3556E-08, 1.4221E-08, 1.5201E-08, 1.6661E-08,      &
     & 1.8385E-08, 2.0585E-08, 2.3674E-08, 2.7928E-08, 3.3901E-08,      &
     & 4.1017E-08, 4.9595E-08, 6.0432E-08, 7.6304E-08, 9.0764E-08,      &
     & 1.0798E-07, 1.2442E-07, 1.4404E-07, 1.6331E-07, 1.8339E-07,      &
     & 2.0445E-07, 2.2288E-07, 2.3083E-07, 2.3196E-07, 2.3919E-07,      &
     & 2.3339E-07, 2.3502E-07, 2.3444E-07, 2.6395E-07, 2.9928E-07,      &
     & 3.0025E-07, 3.0496E-07, 3.1777E-07, 3.4198E-07, 3.4739E-07,      &
     & 3.2696E-07, 3.4100E-07, 3.5405E-07, 3.7774E-07, 3.8285E-07,      &
     & 3.6797E-07, 3.5800E-07, 3.2283E-07, 2.9361E-07, 2.4881E-07,      &
     & 2.0599E-07, 1.7121E-07, 1.3641E-07, 1.1111E-07, 8.9413E-08/
      DATA (S296(I),I=903,952)/                                         &
     & 7.3455E-08, 6.2078E-08, 5.2538E-08, 4.5325E-08, 3.9005E-08,      &
     & 3.4772E-08, 3.1203E-08, 2.8132E-08, 2.5250E-08, 2.2371E-08,      &
     & 2.0131E-08, 1.7992E-08, 1.6076E-08, 1.4222E-08, 1.2490E-08,      &
     & 1.1401E-08, 1.0249E-08, 9.2279E-09, 8.5654E-09, 7.6227E-09,      &
     & 6.9648E-09, 6.2466E-09, 5.7252E-09, 5.3800E-09, 4.6960E-09,      &
     & 4.2194E-09, 3.7746E-09, 3.3813E-09, 3.0656E-09, 2.6885E-09,      &
     & 2.4311E-09, 2.1572E-09, 1.8892E-09, 1.7038E-09, 1.4914E-09,      &
     & 1.3277E-09, 1.1694E-09, 1.0391E-09, 9.2779E-10, 8.3123E-10,      &
     & 7.4968E-10, 6.8385E-10, 6.2915E-10, 5.7784E-10, 5.2838E-10,      &
     & 4.8382E-10, 4.4543E-10, 4.1155E-10, 3.7158E-10, 3.3731E-10/
      DATA (S296(I),I=953,1002)/                                        &
     & 3.0969E-10, 2.8535E-10, 2.6416E-10, 2.4583E-10, 2.2878E-10,      &
     & 2.1379E-10, 2.0073E-10, 1.8907E-10, 1.7866E-10, 1.6936E-10,      &
     & 1.6119E-10, 1.5424E-10, 1.4847E-10, 1.4401E-10, 1.4068E-10,      &
     & 1.3937E-10, 1.3943E-10, 1.4281E-10, 1.4766E-10, 1.5701E-10,      &
     & 1.7079E-10, 1.8691E-10, 2.0081E-10, 2.1740E-10, 2.4847E-10,      &
     & 2.6463E-10, 2.7087E-10, 2.7313E-10, 2.8352E-10, 2.9511E-10,      &
     & 2.8058E-10, 2.7227E-10, 2.7356E-10, 2.8012E-10, 2.8034E-10,      &
     & 2.9031E-10, 3.1030E-10, 3.3745E-10, 3.8152E-10, 4.0622E-10,      &
     & 4.2673E-10, 4.3879E-10, 4.5488E-10, 4.7179E-10, 4.6140E-10,      &
     & 4.6339E-10, 4.6716E-10, 4.7024E-10, 4.7931E-10, 4.8503E-10/
      DATA (S296(I),I=1003,1052)/                                       &
     & 4.9589E-10, 4.9499E-10, 5.0363E-10, 5.3184E-10, 5.6451E-10,      &
     & 6.0932E-10, 6.6469E-10, 7.4076E-10, 8.3605E-10, 9.4898E-10,      &
     & 1.0935E-09, 1.2593E-09, 1.4913E-09, 1.8099E-09, 2.1842E-09,      &
     & 2.7284E-09, 3.2159E-09, 3.7426E-09, 4.5226E-09, 5.3512E-09,      &
     & 6.1787E-09, 6.8237E-09, 7.9421E-09, 9.0002E-09, 9.6841E-09,      &
     & 9.9558E-09, 1.0232E-08, 1.0591E-08, 1.0657E-08, 1.0441E-08,      &
     & 1.0719E-08, 1.1526E-08, 1.2962E-08, 1.4336E-08, 1.6150E-08,      &
     & 1.8417E-08, 2.0725E-08, 2.3426E-08, 2.5619E-08, 2.7828E-08,      &
     & 3.0563E-08, 3.3438E-08, 3.6317E-08, 4.0400E-08, 4.4556E-08,      &
     & 5.0397E-08, 5.3315E-08, 5.9185E-08, 6.5311E-08, 6.9188E-08/
      DATA (S296(I),I=1053,1102)/                                       &
     & 7.7728E-08, 7.9789E-08, 8.6598E-08, 8.7768E-08, 9.1773E-08,      &
     & 9.7533E-08, 1.0007E-07, 1.0650E-07, 1.0992E-07, 1.0864E-07,      &
     & 1.0494E-07, 1.0303E-07, 1.0031E-07, 1.0436E-07, 1.0537E-07,      &
     & 1.1184E-07, 1.2364E-07, 1.3651E-07, 1.4881E-07, 1.4723E-07,      &
     & 1.4118E-07, 1.3371E-07, 1.1902E-07, 1.0007E-07, 7.9628E-08,      &
     & 6.4362E-08, 5.0243E-08, 3.8133E-08, 2.9400E-08, 2.3443E-08,      &
     & 1.9319E-08, 1.6196E-08, 1.4221E-08, 1.2817E-08, 1.1863E-08,      &
     & 1.1383E-08, 1.1221E-08, 1.1574E-08, 1.1661E-08, 1.2157E-08,      &
     & 1.2883E-08, 1.3295E-08, 1.4243E-08, 1.4240E-08, 1.4614E-08,      &
     & 1.4529E-08, 1.4685E-08, 1.4974E-08, 1.4790E-08, 1.4890E-08/
      DATA (S296(I),I=1103,1152)/                                       &
     & 1.4704E-08, 1.4142E-08, 1.3374E-08, 1.2746E-08, 1.2172E-08,      &
     & 1.2336E-08, 1.2546E-08, 1.3065E-08, 1.4090E-08, 1.5215E-08,      &
     & 1.6540E-08, 1.6144E-08, 1.5282E-08, 1.4358E-08, 1.2849E-08,      &
     & 1.0998E-08, 8.6956E-09, 7.0881E-09, 5.5767E-09, 4.2792E-09,      &
     & 3.2233E-09, 2.5020E-09, 1.9985E-09, 1.5834E-09, 1.3015E-09,      &
     & 1.0948E-09, 9.4141E-10, 8.1465E-10, 7.1517E-10, 6.2906E-10,      &
     & 5.5756E-10, 4.9805E-10, 4.3961E-10, 3.9181E-10, 3.5227E-10,      &
     & 3.1670E-10, 2.8667E-10, 2.5745E-10, 2.3212E-10, 2.0948E-10,      &
     & 1.8970E-10, 1.7239E-10, 1.5659E-10, 1.4301E-10, 1.3104E-10,      &
     & 1.2031E-10, 1.1095E-10, 1.0262E-10, 9.5130E-11, 8.8595E-11/
      DATA (S296(I),I=1153,1202)/                                       &
     & 8.2842E-11, 7.7727E-11, 7.3199E-11, 6.9286E-11, 6.5994E-11,      &
     & 6.3316E-11, 6.1244E-11, 5.9669E-11, 5.8843E-11, 5.8832E-11,      &
     & 5.9547E-11, 6.1635E-11, 6.4926E-11, 7.0745E-11, 7.8802E-11,      &
     & 8.6724E-11, 1.0052E-10, 1.1575E-10, 1.3626E-10, 1.5126E-10,      &
     & 1.6751E-10, 1.9239E-10, 2.1748E-10, 2.2654E-10, 2.2902E-10,      &
     & 2.3240E-10, 2.4081E-10, 2.3930E-10, 2.2378E-10, 2.2476E-10,      &
     & 2.2791E-10, 2.4047E-10, 2.5305E-10, 2.8073E-10, 3.1741E-10,      &
     & 3.6592E-10, 4.1495E-10, 4.6565E-10, 5.0990E-10, 5.5607E-10,      &
     & 6.1928E-10, 6.6779E-10, 7.3350E-10, 8.1434E-10, 8.9635E-10,      &
     & 9.9678E-10, 1.1256E-09, 1.2999E-09, 1.4888E-09, 1.7642E-09/
      DATA (S296(I),I=1203,1252)/                                       &
     & 1.9606E-09, 2.2066E-09, 2.4601E-09, 2.7218E-09, 3.0375E-09,      &
     & 3.1591E-09, 3.2852E-09, 3.2464E-09, 3.3046E-09, 3.2710E-09,      &
     & 3.2601E-09, 3.3398E-09, 3.7446E-09, 4.0795E-09, 4.0284E-09,      &
     & 4.0584E-09, 4.1677E-09, 4.5358E-09, 4.4097E-09, 4.2744E-09,      &
     & 4.5449E-09, 4.8147E-09, 5.2656E-09, 5.2476E-09, 5.0275E-09,      &
     & 4.7968E-09, 4.3654E-09, 3.9530E-09, 3.2447E-09, 2.6489E-09,      &
     & 2.1795E-09, 1.7880E-09, 1.4309E-09, 1.1256E-09, 9.1903E-10,      &
     & 7.6533E-10, 6.3989E-10, 5.5496E-10, 4.9581E-10, 4.5722E-10,      &
     & 4.3898E-10, 4.3505E-10, 4.3671E-10, 4.5329E-10, 4.6827E-10,      &
     & 4.9394E-10, 5.1122E-10, 5.1649E-10, 5.0965E-10, 4.9551E-10/
      DATA (S296(I),I=1253,1302)/                                       &
     & 4.8928E-10, 4.7947E-10, 4.7989E-10, 4.9071E-10, 4.8867E-10,      &
     & 4.7260E-10, 4.5756E-10, 4.5400E-10, 4.5993E-10, 4.4042E-10,      &
     & 4.3309E-10, 4.4182E-10, 4.6735E-10, 5.0378E-10, 5.2204E-10,      &
     & 5.0166E-10, 4.6799E-10, 4.3119E-10, 3.8803E-10, 3.3291E-10,      &
     & 2.6289E-10, 2.1029E-10, 1.7011E-10, 1.3345E-10, 1.0224E-10,      &
     & 7.8207E-11, 6.2451E-11, 5.0481E-11, 4.1507E-11, 3.5419E-11,      &
     & 3.0582E-11, 2.6900E-11, 2.3778E-11, 2.1343E-11, 1.9182E-11,      &
     & 1.7162E-11, 1.5391E-11, 1.3877E-11, 1.2619E-11, 1.1450E-11,      &
     & 1.0461E-11, 9.6578E-12, 8.9579E-12, 8.3463E-12, 7.8127E-12,      &
     & 7.3322E-12, 6.9414E-12, 6.6037E-12, 6.3285E-12, 6.1095E-12/
      DATA (S296(I),I=1303,1352)/                                       &
     & 5.9387E-12, 5.8118E-12, 5.7260E-12, 5.6794E-12, 5.6711E-12,      &
     & 5.7003E-12, 5.7670E-12, 5.8717E-12, 6.0151E-12, 6.1984E-12,      &
     & 6.4232E-12, 6.6918E-12, 7.0065E-12, 7.3705E-12, 7.7873E-12,      &
     & 8.2612E-12, 8.7972E-12, 9.4009E-12, 1.0079E-11, 1.0840E-11,      &
     & 1.1692E-11, 1.2648E-11, 1.3723E-11, 1.4935E-11, 1.6313E-11,      &
     & 1.7905E-11, 1.9740E-11, 2.1898E-11, 2.4419E-11, 2.7426E-11,      &
     & 3.0869E-11, 3.4235E-11, 3.7841E-11, 4.1929E-11, 4.6776E-11,      &
     & 5.2123E-11, 5.8497E-11, 6.5294E-11, 7.4038E-11, 8.4793E-11,      &
     & 9.6453E-11, 1.1223E-10, 1.2786E-10, 1.4882E-10, 1.7799E-10,      &
     & 2.0766E-10, 2.4523E-10, 2.8591E-10, 3.3386E-10, 4.0531E-10/
      DATA (S296(I),I=1353,1402)/                                       &
     & 4.7663E-10, 5.4858E-10, 6.3377E-10, 7.1688E-10, 8.4184E-10,      &
     & 9.5144E-10, 1.0481E-09, 1.1356E-09, 1.2339E-09, 1.3396E-09,      &
     & 1.4375E-09, 1.5831E-09, 1.7323E-09, 1.9671E-09, 2.2976E-09,      &
     & 2.6679E-09, 3.0777E-09, 3.4321E-09, 3.8192E-09, 4.2711E-09,      &
     & 4.4903E-09, 4.8931E-09, 5.2253E-09, 5.4040E-09, 5.6387E-09,      &
     & 5.6704E-09, 6.0345E-09, 6.1079E-09, 6.2576E-09, 6.4039E-09,      &
     & 6.3776E-09, 6.1878E-09, 5.8616E-09, 5.7036E-09, 5.5840E-09,      &
     & 5.6905E-09, 5.8931E-09, 6.2478E-09, 6.8291E-09, 7.4528E-09,      &
     & 7.6078E-09, 7.3898E-09, 6.7573E-09, 5.9827E-09, 5.0927E-09,      &
     & 4.0099E-09, 3.1933E-09, 2.4296E-09, 1.8485E-09, 1.4595E-09/
      DATA (S296(I),I=1403,1452)/                                       &
     & 1.2017E-09, 1.0164E-09, 8.7433E-10, 7.7108E-10, 7.0049E-10,      &
     & 6.5291E-10, 6.1477E-10, 5.9254E-10, 5.8150E-10, 5.7591E-10,      &
     & 5.8490E-10, 5.8587E-10, 5.9636E-10, 6.2408E-10, 6.5479E-10,      &
     & 7.0480E-10, 7.2313E-10, 7.5524E-10, 8.0863E-10, 8.3386E-10,      &
     & 9.2342E-10, 9.6754E-10, 1.0293E-09, 1.0895E-09, 1.1330E-09,      &
     & 1.2210E-09, 1.2413E-09, 1.2613E-09, 1.2671E-09, 1.2225E-09,      &
     & 1.1609E-09, 1.0991E-09, 1.0600E-09, 1.0570E-09, 1.0818E-09,      &
     & 1.1421E-09, 1.2270E-09, 1.3370E-09, 1.4742E-09, 1.4946E-09,      &
     & 1.4322E-09, 1.3210E-09, 1.1749E-09, 1.0051E-09, 7.8387E-10,      &
     & 6.1844E-10, 4.6288E-10, 3.4164E-10, 2.5412E-10, 1.9857E-10/
      DATA (S296(I),I=1453,1502)/                                       &
     & 1.5876E-10, 1.2966E-10, 1.0920E-10, 9.4811E-11, 8.3733E-11,      &
     & 7.3906E-11, 6.7259E-11, 6.1146E-11, 5.7119E-11, 5.3546E-11,      &
     & 4.8625E-11, 4.4749E-11, 4.1089E-11, 3.7825E-11, 3.4465E-11,      &
     & 3.1018E-11, 2.8109E-11, 2.5610E-11, 2.2859E-11, 2.0490E-11,      &
     & 1.8133E-11, 1.5835E-11, 1.3949E-11, 1.2295E-11, 1.0799E-11,      &
     & 9.6544E-12, 8.7597E-12, 7.9990E-12, 7.3973E-12, 6.9035E-12,      &
     & 6.4935E-12, 6.1195E-12, 5.8235E-12, 5.5928E-12, 5.4191E-12,      &
     & 5.2993E-12, 5.2338E-12, 5.2272E-12, 5.2923E-12, 5.4252E-12,      &
     & 5.6523E-12, 5.9433E-12, 6.3197E-12, 6.9016E-12, 7.5016E-12,      &
     & 8.2885E-12, 9.4050E-12, 1.0605E-11, 1.2257E-11, 1.3622E-11/
      DATA (S296(I),I=1503,1552)/                                       &
     & 1.5353E-11, 1.7543E-11, 1.9809E-11, 2.2197E-11, 2.4065E-11,      &
     & 2.6777E-11, 2.9751E-11, 3.2543E-11, 3.5536E-11, 3.9942E-11,      &
     & 4.6283E-11, 5.4556E-11, 6.5490E-11, 7.6803E-11, 9.0053E-11,      &
     & 1.0852E-10, 1.2946E-10, 1.4916E-10, 1.7748E-10, 2.0073E-10,      &
     & 2.2485E-10, 2.5114E-10, 2.7715E-10, 3.1319E-10, 3.3305E-10,      &
     & 3.5059E-10, 3.5746E-10, 3.6311E-10, 3.7344E-10, 3.6574E-10,      &
     & 3.7539E-10, 3.9434E-10, 4.3510E-10, 4.3340E-10, 4.2588E-10,      &
     & 4.3977E-10, 4.6062E-10, 4.7687E-10, 4.6457E-10, 4.8578E-10,      &
     & 5.2344E-10, 5.6752E-10, 5.8702E-10, 5.6603E-10, 5.3784E-10,      &
     & 4.9181E-10, 4.3272E-10, 3.5681E-10, 2.8814E-10, 2.3320E-10/
      DATA (S296(I),I=1553,1602)/                                       &
     & 1.8631E-10, 1.4587E-10, 1.1782E-10, 9.8132E-11, 8.2528E-11,      &
     & 6.9174E-11, 6.1056E-11, 5.3459E-11, 4.7116E-11, 4.1878E-11,      &
     & 3.8125E-11, 3.6347E-11, 3.5071E-11, 3.3897E-11, 3.3541E-11,      &
     & 3.3563E-11, 3.5469E-11, 3.8111E-11, 3.8675E-11, 4.1333E-11,      &
     & 4.3475E-11, 4.6476E-11, 4.9761E-11, 5.1380E-11, 5.4135E-11,      &
     & 5.3802E-11, 5.5158E-11, 5.6864E-11, 5.9311E-11, 6.3827E-11,      &
     & 6.7893E-11, 6.8230E-11, 6.6694E-11, 6.6018E-11, 6.4863E-11,      &
     & 6.5893E-11, 6.3813E-11, 6.4741E-11, 6.8630E-11, 7.0255E-11,      &
     & 7.0667E-11, 6.8810E-11, 6.4104E-11, 5.8136E-11, 4.7242E-11,      &
     & 3.7625E-11, 3.1742E-11, 2.5581E-11, 1.8824E-11, 1.3303E-11/
      DATA (S296(I),I=1603,1652)/                                       &
     & 9.6919E-12, 7.5353E-12, 6.0986E-12, 5.0742E-12, 4.3094E-12,      &
     & 3.7190E-12, 3.2520E-12, 2.8756E-12, 2.5680E-12, 2.3139E-12,      &
     & 2.1025E-12, 1.9257E-12, 1.7777E-12, 1.6539E-12, 1.5508E-12,      &
     & 1.4657E-12, 1.3966E-12, 1.3417E-12, 1.2998E-12, 1.2700E-12,      &
     & 1.2514E-12, 1.2437E-12, 1.2463E-12, 1.2592E-12, 1.2823E-12,      &
     & 1.3157E-12, 1.3596E-12, 1.4144E-12, 1.4806E-12, 1.5588E-12,      &
     & 1.6497E-12, 1.7544E-12, 1.8738E-12, 2.0094E-12, 2.1626E-12,      &
     & 2.3354E-12, 2.5297E-12, 2.7483E-12, 2.9941E-12, 3.2708E-12,      &
     & 3.5833E-12, 3.9374E-12, 4.3415E-12, 4.8079E-12, 5.3602E-12,      &
     & 5.9816E-12, 6.7436E-12, 7.6368E-12, 8.6812E-12, 9.8747E-12/
      DATA (S296(I),I=1653,1702)/                                       &
     & 1.1350E-11, 1.3181E-11, 1.5406E-11, 1.7868E-11, 2.0651E-11,      &
     & 2.4504E-11, 2.9184E-11, 3.4159E-11, 3.9979E-11, 4.8704E-11,      &
     & 5.7856E-11, 6.7576E-11, 7.9103E-11, 9.4370E-11, 1.1224E-10,      &
     & 1.3112E-10, 1.5674E-10, 1.8206E-10, 2.0576E-10, 2.3187E-10,      &
     & 2.7005E-10, 3.0055E-10, 3.3423E-10, 3.6956E-10, 3.8737E-10,      &
     & 4.2630E-10, 4.5154E-10, 4.8383E-10, 5.3582E-10, 5.8109E-10,      &
     & 6.3741E-10, 6.3874E-10, 6.3870E-10, 6.5818E-10, 6.5056E-10,      &
     & 6.5291E-10, 6.3159E-10, 6.3984E-10, 6.4549E-10, 6.5444E-10,      &
     & 6.7035E-10, 6.7665E-10, 6.9124E-10, 6.8451E-10, 6.9255E-10,      &
     & 6.9923E-10, 7.0396E-10, 6.7715E-10, 6.0371E-10, 5.3774E-10/
      DATA (S296(I),I=1703,1752)/                                       &
     & 4.6043E-10, 3.7635E-10, 2.9484E-10, 2.2968E-10, 1.8185E-10,      &
     & 1.4191E-10, 1.1471E-10, 9.4790E-11, 7.9613E-11, 6.7989E-11,      &
     & 5.9391E-11, 5.2810E-11, 4.7136E-11, 4.2618E-11, 3.8313E-11,      &
     & 3.4686E-11, 3.1669E-11, 2.9110E-11, 2.6871E-11, 2.5074E-11,      &
     & 2.4368E-11, 2.3925E-11, 2.4067E-11, 2.4336E-11, 2.4704E-11,      &
     & 2.5823E-11, 2.7177E-11, 2.9227E-11, 3.1593E-11, 3.5730E-11,      &
     & 4.0221E-11, 4.3994E-11, 4.8448E-11, 5.3191E-11, 5.8552E-11,      &
     & 6.3458E-11, 6.6335E-11, 7.2457E-11, 7.9091E-11, 8.2234E-11,      &
     & 8.7668E-11, 8.7951E-11, 9.2952E-11, 9.6157E-11, 9.5926E-11,      &
     & 1.0120E-10, 1.0115E-10, 9.9577E-11, 9.6633E-11, 9.2891E-11/
      DATA (S296(I),I=1753,1802)/                                       &
     & 9.3315E-11, 9.5584E-11, 1.0064E-10, 1.0509E-10, 1.1455E-10,      &
     & 1.2443E-10, 1.2963E-10, 1.2632E-10, 1.1308E-10, 1.0186E-10,      &
     & 8.5880E-11, 6.7863E-11, 5.1521E-11, 3.7780E-11, 2.8842E-11,      &
     & 2.2052E-11, 1.7402E-11, 1.4406E-11, 1.1934E-11, 1.0223E-11,      &
     & 8.9544E-12, 7.9088E-12, 7.0675E-12, 6.2222E-12, 5.6051E-12,      &
     & 5.0502E-12, 4.5578E-12, 4.2636E-12, 3.9461E-12, 3.7599E-12,      &
     & 3.5215E-12, 3.2467E-12, 3.0018E-12, 2.6558E-12, 2.3928E-12,      &
     & 2.0707E-12, 1.7575E-12, 1.5114E-12, 1.2941E-12, 1.1004E-12,      &
     & 9.5175E-13, 8.2894E-13, 7.3253E-13, 6.5551E-13, 5.9098E-13,      &
     & 5.3548E-13, 4.8697E-13, 4.4413E-13, 4.0600E-13, 3.7188E-13/
      DATA (S296(I),I=1803,1852)/                                       &
     & 3.4121E-13, 3.1356E-13, 2.8856E-13, 2.6590E-13, 2.4533E-13,      &
     & 2.2663E-13, 2.0960E-13, 1.9407E-13, 1.7990E-13, 1.6695E-13,      &
     & 1.5512E-13, 1.4429E-13, 1.3437E-13, 1.2527E-13, 1.1693E-13,      &
     & 1.0927E-13, 1.0224E-13, 9.5767E-14, 8.9816E-14, 8.4335E-14,      &
     & 7.9285E-14, 7.4626E-14, 7.0325E-14, 6.6352E-14, 6.2676E-14,      &
     & 5.9274E-14, 5.6121E-14, 5.3195E-14, 5.0479E-14, 4.7953E-14,      &
     & 4.5602E-14, 4.3411E-14, 4.1367E-14, 3.9456E-14, 3.7670E-14,      &
     & 3.5996E-14, 3.4427E-14, 3.2952E-14, 3.1566E-14, 3.0261E-14,      &
     & 2.9030E-14, 2.7868E-14, 2.6770E-14, 2.5730E-14, 2.4745E-14,      &
     & 2.3809E-14, 2.2921E-14, 2.2076E-14, 2.1271E-14, 2.0504E-14/
      DATA (S296(I),I=1853,1902)/                                       &
     & 1.9772E-14, 1.9073E-14, 1.8404E-14, 1.7764E-14, 1.7151E-14,      &
     & 1.6564E-14, 1.6000E-14, 1.5459E-14, 1.4939E-14, 1.4439E-14,      &
     & 1.3958E-14, 1.3495E-14, 1.3049E-14, 1.2620E-14, 1.2206E-14,      &
     & 1.1807E-14, 1.1422E-14, 1.1050E-14, 1.0691E-14, 1.0345E-14,      &
     & 1.0010E-14, 9.6870E-15, 9.3747E-15, 9.0727E-15, 8.7808E-15,      &
     & 8.4986E-15, 8.2257E-15, 7.9617E-15, 7.7064E-15, 7.4594E-15,      &
     & 7.2204E-15, 6.9891E-15, 6.7653E-15, 6.5488E-15, 6.3392E-15,      &
     & 6.1363E-15, 5.9399E-15, 5.7499E-15, 5.5659E-15, 5.3878E-15,      &
     & 5.2153E-15, 5.0484E-15, 4.8868E-15, 4.7303E-15, 4.5788E-15,      &
     & 4.4322E-15, 4.2902E-15, 4.1527E-15, 4.0196E-15, 3.8907E-15/
      DATA (S296(I),I=1903,1952)/                                       &
     & 3.7659E-15, 3.6451E-15, 3.5281E-15, 3.4149E-15, 3.3052E-15,      &
     & 3.1991E-15, 3.0963E-15, 2.9967E-15, 2.9004E-15, 2.8071E-15,      &
     & 2.7167E-15, 2.6293E-15, 2.5446E-15, 2.4626E-15, 2.3833E-15,      &
     & 2.3064E-15, 2.2320E-15, 2.1600E-15, 2.0903E-15, 2.0228E-15,      &
     & 1.9574E-15, 1.8942E-15, 1.8329E-15, 1.7736E-15, 1.7163E-15,      &
     & 1.6607E-15, 1.6069E-15, 1.5548E-15, 1.5044E-15, 1.4557E-15,      &
     & 1.4084E-15, 1.3627E-15, 1.3185E-15, 1.2757E-15, 1.2342E-15,      &
     & 1.1941E-15, 1.1552E-15, 1.1177E-15, 1.0813E-15, 1.0461E-15,      &
     & 1.0120E-15, 9.7900E-16, 9.4707E-16, 9.1618E-16, 8.8628E-16,      &
     & 8.5734E-16, 8.2933E-16, 8.0223E-16, 7.7600E-16, 7.5062E-16/
      DATA (S296(I),I=1953,2002)/                                       &
     & 7.2606E-16, 7.0229E-16, 6.7929E-16, 6.5703E-16, 6.3550E-16,      &
     & 6.1466E-16, 5.9449E-16, 5.7498E-16, 5.5610E-16, 5.3783E-16,      &
     & 5.2015E-16, 5.0305E-16, 4.8650E-16, 4.7049E-16, 4.5500E-16,      &
     & 4.4002E-16, 4.2552E-16, 4.1149E-16, 3.9792E-16, 3.8479E-16,      &
     & 3.7209E-16, 3.5981E-16, 3.4792E-16, 3.3642E-16, 3.2530E-16,      &
     & 3.1454E-16, 3.0413E-16, 2.9406E-16, 2.8432E-16, 2.7490E-16,      &
     & 2.6579E-16, 2.5697E-16, 2.4845E-16, 2.4020E-16, 2.3223E-16,      &
     & 2.2451E-16, 2.1705E-16, 2.0984E-16, 2.0286E-16, 1.9611E-16,      &
     & 1.8958E-16, 1.8327E-16, 1.7716E-16, 1.7126E-16, 1.6555E-16,      &
     & 1.6003E-16, 1.5469E-16, 1.4952E-16, 1.4453E-16, 1.3970E-16/
      DATA S296(2003)/1.3503E-16/
      END
