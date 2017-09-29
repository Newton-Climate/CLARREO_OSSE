      BLOCK DATA CPO3BD

!     CPRIME FOR 5 O3 BANDS:

!     COMMONS:
      REAL CPO3
      COMMON/O3/CPO3(447)
      SAVE /O3/

!     LOCAL VARIABLES:
!       I       ARRAY INDEX.
      INTEGER I

!     0-200 CM-1 OZONE BAND:
      DATA (CPO3(I),I=1,41)/                                            &
     & -2.0427, -1.8966, -1.6263, -1.3896, -1.2170, -1.0996, -1.0214,   &
     & -0.9673, -0.9249, -0.8896, -0.8612, -0.8417, -0.8360, -0.8483,   &
     & -0.8785, -0.9273, -0.9932, -1.0720, -1.1639, -1.2662, -1.3771,   &
     & -1.4976, -1.6274, -1.7712, -1.9289, -2.1027, -2.2948, -2.4987,   &
     & -2.7321, -2.9992, -3.3045, -3.6994, -4.1022, -4.6467, -5.1328,   &
     & -5.6481, -6.1634, -6.6787, -7.1940, -7.7093, -8.0000/

!     515-1275 CM-1 OZONE BAND:
      DATA (CPO3(I),I=42,167)/                                          &
     & -7.9274, -7.6418, -7.3562, -7.0706, -6.7850, -6.4994, -6.2138,   &
     & -5.9282, -5.6426, -5.3570, -5.0714, -4.7858, -4.5002, -4.2146,   &
     & -3.9290, -3.6213, -3.3407, -3.0722, -2.8226, -2.5914, -2.3778,   &
     & -2.1823, -2.0057, -1.8456, -1.6991, -1.5659, -1.4436, -1.3323,   &
     & -1.2319, -1.1407, -1.0550, -0.9733, -0.9033, -0.8584, -0.8527,   &
     & -0.8838, -0.9219, -0.9360, -0.9025, -0.8402, -0.7913, -0.7794,   &
     & -0.8123, -0.8750, -0.9484, -1.0206, -1.0864, -1.1520, -1.2202,   &
     & -1.2928, -1.3745, -1.4641, -1.5611, -1.6669, -1.7816, -1.9051,   &
     & -2.0383, -2.1796, -2.3312, -2.4906, -2.6569, -2.8354, -3.0179,   &
     & -3.2121, -3.4106, -3.6208, -3.8332, -4.0584, -4.2854, -4.4979,   &
     & -4.7175, -4.9109, -5.1246, -5.3344, -5.5442, -5.7540, -5.9638,   &
     & -6.1736, -6.3834, -6.5932, -6.8030, -7.0128, -6.9011, -6.2590,   &
     & -5.8119, -5.1603, -4.3327, -3.6849, -3.1253, -2.6304, -2.1903,   &
     & -1.8019, -1.4585, -1.1533, -0.8770, -0.6166, -0.3630, -0.1102,   &
     &  0.1336,  0.3525,  0.5326,  0.6678,  0.7510,  0.7752,  0.7826,   &
     &  0.7874,  0.8006,  0.8241,  0.7614,  0.5662,  0.1949, -0.2770,   &
     & -0.6199, -0.8347, -0.9586, -1.0168, -1.0501, -1.0816, -1.0980,   &
     & -1.0833, -1.0424, -0.9972, -0.9724, -0.9855, -1.0365, -1.1187/
      DATA (CPO3(I),I=168,194)/                                         &
     & -1.2150, -1.3142, -1.4103, -1.4998, -1.5933, -1.6938, -1.8061,   &
     & -1.9332, -2.0737, -2.2279, -2.3966, -2.5787, -2.7755, -2.9855,   &
     & -3.2090, -3.4465, -3.6967, -3.9633, -4.2461, -4.5502, -4.8912,   &
     & -5.2845, -5.7654, -6.4194, -6.9288, -7.4382, -7.9476/

!     1630-2295 CM-1 OZONE BAND:
      DATA (CPO3(I),I=195,320)/                                         &
     & -8.0000, -7.5432, -6.9273, -6.3115, -5.5431, -4.9563, -4.4640,   &
     & -4.0371, -3.6533, -3.3069, -2.9877, -2.7042, -2.4507, -2.2355,   &
     & -2.0651, -1.9477, -1.8705, -1.8422, -1.8235, -1.7782, -1.7367,   &
     & -1.7012, -1.7208, -1.8353, -2.0331, -2.3077, -2.5996, -2.7517,   &
     & -2.7263, -2.6671, -2.6415, -2.6449, -2.6613, -2.6589, -2.6083,   &
     & -2.5250, -2.4529, -2.4157, -2.4298, -2.4906, -2.5823, -2.6873,   &
     & -2.7808, -2.8612, -2.9303, -3.0022, -3.0873, -3.1844, -3.2929,   &
     & -3.4158, -3.5361, -3.6710, -3.8062, -3.9520, -4.1140, -4.2635,   &
     & -4.4395, -4.6138, -4.8372, -5.0837, -5.3302, -5.3665, -5.4358,   &
     & -5.0651, -4.8416, -4.5293, -4.2547, -4.0039, -3.7818, -3.5850,   &
     & -3.4091, -3.2509, -3.0934, -2.9485, -2.8055, -2.6705, -2.5482,   &
     & -2.4362, -2.3380, -2.2486, -2.1645, -2.0834, -2.0035, -1.9081,   &
     & -1.7681, -1.5768, -1.3615, -1.1463, -0.9482, -0.7800, -0.6336,   &
     & -0.5092, -0.4105, -0.3495, -0.3274, -0.3133, -0.3023, -0.2859,   &
     & -0.3055, -0.4374, -0.6972, -1.1064, -1.4904, -1.9687, -2.4498,   &
     & -2.5971, -2.5220, -2.4301, -2.3467, -2.2901, -2.2746, -2.3021,   &
     & -2.3635, -2.4420, -2.5088, -2.5485, -2.5617, -2.5656, -2.5771,   &
     & -2.6134, -2.6822, -2.7885, -2.9379, -3.1200, -3.3260, -3.5464/
      DATA (CPO3(I),I=321,328)/                                         &
     & -3.7736, -4.0311, -4.3651, -4.7794, -5.5152, -6.1240, -7.2193,   &
     & -8.0000/

!     2670-2845 CM-1 OZONE BAND:
      DATA (CPO3(I),I=329,364)/                                         &
     & -7.9721, -7.6118, -7.2515, -6.8913, -6.5310, -6.1707, -5.8105,   &
     & -5.4502, -5.0899, -4.7297, -4.3694, -3.9462, -3.6022, -3.2886,   &
     & -3.0234, -2.7863, -2.5797, -2.4073, -2.2760, -2.1894, -2.1359,   &
     & -2.1160, -2.0808, -2.0151, -1.9666, -1.9409, -1.9868, -2.1450,   &
     & -2.3965, -2.8042, -3.5500, -4.8275, -5.6378, -6.4482, -7.2585,   &
     & -8.0000/

!     2850-3260 CM-1 OZONE BAND:
      DATA (CPO3(I),I=365,447)/                                         &
     & -8.0000, -7.6278, -7.2556, -6.8834, -6.5111, -6.1389, -5.7667,   &
     & -5.3945, -5.0223, -4.6501, -4.2779, -3.9056, -3.5334, -3.3828,   &
     & -3.2452, -3.1411, -3.0403, -2.9428, -2.8436, -2.7573, -2.6853,   &
     & -2.6040, -2.5218, -2.4121, -2.3547, -2.1970, -2.0668, -1.9121,   &
     & -1.7617, -1.6153, -1.4688, -1.4022, -1.3447, -1.2669, -1.1902,   &
     & -1.1805, -1.1707, -1.1609, -1.1609, -1.1805, -1.1999, -1.4214,   &
     & -1.6348, -1.7519, -1.9730, -2.2078, -2.4608, -2.5337, -2.5923,   &
     & -2.6616, -2.6384, -2.6271, -2.6154, -2.5570, -2.4983, -2.4480,   &
     & -2.3890, -2.3663, -2.3431, -2.3314, -2.3200, -2.3200, -2.3314,   &
     & -2.3431, -2.3547, -2.3777, -2.4004, -2.5218, -2.6499, -2.7694,   &
     & -2.9057, -3.0286, -3.1543, -3.3696, -3.6053, -4.1977, -4.7811,   &
     & -5.2933, -5.7554, -6.4542, -7.0239, -7.5937, -8.0000/
      END