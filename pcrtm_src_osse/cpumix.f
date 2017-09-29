      BLOCK DATA CPUMIX

!     C' FOR UNIFORMLY MIXED GASES (CO2, CO, CH4, N2O, O2)

!     COMMONS:

!     /UFMIX/
      REAL CPCO2,CPCO,CPCH4,CPN2O,CPO2
      COMMON/UFMIX/CPCO2(1219),CPCO(173),CPCH4(493),CPN2O(704),CPO2(382)

!     LOCAL VARIABLES:
!       IND      C PRIME INDEX.
      INTEGER IND

!     SAVED COMMONS:
      SAVE /UFMIX/
!=CO2 ====C' FOR    8 BAND MODELS
!=CO2 ====  425-  835
      DATA (CPCO2(IND),IND=   1,  83)/                                  &
     & -9.8495, -9.6484, -9.4472, -9.2461, -9.0449, -8.9544, -8.6127,   &
     & -8.4076, -8.2710, -8.0391, -7.9485, -7.9638, -7.7849, -7.6278,   &
     & -7.1418, -6.7823, -6.3826, -6.0323, -5.7501, -5.5249, -5.3304,   &
     & -5.0105, -4.7703, -4.5714, -4.3919, -4.2974, -4.1370, -3.8761,   &
     & -3.5936, -3.2852, -3.0016, -2.7303, -2.4868, -2.2741, -2.0936,   &
     & -1.9424, -1.8092, -1.6843, -1.5372, -1.3803, -1.2043, -0.9930,   &
     & -0.7724, -0.5509, -0.3465, -0.1785, -0.0470,  0.0449,  0.1114,   &
     &  0.1367,  0.0910,  0.0066, -0.1269, -0.2994, -0.4934, -0.7101,   &
     & -0.9087, -1.1004, -1.2694, -1.4064, -1.5622, -1.6810, -1.7841,   &
     & -1.8973, -2.0274, -2.2079, -2.4264, -2.6763, -2.9312, -3.1896,   &
     & -3.4262, -3.5979, -3.7051, -3.7372, -3.7983, -3.9154, -4.0520,   &
     & -4.2567, -4.4661, -4.6670, -4.9226, -5.2203, -5.5597/
!=CO2 ====  840- 1440
      DATA (CPCO2(IND),IND=  84, 204)/                                  &
     & -5.6403, -5.7039, -5.7674, -5.8310, -5.8948, -5.9503, -6.0217,   &
     & -6.0392, -5.9855, -5.8620, -5.6834, -5.5083, -5.3473, -5.2028,   &
     & -5.0799, -4.9628, -4.8379, -4.7032, -4.5584, -4.4213, -4.3198,   &
     & -4.2786, -4.2843, -4.3099, -4.3210, -4.2769, -4.2229, -4.2179,   &
     & -4.2950, -4.4789, -4.7550, -5.0902, -5.4329, -5.6689, -5.6608,   &
     & -5.4582, -5.1969, -4.9419, -4.7106, -4.5084, -4.3409, -4.2211,   &
     & -4.1563, -4.1259, -4.1108, -4.0803, -4.0211, -3.9824, -4.0053,   &
     & -4.1221, -4.3504, -4.6741, -5.0826, -5.5857, -6.2301, -7.0829,   &
     & -8.1344, -8.8601, -9.0457, -9.1231, -9.0728, -9.1413, -9.1221,   &
     & -9.1882, -9.2752, -9.2237, -9.3604, -9.3058, -9.5455, -9.5567,   &
     & -9.3754, -8.7756, -8.0904, -7.4827, -6.9585, -6.5095, -6.1194,   &
     & -5.7824, -5.4910, -5.2532, -5.0840, -4.9920, -4.9577, -4.9638,   &
     & -4.9741, -4.9555, -4.9466, -4.9774, -5.0719, -5.2558, -5.5213,   &
     & -5.8633, -6.2877, -6.7878, -7.2602, -7.2940, -6.8524, -6.3372,   &
     & -5.8854, -5.5065, -5.2011, -4.9776, -4.8471, -4.7885, -4.7783,   &
     & -4.7815, -4.7538, -4.7228, -4.7259, -4.7860, -4.9231, -5.1270,   &
     & -5.3831, -5.6849, -6.0351, -6.4437, -6.9160, -7.4815, -8.1437,   &
     & -8.9449, -9.8564/
!=CO2 ==== 1805- 2855
      DATA (CPCO2(IND),IND= 205, 330)/                                  &
     & -9.8903, -9.4365, -8.9826, -8.5288, -8.1184, -7.6555, -7.1673,   &
     & -6.7226, -6.3423, -6.0410, -5.8154, -5.6519, -5.5186, -5.3859,   &
     & -5.2279, -5.0238, -4.7865, -4.5343, -4.2846, -4.0560, -3.8717,   &
     & -3.7624, -3.7231, -3.7335, -3.8312, -3.9854, -4.1930, -4.4895,   &
     & -4.7394, -4.8892, -4.9499, -4.9392, -4.9787, -5.1129, -5.3330,   &
     & -5.6093, -5.8862, -6.0581, -6.0274, -5.8356, -5.5989, -5.3738,   &
     & -5.1661, -4.9472, -4.7020, -4.4354, -4.1439, -3.8561, -3.5944,   &
     & -3.3694, -3.2100, -3.1041, -3.0411, -3.0471, -3.1077, -3.2305,   &
     & -3.4274, -3.6115, -3.7542, -3.8666, -3.9338, -4.0079, -4.0962,   &
     & -4.2142, -4.1433, -4.2870, -4.4796, -4.6618, -4.8204, -4.9499,   &
     & -4.9862, -5.0171, -5.0282, -5.0580, -5.0398, -4.9465, -4.7816,   &
     & -4.5538, -4.2975, -4.0286, -3.7528, -3.4715, -3.1899, -2.9041,   &
     & -2.6127, -2.3212, -2.0435, -1.7894, -1.5531, -1.3382, -1.1515,   &
     & -0.9990, -0.8833, -0.8006, -0.7227, -0.6288, -0.4977, -0.3249,   &
     & -0.1349,  0.0576,  0.2487,  0.4386,  0.6260,  0.8081,  0.9681,   &
     &  1.0859,  1.1522,  1.1861,  1.2039,  1.2255,  1.2587,  1.2473,   &
     &  1.1457,  0.9139,  0.5250,  0.0173, -0.5796, -1.3944, -2.3841,   &
     & -2.7244, -2.9264, -3.0689, -3.2120, -3.3353, -3.4510, -3.5566/
      DATA (CPCO2(IND),IND= 331, 415)/                                  &
     & -3.6518, -3.7460, -3.8500, -3.9680, -4.0981, -4.2259, -4.3369,   &
     & -4.4329, -4.5305, -4.6264, -4.7438, -4.8842, -5.0248, -5.1448,   &
     & -5.2371, -5.2781, -5.3299, -5.3766, -5.4233, -5.4699, -5.5166,   &
     & -5.5633, -5.6646, -5.7593, -5.8461, -5.9229, -5.9818, -6.0065,   &
     & -5.9747, -5.8741, -5.7230, -5.5620, -5.4389, -5.3788, -5.3679,   &
     & -5.3827, -5.3837, -5.3460, -5.3186, -5.3394, -5.4320, -5.6095,   &
     & -5.8446, -6.0992, -6.3399, -6.5499, -6.7434, -6.9359, -7.1219,   &
     & -7.2818, -7.3984, -7.4881, -7.5452, -7.5994, -7.6445, -7.6734,   &
     & -7.6422, -7.5057, -7.2650, -6.9975, -6.7749, -6.6398, -6.5875,   &
     & -6.5912, -6.6192, -6.6155, -6.5866, -6.5851, -6.6382, -6.7736,   &
     & -7.0009, -7.2896, -7.6327, -7.9767, -8.2633, -8.4744, -8.5455,   &
     & -8.5813, -8.6025, -8.6459, -8.8948, -9.1436, -9.3925, -9.6413,   &
     & -9.8902/
!=CO2 ==== 3070- 3755
      DATA (CPCO2(IND),IND= 416, 541)/                                  &
     & -9.8006, -9.5049, -9.1947, -8.7254, -8.4410, -8.1781, -8.0182,   &
     & -7.9381, -7.8793, -7.7636, -7.5549, -7.2962, -7.0244, -6.7556,   &
     & -6.4888, -6.2443, -6.0422, -5.9088, -5.8590, -5.8890, -5.9850,   &
     & -6.0949, -6.1164, -6.0207, -5.8592, -5.7110, -5.6328, -5.6369,   &
     & -5.7274, -5.9069, -6.1720, -6.5203, -6.9586, -7.4776, -8.0607,   &
     & -8.5514, -8.7011, -8.4232, -7.9274, -7.6159, -7.3836, -7.1969,   &
     & -7.0523, -6.7685, -6.4022, -6.0354, -5.7125, -5.4659, -5.3088,   &
     & -5.2546, -5.2991, -5.3819, -5.4615, -5.4117, -5.2107, -5.0103,   &
     & -4.8232, -4.7071, -4.6850, -4.7385, -4.8797, -5.1024, -5.4015,   &
     & -5.7758, -6.2225, -6.6681, -6.9127, -6.8919, -6.6972, -6.5012,   &
     & -6.3123, -6.1091, -5.8641, -5.5889, -5.3057, -5.0340, -4.7826,   &
     & -4.5476, -4.3277, -4.1224, -3.9333, -3.7675, -3.6324, -3.5163,   &
     & -3.4043, -3.2744, -3.1180, -2.9557, -2.8254, -2.7359, -2.6721,   &
     & -2.6084, -2.5105, -2.3772, -2.2317, -2.0866, -1.9521, -1.8292,   &
     & -1.7110, -1.5992, -1.4873, -1.3646, -1.2260, -1.0721, -0.9281,   &
     & -0.8379, -0.8123, -0.8261, -0.8483, -0.8305, -0.7792, -0.7626,   &
     & -0.8228, -0.9908, -1.2503, -1.5347, -1.7934, -1.9837, -2.0715,   &
     & -2.0375, -1.8975, -1.6906, -1.4497, -1.2048, -0.9831, -0.8125/
      DATA (CPCO2(IND),IND= 542, 553)/                                  &
     & -0.7157, -0.6707, -0.6532, -0.6297, -0.5706, -0.5263, -0.5489,   &
     & -0.6857, -0.9793, -1.3962, -1.8673, -2.3655/
!=CO2 ==== 3760- 4065
      DATA (CPCO2(IND),IND= 554, 615)/                                  &
     & -3.5436, -4.0424, -4.4084, -4.6848, -4.8663, -4.9516, -4.9790,   &
     & -4.9923, -5.0207, -5.0596, -5.0958, -5.1018, -5.0636, -5.0354,   &
     & -5.0546, -5.1454, -5.3274, -5.5863, -5.8889, -6.1770, -6.3555,   &
     & -6.4096, -6.4371, -6.5112, -6.6680, -6.9183, -7.2418, -7.5827,   &
     & -7.8704, -8.0551, -8.1705, -8.2500, -8.3554, -8.3961, -8.4354,   &
     & -8.3920, -8.2785, -8.0499, -7.7437, -7.4130, -7.1153, -6.8861,   &
     & -6.7422, -6.6786, -6.6774, -6.7053, -6.7090, -6.6794, -6.6055,   &
     & -6.4827, -6.3454, -6.2401, -6.1992, -6.2676, -6.4833, -6.8490,   &
     & -7.4310, -8.4606, -9.7364, -9.8771, -9.8840, -9.9559/
!=CO2 ==== 4530- 5380
      DATA (CPCO2(IND),IND= 616, 741)/                                  &
     & -9.9489, -9.6003, -9.0910, -8.5793, -8.2059, -7.9099, -7.7157,   &
     & -7.6145, -7.5964, -7.5942, -7.5256, -7.3190, -6.9986, -6.6884,   &
     & -6.4102, -6.1769, -5.9882, -5.8421, -5.7499, -5.7201, -5.7189,   &
     & -5.7108, -5.6669, -5.5955, -5.5686, -5.6287, -5.8000, -6.0855,   &
     & -6.4398, -6.7793, -6.9427, -6.9205, -6.8363, -6.7059, -6.5272,   &
     & -6.2903, -6.0085, -5.7224, -5.4722, -5.2772, -5.1501, -5.0768,   &
     & -5.0219, -4.9579, -4.8555, -4.7213, -4.5868, -4.4594, -4.3387,   &
     & -4.2219, -4.1002, -3.9812, -3.8876, -3.8207, -3.7673, -3.7120,   &
     & -3.6223, -3.4912, -3.3444, -3.1983, -3.0732, -3.0262, -3.0078,   &
     & -3.0123, -3.0213, -2.9957, -2.9261, -2.8770, -2.8887, -2.9853,   &
     & -3.1609, -3.3643, -3.5468, -3.6759, -3.7488, -3.7704, -3.7535,   &
     & -3.7113, -3.6368, -3.5277, -3.3812, -3.2020, -3.0043, -2.8020,   &
     & -2.6122, -2.4524, -2.3405, -2.2838, -2.2521, -2.2319, -2.1960,   &
     & -2.1562, -2.1732, -2.2913, -2.5476, -2.9382, -3.3966, -3.8525,   &
     & -4.2541, -4.5682, -4.7376, -4.7524, -4.6733, -4.5170, -4.3123,   &
     & -4.0891, -3.8565, -3.6218, -3.3909, -3.1785, -3.0100, -2.9105,   &
     & -2.8588, -2.8286, -2.7912, -2.7207, -2.6729, -2.6858, -2.7745,   &
     & -2.9414, -3.1445, -3.3617, -3.5954, -3.8508, -4.1739, -4.5122/
      DATA (CPCO2(IND),IND= 742, 786)/                                  &
     & -4.8985, -5.3426, -5.8737, -6.4734, -7.0715, -7.5042, -7.6034,   &
     & -7.5143, -7.4358, -7.4089, -7.3969, -7.3813, -7.3018, -7.1858,   &
     & -7.0633, -6.9962, -6.9905, -7.0319, -7.1331, -7.2054, -7.1856,   &
     & -7.0561, -6.7966, -6.4771, -6.1996, -5.9593, -5.7560, -5.5370,   &
     & -5.2836, -5.0966, -4.9583, -4.9126, -5.0022, -5.1370, -5.3465,   &
     & -5.6279, -5.9364, -6.3695, -6.9602, -7.6823, -8.2701, -8.6427,   &
     & -9.0728, -9.5366, -9.9588/
!=CO2 ==== 5905- 7025
      DATA (CPCO2(IND),IND= 787, 912)/                                  &
     & -9.9871, -9.6762, -9.3358, -8.9954, -8.5140, -8.2066, -7.9742,   &
     & -7.8579, -7.8073, -7.7894, -7.7466, -7.7009, -7.6393, -7.5889,   &
     & -7.5697, -7.5200, -7.3908, -7.1796, -6.9610, -6.7869, -6.6972,   &
     & -6.6735, -6.6775, -6.6495, -6.5292, -6.3435, -6.1371, -5.9268,   &
     & -5.7254, -5.5433, -5.4023, -5.3292, -5.3090, -5.3171, -5.3193,   &
     & -5.2705, -5.2085, -5.1835, -5.2186, -5.3367, -5.5305, -5.7725,   &
     & -6.0228, -6.2150, -6.2857, -6.2634, -6.2250, -6.2234, -6.2616,   &
     & -6.2931, -6.2508, -6.0971, -5.8679, -5.6195, -5.3906, -5.1944,   &
     & -5.0216, -4.8566, -4.6919, -4.5255, -4.3785, -4.2879, -4.2583,   &
     & -4.2636, -4.2768, -4.2484, -4.1853, -4.1586, -4.2079, -4.3651,   &
     & -4.6407, -5.0141, -5.4719, -6.0015, -6.5173, -6.7829, -6.6805,   &
     & -6.4180, -6.0793, -5.7404, -5.4204, -5.1265, -4.8634, -4.6378,   &
     & -4.4559, -4.3360, -4.2752, -4.2461, -4.2257, -4.1768, -4.1068,   &
     & -4.0743, -4.1193, -4.2732, -4.5464, -4.9256, -5.4090, -6.0184,   &
     & -6.7985, -7.7078, -8.3457, -8.5160, -8.6106, -8.8175, -9.1922,   &
     & -9.6775, -9.7423, -9.1980, -8.4120, -7.7499, -7.1685, -6.6817,   &
     & -6.2701, -5.9301, -5.6567, -5.4521, -5.3289, -5.2776, -5.2630,   &
     & -5.2547, -5.2083, -5.1296, -5.0823, -5.0914, -5.1806, -5.3503/
      DATA (CPCO2(IND),IND= 913,1011)/                                  &
     & -5.5600, -5.7877, -5.9936, -6.1720, -6.3801, -6.6371, -6.9964,   &
     & -7.5010, -8.1628, -8.9951, -9.8931,-10.0000,-10.0000,-10.0000,   &
     &-10.0000,-10.0000,-10.0000,-10.0000,-10.0000, -9.4967, -8.9198,   &
     & -8.5081, -8.1255, -7.8286, -7.5478, -7.1487, -6.7853, -6.5537,   &
     & -6.3931, -6.4107, -6.5087, -6.6607, -6.9026, -7.2104, -7.4445,   &
     & -7.6303, -7.6346, -7.4521, -7.2211, -7.0043, -6.7903, -6.5666,   &
     & -6.3499, -6.1534, -5.9988, -5.9033, -5.8760, -5.8693, -5.8277,   &
     & -5.7282, -5.6262, -5.5865, -5.6665, -5.9228, -6.3399, -7.0180,   &
     & -8.4230,-10.0000,-10.0000,-10.0000, -9.4090, -8.8272, -8.3057,   &
     & -7.8885, -7.5044, -7.1560, -6.8292, -6.5250, -6.2461, -5.9904,   &
     & -5.7533, -5.5295, -5.3135, -5.1058, -4.9152, -4.7463, -4.6054,   &
     & -4.4937, -4.3928, -4.2838, -4.1626, -4.0387, -3.9295, -3.8612,   &
     & -3.8501, -3.8647, -3.8625, -3.8099, -3.7351, -3.7179, -3.8549,   &
     & -4.2312, -4.7632, -5.4270, -6.4200, -8.1414, -9.0451, -9.5326,   &
     & -9.8301/
!=CO2 ==== 7395- 7785, 8030- 8335, 9340- 9670
      DATA (CPCO2(IND),IND=1012,1137)/                                  &
     & -9.9472, -9.8274, -8.9797, -8.4298, -7.8906, -7.4477, -7.0750,   &
     & -6.7698, -6.5338, -6.3739, -6.2980, -6.2739, -6.2726, -6.2555,   &
     & -6.1989, -6.1529, -6.1654, -6.2584, -6.4610, -6.7805, -7.2235,   &
     & -7.8191, -8.5850, -9.6084,-10.0000,-10.0000, -9.9199, -9.1093,   &
     & -8.4490, -7.9158, -7.4364, -7.0400, -6.6958, -6.4131, -6.1855,   &
     & -6.0158, -5.9123, -5.8700, -5.8530, -5.8340, -5.7866, -5.7224,   &
     & -5.7048, -5.7653, -5.9281, -6.2234, -6.6646, -7.2957, -8.2799,   &
     & -9.9457,-10.0000,-10.0000,-10.0000,-10.0000,-10.0000,-10.0000,   &
     &-10.0000, -9.2766, -8.6201, -8.0764, -7.6374, -7.2752, -6.9802,   &
     & -6.7578, -6.6163, -6.5546, -6.5392, -6.5397, -6.5132, -6.4531,   &
     & -6.4161, -6.4482, -6.5683, -6.8086, -7.1762, -7.6772, -8.3574,   &
     & -9.2188,-10.0000,-10.0000, -9.5350, -8.9686, -8.5329, -8.1920,   &
     & -7.9237, -7.6797, -7.5039, -7.3667, -7.2856, -7.1969, -7.0745,   &
     & -6.9330, -6.7926, -6.6818, -6.6144, -6.5643, -6.5183, -6.4910,   &
     & -6.4481, -6.3567, -6.2177, -6.0566, -5.9096, -5.7975, -5.7093,   &
     & -5.6165, -5.5127, -5.4124, -5.3426, -5.3061, -5.2648, -5.1864,   &
     & -5.0876, -5.0226, -5.0397, -5.1905, -5.4858, -5.9101, -6.4851,   &
     & -6.7862, -6.5368, -6.2765, -6.0398, -5.8260, -5.6397, -5.4799/
      DATA (CPCO2(IND),IND=1138,1219)/                                  &
     & -5.3438, -5.2274, -5.1411, -5.0917, -5.0473, -4.9820, -4.9114,   &
     & -4.8634, -4.8844, -5.0363, -5.3351, -5.7802, -6.5387, -8.3735,   &
     & -9.9977, -9.7506, -9.1887, -8.6824, -8.3488, -8.0533, -7.8664,   &
     & -7.7346, -7.6934, -7.6674, -7.6268, -7.5451, -7.4677, -7.4520,   &
     & -7.5471, -7.7913, -8.1917, -8.8835,-10.0000,-10.0000,-10.0000,   &
     &-10.0000,-10.0000, -9.7234, -8.9969, -8.5776, -8.1737, -7.8640,   &
     & -7.5729, -7.3186, -7.0973, -6.9131, -6.7782, -6.7073, -6.6768,   &
     & -6.6303, -6.5406, -6.4509, -6.3950, -6.4345, -6.6270, -6.9507,   &
     & -7.5028, -8.6428,-10.0000,-10.0000,-10.0000,-10.0000, -9.5303,   &
     & -8.9369, -8.4952, -8.1465, -7.8567, -7.6177, -7.4249, -7.2876,   &
     & -7.2206, -7.1948, -7.1552, -7.0773, -6.9884, -6.9402, -6.9839,   &
     & -7.1773, -7.4999, -8.0643, -9.1480,-10.0000/
!=CO  ====C' FOR    2 BAND MODEL
!=CO  ====    0-  175
      DATA (CPCO(IND),IND=  1, 36)/                                     &
     & -4.6868, -4.4127, -3.9461, -3.5662, -3.2921, -3.1081, -2.9807,   &
     & -2.8977, -2.8580, -2.8461, -2.8587, -2.9029, -2.9646, -3.0480,   &
     & -3.1589, -3.2836, -3.4277, -3.5993, -3.7963, -4.0164, -4.2799,   &
     & -4.5750, -4.8722, -5.2741, -5.6819, -6.0799, -6.4828, -6.8857,   &
     & -7.2886, -7.6915, -8.0944, -8.4973, -8.9002, -9.3031, -9.7060,   &
     &-10.0000/
!=CO  ==== 1940- 2285, 4040- 4370
      DATA (CPCO(IND),IND= 37,162)/                                     &
     &-10.0000, -9.5312, -8.8977, -8.2642, -7.5767, -6.9972, -6.5408,   &
     & -6.1219, -5.6734, -5.2658, -4.8686, -4.4918, -4.1423, -3.8133,   &
     & -3.4998, -3.2104, -2.9443, -2.7138, -2.5084, -2.3109, -2.1245,   &
     & -1.9387, -1.7608, -1.6054, -1.4733, -1.3594, -1.2540, -1.1480,   &
     & -1.0341, -0.9216, -0.8189, -0.7235, -0.6362, -0.5549, -0.4856,   &
     & -0.4401, -0.4268, -0.4657, -0.5571, -0.6573, -0.7404, -0.7523,   &
     & -0.6601, -0.5380, -0.4211, -0.3367, -0.3167, -0.3320, -0.3753,   &
     & -0.4489, -0.5438, -0.6653, -0.8052, -0.9690, -1.1506, -1.3522,   &
     & -1.5791, -1.8248, -2.1073, -2.4246, -2.7877, -3.2152, -3.7089,   &
     & -4.2832, -4.9518, -5.7251, -6.5319, -7.4879, -9.0885,-10.0000,   &
     &-10.0000, -9.5611, -9.0875, -8.6139, -7.9747, -7.5250, -7.1931,   &
     & -6.8596, -6.5741, -6.2922, -6.0098, -5.7669, -5.5345, -5.3229,   &
     & -5.1461, -4.9882, -4.8493, -4.7239, -4.6064, -4.5009, -4.4071,   &
     & -4.3322, -4.2661, -4.1926, -4.0956, -3.9611, -3.7984, -3.6314,   &
     & -3.4757, -3.3408, -3.2237, -3.1219, -3.0325, -2.9494, -2.8765,   &
     & -2.8117, -2.7531, -2.7023, -2.6635, -2.6440, -2.6550, -2.7225,   &
     & -2.8161, -2.9015, -2.9241, -2.8228, -2.6726, -2.5320, -2.4291,   &
     & -2.3772, -2.3732, -2.3995, -2.4574, -2.5486, -2.6664, -2.8209/
      DATA (CPCO(IND),IND=163,173)/                                     &
     & -3.0129, -3.2516, -3.5482, -3.9165, -4.3714, -4.9326, -5.6394,   &
     & -6.5163, -7.6063, -9.3575,-10.0000/
!=CH4 ====C' FOR    1 BAND MODEL
!=CH4 ==== 1065- 1775, 2345- 3230, 4110- 4690, 5865- 6135
      DATA (CPCH4(IND),IND=  1,126)/                                    &
     &-10.0000, -9.4577, -8.8866, -8.2246, -7.7940, -7.1734, -6.7965,   &
     & -6.5695, -6.1929, -5.9169, -5.7452, -5.4731, -5.3001, -5.1872,   &
     & -4.9672, -4.8474, -4.6939, -4.5210, -4.3377, -4.1346, -3.9322,   &
     & -3.7339, -3.5077, -3.2719, -3.0296, -2.8124, -2.6199, -2.4479,   &
     & -2.2502, -2.0541, -1.8800, -1.7092, -1.5791, -1.4379, -1.2992,   &
     & -1.1735, -1.0510, -0.9646, -0.8779, -0.8002, -0.7574, -0.7356,   &
     & -0.7478, -0.7512, -0.6906, -0.5594, -0.4417, -0.4019, -0.5027,   &
     & -0.7628, -0.9625, -1.0431, -1.0068, -0.8781, -0.7559, -0.6628,   &
     & -0.6128, -0.6118, -0.6575, -0.7620, -0.9217, -1.1264, -1.3660,   &
     & -1.6352, -1.9264, -2.2266, -2.5123, -2.7472, -2.8820, -2.9129,   &
     & -2.9145, -2.8854, -2.8508, -2.8512, -2.8202, -2.8023, -2.8004,   &
     & -2.7800, -2.8175, -2.8413, -2.8943, -2.9876, -3.0688, -3.2424,   &
     & -3.4064, -3.5759, -3.7630, -3.8925, -4.0774, -4.3243, -4.5964,   &
     & -3.8654, -3.0974, -2.5967, -2.2482, -2.1016, -2.1488, -2.3261,   &
     & -2.6448, -3.0446, -3.3958, -3.6510, -3.7049, -3.7240, -3.5992,   &
     & -3.4937, -3.3676, -3.2230, -3.1630, -3.0691, -3.0776, -3.0872,   &
     & -3.0974, -3.1223, -3.1285, -3.1212, -3.1333, -3.1674, -3.1668,   &
     & -3.2433, -3.2398, -3.3135, -3.3975, -3.4427, -3.6434, -3.7528/
      DATA (CPCH4(IND),IND=127,252)/                                    &
     & -3.9466, -4.1940, -4.3362, -4.5539, -4.7410, -4.9155, -5.1345,   &
     & -5.3908, -5.5592, -5.8270, -6.0289, -6.2365, -6.6730, -7.0538,   &
     & -7.6216, -8.5697, -9.8483,-10.0000, -9.3577, -8.5950, -7.8323,   &
     & -7.0696, -6.3069, -5.5442, -5.1501, -4.8853, -4.6900, -4.5262,   &
     & -4.3957, -4.2823, -4.2736, -4.2054, -4.1168, -3.9986, -3.8712,   &
     & -3.8692, -3.8777, -3.8965, -3.9092, -3.8788, -3.7661, -3.6900,   &
     & -3.6239, -3.5597, -3.5193, -3.4906, -3.4415, -3.3730, -3.3579,   &
     & -3.3427, -3.3208, -3.3048, -3.3136, -3.2904, -3.2545, -3.2241,   &
     & -3.1453, -3.0187, -2.9427, -2.8630, -2.8146, -2.8604, -2.8922,   &
     & -2.9650, -2.9959, -2.8920, -2.7989, -2.7028, -2.6506, -2.7285,   &
     & -2.8420, -2.9304, -2.9622, -2.8726, -2.7566, -2.6745, -2.6337,   &
     & -2.6533, -2.6800, -2.7098, -2.7479, -2.6859, -2.6216, -2.5701,   &
     & -2.4683, -2.4426, -2.4463, -2.4194, -2.4578, -2.4894, -2.4639,   &
     & -2.4825, -2.4998, -2.4381, -2.4123, -2.3654, -2.2698, -2.2387,   &
     & -2.2364, -2.2029, -2.1780, -2.1433, -2.0355, -1.9458, -1.8723,   &
     & -1.7936, -1.7639, -1.7782, -1.8022, -1.8115, -1.7818, -1.6986,   &
     & -1.6169, -1.5975, -1.6545, -1.7742, -1.8937, -1.9544, -1.8942,   &
     & -1.7761, -1.6392, -1.5236, -1.4551, -1.4221, -1.4245, -1.4174/
      DATA (CPCH4(IND),IND=253,378)/                                    &
     & -1.4177, -1.3776, -1.3349, -1.2909, -1.2470, -1.2162, -1.1850,   &
     & -1.1677, -1.1449, -1.1229, -1.1031, -1.0795, -1.0687, -1.0692,   &
     & -1.0904, -1.1166, -1.1511, -1.1951, -1.2321, -1.2831, -1.2716,   &
     & -1.1902, -0.9715, -0.6654, -0.4103, -0.3011, -0.5049, -0.8659,   &
     & -1.1777, -1.3847, -1.4359, -1.3908, -1.2992, -1.1923, -1.0951,   &
     & -1.0213, -0.9578, -0.9299, -0.9207, -0.9292, -0.9725, -1.0126,   &
     & -1.0750, -1.1149, -1.1636, -1.2059, -1.2638, -1.3327, -1.4079,   &
     & -1.4983, -1.5711, -1.6872, -1.7870, -1.9266, -2.0774, -2.2119,   &
     & -2.3875, -2.5155, -2.6822, -2.8372, -3.0032, -3.2413, -3.5058,   &
     & -3.9508, -4.5133, -5.3536, -8.0815, -8.9081, -9.8155,-10.0000,   &
     & -7.4757, -5.1602, -4.2454, -3.7640, -3.3256, -3.0103, -2.7726,   &
     & -2.5510, -2.3849, -2.2318, -2.1080, -2.0086, -1.9290, -1.8902,   &
     & -1.8750, -1.8700, -1.8476, -1.7390, -1.5724, -1.4284, -1.3425,   &
     & -1.3791, -1.5132, -1.6508, -1.7283, -1.6684, -1.5432, -1.4447,   &
     & -1.3773, -1.3490, -1.3642, -1.4016, -1.4713, -1.5836, -1.6984,   &
     & -1.8085, -1.8486, -1.7464, -1.6338, -1.5555, -1.5552, -1.6935,   &
     & -1.8165, -1.8417, -1.7697, -1.6346, -1.5589, -1.5466, -1.5604,   &
     & -1.6307, -1.6867, -1.7593, -1.8051, -1.8167, -1.8518, -1.8559/
      DATA (CPCH4(IND),IND=379,493)/                                    &
     & -1.8547, -1.8907, -1.8851, -1.8933, -1.9081, -1.9025, -1.9451,   &
     & -1.9924, -2.0321, -2.0816, -2.1026, -2.1137, -2.1351, -2.1629,   &
     & -2.1876, -2.2340, -2.2960, -2.3747, -2.4970, -2.6244, -2.7641,   &
     & -2.8912, -3.0328, -3.1944, -3.3877, -3.4566, -3.1662, -2.7253,   &
     & -2.3992, -2.2214, -2.2022, -2.3978, -2.7449, -3.2639, -3.9311,   &
     & -4.1470, -3.9351, -3.7471, -3.6245, -3.4791, -3.4710, -3.4210,   &
     & -3.4125, -3.4475, -3.4140, -3.4908, -3.5164, -3.5944, -3.7403,   &
     & -3.8192, -4.0177, -4.1833, -4.3518, -4.6486, -4.8778, -5.2542,   &
     & -5.7834, -6.3451, -7.7212,-10.0000, -9.9134, -7.9181, -6.0815,   &
     & -5.4397, -4.9875, -4.6154, -4.4846, -4.3541, -4.3037, -4.3073,   &
     & -4.2471, -4.2593, -4.1984, -4.1895, -4.1697, -4.1578, -4.1950,   &
     & -4.1878, -4.2299, -4.2209, -4.2646, -4.3123, -4.3911, -4.4588,   &
     & -4.1873, -3.8353, -3.5282, -3.3055, -3.3351, -3.5671, -3.8750,   &
     & -4.2645, -4.4786, -4.4293, -4.3183, -4.1996, -4.0879, -4.0169,   &
     & -3.9787, -3.9536, -3.9454, -3.9283, -3.9166, -3.9152, -3.9336,   &
     & -3.9561, -3.9932, -4.0934, -4.2317, -4.5084, -4.9460, -5.4958,   &
     & -6.5492, -8.5604, -9.6202/
!=N2O ====C' FOR    3 BAND MODEL
!=N2O ====    0-  120
      DATA (CPN2O(IND),IND=  1, 25)/                                    &
     & -2.8003, -2.6628, -2.4313, -2.2579, -2.1700, -2.1702, -2.2490,   &
     & -2.4003, -2.6264, -2.9219, -3.2954, -3.7684, -4.2621, -4.7558,   &
     & -5.2495, -5.7432, -6.2369, -6.7306, -7.2243, -7.7180, -8.2117,   &
     & -8.7054, -9.1991, -9.6928,-10.0000/
!=N2O ====  490-  775,  865-  995, 1065- 1385, 1545- 2040, 2090- 2655
      DATA (CPN2O(IND),IND= 26,151)/                                    &
     & -9.7185, -8.8926, -8.0667, -7.2307, -6.4149, -5.4872, -4.7083,   &
     & -4.0319, -3.4752, -3.0155, -2.6046, -2.2057, -1.8137, -1.4741,   &
     & -1.1914, -0.9603, -0.7923, -0.6629, -0.5849, -0.5402, -0.4975,   &
     & -0.5148, -0.5592, -0.6521, -0.8148, -1.0186, -1.2764, -1.5873,   &
     & -1.9638, -2.3881, -2.8083, -3.2392, -3.6934, -4.0682, -4.1366,   &
     & -3.9423, -3.7143, -3.4975, -3.2602, -3.0976, -2.9815, -2.9153,   &
     & -2.9596, -3.0281, -3.1264, -3.2650, -3.3906, -3.5717, -3.8312,   &
     & -4.1706, -4.6077, -5.1839, -5.9224, -6.9862, -7.6901, -8.3940,   &
     & -9.0979, -9.8018, -9.9154, -9.2271, -8.5388, -7.8504, -7.1621,   &
     & -6.2428, -5.6051, -5.0971, -4.7237, -4.4104, -4.2050, -4.0681,   &
     & -4.0278, -4.0307, -4.0492, -4.0333, -3.9710, -3.9249, -3.9360,   &
     & -4.0316, -4.2317, -4.5414, -4.9787, -5.5623, -6.3335, -7.9968,   &
     & -9.6601, -9.5486, -8.8517, -8.1548, -7.4579, -6.7610, -6.0641,   &
     & -5.3672, -4.6703, -3.6918, -3.0656, -2.5796, -2.1876, -1.8646,   &
     & -1.5919, -1.3587, -1.1684, -1.0286, -0.9470, -0.9271, -0.9442,   &
     & -0.9695, -0.9753, -0.9573, -0.9550, -1.0000, -1.1070, -1.2791,   &
     & -1.4976, -1.7281, -1.9277, -2.0227, -1.9577, -1.7625, -1.5020,   &
     & -1.2186, -0.9270, -0.6326, -0.3429, -0.0768,  0.1500,  0.3215/
      DATA (CPN2O(IND),IND=152,277)/                                    &
     &  0.4104,  0.4385,  0.4288,  0.4185,  0.4570,  0.4972,  0.4987,   &
     &  0.4216,  0.2360, -0.0319, -0.3714, -0.7539, -1.1534, -1.5855,   &
     & -2.0610, -2.6068, -3.2635, -4.1038, -5.2761, -6.1437, -7.0079,   &
     & -7.9440, -8.8801, -9.8162,-10.0000, -9.5951, -9.1305, -8.6659,   &
     & -8.2013, -7.7367, -7.2721, -6.8075, -6.1598, -5.8695, -5.3510,   &
     & -4.9491, -4.6310, -4.3846, -4.0784, -3.7763, -3.5901, -3.4607,   &
     & -3.4386, -3.5481, -3.7014, -3.9310, -4.2251, -4.4593, -4.8210,   &
     & -5.3494, -6.1286, -7.5981,-10.0000,-10.0000,-10.0000,-10.0000,   &
     & -6.3743, -5.5592, -5.0129, -4.6075, -4.3171, -4.0928, -3.7537,   &
     & -3.5406, -3.3869, -3.2913, -3.3633, -3.4932, -3.6924, -4.0074,   &
     & -4.2504, -4.5389, -4.9425, -5.4741, -6.2069, -7.5981,-10.0000,   &
     &-10.0000,-10.0000, -6.9215, -6.0798, -5.1934, -4.6288, -4.1316,   &
     & -3.7322, -3.4089, -3.1573, -2.9573, -2.7298, -2.5615, -2.4382,   &
     & -2.3523, -2.3774, -2.4508, -2.5755, -2.7757, -2.9904, -3.2733,   &
     & -3.6524, -4.1599, -4.7952, -5.7004, -6.8762, -6.9822, -6.2484,   &
     & -5.7613, -5.2586, -4.8674, -4.6633, -4.5332, -4.5158, -4.6593,   &
     & -4.8427, -5.0917, -5.5781, -6.0645, -6.5509, -7.0373, -7.5237,   &
     & -8.0101, -8.4965, -8.9829, -9.4693, -9.9557, -9.7130, -8.6609/
      DATA (CPN2O(IND),IND=278,389)/                                    &
     & -7.6089, -6.5568, -5.0880, -4.4527, -3.9302, -3.4438, -2.9701,   &
     & -2.5423, -2.1616, -1.8076, -1.4763, -1.1580, -0.8445, -0.5455,   &
     & -0.2506,  0.0234,  0.2775,  0.5113,  0.7154,  0.8929,  1.0359,   &
     &  1.1306,  1.1697,  1.1807,  1.1803,  1.1974,  1.2466,  1.2629,   &
     &  1.2068,  1.0472,  0.7695,  0.4083, -0.0244, -0.5477, -1.2202,   &
     & -2.1067, -2.9508, -3.2107, -3.1587, -2.9600, -2.7641, -2.6324,   &
     & -2.5671, -2.5664, -2.6088, -2.6425, -2.6606, -2.6895, -2.7551,   &
     & -2.8837, -3.0884, -3.3746, -3.7078, -4.0975, -4.6272, -5.2484,   &
     &-10.0000,-10.0000,-10.0000, -7.3571, -5.0287, -4.3047, -3.6431,   &
     & -3.1026, -2.6122, -2.1941, -1.8454, -1.5726, -1.3829, -1.2818,   &
     & -1.2505, -1.2579, -1.2731, -1.2502, -1.2092, -1.2044, -1.2577,   &
     & -1.3942, -1.6262, -1.9347, -2.2830, -2.5386, -2.4801, -2.1671,   &
     & -1.8061, -1.4726, -1.1797, -0.9377, -0.7542, -0.6392, -0.5899,   &
     & -0.5743, -0.5669, -0.5339, -0.4745, -0.4471, -0.4779, -0.5877,   &
     & -0.7964, -1.0942, -1.4812, -1.9593, -2.5140, -3.1350, -3.8102,   &
     & -4.5825, -5.5982, -6.4193, -7.2403, -8.0614, -8.8825, -9.7035/
!=N2O ==== 2705- 2865, 3245- 3925, 4260- 4470, 4540- 4785, 4910- 5165
      DATA (CPN2O(IND),IND=390,515)/                                    &
     & -9.8910, -8.9876, -8.0843, -7.1809, -6.1501, -5.3742, -4.7352,   &
     & -4.2051, -3.7525, -3.3562, -2.9916, -2.6649, -2.3872, -2.1499,   &
     & -1.9747, -1.7982, -1.6518, -1.5582, -1.4838, -1.5004, -1.5821,   &
     & -1.6912, -1.8673, -2.0756, -2.3351, -2.7020, -3.1921, -3.8409,   &
     & -4.7085, -5.9588, -6.5829, -8.5585, -9.8584, -9.9723, -9.4215,   &
     & -8.8707, -8.3199, -7.7691, -7.2183, -6.5567, -6.4345, -5.6448,   &
     & -5.0529, -4.4643, -3.9624, -3.5231, -3.1395, -2.8067, -2.5232,   &
     & -2.2858, -2.0820, -1.9049, -1.7554, -1.6485, -1.5959, -1.5838,   &
     & -1.5961, -1.5997, -1.5734, -1.5615, -1.5974, -1.7059, -1.9034,   &
     & -2.1631, -2.4181, -2.5427, -2.4592, -2.2513, -2.0187, -1.7879,   &
     & -1.5612, -1.3399, -1.1265, -0.9226, -0.7379, -0.5790, -0.4573,   &
     & -0.3952, -0.3683, -0.3511, -0.3216, -0.2556, -0.2126, -0.2593,   &
     & -0.4361, -0.7702, -1.2089, -1.7060, -2.2937, -3.1133, -4.4419,   &
     & -6.0119, -6.9457,-10.0000,-10.0000,-10.0000,-10.0000, -7.0394,   &
     & -5.9637, -5.2317, -4.6419, -4.1663, -3.7874, -3.5000, -3.3086,   &
     & -3.2143, -3.1926, -3.2105, -3.2308, -3.1971, -3.1510, -3.1402,   &
     & -3.1969, -3.3477, -3.6005, -3.9534, -4.4117, -4.9729, -5.6009,   &
     & -6.2179, -5.9845, -5.5502, -4.9010, -4.3401, -3.8232, -3.3802/
      DATA (CPN2O(IND),IND=516,641)/                                    &
     & -2.9972, -2.6747, -2.4143, -2.2209, -2.1080, -2.0682, -2.0687,   &
     & -2.0775, -2.0485, -1.9847, -1.9531, -1.9870, -2.1110, -2.3366,   &
     & -2.6293, -2.8922, -2.9474, -2.7627, -2.4999, -2.2554, -2.0537,   &
     & -1.9062, -1.8268, -1.7941, -1.7766, -1.7468, -1.6767, -1.6130,   &
     & -1.6085, -1.6849, -1.8599, -2.1258, -2.4538, -2.8205, -3.2028,   &
     & -3.5988, -4.0691, -4.7117, -5.6320, -6.4806, -7.3731, -8.2602,   &
     & -9.1474,-10.0000,-10.0000, -9.5340, -9.0282, -8.5224, -8.0166,   &
     & -7.5109, -7.0051, -6.4117, -6.0148, -5.4878, -5.1742, -4.8859,   &
     & -4.4873, -4.2249, -4.0285, -3.8669, -3.8247, -3.7652, -3.6521,   &
     & -3.4906, -3.2613, -3.0307, -2.8156, -2.6172, -2.4264, -2.2442,   &
     & -2.0775, -1.9432, -1.8703, -1.8523, -1.8552, -1.8443, -1.7814,   &
     & -1.7104, -1.7043, -1.7952, -2.0205, -2.3968, -2.9374, -3.7689,   &
     & -5.3159, -7.4139, -9.5119, -9.7965, -9.1511, -8.5057, -7.8603,   &
     & -7.2149, -6.5695, -6.2415, -5.5829, -5.0296, -4.5660, -4.1722,   &
     & -3.8364, -3.5551, -3.3398, -3.1970, -3.1363, -3.1232, -3.1257,   &
     & -3.0999, -3.0288, -2.9746, -2.9875, -3.0925, -3.3137, -3.6496,   &
     & -4.0276, -4.1958, -3.9760, -3.6179, -3.2725, -2.9653, -2.6962,   &
     & -2.4677, -2.2828, -2.1547, -2.0949, -2.0763, -2.0606, -2.0142/
      DATA (CPN2O(IND),IND=642,704)/                                    &
     & -1.9239, -1.8618, -1.8813, -2.0099, -2.2825, -2.7071, -3.3277,   &
     & -4.3300, -6.2151, -8.3543,-10.0000, -9.7275, -9.1257, -8.5239,   &
     & -7.9221, -7.3203, -6.7185, -6.6089, -5.8877, -5.4527, -5.0879,   &
     & -4.6598, -4.3806, -4.1830, -4.0426, -4.0175, -4.0178, -3.9811,   &
     & -3.9244, -3.8056, -3.6968, -3.6435, -3.6326, -3.6339, -3.6157,   &
     & -3.5478, -3.4826, -3.4807, -3.5665, -3.7650, -4.0718, -4.3980,   &
     & -4.5075, -4.3358, -4.0765, -3.8674, -3.7221, -3.6588, -3.6429,   &
     & -3.6371, -3.6014, -3.5209, -3.4616, -3.4774, -3.5957, -3.8481,   &
     & -4.2598, -4.8784, -5.8266, -6.7468, -8.1352, -9.2208,-10.0000/
!=O2  ====C' FOR    2 BAND MODEL
!=O2  ====    0-  265
      DATA (CPO2(IND),IND=  1, 54)/                                     &
     & -6.1363, -6.1794, -6.2538, -6.3705, -6.5110, -6.6162, -6.7505,   &
     & -6.7896, -6.8305, -6.8471, -6.8282, -6.8772, -6.8680, -6.9332,   &
     & -6.9511, -7.0048, -7.0662, -7.1043, -7.2055, -7.2443, -7.3520,   &
     & -7.4079, -7.4998, -7.5924, -7.6682, -7.7993, -7.8712, -8.0161,   &
     & -8.1102, -8.2485, -8.3758, -8.4942, -8.6532, -8.7554, -8.9453,   &
     & -9.0665, -9.2631, -9.4387, -9.6325, -9.8757,-10.0628,-10.3761,   &
     &-10.5478,-10.9147,-11.2052,-11.5129,-11.8206,-12.1283,-12.4360,   &
     &-12.7437,-13.0514,-13.3591,-13.6668,-13.9745/
!=O2  ==== 7650- 8080, 9235- 9490,12850-13220,14300-14600,15695-15955
      DATA (CPO2(IND),IND= 55,180)/                                     &
     &-13.9458,-13.7692,-13.5048,-13.1422,-13.0242,-12.6684,-12.3571,   &
     &-12.2428,-11.8492,-11.6427,-11.5173,-11.2108,-11.1584,-11.0196,   &
     &-10.8040,-10.8059,-10.5828,-10.4580,-10.4170,-10.1823,-10.1435,   &
     &-10.0030, -9.8136, -9.7772, -9.5680, -9.4595, -9.3502, -9.1411,   &
     & -9.0476, -8.8628, -8.7051, -8.5838, -8.4282, -8.3271, -8.1958,   &
     & -8.0838, -7.9652, -7.8371, -7.7476, -7.6431, -7.5736, -7.5149,   &
     & -7.4194, -7.2688, -7.0722, -6.8815, -6.7627, -6.8055, -6.9114,   &
     & -6.9936, -7.0519, -7.0597, -7.0680, -7.1242, -7.2088, -7.3265,   &
     & -7.4673, -7.6326, -7.8110, -8.0096, -8.2104, -8.4036, -8.5853,   &
     & -8.7252, -8.8511, -8.9427, -9.0375, -9.1228, -9.2246, -9.3291,   &
     & -9.4436, -9.5716, -9.6951, -9.8408, -9.9759,-10.1489,-10.3027,   &
     &-10.5178,-10.7265,-10.9787,-11.2939,-11.5552,-11.9595,-12.2436,   &
     &-12.6942,-13.2011,-13.8191,-13.9216,-13.7293,-13.5370,-13.3447,   &
     &-13.1523,-12.9600,-12.7677,-12.5754,-12.3830,-12.1907,-11.9948,   &
     &-11.7759,-11.5926,-11.4214,-11.2493,-11.1094,-10.9477,-10.8332,   &
     &-10.7323,-10.6380,-10.5725,-10.4409,-10.2013, -9.8839, -9.6546,   &
     & -9.5053, -9.4638, -9.5526, -9.6558, -9.7430, -9.7958, -9.7896,   &
     & -9.8320, -9.9447,-10.1221,-10.3707,-10.6623,-10.9761,-11.2271/
      DATA (CPO2(IND),IND=181,306)/                                     &
     &-11.4091,-11.4921,-11.6015,-11.6945,-11.8333,-11.9985,-12.1788,   &
     &-12.3822,-12.6605,-13.0796,-13.3528,-13.6463,-13.9398,-13.7034,   &
     &-13.3150,-13.1177,-12.6462,-12.4868,-12.2205,-11.9650,-11.6941,   &
     &-11.4377,-11.2136,-10.9567,-10.7980,-10.5546,-10.3952,-10.2403,   &
     &-10.0491, -9.9226, -9.7871, -9.6557, -9.6106, -9.5142, -9.4763,   &
     & -9.4163, -9.2348, -9.1088, -8.7946, -8.5876, -8.3128, -8.0945,   &
     & -7.9127, -7.7229, -7.5860, -7.4215, -7.2726, -7.1179, -6.9516,   &
     & -6.8075, -6.6413, -6.5043, -6.3519, -6.2112, -6.0839, -5.9337,   &
     & -5.8321, -5.6969, -5.5923, -5.5076, -5.4002, -5.3413, -5.2826,   &
     & -5.2458, -5.2877, -5.3743, -5.4654, -5.5262, -5.4429, -5.2430,   &
     & -5.0284, -4.8464, -4.7534, -4.7825, -4.9462, -5.2290, -5.6440,   &
     & -6.1889, -6.8427, -7.7731, -9.1688, -9.6893,-10.1853,-10.7670,   &
     &-11.4611,-12.3081,-13.1476,-13.8192,-13.5871,-13.2189,-12.9705,   &
     &-12.4825,-12.1301,-11.9430,-11.6636,-11.3197,-11.1678,-10.8967,   &
     &-10.6002,-10.4857,-10.1986, -9.9731, -9.8547, -9.5817, -9.4382,   &
     & -9.3042, -9.0755, -8.9944, -8.8060, -8.6543, -8.5441, -8.3556,   &
     & -8.2557, -8.0959, -7.9717, -7.8453, -7.7076, -7.5910, -7.4567,   &
     & -7.3439, -7.2248, -7.1236, -7.0209, -6.9345, -6.8404, -6.7560/
      DATA (CPO2(IND),IND=307,382)/                                     &
     & -6.6744, -6.5870, -6.5278, -6.4809, -6.5042, -6.5797, -6.6564,   &
     & -6.6939, -6.5912, -6.3776, -6.1438, -6.0062, -6.0469, -6.3081,   &
     & -6.8199, -7.4307, -8.1345, -9.1190,-10.4203,-11.4698,-12.5942,   &
     &-13.5316,-13.8693,-13.9392,-13.6885,-13.4377,-13.1869,-12.9362,   &
     &-12.6854,-12.3720,-12.2852,-11.9331,-11.7575,-11.6297,-11.3290,   &
     &-11.1205,-11.0084,-10.7243,-10.5543,-10.4485,-10.1764,-10.0759,   &
     & -9.9304, -9.7196, -9.6630, -9.4774, -9.3638, -9.2675, -9.1121,   &
     & -9.0368, -8.9025, -8.8028, -8.7012, -8.5909, -8.5121, -8.4141,   &
     & -8.3444, -8.2687, -8.2003, -8.1571, -8.1141, -8.1261, -8.1848,   &
     & -8.2395, -8.2478, -8.0877, -7.7880, -7.5611, -7.4487, -7.4880,   &
     & -7.7644, -8.2142, -8.8765,-10.1091,-12.4483,-13.7228/
      END