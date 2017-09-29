      SUBROUTINE ABCDTA(V)
      COMMON/ABC/ ANH3(2),  ACO2(10), ACO(3),    ACH4(4),   ANO2(3),    &
     &  AN2O(11), AO2(6),   AO3(5),   ASO2(4),   AH2O(14),  ANO,        &
     &  AANH3(2), BBNH3(2), CCNH3(2), AACO2(10), BBCO2(10), CCCO2(10),  &
     &  AACO(3),  BBCO(3),  CCCO(3),  AACH4(4),  BBCH4(4),  CCCH4(4),   &
     &  AANO2(3), BBNO2(3), CCNO2(3), AAN2O(11), BBN2O(11), CCN2O(11),  &
     &  AAO2(6),  BBO2(6),  CCO2(6),  AAO3(5),   BBO3(5),   CCO3(5),    &
     &  AASO2(4), BBSO2(4), CCSO2(4), AAH2O(14), BBH2O(14), CCH2O(14),  &
     &  AANO,     BBNO,     CCNO
      COMMON/AABBCC/AA(11),BB(11),CC(11),IBND(11),QA(11),CPS(11)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /ABC/
      EXTERNAL ABCD

!    MOL
!     1    H2O (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!     2    CO2 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!     3    O3  (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!     4    N2O (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!     5    CO  (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!     6    CH4 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!     7    O2  (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!     8    NO  (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!     9    SO2 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!    10    NO2 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)
!    11    NH3 (ALL REGIONS) (DOUBLE EXPONENTIAL MODELS)

!  ---H2O
      IMOL=1
      IW=-1
      IF(V.GE.     0..AND.V.LE.   345.) IW=17
      IF(V.GE.   350..AND.V.LE.  1000.) IW=18
      IF(V.GE.  1005..AND.V.LE.  1640.) IW=19
      IF(V.GE.  1645..AND.V.LE.  2530.) IW=20
      IF(V.GE.  2535..AND.V.LE.  3420.) IW=21
      IF(V.GE.  3425..AND.V.LE.  4310.) IW=22
      IF(V.GE.  4315..AND.V.LE.  6150.) IW=23
      IF(V.GE.  6155..AND.V.LE.  8000.) IW=24
      IF(V.GE.  8005..AND.V.LE.  9615.) IW=25
      IF(V.GE.  9620..AND.V.LE. 11540.) IW=26
      IF(V.GE. 11545..AND.V.LE. 13070.) IW=27
      IF(V.GE. 13075..AND.V.LE. 14860.) IW=28
      IF(V.GE. 14865..AND.V.LE. 16045.) IW=29
      IF(V.GE. 16340..AND.V.LE. 17860.) IW=30
      IBAND=IW - 16
      IBND(IMOL)=IW
      IF(IW .GT.  0) THEN
           QA(IMOL)= AH2O(IBAND)
           AA(IMOL) =AAH2O(IBAND)
           BB(IMOL) =BBH2O(IBAND)
           CC(IMOL) =CCH2O(IBAND)
      ENDIF
!  ---O3
      IMOL=3
      IW=-1
      IF (V .GE.     0. .AND. V .LE.   200.)  IW=31
      IF (V .GE.   515. .AND. V .LE.  1275.)  IW=32
      IF (V .GE.  1630. .AND. V .LE.  2295.)  IW=33
      IF (V .GE.  2670. .AND. V .LE.  2845.)  IW=34
      IF (V .GE.  2850. .AND. V .LE.  3260.)  IW=35
      IBAND     =IW - 30
      IBND(IMOL)=IW
      IF(IW .GT.  0) THEN
           QA(IMOL)=AO3(IBAND)
           AA(IMOL)=AAO3(IBAND)
           BB(IMOL)=BBO3(IBAND)
           CC(IMOL)=CCO3(IBAND)
      ENDIF
!  ---CO2
      IMOL=2
      IW=-1
      IF (V .GE.   425. .AND. V .LE.   835.)  IW=36
      IF (V .GE.   840. .AND. V .LE.  1440.)  IW=37
      IF (V .GE.  1805. .AND. V .LE.  2855.)  IW=38
      IF (V .GE.  3070. .AND. V .LE.  3755.)  IW=39
      IF (V .GE.  3760. .AND. V .LE.  4065.)  IW=40
      IF (V .GE.  4530. .AND. V .LE.  5380.)  IW=41
      IF (V .GE.  5905. .AND. V .LE.  7025.)  IW=42
      IF((V .GE.  7395. .AND. V .LE.  7785.) .OR.                       &
     &   (V .GE.  8030. .AND. V .LE.  8335.) .OR.                       &
     &   (V .GE.  9340. .AND. V .LE.  9670.)) IW=43
      IBAND=IW - 35
      IBND(IMOL)=IW
      IF(IW .GT.  0) THEN
           QA(IMOL)=ACO2(IBAND)
           AA(IMOL)=AACO2(IBAND)
           BB(IMOL)=BBCO2(IBAND)
           CC(IMOL)=CCCO2(IBAND)
      ENDIF
!  ---CO
      IMOL=5
      IW=-1
      IF (V .GE.     0. .AND. V .LE.   175.) IW=44
      IF((V .GE.  1940. .AND. V .LE.  2285.) .OR.                       &
     &   (V .GE.  4040. .AND. V .LE.  4370.)) IW=45
      IBAND=IW - 43
      IBND(IMOL)=IW
      IF(IW .GT.  0) THEN
           QA(IMOL)=ACO(IBAND)
           AA(IMOL)=AACO(IBAND)
           BB(IMOL)=BBCO(IBAND)
           CC(IMOL)=CCCO(IBAND)
      ENDIF
!  ---CH4
      IMOL=6
      IW=-1
      IF((V .GE.  1065. .AND. V .LE.  1775.) .OR.                       &
     &   (V .GE.  2345. .AND. V .LE.  3230.) .OR.                       &
     &   (V .GE.  4110. .AND. V .LE.  4690.) .OR.                       &
     &   (V .GE.  5865. .AND. V .LE.  6135.))IW=46
      IBAND=IW - 45
      IBND(IMOL)=IW
      IF(IW .GT.  0) THEN
           QA(IMOL)=ACH4(IBAND)
           AA(IMOL)=AACH4(IBAND)
           BB(IMOL)=BBCH4(IBAND)
           CC(IMOL)=CCCH4(IBAND)
      ENDIF
!  ---N2O
      IMOL=4
      IW=-1
      IF (V .GE.     0. .AND. V .LE.   120.)  IW=47
      IF((V .GE.   490. .AND. V .LE.   775.) .OR.                       &
     &   (V .GE.   865. .AND. V .LE.   995.) .OR.                       &
     &   (V .GE.  1065. .AND. V .LE.  1385.) .OR.                       &
     &   (V .GE.  1545. .AND. V .LE.  2040.) .OR.                       &
     &   (V .GE.  2090. .AND. V .LE.  2655.)) IW=48
      IF((V .GE.  2705. .AND. V .LE.  2865.) .OR.                       &
     &   (V .GE.  3245. .AND. V .LE.  3925.) .OR.                       &
     &   (V .GE.  4260. .AND. V .LE.  4470.) .OR.                       &
     &   (V .GE.  4540. .AND. V .LE.  4785.) .OR.                       &
     &   (V .GE.  4910. .AND. V .LE.  5165.)) IW=49
      IBAND=IW - 46
      IBND(IMOL)=IW

      IF(IW .EQ. 49)IBAND=7

!     THIS CORRECTION IS ONLY FOR N2O AS CURRENTLY WRITTEN

      IF(IW .GT.  0) THEN
           QA(IMOL)=AN2O(IBAND)
           AA(IMOL)=AAN2O(IBAND)
           BB(IMOL)=BBN2O(IBAND)
           CC(IMOL)=CCN2O(IBAND)
      ENDIF
!  ---O2
      IMOL=7
      IW=-1
      IF (V .GE.     0. .AND. V .LE.   265.)  IW=50
      IF((V .GE.  7650. .AND. V .LE.  8080.) .OR.                       &
     &   (V .GE.  9235. .AND. V .LE.  9490.) .OR.                       &
     &   (V .GE. 12850. .AND. V .LE. 13220.) .OR.                       &
     &   (V .GE. 14300. .AND. V .LE. 14600.) .OR.                       &
     &   (V .GE. 15695. .AND. V .LE. 15955.)) IW=51
       IF(V .GE. 49600. .AND. V. LE. 52710.)  IW=51
      IBAND=IW - 49
      IBND(IMOL)=IW
      IF(IW .GT.  0) THEN
           QA(IMOL)=AO2(IBAND)
           IF(V .GE. 49600. .AND. V. LE. 52710.)  QA(IMOL) =.4704
           AA(IMOL)=AAO2(IBAND)
           BB(IMOL)=BBO2(IBAND)
           CC(IMOL)=CCO2(IBAND)
      ENDIF
!  ---NH3
      IMOL=11
      IW=-1
      IF (V .GE.     0. .AND. V .LE.   385.)  IW=52
      IF (V .GE.   390. .AND. V .LE.  2150.)  IW=53
      IBAND=IW - 51
      IBND(IMOL)=IW
      IF(IW .GT.  0) THEN
           QA(IMOL)=ANH3(IBAND)
           AA(IMOL)=AANH3(IBAND)
           BB(IMOL)=BBNH3(IBAND)
           CC(IMOL)=CCNH3(IBAND)
      ENDIF
!  ---NO
      IMOL=8
      IW=-1
      IF (V .GE.  1700. .AND. V .LE.  2005.) IW =54
      IBAND=IW - 53
      IBND(IMOL)=IW
      IF(IW .GT.  0) THEN
           QA(IMOL)=ANO
           AA(IMOL)=AANO
           BB(IMOL)=BBNO
           CC(IMOL)=CCNO
      ENDIF
!  ---NO2
      IW=-1
      IMOL=10
      IF((V .GE.   580. .AND. V .LE.   925.) .OR.                       &
     &   (V .GE.  1515. .AND. V .LE.  1695.) .OR.                       &
     &   (V .GE.  2800. .AND. V .LE.  2970.)) IW=55
      IBAND=IW - 54
      IBND(IMOL)=IW
      IF(IW .GT.  0) THEN
           QA(IMOL)=ANO2(IBAND)
           AA(IMOL)=AANO2(IBAND)
           BB(IMOL)=BBNO2(IBAND)
           CC(IMOL)=CCNO2(IBAND)
      ENDIF
!  ---SO2
      IMOL=9
      IW=-1
      IF (V .GE.     0. .AND. V .LE.   185.)  IW=56
      IF((V .GE.   400. .AND. V .LE.   650.) .OR.                       &
     &   (V .GE.   950. .AND. V .LE.  1460.) .OR.                       &
     &   (V .GE.  2415. .AND. V .LE.  2580.)) IW=57
      IBAND=IW - 55
      IBND(IMOL)=IW
      IF(IW .GT.  0) THEN
           QA(IMOL)=ASO2(IBAND)
           AA(IMOL)=AASO2(IBAND)
           BB(IMOL)=BBSO2(IBAND)
           CC(IMOL)=CCSO2(IBAND)
      ENDIF
      RETURN
      END
