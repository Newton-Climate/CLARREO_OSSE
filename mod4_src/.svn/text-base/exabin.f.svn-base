      SUBROUTINE EXABIN(ICH)
      REAL A,A1,A2,ALTB,E1,E2,EC,WRH,X,X1,X2,Y,Y1,Y2,Z1,Z2,ZK
      INTEGER IRH,IM1RH,IREG,IREGC,ITA,ITAS,ITC,IAER,IWAV,ICH(4)

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     LOADS EXTINCTION, ABSORPTION AND ASYMMETRY COEFFICIENTS
!     FOR THE FOUR AEROSOL ALTITUDE REGIONS

!     MODIFIED FOR ASYMMETRY - JAN 1986 (A.E.R. INC.)

!     COMMONS:
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA
      REAL VIS,WSS,WHH,RAINRT
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 &
     &  VIS,WSS,WHH,RAINRT
      COMMON/CARD2D/IREG(14),ALTB(14),IREGC(14)
      INCLUDE 'BASE.h'

!           0-2KM
!             RUREXT=RURAL EXTINCTION   RURABS=RURAL ABSORPTION
!             RURSYM=RURAL ASYMMETRY FACTORS
!             URBEXT=URBAN EXTINCTION   URBABS=URBAN ABSORPTION
!             URBSYM=URBAN ASYMMETRY FACTORS
!             OCNEXT=MARITIME EXTINCTION  OCNABS=MARITIME ABSORPTION
!             OCNSYM=MARITIME ASYMMETRY FACTORS
!             TROEXT=TROPSPHER EXTINCTION  TROABS=TROPOSPHER ABSORPTION
!             TROSYM=TROPSPHERIC ASYMMETRY FACTORS
!             FG1EXT=FOG1 .2KM VIS EXTINCTION  FG1ABS=FOG1 ABSORPTION
!             FG1SYM=FOG1 ASYMMETRY FACTORS
!             FG2EXT=FOG2 .5KM VIS EXTINCTION  FG2ABS=FOG2 ABSORPTION
!             FG2SYM=FOG2 ASYMMETRY FACTORS
!           >2-10KM
!             TROEXT=TROPOSPHER EXTINCTION  TROABS=TROPOSPHER ABSORPTION
!             TROSYM=TROPOSPHERIC ASYMMETRY FACTORS
!           >10-30KM
!             BSTEXT=BACKGROUND STRATOSPHERIC EXTINCTION
!             BSTABS=BACKGROUND STRATOSPHERIC ABSORPTION
!             BSTSYM=BACKGROUND STRATOSPHERIC ASYMMETRY FACTORS
!             AVOEXT=AGED VOLCANIC EXTINCTION
!             AVOABS=AGED VOLCANIC ABSORPTION
!             AVOSYM=AGED VOLCANIC ASYMMETRY FACTORS
!             FVOEXT=FRESH VOLCANIC EXTINCTION
!             FVOABS=FRESH VOLCANIC ABSORPTION
!             FVOSYM=FRESH VOLCANIC ASYMMETRY FACTORS
!           >30-100KM
!             DMEEXT=METEORIC DUST EXTINCTION
!             DMEABS=METEORIC DUST ABSORPTION
!             DMESYM=METEORIC DUST ASYMMETRY FACTORS

!     AEROSOL EXTINCTION AND ABSORPTION DATA:

!     MODIFIED TO INCLUDE ASYMMETRY DATA - JAN 1986 (A.E.R. INC.)
      REAL                   RUREXT,RURABS,RURSYM,URBEXT,URBABS,URBSYM, &
     &  OCNEXT,OCNABS,OCNSYM,TROEXT,TROABS,TROSYM,FG1EXT,FG1ABS,FG1SYM, &
     &  FG2EXT,FG2ABS,FG2SYM,BSTEXT,BSTABS,BSTSYM,AVOEXT,AVOABS,AVOSYM, &
     &  FVOEXT,FVOABS,FVOSYM,DMEEXT,DMEABS,DMESYM
      COMMON/EXTD/RUREXT(NWAVLN,4),RURABS(NWAVLN,4),RURSYM(NWAVLN,4),   &
     &            URBEXT(NWAVLN,4),URBABS(NWAVLN,4),URBSYM(NWAVLN,4),   &
     &            OCNEXT(NWAVLN,4),OCNABS(NWAVLN,4),OCNSYM(NWAVLN,4),   &
     &            TROEXT(NWAVLN,4),TROABS(NWAVLN,4),TROSYM(NWAVLN,4),   &
     &            FG1EXT(NWAVLN),FG1ABS(NWAVLN),FG1SYM(NWAVLN),         &
     &            FG2EXT(NWAVLN),FG2ABS(NWAVLN),FG2SYM(NWAVLN),         &
     &            BSTEXT(NWAVLN),BSTABS(NWAVLN),BSTSYM(NWAVLN),         &
     &            AVOEXT(NWAVLN),AVOABS(NWAVLN),AVOSYM(NWAVLN),         &
     &            FVOEXT(NWAVLN),FVOABS(NWAVLN),FVOSYM(NWAVLN),         &
     &            DMEEXT(NWAVLN),DMEABS(NWAVLN),DMESYM(NWAVLN)
      REAL CCUEXT,CCUABS,CCUSYM,CALEXT,CALABS,CALSYM,                   &
     &  CSTEXT,CSTABS,CSTSYM,CSCEXT,CSCABS,CSCSYM,CNIEXT,CNIABS,CNISYM
      COMMON/CLDDAT/CCUEXT(NWAVLN),CCUABS(NWAVLN),CCUSYM(NWAVLN),       &
     &              CALEXT(NWAVLN),CALABS(NWAVLN),CALSYM(NWAVLN),       &
     &              CSTEXT(NWAVLN),CSTABS(NWAVLN),CSTSYM(NWAVLN),       &
     &              CSCEXT(NWAVLN),CSCABS(NWAVLN),CSCSYM(NWAVLN),       &
     &              CNIEXT(NWAVLN),CNIABS(NWAVLN),CNISYM(NWAVLN)
      REAL CI64XT,CI64AB,CI64G,CIR4XT,CIR4AB,CIR4G
      COMMON /CIRR  / CI64XT(NWAVLN), CI64AB(NWAVLN), CI64G(NWAVLN),    &
     &     CIR4XT(NWAVLN), CIR4AB(NWAVLN), CIR4G(NWAVLN)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /CLDDAT/,/CIRR/,/EXTD/
      EXTERNAL CLDDTA,EXTDTA

!     DATA:
      REAL AFLWC,ASLWC,AVLWC,BSLWC,CULWC,FVLWC,RFLWC,SCLWC,SNLWC,       &
     &  STLWC,MDLWC,RHZONE(4),ELWCR(4),ELWCU(4),ELWCM(4),ELWCT(4)
      DATA RHZONE/0., 70., 80., 99./
      DATA ELWCR/3.517E-04, 3.740E-04, 4.439E-04, 9.529E-04/
      DATA ELWCM/4.675E-04, 6.543E-04, 1.166E-03, 3.154E-03/
      DATA ELWCU/3.102E-04, 3.802E-04, 4.463E-04, 9.745E-04/
      DATA ELWCT/1.735E-04, 1.820E-04, 2.020E-04, 2.408E-04/
      DATA AFLWC/1.295E-02/, RFLWC/1.804E-03/, CULWC/7.683E-03/
      DATA ASLWC/4.509E-03/, STLWC/5.272E-03/, SCLWC/4.177E-03/
      DATA SNLWC/7.518E-03/, BSLWC/1.567E-04/, FVLWC/5.922E-04/
      DATA AVLWC/1.675E-04/, MDLWC/4.775E-04/

!     LOOP OVER AEROSOL REGIONS:
      DO 30 IAER=1,14  !DRF
          AWCCON(IAER)=0
          IF(IREG(IAER).EQ.0 .OR. (ICH(1).EQ.3 .AND. IAER.EQ.1))THEN
              ITA=ICH(IAER)
              ITC=ICH(IAER)-7
              ITAS=ITA
              IF(IREGC(IAER).NE.0)THEN

!                 SECTION TO LOAD EXTINCTION AND ABSORPTION
!                 COEFFICIENTS FOR CLOUD AND OR RAIN MODELS
                  DO 10 IWAV=1,NWAVLN-1
                      IF(ICLD.EQ.2)THEN
                          ABSC(IAER,IWAV)=CALABS(IWAV)
                          EXTC(IAER,IWAV)=CALEXT(IWAV)
                          ASYM(IAER,IWAV)=CALSYM(IWAV)
                          IF(IWAV.EQ.1)AWCCON(IAER)=ASLWC
                      ELSEIF(ICLD.EQ.3 .OR. ICLD.EQ.6)THEN
                          ABSC(IAER,IWAV)=CSTABS(IWAV)
                          EXTC(IAER,IWAV)=CSTEXT(IWAV)
                          ASYM(IAER,IWAV)=CSTSYM(IWAV)
                          IF(IWAV.EQ.1)AWCCON(IAER)=STLWC
                      ELSEIF(ICLD.EQ.4)THEN
                          ABSC(IAER,IWAV)=CSCABS(IWAV)
                          EXTC(IAER,IWAV)=CSCEXT(IWAV)
                          ASYM(IAER,IWAV)=CSCSYM(IWAV)
                          IF(IWAV.EQ.1)AWCCON(IAER)=SCLWC
                      ELSEIF(ICLD.EQ.5 .OR. ICLD.EQ.7                   &
     &                                 .OR. ICLD.EQ.8)THEN
                          ABSC(IAER,IWAV)=CNIABS(IWAV)
                          EXTC(IAER,IWAV)=CNIEXT(IWAV)
                          ASYM(IAER,IWAV)=CNISYM(IWAV)
                          IF(IWAV.EQ.1)AWCCON(IAER)=SNLWC
                      ELSE
                          ABSC(IAER,IWAV)=CCUABS(IWAV)
                          EXTC(IAER,IWAV)=CCUEXT(IWAV)
                          ASYM(IAER,IWAV)=CCUSYM(IWAV)
                          IF(IWAV.EQ.1)AWCCON(IAER)=CULWC
                      ENDIF
   10             CONTINUE
              ELSE
                  WRH=W(15)
                  IF(ICH(IAER).EQ.6 .AND. IAER.NE.1)WRH=70.

!                 THIS DOES NOT ALLOW TROP RH DEPENDENT ABOVE EH(7,I)
!                 DEFAULTS TO TROPOSPHERIC AT 70. PERCENT
                  IF(WRH.LT.RHZONE(2))THEN
                      IM1RH=1
                      IRH=2
                  ELSEIF(WRH.LT.RHZONE(3))THEN
                      IM1RH=2
                      IRH=3
                  ELSE
                      IM1RH=3
                      IRH=4
                  ENDIF
                  IF(WRH.GT.0. .AND. WRH.LT.99.)X=ALOG(100.-WRH)
                  X1=ALOG(100.-RHZONE(IM1RH))
                  X2=ALOG(100.-RHZONE(IRH))
                  IF(WRH.GE.99.)X=X2
                  IF(WRH.LE.0.)X=X1
                  DO 20 IWAV=1,NWAVLN-1
                      ITA=ITAS
                      IF(ITA.NE.3 .OR. IAER.NE.1)THEN
                          ABSC(IAER,IWAV)=0.
                          EXTC(IAER,IWAV)=0.
                          ASYM(IAER,IWAV)=0.
                          IF(ITA.GT.6)THEN
                              IF(ITA.GT.19)THEN
                                  ABSC(IAER,IWAV)=DMEABS(IWAV)
                                  EXTC(IAER,IWAV)=DMEEXT(IWAV)
                                  ASYM(IAER,IWAV)=DMESYM(IWAV)
                                  IF(IWAV.EQ.1)AWCCON(IAER)=MDLWC
                                  GOTO20
                              ENDIF
                              IF(ITC.GE.1)THEN
                                  IF(ITC.EQ.2)THEN
                                      ABSC(IAER,IWAV)=FG2ABS(IWAV)
                                      EXTC(IAER,IWAV)=FG2EXT(IWAV)
                                      ASYM(IAER,IWAV)=FG2SYM(IWAV)
                                      IF(IWAV.EQ.1)AWCCON(IAER)=RFLWC
                                  ELSEIF(ITC.EQ.4 .OR. ITC.EQ.9         &
     &                                            .OR. ITC.EQ.10)THEN
                                      ABSC(IAER,IWAV)=BSTABS(IWAV)
                                      EXTC(IAER,IWAV)=BSTEXT(IWAV)
                                      ASYM(IAER,IWAV)=BSTSYM(IWAV)
                                      IF(IWAV.EQ.1)AWCCON(IAER)=BSLWC
                                  ELSEIF(ITC.EQ.5 .OR. ITC.EQ.7)THEN
                                      ABSC(IAER,IWAV)=AVOABS(IWAV)
                                      EXTC(IAER,IWAV)=AVOEXT(IWAV)
                                      ASYM(IAER,IWAV)=AVOSYM(IWAV)
                                      IF(IWAV.EQ.1)AWCCON(IAER)=AVLWC
                                  ELSEIF(ITC.EQ.6 .OR. ITC.EQ.8         &
     &                                           .OR. ITC.EQ.11)THEN
                                      ABSC(IAER,IWAV)=FVOABS(IWAV)
                                      EXTC(IAER,IWAV)=FVOEXT(IWAV)
                                      ASYM(IAER,IWAV)=FVOSYM(IWAV)
                                      IF(IWAV.EQ.1)AWCCON(IAER)=FVLWC
                                  ELSEIF(ITC.EQ.12)THEN
                                      ABSC(IAER,IWAV)=DMEABS(IWAV)
                                      EXTC(IAER,IWAV)=DMEEXT(IWAV)
                                      ASYM(IAER,IWAV)=DMESYM(IWAV)
                                      IF(IWAV.EQ.1)AWCCON(IAER)=MDLWC
                                      GOTO20
                                  ELSEIF(ITC.NE.3)THEN
                                      ABSC(IAER,IWAV)=FG1ABS(IWAV)
                                      EXTC(IAER,IWAV)=FG1EXT(IWAV)
                                      ASYM(IAER,IWAV)=FG1SYM(IWAV)
                                      IF(IWAV.EQ.1)AWCCON(IAER)=AFLWC
                                  ENDIF
                              ENDIF
                              GOTO20
                          ELSEIF(ITA.LE.0)THEN
                              GOTO20
                          ENDIF
                      ENDIF

!                     NAVY MARITIME AEROSOL BECOMES MARINE IN MICROWAVE
!                     IWAV=743 for 50.00 microns.
                      IF(IWAV.GE.743 .AND. ITA.EQ.3)ITA=4

!                     RH DEPENDENT AEROSOLS
                      IF(ITA.EQ.3 .AND. IAER.EQ.1)THEN
                          IF(IAER.EQ.1)THEN
                              A2=ALOG(OCNSYM(IWAV,IRH))
                              A1=ALOG(OCNSYM(IWAV,IM1RH))
                              A=A1+(A2-A1)*(X-X1)/(X2-X1)
                              ASYM(IAER,IWAV)=EXP(A)
                              E2=ALOG(ELWCM(IRH))
                              E1=ALOG(ELWCM(IM1RH))
                              GOTO20
                          ENDIF
                      ELSEIF(ITA.EQ.5)THEN
                          Y2=ALOG(URBEXT(IWAV,IRH))
                          Y1=ALOG(URBEXT(IWAV,IM1RH))
                          Z2=ALOG(URBABS(IWAV,IRH))
                          Z1=ALOG(URBABS(IWAV,IM1RH))
                          A2=ALOG(URBSYM(IWAV,IRH))
                          A1=ALOG(URBSYM(IWAV,IM1RH))
                          E2=ALOG(ELWCU(IRH))
                          E1=ALOG(ELWCU(IM1RH))
                      ELSEIF(ITA.EQ.6)THEN
                          Y2=ALOG(TROEXT(IWAV,IRH))
                          Y1=ALOG(TROEXT(IWAV,IM1RH))
                          Z2=ALOG(TROABS(IWAV,IRH))
                          Z1=ALOG(TROABS(IWAV,IM1RH))
                          A2=ALOG(TROSYM(IWAV,IRH))
                          A1=ALOG(TROSYM(IWAV,IM1RH))
                          E2=ALOG(ELWCT(IRH))
                          E1=ALOG(ELWCT(IM1RH))
                      ELSEIF(ITA.EQ.3 .OR. ITA.EQ.4)THEN
                          Y2=ALOG(OCNEXT(IWAV,IRH))
                          Y1=ALOG(OCNEXT(IWAV,IM1RH))
                          Z2=ALOG(OCNABS(IWAV,IRH))
                          Z1=ALOG(OCNABS(IWAV,IM1RH))
                          A2=ALOG(OCNSYM(IWAV,IRH))
                          A1=ALOG(OCNSYM(IWAV,IM1RH))
                          E2=ALOG(ELWCM(IRH))
                          E1=ALOG(ELWCM(IM1RH))
                      ELSE
                          Y2=ALOG(RUREXT(IWAV,IRH))
                          Y1=ALOG(RUREXT(IWAV,IM1RH))
                          Z2=ALOG(RURABS(IWAV,IRH))
                          Z1=ALOG(RURABS(IWAV,IM1RH))
                          A2=ALOG(RURSYM(IWAV,IRH))
                          A1=ALOG(RURSYM(IWAV,IM1RH))
                          E2=ALOG(ELWCR(IRH))
                          E1=ALOG(ELWCR(IM1RH))
                      ENDIF
                      Y=Y1+(Y2-Y1)*(X-X1)/(X2-X1)
                      ZK=Z1+(Z2-Z1)*(X-X1)/(X2-X1)
                      A=A1+(A2-A1)*(X-X1)/(X2-X1)
                      ABSC(IAER,IWAV)=EXP(ZK)
                      EXTC(IAER,IWAV)=EXP(Y)
                      ASYM(IAER,IWAV)=EXP(A)
                      IF(IWAV.EQ.1)THEN
                          EC=E1+(E2-E1)*(X-X1)/(X2-X1)
                          AWCCON(IAER)=EXP(EC)
                      ENDIF
   20             CONTINUE
              ENDIF
          ENDIF
   30 CONTINUE
      DO 40 IWAV=1,NWAVLN
          IF(ICLD.EQ.18)THEN
              EXTC(5,IWAV)=CI64XT(IWAV)
              ABSC(5,IWAV)=CI64AB(IWAV)
              ASYM(5,IWAV)=CI64G(IWAV)
              AWCCON(5)=5.811E-2
          ELSEIF(ICLD.EQ.19)THEN
             EXTC(5,IWAV)=CIR4XT(IWAV)
             ABSC(5,IWAV)=CIR4AB(IWAV)
             ASYM(5,IWAV)=CIR4G(IWAV)
             AWCCON(5)=3.446E-3
          ELSE
              ABSC(5,IWAV)=0.
              EXTC(5,IWAV)=0.
              ASYM(5,IWAV)=0.
              AWCCON(5)=0.
          ENDIF
   40 CONTINUE
      RETURN
      END
