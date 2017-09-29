      SUBROUTINE CD2C(LMODEL,LAPLUS,LENDAT,ISEED,MARIC1,MARK,JPRT,ICH,  &
     &  H2OCOL,O3COL)

!     PROCESS CARD2C INPUTS:

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
!       LMODEL   FLAG, .TRUE. IF MODEL ATMOSPHERE IS USED.
!       LAPLUS   LOGICAL FLAG FOR THE AEROSOL A+ OPTION.
!       LENDAT   LENGTH OF DATDIR STRING.
!       ISEED    RANDOM NUMBER SEED.
!       H2OCOL   WATER COLUMN [KM GM /M3].
!       O3COL    OZONE COLUMN [KM GM /M3].
      LOGICAL LMODEL,LAPLUS
      INTEGER LENDAT,ISEED,MARIC1,MARK,JPRT,ICH(4)
      REAL H2OCOL,O3COL

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'YPROP.h'
      INCLUDE 'YPROPC.h'

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT

!     /CARD1/
!       MODEL    ATMOSPHERIC PROFILE MODEL NUMBER.
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,MODTRN

!     /CARD1A/
      INTEGER M4,M5,M6,MDEF,IRD1,IRD2
      COMMON/CARD1A/M4,M5,M6,MDEF,IRD1,IRD2

!     /TITL/
      CHARACTER HHAZE(16)*20,HSEASN(2)*20,HVULCN(8)*20,                 &
     &  HMET(2)*20,HMODEL(0:8)*20
      COMMON/TITL/HHAZE,HSEASN,HVULCN,HMET,HMODEL

!     /JM1A1/
!       DISSTR   CHARACTER STRING USED TO READ IN DISORT LOGICALS.
!       H2OSTR   VERTICAL WATER COLUMN CHARACTER STRING (IF THE
!                FIRST NON-BLANK CHARACTER IS A "G", THE WATER
!                COLUMN IN GM/CM2 FOLLOWS "G"; IF THE FIRST
!                NON-BLANK CHARACTER IS AN "A", THE WATER COLUMN
!                IN ATM-CM AT 273.15K FOLLOWS "A"; OTHERWISE THE
!                STRING CONTAINS A WATER COLUMN SCALING FACTOR).
!       O3STR    VERTICAL OZONE COLUMN CHARACTER STRING (IF THE
!                FIRST NON-BLANK CHARACTER IS A "G", THE OZONE
!                COLUMN IN GM/CM2 FOLLOWS "G"; IF THE FIRST
!                NON-BLANK CHARACTER IS AN "A", THE OZONE COLUMN
!                IN ATM-CM AT 273.15K FOLLOWS "A"; OTHERWISE THE
!                STRING CONTAINS AN OZONE COLUMN SCALING FACTOR).
!       SUNFL2   TOP-OF-ATMOSPHERE SPECTRAL SOLAR IRRADIANCES FILE NAME.
!       BMROOT   PREFIX OF MOLECULAR BAND MODEL PARAMETERS FILE.
!       FILTNM   NAME OF FILTER RESPONSE FUNCTION FILE.
!       H2OAER   FLAG, TRUE IF DEFAULT AEROSOL PROPERTIES ARE REVISED
!                BASED ON WATER COLUMN SCALING.
!       DATDIR   NAME OF THE MODTRAN DATA DIRECTORY.
      CHARACTER DISSTR*3,H2OSTR*10,O3STR*10,SUNFL2*(NAMLEN),            &
     &  FILTNM*(NAMLEN),BMROOT*(NAMLEN),H2OAER*1,DATDIR*(NAMLEN-12)
      COMMON/JM1A1/DISSTR,H2OSTR,O3STR,SUNFL2,                          &
     &  FILTNM,BMROOT,H2OAER,DATDIR

!     /CARD2/
!       ICLD     CLOUD MODEL NUMBER.
      INTEGER IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA
      REAL VIS,WSS,WHH,RAINRT
      COMMON/CARD2/IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,                 &
     &  VIS,WSS,WHH,RAINRT

!     /JM2/
      CHARACTER APLUS*2,CNOVAM*1,ARUSS*3
      COMMON/JM2/APLUS,CNOVAM,ARUSS

!     /GRAUND/
      REAL GNDALT
      COMMON/GRAUND/GNDALT

!     /CARD3/
!       H1       OBSERVER (SENSOR) ALTITUDE [KM].
!       H2       FINAL (TARGET) ALTITUDE [KM].
!       ANGLE    ZENITH ANGLE FROM H1 TO H2 [DEG].
!       HRANGE   DISTANCE FROM H1 TO H2 [KM].
!       BETA     EARTH CENTER ANGLE BETWEEN H1 AND H2 [DEG].
!       REE      RADIUS OF THE EARTH [KM].
!       LENN     PATH LENGTH SWITCH (0=SHORT, 1=LONG).
      INTEGER LENN
      REAL H1,H2,ANGLE,HRANGE,BETA,REE
      COMMON/CARD3/H1,H2,ANGLE,HRANGE,BETA,REE,LENN

!     /COMNOV/
!       LNOVAM   LOGICAL FLAG, .TRUE. IF NOVAM AEROSOLS ARE USED.
      LOGICAL LNOVAM
      REAL EXTNOV(MNOV,MXWVLN),ABSNOV(MNOV,MXWVLN),                     &
     &  ASMNOV(MNOV,MXWVLN),WLNOV(MXWVLN)
      INTEGER NNOV,NWLNOV
      COMMON/COMNOV/LNOVAM,EXTNOV,ABSNOV,ASMNOV,WLNOV,NNOV,NWLNOV

!     /NOVDAT/
      INTEGER NLNOV
      REAL ALTNOV,RHNOV,PNOV,TNOV,DENNOV
      COMMON/NOVDAT/NLNOV,ALTNOV(MLNOV),RHNOV(MLNOV),                   &
     &  PNOV(MLNOV),TNOV(MLNOV),DENNOV(MNOV,MLNOV)

!     /CARD2A/
      INTEGER NCRALT,NCRSPC
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      &
     &  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP

!     /JM2APLUS/
      REAL ZAER11,ZAER12,SCALE1,ZAER21,ZAER22,SCALE2,                   &
     &  ZAER31,ZAER32,SCALE3,ZAER41,ZAER42,SCALE4
      COMMON/JM2APLUS/ZAER11,ZAER12,SCALE1,ZAER21,ZAER22,SCALE2,        &
     &  ZAER31,ZAER32,SCALE3,ZAER41,ZAER42,SCALE4

!     SAVED COMMONS:
      SAVE /TITL/

!     LOCAL VARIABLES:
!       ICLDV    SAVED VALUE OF ICLD.
!       RAINSV   SAVED VALUE OF RAINRT.
!       CPROB    PROBABILITY OF CLOUD OCCURRING [%].
!       LEXIST   LOGICAL FLAG, .TRUE. IF FILE EXISTS.
!       LOPEN    LOGICAL FLAG, .TRUE. IF FILE IS OPEN.
      CHARACTER YNAME(MMOLY)*10,YBM(2)*11,CDUM*8
      REAL DUM,RAINSV,CPROB
      LOGICAL LEXIST,LOPEN
      INTEGER FUNIT,IBMFIL,NBMFIL,IDUM,NUNIT,IMOLY,JMOLY,               &
     &  ICLDSV,IHMET,IHVUL

	WRITE(*,*) 'into cd2c'
!     The number of user-defined species is initialized to zero.
      NMOLY=0
      NMOLYS=0
      NMOLT=NMOLXT+NMOLY

      IF(LMODEL .OR. IM.NE.1)THEN
          IRD1=0
          IRD2=0
      ELSE

!         CARD 2C:  USER SUPPLIED ATMOSPHERIC PROFILE
          IF(LJMASS)THEN
              CALL INITCARD('CARD2C')
          ELSE
              READ(IRD,'(3I5,A20,F10.0,I5)')                            &
     &          ML,IRD1,IRD2,HMODEL(MODEL),REE,NMOLYC
              IF(REE.LE.0.)THEN
                  IF(M1.EQ.1)THEN
                      REE=6378.39
                  ELSEIF(M1.EQ.4)THEN
                      REE=6356.91
                  ELSEIF(M1.EQ.5)THEN
                      REE=6356.91
                  ELSE
                      REE=6371.23
                  ENDIF
              ENDIF
              WRITE(IPR,'(/A,3I5,A20,F10.0,I5)')' CARD 2C *****',       &
     &          ML,IRD1,IRD2,HMODEL(MODEL),REE,NMOLYC
              IF(NMOLYC.GT.0)THEN

!                 READ IN "Y-SPECIES":
                  READ(IRD,'((8A10))')(YNAME(IMOLY),IMOLY=1,NMOLYC)
                  WRITE(IPR,'((14X,8A10))')(YNAME(IMOLY),IMOLY=1,NMOLYC)
                  NMOLYS=0
                  DO IMOLY=1,NMOLYC
                      YBM(2)='           '
                      CALL YFILE(IPR,YNAME(IMOLY),YBM)
                      NBMFIL=1
                      IF(YNAME(IMOLY)(1:1).NE.'-')THEN

                          NMOLY=NMOLY+1
                          LMOLY(NMOLY)=IMOLY
                          LMOLYS(NMOLY)=0
                          CALL MVRGT(YNAME(IMOLY))
                          CNAMEY(NMOLY)=YNAME(IMOLY)(3:10)

!                         YBM(2), defined in routine YFILE, is not
!                         blank.  This means two band model parameter
!                         files per species (CENTER+TAIL).
                          IF(YBM(2).NE.'           ')THEN
                              NBMFIL=2
                              NMOLYS=NMOLYS+1
                              LMOLYS(NMOLY)=IMOLY
                          ENDIF
                          DO IBMFIL=1,NBMFIL
                              INQUIRE                                   &
     &                          (FILE=DATDIR(1:LENDAT)//YBM(IBMFIL),    &
     &                          EXIST=LEXIST,OPENED=LOPEN,NUMBER=FUNIT)
                              IF(LOPEN)CLOSE(FUNIT)
                              IF(LEXIST)THEN
                                  ITBY(NMOLY,IBMFIL)=NUNIT()
                                  OPEN(ITBY(NMOLY,IBMFIL),              &
     &                              FILE=DATDIR(1:LENDAT)//YBM(IBMFIL), &
     &                              STATUS='OLD')

!                                 For full band model files, read
!                                 molecular weight which is in both
!                                 tail and center files.  So must
!                                 read it twice but it is redundant.
                                  READ(ITBY(NMOLY,IBMFIL),*)AMWTY(NMOLY)

!                                 Compute DOPY (also done twice also):
                                  DOPY(NMOLY)                           &
     &                              =5.91861E-6/SQRT(AMWTY(NMOLY))
                              ELSE
                                  WRITE(*,'(/(A))')'Routine DRIVER'//   &
     &                              ' error:  BANDMODEL FILE ',         &
     &                              DATDIR(1:LENDAT)//YBM(IBMFIL),      &
     &                              ' does not exist.  Create it or'//  &
     &                              ' delete species with a minus',     &
     &                              ' sign in tape5.'
                                  STOP 'Band Model File does not exist.'
                              ENDIF
                          ENDDO
                      ENDIF
                  ENDDO
                  NMOLT=NMOLXT+NMOLY
              ENDIF
          ENDIF

!         REORDER SPECIES SO "NMOLYS" SPECIES ARE FIRST FOLLOWED
!         BY THE NON-STAR SPECIES.  LMOLYS(JMOLY)=0 IF JMOLY IS
!         A NON-STAR SPECIES, AND 0 < JMOLY <= NMOLY FOR STAR
!         SPECIES, LMOLYS(JMOLY) DENOTES POSITION.
          JMOLY=1
          DO IMOLY=1,NMOLY
              IF(LMOLYS(IMOLY).NE.0)THEN

                  IDUM=LMOLY(JMOLY)
                  LMOLY(JMOLY)=LMOLY(IMOLY)
                  LMOLY(IMOLY)=IDUM

!                 SWAP ITBY, CNAMEY, DOPY, AMWTY.
                  DUM=DOPY(JMOLY)
                  DOPY(JMOLY)=DOPY(IMOLY)
                  DOPY(IMOLY)=DUM

                  DUM=AMWTY(JMOLY)
                  AMWTY(JMOLY)=AMWTY(IMOLY)
                  AMWTY(IMOLY)=DUM

                  IDUM=ITBY(JMOLY,1)
                  ITBY(JMOLY,1)=ITBY(IMOLY,1)
                  ITBY(IMOLY,1)=IDUM

                  IDUM=ITBY(JMOLY,2)
                  ITBY(JMOLY,2)=ITBY(IMOLY,2)
                  ITBY(IMOLY,2)=IDUM

                  CDUM=CNAMEY(JMOLY)
                  CNAMEY(JMOLY)=CNAMEY(IMOLY)
                  CNAMEY(IMOLY)=CDUM

                  JMOLY=JMOLY+1
              ENDIF
          ENDDO
          IF(LAPLUS .AND. (IRD2.EQ.1.OR.IRD2.EQ.2))THEN
              WRITE(IPR,'(2A,/(10X,A))')' WARNING: ',                   &
     &          ' WHEN APLUS ="A+", IRD2 CANNOT BE 1 OR 2 (FOR',        &
     &          ' READING AEROSOL PROFILES WITH MODEL=0/7).',           &
     &          ' *****  APLUS OPTION WILL BE IGNORED  *****'
              LAPLUS=.FALSE.
              APLUS='  '
          ENDIF
          IF(IVSA.EQ.1)CALL RDNSM
      ENDIF
      MARIC1=0
      MARK=0
      JPRT=0

      IF(ICLD.GE.1 .AND. ICLD.LE.10)THEN

!         CLOUD/RAIN MODELS 1-10 ARE NOW SET UP IN ROUTINE CRDRIV, NOT
!         ROUTINE AERNSM; TEMPORARILY SET ICLD AND RAINRT TO ZERO.
          ICLDSV=ICLD
          RAINSV=RAINRT
          ICLD=0
          RAINRT=0.
          CALL AERNSM(JPRT,GNDALT,MARIC1,MARK,ICH,LMODEL)
          ICLD=ICLDSV
          RAINRT=RAINSV
          CALL CRDRIV
      ELSE
          CALL AERNSM(JPRT,GNDALT,MARIC1,MARK,ICH,LMODEL)
      ENDIF

!     IF NOVAM IS CALLED MERGE NOVAM LAYERS WITH OTHER.
      IF(LNOVAM)THEN
          NLNOV=2*(NNOV+1)
          CALL NOVMRG(ALTNOV,RHNOV,PNOV,TNOV,DENNOV,NLNOV)
      ENDIF

!     AERNSM IS USED AS BEFORE BUT NEW A+ SCHEME NEEDS EXTRA
!     HANDLING.  THE NEW SCHEME GOES INTO EFFECT AFTER AERNSM.
      IF(LAPLUS)CALL APRFNU(LMODEL,IHAZE,ZAER11,ZAER12,SCALE1,          &
     &  ZAER21,ZAER22,SCALE2,ZAER31,ZAER32,SCALE3,                      &
     &  ZAER41,ZAER42,SCALE4,GNDALT)
      IF(ICLD.EQ.20)THEN

!         SET UP CIRRUS MODEL
          CALL CIRRUS(CTHIK,CALT,ISEED,CPROB,CEXT)
          IF(.NOT.LJMASS)THEN
              WRITE(IPR,'(15X,A)')                                      &
     &          'CIRRUS ATTENUATION INCLUDED (N O A A CIRRUS)'
              IF(ISEED.EQ.0)THEN
                  WRITE(IPR,'((15X,2A,F10.5,A))')' CIRRUS THICKNESS ',  &
     &              'DEFAULTED TO MEAN VALUE OF',CTHIK,'KM','CIRRUS ',  &
     &              'BASE ALTITUDE DEFAULTED TO MEAN VALUE OF',CALT,'KM'
              ELSE
                  IF(CTHIK.NE.0.)THEN
                      WRITE(IPR,'(15X,2A,F10.5,A)')' CIRRUS THICKNESS', &
     &                  ' USER DETERMINED TO BE',CTHIK,'KM'

                  ELSE
                      WRITE(IPR,'(15X,2A,F10.5,A)')'CIRRUS ATTENUATION',&
     &                  ' STATISTICALLY DETERMINED TO BE',CTHIK,'KM'
                  ENDIF
                  IF(CALT.NE.0)THEN
                      WRITE(IPR,'(15X,2A,F10.5,A)')'CIRRUS BASE',       &
     &                  ' ALTITUDE USER DETERMINED TO BE',CALT,'KM'
                  ELSE
                      WRITE(IPR,'(15X,2A,F10.5,A)')'CIRRUS BASE ALTI',  &
     &                  'TUDE STATISTICALLY DETERMINED TO BE',CALT,'KM'
                  ENDIF
              ENDIF
              WRITE(IPR,'(15X,A,F7.1,A)')                               &
     &          'PROBABILITY OF CLOUD OCCURRING IS',CPROB,' PERCENT'
          ENDIF
      ENDIF

!     SCALE OZONE AND/OR WATER PROFILE:
      IF(H2OSTR.NE.'          ' .OR. O3STR.NE.'          ')             &
     &  CALL SCLCOL(H2OSTR,O3STR,H2OCOL,O3COL)

      IF(LJMASS)THEN

!         SET DEFAULT AEROSOLS:
          IF(IHAZE.GT.0 .AND. JPRT.NE.0)THEN
              IF(ISEASN.EQ.0)ISEASN=1
              IF(IVULCN.LE.0)IVULCN=1
          ENDIF
          RETURN
      ENDIF

      IF(IHAZE.LE.0 .OR. JPRT.EQ.0)RETURN
      IF(ISEASN.EQ.0)ISEASN=1
      IF(IVULCN.LE.1)THEN
          IHMET=1
          IVULCN=1
      ELSE
          IHMET=2
      ENDIF
      IF(IVULCN.EQ.6 .OR. IVULCN.EQ.7)THEN
          IHVUL=11
      ELSEIF(IVULCN.EQ.8)THEN
          IHVUL=13
      ELSE
          IHVUL=IVULCN+10
      ENDIF
      WRITE(IPR,'(/A,3(/5X,A),F10.2,A,/(4(5X,A)))')' AEROSOL MODEL',    &
     &  'REGIME                      AEROSOL TYPE             '//       &
     &     'PROFILE                  SEASON',                           &
     &  '-----------------------     --------------------     '//       &
     &     '--------------------     --------------------',             &
     &  'BOUNDARY LAYER (0-2KM)      '//HHAZE(IHAZE),                   &
     &     VIS,' KM METEOROLOGICAL_RANGE_AT_SEA_LEVEL',                 &
     &  'TROPOSPHERE  (2-10KM)  ',HHAZE(6),HHAZE(6),HSEASN(ISEASN),     &
     &  'STRATOSPHERE (10-30KM) ',HHAZE(IHVUL),HVULCN(IVULCN),          &
     &  HSEASN(ISEASN),'UPPER ATMOS (30-100KM) ',HHAZE(16),HMET(IHMET)

!     RETURN TO DRIVER:
      RETURN
      END
