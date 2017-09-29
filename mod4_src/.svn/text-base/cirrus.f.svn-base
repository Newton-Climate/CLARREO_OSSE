      SUBROUTINE CIRRUS(CTHIK,CALT,ISEED,CPROB,CEXT)
!*********************************************************************
!*  ROUTINE TO GENERATE ALTITUDE PROFILES OF CIRRUS DENSITY         **
!*  PROGRAMMED BY   M.J. POST                                       **
!*                  R.A. RICHTER        NOAA/WPL                    **
!*                                      BOULDER, COLORADO           **
!*                                      01/27/1981                  **
!*                                                                  **
!*  INPUTS!                                                         **
!*           CTHIK    -  CIRRUS THICKNESS (KM)                      **
!*                       0 = USE THICKNESS STATISTICS               **
!*                       .NE. 0 = USER DEFINES THICKNESS            **
!*                                                                  **
!*           CALT     -  CIRRUS BASE ALTITUDE (KM)                  **
!*                       0 = USE CALCULATED VALUE                   **
!*                       .NE. 0 = USER DEFINES BASE ALTITUDE        **
!*                                                                  **
!*           ICIR     -  CIRRUS PRESENCE FLAG                       **
!*                       0 = NO CIRRUS                              **
!*                       .NE. 0 = USE CIRRUS PROFILE                **
!*                                                                  **
!*           MODEL    -  ATMOSPHERIC MODEL                          **
!*                       1-5  AS IN MAIN PROGRAM                    **
!*                       MODEL = 0,6,7 NOT USED SET TO 2            **
!*                                                                  **
!*           ISEED    -  RANDOM NUMBER INITIALIZATION FLAG.         **
!*                       0 = USE DEFAULT MEAN VALUES FOR CIRRUS     **
!*                       .NE. 0 = INITIAL VALUE OF SEED FOR RANF    **
!*                       FUNCTION. CHANGE SEED VALUE EACH RUN FOR   **
!*                       DIFFERENT RANDOM NUMBER SEQUENCES. THIS    **
!*                       PROVIDES FOR STATISTICAL DETERMINATION     **
!*                       OF CIRRUS BASE ALTITUDE AND THICKNESS.     **
!*                                                                  **
!*  OUTPUTS!                                                        **
!*         CTHIK        -  CIRRUS THICKNESS (KM)                    **
!*         CALT         -  CIRRUS BASE ALTITUDE (KM)                **
!*         DENSTY(16,I) -  ARRAY, ALTITUDE PROFILE OF CIRRUS DENSITY**
!*         CPROB        -  CIRRUS PROBABILITY                       **
!*                                                                  **
!*********************************************************************

!     /CARD1/
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,MODTRN
      INCLUDE 'PARAMS.h'

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT

!     /MODEL/
!       ZM       PROFILE BOUNDARY ALTITUDES [KM].
!       PM       PROFILE BOUNDARY PRESSURES [MBAR].
!       TM       PROFILE BOUNDARY TEMPERATURES [K].
!       RFNDX    PROFILE BOUNDARY REFRACTIVITIES.
!       DENSTY   PROFILE BOUNDARY DENSITIES [UNITS DEPEND ON SPECIES].
!       LRHSET   FLAG, TRUE IF RELATIVE HUMIDITY CANNOT BE SCALED.
      REAL ZM,PM,TM,RFNDX,DENSTY
      LOGICAL LRHSET
      COMMON/MODEL/ZM(LAYDIM),PM(LAYDIM),TM(LAYDIM),                    &
     &  RFNDX(LAYDIM),DENSTY(MEXT,LAYDIM),LRHSET(LAYDIM)
      DIMENSION CBASE(5,2),TSTAT(11),PTAB(5),CAMEAN(5)
      DIMENSION CBASE1(5),CBASE2(5)
      EQUIVALENCE (CBASE1(1),CBASE(1,1)),(CBASE2(1),CBASE(1,2))

      DATA  CAMEAN           / 11.0, 10.0, 8.0, 7.0, 5.0 /
      DATA  PTAB           / 0.8, 0.4, 0.5, 0.45, 0.4/
      DATA  CBASE1            / 7.5, 7.3, 4.5, 4.5, 2.5 /
      DATA  CBASE2            /16.5,13.5,14.0, 9.5,10.0 /
      DATA  TSTAT             / 0.0,.291,.509,.655,.764,.837,.892,      &
     & 0.928, 0.960, 0.982, 1.00 /

!  SET CIRRUS PROBABILITY AND PROFILE TO ALL ZEROES

      CPROB = 0.0
      MDL = MODEL

      DO 10 I=1,ML
          DENSTY(16,I)=0.
   10 CONTINUE

!  CHECK IF USER WANTS TO USE A THICKNESS VALUE HE PROVIDES, CALCULATE
!  A STATISTICAL THICKNESS, OR USE A MEAN THICKNESS (ISEED = 0).
!  DEFAULTED MEAN CIRRUS THICKNESS IS 1.0 KM.

      IF ( CTHIK .GT. 0.0 ) GO TO 25
      IF ( ISEED .NE. 0 ) GO TO 15
      CTHIK = 1.0
      GO TO 25

!  CALCULATE CLOUD THICKNESS USING LOWTRAN CIRRUS THICKNESS STATISTICS
!  NOTE - THIS ROUTINE USES A UNIFORM RANDOM NUMBER GENERATOR
!  FUNCTION (RANF) WHICH RETURNS A NUMBER BETWEEN 0 AND 1.
!  THIS FEATURE IS MACHINE DEPENDENT

   15 CONTINUE
      URN = RANDOM()
      DO 20 I = 1, 10
         IF (URN .GE. TSTAT(I) .AND. URN .LT. TSTAT(I+1)) CTHIK = I-1
   20 CONTINUE
      CTHIK = CTHIK / 2.0  +  RANDOM() / 2.0

!  DENCIR IS CIRRUS DENSITY IN KM-1

25    IF(CEXT .GT. 0.) THEN
           DENCIR = CEXT / 2.
      ELSE
           DENCIR = 0.07 * CTHIK
      ENDIF

!  BASE HEIGHT CALCULATIONS

      IF ( MODEL .LT. 1  .OR.  MODEL .GT. 5 ) MDL = 2
      CPROB = 100.0 * PTAB(MDL)

      HMAX = CBASE(MDL,2) - CTHIK
      BRANGE = HMAX - CBASE(MDL,1)
      IF ( CALT .GT. 0.0 ) GO TO 27
      IF ( ISEED .NE. 0 ) GO TO 26
      CALT = CAMEAN(MDL)
      GO TO 27
   26 CALT = BRANGE * RANDOM()+ CBASE(MDL,1)

!  PUT CIRRUS DENSITY IN CORRECT ALTITUDE BINS. IF MODEL = 7,
!  INTERPOLATE EH(16,I) FOR NON-STANDARD ALTITUDE BOUNDARIES.

   27 IF(MODEL .EQ. 7) GO TO 60
      IV1=INT(CALT )
      IV2=INT(CALT+CTHIK )
      DO 30 I = 2, 16
         IF(I .GE. IV1 .AND. I .LE. IV2) DENSTY(16,I+1) =  DENCIR
   30 CONTINUE

!  ADJUST FIRST AND LAST CIRRUS LEVEL IF CLOUD DOES NOT ENTIRELY
!  FILL EACH LEVEL.

      IHGT1 = INT( CALT )
      IHGT2 = INT( CALT + CTHIK)
      IF( IHGT1 . NE . IHGT2 ) GO TO 35
      DENSTY(16,IHGT1+1) = DENSTY( 16,IHGT1+1)*CTHIK
      RETURN
   35 PCT1  = 1.0 - ( CALT - IHGT1 )
      DENSTY(16,IHGT1+1) = DENSTY(16,IHGT1+1) * PCT1
      PCT2 =  ( CALT + CTHIK) - IHGT2
      DENSTY(16,IHGT2+1) = DENSTY(16,IHGT2+1) * PCT2
      RETURN

!  INTERPOLATE DENSTY(16,I) FOR USER SUPPLIED ALTITUDE BOUNDARIES

   60 TOP = CALT + CTHIK
      BOTTOM = CALT
      IF(TOP.LT.ZM(1) .OR. BOTTOM.GT.ZM(ML))RETURN
      IML = ML - 1
      DO 70 I=1,IML
         ZMIN=ZM(I)
         ZMAX=ZM(I+1)
         DENOM = ZMAX - ZMIN
         IF(BOTTOM.LE.ZMIN .AND. TOP.GE.ZMAX) DENSTY(16,I) = DENCIR
         IF(BOTTOM.GE.ZMIN .AND. TOP.LT.ZMAX)                           &
     &        DENSTY(16,I) = DENCIR * CTHIK/DENOM
         IF(BOTTOM.GE.ZMIN .AND. TOP.GE.ZMAX .AND. BOTTOM.LT.ZMAX)      &
     &        DENSTY(16,I) = DENCIR * (ZMAX - BOTTOM)/ DENOM
         IF(BOTTOM.LT.ZMIN .AND. TOP.LE.ZMAX .AND.TOP.GT.ZMIN)          &
     &        DENSTY(16,I) = DENCIR * (TOP - ZMIN) / DENOM
   70 CONTINUE
      RETURN
      END
