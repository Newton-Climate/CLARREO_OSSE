      SUBROUTINE BMDATA(IV,IKMX,IMSMX)

!     BMDATA (CALLED BY TRANS) MAKES THE INITIAL BAND MODEL TAPE READ
!     AND CALCULATES WAVENUMBER-INDEPENDENT PARAMETERS FOR USE BY BMOD

!     ARGUMENTS:
!       IV     SPECTRAL FREQUENCY [CM-1].
!       IKMX   CURRENT NUMBER OF LAYER BOUNDARIES.
!       IMSMX  NUMBER OF LAYER BOUNDARIES FOR MULTIPLE SCATTERING PATH.
      INTEGER IV,IKMX,IMSMX

!     PARAMETER MEXT DENOTES THE NUMBER OF MODTRAN "SPECIES".
!     THIS INCLUDES THE 12 ORIGINAL BAND MODEL PARAMETER MOLECULES
!     PLUS A HOST OF OTHER ABSORPTION AND/OR SCATTERING SOURCES.
      INCLUDE 'PARAMS.h'

!     TRANS VARIABLES
      INTEGER IBINX,IMOLX,IALFX
      REAL SDZX,ODZX
      COMMON/BMDCMX/SDZX(MXTEMP),ODZX(MXTEMP),IBINX,IMOLX,IALFX
      REAL TBBYSX
      COMMON/SOLSX/TBBYSX(LAYTHR,NMOLX)
      INCLUDE 'BASE.h'
      INCLUDE 'IFIL.h'
      INCLUDE 'SOLS.h'

!     /CARD1/
      INTEGER MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT
      LOGICAL MODTRN
      COMMON/CARD1/MODEL,ITYPE,IEMSCT,M1,M2,M3,IM,NOPRNT,MODTRN

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT
      INCLUDE 'BMHEAD.h'
      INCLUDE 'BMDAT.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
      INTEGER IT,MSMAX,IKHI,MSOFF,IK,IKOFF,J,K
      REAL FREQ,TT

!     DATA:
!       TZERO    FREEZING POINT TEMPERATURE OF WATER [KELVIN].
      REAL TZERO
      DATA TZERO/273.15/

!     REWIND THE FORMATTED BAND MODEL DATA FILE, UNIT ITBX.
      REWIND(ITBX)

!     PERFORM THE STANDARD LOWTRAN CALCULATION IF NO BAND MODEL DATA
!     EXISTS FOR THE CHOSEN FREQUENCY RANGE.
      IF(IV.GT.MXFREQ)THEN
          WRITE(IPR,'(//2A)')' (THE BAND MODEL TAPE DOES NOT',          &
     &      ' HAVE DATA IN THE REQUESTED WAVENUMBER RANGE.)'
          RETURN
      ENDIF

!     FIND 1ST RECORD (IP) IN BAND MODEL FILE WHERE FREQUENCY IV OCCURS.
      CALL GTSTRT(IV,IP,ITB,LSTREC)

!     READ THE FIRST RECORD
      READ(ITB,REC=IP)                                                  &
     &  IBIN,IMOL,(SDZ(I),I=1,MXTEMP),IALF,(ODZ(I),I=1,MXTEMP)

!     READ THE BAND MODEL PARAMETERS (ABSORPTION CROSS-SECTIONS ONLY)
!     AND POSITION AT THE PROPER PLACE OF EACH FILE OF EACH X-SPECIES.
   10 CONTINUE
      READ(ITBX,*,END=20)FREQ,IMOLX,(SDZX(IT),IT=1,NTEMP)
      IBINX=INT(FREQ+0.01)
      READ(ITBX,*)IALFX,(ODZX(IT),IT=1,NTEMP)
      IF(IV.GT.IBINX)GOTO10
   20 CONTINUE

!     SET MSMAX TO LAYTWO FOR MULTIPLE SCATTERING
      MSMAX=ABS(IMULT)*LAYTWO

!     SET TEMPERATURE INTERPOLATION INDICES FOR EACH LAYER
      IKHI=IKMX
      DO 60 MSOFF=0,MSMAX,LAYTWO
          DO 50 IK=1,IKHI
              IKOFF=IK+MSOFF
              TT=TBBY(IKOFF)
              IF(TT.LE.TBAND(1))THEN
                  JJ(IKOFF)=2
                  FF(IKOFF)=1.
              ELSEIF(TT.GE.TBAND(NTEMP))THEN
                  JJ(IKOFF)=NTEMP
                  FF(IKOFF)=0.
              ELSE
                  DO 30 J=2,NTEMP
                      IF(TT.LE.TBAND(J))GOTO40
   30             CONTINUE
   40             JJ(IKOFF)=J
                  FF(IKOFF)=(TBAND(J)-TT)/(TBAND(J)-TBAND(J-1))
              ENDIF

!         SET TEMPERATURE SCALING PARAMETERS
!           T5      LAYER TEMPERATURE DIVIDED BY 273.15K RAISED TO 0.5
!           PTM75   LAYER TEMPERATURE DIVIDED BY 273.15K RAISED TO -0.75
!                   TIMES THE LAYER PRESSURE IN ATMOSPHERES.
          T5(IKOFF)=SQRT(TT/TZERO)
          PTM75(IKOFF)=PATM(IKOFF)*(TZERO/TT)**.75
   50     CONTINUE
          IKHI=IMSMX
   60  CONTINUE

!     IF NO SINGLE SCATTER SOLAR, RETURN
      IF(IEMSCT.NE.2)RETURN

!     SET TEMPERATURE INTERPOLATION INDICES FOR SOLAR LAYERS
      IKHI=IKMX+1
      DO 110 MSOFF=0,MSMAX,LAYTWO
          DO 100 IK=1,IKHI
              IKOFF=IK+MSOFF

!             SKIP SET UP IF THE SUN IS IN THE SHADE.
              IF(WPATHS(IKOFF,36).LT.0.)GOTO100
              DO 90 K=1,NMOLXT
                  IF(K.LE.NMOL)THEN
                      TT=TBBYS(IKOFF,K)
                  ELSE
                      TT=TBBYSX(IKOFF,K-NMOL)
                  ENDIF
                  IF(TT.LE.TBAND(1))THEN
                      JJS(IKOFF,K)=2
                      FFS(IKOFF,K)=1.
                  ELSEIF(TT.GE.TBAND(NTEMP))THEN
                      JJS(IKOFF,K)=NTEMP
                      FFS(IKOFF,K)=0.
                  ELSE
                      DO 70 J=2,NTEMP
                          IF(TT.LE.TBAND(J))GOTO80
   70                 CONTINUE
   80                 JJS(IKOFF,K)=J
                      FFS(IKOFF,K)=(TBAND(J)-TT)/(TBAND(J)-TBAND(J-1))
                  ENDIF
                  IF(K.GT.NMOL)GOTO90

!                 SET TEMPERATURE SCALING PARAMETERS FOR SOLAR PATHS.
!                   T5S     SOLAR PATH TEMPERATURE DIVIDED
!                           BY 273.15K RAISED TO 0.5
!                   PTM75S  SOLAR PATH TEMPERATURE DIVIDED BY
!                           273.15K RAISED TO -0.75 TIMES SOLAR
!                           PATH PRESSURE IN ATMOSPHERES.
                  T5S(IKOFF,K)=SQRT(TT/TZERO)
                  PTM75S(IKOFF,K)=PATMS(IKOFF,K)*(TZERO/TT)**.75
   90         CONTINUE
  100     CONTINUE
          IKHI=IMSMX+1
  110 CONTINUE
      RETURN
      END