      SUBROUTINE CRUSPC(REFWAV,REFDEP,                                  &
     &      ICLD0,ICIR0,XNRMWD,XNRMIP,CWDCOL,CIPCOL,CTHICK,WVL_CLD,     &
     &      EXT_LIQ_CLD,SSA_LIQ_CLD,ASYM_LIQ_CLD,                       &
     &      EXT_ICE_CLD,SSA_ICE_CLD,ASYM_ICE_CLD,NUM_LEVS,NUM_AC_WVL)

!     THIS ROUTINE DEFINES USER-DEFINED CLOUD SPECTRAL DATA.
!     EXTC(6,IWAV)   WATER DROPLET EXTINCTION COEFS [KM-1 M3/GM]
!     ABSC(6,IWAV)   WATER DROPLET ABSORPTION COEFS [KM-1 M3/GM]
!     ASYM(6,IWAV)   WATER DROPLET HENYEY-GREENSTEIN ASYMMETRY FACTORS
!     EXTC(7,IWAV)   ICE PARTICLE EXTINCTION COEFS [KM-1 M3/GM]
!     ABSC(7,IWAV)   ICE PARTICLE ABSORPTION COEFS [KM-1 M3/GM]
!     ASYM(7,IWAV)   ICE PARTICLE HENYEY-GREENSTEIN ASYMMETRY FACTORS

!     DECLARE INPUTS
!     REFWAV   REFERENCE WAVELENGTH [MICRONS]
!     REFDEP   REFERENCE VERTICAL OPTICAL DEPTH AT REFWAV
!     ICLD0    INDEX OF CLOUD WATER DROPLET MODEL
!     ICIR0    INDEX OF CLOUD ICE PARTICLE MODEL
!     XNRMWD   DEFAULT WATER DROPLET MODEL EXTINCTION
!              COEFFICIENT AT 0.55 MICRONS [KM-1 M3/GM]
!     XNRMIP   DEFAULT ICE PARTICLE MODEL EXTINCTION
!              COEFFICIENT AT 0.55 MICRONS [KM-1 M3/GM]
!     CWDCOL   CLOUD WATER DROPLET VERTICAL COLUMN DENSITY [KM GM/M3]
!     CIPCOL   CLOUD ICE PARTICLE VERTICAL COLUMN DENSITY [KM GM/M3]
!     CTHICK   CLOUD VERTICAL THICKNESS [KM]
      INTEGER ICLD0,ICIR0
      REAL REFWAV,REFDEP,XNRMWD,XNRMIP,CWDCOL,CIPCOL,CTHICK
      REAL WVL_CLD(24),EXT_LIQ_CLD(24)
      REAL SSA_LIQ_CLD(24),ASYM_LIQ_CLD(24)
      REAL EXT_ICE_CLD(24),SSA_ICE_CLD(24),ASYM_ICE_CLD(24)

!     PARAMETERS:
      INTEGER NCLDS,NCIRS
      PARAMETER(NCLDS=5,NCIRS=2)
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:
      INCLUDE 'BASE.h'
      INCLUDE 'IFIL.h'
      REAL VX0
      COMMON/EXTWAV/VX0(NWAVLN)
      REAL CLDSPC
      COMMON/CLDDAT/CLDSPC(NWAVLN,3,NCLDS)
      REAL CIRSPC
      COMMON/CIRR/CIRSPC(NWAVLN,3,NCIRS)
      INTEGER NCRALT,NCRSPC
      REAL CTHIK,CALT,CEXT,CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP
      COMMON/CARD2A/CTHIK,CALT,CEXT,NCRALT,NCRSPC,                      &
     &  CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /CLDDAT/,/CIRR/,/EXTWAV/
      EXTERNAL DEVCBD,CLDDTA,EXTDTA

!     LIST LOCAL VARIABLES AND ARRAYS
      INTEGER NWAVM1,I0LO,I0HI,I,N,IWAV,IWAVM1
      REAL WAVOLD,FACTOR,X,EXTWD,EXTIP,RATDEP

!     CHECK THAT NCRSPC IS NOT TOO LARGE
      IF(NCRSPC.GT.MXWVLN)THEN
          WRITE(IPR,'(/2A,I3,/14X,A,I3,A)')' FATAL ERROR:  THE INPUT',  &
     &      ' NUMBER OF CLOUD SPECTRAL DATA POINTS (NCRSPC =',NCRSPC,   &
     &      ' EXCEEDS THE MAXIMUM NUMBER (PARAMETER MXWVLN =',MXWVLN,   &
     &      ' IN PARAMETER.LST).'
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'INPUT NUMBER OF CLOUD SPECTRAL DATA POINTS IS TOO LARGE'
      ENDIF

!     ANNOUNCE READING OF CLOUD SPECTRAL DATA:
      IF(.NOT.LJMASS)WRITE(IPR,'(/2A,I3,A)')' USER-DEFINED CLOUD',      &
     &  ' SPECTRAL DATA AT',NCRSPC,' SPECTRAL POINTS.'

!     LOOP OVER SPECTRAL DATA:
      NWAVM1=NWAVLN-1
      I0LO=1
      WAVOLD=-.1
      DO I=1,NCRSPC

!         READ SPECTRAL DATA:
          IF(LJMASS) THEN
              CALL INITCARD( 'CARD2E2' )
          ELSE
!DRF              READ(IRD,'(7F10.0)',ERR=30)WAVLEN(I),EXTC(6,I),ABSC(6,I), &
!DRF     &          ASYM(6,I),EXTC(7,I),ABSC(7,I),ASYM(7,I)
             WAVLEN(I) = WVL_CLD(I)
             EXTC(6,I) = EXT_LIQ_CLD(I)
             ABSC(6,I) = SSA_LIQ_CLD(I)
             ASYM(6,I) = ASYM_LIQ_CLD(I)
             EXTC(7,I) = EXT_ICE_CLD(I)
             ABSC(7,I) = SSA_ICE_CLD(I)
             ASYM(7,I) = ASYM_ICE_CLD(I)
          ENDIF

!         CHECK WAVELENGTH:
          IF(WAVLEN(I).LE.WAVOLD)THEN
              WRITE(IPR,'(/3A,/14X,2A,/(14X,I6,F10.6))')' FATAL',       &
     &          ' ERROR:  CLOUD SPECTRAL DATA MUST BE READ IN',         &
     &          ' INCREASING WAVELENGTH ORDER',' THE WAVELENGTHS',      &
     &          ' ENCOUNTERED THUS FAR ARE:',(N,WAVLEN(N),N=1,I)
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP 'NON-INCREASING WAVELENGTH IN CLOUD SPECTRAL TABLE'
          ENDIF
          WAVOLD=WAVLEN(I)

!         SET UP INTERPOLATION
          IF(WAVLEN(I).LT.VX0(1))THEN
              I0HI=1
              FACTOR=0.
          ELSE
              DO I0HI=I0LO+1,NWAVM1
                  IF(WAVLEN(I).LE.VX0(I0HI))THEN
                      FACTOR=(WAVLEN(I)-VX0(I0LO))/(VX0(I0HI)-VX0(I0LO))
                      GOTO 10
                  ENDIF
                  I0LO=I0HI
              ENDDO
              I0HI=NWAVM1
              FACTOR=0.
          ENDIF
   10     CONTINUE

!         CHECK WATER DROPLET SPECTRAL DATA

!         IF THE INPUT SPECTRAL EXTINCTION COEFFICIENT IS NEGATIVE, IT
!         IS REPLACED BY THE WAVELENGTH-INTERPOLATED CLOUD MODEL VALUE.
!         IF THE INPUT SPECTRAL ABSORPTION COEFFICIENT IS LESS THAN -1.
!         OR IF IT EXCEEDS THE EXTINCTION COEFFICIENT, IT IS REPLACED
!         BY THE VALUE WHICH YIELDS THE CLOUD MODEL SINGLE SCATTERING
!         ALBEDO.  IF THE INPUT SPECTRAL ABSORPTION COEFFICIENT IS
!         NEGATIVE BUT NOT LESS THAN -1., THAN THE INPUT VALUE IS TAKEN
!         TO BE THE SINGLE SCATTERING ALBEDO MINUS ONE, -(1-OMEGA).
          IF(EXTC(6,I).LT.0.)THEN
              EXTC(6,I)=XNRMWD*(CLDSPC(I0LO,1,ICLD0)                    &
     &          +FACTOR*(CLDSPC(I0HI,1,ICLD0)-CLDSPC(I0LO,1,ICLD0)))
              IF(ABSC(6,I).LT.-1. .OR. ABSC(6,I).GT.EXTC(6,I))THEN
                  ABSC(6,I)=XNRMWD*(CLDSPC(I0LO,2,ICLD0)                &
     &              +FACTOR*(CLDSPC(I0HI,2,ICLD0)-CLDSPC(I0LO,2,ICLD0)))
              ELSEIF(ABSC(6,I).LT.0.)THEN
                  ABSC(6,I)=-ABSC(6,I)*EXTC(6,I)
              ENDIF
          ELSEIF(ABSC(6,I).LT.-1. .OR. ABSC(6,I).GT.EXTC(6,I))THEN
              X=CLDSPC(I0LO,1,ICLD0)                                    &
     &          +FACTOR*(CLDSPC(I0HI,1,ICLD0)-CLDSPC(I0LO,1,ICLD0))
              ABSC(6,I)=0.
              IF(X.GT.0.)ABSC(6,I)=EXTC(6,I)*(CLDSPC(I0LO,2,ICLD0)      &
     &          +FACTOR*(CLDSPC(I0HI,2,ICLD0)-CLDSPC(I0LO,2,ICLD0)))/X
          ELSEIF(ABSC(6,I).LT.0.)THEN
              ABSC(6,I)=-ABSC(6,I)*EXTC(6,I)
          ENDIF

!         IF CARD2A INPUT IS AN ACCEPTABLE VALUE (-1 < ASYMWD < 1) FOR
!         THE WATER DROPLET HENYEY-GREENSTEIN PHASE FUNCTION ASYMMETRY
!         FACTOR, THE CONSTANT VALUE IS USED AT ALL WAVELENGTHS.
!         IF THE ASYMWD IS OUT OF RANGE, THE USER-DEFINED SPECTRAL
!         ASYMMETRY FACTOR IS USED UNLESS IT IS ALSO OUT OF RANGE.
!         IF BOTH ASYMMETRY FACTOR INPUTS ARE OUT OF RANGE, THE
!         WAVELENGTH-INTERPOLATED CLOUD MODEL VALUE IS USED.
          IF(ABS(ASYMWD).LT.1.)THEN
              ASYM(6,I)=ASYMWD
          ELSEIF(ABS(ASYM(6,I)).GE.1.)THEN
              ASYM(6,I)=CLDSPC(I0LO,3,ICLD0)                            &
     &          +FACTOR*(CLDSPC(I0HI,3,ICLD0)-CLDSPC(I0LO,3,ICLD0))
          ENDIF

!         CHECK ICE PARTICLE SPECTRAL DATA

!         IF THE INPUT SPECTRAL EXTINCTION COEFFICIENT IS NEGATIVE, IT
!         IS REPLACED BY THE WAVELENGTH-INTERPOLATED CIRRUS MODEL VALUE.
!         IF THE INPUT SPECTRAL ABSORPTION COEFFICIENT IS NEGATIVE OR
!         IF IT EXCEEDS THE EXTINCTION COEFFICIENT, IT IS REPLACED BY
!         THE VALUE WHICH YIELDS THE CIRRUS MODEL SINGLE SCATTERING
!         ALBEDO.  IF THE INPUT SPECTRAL ABSORPTION COEFFICIENT IS
!         NEGATIVE BUT NOT LESS THAN -1., THAN THE INPUT VALUE IS TAKEN
!         TO BE THE SINGLE SCATTERING ALBEDO MINUS ONE, -(1-OMEGA).
          IF(EXTC(7,I).LT.0.)THEN
              EXTC(7,I)=XNRMIP*(CIRSPC(I0LO,1,ICIR0)                    &
     &          +FACTOR*(CIRSPC(I0HI,1,ICIR0)-CIRSPC(I0LO,1,ICIR0)))
              IF(ABSC(7,I).LT.-1. .OR. ABSC(7,I).GT.EXTC(7,I))THEN
                  ABSC(7,I)=XNRMIP*(CIRSPC(I0LO,2,ICIR0)                &
     &              +FACTOR*(CIRSPC(I0HI,2,ICIR0)-CIRSPC(I0LO,2,ICIR0)))
              ELSEIF(ABSC(7,I).LT.0.)THEN
                  ABSC(7,I)=-ABSC(7,I)*EXTC(7,I)
              ENDIF
          ELSEIF(ABSC(7,I).LT.-1. .OR. ABSC(7,I).GT.EXTC(7,I))THEN
              X=CIRSPC(I0LO,1,ICIR0)                                    &
     &          +FACTOR*(CIRSPC(I0HI,1,ICIR0)-CIRSPC(I0LO,1,ICIR0))
              ABSC(7,I)=0.
              IF(X.GT.0.)ABSC(7,I)=EXTC(7,I)*(CIRSPC(I0LO,2,ICIR0)      &
     &          +FACTOR*(CIRSPC(I0HI,2,ICIR0)-CIRSPC(I0LO,2,ICIR0)))/X
          ELSEIF(ABSC(7,I).LT.0.)THEN
              ABSC(7,I)=-ABSC(7,I)*EXTC(7,I)
          ENDIF

!         IF CARD2A INPUT IS AN ACCEPTABLE VALUE (-1 < ASYMIP < 1) FOR
!         THE ICE PARTICLE HENYEY-GREENSTEIN PHASE FUNCTION ASYMMETRY
!         FACTOR, THE CONSTANT VALUE IS USED AT ALL WAVELENGTHS.
!         IF THE ASYMIP IS OUT OF RANGE, THE USER-DEFINED SPECTRAL
!         ASYMMETRY FACTOR IS USED UNLESS IT IS ALSO OUT OF RANGE.
!         IF BOTH ASYMMETRY FACTOR INPUTS ARE OUT OF RANGE, THE
!         WAVELENGTH-INTERPOLATED CIRRUS MODEL VALUE IS USED.
          IF(ABS(ASYMIP).LT.1.)THEN
              ASYM(7,I)=ASYMIP
          ELSEIF(ABS(ASYM(7,I)).GE.1.)THEN
              ASYM(7,I)=CIRSPC(I0LO,3,ICIR0)                            &
     &          +FACTOR*(CIRSPC(I0HI,3,ICIR0)-CIRSPC(I0LO,3,ICIR0))
          ENDIF
      ENDDO

!     RENORMALIZE SPECTRAL DATA IF REFERENCE OPTICAL DEPTH WAS INPUT
      IF(REFDEP.GT.0.)THEN

!         DETERMINE BRACKETING WAVELENGTHS FOR MODEL DATA
          IF(REFWAV.LT.WAVLEN(1) .OR. REFWAV.GT.WAVLEN(NCRSPC))THEN
              WRITE(IPR,'(/A,F10.6,A,/14X,2A,/14X,A,2(F10.6,A))')       &
     &          ' FATAL ERROR:  THE WAVELENGTH (',REFWAV,               &
     &          ' MICRONS) USED TO DEFINE THE CLOUD',                   &
     &          ' VERTICAL OPTICAL DEPTH IS OUTSIDE THE RANGE',         &
     &          ' OF THE USER-DEFINED',' CLOUD SPECTRAL DATA (',        &
     &          WAVLEN(1),' TO',WAVLEN(NCRSPC),' MICRONS).'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP 'INPUT SPECTRAL CLOUD DEPTH OUTSIDE SPECTRAL RANGE'
          ENDIF
          IWAVM1=1
          DO IWAV=2,NCRSPC-1
              IF(REFWAV.LE.WAVLEN(IWAV))GOTO 20
              IWAVM1=IWAV
          ENDDO
          IWAV=NCRSPC
   20     CONTINUE

!         DETERMINE DEFAULT EXTINCTION COEFFICIENTS AT REFWAV
          FACTOR=(REFWAV-WAVLEN(IWAVM1))/(WAVLEN(IWAV)-WAVLEN(IWAVM1))
          EXTWD=EXTC(6,IWAVM1)+FACTOR*(EXTC(6,IWAV)-EXTC(6,IWAVM1))
          EXTIP=EXTC(7,IWAVM1)+FACTOR*(EXTC(7,IWAV)-EXTC(7,IWAVM1))

!         DETERMINE RATIO OF INPUT TO CURRENT CLOUD DEPTH
          RATDEP=REFDEP/(EXTWD*CWDCOL+EXTIP*CIPCOL)

!         SCALE THE CLOUD PARTICLE SPECTRAL DATA
          DO IWAV=1,NCRSPC
              EXTC(6,IWAV)=RATDEP*EXTC(6,IWAV)
              ABSC(6,IWAV)=RATDEP*ABSC(6,IWAV)
              EXTC(7,IWAV)=RATDEP*EXTC(7,IWAV)
              ABSC(7,IWAV)=RATDEP*ABSC(7,IWAV)
          ENDDO
      ENDIF

!     RETURN TO CRUSPC IF SPECTRAL DATA IS NOT TO BE OUTPUT.
      !DRF IF(NPR.GE.0 .OR. LJMASS)RETURN

!     WRITE SPECTRAL DATA HEADER
      WRITE(IPR,'(A,/A,//54X,A,34X,A,/38X,A,2X,A,/(3A))')'1',           &
     &  ' CLOUD SPECTRAL DATA','WATER DROPLETS','ICE PARTICLES',        &
     &  '----------------------------------------------',               &
     &  '----------------------------------------------',               &
     &  ' IWAV   WAVLEN       FREQ   VERT EXT',                         &
     &  '  EXT COEF  ABS COEF  SCT COEF     ASYM  SCT ALB',             &
     &  '  EXT COEF  ABS COEF  SCT COEF     ASYM  SCT ALB',             &
     &  '      (MICRON)     (CM-1)     (KM-1)',                         &
     &  '  (        KM-1 M3/GM        )                  ',             &
     &  '  (        KM-1 M3/GM        )'

!     WRITE SPECTRAL DATA
      WRITE(IPR,'((I4,F10.4,F11.3,F11.5,2(3F10.5,2F9.5)))')             &
     &  (IWAV,WAVLEN(IWAV),10000./WAVLEN(IWAV),                         &
     &  (EXTC(6,IWAV)*CWDCOL+EXTC(7,IWAV)*CIPCOL)/CTHICK,               &
     &  EXTC(6,IWAV),ABSC(6,IWAV),EXTC(6,IWAV)-ABSC(6,IWAV),            &
     &  ASYM(6,IWAV),1-ABSC(6,IWAV)/(EXTC(6,IWAV)+1.E-20),              &
     &  EXTC(7,IWAV),ABSC(7,IWAV),EXTC(7,IWAV)-ABSC(7,IWAV),            &
     &  ASYM(7,IWAV),1-ABSC(7,IWAV)/(EXTC(7,IWAV)+1.E-20),IWAV=1,NCRSPC)
      WRITE(IPR,'(/A,//)')' END OF CLOUD PARTICLE SPECTRAL DATA'

!     RETURN TO ROUTINE CRSPEC
      RETURN

!     FATAL ERROR READING CLOUD SPECTRAL DATA
   30 CONTINUE
      WRITE(IPR,'(/A,I3,A)')' FATAL ERROR:  UNABLE TO READ LINE',       &
     &  I,' OF CLOUD SPECTRAL DATA.'
      IF(I.GT.1)WRITE(IPR,'(/A,I3,A,/(5X,7F10.6))')                     &
     &  ' THE FIRST',I-1,' LINES OF DATA ARE:',(WAVLEN(N),EXTC(6,N),    &
     &  ABSC(6,N),ASYM(6,N),EXTC(7,N),ABSC(7,N),ASYM(7,N),N=1,I-1)
      IF(LJMASS)CALL WRTBUF(FATAL)
      STOP 'PROBLEM READING CLOUD SPECTRAL DATA'
      END
