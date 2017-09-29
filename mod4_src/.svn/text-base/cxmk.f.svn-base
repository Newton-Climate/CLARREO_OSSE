      LOGICAL FUNCTION GTCXMK(NPARAM,BRDFUN,ISURF,GSURF,BRDF_PARAM)

!     GTCXMK READS IN BRDF PARAMETERS AND CALCULATES REQUIRED
!     REFLECTANCE INTEGRALS USING THE BRDF FUNCTION, BRDFUN.
!     IF SUCCESSFUL, GTBRDF RETURNS A VALUE OF TRUE.

!     INPUT ARGUMENTS:
!       NPARAM  NUMBER OF PARAMETERS IN BRDF REPRESENTATION.
!       BRDFUN  BIDIRECTIONAL REFLECTIVITY DISTRIBUTION FUNCTION.
!       ISURF   SURFACE TYPE INDEX (0=TARGET, 1=AREA-AVERAGED GROUND).
!       GSURF   LOGICAL FLAG, .TRUE. IF CURRENT SURFACE IS GROUND.
      INTEGER NPARAM,ISURF
      REAL BRDFUN
      LOGICAL GSURF
      REAL BRDF_PARAM(7,3)
      EXTERNAL BRDFUN

!     PARAMETERS:
!       MPARAM  MAXIMUM NUMBER OF PARAMETERS IN BRDF REPRESENTATION.
!       NGAUS   NUMBER OF QUADRATURE POINTS FOR INTEGRATING THE BRDF.
      INTEGER MPARAM,NGAUS
      PARAMETER(MPARAM=5,NGAUS=50)
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /SRFACE/
!       NWVSRF  NUMBER OF WAVELENGTH GRID POINTS.
!       WVSURF  WAVELENGTH GRID FOR TARGET-PIXEL (1) AND
!               AREA-AVERAGE (2) SURFACES [MICRONS].
!       SALB    AREA-AVERAGE GROUND SURFACE ALBEDO SPECTRAL ARRAY.
!       GDIRRF  AREA-AVERAGE GROUND SURFACE DIRECTIONAL REFLECTIVITY
!               AT SOLAR ANGLE.
!       HDIR    TARGET-PIXEL HEMISPHERE DIRECTIONAL REFLECTANCES
!               AT VIEWING ANGLE.
!       BRDF    TARGET-PIXEL BIDIRECTIONAL REFLECTANCE DISTRIBUTION
!               FUNCTION AT VIEWING, SUN AND RELATIVE AZIMUTH ANGLES.
!       GDIREM  AREA-AVERAGED GROUND DIRECTIONAL EMISSIVITIES
!               AT VIEWING AND GAUSSIAN QUADRATURE ANGLES.
!       GNDMOM  AREA-AVERAGED GROUND BRDF AZIMUTH FOURIER MOMENTS.
      INTEGER NWVSRF
      REAL WVSURF,SALB,GDIRRF,HDIR,BRDF,GDIREM,GNDMOM
      COMMON/SRFACE/NWVSRF(2),WVSURF(2,MWVSRF),SALB(MWVSRF),            &
     &  GDIRRF(MWVSRF),HDIR(MWVSRF),BRDF(MWVSRF),                       &
     &  GDIREM(0:MI,1:MWVSRF),GNDMOM(0:MI,0:MI,0:MAZ,1:MWVSRF)

!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI.
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

!     /ANGSRF/
!       CVWSRF  COSINE OF THE VIEW ZENITH ANGLE FROM THE SURFACE [RAD].
!       CSNSRF  COSINE OF THE SOLAR (LUNAR) ZENITH AT SURFACE [RAD].
!       AZMSRF  RELATIVE AZIMUTH ANGLE (SUN - SENSOR AT SURFACE) [RAD].
!       UMU1    COSINE OF THE PATH NADIR ANGLE.
!               (AT H1 IF IMULT=1; AT OR "NEAR" H2 IF IMULT=-1)
!       UMU0    COSINE OF THE SOLAR ZENITH ANGLE.
!               (AT H1 IF IMULT=1; AT OR "NEAR" H2 IF IMULT=-1)
!       PHI1    RELATIVE AZIMUTH ANGLE (SUN - LOS PATH AT SENSOR) [DEG].
!               (AT H1 IF IMULT=1; AT OR "NEAR" H2 IF IMULT=-1)
!       CMU     COSINE OF THE NADIR ANGLES USED IN DISORT.
      REAL CVWSRF,CSNSRF,AZMSRF,UMU1,UMU0,PHI1,CMU
      COMMON/ANGSRF/CVWSRF,CSNSRF,AZMSRF,UMU1,UMU0,PHI1,CMU(MI)

!     /DISRT/
!       DIS      LOGICAL FLAG, TRUE FOR DISORT MULTIPLE SCATTERING.
!       DISAZM   LOGICAL FLAG, TRUE FOR DISORT WITH AZIMUTH DEPENDENCE.
!       DISALB   LOGICAL FLAG, TRUE FOR DISORT SPHERICAL ALBEDO OPTION.
!       LDISCL   LOGICAL FLAG, TRUE FOR ISAACS SCALED TO DISORT.
!       NSTR     NUMBER OF DISCRETE ORDINATE STREAMS.
!       NAZ      NUMBER OF DISORT AZIMUTH COMPONENTS.
!       N2GAUS   ORDER OF DOUBLE-GAUSS QUADRATURES.
      LOGICAL DIS,DISAZM,DISALB,LDISCL
      INTEGER NSTR,NAZ,N2GAUS
      COMMON/DISRT/DIS,DISAZM,DISALB,LDISCL,NSTR,NAZ,N2GAUS

!     /JM4B2/
!       SURFZN  ZENITH ANGLE OF SURFACE NORMAL [DEG].
!       SURFAZ  SURFACE NORMAL TO VIEW DIRECTION RELATIVE
!               AZIMUTH ANGLE DEFINED AT THE SURFACE [DEG].
      REAL SURFZN,SURFAZ
      COMMON/JM4B2/SURFZN,SURFAZ

!     /JM4B3/
!       PARAMS  PARAMETERS IN BRDF REPRESENTATION.
      REAL PARAMS
      COMMON/JM4B3/PARAMS(MPARAM)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:

!     LOCAL VARIABLES:
!       IOTEST  RESULT OF IOSTAT TEST IN READ.
!       IWVSRF  CURRENT SPECTRAL WAVELENGTH INDEX.
!       IPARAM  LOOP INDEX FOR PARAMETERS IN BRDF REPRESENTATION.
!       IGAUS   LOOP INDEX FOR ANGLE COSINE QUADRATURE POINTS ON (0,1).
!       JGAUS   LOOP INDEX FOR ANGLE COSINE QUADRATURE POINTS ON (0,1).
!       KGAUS   LOOP INDEX FOR ANGLE COSINE QUADRATURE POINTS ON (0,PI).
!       I2GAUS  LOOP INDEX FOR DISORT QUADRATURE POINTS.
!       J2GAUS  LOOP INDEX FOR DISORT QUADRATURE POINTS.
!       IAZ     DISORT AZIMUTH COMPONENT LOOP INDEX.
!       J2GMAX  I2GAUS MINUS 1 (USED TO EXPLOIT BRDF RECIPROCITY).
!       DREF    DIRECTIONAL REFLECTIVITY.
!       PHISUM  BRDF INTEGRATED OVER RELATIVE AZIMUTH ANGLE.
!       GMUI    IGAUS VALUE OF GMU ARRAY.
!       GMUJ    JGAUS VALUE OF GMU ARRAY.
!       CMUI    I2GAUS VALUE OF CMU ARRAY.
!       GWT2    QUADRATURE WEIGHTS ON (0,1) DOUBLED.
      INTEGER IOTEST,IWVSRF,IPARAM,IGAUS,JGAUS,KGAUS,                   &
     &  I2GAUS,J2GAUS,IAZ,J2GMAX
      REAL DREF,PHISUM,GMUI,GMUJ,CMUI,GWT2
      LOGICAL DRF_FLAG
      LOGICAL DEBUG_FLAG

!     LOCAL ARRAYS:
!       GMU     QUADRATURE POINTS ON (0,1).
!       GWT     QUADRATURE WEIGHTS ON (0,1).
!       GMUWT   PRODUCT OF QUADRATURE POINTS AND WEIGHTS ON (0,1).
!       GPHI    QUADRATURE POINTS ON (0,PI).
!       GPHIWT  QUADRATURE-WEIGHTED FOURIER COMPONENT COSINE TERMS.
      REAL GMU(NGAUS),GWT(NGAUS),GMUWT(NGAUS),                          &
     &  GPHI(NGAUS),GPHIWT(0:MAZ,1:NGAUS)
      SAVE GMU,GWT,GMUWT,GPHI,GPHIWT

!     DATA:
!       PASS1   LOGICAL FLAG, .TRUE. FOR FIRST PASS THROUGH.
!       RECIP   LOGICAL FLAG, .TRUE. IF BRDF OBEYS RECIPROCITY LAW.
      LOGICAL PASS1,RECIP
      SAVE PASS1,RECIP
      DATA PASS1,RECIP/.TRUE.,.TRUE./

      DRF_FLAG = .TRUE.  
      DEBUG_FLAG = .TRUE.

!     CHECK NUMBER OF BRDF PARAMETERS:
      IF(NPARAM.GT.MPARAM)THEN
         IF(DEBUG_FLAG) THEN
          WRITE(IPR,'(/2A,/17X,A,I3)')' Error in GTBRDF:  Number',      &
     &      ' of parameters in BRDF representation is too large.',      &
     &      ' Please increase parameter MPARAM to',NPARAM
         ENDIF
          GTBRDF=.FALSE.
          RETURN
      ENDIF

!     READ IN NUMBER OF WAVELENGTH GRID POINTS AND SURFACE NORMAL:
      IF(LJMASS)THEN
          CALL INITCARD('CARD4B2')
      ELSE
         IF(.NOT.DRF_FLAG) THEN
          READ(IRD,*,IOSTAT=IOTEST)NWVSRF(ISURF),SURFZN,SURFAZ
         ELSE
            NWVSRF = 7.0
            SURFZN = 0.0
            SURVAZ = 0.0
         ENDIF
      ENDIF
      IF(IOTEST.NE.0)THEN
         IF(DEBUG_FLAG) THEN
          WRITE(IPR,'(/2A,I2)')' Error in GTBRDF:  Unable to',          &
     &      ' read number of wavelengths for surface',ISURF
         ENDIF
          GTBRDF=.FALSE.
          RETURN
      ELSEIF(NWVSRF(ISURF).LE.0)THEN
         IF(DEBUG_FLAG) THEN
          WRITE(IPR,'(/2A,I8)')' Error in GTBRDF:  Input number of',    &
     &      ' wavelength points for surface BRDF is',NWVSRF(ISURF)
         ENDIF
          GTBRDF=.FALSE.
          RETURN
      ELSEIF(NWVSRF(ISURF).GT.MWVSRF)THEN
         IF(DEBUG_FLAG) THEN
          WRITE(IPR,'(/2A,/17X,A,I3)')' Error in GTBRDF:  Input',       &
     &      ' number of wavelength points for surface BRDF is too',     &
     &      ' large.  Increase parameter MWVSRF to',NWVSRF(ISURF)
         ENDIF
          GTBRDF=.FALSE.
          RETURN
      ELSEIF(SURFZN.NE.0.)THEN
         IF(DEBUG_FLAG) THEN
          WRITE(IPR,'(/2A,/17X,A)')' Error in GTBRDF: ',                &
     &      ' Currently, only surface normal zenith',                   &
     &      ' angles of 0. degrees are allowed.'
         ENDIF
      ENDIF

!     DEFINE QUADRATURE POINTS
      IF(PASS1)THEN
          PASS1=.FALSE.
          CALL QGAUSN(NGAUS,GMU,GWT)
          DO IGAUS=1,NGAUS

!             COMBINE GMU AND GWT FOR THE COSINE ZENITH
!             INTEGRALS.  FOR THE AZIMUTH INTEGRALS,
!             SCALE GMU BY PI (ASSUME AZIMUTH SYMMETRY).
              GWT2=2*GWT(IGAUS)
              GMUWT(IGAUS)=GWT2*GMU(IGAUS)
              GPHI(IGAUS)=PI*GMU(IGAUS)
              IF(DIS)THEN

!                 FOR DISORT, PHI=0 CORRESPONDS TO THE TARGET PIXEL AND
!                 THE SUN BOTH BEING IN THE FORWARD HALF-PLANE.  THE
!                 BRDF IS DEFINED WITH PHI=0 CORRESPONDING TO BOTH THE
!                 SENSOR AND SUN BEING IN THE SAME HALF-PLANE.  TO
!                 ACCOUNT FOR THIS DIFFERENCE, A FACTOR OF (-1)**IAZ IS
!                 INCORPORATED INTO THE BRDF AZIMUTH FOURIER COMPONENTS:
                  GPHIWT(0,IGAUS)=GWT(IGAUS)
                  GPHIWT(1,IGAUS)=-GWT2*COS(GPHI(IGAUS))
                  DO IAZ=2,NAZ
                      GPHIWT(IAZ,IGAUS)=GWT2*COS(IAZ*GPHI(IGAUS))
                      GWT2=-GWT2
                  ENDDO
              ENDIF
          ENDDO
      ENDIF

!     BRDF TABLE HEADER:
      IF(.NOT.LJMASS)THEN
         IF(DEBUG_FLAG) THEN
          WRITE(IPR,'(/A,//A12,8(A12,I2))')' BRDF DATA TABLE:',         &
     &      ' WAVLEN (um)',('   PARAMETER',IPARAM,IPARAM=1,NPARAM)
          WRITE(IPR,'(A12,8A14)')                                       &
     &      ' -----------',('   -----------',IPARAM=1,NPARAM)
         ENDIF
      ENDIF

!     LOOP OVER WAVELENGTHS:
      DO IWVSRF=1,NWVSRF(ISURF)

!         READ IN SPECTRAL PARAMETERS:
          IF(LJMASS)THEN
              CALL INITCARD('CARD4B3')
          ELSE
             IF(.NOT.DRF_FLAG) THEN
              READ(IRD,*,IOSTAT=IOTEST)                                 &
     &          WVSURF(ISURF,IWVSRF),(PARAMS(IPARAM),IPARAM=1,NPARAM)
             ELSE
                IF(IWVSRF.EQ.1) WVSURF(ISURF,IWVSRF) = 0.47
                IF(IWVSRF.EQ.2) WVSURF(ISURF,IWVSRF) = 0.56
                IF(IWVSRF.EQ.3) WVSURF(ISURF,IWVSRF) = 0.64
                IF(IWVSRF.EQ.4) WVSURF(ISURF,IWVSRF) = 0.86
                IF(IWVSRF.EQ.5) WVSURF(ISURF,IWVSRF) = 1.24
                IF(IWVSRF.EQ.6) WVSURF(ISURF,IWVSRF) = 1.64
                IF(IWVSRF.EQ.7) WVSURF(ISURF,IWVSRF) = 2.13
                !WRITE(*,*) 'BRDF_PARAM = ',BRDF_PARAM
                PARAMS(1) = BRDF_PARAM(IWVSRF,1)
                PARAMS(2) = BRDF_PARAM(IWVSRF,2)
                PARAMS(3) = BRDF_PARAM(IWVSRF,3)
                PARAMS(4) = 0.0
                PARAMS(5) = 0.0
             ENDIF
          ENDIF
          IF(IOTEST.NE.0)THEN
             IF(DEBUG_FLAG) THEN
              WRITE(IPR,'(/2A,I2,A,/17X,A,I3)')' Error in GTBRDF: ',    &
     &          ' Unable to read surface',ISURF,' BRDF',                &
     &          ' parameters for spectral point',IWVSRF
             ENDIF
              GTBRDF=.FALSE.
              RETURN
          ENDIF
          IF(DEBUG_FLAG) THEN
          IF(.NOT.LJMASS)WRITE(IPR,'(F12.5,8F14.6)')                    &
     &      WVSURF(ISURF,IWVSRF),(PARAMS(IPARAM),IPARAM=1,NPARAM)
          ENDIF
          IF(IWVSRF.GT.1)THEN
              IF(WVSURF(ISURF,IWVSRF).LE.WVSURF(ISURF,IWVSRF-1))THEN
                 IF(DEBUG_FLAG) THEN
                  WRITE(IPR,'(/2A)')' Error in GTBRDF: ',               &
     &              ' BRDF surface wavelengths out of order.'
                 ENDIF
                  GTBRDF=.FALSE.
                  RETURN
              ENDIF
          ELSEIF(WVSURF(ISURF,1).LT.0.)THEN
             IF(DEBUG_FLAG) THEN
              WRITE(IPR,'(/2A)')' Error in GTBRDF: ',                   &
     &          ' BRDF surface wavelength less than 0.'
             ENDIF
              GTBRDF=.FALSE.
              RETURN
          ENDIF

!         SURFACE ALBEDO:
          SALB(IWVSRF)=0.
          DO IGAUS=1,NGAUS
              GMUI=GMU(IGAUS)
              DREF=0.
              DO JGAUS=1,NGAUS
                  GMUJ=GMU(JGAUS)
                  PHISUM=0.
                  DO KGAUS=1,NGAUS
                      PHISUM=PHISUM+GWT(KGAUS)                          &
     &                  *MAX(BRDFUN(GMUI,GMUJ,GPHI(KGAUS),PARAMS),0.)
                  ENDDO
                  DREF=DREF+GMUWT(JGAUS)*PHISUM
              ENDDO
              SALB(IWVSRF)=SALB(IWVSRF)+GMUWT(IGAUS)*DREF
          ENDDO

!         TARGET-PIXEL SURFACE:
          IF(ISURF.EQ.1)THEN

!             CHECK SURFACE ALBEDO:
              IF(SALB(IWVSRF).GT.1.)THEN
                 IF(DEBUG_FLAG) THEN
                  WRITE(IPR,'(/2A,/(22X,A,F12.5,A))')' Error in',       &
     &              ' GTBRDF:  The target-pixel surface albedo',        &
     &              ' at',WVSURF(1,IWVSRF),' Microns exceeds 1.'
                 ENDIF
                  GTBRDF=.FALSE.
                  RETURN
              ENDIF

!             BRDF:
              BRDF(IWVSRF)=BRDFUN(CVWSRF,CSNSRF,AZMSRF,PARAMS)
              IF(BRDF(IWVSRF).LT.0.)THEN
                 IF(DEBUG_FLAG) THEN
                  WRITE(IPR,'(/2A,F12.5,A,2(/23X,A),F12.5,A)')          &
     &              ' Warning from GTBRDF:  The target-pixel',          &
     &              ' BRDF value at',WVSURF(1,IWVSRF),' Microns',       &
     &              '(for the current sun and view angles) is negative',&
     &              '(=',BRDF(IWVSRF),                                  &
     &              ').  The value has been reset to 0.'
                 ENDIF
                  BRDF(IWVSRF)=0.
              ENDIF

!             HEMISPHERE DIRECTIONAL REFLECTIVITY:
              HDIR(IWVSRF)=0.
              DO JGAUS=1,NGAUS
                  GMUJ=GMU(JGAUS)
                  PHISUM=0.
                  DO KGAUS=1,NGAUS
                      PHISUM=PHISUM+GWT(KGAUS)                          &
     &                  *MAX(0.,BRDFUN(CVWSRF,GMUJ,GPHI(KGAUS),PARAMS))
                  ENDDO
                  HDIR(IWVSRF)=HDIR(IWVSRF)+GMUWT(JGAUS)*PHISUM
              ENDDO
              IF(HDIR(IWVSRF).GT.1.)THEN
                 IF(DEBUG_FLAG) THEN
                  WRITE(IPR,'(/2A,/22X,2(A,F12.5),A,/(A,F12.5))')       &
     &              ' Warning from GTBRDF:  The target-pixel',          &
     &              ' hemisphere directional reflectance',              &
     &              ' exceeds 1 (=',HDIR(IWVSRF),') at',                &
     &              WVSURF(1,IWVSRF),' Microns when',                   &
     &              ' the cosine of the view angle equals',CVWSRF,      &
     &              ' The value has been reset to 1.'
                 ENDIF
                  HDIR(IWVSRF)=1.
              ENDIF
          ENDIF

!         AREA-AVERAGED GROUND SURFACE:
          IF(GSURF)THEN

!             CHECK SURFACE ALBEDO:
              IF(SALB(IWVSRF).GT.1.)THEN
                 IF(DEBUG_FLAG) THEN
                  WRITE(IPR,'(/2A,/(22X,A,F12.5,A))')' Error in',       &
     &              ' GTBRDF:  The area-averaged ground surface albedo',&
     &              ' at',WVSURF(ISURF,IWVSRF),' Microns exceeds 1.'
                 ENDIF
                  GTBRDF=.FALSE.
                  RETURN
              ENDIF
              IF(DIS)THEN

!                 DIRECTIONAL EMISSIVITY:
                  GDIREM(0,IWVSRF)=1.-SALB(IWVSRF)

!                 LOOP OVER DISORT COMPUTATIONAL ANGLES:
                  DO I2GAUS=1,N2GAUS
                      GDIREM(I2GAUS,IWVSRF)=1.-SALB(IWVSRF)
                  ENDDO

!                 VIEW - SUN BRDF AZIMUTH FOURIER MOMENTS:
                  IF(UMU0.GT.0. .AND. UMU1.GT.0.)THEN
                      CALL DEFMOM(UMU1,UMU0,0,0,                        &
     &                  IWVSRF,BRDFUN,PARAMS,NGAUS,NAZ,GPHI,GPHIWT)
                  ELSE
                      DO IAZ=0,NAZ
                          GNDMOM(0,0,IAZ,IWVSRF)=0.
                      ENDDO
                  ENDIF

!                 VIEW - QUADRATURES BRDF AZIMUTH FOURIER MOMENTS:
                  IF(UMU1.GT.0.)THEN
                      DO J2GAUS=1,N2GAUS
                          CALL DEFMOM(UMU1,CMU(J2GAUS),0,J2GAUS,        &
     &                      IWVSRF,BRDFUN,PARAMS,NGAUS,NAZ,GPHI,GPHIWT)
                      ENDDO
                  ELSE
                      DO J2GAUS=1,N2GAUS
                          DO IAZ=0,NAZ
                              GNDMOM(0,J2GAUS,IAZ,IWVSRF)=0.
                          ENDDO
                      ENDDO
                  ENDIF

!                 QUADRATURES - SUN BRDF AZIMUTH FOURIER MOMENTS:
                  IF(UMU0.GT.0.)THEN
                      DO I2GAUS=1,N2GAUS
                          CALL DEFMOM(CMU(I2GAUS),UMU0,I2GAUS,0,        &
     &                      IWVSRF,BRDFUN,PARAMS,NGAUS,NAZ,GPHI,GPHIWT)
                      ENDDO
                  ELSE
                      DO I2GAUS=1,N2GAUS
                          DO IAZ=0,NAZ
                              GNDMOM(I2GAUS,0,IAZ,IWVSRF)=0.
                          ENDDO
                      ENDDO
                  ENDIF

!                 QUADRATURES - QUADRATURES BRDF AZM FOURIER MOMENTS:
                  IF(RECIP)THEN

!                     EXPLOIT BRDF RECIPROCITY:
                      CALL DEFMOM(CMU(1),CMU(1),1,1,                    &
     &                  IWVSRF,BRDFUN,PARAMS,NGAUS,NAZ,GPHI,GPHIWT)
                      J2GMAX=1
                      DO I2GAUS=2,N2GAUS
                          CMUI=CMU(I2GAUS)
                          DO J2GAUS=1,J2GMAX
                              CALL DEFMOM(CMUI,CMU(J2GAUS),             &
     &                          I2GAUS,J2GAUS,IWVSRF,BRDFUN,            &
     &                          PARAMS,NGAUS,NAZ,GPHI,GPHIWT)
                              DO IAZ=0,NAZ
                                  GNDMOM(J2GAUS,I2GAUS,IAZ,IWVSRF)      &
     &                              =GNDMOM(I2GAUS,J2GAUS,IAZ,IWVSRF)
                              ENDDO
                          ENDDO
                          CALL DEFMOM(CMUI,CMUI,I2GAUS,I2GAUS,          &
     &                      IWVSRF,BRDFUN,PARAMS,NGAUS,NAZ,GPHI,GPHIWT)
                          J2GMAX=I2GAUS
                      ENDDO
                  ELSE

!                     IGNORE BRDF RECIPROCITY:
                      DO I2GAUS=1,N2GAUS
                          CMUI=CMU(I2GAUS)
                          DO J2GAUS=1,N2GAUS
                              CALL DEFMOM(CMUI,CMU(J2GAUS),             &
     &                          I2GAUS,J2GAUS,IWVSRF,BRDFUN,            &
     &                          PARAMS,NGAUS,NAZ,GPHI,GPHIWT)
                          ENDDO
                      ENDDO
                  ENDIF
              ENDIF
              IF(.NOT.DIS .OR. LDISCL)THEN

!                 SOLAR ANGLE DIRECTIONAL REFLECTIVITY:
                  GDIRRF(IWVSRF)=0.
                  IF(UMU0.GT.0.)THEN
                      DO JGAUS=1,NGAUS
                          GMUJ=GMU(JGAUS)
                          PHISUM=0.
                          DO KGAUS=1,NGAUS
                              PHISUM=PHISUM+GWT(KGAUS)*MAX(0.,          &
     &                          BRDFUN(UMU1,GMUJ,GPHI(KGAUS),PARAMS))
                          ENDDO
                          GDIRRF(IWVSRF)                                &
     &                      =GDIRRF(IWVSRF)+GMUWT(JGAUS)*PHISUM
                      ENDDO
                      IF(GDIRRF(IWVSRF).LT.0.)THEN
                         IF(DEBUG_FLAG) THEN
                          WRITE(IPR,                                    &
     &                      '(/2A,/22X,2(A,F12.5),A,/(22X,A,F12.5))')   &
     &                      ' Warning from GTBRDF:  The area-averaged', &
     &                      ' ground directional reflectivity is',      &
     &                      ' negative (=',GDIRRF(IWVSRF),              &
     &                      ') at',WVSURF(ISURF,IWVSRF),' Microns',     &
     &                      ' when the cosine of the sun angle equals', &
     &                      UMU0,' The value has been reset to 0.'
                          ENDIF
                          GDIRRF(IWVSRF)=0.
                      ELSEIF(GDIRRF(IWVSRF).GT.1.)THEN
                         IF(DEBUG_FLAG) THEN
                          WRITE(IPR,                                    &
     &                      '(/2A,/22X,2(A,F12.5),A,/(22X,A,F12.5))')   &
     &                      ' Warning from GTBRDF:  The area-averaged', &
     &                      ' ground directional reflectivity exceeds', &
     &                      ' one (=',GDIRRF(IWVSRF),                   &
     &                      ') at',WVSURF(ISURF,IWVSRF),' Microns',     &
     &                      ' when the cosine of the sun angle equals', &
     &                      UMU0,' The value has been reset to 1.'
                          ENDIF
                          GDIRRF(IWVSRF)=1.
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF
      ENDDO
      GTBRDF=.TRUE.
      RETURN
      END
