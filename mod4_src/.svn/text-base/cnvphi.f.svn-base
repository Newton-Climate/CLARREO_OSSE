      SUBROUTINE CNVPHI(IERROR,ITYPE,H1,H2,ANGLE,RANGE,BETA,LENN,PHI)

!     THIS ROUTINE CONVERTS TARGET-BASED LINE-OF-SIGHT INPUTS INTO
!     MODTRAN STANDARD OBSERVER-BASED INPUTS.  MORE SPECIFICALLY, THREE
!     TARGET-BASED LINE-OF-SIGHT INPUT OPTIONS ARE NOW AVAILABLE:

!       CASE     INPUTS
!       ----     --------------------
!        2E      H2, PHI, H1 and LENN
!        2F      H2, PHI and RANGE
!        3C      H2, PHI and H1=SPACE

!     BOTH CASE 2E AND CASE 3C ARE CONVERTED TO CASE 2A (H1, H2, ANGLE
!     AND LENN), AND CASE 2F IS CONVERTED TO CASE 2C (H1, H2 AND RANGE)

!     DECLARE ROUTINE ARGUMENTS:
!       ITYPE   PATH TYPE LABEL (1=HORIZONTAL, 2=SLANT
!               BETWEEN ALTITUDES and 3=SLANT TO SPACE).
!       H1      OBSERVER (SENSOR) ALTITUDE [KM].
!       H2      TARGET (FINAL) ALTITUDE [KM].
!       ANGLE   ZENITH ANGLE AT OBSERVER TOWARDS TARGET [DEG].
!       RANGE   OBSERVER TO TARGET SLANT RANGE [KM].
!       BETA    OBSERVER TO TARGET EARTH CENTER ANGLE [DEG].
!       RE      RADIUS OF THE EARTH [KM].
!       LENN    PATH LENGTH SWITCH (0=SHORT PATH and
!               1=LONG PATH THROUGH TANGENT HEIGHT).
!       PHI     ZENITH ANGLE AT TARGER TOWARDS OBSERVER [DEG].
      INTEGER ITYPE,LENN
      REAL H1,H2,ANGLE,RANGE,BETA,PHI

!     INCLUDED FILES:
      INCLUDE 'ERROR.h'

!     LIST COMMONS:
      INCLUDE 'IFIL.h'

!     COMMON /GRAUND/:
!       GNDALT   GROUND ALTITUDE [KM].
      REAL GNDALT
      COMMON/GRAUND/GNDALT

!     COMMON /PARMTR/:
!       RE       EARTH RADIUS [KM].
!       ZMAX     TOP OF THE ATMOSPHERE [KM].
!       IPATH    SWITCH USED TO DISTINGUISH LOS FROM SUN/MOON PATHS.
      REAL RE,ZMAX
      INTEGER IPATH
      COMMON/PARMTR/RE,ZMAX,IPATH

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     DECLARE LOCAL VARIABLES:
!       IERROR   ERROR FLAG
!       LENNSV   SAVED VALUE OF INPUT LENGTH SWITCH, LENN.
!       DPH1     DOUBLE PRECISION H1 [KM].
!       DPH2     DOUBLE PRECISION H2 [KM].
!       DPANG    DOUBLE PRECISION ANGLE [DEG].
!       DPRANG   DOUBLE PRECISION RANGE [KM].
!       DPBETA   DOUBLE PRECISION BETA [DEG].
!       DPPHI    DOUBLE PRECISION PHI [DEG].
!       DPHMIN   MINIMUM PATH ALTITUDE [KM].
      INTEGER IERROR,LENNSV
      DOUBLE PRECISION DPH1,DPH2,DPANG,DPRANG,DPBETA,DPPHI,DPHMIN

!     LIST DATA:
!       ITER     ITERATION COUNTER USED BY ROUTINE DPFNMN.
!       LWARN    WARNING FLAG FOR CALL TO DPFNMN.
      INTEGER ITER
      LOGICAL LWARN
      DATA ITER/0/,LWARN/.TRUE./

!     DEFINE DOUBLE PRECISION PARAMETERS:
      DPH2=DBLE(H2)
      DPPHI=DBLE(PHI)

!     INITIALIZE ERROR FLAG
      IERROR=0

!     BRANCH BASED ON CASE.
      IF(ITYPE.EQ.3)THEN

!         CASE 3C:  H2, PHI and H1=SPACE.
          DPH1=DBLE(ZMAX)
          LENN=0
          IF(PHI.GT.90.)LENN=1
          CALL DPFNMN(DPH2,DPPHI,DPH1,LENN,                             &
     &      ITER,DPHMIN,DPANG,IERROR,LWARN)
          ITYPE=2
          H1=REAL(DPH1)
          ANGLE=REAL(DPANG)
          RANGE=0.
          IF(.NOT.LJMASS)WRITE(IPR,                                     &
     &      '(2(/A,T50,A),3(/A,F12.5,A,T50,A,F12.5,A),/T50,A,I6)')      &
     &      ' CASE 3C                    CONVERTED TO','CASE 2A',       &
     &      ' -------',                      '-------',                 &
     &      ' H2   ',H2,   ' KM',            'H1   ',H1,   ' KM',       &
     &      ' PHI  ',PHI,  ' DEG',           'ANGLE',ANGLE,' DEG',      &
     &      ' H1   ',ZMAX, ' KM',            'H2   ',H2,   ' KM',       &
     &                                       'LENN ',LENN
      ELSEIF(RANGE.GT.0.)THEN

!         CASE 2F:  H2, PHI and RANGE.
          DPRANG=DBLE(RANGE)
          CALL NEWH2(DPH2,DPH1,DPPHI,DPRANG,DPBETA,LENN,DPHMIN,DPANG)
          H1=REAL(DPH1)
          ANGLE=0.
          IF(.NOT.LJMASS)WRITE(IPR,                                     &
     &      '(2(/A,T50,A),/(A,F12.5,A,T50,A,F12.5,A))')                 &
     &      ' CASE 2F                    CONVERTED TO','CASE 2C',       &
     &      ' -------',                      '-------',                 &
     &      ' H2   ',H2,   ' KM',            'H1   ',H1,   ' KM',       &
     &      ' PHI  ',PHI,  ' DEG',           'RANGE',RANGE,' KM',       &
     &      ' RANGE',RANGE,' KM',            'H2   ',H2,   ' KM'
      ELSE

!         CASE 2E:  H2, PHI, H1 and LENN.
          ANGLE=0.
          LENNSV=LENN
          IF(H2.GE.H1 .AND. PHI.LE.90.)THEN
              WRITE(IPR,'(/2A)')' ERROR in CNVPHI: ',                   &
     &          ' UPWARD PATH BEGINS ABOVE FINAL ALTITUDE.'
              IERROR=1
          ELSEIF(H2.LE.GNDALT .AND. PHI.GT.90.)THEN
              WRITE(IPR,'(/2A)')' ERROR in CNVPHI: ',                   &
     &          ' DOWNWARD PATH BEGINS AT OR BELOW THE EARTH SURFACE.'
              IERROR=1
          ELSE
              IF(.NOT.LJMASS .AND. H1.LT.H2 .AND. PHI.GT.90.)           &
     &          WRITE(IPR,'(/2A,I6)')                                   &
     &          ' EITHER A SHORT PATH (LENN=0) OR A LONG PATH THROUGH', &
     &          ' A TANGENT HEIGHT (LENN=1) IS POSSIBLE:  LENN =',LENN
              DPH1=DBLE(H1)
              CALL DPFNMN(DPH2,DPPHI,DPH1,LENN,                         &
     &          ITER,DPHMIN,DPANG,IERROR,LWARN)
              IF(H1.NE.REAL(DPH1))THEN
                  WRITE(IPR,'(/2A)')' ERROR in CNVPHI:  SLANT PATH',    &
     &              ' INTERSECTS THE EARTH AND CANNOT REACH H1.'
                  IERROR=1
              ENDIF
              ANGLE=REAL(DPANG)
          ENDIF
          IF(.NOT.LJMASS)WRITE(IPR,                                     &
     &      '(2(/A,T50,A),3(/A,F12.5,A,T50,A,F12.5,A),/A,I6,T50,A,I6)') &
     &      ' CASE 2E                    CONVERTED TO','CASE 2A',       &
     &      ' -------',                      '-------',                 &
     &      ' H2   ',H2,   ' KM',            'H1   ',H1,   ' KM',       &
     &      ' PHI  ',PHI,  ' DEG',           'ANGLE',ANGLE,' DEG',      &
     &      ' H1   ',H1,   ' KM',            'H2   ',H2,   ' KM',       &
     &      ' LENN ',LENNSV,                 'LENN ',LENN
      ENDIF
      BETA=0.

!     RETURN TO DRIVER.
      RETURN
      END
