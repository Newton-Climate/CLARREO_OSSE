      SUBROUTINE GEOINP(SPH1,SPH2,SPANGL,SPRANG,SPBETA,                 &
     &  ITYPE,LENN,SPHMIN,SPPHI,IERROR,ISLCT,LOGLOS,SPH2MX)

!     PARAMETERS:
      INCLUDE 'ERROR.h'

!     THIS ROUTINE INTERPRETS THE ALLOWABLE COMBINATIONS OF INPUT PATH
!     PARAMETERS INTO THE STANDARD SET:  H1, H2, ANGLE, PHI, HMIN AND
!     LENN.  THE ALLOWABLE COMBINATIONS OF INPUT PARAMETERS ARE

!          FOR ITYPE = 2  (SLANT PATH H1 TO H2)
!              A.  H1, H2 AND ANGLE;
!              B.  H1, ANGLE AND RANGE;
!              C.  H1, H2 AND RANGE; OR
!              D.  H1, H2, AND BETA.

!          FOR ITYPE = 3  (SLANT PATH H1 TO SPACE/GROUND)
!              A. H1 AND ANGLE; OR
!              B. H1 AND HMIN (INPUT AS H2).

!     THIS ROUTINE ALSO DETECTS BAD INPUT (IMPOSSIBLE GEOMETRY) AND
!     ITYPE = 2 CASES WHICH INTERSECT THE EARTH, AND RETURNS THESE
!     CASES WITH ERROR FLAGS.

!     DECLARE ARGUMENTS:
!       SPH1     ALTITUDE OF OBSERVER [KM].
!       SPH2     FINAL OR TANGENT ALTITUDE [KM].
!       SPANGL   PATH ZENITH ANGLE FROM OBSERVER TOWARDS SPH2 [DEG].
!       SPRANG   SLANT PATH RANGE [KM].
!       SPBETA   SLANT PATH EARTH CENTER ANGLE [DEG].
!       ITYPE    FLAG FOR GEOMETRY TYPE.
!       LENN     LENGTH FLAG USED FOR CASE 2A WHEN H2<H1 AND ANGLE>90
!                (=0 FOR SHORT PATH, =1 FOR PATH THROUGH TANGENT POINT).
!       SPHMIN   SLANT PATH MINIMUM ALTITUDE [KM].
!       SPPHI    PATH ZENITH ANGLE FROM FINAL ALTITUDE TOWARDS
!                OBSERVER [DEG].
!       IERROR   ERROR FLAG (=0 FOR ACCEPTABLE GEOMETRY INPUTS).
!       ISLCT    ????
!       LOGLOS   LOGICAL FLAG(.TRUE. FOR ITYPE=2 WITH NO SUN).
!       SPH2MX   VALUE OF H2 PRIOR TO CALLING ROUTINE REDUCE.
      REAL SPH1,SPH2,SPANGL,SPRANG,SPBETA,SPHMIN,SPPHI,SPH2MX
      INTEGER ITYPE,LENN,IERROR,ISLCT
      LOGICAL LOGLOS

!     LIST COMMONS:
      INCLUDE 'IFIL.h'
      REAL RE,ZMAX
      INTEGER IPATH
      COMMON/PARMTR/RE,ZMAX,IPATH
      REAL GNDALT
      COMMON/GRAUND/GNDALT

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     DECLARE LOCAL VARIABLES:
      DOUBLE PRECISION H1,H2,ANGLE,RANGE,BETA,HMIN,PHI,H2ST,STORE

!     DEFINE DATA
!       ITER    COUNTER USED IN CALL TO DPFNMN.
!       LWARN   WARNING FLAG USED IN DPFNMN.
      INTEGER ITER
      LOGICAL LWARN
      DATA ITER/0/,LWARN/.TRUE./

!     DEFINE DOUBLE PRECISION GEOMETRY PARAMETERS:
      H1=DBLE(SPH1)
      H2=DBLE(SPH2)
      ANGLE=DBLE(SPANGL)
      RANGE=DBLE(SPRANG)
      BETA=DBLE(SPBETA)

!     BRANCH BASED ON PATH TYPE.
      IF(ITYPE.EQ.3)THEN

!         SLANT PATH TO SPACE.
          IF(H2.EQ.0.)THEN

!             CASE 3A:  H1 AND ANGLE.
              IF(.NOT.LJMASS .AND. NPR.LT.1)WRITE(IPR,'(//A)')          &
     &          ' CASE 3A:  GIVEN H1,H2=SPACE,ANGLE'
              H2=DBLE(ZMAX)
              IF(ANGLE.GT.90.)LENN=1
              CALL DPFNMN(H1,ANGLE,H2,LENN,ITER,HMIN,PHI,IERROR,LWARN)
          ELSE

!             CASE 3B:  H1 AND HMIN
              IF(.NOT.LJMASS .AND. NPR.LT.1)WRITE(IPR,'(//A)')          &
     &          ' CASE 3B:  GIVEN H1, HMIN, H2=SPACE'
              HMIN=H2
              H2=DBLE(ZMAX)
              IF(H1.LT.HMIN)THEN
                  WRITE(IPR,'(/2A,//10X,2(A,F13.6))')' GEOINP,',        &
     &              ' CASE 3B (H1,HMIN,SPACE):  ERROR IN INPUT DATA',   &
     &              'H1 =',H1,'    IS LESS THAN HMIN =',HMIN
                  IERROR=1
              ELSE
                  STORE=90.D0
                  CALL DPFNMN(HMIN,STORE,H1,LENN,                       &
     &              ITER,HMIN,ANGLE,IERROR,LWARN)
                  STORE=90.D0
                  CALL DPFNMN(HMIN,STORE,H2,LENN,                       &
     &              ITER,HMIN,PHI,IERROR,LWARN)
                  LENN=1
                  IF(HMIN.EQ.H1)LENN=0
              ENDIF
          ENDIF
      ELSEIF(ITYPE.EQ.2)THEN

!         SLANT PATH TO BETWEEN ALTITUDES.
          IF((LOGLOS .AND. ISLCT.EQ.21) .OR.                            &
     &      (.NOT.LOGLOS .AND. RANGE.LE.0. .AND. BETA.LE.0.))THEN

!             CASE 2A:  H1, H2, ANGLE
              IF(.NOT.LJMASS .AND. NPR.LT.1)WRITE(IPR,'(//A)')          &
     &          ' CASE 2A:  GIVEN H1, H2 AND ANGLE'
              IF(H1.GE.H2 .AND. ANGLE.LE.90.)THEN
                  WRITE(IPR,'(/2A,//(10X,2(A,F13.6),A))')' GEOINP,',    &
     &              ' CASE 2A (H1,H2,ANGLE):  ERROR IN INPUT DATA',     &
     &              'H1 (=',H1,'KM) IS GREATER THAN OR EQUAL TO H2 (= ',&
     &              H2,'KM)','AND ANGLE (=',ANGLE,                      &
     &              'DEG) DOES NOT EXCEED 90DEG.'
                  IERROR=1
              ELSEIF(H1.LE.GNDALT .AND. ANGLE.GT.90.)THEN
                  WRITE(IPR,'(/2A)')                                    &
     &              ' GEOINP, ITYPE = 2: SLANT PATH INTERSECTS',        &
     &              ' THE EARTH OR GNDALT, AND CANNOT REACH H2.'
                  IERROR=1
              ELSE
                  IF(.NOT.LJMASS .AND. H2.LT.H1 .AND. ANGLE.GT.90.      &
     &              .AND. NPR.LT.1)                                     &
     &              WRITE(IPR,'(//3A,I3)')' EITHER A SHORT PATH',       &
     &              ' (LENN=0) OR A LONG PATH THROUGH A TANGENT',       &
     &              ' HEIGHT (LENN=1) IS POSSIBLE: LENN = ',LENN
                  H2ST=H2
                  CALL DPFNMN(H1,ANGLE,H2,LENN,                         &
     &              ITER,HMIN,PHI,IERROR,LWARN)
                  IF(H2.NE.H2ST)THEN
                      WRITE(IPR,'(/2A)')                                &
     &                  ' GEOINP, ITYPE = 2: SLANT PATH INTERSECTS',    &
     &                  ' THE EARTH OR GNDALT, AND CANNOT REACH H2.'
                      IERROR=1
                  ENDIF
              ENDIF
          ELSEIF((LOGLOS .AND. ISLCT.EQ.24) .OR.                        &
     &      (.NOT.LOGLOS .AND. BETA.GT.0.))THEN

!             CASE 2D:  H1, H2, BETA
              CALL FDBETA(H1,H2,BETA,ANGLE,PHI,LENN,HMIN,IERROR)
          ELSEIF((LOGLOS .AND. ISLCT.EQ.22) .OR.                        &
     &      (.NOT.LOGLOS .AND. ANGLE.GT.0.))THEN

!             CASE 2B:  H1, ANGLE, RANGE
              CALL NEWH2(H1,H2,ANGLE,RANGE,BETA,LENN,HMIN,PHI)
              IF(ANGLE.GT.90. .AND.PHI.GT.90.)LENN=1
              CALL DPFNMN(H1,ANGLE,H2,LENN,ITER,HMIN,PHI,IERROR,LWARN)
          ELSE

!             CASE 2C:  H1, H2, RANGE
              CALL FTRANG(H1,H2,RANGE,ANGLE,PHI,LENN,HMIN,IERROR)
          ENDIF
      ELSE
          WRITE(IPR,'(/2A,I10)')' GEOINP:  ERROR IN INPUT DATA,',       &
     &      ' ITYPE NOT EQUAL TO 2 OR 3.   ITYPE =',ITYPE
          IERROR=1
      ENDIF

!     TEST IERROR AND RECHECK LENN
      SPH2MX=REAL(H2)
      IF(IERROR.EQ.0)THEN
          LENN=0
          IF(ABS(HMIN-MIN(H1,H2)).GT..00005)LENN=1

!         REDUCE PATH END POINTS ABOVE ZMAX TO ZMAX
          IF(HMIN.GE.ZMAX)THEN
              WRITE(IPR,'(/2A,//4(4X,A,F11.5))')                        &
     &          ' GEOINP:  THE ENTIRE PATH LIES ABOVE ZMAX,',           &
     &          ' THE TOP OF THE ATMOSPHERIC PROFILE',                  &
     &          'ZMAX =',ZMAX,'H1 =',H1,'H2 =',H2,'HMIN =',HMIN
              IERROR=1
          ELSE
              IF(H1.GT.ZMAX .OR. H2.GT.ZMAX)CALL REDUCE(H1,H2,ANGLE,PHI)
              IF(.NOT.LJMASS .AND. NPR.LT.1)                            &
     &          WRITE(IPR,'(//A,/5(/10X,A,F11.5,A),                     &
     &          /10X,A,I11)')' SLANT PATH PARAMETERS IN STANDARD FORM', &
     &          'H1      =',H1   ,' KM' ,'H2      =',H2 ,' KM' ,        &
     &          'ANGLE   =',ANGLE,' DEG','PHI     =',PHI,' DEG',        &
     &          'HMIN    =',HMIN ,' KM' ,'LENN    =',LENN
          ENDIF
      ENDIF

!     DEFINE SINGLE PRECISION VALUES AND RETURN
      SPH1=REAL(H1)
      SPH2=REAL(H2)
      SPANGL=REAL(ANGLE)
      SPRANG=REAL(RANGE)
      SPBETA=REAL(BETA)
      SPHMIN=REAL(HMIN)
      SPPHI=REAL(PHI)
      RETURN
      END
