      SUBROUTINE FDBETA(REARTH,H1ALT,H2ALT,BETA,OBSZEN,BCKZEN,          &
     &  LENN,HMIN,IERROR)

!     GIVEN H1ALT, H2ALT AND BETA (THE EARTH CENTERED ANGLE),
!     FDBETA CALCULATES THE ZENITH ANGLE AT H1ALT (OBSZEN)
!     AND AT H2ALT (BCKZEN) BASED ON A NEWTON-RAPHSON METHOD.
      IMPLICIT NONE

!     PARAMETERS:
      DOUBLE PRECISION TOLRNC
      INTEGER ITERMX
      INCLUDE 'PARAMS.h'
      PARAMETER(TOLRNC=.5D-4,ITERMX=35)

!     INPUT ARGUMENTS:
!       REARTH   RADIUS OF THE EARTH [KM].
!       H1ALT    OBSERVER (SENSOR) ALTITUDE [KM].
!       H2ALT    FINAL (TARGET) ALTITUDE [KM].
      DOUBLE PRECISION REARTH,H1ALT,H2ALT,BETA

!     OUTPUT ARGUMENTS:
!       OBSZEN   ZENITH ANGLE FROM H1ALT TO H2ALT [DEG].
!       BCKZEN   ZENITH AT TARGET (H2ALT) TO OBSERVER (H1ALT) [DEG].
!       LENN     PATH LENGTH SWITCH (0=SHORT, 1=LONG).
!       HMIN     PATH MINIMUM ALTITUDE.
      DOUBLE PRECISION OBSZEN,BCKZEN,HMIN
      INTEGER LENN,IERROR

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /GRND/
!       GNDALT   GROUND ALTITUDE, ABOVE SEA LEVEL [KM].
      DOUBLE PRECISION GNDALT
      COMMON/GRND/GNDALT

!     /DCNSTN/
!       DRIGHT   SMALLEST DOUBLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
!       DPDEG    NUMBER OF DEGREES IN ONE RADIAN IN DOUBLE PRECISION
      DOUBLE PRECISION DRIGHT,DPDEG
      COMMON/DCNSTN/DRIGHT,DPDEG

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
      INTEGER ITER
      DOUBLE PRECISION HA,HB,BETA1,HRANGE,BEND,ANGLE1,RATIOA,RATIOB,    &
     &  STORE,DENOM,DERIV,APREV,AWGHTD,DBPREV,DBETA
      IF(BETA.EQ.0.)THEN
          LENN=0
          OBSZEN=0.D0
          HMIN=H1ALT
          IF(H1ALT.GT.H2ALT)THEN
              OBSZEN=180.D0
              HMIN=H2ALT
          ENDIF
          BCKZEN=180-OBSZEN
      ENDIF
      IF(H1ALT.GT.H2ALT)THEN
          HA=H2ALT
          HB=H1ALT
      ELSE
          HA=H1ALT
          HB=H2ALT
      ENDIF

!***  GUESS AT OBSZEN,INTEGRATETO FIND BETA, TEST FOR CONVERGENCE, AND
!***  ITERATE FIRST GUESS AT OBSZEN: USE GEOMETRIC SOLN (NO REFRACTION)
      IF(.NOT.LJMASS)WRITE(IPR,'(/3(//A),/A,/)')                        &
     &  ' CASE 2D: GIVEN H1ALT, H2ALT,  BETA:',                         &
     &  ' ITERATE AROUND OBSZEN UNTIL BETA CONVERGES',                  &
     &  ' ITER     OBSZEN     BETA     DBETA'//                         &
     &  '    HRANGE      HMIN    BCKZEN  BENDING',                      &
     &  '          (DEG)     (DEG)     (DEG)'//                         &
     &  '      (KM)      (KM)     (DEG)    (DEG)'

!     CALCULATE ANGLE1, A GUESS VALUE FOR OBSZEN:
      RATIOA=(HB-HA)/(REARTH+HA)
      RATIOB=(HB-HA)/(REARTH+HB)
      STORE=2*SIN(BETA/(2*DPDEG))**2
      DENOM=RATIOB-STORE
      ANGLE1=90.D0
      IF(DENOM.NE.0.)ANGLE1=DPDEG*ATAN(SIN(BETA/DPDEG)/DENOM)
      IF(ANGLE1.LT.0.)ANGLE1=ANGLE1+180.D0

!     CALCULATE THE DERIVATIVE D(OBSZEN)/D(BETA):
      DERIV=(RATIOA+STORE)/(RATIOA*RATIOB+2*STORE)

!     DERIV TENDS TO OVERSHOOT VALUE; FUDGE=.6 SPEEDS UP CONVERGENCE:
      DERIV=.6D0*DERIV

!     BEGIN ITERATIVE PROCEDURE:
      ITER=0
 10   ITER=ITER+1
      IF(.NOT.LJMASS .AND. ITER.GT.ITERMX)THEN
          WRITE(IPR,'(/A,/(/6X,3(4X,A,F16.9),4X,A,I5,//10X,A))')        &
     &      ' FDBETA CASE 2D (H1ALT,H2ALT,BETA): Solution did not'      &
     &      //' converge.','H1ALT =',H1ALT,'H2ALT =',H2ALT,             &
     &      'BETA =',BETA,'ITERATIONS =',ITER,                          &
     &      'LAST ITERATION','OBSZEN =',ANGLE1,'BETA =',BETA1
          IERROR=1
          RETURN
      ENDIF

!     DETERMINE BETA1, THE BETA CORRESPONDING TO ANGLE1:
      CALL FINDMN(REARTH,HA,ANGLE1,HB,LENN,ITER,                        &
     &  HMIN,BCKZEN,IERROR,.TRUE.)
      CALL RFPATH(REARTH,HA,HB,ANGLE1,LENN,HMIN,.FALSE.,                &
     &  BCKZEN,BETA1,HRANGE,BEND)
      DBETA=BETA1-BETA
      IF(.NOT.LJMASS)WRITE(IPR,'(I5,3F10.4,2F10.3,2F10.4)')             &
     &  ITER,ANGLE1,BETA1,-DBETA,HRANGE,HMIN,BCKZEN,BEND

!     CHECK FOR CONVERGENCE:
      IF(ABS(DBETA).LT.TOLRNC)THEN
          IF(.NOT.LJMASS .AND. HMIN.LT.GNDALT)THEN
              WRITE(IPR,'(/A,//9X,A)')                                  &
     &          ' FDBETA, CASE 2D(H1ALT,H2ALT,BETA): REFRACTED TANGENT' &
     &          //' HEIGHT IS LESS THAN ZERO-PATH INTERSECTS THE EARTH',&
     &          ' BETA IS TOO LARGE FOR THIS H1ALT AND H2ALT'
              IERROR=1
          ELSEIF(H1ALT.LE.H2ALT)THEN
              BETA=BETA1
              OBSZEN=ANGLE1
          ELSE
              BETA=BETA1
              OBSZEN=BCKZEN
              BCKZEN=ANGLE1
          ENDIF
          RETURN
      ENDIF

      IF(ITER.GT.5)THEN
          IF(DBPREV*DBETA.LT.0.D0)THEN
              AWGHTD=(ANGLE1*ABS(DBETA)+APREV*ABS(DBPREV))              &
     &          /(ABS(DBETA)+ABS(DBPREV))
              APREV=ANGLE1
              ANGLE1=AWGHTD
              DBPREV=DBETA
              GOTO 10
          ENDIF
      ENDIF
      APREV=ANGLE1
      DBPREV=DBETA
      ANGLE1=ANGLE1-DERIV*DBPREV
      GOTO 10
      END
