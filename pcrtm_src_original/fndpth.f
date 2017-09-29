      SUBROUTINE FNDPTH(CPATH,REARTH,ZMAX,H1ALT,HTAN,H2ALT,             &
     &  RANGEI,BETA,LENN,ANGLE,PHI)

!     FNDPTH DETERMINES H2ALT, BETA AND LENN GIVEN H1ALT, ANGLE, HRANGE,
!     HTAN, AND CPATH (THE LOWTRAN 6 MANUAL, AFGL-TR-83-0187, PP 14-17).
      IMPLICIT NONE

!     PARAMETERS:
!       DR    SHORT SLANT RANGE PARAMETER [KM].
      DOUBLE PRECISION DR
      PARAMETER(DR=.005D0)
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
!       REARTH   RADIUS OF THE EARTH [KM].
!       ZMAX     MAXIMUM ATMOSPHERIC PROFILE ALTITUDE [KM].
!       CPATH   REFRACTIVE PATH CONSTANT [KM].
!       H1ALT   OBSERVER ALTITUDE [KM].
!       HTAN    SLANT PATH TANGENT ALTITUDE [KM].
!       H2ALT   FINAL ALTITUDE [KM].
!       RANGEI  SLANT PATH RANGE [KM].
!       BETA    EACH CENTER ANGLE [DEG].
!       LENN    LENGTH SWITCH (=0 FOR SHORT PATHS,
!               =1 FOR PATHS THROUGH TANGENT POINT WITH H2ALT < H1ALT).
!       ANGLE   PATH ZENITH ANGLE AT H1ALT TOWARDS H2ALT [DEG].
!       PHI     PATH ZENITH ANGLE AT H2ALT TOWARDS H1ALT [DEG].
      INTEGER LENN
      DOUBLE PRECISION CPATH,REARTH,ZMAX,                               &
     &  H1ALT,HTAN,H2ALT,RANGEI,BETA,ANGLE,PHI

!     COMMONS:

!     /GRND/
!       GNDALT   GROUND ALTITUDE, ABOVE SEA LEVEL [KM].
      DOUBLE PRECISION GNDALT
      COMMON/GRND/GNDALT

!     /DCNSTN/
!       DRIGHT   SMALLEST DOUBLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
!       DPDEG    NUMBER OF DEGREES IN ONE RADIAN IN DOUBLE PRECISION
      DOUBLE PRECISION DRIGHT,DPDEG
      COMMON/DCNSTN/DRIGHT,DPDEG

!     LOCAL VARIABLES:
!       CAPRJ IS FOR CAPITAL R WITH SUBSCRIPT J
!       PNTGRN IS THE INTEGRAND OF EQUATION 21.
      INTEGER I,ILO
      DOUBLE PRECISION R1,R2,R,RPLDR,RANGEO,PNTGRN,STHTSV,RX,RXRAT,DX,  &
     &  CAPRJ,XJ,XJPL1,CTHETA,STHETA,DBETA,Z,DZ,BASE,PERP,DRANGE,DIFF

!     (RANGEI.LT.DR) SHOULD NOT HAPPEN; SO THIS CHECK IS REDUNDANT.
      IF(RANGEI.LT.DR)THEN
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP ' STOPPED IN FNDPTH'
      ENDIF
      RANGEO=0.D0
      BETA=0.D0
      IF(ANGLE.GT.90.D0)THEN
          ILO=1
          R1=REARTH+H1ALT
          R2=REARTH+HTAN
      ELSE
          ILO=2
      ENDIF
      DO I=ILO,2
          IF(I.EQ.2)THEN
              IF(HTAN.LT.GNDALT+.001D0 .AND. ANGLE.GT.90.D0)GOTO 20

!             IF HTAN IS NEAR 0, THEN YOU ARE ABOUT TO HIT THE EARTH:
              R2=REARTH+ZMAX
              IF(ANGLE.LE.90.D0)THEN
                  R1=REARTH+H1ALT
              ELSE
                  R1=REARTH+HTAN
              ENDIF
          ENDIF
          IF(R2.LT.R1)THEN
              DZ=-DR
          ELSE
              DZ=DR
          ENDIF
          R=R1
   10     CONTINUE
              Z=R-REARTH
              CALL IRFXN(Z,RX,RXRAT)
              STHETA=CPATH/(RX*R)
              IF(STHETA.GT. DBLE(1.))STHETA= DBLE(1.)
              IF(STHETA.LT.-DBLE(1.))STHETA=-DBLE(1.)
              STHTSV=STHETA
              CTHETA=SQRT(DBLE(1.)-STHETA**2)

!             IF(R1.GT.R2)THEN CTHETA IS NEGATIVE BECAUSE THETA.GT.9
              IF(R1.GT.R2)CTHETA=-CTHETA
              XJ=R*CTHETA
              CAPRJ=-R*RXRAT
              PNTGRN=1/(1-CAPRJ*STHETA*STHETA)
              RPLDR=R+DZ
              Z=RPLDR-REARTH
              CALL IRFXN(Z,RX,RXRAT)
              STHETA=CPATH/(RX*RPLDR)
              CTHETA=SQRT(1-STHETA**2)
              IF(R1.GT.R2)CTHETA=-CTHETA
              XJPL1=RPLDR*CTHETA
              DX=XJPL1-XJ
              DRANGE=PNTGRN*DX
              RANGEO=RANGEO+DRANGE
              DBETA=(((STHTSV+STHETA)/2)*(PNTGRN*DX))/(R-DZ/2)
              BETA=BETA+DBETA
              IF(RANGEO.GE.RANGEI)THEN
                  DIFF=(RANGEI-(RANGEO-DRANGE))
                  H2ALT=R-REARTH+(DZ/DRANGE)*DIFF
                  BETA=BETA*DPDEG
                  IF(I.EQ.2)THEN
                      LENN=1
                      IF(ANGLE.LE.90.D0)LENN=0
                      IF(H2ALT.LT.HTAN)THEN

!                         THIS WILL BE THE CASE IF I=2, AND YOU HAVE
!                         GONE THROUGH THE R-LOOP BARELY (ONLY) ONCE.
                          H2ALT=HTAN
                          LENN=0
                      ENDIF
                  ELSE
                      LENN=0
                  ENDIF

!                 CORRECTION FOR VERY SHORT PATHS, 5 KM OR LESS
                  IF(RANGEI.LT.5.D0 .AND. RANGEO/RANGEI.GT.1.05D0)THEN

!                     CALCULATE BETA BY STARIGHT LINE GEOMETRY.
                      PERP=SIN(ANGLE/DPDEG)*RANGEI
                      BASE=COS(ANGLE/DPDEG)*RANGEI+REARTH+H1ALT
                      BETA=ATAN(PERP/BASE)*DPDEG
                      RANGEO=RANGEI
                      H2ALT=BASE-REARTH
                  ENDIF
                  PHI=DBLE(180.)-ACOS(CTHETA)*DPDEG
                  RETURN
              ENDIF
              R=R+DZ
              IF((DZ.GT.0. .AND. R.LE.R2-DZ) .OR.                       &
     &           (DZ.LT.0. .AND. R.GE.R2-DZ))GOTO 10
      ENDDO

!     COMES HERE IF YOU HAVE REACHED ZMAX, BUT YOUR RANGEI IS STILL
!     NOT EQUAL TO OUTPUT VALUE.  IN THIS CASE DO THE FOLLOWING.
   20 CONTINUE
      RANGEI=RANGEO
      H2ALT=ZMAX
      IF(ANGLE.LE.90.D0)THEN
          LENN=0
      ELSE
          LENN=1
      ENDIF
      IF(HTAN.LT.GNDALT+.001D0 .AND. ANGLE.GT.90.D0)THEN

!         YOU HAVE HIT THE EARTH IF YOU ARE AT THIS POINT OF THE CODE
          LENN=0
          H2ALT=GNDALT
      ENDIF
      BETA=BETA*DPDEG
      PHI=180-ACOS(CTHETA)*DPDEG
      RETURN
      END