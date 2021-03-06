      REAL FUNCTION SOURCE(SINIT,IVX,ISOURC,IDAY,ANGLEM,SUNFUN)

!     THIS ROUTINE RETURNS SOURCE IRRADIANCE [W CM-2 / CM-1].
!     CORRECTIONS ARE MADE FOR THE SUN'S ELLIPTIC ORBIT.  IF THE
!     SOURCE IS THE MOON RATHER THAN THE EARTH, THE PHASE ANGLE
!     BETWEEN THE SUN, MOON AND EARTH IS TAKEN IN ACCOUNT.

!     INPUT ARBUMENTS:
!       SINIT    FLAG, TRUE FOR INITIAL CALL TO SOURCE FROM TRANS.
!       IVX      FREQUENCY [CM-1].
!       ISOURC   SOURCE FLAG [0=SUN AND 1=MOON].
!       IDAY     DAY OF YEAR [0-366, DEFAULT (0) IS DAY 91].
!       ANGLEM   LUNAR PHASE ANGLE [0 TO 180 DEG].
!       SUNFUN   TOP-OF-ATMOSPHERE SOLAR IRRADIANCE
!                FUNCTION [W CM-2 / CM-1].
      LOGICAL SINIT
      INTEGER IVX,ISOURC,IDAY
      REAL ANGLEM,SUNFUN
      EXTERNAL SUNFUN

!     PARAMETERS:

      INCLUDE 'ERROR.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
      INTEGER I,IM1
      REAL V

!     SAVED VARIABLES:
      REAL FPHS,FALB,FORBIT
      SAVE FPHS,FALB,FORBIT

!     DATA:
      INTEGER NDAY(13)
      REAL RAT(13),PHS(0:16),ALB(29)
      SAVE NDAY,RAT,PHS,ALB
      DATA NDAY/1,32,60,91,121,152,182,213,244,274,305,335,366/
      DATA RAT/1.034,1.030,1.019,1.001,.985,.972,.967,.971,.982,        &
     &    .998,1.015,1.029,1.034/
      DATA PHS/100.,73.2,57.8,42.3,32.0,23.3,16.7,12.4,8.7,6.7,         &
     &    4.7,3.6,2.4,1.2,0.9,0.4,.002/
      DATA ALB/.001,.01,.03,.075,.1,.13,.155,.17,.178,.185,.2,.211,     &
     &    .231,.25,.275,.289,.285,.287,.3,.29,.3,.31,.313,.319,.329     &
     &    ,.339,.345,.350,.4/

!     CHECK FOR LUNAR SOURCE.
      IF(ISOURC.EQ.1)THEN

!         CHECK FOR INITIAL CALL.
          IF(SINIT)THEN

!             LUNAR PHASE ANGLE FACTOR
              FPHS=0.
              I=INT(ANGLEM/10)
              IF(I.LT.16)FPHS=PHS(I)+(ANGLEM/10-I)*(PHS(I+1)-PHS(I))
          ENDIF

!         GEOMETRICAL ALBEDO OF THE MOON
          FALB=.4
          IF(IVX.GT.3571)THEN
              V=10000./IVX
              I=INT(10*V)
              FALB=ALB(I)+(ALB(I+1)-ALB(I))*(10*V-I)
          ELSEIF(IVX.GT.2000)THEN
              V=10000./IVX
              FALB=ALB(28)+(ALB(29)-ALB(28))*(V-2.8)/2.2
          ENDIF
      ENDIF

!     CHECK FOR INITIAL CALL.
      IF(SINIT)THEN
          SINIT=.FALSE.

!         SUN ELLIPTIC ORBIT FACTOR
          IF(IDAY.LT.NDAY(1) .OR. IDAY.GT.NDAY(13))THEN
              FORBIT=1.
          ELSEIF(IDAY.EQ.1)THEN
              FORBIT=RAT(1)
          ELSE
              IM1=1
              DO 10 I=2,12
                  IF(IDAY.LE.NDAY(I))GOTO20
   10         IM1=I
              I=13
   20         FORBIT=RAT(IM1)+(IDAY-NDAY(IM1))*                         &
     &          (RAT(I)-RAT(IM1))/(NDAY(I)-NDAY(I-1))
          ENDIF

!         JMASS TREATS FOLLOWING MESSAGE AS WARNING:
          IF(.NOT.LJMASS)WRITE(IPR,'(/A,I4,A,F10.5)')                   &
     &      ' SUN ELLIPTIC ORBIT FACTOR FOR DAY',IDAY,':',FORBIT
      ENDIF

!     SOLAR INTENSITY [W CM-2 / CM-1]
      SOURCE=FORBIT*SUNFUN(IVX)

!     LUNAR INTENSITY [W CM-2 / CM-1]
      IF(ISOURC.EQ.1)SOURCE=2.04472E-7*FPHS*FALB*SOURCE
      RETURN
      END
      REAL FUNCTION GETS01(IVX)

!     RETURNS TOP-OF-ATMOSPHERE SOLAR IRRADIANCE VALUES TABULATED
!     EACH WAVENUMBER FROM 0 TO 50,000 CM-1.

!     INPUT ARGUMENT:
!       IVX      SPECTRAL FREQUENCY [CM-1].
      INTEGER IVX

!     PARAMETERS:
      INTEGER MAXSUN
      PARAMETER(MAXSUN=50000)

!     COMMONS:
      REAL SUN
      COMMON/SOLAR1/SUN(0:MAXSUN)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /SOLAR1/
      EXTERNAL SUNBD

!     RETURN VALUE:
      GETS01=SUN(IVX)
      RETURN
      END
      REAL FUNCTION GETS15(IVX)

!     RETURNS TOP-OF-ATMOSPHERE SOLAR IRRADIANCE VALUES TABULATED
!     EACH 15 CM-1 FROM 0 TO 49,995 CM-1.

!     INPUT ARGUMENT:
!       IVX      SPECTRAL FREQUENCY [CM-1].
      INTEGER IVX

!     PARAMETERS:
      INTEGER MSUN15
      PARAMETER(MSUN15=3333)

!     COMMONS:
      REAL SUN15
      COMMON/SOL15/SUN15(0:MSUN15)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /SOL15/
      EXTERNAL S15BD

!     RETURN VALUE:
      GETS15=SUN15(IVX/15)
      RETURN
      END
