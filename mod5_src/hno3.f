      REAL FUNCTION HNO3(V)

!     HNO3 STATISTICAL BAND PARAMETERS [CM-1/ATM].
      IMPLICIT NONE

!     INPUT ARGUMENTS:
      REAL V

!     LOCAL VARIABLES
      INTEGER I

!     DATA:
      REAL H1(15),H2(16),H3(13)
      SAVE H1,H2,H3

!     ARRAY H1 CONTAINS HNO3 ABS COEF(CM-1 / ATM) FROM  850 TO 920 CM-1
      DATA H1/2.197,3.911,6.154,8.150,9.217,9.461,11.56,11.10,          &
     &  11.17,12.40,10.49,7.509,6.136,4.899,2.866/

!     ARRAY H2 CONTAINS HNO3 ABS COEF(CM-1 / ATM) FROM 1275 TO1350 CM-1
      DATA H2/2.828,4.611,6.755,8.759,10.51,13.74,18.00,21.51,          &
     &  23.09,21.68,21.32,16.82,16.42,17.87,14.86,8.716/

!     ARRAY H3 CONTAINS HNO3 ABS COEF(CM-1 / ATM) FROM 1675 TO1735 CM-1
      DATA H3/5.003,8.803,14.12,19.83,23.31,23.58,23.22,21.09,          &
     &  26.99,25.84,24.79,17.68,9.420/
      IF(V.GE.850. .AND. V.LE.920.)THEN
          I=INT((V-845.)/5.)
          HNO3=H1(I)
      ELSEIF(V.GE.1275. .AND. V.LE.1350.)THEN
          I=INT((V-1270.)/5.)
          HNO3=H2(I)
      ELSEIF(V.GE.1675. .AND. V.LE.1735.)THEN
          I=INT((V-1670.)/5.)
          HNO3=H3(I)
      ELSE
          HNO3=0.
      ENDIF
      RETURN
      END
