      SUBROUTINE O3UV(V,O3T0)
      COMMON/O3UVF/V1,DV,NPT,S(133)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /O3UVF/
      EXTERNAL O3UVFD

!     INTERPOLATION FOR O3 CONTINUUM WITH LOWTRAN
      I=INT((V-V1)/DV+1.00001)
      IF(I.LT.1 .OR. I.GT.NPT)THEN
          O3T0=0.
      ELSE
          VR=(I-1)*DV+V1
          IF(VR.GT.V+.1 .OR. VR.LT.V-.1)THEN
              IF(I.EQ.NPT)I=NPT-1
              O3T0=(S(I+1)-S(I))*(V-VR)/DV+S(I)
          ELSE
              O3T0=S(I)
          ENDIF
      ENDIF
      RETURN
      END
