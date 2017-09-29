      FUNCTION   GAMFOG(MR,FREQ,T,RHO)

!        COMPUTES ATTENUATION OF EQUIVALENT LIQUID WATER CONTENT
!       IN CLOUDS OR FOG IN DB/KM
!       CONVERTED TO NEPERS BY NEW CONSTANT 1.885

!        FREQ = WAVENUMBER (INVERSE CM)
!        T    = TEMPERATURE (DEGREES KELVIN)
!        RHO  = EQUIVALENT LIQUID CONTENT  (G/CUBIC METER)
!      CINDEX=COMPLEX DIELECTRIC CONSTANT M  FROM INDEX
!      WAVL = WAVELENGTH IN CM

      COMPLEX CINDEX
      IF(RHO.GT.0.) GO TO 2
      GAMFOG=0.
      RETURN
    2 CONTINUE
      KEY=1
      IF(MR. GE. 5) KEY = 0
      WAVL=1.0/FREQ
      TC=T-273.2
!CC
!CC    CHANGE TEMP SO THAT MINIMUM IS -20.0 CENT.
!CC
      IF(TC.LT.-20.0) TC=-20.0
      CALL INDX (WAVL,TC,KEY,REIL,AIMAK)
      CINDEX=CMPLX(REIL,AIMAK)
!CC
!CC   ATTENUATION = 6.0*PI*FREQ*RHO*IMAG(-K)
!CC    6.0*PI/10. = 1.885 (THE FACTOR OF 10 IS FOR UNITS CONVERSION)
!CC
!     GAMFOG=8.1888*FREQ*RHO*AIMAG( -  (CINDEX**2-1)/(CINDEX**2+2))
      GAMFOG=1.885 *FREQ*RHO*AIMAG( -  (CINDEX**2-1)/(CINDEX**2+2))
      RETURN
      END