      FUNCTION   AB(WAVL,A,CEN,B,C)
!CC
!CC    DESCRIBES THE IMAGINARY PART OF THE DIELECTRIC CONSTANT
!CC
      AB=-A*EXP(-ABS((LOG10(10000.*WAVL/CEN)/B))**C)
      RETURN
      END
