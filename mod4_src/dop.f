      FUNCTION   DOP(WAVL,A,CEN1,B,C,CEN2,D,E,CEN3,F,G)
!CC
!CC    DESCRIBES THE REAL PART OF THE DIELECTRIC CONSTANT
!CC
      V=1./WAVL
      V2=V*V
      H1=CEN1**2-V2
      H2=CEN2**2-V2
      H3=CEN3**2-V2
      DOP=SQRT(A+B*H1/(H1*H1+C*V2)+D*H2/(H2*H2+E*V2)+F*H3/(H3*H3+G*V2))
      RETURN
      END