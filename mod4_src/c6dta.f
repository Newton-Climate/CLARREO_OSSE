      REAL FUNCTION C6DTA(V)

!     RETURNS MOLECULAR RAYLEIGH SCATTERING COEFFICIENT [KM-1 @ 273K
!     AND 1 ATM] USING 1999 APPROXIMATION FROM ERIC P. SHETTLE (NRL).

!     INPUT ARGUMENTS:
!       V       SPECTRAL FREQUENCY [CM-1].
      REAL V

!     PARAMETERS:
!       A0-A8   POLYNOMIAL EXPANSION COEFFICIENTS:
      REAL A0,A2,A4,A6,A8
      PARAMETER(A0=+9.38508E+18,A2=-1.09977E+09,                        &
     &  A4=-3.45167E-02,A6=1.07430E-11,A8=-3.69648E-21)

!     LOCAL VARIABLES:
!       V2      SPECTRAL FREQUENCY SQUARED [CM-2].
      REAL V2

!     SCATTERING COEFFICIENT:
      V2=V*V
      C6DTA=V2*V2/(A0+V2*(A2+V2*(A4+V2*(A6+V2*A8))))
      RETURN
      END