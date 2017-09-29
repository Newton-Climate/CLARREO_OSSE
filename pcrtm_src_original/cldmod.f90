      MODULE CLDMOD

!         CLDMOD CONTAINS WATER CLOUD DATA FROM AN EXTERNAL FILE.
!           CLDTIT   WATER CLOUD TYPE OR NAME.
!           NCLDAN   NUMBER OF WATER CLOUD SCATTERING PHASE
!                    FUNCTION ANGULAR GRID POINTS.
!           NCLDLG   NUMBER OF WATER CLOUD SCATTERING PHASE FUNCTION
!                    LEGENDRE POLYNOMIAL EXPANSION COEFFICIENTS.
!           NCLDWV   NUMBER OF WATER CLOUD SPECTRAL GRID POINTS.
!           ICLDAN   LOOP INDEX FOR WATER CLOUD SCATTERING PHASE
!                    FUNCTION ANGULAR GRID POINTS.
!           ICLDLG   LOOP INDEX FOR WATER CLOUD SCATTERING PHASE
!                    FUNCTION LEGENDRE POLYNOMIAL EXPANSION COEFFICIENTS.
!           ICLDWV   LOOP INDEX FOR WATER CLOUD SPECTRAL GRID POINTS.
!           CLDANG   WATER CLOUD SCATTERING PHASE FUNCTION ANGULAR
!                    GRID [DEGREES].
!           CLDCOS   WATER CLOUD SCATTERING PHASE FUNCTION COSINE GRID.
!           CLDWAV   WATER CLOUD SCATTERING PHASE FUNCTION SPECTRAL
!                    GRID [MICRONS].
!           CLDEXT   WATER CLOUD EXTINCTION COEFFICIENTS [KM-1 M3 / G].
!           CLDABS   WATER CLOUD ABSORPTION COEFFICIENTS [KM-1 M3 / G].
!           CLDPF    WATER CLOUD PHASE FUNCTION [SR-1, 4 pi NORMALIZED].
!           CLDLEG   WATER CLOUD SCATTERING PHASE FUNCTION LEGENDRE
!                    POLYNOMIAL EXPANSION COEFFICIENTS.
!           CLDPFV   CIRRUS CLOUD PHASE FUNCTION INTERPOLATED TO A
!                    SPECTRAL WAVELENGTH [SR-1, 4 pi NORMALIZED].
          CHARACTER CLDTIT*80
          INTEGER NCLDAN, NCLDLG, NCLDWV, ICLDAN, ICLDLG, ICLDWV
          REAL, ALLOCATABLE :: CLDANG(:), CLDCOS(:), CLDWAV(:),         &
     &      CLDEXT(:), CLDABS(:), CLDPF(:,:), CLDLEG(:,:), CLDPFV(:)
      END MODULE CLDMOD
