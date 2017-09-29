#include <misc.h>
#include <params.h>

module cloud_optics

  !Purpose:
  !   get liquid and ice cloud optics information
  !   liquid: Slingo, 1989 Journal of the Atmospheric Sciences
  !   ice: Ebert & Curry, 1992 JGR
  !   Also will do cloud-overlap approximation
  !   Author: D. Feldman

  use ppgrid
  use pmgrid
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  save

contains

  subroutine get_ext_ssa_asym_liq(r_eff,wvl_out,ext_liq_out,ssa_liq_out,asym_liq_out)

    !Will return extinction, absorption, and asymmetry
    !for liquid clouds of a given effective radius
    ! Author: D. Feldman
    !
    ! Inputs 
    ! r_eff = effective radius (in um)
    !
    ! Outputs 
    ! wvl_out = wavelength (in um) possibly in vector form
    ! ext_liq_out = extinction coefficient (km^-1 m^3/g)
    ! ssa_liq_out = single-scattering albedo - 1 (unitless, and negative)
    ! asym_liq_out = liquid water droplet henyey-greenstein scattering phase function asymmetry factor


    use ioFileMod, only: getfil
    implicit none

    real(r8), intent(in) :: r_eff
    real(r8), intent(out) :: wvl_out(24)
    real(r8), intent(out) :: ext_liq_out(24)
    real(r8), intent(out) :: ssa_liq_out(24)
    real(r8), intent(out) :: asym_liq_out(24)

    !local variables for loops
    integer :: i

    !Variables from Table I in Slingo, 1989
    real(r8) :: wvl_low(24) 
    real(r8) :: wvl_high(24) 
    real(r8) :: wvl_mid(24) 
    real(r8) :: a_i(24) 
    real(r8) :: b_i(24) 
    real(r8) :: c_i(24) 
    real(r8) :: d_i(24) 
    real(r8) :: e_i(24) 
    real(r8) :: f_i(24) 

    !Data from Table I of Slingo 1989
    wvl_low(1) = 0.25
    wvl_low(2) = 0.3
    wvl_low(3) = 0.33
    wvl_low(4) = 0.36
    wvl_low(5) = 0.4
    wvl_low(6) = 0.44
    wvl_low(7) = 0.48
    wvl_low(8) = 0.52
    wvl_low(9) = 0.57
    wvl_low(10) = 0.64
    wvl_low(11) = 0.69
    wvl_low(12) = 0.75
    wvl_low(13) = 0.78
    wvl_low(14) = 0.87
    wvl_low(15) = 1.00
    wvl_low(16) = 1.1
    wvl_low(17) = 1.19
    wvl_low(18) = 1.28
    wvl_low(19) = 1.53
    wvl_low(20) = 1.64
    wvl_low(21) = 2.13
    wvl_low(22) = 2.38
    wvl_low(23) = 2.91
    wvl_low(24) = 3.42
    
    wvl_high(1) = 0.3
    wvl_high(2) = 0.33
    wvl_high(3) = 0.36
    wvl_high(4) = 0.4
    wvl_high(5) = 0.44
    wvl_high(6) = 0.48
    wvl_high(7) = 0.52
    wvl_high(8) = 0.57
    wvl_high(9) = 0.64
    wvl_high(10) = 0.69
    wvl_high(11) = 0.75
    wvl_high(12) = 0.78
    wvl_high(13) = 0.87
    wvl_high(14) = 1.00
    wvl_high(15) = 1.1
    wvl_high(16) = 1.19
    wvl_high(17) = 1.28
    wvl_high(18) = 1.53
    wvl_high(19) = 1.64
    wvl_high(20) = 2.13
    wvl_high(21) = 2.38
    wvl_high(22) = 2.91
    wvl_high(23) = 3.42
    wvl_high(24) = 4.0

    wvl_mid(1) = 0.275
    wvl_mid(2) = 0.315
    wvl_mid(3) = 0.345
    wvl_mid(4) = 0.38
    wvl_mid(5) = 0.42
    wvl_mid(6) = 0.46
    wvl_mid(7) = 0.5
    wvl_mid(8) = 0.545
    wvl_mid(9) = 0.605
    wvl_mid(10) = 0.665
    wvl_mid(11) = 0.72
    wvl_mid(12) = 0.765
    wvl_mid(13) = 0.825
    wvl_mid(14) = 0.935
    wvl_mid(15) = 1.05
    wvl_mid(16) = 1.145
    wvl_mid(17) = 1.235
    wvl_mid(18) = 1.405
    wvl_mid(19) = 1.585
    wvl_mid(20) = 1.9
    wvl_mid(21) = 2.255
    wvl_mid(22) = 2.5
    wvl_mid(23) = 3.15
    wvl_mid(24) = 3.7

    a_i(1) = 3.094
    a_i(2) = 2.944
    a_i(3) = 3.308
    a_i(4) = 2.801
    a_i(5) = 2.668
    a_i(6) = 2.698
    a_i(7) = 2.672
    a_i(8) = 2.838
    a_i(9) = 2.831
    a_i(10) = 2.895
    a_i(11) = 3.115
    a_i(12) = 2.65
    a_i(13) = 2.622
    a_i(14) = 2.497
    a_i(15) = 2.632
    a_i(16) = 2.589
    a_i(17) = 2.551
    a_i(18) = 2.463
    a_i(19) = 2.237
    a_i(20) = 1.97
    a_i(21) = 1.85
    a_i(22) = 1.579
    a_i(23) = 1.95
    a_i(24) = -1.023

    b_i(1) = 1.252
    b_i(2) = 1.27
    b_i(3) = 1.246
    b_i(4) = 1.293
    b_i(5) = 1.307
    b_i(6) = 1.315
    b_i(7) = 1.32
    b_i(8) = 1.3
    b_i(9) = 1.317
    b_i(10) = 1.315
    b_i(11) = 1.244
    b_i(12) = 1.349
    b_i(13) = 1.362
    b_i(14) = 1.376
    b_i(15) = 1.365
    b_i(16) = 1.385
    b_i(17) = 1.401
    b_i(18) = 1.42
    b_i(19) = 1.452
    b_i(20) = 1.501
    b_i(21) = 1.556
    b_i(22) = 1.611
    b_i(23) = 1.54
    b_i(24) = 1.933

    c_i(1) = 7.9e-7
    c_i(2) = -6.5e-7
    c_i(3) = -3.0e-7
    c_i(4) = 1.0e-6
    c_i(5) = 0.
    c_i(6) = 1.0e-6
    c_i(7) = 0.
    c_i(8) = 0.
    c_i(9) = -1.2e-6
    c_i(10) = -1.2e-7
    c_i(11) = -2.7e-7
    c_i(12) = 2.3e-6
    c_i(13) = 3.3e-6
    c_i(14) = 9.8e-6
    c_i(15) = -4.6e-5
    c_i(16) = -2.8e-5
    c_i(17) = 6.2e-5
    c_i(18) = 2.4e-4
    c_i(19) = 1.2e-4
    c_i(20) = 1.2e-3
    c_i(21) = 1.9e-4
    c_i(22) = 1.23e-1
    c_i(23) = 4.49e-1
    c_i(24) = 2.5e-2

    d_i(1) = 3.69e-7
    d_i(2) = 4.33e-7
    d_i(3) = 2.36e-7
    d_i(4) = 0.
    d_i(5) = 0.
    d_i(6) = 0.
    d_i(7) = 0.
    d_i(8) = 0.
    d_i(9) = 4.0e-7
    d_i(10) = 4.4e-7
    d_i(11) = 1.4e-6
    d_i(12) = 1.7e-6
    d_i(13) = 2.8e-6
    d_i(14) = 2.1e-5
    d_i(15) = 5.0e-5
    d_i(16) = 8.0e-5
    d_i(17) = 2.6e-4
    d_i(18) = 8.56e-4
    d_i(19) = 6.67e-4
    d_i(20) = 2.16e-3
    d_i(21) = 2.54e-3
    d_i(22) = 9.35e-3
    d_i(23) = 1.54e-3
    d_i(24) = 1.22e-2

    e_i(1) = 0.844
    e_i(2) = 0.841
    e_i(3) = 0.839
    e_i(4) = 0.836
    e_i(5) = 0.840
    e_i(6) = 0.820
    e_i(7) = 0.828
    e_i(8) = 0.825
    e_i(9) = 0.828
    e_i(10) = 0.818
    e_i(11) = 0.804
    e_i(12) = 0.809
    e_i(13) = 0.806
    e_i(14) = 0.783
    e_i(15) = 0.784
    e_i(16) = 0.780
    e_i(17) = 0.773
    e_i(18) = 0.754
    e_i(19) = 0.749
    e_i(20) = 0.740
    e_i(21) = 0.769
    e_i(22) = 0.851
    e_i(23) = 0.831
    e_i(24) = 0.726

    f_i(1) = 1.558
    f_i(2) = 1.68
    f_i(3) = 1.946
    f_i(4) = 2.153
    f_i(5) = 1.881
    f_i(6) = 3.004
    f_i(7) = 2.467
    f_i(8) = 2.776
    f_i(9) = 2.492
    f_i(10) = 2.989
    f_i(11) = 3.520
    f_i(12) = 3.387
    f_i(13) = 3.355
    f_i(14) = 5.035
    f_i(15) = 4.745
    f_i(16) = 4.989
    f_i(17) = 5.405
    f_i(18) = 6.555
    f_i(19) = 6.931
    f_i(20) = 7.469
    f_i(21) = 5.171
    f_i(22) = 2.814
    f_i(23) = 6.102
    f_i(24) = 6.652

    do i=1,24
       wvl_out(i) = wvl_mid(i)
       ! in m^2/g
       ext_liq_out(i) = 0.01*a_i(i) + b_i(i)/r_eff  
       ! convert to km^-1*m^3/g
       ext_liq_out(i) = ext_liq_out(i)*1000.0
       ssa_liq_out(i) = -1.0*(c_i(i)+d_i(i)*r_eff)
       asym_liq_out(i) = e_i(i)+0.001*f_i(i)*r_eff
    end do

  end subroutine  get_ext_ssa_asym_liq


  subroutine get_ext_ssa_asym_ice(r_eff,wvl_in,ext_ice_out,ssa_ice_out,asym_ice_out)

    !Will return extinction, absorption, and asymmetry
    !for ice clouds of a given effective radius
    ! based on data from Ebert & Curry, JGR 1992
    ! Author: D. Feldman
    ! 
    ! r_eff = effective radius (in um)
    ! wvl_in = wavelength grid (in um)
    ! ext_ice_out = extinction coefficient (km^-1 m^3/g)
    ! ssa_ice_out = single-scattering albedo - 1 (unitless, and negative)
    ! asym_ice_out = liquid water droplet henyey-greenstein scattering phase function asymmetry factor


    use ioFileMod, only: getfil
    implicit none

    real(r8), intent(in) :: r_eff
    real(r8), intent(in) :: wvl_in(24)
    real(r8), intent(out) :: ext_ice_out(24)
    real(r8), intent(out) :: ssa_ice_out(24)
    real(r8), intent(out) :: asym_ice_out(24)

    !local variables for loops
    integer :: i,j
    integer :: liq_ice_index

    real(r8) :: wvl_mid(5) 
    real(r8) :: a_i(5) 
    real(r8) :: b_i(5) 
    real(r8) :: c_i(5) 
    real(r8) :: d_i(5) 
    real(r8) :: e_i(5) 
    real(r8) :: f_i(5) 

    !Data from Table II of Ebert & Curry 1992
    wvl_mid(1) = 0.475
    wvl_mid(2) = 1.0
    wvl_mid(3) = 1.6
    wvl_mid(4) = 2.2
    wvl_mid(5) = 3.0

    a_i(1) = 3.448e-3
    a_i(2) = 3.448e-3
    a_i(3) = 3.448e-3
    a_i(4) = 3.448e-3
    a_i(5) = 3.448e-3

    b_i(1) = 2.431
    b_i(2) = 2.431
    b_i(3) = 2.431
    b_i(4) = 2.431
    b_i(5) = 2.431

    c_i(1) = 0.00001
    c_i(2) = 0.00011
    c_i(3) = 0.01240
    c_i(4) = 0.03779
    c_i(5) = 0.46658

    d_i(1) = 0.0
    d_i(2) = 1.405e-5
    d_i(3) = 6.867e-4
    d_i(4) = 1.284e-3
    d_i(5) = 2.05e-5

    e_i(1) = 0.7661
    e_i(2) = 0.7730
    e_i(3) = 0.7865
    e_i(4) = 0.8172
    e_i(5) = 0.9595

    f_i(1) = 5.851e-4
    f_i(2) = 5.665e-4
    f_i(3) = 7.204e-4
    f_i(4) = 7.463e-4
    f_i(5) = 1.076e-4

    do i=1,24
       liq_ice_index = minloc(abs(wvl_in(i)-wvl_mid),1)

       ! in m^2/g
       ext_ice_out(i) = a_i(liq_ice_index) + b_i(liq_ice_index)/r_eff 
       ! convert to km^-1*m^3/g
       ext_ice_out(i) = ext_ice_out(i)*1000.0
       ssa_ice_out(i) = -1.0*(c_i(liq_ice_index)+d_i(liq_ice_index)*r_eff)
       asym_ice_out(i) = e_i(liq_ice_index)+f_i(liq_ice_index)*r_eff
    end do

  end subroutine get_ext_ssa_asym_ice

  subroutine get_ext_ssa_asym_liq_lw(r_eff,wvl_in,ext_liq_out,ssa_liq_out,asym_liq_out)
    !Will return extinction, absorption, and asymmetry for longwave RT calculations
    !for ice clouds of a given effective radius
    ! Refer to NCAR CCSM3.0 description NCAR/TN-464+STR section 4.9.5
    ! Author: D. Feldman
    ! 
    ! r_eff = effective radius (in um)
    ! wvl_in = wavelength grid (in um)
    ! ext_liq_out = extinction coefficient (km^-1 m^3/g)
    ! ssa_liq_out = single-scattering albedo - 1 (unitless, and negative)
    ! asym_liq_out = liquid water droplet henyey-greenstein scattering phase function asymmetry factor


    use ioFileMod, only: getfil
    implicit none

    real(r8), intent(in) :: r_eff
    real(r8), intent(in) :: wvl_in(24)
    real(r8), intent(out) :: ext_liq_out(24)
    real(r8), intent(out) :: ssa_liq_out(24)
    real(r8), intent(out) :: asym_liq_out(24)

    !local variables for loops
    integer :: i

    do i=1,24
       !ext_liq_out(i) = 0.090361*1000. !converts m^2/g to km^-1*m^3/g
       ext_liq_out(i) = 0.000001*1000.
       ssa_liq_out(i) = -0.9999
       asym_liq_out(i) = 0.0 
    end do

  end subroutine get_ext_ssa_asym_liq_lw

  subroutine get_ext_ssa_asym_ice_lw(r_eff,wvl_in,ext_ice_out,ssa_ice_out,asym_ice_out)
    !Will return extinction, absorption, and asymmetry for longwave RT calculations
    !for ice clouds of a given effective radius
    ! based on the Ebert and Curry formulation as adapted for the CAM RT code
    ! Refer to NCAR CCSM3.0 description NCAR/TN-464+STR
    ! Author: D. Feldman
    ! 
    ! r_eff = effective radius (in um)
    ! wvl_in = wavelength grid (in um)
    ! ext_ice_out = extinction coefficient (km^-1 m^3/g)
    ! ssa_ice_out = single-scattering albedo - 1 (unitless, and negative)
    ! asym_ice_out = liquid water droplet henyey-greenstein scattering phase function asymmetry factor


    use ioFileMod, only: getfil
    implicit none

    real(r8), intent(in) :: r_eff
    real(r8), intent(in) :: wvl_in(24)
    real(r8), intent(out) :: ext_ice_out(24)
    real(r8), intent(out) :: ssa_ice_out(24)
    real(r8), intent(out) :: asym_ice_out(24)

    !local variables for loops
    integer :: i

    do i=1,24
       ext_ice_out(i) = (0.005 + 1./r_eff)*1000  !converts m^2/g to km^-1*m^3/g
       ssa_ice_out(i) = -0.9999
       asym_ice_out(i) = 0.0 
    end do

  end subroutine get_ext_ssa_asym_ice_lw

  
   subroutine cldovrlap(pint    ,cld     ,nmxrgn  ,pmxrgn  ) 
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Partitions each column into regions with clouds in neighboring layers.
! This information is used to implement maximum overlap in these regions
! with random overlap between them.
! On output,
!    nmxrgn contains the number of regions in each column
!    pmxrgn contains the interface pressures for the lower boundaries of
!           each region! 
! Method: 
! 
! Author: W. Collins
!         modified by D. Feldman, U.C. Berkeley 06-10-09
!-----------------------------------------------------------------------

    implicit none
!
! Input arguments
!

    real(r8), intent(in) :: pint(pverp)   ! Interface pressure
    real(r8), intent(in) :: cld(pver)     ! Fractional cloud cover
!
! Output arguments
!
    real(r8), intent(out) :: pmxrgn(pverp)! Maximum values of pressure for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pressure for
!    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc
    integer nmxrgn                    ! Number of maximally overlapped regions
!
!---------------------------Local variables-----------------------------
!
    integer i                    ! Longitude index
    integer k                    ! Level index
    integer n                    ! Max-overlap region counter

    real(r8) pnm(pverp)    ! Interface pressure

    logical cld_found            ! Flag for detection of cloud
    logical cld_layer(pver)      ! Flag for cloud in layer
!
!------------------------------------------------------------------------
!

       cld_found = .false.
       cld_layer(:) = cld(:) > 0.0_r8  !Threshold for whether a cloud exists
       pmxrgn(:) = 0.0
       pnm(:)=pint(:)*10.
       n = 1
       do k = 1, pver
          if (cld_layer(k) .and.  .not. cld_found) then
             cld_found = .true.
          else if ( .not. cld_layer(k) .and. cld_found) then
             cld_found = .false.
             if (count(cld_layer(k:pver)) == 0) then
                exit
             endif
             pmxrgn(n) = pnm(k)
             n = n + 1
          endif
       end do
       pmxrgn(n) = pnm(pverp)
       nmxrgn = n

    return
  end subroutine cldovrlap


  subroutine denfac_copy(ACOEF,B2COEF,RATLOG,denfac_out)

!     This subroutine is a copy of Modtran's exponential integral evaluator segment.f DENFAC
!     DENFAC EVALUATES THE INTEGRAL
!          1
!          /       - RATLOG (B2COEF X + 2 ACOEF) X
!          |   EXP                                  DX
!          /
!          0
!
!     THE PRODUCT B2COEF * RATLOG IS POSITIVE UNLESS THERE IS A
!     DENSITY INVERSION.  THE INTEGRAL CAN BE EXPRESSED IN TERMS OF
!     ERROR FUNCTIONS OF REAL ARGUMENTS AS LONG AS THE PRODUCT IS
!     POSITIVE.  IF THE PRODUCT IS NEGATIVE, DAWSON INTEGRALS RESULT.
!     THE INTEGRAL IS EVALUATED ONLY IF B2COEF * RATLOG IS GREATER
!     THAN -0.25;  OTHERWISE, THE QUADRATIC TERM IN THE EXPANSION,
!     B2COEF, IS SET TO ZERO AND A WARNING IS PRINTED IN THE LOG FILE.
!     (THE DAWSON INTEGRAL COULD BE EVAULATED BY INTERPOLATING TABLE
!     7.5 OF ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS)
    implicit none

!     INPUT ARGUMENTS:
!       ACOEF    COEFFICIENT OF LINEAR EXPONENT TERMS.
!       B2COEF   COEFFICIENT OF QUADRATIC EXPONENT TERMS.
!       RATLOG   EXPONENTIAL ARGMENT COEFFICIENT.
      real(r8) ACOEF,B2COEF,RATLOG

      real(r8) denfac_out

!     LOCAL VARIABLES:
!       A        HALF THE LINEAR EXPONENT TERM.
!       TWOA     THE LINE EXPONENT TERM.
!       A2       A SUQARED.
!       E        THE EXPONENTIAL OF NEGATIVE TWOA.
!       B2       THE QUADRATIC EXPONENT TERM.
!       B4       B2 SQUARED.
!       B6       B2 CUBED.
!       B        SQUARE ROOT OF B2.
!       APB2     A PLUS B2.
!       RAT      RATIO OF B2 TO TWOA SQUARED.
      real(r8) A,TWOA,A2,E,B2,B4,B6,B,APB2,RAT
      real(r8) facfnc1,facfnc2
      real(r8) RTPI

!     TEST INTEGRAL DOMAIN:
      RTPI = 1.7724538509055

      A=ACOEF*RATLOG
      B2=B2COEF*RATLOG
      IF(B2.LT.-.25)THEN
         B2=0.
      ENDIF

      IF(B2.LT.-.0001)THEN
          IF(ABS(A).LE..25)THEN
              B4=B2*B2
              B6=B2*B4
              A2=A*A
              denfac_out=1-A                                                &
               +(2*A2-B2                                               &
               +(A*(3*B2-2*A2)                                         &
               +(3*B4-4*A2*(3*B2-A2)                                   &
               -(A*(15*B4-2*A2*(10*B2-A2))                             &
               +(15*B6-A2*(90*B4-A2*(60*B2-8*A2))                      &
               -(A*(105*B6-A2*(210*B4-A2*(84*B2-8*A2)))                &
                                                 )/4)/7)/3)/5)/2)/3
          ELSE
              TWOA=2*A
              E=EXP(-TWOA)
              RAT=B2/TWOA**2
              denfac_out=((1-E)                                             &
               -RAT*(2-E*(2+TWOA*(2+TWOA))                             &
               -RAT*(12-E*(12+TWOA*(12+TWOA*(6+TWOA*(2+A))))           &
               -RAT*(120-E*(120+TWOA*(120+TWOA*(60+TWOA*(20+TWOA*(5    &
                                     +TWOA*(1+TWOA/6))))))))))/TWOA
          ENDIF
      ELSEIF(B2.LE..0001)THEN
          IF(ABS(A).LE..005)THEN
              denfac_out=1-(A-(A*A*(2-A)-B2*(1-1.5*A))/3)
          ELSE
              TWOA=2*A
              E=EXP(-TWOA)
              denfac_out=((1-E)-(2-E*(2+TWOA*(2+TWOA)))*B2/TWOA**2)/TWOA
          ENDIF
      ELSE
        APB2=A+B2
        B=SQRT(B2)
        IF(A.GE.0.D0)THEN
            call facfnc_copy(B,A,facfnc1)
            call facfnc_copy(B,APB2,facfnc2)

            denfac_out= facfnc1-EXP(-A-APB2)*facfnc2
        ELSEIF(APB2.LE.0.D0)THEN
            call facfnc_copy(B,-APB2,facfnc1)
            call facfnc_copy(B,-A,facfnc2)
  
            denfac_out=EXP(-A-APB2)*facfnc1-facfnc2
        ELSE
            call facfnc_copy(B,APB2,facfnc1)
            call facfnc_copy(B,-A,facfnc2)
  
            denfac_out=RTPI*EXP((A/B)**2)/B                                   &
             -EXP(-A-APB2)*facfnc1-facfnc2
        ENDIF
      ENDIF
   return

  end subroutine denfac_copy

  subroutine facfnc_copy(B,Z,facfnc_out)
!----------------------------------------------------------------------- 
! 

  real(r8) B,T,Z,DENOM,facfnc_out

       IF (Z.GT.1.883*B) THEN
          T = (B/Z)**2
          facfnc_out =(8+T*(140+T*(690+T*(975+T*192))))/                     &
                    (Z*(16+T*(288+T*(1512+T*(2520+945*T)))))
       ELSE
          DENOM=B+0.3275911*Z
          T=B/DENOM
          facfnc_out = (0.225836846+T*(-0.252128668+T*(1.259695129+T*(-1.287822453+T*0.940646070))))
       ENDIF
       RETURN

  end subroutine facfnc_copy




  subroutine generate_stochastic_clouds(nlayers, icld, changeSeed,pmid, cldfrac, clwp, ciwp,  &
                               num_permute,cld_stoch,clwp_stoch,ciwp_stoch,iscloudy_matrix) 

!----------------------------------------------------------------------- 
! 
! Purpose: 
!   This subroutine will yield random vertical profiles of liquid and ice cloud water
!   content that is consistent with a maximum-random cloud-overlap approximation given
!   cloud fraction & input profiles of cloud fraction and liquid/ice water content.
!----------------------------------------------------------------------
! Modified: D. Feldman 06-10-09
!----------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------

  !----------------------------------------------------------------------------------------------------------------
  ! ---------------------
  ! Contact: Cecile Hannay (hannay@ucar.edu)
  ! 
  ! Modifications: Generalized for use with RRTMG and added Mersenne Twister as the default
  !   random number generator, which can be changed to the optional kissvec random number generator
  !   with flag 'irnd' below . Some extra functionality has been commented or removed.  
  !   Michael J. Iacono, AER, Inc., February 2007
  !
  ! Given a profile of cloud fraction, cloud water and cloud ice, we produce a set of subcolumns.
  ! Each layer within each subcolumn is homogeneous, with cloud fraction equal to zero or one 
  ! and uniform cloud liquid and cloud ice concentration.
  ! The ensemble as a whole reproduces the probability function of cloud liquid and ice within each layer 
  ! and obeys an overlap assumption in the vertical.   
  ! 
  ! Overlap assumption:
  !  The cloud are consistent with 4 overlap assumptions: random, maximum, maximum-random and exponential. 
  !  The default option is maximum-random (option 3)
  !  The options are: 1=random overlap, 2=max/random, 3=maximum overlap, 4=exponential overlap
  !  This is set with the variable "overlap" 
  !mji - Exponential overlap option (overlap=4) has been deactivated in this version
  !  The exponential overlap uses also a length scale, Zo. (real,    parameter  :: Zo = 2500. ) 
  ! 
  ! Seed:
  !  If the stochastic cloud generator is called several times during the same timestep, 
  !  one should change the seed between the call to insure that the subcolumns are different.
  !  This is done by changing the argument 'changeSeed'
  !  For example, if one wants to create a set of columns for the shortwave and another set for the longwave ,
  !  use 'changeSeed = 1' for the first call and'changeSeed = 2' for the second call 
  !
  ! PDF assumption:
  !  We can use arbitrary complicated PDFS. 
  !  In the present version, we produce homogeneuous clouds (the simplest case).  
  !  Future developments include using the PDF scheme of Ben Johnson. 
  !
  ! History file:
  !  Option to add diagnostics variables in the history file. (using FINCL in the namelist)
  !  nsubcol = number of subcolumns
  !  overlap = overlap type (1-3)
  !  Zo = length scale 
  !  CLOUD_S = mean of the subcolumn cloud fraction ('_S" means Stochastic)
  !  CLDLIQ_S = mean of the subcolumn cloud water
  !  CLDICE_S = mean of the subcolumn cloud ice 
  !
  ! Note:
  !   Here: we force that the cloud condensate to be consistent with the cloud fraction 
  !   i.e we only have cloud condensate when the cell is cloudy. 
  !   In CAM: The cloud condensate and the cloud fraction are obtained from 2 different equations 
  !   and the 2 quantities can be inconsistent (i.e. CAM can produce cloud fraction 
  !   without cloud condensate or the opposite).
  !---------------------------------------------------------------------------------------------------------------

      use mcica_random_numbers
! The Mersenne Twister random number engine
      use MersenneTwister, only: randomNumberSequence, &   
                                 new_RandomNumberSequence, getRandomReal

      type(randomNumberSequence) :: randomNumbers

! -- Arguments

      integer, intent(in) :: nlayers         ! number of layers
      integer, intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer, intent(in) :: changeSeed      ! allows permuting seed
      integer, intent(in) :: num_permute     ! number of permutations

! Column state (cloud fraction, cloud water, cloud ice) + variables needed to read physics state 
      real(r8), intent(in) :: pmid(pver)            ! layer pressure (Pa)
                                                        !    Dimensions: (nlayers)
      real(r8), intent(in) :: cldfrac(pver)             ! cloud fraction 
                                                        !    Dimensions: (nlayers)
      real(r8), intent(in) :: clwp(pver)            ! cloud liquid water path (g/m2)
                                                        !    Dimensions: (nlayers)
      real(r8), intent(in) :: ciwp(pver)            ! cloud ice water path (g/m2)
                                                        !    Dimensions: (nlayers)
!      real(r8), intent(in) :: tauc(1,pver)          ! cloud optical depth (non-delta scaled)
!                                                        !    Dimensions: (nbndsw,nlayers)
!      real(r8), intent(in) :: ssac(1,pver)          ! cloud single scattering albedo (non-delta scaled)
!                                                        !    Dimensions: (nbndsw,nlayers)
!      real(r8), intent(in) :: asmc(1,pver)          ! cloud asymmetry parameter (non-delta scaled)
                                                        !    Dimensions: (nbndsw,nlayers)

      real(r8), intent(out) :: cld_stoch(num_permute,pver)    ! subcolumn cloud fraction 
                                                        !    Dimensions: (ngptsw,nlayers)
      real(r8), intent(out) :: clwp_stoch(num_permute,pver)   ! subcolumn cloud liquid water path
                                                        !    Dimensions: (ngptsw,nlayers)
      real(r8), intent(out) :: ciwp_stoch(num_permute,pver)   ! subcolumn cloud ice water path
                                                        !    Dimensions: (ngptsw,nlayers)
      integer, intent(out) :: iscloudy_matrix(num_permute,pver) !0 if no cloud, 1 if cloud

!      real(r8), intent(out) :: tauc_stoch(num_permute,pver)   ! subcolumn cloud optical depth
!                                                        !    Dimensions: (ngptsw,nlayers)
!      real(r8), intent(out) :: ssac_stoch(num_permute,pver)   ! subcolumn cloud single scattering albedo
!                                                        !    Dimensions: (ngptsw,nlayers)
!      real(r8), intent(out) :: asmc_stoch(num_permute,nlayers)   ! subcolumn cloud asymmetry parameter
!                                                        !    Dimensions: (ngptsw,nlayers)

! -- Local variables
      !integer, parameter :: nsubcol = ngptsw ! number of sub-columns (g-point intervals)
      real(r8) :: cldf(nlayers)                  ! cloud fraction 
    
! Mean over the subcolumns (cloud fraction, cloud water , cloud ice) - inactive
!      real(kind=jprb) :: mean_cld_stoch(nlayers)       ! cloud fraction 
!      real(kind=jprb) :: mean_clwp_stoch(nlayers)      ! cloud water
!      real(kind=jprb) :: mean_ciwp_stoch(nlayers)      ! cloud ice
!      real(kind=jprb) :: mean_tauc_stoch(nlayers)      ! cloud optical depth
!      real(kind=jprb) :: mean_ssac_stoch(nlayers)      ! cloud single scattering albedo
!      real(kind=jprb) :: mean_asmc_stoch(nlayers)      ! cloud asymmetry parameter

! Set overlap
      integer :: overlap                     ! 1 = random overlap, 2 = maximum/random,
                                                        ! 3 = maximum overlap, 
!      real(kind=jprb), parameter  :: Zo = 2500._jprb   ! length scale (m) 
!      real(kind=jprb) :: zm(nlayers)                   ! Height of midpoints (above surface)
!      real(kind=jprb), dimension(nlayers) :: alpha=0.0_jprb ! overlap parameter  

! Constants (min value for cloud fraction and cloud water and ice)
      real(r8), parameter :: cldmin = 1.0e-2! min cloud fraction
!      real(kind=jprb), parameter :: qmin   = 1.0e-10_jprb   ! min cloud water and cloud ice (not used)

! Variables related to random number and seed 
      integer :: irnd                        ! flag for random number generator
                                                        !  0 = kissvec
                                                        !  1 = Mersenne Twister

      real(r8), dimension(num_permute, nlayers) :: CDF, CDF2 ! random numbers
      integer :: seed1, seed2, seed3, seed4          ! seed to create random number (kissvec)
      real(r8) :: rand_num                       ! random number (kissvec)
      integer :: iseed                       ! seed to create random number (Mersenne Twister)
      real(r8) :: rand_num_mt                    ! random number (Mersenne Twister)

! Flag to identify cloud fraction in subcolumns
      logical,  dimension(num_permute, nlayers) :: iscloudy ! flag that says whether a gridbox is cloudy

! Indices
      integer(r8) :: ilev, isubcol, i,j, n         ! indices

!------------------------------------------------------------------------------------------ 
! Set randum number generator to use (0 = kissvec; 1 = mersennetwister)
!      irnd = 0
      irnd = 1

! Pass input cloud overlap setting to local variable
    overlap = icld

! ensure that cloud fractions are in bounds 
! to avoid to get ql_stoch getting to big (ql_stoch = ql/cld).
    cldf(:) = cldfrac(:)
    where (cldf(:) < cldmin)
        cldf(:) = 0.
    end where

! ----- Create seed  --------
   
! Advance randum number generator by changeseed values
      if (irnd.eq.0) then   
! For kissvec, create a seed that depends on the state of the columns. Maybe not the best way, but it works.  
! Must use pmid from bottom four layers. 
         if (pmid(1).lt.pmid(2)) then
            stop 'MCICA_SUBCOL: KISSVEC SEED GENERATOR REQUIRES PMID FROM BOTTOM FOUR LAYERS.'
         endif
         seed1 = (pmid(1) - int(pmid(1)))  * 1000000000
         seed2 = (pmid(2) - int(pmid(2)))  * 1000000000
         seed3 = (pmid(3) - int(pmid(3)))  * 1000000000
         seed4 = (pmid(4) - int(pmid(4)))  * 1000000000
         do i=1,changeSeed
            call kissvec(seed1, seed2, seed3, seed4, rand_num)
         enddo
      elseif (irnd.eq.1) then
         randomNumbers = new_RandomNumberSequence(seed = changeSeed)
      endif 


! ------ Apply overlap assumption --------

! generate the random numbers  

      select case (overlap)

      case(1) 
! Random overlap
! i) pick a random value at every level
  
         if (irnd.eq.0) then 
            do isubcol = 1,num_permute
               do ilev = 1,nlayers
                  call kissvec(seed1, seed2, seed3, seed4, rand_num)  ! we get different random number for each level
                  CDF(isubcol, ilev) = rand_num
               enddo
            enddo
         elseif (irnd.eq.1) then
            do isubcol = 1, num_permute
               do ilev = 1, nlayers
                  rand_num_mt = getRandomReal(randomNumbers)
                  CDF(isubcol,ilev) = rand_num_mt
               enddo
             enddo
         endif

      case(2) 
! Maximum-Random overlap
! i) pick  a random number for top layer.
! ii) walk down the column: 
!    - if the layer above is cloudy, we use the same random number than in the layer above
!    - if the layer above is clear, we use a new random number 

         if (irnd.eq.0) then 
            do isubcol = 1,num_permute
               do ilev = 1,nlayers
                  call kissvec(seed1, seed2, seed3, seed4, rand_num) 
                  CDF(isubcol, ilev) = rand_num
               enddo
            enddo
         elseif (irnd.eq.1) then
            do isubcol = 1, num_permute
               do ilev = 1, nlayers
                  rand_num_mt = getRandomReal(randomNumbers)
                  CDF(isubcol,ilev) = rand_num_mt
               enddo
             enddo
         endif

         do ilev = 2,nlayers
            where (CDF(:, ilev-1) > spread(1.- cldf(ilev-1), dim=1, nCopies=num_permute) )
               CDF(:,ilev) = CDF(:,ilev-1) 
            elsewhere
               CDF(:,ilev) = CDF(:,ilev) *  spread(1.- cldf(ilev-1), dim=1, nCopies=num_permute) 
            end where
         enddo

    case(3) 
! Maximum overlap
! i) pick the same random number at every level  

         if (irnd.eq.0) then 
            do isubcol = 1,num_permute
               call kissvec(seed1, seed2, seed3, seed4, rand_num) 
               do ilev = 1,nlayers
                  CDF(isubcol, ilev) = rand_num
               enddo
            enddo
         elseif (irnd.eq.1) then
            do isubcol = 1, num_permute
                  rand_num_mt = getRandomReal(randomNumbers)
               do ilev = 1, nlayers
                  CDF(isubcol,ilev) = rand_num_mt
               enddo
             enddo
         endif

      end select

 
! -- generate subcolumns for homogeneous clouds -----

      do ilev = 1, nlayers
         iscloudy(:,ilev) = (CDF(:,ilev) >= 1.- spread(cldf(ilev), dim=1, nCopies=num_permute) )
         
         do j=1,num_permute
           if (iscloudy(j,ilev)) then
             iscloudy_matrix(j,ilev) = 1
           else
             iscloudy_matrix(j,ilev) = 0
           endif
         end do
      enddo

! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0

      do ilev = 1, nlayers
         where (iscloudy(:,ilev) )
            cld_stoch(:,ilev) = 1.
         elsewhere (.not. iscloudy(:,ilev) )
            cld_stoch(:,ilev) = 0.
         end where
      enddo

! where there is a cloud, set the subcolumn cloud properties;
! cloud-averaged and not grid-averaged cloud properties are assumed for single-column input;
! do NOT divide by cldf here to convert grid-averaged to cloud-averaged quantities

      do ilev = 1, nlayers
         where ( iscloudy(:,ilev) .and. (spread(cldf(ilev), dim=1, nCopies=num_permute) > 0.) )
            clwp_stoch(:,ilev) = spread(clwp(ilev), dim=1, nCopies=num_permute)
            ciwp_stoch(:,ilev) = spread(ciwp(ilev), dim=1, nCopies=num_permute)
         elsewhere
            clwp_stoch(:,ilev) = 0.
            ciwp_stoch(:,ilev) = 0.
         end where
      enddo
!      do ilev = 1, nlayers
!         do isubcol = 1, num_permute
!            if ( iscloudy(isubcol,ilev) .and. (cldf(ilev) > 0.) ) then
!               n = 1
!               tauc_stoch(isubcol,ilev) = tauc(n,ilev)
!               ssac_stoch(isubcol,ilev) = ssac(n,ilev)
!               asmc_stoch(isubcol,ilev) = asmc(n,ilev)
!            else
!               tauc_stoch(isubcol,ilev) = 0.
!               ssac_stoch(isubcol,ilev) = 1.
!               asmc_stoch(isubcol,ilev) = 1.
!            endif
!         enddo
!      enddo

      end subroutine generate_stochastic_clouds


!-------------------------------------------------------------------------------------------------- 
      subroutine kissvec(seed1,seed2,seed3,seed4,ran_arr)
!-------------------------------------------------------------------------------------------------- 

! public domain code
! made available from http://www.fortran.com/
! downloaded by pjr on 03/16/04
! converted to vector form, functions inlined by pjr,mvr on 05/10/2004

! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123; 
!
      real(r8), intent(inout)  :: ran_arr
      integer, intent(inout) :: seed1,seed2,seed3,seed4
!      integer(kind=jpim) :: i,sz,kiss
      integer :: kiss
      integer :: m, k, n

! inline function 
      m(k, n) = ieor (k, ishft (k, n) )

!      sz = size(ran_arr)
!      do i = 1, sz
         seed1 = 69069 * seed1 + 1327217885
         seed2 = m (m (m (seed2, 13), - 17), 5)
         seed3 = 18000 * iand (seed3, 65535) + ishft (seed3, - 16)
         seed4 = 30903 * iand (seed4, 65535) + ishft (seed4, - 16)
         kiss = seed1 + seed2 + ishft (seed3, 16) + seed4
         ran_arr = kiss*2.328306e-10+ 0.5
!     end do
    
      end subroutine kissvec


!  end subroutine generate_stochastic_clouds 

end module cloud_optics 
