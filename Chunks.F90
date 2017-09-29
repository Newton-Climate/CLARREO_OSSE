module Chunks

  use shr_kind_mod, only: r8 => shr_kind_r8                                   
  use ppgrid
  use pmgrid
  use RadInput
  use prescribed_aerosols, only: naer_all

  implicit none
  save
!
! Fields input from history file
!
! 3D
!
  real(r8) :: c_allaer(plon,pver,plat,naer_all) ! Aerosol mixing ratio
  real(r8) :: c_cfc11(plon,pver,plat)     ! CFC11 mixing ratio
  real(r8) :: c_cfc12(plon,pver,plat)     ! CFC12 mixing ratio
  real(r8) :: c_ch4(plon,pver,plat)       ! CH4 mixing ratio
  real(r8) :: c_cloud(plon,pver,plat)     ! cloud amount
  real(r8) :: c_emis(plon,pver,plat)      ! Cloud longwave emissivity
  real(r8) :: c_icldiwp(plon,pver,plat)   ! in-cloud ice water path
  real(r8) :: c_icldlwp(plon,pver,plat)   ! in-cloud total water path
  real(r8) :: c_n2o(plon,pver,plat)       ! N2O mixing ratio
  real(r8) :: c_o3vmr(plon,pver,plat)     ! Ozone volume mixing ratio
  real(r8) :: c_q(plon,pver,plat)         ! Water vapor mixing ratio
  real(r8) :: c_rh(plon,pver,plat)         ! Water vapor mixing ratio
  real(r8) :: c_rel(plon,pver,plat)       ! Liquid cloud particle eff. radius
  real(r8) :: c_rei(plon,pver,plat)       ! Ice effective drop size (microns)
  real(r8) :: c_t(plon,pver,plat)         ! Atmospheric temperature
  real(r8) :: c_t_cld(plon,pver,plat)     ! Atmospheric temperature for clouds
!
! 2D
!
  real(r8) :: c_asdir(plon,plat)          ! Albedo: shortwave, direct
  real(r8) :: c_asdif(plon,plat)          ! Albedo: shortwave, diffuse
  real(r8) :: c_aldir(plon,plat)          ! Albedo: longwave, direct
  real(r8) :: c_aldif(plon,plat)          ! Albedo: longwave, diffuse
  real(r8) :: c_icefrac(plon,plat)        ! Fractional ice amount
  real(r8) :: c_landfrac(plon,plat)       ! Fractional land amount
  real(r8) :: c_landmcos(plon,plat)       ! landm*cos(lat) field for rel/rei
  real(r8) :: c_ps(plon,plat)             ! Surface pressure
  real(r8) :: c_snowh(plon,plat)          ! Snow height (equivalent water depth)
  real(r8) :: c_lwup(plon,plat)           ! Upwelling flux from surface
  real(r8) :: c_ts(plon,plat)             ! Radiative surface temperature (DRF)
  real(r8) :: c_windspeed(plon,plat)      ! surface wind speed (m/s) (DRF)
!
! 1D
!
  real(r8) :: c_lat(plat)                 ! latitude
  real(r8) :: c_lon(plon)                 ! longitude
!
!
CONTAINS

  subroutine transpose_input(build_aermmr, build_trace, build_emis, &
                             build_re, build_ozone)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Transpose input data onto internal arrays with parallelized dimension
!   (latitude) as the slowest varying dimension
! 
! Method: 
! Standard F90 transpose and array copy calls
! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------

    use RadInput

    implicit none

!
! Input arguments
!
    character(len=16), intent(in) :: build_aermmr  ! Build AERMMR internally
    logical, intent(in) :: build_trace       ! Build CFCs,CH4,& N2O internally
    logical, intent(in) :: build_emis        ! Build EMIS internally
    logical, intent(in) :: build_re          ! Build RE internally
    logical, intent(in) :: build_ozone       ! Build O3 internally 
!
! Local variables
!
    integer ilon                            ! Longitude index
    integer iaer                            ! Aerosol index

!
! Do Transpose
!
    do ilon = 1, plon
!
! 3D fields
!
       if (build_aermmr == 'NONE') then
          do iaer = 1, naer_all
             c_allaer(ilon,:,:,iaer) =    transpose(inp_allaer(ilon,:,:,iaer))
          end do
       endif
       if (build_aermmr == 'IPCC') then
!
! inp_aermass contains sulfates 
! Determine index for sulfates
!
          iaer = minval(minloc(abs(aerosol_index - idxSUL)))
          c_allaer(ilon,:,:,iaer) =    transpose(inp_aermass(ilon,:,:))
       endif
       if (.not. build_trace) then
          c_cfc11(ilon,:,:) =     transpose(inp_cfc11(ilon,:,:))
          c_cfc12(ilon,:,:) =     transpose(inp_cfc12(ilon,:,:))
          c_ch4(ilon,:,:) =       transpose(inp_ch4(ilon,:,:))
          c_n2o(ilon,:,:) =       transpose(inp_n2o(ilon,:,:))
       endif
       if (.not. build_emis) then
          c_emis(ilon,:,:) =      transpose(inp_emis(ilon,:,:))
       endif
       if (.not. build_re) then
          c_rel(ilon,:,:) =       transpose(inp_rel(ilon,:,:))
          c_rei(ilon,:,:) =       transpose(inp_rei(ilon,:,:))
       endif
       if (.not. build_ozone) then
          c_o3vmr(ilon,:,:) =     transpose(inp_o3vmr(ilon,:,:))
       endif

       c_cloud(ilon,:,:) =     transpose(inp_cloud(ilon,:,:))
       c_icldiwp(ilon,:,:) =   transpose(inp_icldiwp(ilon,:,:))
       c_icldlwp(ilon,:,:) =   transpose(inp_icldlwp(ilon,:,:))
       c_q(ilon,:,:) =         transpose(inp_q(ilon,:,:))
       c_rh(ilon,:,:) =        transpose(inp_rh(ilon,:,:))/100.
       c_t(ilon,:,:) =         transpose(inp_t(ilon,:,:))
       c_t_cld(ilon,:,:) =     transpose(inp_t_cld(ilon,:,:))
    end do
!
! 2D fields
!
    c_asdir =       inp_asdir
    c_asdif =       inp_asdif
    c_aldir =       inp_aldir
    c_aldif =       inp_aldif
    !if (build_re) then !DRF
       c_landfrac =    inp_landfrac
       !c_landmcos =    inp_landmcos
       c_icefrac  =    inp_icefrac
       c_snowh    =    inp_snowh
    !endif !DRF
    c_ps =          inp_ps
    c_lwup =        inp_lwup       
    c_ts   =        inp_ts
    c_windspeed =   inp_windspeed
!
! 1D fields
!
    c_lat = inp_lat
    c_lon = inp_lon

    return

  end subroutine transpose_input

end module Chunks
