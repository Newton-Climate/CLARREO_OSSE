#include <misc.h>
#include <params.h>

subroutine radctl(lchnk   ,ncol    , in2o,ich4,if11,if12,       &
                  lwup    ,ts, emis    ,          &
                  pmid    ,pint    ,pmln    ,piln    ,t       , &
                  qm1     ,cld     ,cicewp  ,cliqwp  ,coszrs  , jday_input,&
                  asdir   ,asdif   ,aldir   ,aldif   ,pmxrgn  , &
                  nmxrgn  ,fsns    ,fsnt    ,fsnirtoa, flns    ,flnt    , &
                  qrs     ,qrl     ,flwds   ,rel     ,rei     , &
                  sols    ,soll    ,solsd   ,solld   , &
#ifndef OFFLINE
                  landfrac,zm      ,state, fsds)
#else
                  build_ozone, state   ,aerosol ,calday  ,o3vmr   ,&
                  flnsc   ,flntc   ,fsds    , &
                  fsdsc   ,fsnsc   ,fsntc   ,fsnirtoac, solin   , &
                  fln200  ,fln200c ,fsn200  ,fsn200c , &
                  fln     ,fsn     ,aertauout, sol_ann_mean, no_o2_abs, pbr,pnm,eccf,o3mmr, &
		  esat, qsat, rh, latval,lonval,co2vmr, &
                  wavelength_lres, radiance_lres_clr,radiance_lres_all,wavelength_hres, &
                  radiance_hres_clr,radiance_hres_all,solar_flux,land_flag,icefrac,fsno_read,brdf_param,emis_array, &
                  ocean_reflectance,snow_reflectance, solar_az,&
                  modis_asdir,modis_asdif,modis_aldir,modis_aldif,phase_functions,phase_rhs, &
                  phase_wvls,phase_angles,diffuse_flux_clr,diffuse_flux_all,solar_zenith,ipath_in,mcdate, &
                  aod, tau_sulf,tau_dust,tau_carb,tau_sslt, &
                  cldtau_ice_sw,cldtau_liq_sw, cldtau_lw,&
                  bb_updiffuse_clr, bb_updiffuse_all, &
                  bb_dndiffuse_clr, bb_dndiffuse_all, &
                  bb_dndirect_clr, bb_dndirect_all, &
                  flx_updiffuse_clr, flx_updiffuse_all, &
                  flx_dndiffuse_clr, flx_dndiffuse_all, &
                  flx_dndirect_clr, flx_dndirect_all, &
                  lw_flag,windspeed,latitude_index,my_longitude_start,my_longitude_end)
#endif

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Driver for radiation computation.
! 
! Method: 
! Radiation uses cgs units, so conversions must be done from
! model fields to radiation fields.
!
! Author: CCM1,  CMS Contact: J. Truesdale
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use pspect
   use commap
   use history, only: outfld
   use constituents, only: ppcnst, cnst_get_ind
#ifndef OFFLINE
   use prescribed_aerosols, only: get_aerosol, naer_all, aerosol_diagnostics, &
      aerosol_indirect, get_rf_scales, get_int_scales, radforce, idxVOLC
#else
   !use prescribed_aerosols, only: naer_all, radforce, idxVOLC
   use prescribed_aerosols, only: naer_all, radforce, idxBG, idxSUL, idxSSLT, idxOCPHO, idxBCPHO, idxOCPHI, idxBCPHI, &
     idxDUSTfirst, numDUST, idxVOLC, naer_all
#endif

   use physics_types, only: physics_state
   use wv_saturation, only: aqsat
   use chemistry,    only: trace_gas
   use physconst, only: cpair, epsilo
   use aer_optics, only: idxVIS
   use cloud_optics
#ifndef OFFLINE
   use aerosol_intr, only: set_aerosol_from_prognostics
#endif
   use aer_optics, only: nrh, ndstsz, ksul, wsul, gsul, &
     ksslt, wsslt, gsslt, kcphil, wcphil, gcphil, kcphob, wcphob, gcphob, &
     kcb, wcb, gcb, kdst, wdst, gdst, kbg, wbg, gbg, kvolc, wvolc, gvolc
   use write_modtran_tp5
   use ZJIN_OCEAN_WGTLUT, only:zjin_wgtlut 

   implicit none

#include <ptrrgrid.h>
#include <comctl.h>
#include <comsol.h>
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier (latitude index DRF)
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: in2o, ich4, if11, if12 ! indexes of gases in constituent array DRF

   integer nspint            ! Num of spctrl intervals across solar spectrum
   integer naer_groups       ! Num of aerosol groups for optical diagnostics
   parameter ( nspint = 19 )
   parameter ( naer_groups = 7 )    ! current groupings are sul, sslt, all carbons, all dust, background, and all aerosols


   real(r8), intent(in) :: lwup(pcols)          ! Longwave up flux at surface
   real(r8), intent(in) :: ts(pcols)            ! Radiative surface temperature (DRF)
   real(r8), intent(in) :: emis(pcols,pver)     ! Cloud emissivity
   real(r8), intent(in) :: pmid(pcols,pver)     ! Model level pressures
   real(r8), intent(in) :: pint(pcols,pverp)    ! Model interface pressures
   real(r8), intent(in) :: pmln(pcols,pver)     ! Natural log of pmid
   real(r8), intent(in) :: rel(pcols,pver)      ! liquid effective drop size (microns)
   real(r8), intent(in) :: rei(pcols,pver)      ! ice effective drop size (microns)
   real(r8), intent(in) :: piln(pcols,pverp)    ! Natural log of pint
   real(r8), intent(in) :: t(pcols,pver)        ! Model level temperatures
   real(r8), intent(in) :: qm1(pcols,pver,ppcnst) ! Specific humidity and tracers
   real(r8), intent(in) :: cld(pcols,pver)      ! Fractional cloud cover
   real(r8), intent(in) :: cicewp(pcols,pver)   ! in-cloud cloud ice water path
   real(r8), intent(in) :: cliqwp(pcols,pver)   ! in-cloud cloud liquid water path
   real(r8), intent(in) :: coszrs(pcols)        ! Cosine solar zenith angle
   real(r8), intent(in) :: solar_az(pcols)      ! Solar azimuth angle (in degrees) from Wiscombe calculator (DRF)
   integer, intent(in) :: jday_input            ! Julian day input
   real(r8), intent(in) :: asdir(pcols)         ! albedo shortwave direct
   real(r8), intent(in) :: asdif(pcols)         ! albedo shortwave diffuse
   real(r8), intent(in) :: aldir(pcols)         ! albedo longwave direct
   real(r8), intent(in) :: aldif(pcols)         ! albedo longwave diffuse
#ifndef OFFLINE
   real(r8), intent(in) :: landfrac(pcols)      ! land fraction
   real(r8), intent(in) :: zm(pcols,pver)       ! Height of midpoints (above surface)
#else
   real(r8), intent(in) :: calday               ! current calendar day
   real(r8), intent(in) :: o3vmr(pcols,pverr)   ! Ozone volume mixing ratio
   logical,  intent(in) :: build_ozone          ! Flag for reconstruction of O3
   logical,  intent(in) :: sol_ann_mean         ! Flag for setting eccf = 1.0
   logical,  intent(in) :: no_o2_abs            ! Flag to set O2 SW abs to 0
#endif
   real(r8), intent (in) ::  pbr(pcols,pverr)     ! Model mid-level pressures (dynes/cm2) (DRF)
   real(r8), intent (in) ::  pnm(pcols,pverrp)    ! Model interface pressures (dynes/cm2) (DRF)
   real(r8), intent (in) :: eccf                  ! Earth/sun distance factor (DRF)
   real(r8), intent (in) :: o3mmr(pcols,pverr)   ! Ozone mass mixing ratio
   real(r8), intent (in) ::  esat(pcols,pverr)    ! saturation vapor pressure (DRF)
   real(r8), intent (in) ::  qsat(pcols,pverr)    ! saturation specific humidity (DRF)
   real(r8), intent (in) ::  rh(pcols,pverr)      ! level relative humidity (fraction) (DRF)
   real(r8), intent (in) :: latval                ! input latitude in degrees (DRF)
   real(r8), intent (in) :: lonval(pcols)         ! input longitude in degrees (DRF)
   real(r8), intent (in) :: co2vmr                ! input co2vmr (DRF)
   logical, intent (in) :: land_flag(pcols)       ! input flag == true if land (DRF)
   real(r8), intent (in) :: icefrac(pcols)        ! inputfrac of ice (DRF)
   real(r8), intent(in) :: fsno_read(pcols)        ! snow fraction
   real(r8), intent (in) :: brdf_param(pcols,7,3)  ! input brdf properties (DRF)
   real(r8), intent (in) :: emis_array(pcols,6)    ! input surface emissivity from MODIS (DRF)
   real(r8), intent (in) :: ocean_reflectance(24,2) ! input reflectance of ocean (DRF)
   real(r8), intent (in) :: snow_reflectance(969,4) ! input snow reflectance (DRF)
   real(r8), intent (in) :: modis_asdir(pcols)      !MODIS direct SW albedo
   real(r8), intent (in) :: modis_asdif(pcols)      !MODIS diffuse SW albedo
   real(r8), intent (in) :: modis_aldir(pcols)      !MODIS direct LW albedo
   real(r8), intent (in) :: modis_aldif(pcols)      !MODIS diffuse LW albedo
   real(r8), intent (in) :: phase_functions(50,15,8,4) ! input modtran phase functions (DRF)
   real(r8), intent (in) :: phase_rhs(8) !input relative humidity values for phase functions
   real(r8), intent (in) :: phase_wvls(15) !phase function wavelengths (in um)
   real(r8), intent (in) :: phase_angles(50) !phase function angles (in radians)
   character, intent(in) :: ipath_in(120)        ! path to input data
   integer, intent(in) :: mcdate   !date from input_file

   type(physics_state), intent(in) :: state     
   real(r8), intent(inout) :: pmxrgn(pcols,pverp) ! Maximum values of pmid for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pmid for
!    1st region, pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc
   integer, intent(inout) :: nmxrgn(pcols)     ! Number of maximally overlapped regions
   integer, intent(in)    :: latitude_index    ! The latitude index we are presently on in Main
   integer, intent(in)    :: my_longitude_start    ! The latitude index we are presently on in Main
   integer, intent(in)    :: my_longitude_end    ! The latitude index we are presently on in Main
   real(r8), intent(in)   :: windspeed(pcols)    !surface wind-speed derived from cam history file (DRF)
   logical, intent(in)    :: lw_flag             ! true for LW calculations, false for SW calculations



    !BEGIN DRF ADDITIONS
    real real_windspeed
    logical ice_flag_local
    logical no_aerosol_flag  !set to true to turn aerosols off
    logical no_sulf_flag,no_dust_flag,no_carb_flag,no_sslt_flag  !Flags for individual aerosol tests
    logical ocean_flag
    logical snow_flag_local(pcols)        ! input flag == true if snow (DRF)
    real(r8) pflx(pcols,0:pverp)  !interface pressures, including extra layer !DRF
    real(r8) gravit !DRF
    real*4 :: full_atm(pver,28)             !DRF 
    real(r8) :: dummy3(pver)         !DRF
    real(r8) :: dummy4(pver)         !DRF

    !scheme for offsetting Modtran's interpolation
    real(r8) :: total_modtran_extinction
    real(r8) :: scale_fac !for transfer from 495 to 550 nm
    real(r8) :: zp_beg, zp_end
    real(r8) :: halfdr,ecdbeg,ecdend,ecdavg,acoef,abeg,aend,b2beg,b2end,ratlog,denln,denfac_beg,denfac_end
    real(r8) :: aerosol_interpolated(pver)

    real(r8) :: angstrom_exponent     !DRF
    real*4 :: g
    real(r8) :: wvl_ratio !DRF for angstrom exponent
    real*4 aer_wavelen(42)              !DRF
    real*4 cthik_input,cext_input    !DRF 
    integer cam_band_mapping(42)     !DRF
    integer num_levs                 !DRF
    integer num_hres,num_lres        !DRF
    integer aerosol_group_trans(12)  !DRF
    real(r8) :: pmxrgnrf(pcols,pverp)             ! temporary copy of pmxrgn
    integer  :: nmxrgnrf(pcols)     ! temporary copy of nmxrgn
    integer dummy2 !DRF
    integer brdf_len !DRF
    integer spectral_albedo_len
    integer num_ac_wvl !DRF
    integer clr_flx_index !DRF
    real(r8) :: kaerosol(12,nspint)  !DRF extinction matrix
    real(r8) :: waerosol(12,nspint)   !DRF single-scattering albedo matrix
    real(r8) :: gaerosol(12,nspint)   !DRF asymmetry matrix
    real(r8) :: cld_local(pcols,pver)   ! cloud fraction (for manipulating)
    real(r8) :: emis_local(pcols,pver)   ! cloud fraction (for manipulating)
    real(r8) :: cicewp_local(pcols,pver)   ! in-cloud cloud ice water path
    real(r8) :: cliqwp_local(pcols,pver)   ! in-cloud cloud liquid water path
    real(r8) :: t_local(pcols,pver)        ! temperature (copy of 't' but editable) (DRF)
    real(r8) :: h2o_local(pcols,pver)      ! h2o (copy of 'h2o' but editable) (DRF)
    real(r8) :: o3_local(pcols,pver)       ! o3 (copy of 'o3vmr' but editable) (DRF)
    real(r8) :: ch4_local(pcols,pver)      ! ch4 (copy of 'ch4' but editable) (DRF)
    real(r8) :: n2o_local(pcols,pver)      ! n2o (copy of 'n2o' but editable)
    real(r8) :: cfc11_local(pcols,pver)    ! cfc11 (copy of 'cfc11' but editable) (DRF)
    real(r8) :: cfc12_local(pcols,pver)    ! cfc12 (copy of 'cfc12' but editable) (DRF)
    real(r8) :: rei_local(pcols,pver)      ! ice effective drop size (microns)
    real(r8) :: liq_ice_tolerance        ! value for separating liquid and ice clouds (DRF)
    real(r8) :: liq_ice_diff             ! separation between cliqwp and icewp (DRF)
    integer :: max_cld_index             !highest vertical layer with clouds present
    integer :: num_cld_layer             !number of cloud layers (for adding an optically thin cloud)
    logical :: top_cld_ice                !true of highest cloud is ice, false if liquid
    real(r8) :: lowest_t                  !lowest level of temperature that is not undefined (DRF)
    real(r8) :: lowest_h2o                !lowest level of h2o that is not undefined (DRF)
    real(r8) :: lowest_o3                !lowest level of o3 that is not undefined (DRF)
    real(r8) :: lowest_ch4                !lowest level of ch4 that is not undefined (DRF)
    real(r8) :: lowest_n2o                !lowest level of n2o that is not undefined (DRF)
    real(r8) :: lowest_cfc11                !lowest level of cfc11 that is not undefined (DRF)
    real(r8) :: lowest_cfc12                !lowest level of cfc12 that is not undefined (DRF)
    real(r8) :: surface_offset		!If hypsometric equation leads to altitude<0, offset all altitudes by this amount (DRF)

    real*4 brdf_wvl(15) !DRF
    real*4 brdf_param2(15,3) !For entry into driver.f !DRF
    real*4 spec_emissivity_wvl(6) !DRF (must by monotonic increasing in wavelength (um))
    real*4 spec_emissivity_vals(6) !DRF
    integer spec_emissivity_len

    real*4, allocatable :: tmp_rad_lres_clr_noice(:)
    real*4, allocatable :: tmp_rad_lres_clr_ice(:)
    real*4, allocatable :: tmp_rad_hres_clr_noice(:)
    real*4, allocatable :: tmp_rad_hres_clr_ice(:)
    real*4, allocatable :: tmp_rad_lres_all_noice(:)
    real*4, allocatable :: tmp_rad_lres_all_ice(:)
    real*4, allocatable :: tmp_rad_hres_all_noice(:)
    real*4, allocatable :: tmp_rad_hres_all_ice(:)

    real*4, allocatable :: tmp_bbupdif_clr_noice(:)
    real*4, allocatable :: tmp_bbupdif_clr_ice(:)
    real*4, allocatable :: tmp_bbupdif_all_noice(:)
    real*4, allocatable :: tmp_bbupdif_all_ice(:)

    real*4, allocatable :: tmp_bbdndif_clr_noice(:)
    real*4, allocatable :: tmp_bbdndif_clr_ice(:)
    real*4, allocatable :: tmp_bbdndif_all_noice(:)
    real*4, allocatable :: tmp_bbdndif_all_ice(:)

    real*4, allocatable :: tmp_bbdndir_clr_noice(:)
    real*4, allocatable :: tmp_bbdndir_clr_ice(:)
    real*4, allocatable :: tmp_bbdndir_all_noice(:)
    real*4, allocatable :: tmp_bbdndir_all_ice(:)

    real*4, allocatable :: tmp_flxupdif_clr_noice(:,:)
    real*4, allocatable :: tmp_flxupdif_clr_ice(:,:)
    real*4, allocatable :: tmp_flxupdif_all_noice(:,:)
    real*4, allocatable :: tmp_flxupdif_all_ice(:,:)

    real*4, allocatable :: tmp_flxdndif_clr_noice(:,:)
    real*4, allocatable :: tmp_flxdndif_clr_ice(:,:)
    real*4, allocatable :: tmp_flxdndif_all_noice(:,:)
    real*4, allocatable :: tmp_flxdndif_all_ice(:,:)

    real*4, allocatable :: tmp_flxdndir_clr_noice(:,:)
    real*4, allocatable :: tmp_flxdndir_clr_ice(:,:)
    real*4, allocatable :: tmp_flxdndir_all_noice(:,:)
    real*4, allocatable :: tmp_flxdndir_all_ice(:,:)

    real*4, allocatable :: tmp_difflux_clr_noice(:)
    real*4, allocatable :: tmp_difflux_clr_ice(:)
    real*4, allocatable :: tmp_difflux_all_noice(:)
    real*4, allocatable :: tmp_difflux_all_ice(:)

    real(r8), allocatable :: brdf_wvl_dble(:) !DRF for calculating modis_asdir_local, etc
    real(r8), allocatable :: brdf_param_dble(:,:) !DRF for calculating modis_asdir_local, etc
    real*4 :: pfs(50,15,4) !DRF phase functions at appropriate relative humidity
    real*4 :: pfs_sum(50,15,4) !DRF phase functions at appropriate relative humidity
    real(r8) :: max_rh(4) !DRF
    integer krh               ! relative humidity bin index
    real(r8) wrh              ! weight for linear interpolation between lut points
    real(r8) :: rhtrunc       ! rh, truncated for the purposes of extrapolating
    real(r8) :: rh_vec(1000) !DRF
    integer :: max_rh_index !DRF
    real*4 :: tau_weight(42)  !used in o_aerosol calculation
    real*4 :: tau_total(42)   !used in o_aerosol calculation
    real*4 :: w_tau_weight(42)  !used in o_aerosol for ssa weighting
    real*4 :: g_tau_weight(42)  !used in o_aerosol for asymmetry weighting
    real*4 :: g_tau_total(42) ! used in o_aerosol for asymmetry weighting
    real*4 :: w_tau_total(42)   !used in o_aerosol for ssa weighting
    real*4 :: w_tau_weight2  !used for phase-function scaling 
    real*4 :: w_tau_total2(42)   !used phase-function scaling
    integer :: pfs_rh_index    !index for matching phase function relative humidity
    integer :: pfs_wvl_index    !index for matching phase function wavelength
    real*4 :: dummy_dust !DRF
    real*4 :: dummy_carb,dummy_soot2 !DRF
    real*4 :: tau_aerosol_species(4) !DRF
    real(r8) :: dummy_sslt
    real(r8) :: dummy_sulf
    character :: dummy_char !DRF
    character(len=23) :: fname_spectrum ! path to write to netcdf file (DRF)
    integer :: fname_index_first !DRF
    integer :: fname_index_last !DRF
    logical :: fname_flag !DRF
    integer :: fname_length !DRF
    integer :: fname_length_cld !DRF
    integer :: mcdate_local !DRF
    real*4  :: latval_real                           ! input latitude in degrees (DRF)
    real*4  :: lonval_real                           ! input longitude in degrees (DRF)
    real*4 :: co2vmr_real                ! input co2vmr (DRF)
    real*4 :: gndalt_real                ! ground altitude (for submission to driver) DRF
    real*4 :: surftemp_real                ! surface temp (for submission to driver) DRF
    real*4 :: gmt_real                ! gmt time (for submission to driver) DRF
    real*4 :: zen_real                ! zenith_angle (for submission to driver) DRF
    real*4 :: cos_zen_real            ! cosine zenith_angle (for submission to driver) DRF
    real*4 :: v1_real                ! v1 in wavelength (for submission to driver) DRF
    real*4 :: v2_real                ! v2 in wavelengt (for submission to driver) DRF
    real*4 wvl_cld_real(24)
    real*4 liq_ext_avg_real(24),liq_ssa_avg_real(24),liq_asym_avg_real(24)
    real*4 ice_ext_avg_real(24),ice_ssa_avg_real(24),ice_asym_avg_real(24)
    real*4 :: phase_wvls_real(15)        ! phase function wavelengths (for submission to driver) DRF
    real*4 :: phase_angles_real(50)        ! phase function angles (degrees) (for submission to driver) DRF
    real*4 :: z_cld_val_real(pver-1)        ! cloud altitude values (for submission to driver) DRF
    real*4 :: dv_input_real               ! output spectral resolution (for submission to driver DRF
    real*4 :: psipo_real                  ! solar azimuth angle that modtran requires (DRF)
    real*4 :: declination_angle           ! solar declination angle for azimuth calculation (DRF)
    real*4 :: elevation_angle             ! solar elevation angle for azimuth calculation (DRF)
    real*4 :: pi_real                     ! value of pi for azimuth calculation (DRF)
    real*4 :: fuspecalb(15)               !fu-liou spectral albedo
    real*4 :: fu_wvl(15)               !fu-liou spectral wavelengths
    real*4 :: total_od                 !cloud+aerosol optical depth

    !For modis asdir,asdif, aldir,aldif
    real(r8) :: brdf_mult_vis(7) !multiplying factor for wavelength-integration over visible (DRF)
    real(r8) :: brdf_mult_nir(7) !multiplying factor for wavelength-integration over near-ir (DRF)
    real(r8) :: spec_alb_ws(7)  !white-sky spectral albedo (DRF)
    real(r8) :: spec_alb_bs(7)  !black-sky spectral albedo (DRF)
    real(r8) :: gws_iso  ! multiplication factor for white-sky albedo isotropic factor (DRF)
    real(r8) :: gws_vol  ! multiplication factor for white-sky albedo volumetric factor (DRF)
    real(r8) :: gws_geo  ! multiplication factor for white-sky albedo geometric factor (DRF)
    real(r8) :: gbs_iso(3)  ! multiplication factor for black-sky albedo isotropic factor (DRF)
    real(r8) :: gbs_vol(3)  ! multiplication factor for black-sky albedo volumetric factor (DRF)
    real(r8) :: gbs_geo(3)  ! multiplication factor for black-sky albedo geometric factor (DRF)
    real(r8) :: solar_zen_rad
    real(r8) :: modis_asdir_local(pcols) !copy of modis_asdir, modified with brdf vals (DRF)
    real(r8) :: modis_asdif_local(pcols) !copy of modis_asdif, modified with brdf vals (DRF)
    real(r8) :: modis_aldir_local(pcols) !copy of modis_aldir, modified with brdf vals (DRF)
    real(r8) :: modis_aldif_local(pcols) !copy of modis_aldif, modified with brdf vals (DRF)

    real(r8) :: modis_asdir_ocean !copy of modis_asdir, modified with brdf vals (DRF)
    real(r8) :: modis_asdif_ocean !copy of modis_asdif, modified with brdf vals (DRF)
    real(r8) :: modis_aldir_ocean !copy of modis_aldir, modified with brdf vals (DRF)
    real(r8) :: modis_aldif_ocean !copy of modis_aldif, modified with brdf vals (DRF)

!
! Output solar arguments
!
   real(r8), intent(out) :: fsns(pcols)          ! Surface absorbed solar flux
   real(r8), intent(out) :: fsnt(pcols)          ! Net column abs solar flux at model top
   real(r8), intent(out) :: flns(pcols)          ! Srf longwave cooling (up-down) flux
   real(r8), intent(out) :: flnt(pcols)          ! Net outgoing lw flux at model top
   real(r8), intent(out) :: sols(pcols)          ! Downward solar rad onto surface (sw direct)
   real(r8), intent(out) :: soll(pcols)          ! Downward solar rad onto surface (lw direct)
   real(r8), intent(out) :: solsd(pcols)         ! Downward solar rad onto surface (sw diffuse)
   real(r8), intent(out) :: solld(pcols)         ! Downward solar rad onto surface (lw diffuse)
   real(r8), intent(out) :: qrs(pcols,pver)      ! Solar heating rate
   real(r8), intent(out) :: fsds(pcols)          ! Flux Shortwave Downwelling Surface
   real(r8), intent(out) :: aertauout(pcols,naer_groups)   ! aerosol optical depth from radcswmx.F90 (DRF)
   real(r8), intent(out) :: fsnirtoa(pcols)      ! NIR TOA flux all-sky (DRF)
   real(r8), intent(out) :: fsnirtoac(pcols)      ! NIR TOA flux clear-sky (DRF)

#ifdef OFFLINE
   real(r8), intent(out) :: fln(pcols,pverp)     ! Net LW fluxes at interfaces
   real(r8), intent(out) :: flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
   real(r8), intent(out) :: flntc(pcols)         ! Clear sky lw flux at model top!
   real(r8), intent(out) :: fsdsc(pcols)         ! Clear sky surface downwelling solar flux
   real(r8), intent(out) :: fsnsc(pcols)         ! Clear sky surface abs solar flux
   real(r8), intent(out) :: fsn(pcols,pverp)     ! Net SW fluxes at interfaces
   real(r8), intent(out) :: fsntc(pcols)         ! Clear sky total column abs solar flux

   real(r8), intent(out) :: solin(pcols)         ! Solar incident flux!-WDC
#endif
!
! Output longwave arguments
!
   real(r8), intent(out) :: qrl(pcols,pver)      ! Longwave cooling rate
   real(r8), intent(out) :: flwds(pcols)         ! Surface down longwave flux
   real*4, intent(out) :: wavelength_lres(wvlng2)  ! Wavelength from Modtran (DRF)
   real*4, intent(out) :: wavelength_hres(wvlng_hres2)  ! Wavelength from Modtran (DRF)
   real*4, intent(out) :: radiance_lres_clr(pcols,wvlng2)   ! Radiance from Modtran (DRF)
   real*4, intent(out) :: radiance_lres_all(pcols,wvlng2)   ! Radiance from Modtran (DRF)
   real*4, intent(out) :: radiance_hres_clr(pcols,wvlng_hres2)   ! Radiance from Modtran (DRF)
   real*4, intent(out) :: radiance_hres_all(pcols,wvlng_hres2)   ! Radiance from Modtran (DRF)
   real*4, intent(out) :: solar_flux(pcols,wvlng2)   ! TOA solar flux from Modtran (DRF)
   real*4, intent(out) :: diffuse_flux_clr(pcols,wvlng2)   !SW diffuse flux from Modtran (DRF)
   real*4, intent(out) :: diffuse_flux_all(pcols,wvlng2)   !SW diffuse flux from Modtran (DRF)
   real*4, intent(out) :: solar_zenith(pcols)   !Solar zenith angle from Modtran (DRF)
   real*4, intent(out) :: aod(pcols)   !Solar zenith angle from Modtran (DRF)
   real*4, intent(out) :: tau_sulf(pcols)   !Solar zenith angle from Modtran (DRF)
   real*4, intent(out) :: tau_dust(pcols)   !Solar zenith angle from Modtran (DRF)
   real*4, intent(out) :: tau_carb(pcols)   !Solar zenith angle from Modtran (DRF)
   real*4, intent(out) :: tau_sslt(pcols)   !Solar zenith angle from Modtran (DRF)
   real(r8), intent(out) :: cldtau_ice_sw(pcols) !ice cloud optical depth
   real(r8), intent(out) :: cldtau_liq_sw(pcols) !ice cloud optical depth
   real(r8), intent(out) :: cldtau_lw(pcols) !cloud optical depth for LW calcs
   real*4, intent(out) :: bb_updiffuse_clr(pcols,pver)   !broadband upwelling diffuse flux from MODTRAN (DRF)
   real*4, intent(out) :: bb_updiffuse_all(pcols,pver)   !broadband upwelling diffuse flux from MODTRAN (DRF)
   real*4, intent(out) :: bb_dndiffuse_clr(pcols,pver)   !broadband downwelling diffuse flux from MODTRAN (DRF)
   real*4, intent(out) :: bb_dndiffuse_all(pcols,pver)   !broadband downwelling diffuse flux from MODTRAN (DRF)
   real*4, intent(out) :: bb_dndirect_clr(pcols,pver)   !broadband downwelling direct flux from MODTRAN (DRF)
   real*4, intent(out) :: bb_dndirect_all(pcols,pver)   !broadband downwelling direct flux from MODTRAN (DRF)
 
   real, intent(out) :: flx_updiffuse_clr(pcols,pver,num_bands)   !vis/nir upwelling diffuse clear-sky flux from MODTRAN (DRF)
   real, intent(out) :: flx_updiffuse_all(pcols,pver,num_bands)   !vis/nir upwelling diffuse all-sky flux from MODTRAN (DRF)
   real, intent(out) :: flx_dndiffuse_clr(pcols,pver,num_bands)   !vis/nir downwelling diffuse clear-sky flux from MODTRAN (DRF)
   real, intent(out) :: flx_dndiffuse_all(pcols,pver,num_bands)   !vis/nir downwelling diffuse all-sky flux from MODTRAN (DRF)
   real, intent(out) :: flx_dndirect_clr(pcols,pver,num_bands)   !vis/nir downwelling direct clear-sky flux from MODTRAN (DRF)
   real, intent(out) :: flx_dndirect_all(pcols,pver,num_bands)   !vis/nir downwelling direct all-sky flux from MODTRAN (DRF)

!
!---------------------------Local variables-----------------------------
!
   integer i, j, k, jj, kk,mm              ! index
   real*4 ii
!   integer :: in2o, ich4, if11, if12 ! indexes of gases in constituent array

#ifndef OFFLINE
   real(r8) solin(pcols)         ! Solar incident flux
!  real(r8) fsds(pcols)          ! Flux Shortwave Downwelling Surface
#endif
   real(r8) fsntoa(pcols)        ! Net solar flux at TOA
   real(r8) fsntoac(pcols)       ! Clear sky net solar flux at TOA
   real(r8) fsnirt(pcols)        ! Near-IR flux absorbed at toa
   real(r8) fsnrtc(pcols)        ! Clear sky near-IR flux absorbed at toa
   real(r8) fsnirtsq(pcols)      ! Near-IR flux absorbed at toa >= 0.7 microns
#ifndef OFFLINE
   real(r8) fsntc(pcols)         ! Clear sky total column abs solar flux
   real(r8) fsnsc(pcols)         ! Clear sky surface abs solar flux
   real(r8) fsdsc(pcols)         ! Clear sky surface downwelling solar flux
#endif
   real(r8) flut(pcols)          ! Upward flux at top of model
   real(r8) lwcf(pcols)          ! longwave cloud forcing
   real(r8) swcf(pcols)          ! shortwave cloud forcing
   real(r8) flutc(pcols)         ! Upward Clear Sky flux at top of model
#ifndef OFFLINE
   real(r8) flntc(pcols)         ! Clear sky lw flux at model top
   real(r8) flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
#endif
   real(r8) ftem(pcols,pver)     ! temporary array for outfld
   real(r8) fln200(pcols)        ! net longwave flux interpolated to 200 mb
   real(r8) fln200c(pcols)       ! net clearsky longwave flux interpolated to 200 mb
#ifndef OFFLINE
   real(r8) fsn(pcols,pverp)     ! net shortwave flux
#endif
   real(r8) fsnc(pcols,pverp)    ! net clear-sky shortwave flux
   real(r8) fsn200(pcols)        ! fns interpolated to 200 mb
   real(r8) fsn200c(pcols)       ! fsnc interpolated to 200 mb
#ifndef OFFLINE
   real(r8) fln(pcols,pverp)     ! net longwave flux
#endif
   real(r8) flnc(pcols,pverp)    ! net clear-sky longwave flux

!DRF   real(r8) pbr(pcols,pverr)     ! Model mid-level pressures (dynes/cm2)
!DRF   real(r8) pnm(pcols,pverrp)    ! Model interface pressures (dynes/cm2)
#ifndef OFFLINE
   real(r8) o3vmr(pcols,pverr)   ! Ozone volume mixing ratio
#endif
!DRF   real(r8) o3mmr(pcols,pverr)   ! Ozone mass mixing ratio
!DRF   real(r8) eccf                 ! Earth/sun distance factor
   real(r8) n2o(pcols,pver)      ! nitrous oxide mass mixing ratio
   real(r8) ch4(pcols,pver)      ! methane mass mixing ratio
   real(r8) cfc11(pcols,pver)    ! cfc11 mass mixing ratio
   real(r8) cfc12(pcols,pver)    ! cfc12 mass mixing ratio
!DRF   real(r8) rh(pcols,pverr)      ! level relative humidity (fraction)
   real(r8) lwupcgs(pcols)       ! Upward longwave flux in cgs units

!DRF   real(r8) esat(pcols,pverr)    ! saturation vapor pressure
!DRF   real(r8) qsat(pcols,pverr)    ! saturation specific humidity

   real(r8) :: frc_day(pcols) ! = 1 for daylight, =0 for night colums
   real(r8) :: aertau(pcols,nspint,naer_groups) ! Aerosol column optical depth
   real(r8) :: aerssa(pcols,nspint,naer_groups) ! Aerosol column averaged single scattering albedo
   real(r8) :: aerasm(pcols,nspint,naer_groups) ! Aerosol column averaged asymmetry parameter
   real(r8) :: aerfwd(pcols,nspint,naer_groups) ! Aerosol column averaged forward scattering

   real(r8) aerosol(pcols, pver, naer_all) ! aerosol mass mixing ratios
   real(r8) scales(naer_all)               ! scaling factors for aerosols
   integer time_array_0(8)
   integer time_array_1(8)
   integer rate,start_time,end_time

   real*4 m_aerosol(pver,4)             !profile of total aerosol mass (in km^-1 per layer)
   real*4 o_aerosol(42,12)           !absorption, extinction, ssa spectral data for aerosols
   real*4 min_cld_val                !min val for cloud determination
   real(r8) :: m_cloud_liq(pver)              !profile of liquid water cloud DRF
   real*4 :: m_cloud_liq_real(pver)              !profile of liquid water cloud DRF
   real(r8) :: m_cloud_ice(pver)              !profile of ice water cloud
   real*4 :: m_cloud_ice_real(pver)              !profile of liquid water cloud DRF
   real(r8) flip_cldfrac(pver)             !profile of cloud fraction
   logical a_flag,c_flag,dummy_flag,debug_flag,modtran_flag         !flags for aerosols and clouds
   character(len=120) :: fname !to be sent to modtran
   character(len=120) :: fname_cld !to be sent to modtran
   character* (1) fval_1  !DRF for modtran temp file naming
   character* (2) fval_2  !ditto
   character* (3) fval_3  !ditto
   character* (1) filler_period !ditto
   real tau_rayleigh   !DRF rayleigh scattering optical depth at 0.55 um
   real(r8) aerosol_mult !DRF scaling for aerosols
   real*4 dummy_flx(pver,wvlng2)
   character tp5_file(1000)*110
   logical tp5_flag
 
   !Variables for trapezoidal rule
   real*4 dummy_trapval   !variable used to calculate trapezoidal rule
   integer vis_nir_index  !separation index for vis/nir (700 nm)
   real*4 wvl_vis_nir     !separation value for vis/nir (700 nm)
   integer end_nir_index  !max index greater than 0 

   real(r8), allocatable :: out_cldfrac(:,:,:) !DRF for cloud overlap scheme
   real(r8), allocatable :: out_cldliqprof(:,:,:) !DRF for cloud overlap scheme
   real(r8), allocatable :: out_cldiceprof(:,:,:) !DRF for cloud overlap scheme
   
   real, allocatable :: tmp_radiance_lres(:,:)  !Temporary arrays for cloud-overlap approximation averaging
   real, allocatable :: tmp_radiance_hres(:,:)
   real, allocatable :: tmp_diffuse_flux(:,:)
   real, allocatable :: tmp_bb_updiffuse(:,:)
   real, allocatable :: tmp_bb_dndiffuse(:,:)
   real, allocatable :: tmp_bb_dndirect(:,:)
   real, allocatable :: tmp_flx_updiffuse(:,:,:)
   real, allocatable :: tmp_flx_dndiffuse(:,:,:)
   real, allocatable :: tmp_flx_dndirect(:,:,:)
   real ocean_alb  !for debugging

   !Variables for cloud optical properties
   real(r8) dummy_ext(24),dummy_ssa(24),dummy_asym(24),dummy_sum
   real(r8) wvl_cld(24),ext_liq_cld(24,pver),ssa_liq_cld(24,pver),asym_liq_cld(24,pver)
   real(r8) ext_ice_cld(24,pver),ssa_ice_cld(24,pver),asym_ice_cld(24,pver)
   real(r8) liq_ext_avg(24),liq_ssa_avg(24),liq_asym_avg(24)
   real(r8) ice_ext_avg(24),ice_ssa_avg(24),ice_asym_avg(24)

   !for cloud overlap approximation
   integer num_profs !DRF Number of times to call cloud overlap maximum-random overlap profile generator
   integer :: iscloudy_matrix(16,pver) !0 if no cloud, 1 if cloud, set to num_profs
   integer :: cloud_overlap_approx
   integer :: changeSeed 
   real(r8) :: press_levs(pver,pcols)
   integer cloud_prof_id(16)  !array of values ranging from 1 to 2^16
   integer weights(16,pcols) 
   logical :: calc_or_nocalc(16,pcols)
   integer :: num_modtran_calls_per_grid(pcols)

! 
! Maximum number of configurations to include in solution
! 
   integer nconfgmax
   parameter (nconfgmax = 15)
! 
! Cloud configurations
! 
   real(r8) wgtv_v(pcols,nconfgmax)   ! Fractional area for cloud configurations
   real(r8) totwgt_v(pcols)           ! Sum of wgtv_v = total fractional area
   integer ccon_v(pcols,pver,nconfgmax)!Binary cloud configurations (0=no,1=yes)
   integer nconfig_v(pcols)           ! # of cloud configurations in CAM columns
   logical cldovrlp


!DRF
   !set to false if going to use aerosols
   modtran_flag = .true.
   debug_flag = .false.
   no_aerosol_flag = .false.  !doesn't work, must be set to false
   
   no_sulf_flag = .false.
   no_dust_flag = .false.
   no_carb_flag = .false.
   no_sslt_flag = .false.
   !if (debug_flag) then
      no_sulf_flag = .true.
      no_dust_flag = .true.
      no_carb_flag = .true.
      no_sslt_flag = .true.
   !endif 

   if (no_aerosol_flag) then
      do i=1,ncol
         do j=1,pver
            do k=1,12
               aerosol(i,j,k) = 0.
            end do   
         end do
      end do
   endif

   if (no_sulf_flag) then
      do j=1,ncol
        do kk=1,pver
           aerosol(j,kk,1) = 0.
        end do
      end do
   endif 

   if (no_dust_flag) then
      do j=1,ncol
        do kk=1,pver
           aerosol(j,kk,idxDUSTfirst) = 0.
           aerosol(j,kk,idxDUSTfirst+1) = 0.
           aerosol(j,kk,idxDUSTfirst+2) = 0.
           aerosol(j,kk,idxDUSTfirst+3) = 0.
        end do
      end do
   endif 
 
   if (no_carb_flag) then
      do j=1,ncol
        do kk=1,pver
           aerosol(j,kk,7) = 0.
           aerosol(j,kk,8) = 0.
           aerosol(j,kk,9) = 0.
           aerosol(j,kk,10) = 0.
        end do
      end do
   endif

   if (no_sslt_flag) then
      do j=1,ncol
        do kk=1,pver
           aerosol(j,kk,2) = 0.
        end do
      end do
   endif

!!! one-layer aerosol
!   do j=1,ncol
!     do kk=1,pver
!        aerosol(j,kk,2) = 0.
!     end do
!     aerosol(j,20,2) = 0.001
!   end do

   !!!


   do i=1,pver
     do j=1,ncol
        if (rei(j,i).gt.0) then
           rei_local(j,i) = rei(j,i)
        else
           rei_local(j,i) = 5.
        endif
     end do
   end do

   do i=1,1000
     rh_vec(i) = 1.0_r8 / 1000 * (i-1)
   end do

   num_profs = 16
   cloud_overlap_approx = 2  !maximum-random overlap in generate_stochastic clouds
   changeSeed = 0
   gravit = 9.80616
   clr_flx_index = 1
   cthik_input = -9.
   cext_input = -9.
!

!
! Solar radiation computation
!

!write(*,*) "PNM at start = ",pnm

   if (dosw) then

!
! calculate heating with aerosols
!


      if (radforce) then

         pmxrgnrf = pmxrgn
         nmxrgnrf = nmxrgn

#ifndef OFFLINE
         call get_rf_scales(scales)

         call get_aerosol(lchnk, pint, aerosol, scales)

         ! overwrite with prognostics aerosols
         call set_aerosol_from_prognostics (state, aerosol)

         call aerosol_indirect(ncol,lchnk,landfrac,pmid,t,qm1,cld,zm,rel)
#endif   
         call t_startf('radcswmx_rf')


!         call radcswmx(lchnk   ,ncol ,                            &
!                    pnm     ,pbr     ,qm1     ,rh      ,o3mmr   , &
!                    aerosol ,cld     ,cicewp  ,cliqwp  ,rel     , &
!                    rei     ,eccf    ,coszrs  ,scon    ,solin   , &
!                    asdir   ,asdif   ,aldir   ,aldif   ,nmxrgnrf, &
!                    pmxrgnrf,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
!                    fsntoac ,fsnirt  ,fsnrtc  ,fsnirtsq,fsns    , &
!                    fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
!                    solsd   ,solld   ,frc_day ,                   &
!                    aertau  ,aerssa  ,aerasm  ,aerfwd  ,fsn     , &
!#ifndef OFFLINE
!                    fsnc    )
!#else
!                    fsnc    ,no_o2_abs)
!#endif

         call t_stopf('radcswmx_rf')

!
! Convert units of shortwave fields needed by rest of model from CGS to MKS
!

            do i = 1, ncol
            solin(i) = solin(i)*1.e-3
            fsnt(i)  = fsnt(i) *1.e-3
            fsns(i)  = fsns(i) *1.e-3
            fsntc(i) = fsntc(i)*1.e-3
            fsnsc(i) = fsnsc(i)*1.e-3
         end do
         ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair

!
! Dump shortwave radiation information to history tape buffer (diagnostics)
!
         call outfld('QRS_RF  ',ftem  ,pcols,lchnk)
         call outfld('FSNT_RF ',fsnt  ,pcols,lchnk)
         call outfld('FSNS_RF ',fsns  ,pcols,lchnk)
         call outfld('FSNTC_RF',fsntc ,pcols,lchnk)
         call outfld('FSNSC_RF',fsnsc ,pcols,lchnk)
 
      endif ! if (radforce)

#ifndef OFFLINE
      call get_int_scales(scales)
 
      call get_aerosol(lchnk, pint, aerosol, scales)

      ! overwrite with prognostics aerosols
      call set_aerosol_from_prognostics (state, aerosol)

      call aerosol_indirect(ncol,lchnk,landfrac,pmid,t,qm1,cld,zm,rel)
#endif

      call t_startf('radcswmx')

!DRF variables
      write(*,*) 'assigning cliqwp',ncol, pver

      !Test different cloud configurations
      do i=1,ncol
         do k=1,pver
           cliqwp_local(i,k) = cliqwp(i,k)
           cicewp_local(i,k) = cicewp(i,k)
           cld_local(i,k) = cld(i,k)
           emis_local(i,k) = emis(i,k)
           t_local(i,k) = t(i,k)
           h2o_local(i,k) = qm1(i,k,1)
           o3_local(i,k) = o3vmr(i,k)
           n2o_local(i,k) = qm1(i,k,in2o)
           ch4_local(i,k) = qm1(i,k,ich4)
           cfc11_local(i,k) = qm1(i,k,if11)
           cfc12_local(i,k) = qm1(i,k,if12)

           !if (maxval(t_local(i,:)).lt.400.0_r8 .and. maxval(h2o_local(i,:)).lt.400.0_r8) then
           !    lowest_t = t(i,k)
           !endif
!BEGIN 3/28/16 test
           !if (t_local(i,k).lt.400.0_r8 .and.h2o_local(i,k).lt.400.0_r8 .and. maxval(t_local(i,:)).lt.400.0_r8 .and. maxval(h2o_local(i,:)).lt.400.0_r8) then
           !   lowest_t = t(i,k)
           !   lowest_h2o = h2o_local(i,k)
           !   lowest_o3 = o3_local(i,k)
           !   lowest_ch4 = ch4_local(i,k)
           !   lowest_n2o = n2o_local(i,k)
           !   lowest_cfc11 = cfc11_local(i,k)
           !   lowest_cfc12 = cfc12_local(i,k)
           !else
           !   !write(*,*) 'replacing t',lowest_t,h2o_local(i,1),t_local(i,1:k-1)
           !   t_local(i,k) = lowest_t
           !   h2o_local(i,k) = lowest_h2o
           !   o3_local(i,k) = lowest_o3 ! 0.0000000001_r8    !0.0_r8
           !   n2o_local(i,k) = lowest_n2o !0.0000000001_r8   !0.0_r8
           !   ch4_local(i,k) = lowest_ch4 !0.0000000001_r8   !0.0_r8
           !   cfc11_local(i,k) = lowest_cfc11 !0.0000000001_r8 !0.0_r8
           !   cfc12_local(i,k) = lowest_cfc12 !0.0000000001_r8 !0.0_r8
           !endif

           !if (maxval(t_local(i,:)).gt.400.0_r8 .or. maxval(h2o_local(i,:)).gt.400.0_r8) then
           !   do kk=1,pver
           !      if (t_local(i,kk).gt.400.0_r8) then
           !          t_local(i,kk) = lowest_t
           !          h2o_local(i,kk) = lowest_h2o    !0.0_r8
           !          o3_local(i,kk) = lowest_o3 !0.0000000001_r8    !0.0_r8
           !          n2o_local(i,kk) = lowest_n2o !0.0000000001_r8   !0.0_r8
           !          ch4_local(i,kk) = lowest_ch4 !0.0000000001_r8   !0.0_r8
           !          cfc11_local(i,kk) = lowest_cfc11 !0.0000000001_r8 !0.0_r8
           !          cfc12_local(i,kk) = lowest_cfc12 !0.0000000001_r8 !0.0_r8
           !       endif
           !  end do
           !endif
!END 3/28/16 test

           !if (i.eq.144 .and. lchnk.eq.144) then
           !    write(*,*) 't_local at 144 = ',t_local(i,:)
           !    write(*,*) 'h2o_local at 144 = ',h2o_local(i,:)
           !    write(*,*) 'lowest_t = ',lowest_t
           !endif

     if (i.eq.166 .and. lchnk.eq.166) then
     !write(*,*) 'brdf_param1 = ',brdf_param(i,:,1)
     !write(*,*) 'brdf_param2 = ',brdf_param(i,:,2)
     !write(*,*) 'brdf_param3 = ',brdf_param(i,:,3)
     !write(*,*) 'land_flag = ',land_flag(:)
     !write(*,*) 'landfrac_166 = ',land_flag(166)
     !write(*,*) 't = ',t_local(1,:)
     !write(*,*) 'qm1 = ',qm1(1,1,1)
     !write(*,*) 'o3vmr = ',o3vmr(1,:)
     !write(*,*) 'cld = ',cld(1,:)
     !write(*,*) 'emis = ',emis_local(i,:)
     !write(*,*) 'rei_local = ',rei_local(i,:)
     !write(*,*) 'cliqwp_local = ',cliqwp_local(i,:)
     !write(*,*) 'cicewp_local = ',cicewp_local(i,:)
     endif


           !write(*,*) 'max cliqwp = ',maxval(cliqwp_local)
           !write(*,*) 'min cliqwp = ',maxval(cliqwp_local)

           !write(*,*) 'max cicewp = ',maxval(cicewp_local)
           !write(*,*) 'min cicewp = ',maxval(cicewp_local)


           if (isnan(cliqwp_local(i,k))) then
              cliqwp_local(i,k) = 0.
              write(*,*) 'nan cliqwp',cliqwp(i,k),i,k
           endif  
           if (isnan(cicewp_local(i,k))) then
              cicewp_local(i,k) = 0.
              write(*,*) 'nan cicewp',cicewp(i,k),i,k
           endif  
       
          if (cicewp_local(i,k).lt.0.0_r8 .or. cliqwp_local(i,k).lt.0.0_r8) then
             write(*,*) 'liq cloud lt zero', i, k, cicewp_local(i,k), cliqwp_local(i,k)
             if (cicewp_local(i,k).lt.0.0_r8) then
                cicewp_local(i,k) = 0.
             endif
             if (cliqwp_local(i,k).lt.0.0_r8) then
                cliqwp_local(i,k) = 0.
             endif
          endif

         end do 

         !cliqwp_local(i,22) = 400.  !test liquid clouds
         !cld_local(i,22) = 1.0     !test liquid clouds
         !cliqwp_local(i,21) = 0.01  !test liquid clouds
         !cld_local(i,21) = 1.0     !test liquid clouds
         !t_local(i,21) = t_local(i,22) !make layers isothermal
         !t_local(i,20) = t_local(i,22) !make layers isothermal

         !emis_local(i,22) = 1. - exp(-1.66*0.090361*cliqwp_local(i,22))
         !emis_local(i,21) = 1. - exp(-1.66*0.090361*cliqwp_local(i,21))

         !parameter (kabsl = 0.090361)
         ! kabsi = 0.005 + 1./rei(i,k)
         ! kabs = kabsl*(1.-fice(i,k)) + kabsi*fice(i,k)
         ! emis(i,k) = 1. - exp(-1.66*kabs*clwp(i,k))
         
         !cicewp_local(i,15) = 10.  !test ice clouds
         !cld_local(i,15) = 1.0     !test ice clouds
         !cicewp_local(i,14) = 0.01  !test ice clouds
         !cld_local(i,14) = 1.0     !test ice clouds

         !Calculate cloud optical depth in the LW
         dummy_sum = 0.
         do k=1,pver
            if (rei_local(i,k).gt.0.) then
              dummy_sum = dummy_sum + 1.66*(0.090361*cliqwp_local(i,k) + (0.005+1./rei_local(i,k))*cicewp_local(i,k))
            else
              dummy_sum = dummy_sum + 1.66*(0.090361*cliqwp_local(i,k) + (0.005+1./41.)*cicewp_local(i,k))
            endif
         end do
         cldtau_lw(i) = dummy_sum
      end do

      !modis albedo modifications

      !for white-sky albedo integration
      gws_iso = 1.
      gws_vol = 0.189184
      gws_geo = -1.377622

      !for black-sky albedo integration
      gbs_iso(1) = 1. 
      gbs_iso(2) = 0.
      gbs_iso(3) = 0.
      gbs_vol(1) = -0.007574
      gbs_vol(2) = -0.070987
      gbs_vol(3) = 0.307588
      gbs_geo(1) = -1.284909
      gbs_geo(2) = -0.166314
      gbs_geo(3) = 0.041840

      !brdf_mult_vis(1) = 0.4364
      !brdf_mult_vis(2) = 0.2366
      !brdf_mult_vis(3) = 0.3265 
      !brdf_mult_vis(1) = 0.5096
      !brdf_mult_vis(2) = 0.2437
      !brdf_mult_vis(3) = 0.2465


      !brdf_mult_nir(4) =  0.5447
      !brdf_mult_nir(5) =  0.1363
      !brdf_mult_nir(6) =  0.0469
      !brdf_mult_nir(7) =  0.2536

      !brdf_mult_nir(4) =  0.5318
      !brdf_mult_nir(5) =  0.2790
      !brdf_mult_nir(6) =  0.1539
      !brdf_mult_nir(7) =  0.1446

      !   brdf_mult_nir(4) =  0.5271
      !   brdf_mult_nir(5) =  0.1795
      !   brdf_mult_nir(6) =  0.0000
      !   brdf_mult_nir(7) =  0.2755


      do i=1,ncol
         modis_asdir_local(i) = modis_asdir(i)
         modis_asdif_local(i) = modis_asdif(i)
         modis_aldir_local(i) = modis_aldir(i)
         modis_aldif_local(i) = modis_aldif(i)

 !       if (j.eq.210 .and. lchnk.eq.132) then
	 !write(*,*) 'asdir = ',maxval(modis_asdir),minval(modis_asdir)
	 !write(*,*) 'asdif = ',maxval(modis_asdif),minval(modis_asdif)
	 !write(*,*) 'aldir = ',maxval(modis_asdir),minval(modis_asdir)
	 !write(*,*) 'aldif = ',maxval(modis_asdif),minval(modis_asdif)


         dummy_flag = .false.
         solar_zen_rad = acos(coszrs(i))

         if (land_flag(i)) then
            !write(*,*) 'land, either snowy or snow-free'

            brdf_mult_vis(1) = (1.-fsno_read(i))*0.4364+fsno_read(i)*0.5096
            brdf_mult_vis(2) = (1.-fsno_read(i))*0.2366+fsno_read(i)*0.2437
            brdf_mult_vis(3) = (1.-fsno_read(i))*0.3265+fsno_read(i)*0.2465
            !brdf_mult_vis(1) = 0.5096
            !brdf_mult_vis(2) = 0.2437
            !brdf_mult_vis(3) = 0.2465
            brdf_mult_vis(4) = 0.
            brdf_mult_vis(5) = 0.
            brdf_mult_vis(6) = 0.
            brdf_mult_vis(7) = 0.

            brdf_mult_nir(1) = 0.
            brdf_mult_nir(2) = 0.
            brdf_mult_nir(3) = 0.
            brdf_mult_nir(4) = (1.-fsno_read(i))*0.5271+fsno_read(i)*0.5318
            brdf_mult_nir(5) = (1.-fsno_read(i))*0.1795+fsno_read(i)*0.2790
            brdf_mult_nir(6) = (1.-fsno_read(i))*0.0 + fsno_read(i)*0.1539
            brdf_mult_nir(7) = (1.-fsno_read(i))*0.2755+fsno_read(i)*0.1446
            !brdf_mult_nir(4) =  0.5318
            !brdf_mult_nir(5) =  0.2790
            !brdf_mult_nir(6) =  0.1539
            !brdf_mult_nir(7) =  0.1446

            allocate(brdf_param_dble(7,3))
            do ii=1,7
               brdf_param_dble(ii,1) = brdf_param(i,ii,1)/1000.
               brdf_param_dble(ii,2) = brdf_param(i,ii,3)/1000.
               brdf_param_dble(ii,3) = brdf_param(i,ii,2)/1000.
            end do

            do ii=1,7
               spec_alb_ws(ii) = 0.
               spec_alb_ws(ii) = spec_alb_ws(ii)+gws_iso*brdf_param_dble(ii,1)
               spec_alb_ws(ii) = spec_alb_ws(ii)+gws_vol*brdf_param_dble(ii,2)
               spec_alb_ws(ii) = spec_alb_ws(ii)+gws_geo*brdf_param_dble(ii,3)
   
               spec_alb_bs(ii) = 0.
               spec_alb_bs(ii) = spec_alb_bs(ii) + brdf_param_dble(ii,1)* &
                    (gbs_iso(1) + gbs_iso(2)*solar_zen_rad**2 + gbs_iso(3)*solar_zen_rad**3)
               spec_alb_bs(ii) = spec_alb_bs(ii)+ brdf_param_dble(ii,2)* &
                    (gbs_vol(1) + gbs_vol(2)*solar_zen_rad**2 + gbs_vol(3)*solar_zen_rad**3)
               spec_alb_bs(ii) = spec_alb_bs(ii)+ brdf_param_dble(ii,3)* &
                    (gbs_geo(1) + gbs_geo(2)*solar_zen_rad**2 + gbs_geo(3)*solar_zen_rad**3)
            end do
            modis_asdir_local(i) = 0.
            modis_asdif_local(i) = 0.
            modis_aldir_local(i) = 0.
            modis_aldif_local(i) = 0.

            do ii=1,7
               modis_asdir_local(i) = modis_asdir_local(i) + spec_alb_bs(ii)*brdf_mult_vis(ii) 
               modis_asdif_local(i) = modis_asdif_local(i) + spec_alb_ws(ii)*brdf_mult_vis(ii) 
               modis_aldir_local(i) = modis_aldir_local(i) + spec_alb_bs(ii)*brdf_mult_nir(ii) 
               modis_aldif_local(i) = modis_aldif_local(i) + spec_alb_ws(ii)*brdf_mult_nir(ii) 
            end do

            if (debug_flag) then
                write(*,*) 'i = ',i
                write(*,*) 'solar_zen_rad = ',solar_zen_rad
                write(*,*) 'brdf_param2_1 = ',brdf_param_dble(:,1) 
                write(*,*) 'brdf_param2_2 = ',brdf_param_dble(:,2) 
                write(*,*) 'brdf_param2_3 = ',brdf_param_dble(:,3) 
                write(*,*) 'modis_asdir = ',modis_asdir_local(i)
                write(*,*) 'modis_asdif = ',modis_asdif_local(i)
                write(*,*) 'modis_aldir = ',modis_aldir_local(i)
                write(*,*) 'modis_aldif = ',modis_aldif_local(i)
            endif

            deallocate(brdf_param_dble)
         endif

         ice_flag_local = .false. 
         if (icefrac(i).ge.0.001) then
            !write(*,*) 'icy'
            ice_flag_local = .true.

            allocate(brdf_param_dble(7,3))
            allocate(brdf_wvl_dble(7))

            brdf_mult_vis(1) = 0.5096
            brdf_mult_vis(2) = 0.2437
            brdf_mult_vis(3) = 0.2465
            brdf_mult_vis(4) = 0.
            brdf_mult_vis(5) = 0.
            brdf_mult_vis(6) = 0.
            brdf_mult_vis(7) = 0.

            brdf_mult_nir(1) = 0.
            brdf_mult_nir(2) = 0.
            brdf_mult_nir(3) = 0.
            brdf_mult_nir(4) = 0.5318
            brdf_mult_nir(5) = 0.2790
            brdf_mult_nir(6) = 0.1539
            brdf_mult_nir(7) = 0.1446

            !integrate up ice-brdf
            brdf_param_dble(1,1) = 0.76
            brdf_param_dble(2,1) = 0.74
            brdf_param_dble(3,1) = 0.72
            brdf_param_dble(4,1) = 0.6
            brdf_param_dble(5,1) = 0.25
            brdf_param_dble(6,1) = 0.07
            brdf_param_dble(7,1) = 0.06

            do jj=1,7
               brdf_param_dble(jj,2) = 0. !vol
               brdf_param_dble(jj,3) = 0. !geo 
            end do

            do ii=1,7
                spec_alb_ws(ii) = 0.
                spec_alb_ws(ii) = spec_alb_ws(ii)+gws_iso*brdf_param_dble(ii,1)
                spec_alb_ws(ii) = spec_alb_ws(ii)+gws_vol*brdf_param_dble(ii,2)
                spec_alb_ws(ii) = spec_alb_ws(ii)+gws_geo*brdf_param_dble(ii,3)
   
                spec_alb_bs(ii) = 0.
                spec_alb_bs(ii) = spec_alb_bs(ii)+ brdf_param_dble(ii,1)* &
                     (gbs_iso(1) + gbs_iso(2)*solar_zen_rad**2 + gbs_iso(3)*solar_zen_rad**3)
                spec_alb_bs(ii) = spec_alb_bs(ii)+ brdf_param_dble(ii,2)* &
                     (gbs_vol(1) + gbs_vol(2)*solar_zen_rad**2 + gbs_vol(3)*solar_zen_rad**3)
                spec_alb_bs(ii) = spec_alb_bs(ii)+ brdf_param_dble(ii,3)* &
                     (gbs_geo(1) + gbs_geo(2)*solar_zen_rad**2 + gbs_geo(3)*solar_zen_rad**3)
            end do
            modis_asdir_local(i) = 0.
            modis_asdif_local(i) = 0.
            modis_aldir_local(i) = 0.
            modis_aldif_local(i) = 0.

            do ii=1,7
               modis_asdir_local(i) = modis_asdir_local(i) + spec_alb_bs(ii)*brdf_mult_vis(ii) 
               modis_asdif_local(i) = modis_asdif_local(i) + spec_alb_ws(ii)*brdf_mult_vis(ii) 
               modis_aldir_local(i) = modis_aldir_local(i) + spec_alb_bs(ii)*brdf_mult_nir(ii) 
               modis_aldif_local(i) = modis_aldif_local(i) + spec_alb_ws(ii)*brdf_mult_nir(ii) 
            end do

            !Recombine with ocean asdir,aldir,asdif,aldif and scale by icefrac
            call albocean(coszrs(i),modis_asdir_ocean,modis_aldir_ocean, &
                       modis_asdif_ocean,modis_aldif_ocean)
         
            modis_asdir_local(i) = icefrac(i)*modis_asdir_local(i)+(1.-icefrac(i))*modis_asdir_ocean
            modis_asdif_local(i) = icefrac(i)*modis_asdif_local(i)+(1.-icefrac(i))*modis_asdif_ocean
            modis_aldir_local(i) = icefrac(i)*modis_aldir_local(i)+(1.-icefrac(i))*modis_aldir_ocean
            modis_aldif_local(i) = icefrac(i)*modis_aldif_local(i)+(1.-icefrac(i))*modis_aldif_ocean

            if (debug_flag) then
                write(*,*) 'brdf_param2_1 = ',brdf_param_dble(:,1) 
                write(*,*) 'brdf_param2_2 = ',brdf_param_dble(:,2) 
                write(*,*) 'brdf_param2_3 = ',brdf_param_dble(:,3) 
                write(*,*) 'icefrac = ',icefrac(i)
                write(*,*) 'modis_asdir = ',modis_asdir_ocean
                write(*,*) 'modis_asdif = ',modis_asdif_ocean
                write(*,*) 'modis_aldir = ',modis_aldir_ocean
                write(*,*) 'modis_aldif = ',modis_aldif_ocean
                write(*,*) 'modis_asdir = ',modis_asdir_local(i)
                write(*,*) 'modis_asdif = ',modis_asdif_local(i)
                write(*,*) 'modis_aldir = ',modis_aldir_local(i)
                write(*,*) 'modis_aldif = ',modis_aldif_local(i)
                write(*,*) 'solar_zen = ',solar_zen_rad
            endif
            deallocate(brdf_param_dble)
            deallocate(brdf_wvl_dble)
         endif
      end do





      cldovrlp = .true.
      call radcswmx_cldovrlp(lchnk   ,ncol    ,                            &
                    pnm     ,pbr     ,qm1     ,rh      ,o3mmr   , &
                    aerosol ,cld_local ,cicewp_local  ,cliqwp_local  ,rel     , &
                    rei_local ,eccf    ,coszrs  ,scon    ,solin   , &
                    modis_asdir_local,modis_asdif_local,modis_aldir_local,modis_aldif_local,nmxrgn, &
                    pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                    fsntoac ,fsnirt  ,fsnrtc  ,fsnirtsq,fsns    , &
                    fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                    solsd   ,solld   ,frc_day ,                   &
                    aertau  ,aerssa  ,aerasm  ,aerfwd  ,fsn     , &
                    cldtau_ice_sw, cldtau_liq_sw                , &
                    wgtv_v  ,totwgt_v,ccon_v  ,nconfig_v,         &
#ifndef OFFLINE
                    fsnc,cldovrlp    )
#else
                    fsnc    ,no_o2_abs,cldovrlp)
#endif


      call radcswmx(lchnk   ,ncol    ,                            &
                    pnm     ,pbr     ,qm1     ,rh      ,o3mmr   , &
                    aerosol ,cld_local ,cicewp_local  ,cliqwp_local  ,rel     , &
                    rei_local ,eccf    ,coszrs  ,scon    ,solin   , &
                    modis_asdir_local,modis_asdif_local,modis_aldir_local,modis_aldif_local,nmxrgn, &
                    pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                    fsntoac ,fsnirt  ,fsnrtc  ,fsnirtsq,fsns    , &
                    fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                    solsd   ,solld   ,frc_day ,                   &
                    aertau  ,aerssa  ,aerasm  ,aerfwd  ,fsn     , &
                    cldtau_ice_sw, cldtau_liq_sw                , &
                    wgtv_v  ,totwgt_v,ccon_v  ,nconfig_v,         &
#ifndef OFFLINE
                    fsnc    )
#else
                    fsnc    ,no_o2_abs)
#endif

!	if (lchnk.eq.62) then
!		write(*,*) 'ncol = ',ncol
!		write(*,*) 'naer_all = ',naer_all
!		write(*,*) 'pnm = ',pnm(211,:)
!		write(*,*) 'pint = ',pint(211,:)
!		write(*,*) 'ppb = ',pbr(211,:)
!		write(*,*) 'rh in radctl= ',rh(211,:)
!		write(*,*) 'o3mmr in radctl= ',o3mmr(211,:)
!		write(*,*) 'aerosol 1 in radctl= ',aerosol(211,:,1)
!		write(*,*) 'aerosol 2 in radctl= ',aerosol(211,:,2)
!		write(*,*) 'aerosol 3 in radctl= ',aerosol(211,:,3)
!		write(*,*) 'aerosol 4 in radctl= ',aerosol(211,:,4)
!		write(*,*) 'aerosol 5 in radctl= ',aerosol(211,:,5)
!		write(*,*) 'aerosol 6 in radctl= ',aerosol(211,:,6)
!		write(*,*) 'cld_local = ',cld_local(211,:)
!		write(*,*) 'cloud ice = ',cicewp_local(211,:)
!		write(*,*) 'cloud liq = ',cliqwp_local(211,:)
!		write(*,*) 'rel = ',rel(211,:)
!		write(*,*) 'rei = ',rei(211,:)
!		write(*,*) 'h2o = ',qm1(211,:,1)
!		write(*,*) 'o3 = ',o3vmr(211,:)
!		write(*,*) 'n2o = ',qm1(211,:,in2o)
!		write(*,*) 'ch4 = ',qm1(211,:,ich4)
!		write(*,*) 'cfc11 = ',qm1(211,:,if11)
!		write(*,*) 'cfc12 = ',qm1(211,:,if12)
!
!		write(*,*) 'scon = ',scon
!		write(*,*) 'solin = ',solin(211)
!		write(*,*) 'eccf = ',eccf
!		write(*,*) 'doszrs = ',coszrs(17)
!		write(*,*) 'modis_asdir = ',modis_asdir(211)
!		write(*,*) 'modis_asdif = ',modis_asdif(211)
!		write(*,*) 'modis_aldir = ',modis_aldir(211)
!		write(*,*) 'modis_aldif = ',modis_aldif(211)

!		write(*,*) 'pmxrgn = ',pmxrgn(211,:)
!		write(*,*) 'nmxrgn = ',nmxrgn(211)

!		write(*,*) 'fsnt = ',fsnt(211)
!		write(*,*) 'fsntc = ',fsntc(211)
!		stop 'radctl in 1260'
!	endif
        
     do ii=1,pcols
       if (cldtau_liq_sw(ii).lt.1.) then
         !write(*,*) 'cldtau_liq zeros at ',i,latitude_index
         !write(*,*) 'cld liq = ',cliqwp_local(i,:)
       endif
     end do

     !DRF
     !modify cliqwp_local
     !if (lw_flag) then 
     !   do i=1,pcols
     !     do j=1,pver
     !        liq_ice_diff = abs(cliqwp_local(i,j) - cicewp_local(i,j))
     !        liq_ice_tolerance = 0.0001
     !        if (cliqwp_local(i,j).gt.0. .and. cicewp_local(i,j).gt.0. .and. liq_ice_diff.lt.liq_ice_tolerance) then
     !            cliqwp_local(i,j) = 0.
     !        endif
     !     end do
     !   end do
     !endif 
     !DRF

     !DRF for assigning aerosols
     aertauout(:,:) = aertau(:,8,:)
     !write(*,*) 'AFTER RADCSWMX',lchnk

     !DRF


!DRF set qrs, and others to 0.0

      !do i=1,pcols
      !  do  j=1,pver
      !    qrs(i,j) = 0.0
      !    qrl(i,j) = 0.0
      !  end do
      !  sols(i) = 0.0
      !  solsd(i) = 0.0
      !  soll(i) = 0.0
      !  solld(i) = 0.0
      !end do


!DRF

      call t_stopf('radcswmx')


! -- tls ---------------------------------------------------------------2

!  Output net fluxes at 200 mb

      call vertinterp(ncol, pcols, pverp, pint, 20000._r8, fsnc, fsn200c)
      call vertinterp(ncol, pcols, pverp, pint, 20000._r8, fsn, fsn200)

!
! Convert units of shortwave fields needed by rest of model from CGS to MKS
!
      do i=1,ncol
         solin(i) = solin(i)*1.e-3
         fsds(i)  = fsds(i)*1.e-3
         fsnirt(i)= fsnirt(i)*1.e-3
         fsnirtoa(i) = fsnirt(i) !DRF 
         fsnrtc(i)= fsnrtc(i)*1.e-3
         fsnirtoac(i) = fsnrtc(i)  !DRF
         fsnirtsq(i)= fsnirtsq(i)*1.e-3
         fsnt(i)  = fsnt(i) *1.e-3
         fsns(i)  = fsns(i) *1.e-3
         fsntc(i) = fsntc(i)*1.e-3
         fsnsc(i) = fsnsc(i)*1.e-3
         fsdsc(i) = fsdsc(i)*1.e-3
         fsntoa(i)=fsntoa(i)*1.e-3
         fsntoac(i)=fsntoac(i)*1.e-3
         fsn200(i)  = fsn200(i)*1.e-3
         fsn200c(i) = fsn200c(i)*1.e-3
      end do
      ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair

!
! Dump shortwave radiation information to history tape buffer (diagnostics)
!

      call outfld('frc_day ', frc_day, pcols, lchnk)
      call outfld('SULOD_v ', aertau(:,idxVIS,1) ,pcols,lchnk)
      call outfld('SSLTOD_v', aertau(:,idxVIS,2) ,pcols,lchnk)
      call outfld('CAROD_v ', aertau(:,idxVIS,3) ,pcols,lchnk)
      call outfld('DUSTOD_v', aertau(:,idxVIS,4) ,pcols,lchnk)
      call outfld('BGOD_v  ', aertau(:,idxVIS,5) ,pcols,lchnk)
      call outfld('VOLCOD_v', aertau(:,idxVIS,6) ,pcols,lchnk)
      call outfld('AEROD_v ', aertau(:,idxVIS,7) ,pcols,lchnk)
      call outfld('AERSSA_v', aerssa(:,idxVIS,7) ,pcols,lchnk)
      call outfld('AERASM_v', aerasm(:,idxVIS,7) ,pcols,lchnk)
      call outfld('AERFWD_v', aerfwd(:,idxVIS,7) ,pcols,lchnk)
#ifndef OFFLINE
      call aerosol_diagnostics (state, aerosol)
#endif

      call outfld('QRS     ',ftem  ,pcols,lchnk)
      call outfld('SOLIN   ',solin ,pcols,lchnk)
      call outfld('FSDS    ',fsds  ,pcols,lchnk)
      call outfld('FSNIRTOA',fsnirtoa,pcols,lchnk)
      call outfld('FSNRTOAC',fsnirtoac,pcols,lchnk)
      call outfld('FSNRTOAS',fsnirtsq,pcols,lchnk)
      call outfld('FSNT    ',fsnt  ,pcols,lchnk)
      call outfld('FSNS    ',fsns  ,pcols,lchnk)
      call outfld('FSNTC   ',fsntc ,pcols,lchnk)
      call outfld('FSNSC   ',fsnsc ,pcols,lchnk)
      call outfld('FSDSC   ',fsdsc ,pcols,lchnk)
      call outfld('FSNTOA  ',fsntoa,pcols,lchnk)
      call outfld('FSNTOAC ',fsntoac,pcols,lchnk)
      call outfld('SOLS    ',sols  ,pcols,lchnk)
      call outfld('SOLL    ',soll  ,pcols,lchnk)
      call outfld('SOLSD   ',solsd ,pcols,lchnk)
      call outfld('SOLLD   ',solld ,pcols,lchnk)
      call outfld('FSN200  ',fsn200,pcols,lchnk)
      call outfld('FSN200C ',fsn200c,pcols,lchnk)

   end if !dosw
!
! Longwave radiation computation
!
   if (dolw) then
!
! Convert upward longwave flux units to CGS
!
      do i=1,ncol
         !write(*,*) 'ts = ',ts(i)
         lwupcgs(i) = 0.000056704*ts(i)**(4.0_r8)
      end do
!
! Do longwave computation. If not implementing greenhouse gas code then
! first specify trace gas mixing ratios. If greenhouse gas code then:
!  o ixtrcg   => indx of advected n2o tracer
!  o ixtrcg+1 => indx of advected ch4 tracer
!  o ixtrcg+2 => indx of advected cfc11 tracer
!  o ixtrcg+3 => indx of advected cfc12 tracer
!
      if (trace_gas) then

         ii = 1.0
         num_levs = pverr

         call radclwmx(lchnk   ,ncol    ,                            &
                       lwupcgs ,t_local ,h2o_local       ,o3_local , &
                       pbr     ,pnm     ,pmln    ,piln    ,          &
                       n2o_local    ,ch4_local    ,                  &
                       cfc11_local    ,cfc12_local    ,              &
                       cld_local ,emis_local ,pmxrgn  ,nmxrgn  ,qrl, &
                       flns    ,flnt    ,flnsc   ,flntc   ,flwds   , &
                       flut    ,flutc   ,                            &
                       aerosol(:,:,idxVOLC)      ,fln     , flnc   )
         !   write(*,*) 'isnan flnsc'
         !write(*,*) 'lchnk = ',lchnk
         !write(*,*) 'ncol = ',ncol
         !write(*,*) 'lwupcgs = ',lwupcgs(1)
         !write(*,*) 't = ',t_local(1,:)
         !write(*,*) 'qm1 = ',h2o_local(1,:)
         !write(*,*) 'o3vmr = ',o3_local(1,:)
         !write(*,*) 'pbr = ',pbr(1,:)
         !write(*,*) 'pnm = ',pnm(1,:)
         !write(*,*) 'pmln = ',pmln(1,:)
         !write(*,*) 'piln = ',piln(1,:)
         !write(*,*) 'qm1_n2o = ',qm1(1,:,in2o)
         !write(*,*) 'qm1_ch4 = ',qm1(1,:,ich4)
         !write(*,*) 'qm1_if11 = ',qm1(1,:,if11)
         !write(*,*) 'cld = ',cld_local(1,:)
         !write(*,*) 'emis = ',emis_local(1,:)
         !write(*,*) 'pmxrgn = ',pmxrgn(1,:)
         !write(*,*) 'nmxrgn = ',nmxrgn(1)
         !write(*,*) 'flns = ',flns(1)
         !write(*,*) 'flnt= ',flnt(1)
         !write(*,*) 'flwds=',flwds(1)
         !write(*,*) 'aerosol = ',aerosol(1,:,idxVOLC)
         !stop 'prob with flns'

         !write(*,*) 'flnt = ',flnt(1)
         !write(*,*) 'flntc = ',flntc(1)
         !write(*,*) 'flns = ',flns(1)
         !write(*,*) 'flwds = ',flwds(1)
         !write(*,*) 'flut = ',flut(1)
         !write(*,*) 'flutc = ',flutc(1)
         !write(*,*) 'fln = ',fln(1,:)
         !write(*,*) 'flnc = ',flnc(1,:)
         !stop 'error radctl line 1176'


         call t_stopf("radclwmx")
      else
#ifndef OFFLINE
         call trcmix(lchnk   ,ncol    , &
                     pmid    ,n2o     ,ch4     ,                     &
                     cfc11   ,cfc12   )

         call t_startf("radclwmx")
         write(*,*) 'other radclwmx'
         call radclwmx(lchnk     ,ncol    ,                            &
                       lwupcgs   ,t_local    ,qm1(1,1,1)    ,o3vmr ,   &
                       pbr       ,pnm     ,pmln    ,piln    ,          &
                       n2o       ,ch4     ,cfc11   ,cfc12   ,          &
                       cld       ,emis    ,pmxrgn  ,nmxrgn  ,qrl     , &
                       flns      ,flnt    ,flnsc   ,flntc   ,flwds   , &
                       flut      ,flutc   ,                            &
                       aerosol(:,:,idxVOLC)        ,fln     ,flnc    )
         call t_stopf("radclwmx")
#endif
      endif

      
      !For cloud-overlap approximation purposes
      !allocate(out_cldfrac(num_profs,num_levs,pcols))
      !allocate(out_cldliqprof(num_profs,num_levs,pcols))
      !allocate(out_cldiceprof(num_profs,num_levs,pcols))
 
      min_cld_val = 1.e-5    !DRF modification

      do j=1,pcols

          pflx(j,0) = 0._r8
          do i=1,pverp
            pflx(j,i) = pint(j,i)
          end do

          num_modtran_calls_per_grid(j) = 0

         !First determine if clouds are present
          if (maxval(cld_local(j,:)).gt.1.0e-2) then 
             c_flag = .true.
          else
             c_flag = .false.
          endif
          ocean_flag = .false.

          if (c_flag) then

             do i=1,pver
               press_levs(i,j) = (pnm(j,pver+1-i)*1.0e-3)
             end do
          
             !Flip cloud-fraction so that it refers to layers of descending pressure
             do jj=1,pver
               flip_cldfrac(jj) = cld_local(j,1+pver-jj)
             end do

             !Convert cloud liquid and ice water path profiles to g/m^3
             !Reverse direction, lowest indices start on the lowest height

             do i=1,num_levs
                m_cloud_liq(i) = cliqwp_local(j,pver+1-i)/(7000.*dlog(pnm(j,pver+2-i)/pnm(j,1+pver-i)))
                m_cloud_ice(i) = cicewp_local(j,pver+1-i)/(7000.*dlog(pnm(j,pver+2-i)/pnm(j,1+pver-i)))
             end do

                   !write(*,*) "m_cloud_liq= ",m_cloud_liq
                   !write(*,*) "m_cloud_ice= ",m_cloud_ice
                   !write(*,*) "cicewp_local=",cicewp_local(j,:)
                   !stop

             !call generate_stochastic_clouds(num_levs,cloud_overlap_approx,changeSeed,press_levs(:,j),flip_cldfrac, &
             !          m_cloud_liq,m_cloud_ice,num_profs, &
             !          out_cldfrac(:,:,j),out_cldliqprof(:,:,j),out_cldiceprof(:,:,j),iscloudy_matrix )

             !Short cuts for not necessarily doing all of the cloud calculations, important to see if calculation has already been done

             !First reset values
             do jj=1,num_profs
                cloud_prof_id(jj) = 0
                weights(jj,j) = 0
                calc_or_nocalc(jj,j) = .false.
             end do          

             do jj=1,num_profs
                !calculate cloud id for each profile
                do kk=1,num_levs
                   if (iscloudy_matrix(jj,kk)) then
                      cloud_prof_id(jj) = cloud_prof_id(jj) + 2**kk
                   endif     
                end do 
             end do

             do jj=1,num_profs 
                calc_or_nocalc(jj,j) = .true.
                weights(jj,j) = 1
             end do

             do jj=1,num_profs
                !check the array to see if more than one of each cloud id exists
                do kk=jj+1,num_profs
                   if (weights(kk,j).gt.0 .and. cloud_prof_id(kk).eq.cloud_prof_id(jj)) then
                        calc_or_nocalc(kk,j) = .false.
                        weights(kk,j) = 0
                        weights(jj,j) = weights(jj,j)+1
                   endif
                end do
             end do   

             num_modtran_calls_per_grid(j) = 0
             do jj=1,num_profs
               if (weights(jj,j).gt.0) then
                  num_modtran_calls_per_grid(j) = num_modtran_calls_per_grid(j) + 1
               end if
             end do

             !write(*,*) "done with cloud generation"

             !1 calculation needed for clear-sky
          else
             num_modtran_calls_per_grid(j) = 1
          endif

      end do
  

      !loop through the Modtran calls
      do j= my_longitude_start, my_longitude_end !pcols
         ii = 6.0 + real(jj)*12.0/99.0
	 surface_offset = 0.0
         do i=1,pver
            full_atm(i,1) = -7.0*log(real(pnm(j,pver+1-i))/(1.0135*1.0e6)) 
            if (i.eq.1 .and. full_atm(i,1).lt.0.0) then
		surface_offset = abs(full_atm(i,1))
            end if
	    if (surface_offset.gt.0.0) then
		write(*,*) 'surf_offset = ',surface_offset,j,lchnk
	       full_atm(i,1) = full_atm(i,1) + surface_offset
            end if
            full_atm(i,2) = real(pnm(j,pver+1-i)*1.0e-3)
            if (i.eq.1) then
              full_atm(i,3) = real(ts(j))
            else
              full_atm(i,3) = real(t_local(j,pver+1-i))
            endif
            full_atm(i,4) = real(qm1(j,pver+1-i,1)*28.0/18.0*1.0e6) !H2O
            full_atm(i,5) = real(co2vmr*1.0e6)
            full_atm(i,6) = real(o3vmr(j,pver+1-i)*1.0e6)
            full_atm(i,7) = real(qm1(j,pver+1-i,in2o)*28.0/44.0*1.0e6) !N2O
            full_atm(i,8) = 0. !CO
            full_atm(i,9) = real(qm1(j,pver+1-i,ich4)*28.0/16.0*1.0e6) !CH4
!write(*,*) 'CH4 = ',full_atm(i,9)
!stop
            full_atm(i,10) = 0.
            full_atm(i,11) = 0.
            full_atm(i,12) = 0. 
            full_atm(i,13) = 0. 
            full_atm(i,14) = 0.
            full_atm(i,15) = 0.
            full_atm(i,16) = real(qm1(j,pver+1-i,if11)*28.0/137.0*1.0e6) !CFC-11
            full_atm(i,17) = real(qm1(j,pver+1-i,if12)*28.0/121.0*1.0e6) !CFC-12
            full_atm(i,18) = 0.
            full_atm(i,19) = 0.
            full_atm(i,20) = 0.
            full_atm(i,21) = 0.
            full_atm(i,22) = 0.
            full_atm(i,23) = 0.
            full_atm(i,24) = 0.
            full_atm(i,25) = 0.
            full_atm(i,26) = 0.
            full_atm(i,27) = 0.
            full_atm(i,28) = 0.
         end do
         
         !Do aerosols too
         
         !First 9 entries are read straight from AerosolOptics file
         aer_wavelen(1) = 0.2225
         aer_wavelen(2) = 0.255
         aer_wavelen(3) = 0.27
         aer_wavelen(4) = 0.28
         aer_wavelen(5) = 0.29
         aer_wavelen(6) = 0.3
         aer_wavelen(7) = 0.3275
         aer_wavelen(8) = 0.495
         aer_wavelen(9) = 0.67
         
         !Now map between pseudo bands and real bands
         aer_wavelen(10) = 0.701
         aer_wavelen(11) = 0.7095
         aer_wavelen(12) = 0.7245
         aer_wavelen(13) = 0.7725
         aer_wavelen(14) = 0.8210
         aer_wavelen(15) = 0.8615
         aer_wavelen(16) = 0.9110
         aer_wavelen(17) = 0.9445
         aer_wavelen(18) = 0.9725
         aer_wavelen(19) = 1.0380
         aer_wavelen(20) = 1.1010
         aer_wavelen(21) = 1.1345
         aer_wavelen(22) = 1.1830
         aer_wavelen(23) = 1.2515
         aer_wavelen(24) = 1.3100
         aer_wavelen(25) = 1.3320
         aer_wavelen(26) = 1.3510
         aer_wavelen(27) = 1.3860
         aer_wavelen(28) = 1.4325
         aer_wavelen(29) = 1.4740
         aer_wavelen(30) = 1.5090
         aer_wavelen(31) = 1.6260
         aer_wavelen(32) = 1.7475
         aer_wavelen(33) = 1.7820
         aer_wavelen(34) = 1.8120
         aer_wavelen(35) = 1.8710
         aer_wavelen(36) = 1.9345
         aer_wavelen(37) = 1.9720
         aer_wavelen(38) = 2.0180
         aer_wavelen(39) = 2.1835
         aer_wavelen(40) = 2.3620
         aer_wavelen(41) = 2.4330
         aer_wavelen(42) = 2.4860
         
         cam_band_mapping(1) = 1
         cam_band_mapping(2) = 2
         cam_band_mapping(3) = 3
         cam_band_mapping(4) = 5
         cam_band_mapping(5) = 5
         cam_band_mapping(6) = 6
         cam_band_mapping(7) = 7
         cam_band_mapping(8) = 8
         cam_band_mapping(9) = 9
         cam_band_mapping(10) = 11
         cam_band_mapping(11) = 10
         cam_band_mapping(12) = 11
         cam_band_mapping(13) = 10
         cam_band_mapping(14) = 11
         cam_band_mapping(15) = 10
         cam_band_mapping(16) = 11
         cam_band_mapping(17) = 12
         cam_band_mapping(18) = 11
         cam_band_mapping(19) = 10
         cam_band_mapping(20) = 11
         cam_band_mapping(21) = 12
         cam_band_mapping(22) = 11
         cam_band_mapping(23) = 10
         cam_band_mapping(24) = 11
         cam_band_mapping(25) = 12
         cam_band_mapping(26) = 13
         cam_band_mapping(27) = 14
         cam_band_mapping(28) = 13
         cam_band_mapping(29) = 12
         cam_band_mapping(30) = 11
         cam_band_mapping(31) = 10
         cam_band_mapping(32) = 11
         cam_band_mapping(33) = 12
         cam_band_mapping(34) = 13
         cam_band_mapping(35) = 14
         cam_band_mapping(36) = 13
         cam_band_mapping(37) = 12
         cam_band_mapping(38) = 11
         cam_band_mapping(39) = 10
         cam_band_mapping(40) = 11
         cam_band_mapping(41) = 12
         cam_band_mapping(42) = 13
         
         !Test out different aerosols            
         ! 1. MBCPHI_V
         ! 2. MBCPHO_V
         ! 3. MBG_V
         ! 4. MDUST1_V
         ! 5. MDUST2_V
         ! 6. MDUST3_V
         ! 7. MDUST4_V
         ! 8. MOCPHI_V
         ! 9. MOCPHO_V
         ! 10. MSSLT_V
         ! 11. MSU_V
         ! 12. MVOLC
         ! Modtran species: sulfate, dust, carb, sea-salt
         
         do kk=1,num_levs
           m_aerosol(kk,1) = 0.0
           m_aerosol(kk,2) = 0.0
           m_aerosol(kk,3) = 0.0
           m_aerosol(kk,4) = 0.0
           dummy3(kk) = 0.0
         end do
         do jj=1,42
           do kk=1,12
              o_aerosol(jj,kk) = 0.0
           end do
         end do
         do jj=1,50
            do kk=1,15
              do mm=1,4
                 pfs(jj,kk,mm) = 1.
                 pfs_sum(jj,kk,mm) = 0.
              end do
            end do
         end do

         aerosol_mult = 1.0e6

!!!!!!!!!!!!!!!
!Calculate Rayleigh optical depth at 0.55 um
         tau_rayleigh = 0.0956* full_atm(1,2)/1013.5
         !write(*,*) "tau_rayleigh = ",tau_rayleigh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Sulfate (aerosol index 1), converted to units of g/m^3 from mass below level
         ! units of aerosol are kg aerosol / kg air, k_ext units are m^2/g, need km^-1
         ! m_aerosol = ksul * aerosol * density, where density = p/R/T
         do kk=1,num_levs
            !k*u_path/dz where k = ext. coeff, u_path = aerosol*dp/g, and  dz is layer separation
            !max_rh_index = minloc(abs(rh_vec - rh(j,kk)),1)

            rhtrunc = rh(j,kk)
            rhtrunc = min(rh(j,kk),1._r8)
            krh = min(floor( rhtrunc * nrh ) + 1, nrh - 1)
            wrh = rhtrunc * nrh - krh

            !dummy_sulf = ksul(krh + 1, 9) * (wrh + 1) - ksul(krh, 9) * wrh
            !write(*,*) 'ksul9 = ',dummy_sulf

            dummy_sulf = ksul(krh + 1, 8) * (wrh + 1) - ksul(krh, 8) * wrh
            !write(*,*) 'ksul8 = ',dummy_sulf

            dummy3(kk) = dummy_sulf*1.e4*aerosol(j,kk,1)* &
                (pnm(j,kk)-pnm(j,kk+1))/980.616/(7.0*dlog(pnm(j,kk+1)/pnm(j,kk)))
         end do
         !dummy3(9) = 0.1 !DRF
         !dummy3(10) = 0.1 !DRF

         do kk=1,num_levs
            m_aerosol(kk,1) = real(abs(dummy3(1+num_levs-kk)))
            if (no_aerosol_flag .or. no_sulf_flag) then
               m_aerosol(kk,1) = 0.0
            endif 
         end do

         tau_aerosol_species(1) = 0.0
         do kk=1,num_levs
            tau_aerosol_species(1) = tau_aerosol_species(1) + & 
                m_aerosol(kk,1)*real((7.0*dlog(pnm(j,2+num_levs-kk)/pnm(j,1+num_levs-kk))))
         end do
         tau_sulf(j) = tau_aerosol_species(1)

         total_modtran_extinction = 0.
         do kk=1,num_levs-1
            halfdr = dble(full_atm(1+num_levs-kk,1)-full_atm(num_levs-kk,1))/2.
            zp_beg = (-7.0*dlog(pnm(j,kk+1)/1013500.)) !lower level
            zp_end = (-7.0*dlog(pnm(j,kk)/1013500.)) !upper level
            ecdbeg = 6371.23+zp_beg
            ecdend = 6371.23+zp_end
            ecdavg = 6371.23+(zp_beg+zp_end)/2.
            acoef = halfdr/(2.*(zp_end-zp_beg))
            abeg = acoef
            aend = -1.*abeg
            b2beg = (halfdr+ecdavg)*(0.25-abeg)/ecdbeg
            b2end = (halfdr+ecdavg)*(0.25+aend)/ecdend

            if (abs(dummy3(kk)).gt.0) then
               denln = dlog(dummy3(kk+1)/dummy3(kk))
            endif
            call denfac_copy(abeg,b2beg,denln,denfac_beg)
            call denfac_copy(aend,b2end,denln,denfac_end)
 
            aerosol_interpolated(kk) = abs(halfdr*(denfac_beg*dummy3(kk+1)+denfac_end*dummy3(kk)))
            total_modtran_extinction = total_modtran_extinction  + aerosol_interpolated(kk)
         end do
         total_modtran_extinction = abs(total_modtran_extinction)

         do kk=1,num_levs
            if(total_modtran_extinction.gt.0.) then
              m_aerosol(kk,1) = m_aerosol(kk,1)*tau_aerosol_species(1)/real(total_modtran_extinction)
            endif
         end do

         if (tau_sulf(j).gt.0.001*tau_rayleigh) then

           do kk=1,42
              o_aerosol(kk,1) = real(aertau(j,cam_band_mapping(kk),1)/aertau(j,8,1))
              o_aerosol(kk,2) = real(aertau(j,cam_band_mapping(kk),1)*(1.-aerssa(j,cam_band_mapping(kk),1))/aertau(j,8,1))
              o_aerosol(kk,3) = real(aerasm(j,cam_band_mapping(kk),1))
           end do
           do kk=1,15
              pfs_wvl_index = minloc(abs(aer_wavelen-real(phase_wvls(kk))),1)
              do mm=1,50
                 g = aerasm(j,cam_band_mapping(pfs_wvl_index),1)
                 pfs(mm,kk,1) = real((1.-g**2)/(1+g**2-2*g*cos(phase_angles(mm)))**1.5) !Henyey-Greenstein phase function
              end do
           end do 
        endif
         
           !   write(*,*) 'aertau2 = ',aertau(j,:,1)
           !   write(*,*) 'aerssa2 = ',aerssa(j,:,1)
           !   write(*,*) 'aerasm2 = ',aerasm(j,:,1)

           !   write(*,*) 'mod_ext = ',o_aerosol(:,1)
           !   write(*,*) 'mod_abs = ',o_aerosol(:,2)
           !   write(*,*) 'mod_asym = ',o_aerosol(:,3)

           !   write(*,*) 'tau_weight = ',tau_weight(kk)
           !   write(*,*) 'w_tau_weight = ',w_tau_weight(kk)
           !   write(*,*) 'g_tau_weight = ',g_tau_weight(kk)

           !   stop 'error radct.F90 line 1419'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!

         !Dust (aerosol indices 3-6) converted to units of g/m^3
         do kk=1,num_levs
            m_aerosol(kk,2) = 0.0
            o_aerosol(kk,4) = 0.0
            o_aerosol(kk,5) = 0.0
            o_aerosol(kk,6) = 0.0
            dummy3(kk) = 0.
         end do
         do kk=1,num_levs
               !Create effective dust profiles in units of km^-1
            dummy3(kk) = dummy3(kk) + (kdst(1,8)*1.0e4*aerosol(j,kk,idxDUSTfirst)* &
                (pnm(j,kk)-pnm(j,kk+1))/980.616)
            dummy3(kk) = dummy3(kk) + (kdst(2,8)*1.0e4*aerosol(j,kk,idxDUSTfirst+1)* &
                (pnm(j,kk)-pnm(j,kk+1))/980.616)
            dummy3(kk) = dummy3(kk) + (kdst(3,8)*1.0e4*aerosol(j,kk,idxDUSTfirst+2)* &
                (pnm(j,kk)-pnm(j,kk+1))/980.616)
            dummy3(kk) = dummy3(kk) + (kdst(4,8)*1.0e4*aerosol(j,kk,idxDUSTfirst+3)* &
                (pnm(j,kk)-pnm(j,kk+1))/980.616)
            dummy3(kk) = dummy3(kk)/(7.0*dlog(pnm(j,kk+1)/pnm(j,kk)))
         end do

         !dummy3(9) = 0.2 !DRF
         !dummy3(10) = 0.2 !DRF

         do kk=1,num_levs
            m_aerosol(kk,2) = real(abs(dummy3(1+num_levs-kk)))

            !scale aerosol by angstrom exponent to get 550 nm extinction
            if (no_aerosol_flag .or. no_dust_flag) then
               m_aerosol(kk,2) = 0.0
            endif
         end do

         !get effective profile extinction for modtran level-layer interpolation

         tau_aerosol_species(2) = 0.0
         do kk=1,num_levs
            tau_aerosol_species(2) = tau_aerosol_species(2) + & 
                m_aerosol(kk,2)*real((7.0*dlog(pnm(j,2+num_levs-kk)/pnm(j,1+num_levs-kk))))
         end do
         tau_dust(j) = tau_aerosol_species(2)

         total_modtran_extinction = 0.
         do kk=1,num_levs-1
            halfdr = dble(full_atm(1+num_levs-kk,1)-full_atm(num_levs-kk,1))/2.
            zp_beg = (-7.0*dlog(pnm(j,kk+1)/1013500.)) !lower level
            zp_end = (-7.0*dlog(pnm(j,kk)/1013500.)) !upper level
            ecdbeg = 6371.23+zp_beg
            ecdend = 6371.23+zp_end
            ecdavg = 6371.23+(zp_beg+zp_end)/2.
            acoef = halfdr/(2.*(zp_end-zp_beg))
            abeg = acoef
            aend = -1.*abeg
            b2beg = (halfdr+ecdavg)*(0.25-abeg)/ecdbeg
            b2end = (halfdr+ecdavg)*(0.25+aend)/ecdend

            if (abs(dummy3(kk)).gt.0) then
               denln = dlog(dummy3(kk+1)/dummy3(kk))
            endif
            call denfac_copy(abeg,b2beg,denln,denfac_beg)
            call denfac_copy(aend,b2end,denln,denfac_end)
 
            aerosol_interpolated(kk) = abs(halfdr*(denfac_beg*dummy3(kk+1)+denfac_end*dummy3(kk)))
            total_modtran_extinction = total_modtran_extinction  + aerosol_interpolated(kk)
         end do
         total_modtran_extinction = abs(total_modtran_extinction)
         !write(*,*) 'total_modtran_extinction = ',total_modtran_extinction 

         scale_fac = 1.

         do kk=1,num_levs
            if(total_modtran_extinction.gt.0.) then
              m_aerosol(kk,2) = m_aerosol(kk,2)*tau_aerosol_species(2)/real(total_modtran_extinction)*real(scale_fac)
            endif
         end do

         if (tau_aerosol_species(2).gt.0.001*tau_rayleigh) then
             do kk=1,42
               if (aertau(j,8,4).gt.0.) then
                 o_aerosol(kk,4) = real(aertau(j,cam_band_mapping(kk),4)/(aertau(j,8,4)*scale_fac))
                 o_aerosol(kk,5) = real(aertau(j,cam_band_mapping(kk),4)* &
                                        (1.-aerssa(j,cam_band_mapping(kk),4))/(scale_fac*aertau(j,8,4)))
                 o_aerosol(kk,6) = real(aerasm(j,cam_band_mapping(kk),4))
               endif
             end do
             do kk=1,15
                pfs_wvl_index = minloc(abs(aer_wavelen-real(phase_wvls(kk))),1)
                do mm=1,50
                   g = aerasm(j,cam_band_mapping(pfs_wvl_index),4)
                   pfs(mm,kk,2) = real((1.-g**2)/(1+g**2-2*g*cos(phase_angles(mm)))**1.5) !Henyey-Greenstein phase function
                end do
             end do

              !write(*,*) 'mod_ext = ',o_aerosol(:,4)
              !write(*,*) 'mod_abs = ',o_aerosol(:,5)
              !write(*,*) 'mod_asym = ',o_aerosol(:,6)

              !write(*,*) 'pfs(2,:) = ',pfs(:,1,2)
              !write(*,*) 'phase_angles = ',phase_angles
              !stop 'error radct.F90 line 1504'

         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Soot (aerosol indices 7-10) converted to units of g/m^3
         
         do kk=1,num_levs
            m_aerosol(kk,3) = 0.0
            o_aerosol(kk,7) = 0.0
            o_aerosol(kk,8) = 0.0
            o_aerosol(kk,9) = 0.0
            dummy3(kk) = 0.0
         end do 
         !Get RH max first
         do kk=1,num_levs
            max_rh_index = minloc(abs(rh_vec - rh(j,kk)),1)
            dummy3(kk) = dummy3(kk) + (kcphob(8)*1.0e4*aerosol(j,kk,7)* &
                (pnm(j,kk)-pnm(j,kk+1))/980.616)
            dummy3(kk) = dummy3(kk) + (kcb(8)*1.0e4*aerosol(j,kk,8)* &
                (pnm(j,kk)-pnm(j,kk+1))/980.616)
            dummy3(kk) = dummy3(kk) + (kcphil(max_rh_index,8)*1.0e4*aerosol(j,kk,9)* &
                (pnm(j,kk)-pnm(j,kk+1))/980.616)
            dummy3(kk) = dummy3(kk) + (kcb(8)*1.0e4*aerosol(j,kk,10)* &
                (pnm(j,kk)-pnm(j,kk+1))/980.616)
            dummy3(kk) = dummy3(kk)/(7.0*dlog(pnm(j,kk+1)/pnm(j,kk)))
         end do
         !dummy3(9) = 0.3 !DRF
         !dummy3(10) = 0.3 !DRF

         do kk=1,num_levs
            m_aerosol(kk,3) = real(abs(dummy3(1+num_levs-kk)))
            if (no_aerosol_flag .or. no_carb_flag) then
               m_aerosol(kk,3) =  0.0
            endif
         end do

         tau_aerosol_species(3) = 0.0
         do kk=1,pver-1
            tau_aerosol_species(3) = tau_aerosol_species(3) + & 
                m_aerosol(kk,3)*real((7.0*dlog(pnm(j,2+num_levs-kk)/pnm(j,1+num_levs-kk))))
         end do
         tau_carb(j) = tau_aerosol_species(3)

         total_modtran_extinction = 0.
         do kk=1,num_levs-1
            halfdr = dble(full_atm(1+num_levs-kk,1)-full_atm(num_levs-kk,1))/2.
            zp_beg = (-7.0*dlog(pnm(j,kk+1)/1013500.)) !lower level
            zp_end = (-7.0*dlog(pnm(j,kk)/1013500.)) !upper level
            ecdbeg = 6371.23+zp_beg
            ecdend = 6371.23+zp_end
            ecdavg = 6371.23+(zp_beg+zp_end)/2.
            acoef = halfdr/(2.*(zp_end-zp_beg))
            abeg = acoef
            aend = -1.*abeg
            b2beg = (halfdr+ecdavg)*(0.25-abeg)/ecdbeg
            b2end = (halfdr+ecdavg)*(0.25+aend)/ecdend

            if (abs(dummy3(kk)).gt.0) then
               denln = dlog(dummy3(kk+1)/dummy3(kk))
            endif
            call denfac_copy(abeg,b2beg,denln,denfac_beg)
            call denfac_copy(aend,b2end,denln,denfac_end)
 
            aerosol_interpolated(kk) = abs(halfdr*(denfac_beg*dummy3(kk+1)+denfac_end*dummy3(kk)))
            total_modtran_extinction = total_modtran_extinction  + aerosol_interpolated(kk)
         end do
         total_modtran_extinction = abs(total_modtran_extinction)
         do kk=1,num_levs
            if(total_modtran_extinction.gt.0.) then
               m_aerosol(kk,3) = m_aerosol(kk,3)*tau_aerosol_species(3)/real(total_modtran_extinction)
            endif
         end do
         
         if (tau_aerosol_species(3).gt.0.001*tau_rayleigh) then
               do kk=1,42
                 if (aertau(j,8,3).gt.0.) then
                   o_aerosol(kk,7) = real(aertau(j,cam_band_mapping(kk),3)/aertau(j,8,3))
                   o_aerosol(kk,8) = real(aertau(j,cam_band_mapping(kk),3)*(1.-aerssa(j,cam_band_mapping(kk),3))/aertau(j,8,3))
                 endif
                 o_aerosol(kk,9) = real(aerasm(j,cam_band_mapping(kk),3))
               end do
              do kk=1,15
                pfs_wvl_index = minloc(abs(aer_wavelen-real(phase_wvls(kk))),1)

                do mm=1,50
                   g = aerasm(j,cam_band_mapping(pfs_wvl_index),3)
                   pfs(mm,kk,3) = real((1.-g**2)/(1+g**2-2*g*cos(phase_angles(mm)))**1.5) !Henyey-Greenstein phase function
                end do
              end do
           !write(*,*) 'o_aerosol(:,7) = ',o_aerosol(:,7)
           !write(*,*) 'o_aerosol(:,8) = ',o_aerosol(:,8)
           !write(*,*) 'o_aerosol(:,9) = ',o_aerosol(:,9)
           !stop 'error radctl line 1569'
         endif
 
!!!!!!!!!!!!!!!!!!!!!!!!
!Sea-salt

         !Sea-salt (aerosol index 2), converted to units of g/m^3 from mass below level

         do kk=1,num_levs
            m_aerosol(kk,4) = 0.0
            o_aerosol(kk,10) = 0.0
            o_aerosol(kk,11) = 0.0
            o_aerosol(kk,12) = 0.0
            dummy3(kk) = 0.0
         end do

         do kk=1,num_levs
            rhtrunc = rh(j,kk)
            rhtrunc = min(rh(j,kk),1._r8)
            krh = min(floor( rhtrunc * nrh ) + 1, nrh - 1)
            wrh = rhtrunc * nrh - krh

            !dummy_sslt = ksslt(krh + 1, 9) * (wrh + 1) - ksslt(krh, 9) * wrh
            !write(*,*) 'ksslt9 = ',dummy_sslt

            dummy_sslt = ksslt(krh + 1, 8) * (wrh + 1) - ksslt(krh, 8) * wrh
            !write(*,*) 'ksslt8= ',dummy_sslt

            dummy_sslt = ksslt(krh + 1, 8) * (wrh + 1) - ksslt(krh, 8) * wrh
            dummy3(kk) = dummy_sslt*1.e4*aerosol(j,kk,2)* &
                (pnm(j,kk)-pnm(j,kk+1))/980.616/(7.0*dlog(pnm(j,kk+1)/pnm(j,kk)))
         end do

         do kk=1,num_levs
            m_aerosol(kk,4) = real(abs(dummy3(1+num_levs-kk)))
            if (no_aerosol_flag .or. no_sslt_flag) then
               m_aerosol(kk,4) = 0.0
            endif
         end do
 
         tau_aerosol_species(4) = 0.0
         do kk=1,num_levs
            tau_aerosol_species(4) = tau_aerosol_species(4) + & 
                m_aerosol(kk,4)*real((7.0*dlog(pnm(j,2+num_levs-kk)/pnm(j,1+num_levs-kk))))
         end do
         tau_sslt(j) = tau_aerosol_species(4)

         total_modtran_extinction = 0.
         do kk=1,num_levs-1
            halfdr = dble(full_atm(1+num_levs-kk,1)-full_atm(num_levs-kk,1))/2.
            zp_beg = (-7.0*dlog(pnm(j,kk+1)/1013500.)) !lower level
            zp_end = (-7.0*dlog(pnm(j,kk)/1013500.)) !upper level
            ecdbeg = 6371.23+zp_beg
            ecdend = 6371.23+zp_end
            ecdavg = 6371.23+(zp_beg+zp_end)/2.
            acoef = halfdr/(2.*(zp_end-zp_beg))
            abeg = acoef
            aend = -1.*abeg
            b2beg = (halfdr+ecdavg)*(0.25-abeg)/ecdbeg
            b2end = (halfdr+ecdavg)*(0.25+aend)/ecdend

            if (abs(dummy3(kk)).gt.0) then
               denln = dlog(dummy3(kk+1)/dummy3(kk))
            endif
            call denfac_copy(abeg,b2beg,denln,denfac_beg)
            call denfac_copy(aend,b2end,denln,denfac_end)
 
            aerosol_interpolated(kk) = abs(halfdr*(denfac_beg*dummy3(kk+1)+denfac_end*dummy3(kk)))
            total_modtran_extinction = total_modtran_extinction  + aerosol_interpolated(kk)
         end do

         total_modtran_extinction = abs(total_modtran_extinction)
         !write(*,*) 'total_ext =',total_modtran_extinction
         !write(*,*) 'tau = ',tau_aerosol_species
         !write(*,*) 'm_aerosol(:,4) = ',m_aerosol(:,4)
         do kk=1,num_levs
            if(total_modtran_extinction.gt.0.) then
               m_aerosol(kk,4) = m_aerosol(kk,4)*tau_aerosol_species(4)/real(total_modtran_extinction)
            endif
         end do

         do kk=1,42
            tau_weight(kk) = 0.
            tau_total(kk) = 0.
            w_tau_weight(kk) = 0.
            w_tau_total(kk) = 0.
            w_tau_total2(kk) = 0.
         end do
         if (tau_aerosol_species(4).gt.0.) then

           do kk=1,42
              o_aerosol(kk,10) = real(aertau(j,cam_band_mapping(kk),2)/aertau(j,8,2))
              o_aerosol(kk,11) = real(aertau(j,cam_band_mapping(kk),2)*(1.-aerssa(j,cam_band_mapping(kk),2))/aertau(j,8,2))
              o_aerosol(kk,12) = real(aerasm(j,cam_band_mapping(kk),2))
           end do

           do kk=1,15
              pfs_wvl_index = minloc(abs(aer_wavelen-real(phase_wvls(kk))),1)
              do mm=1,50
                  g = aerasm(j,cam_band_mapping(pfs_wvl_index),2)
                  pfs(mm,kk,4) = real((1.-g**2)/(1+g**2-2*g*cos(phase_angles(mm)))**1.5) !Henyey-Greenstein phase function
              end do
           end do

           !write(*,*) 'o_aerosol(:,10) = ',o_aerosol(:,10)
           !write(*,*) 'o_aerosol(:,11) = ',o_aerosol(:,11)
           !write(*,*) 'o_aerosol(:,12) = ',o_aerosol(:,12)

           !write(*,*) 'tau = ',tau_sslt(j)
           !write(*,*) 'pfs(1) = ',pfs(:,1,4)
           !write(*,*) 'pfs(2) = ',pfs(:,2,4)
           !write(*,*) 'pfs(3) = ',pfs(:,3,4)

           !stop 'error radctl line 1883'
         endif

         !write(*,*) "tau_sulfate = ",tau_aerosol_species(1)
         !write(*,*) "tau_dust = ",tau_aerosol_species(2)
         !write(*,*) "tau_carb = ",tau_aerosol_species(3)
         !write(*,*) "tau_sslt = ",tau_aerosol_species(4)
         !write(*,*) "tau_total = ",tau_aerosol_species(1)+tau_aerosol_species(2)+ &
         !                         tau_aerosol_species(3)+tau_aerosol_species(4)         
         aod(j) = tau_aerosol_species(1)+tau_aerosol_species(2)+tau_aerosol_species(3)+tau_aerosol_species(4)         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Get data from phase_functions netcdf file, tranlate to matrix(50,15,4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Get cloud optical properties

            if (lw_flag) then
               do i=1,24
                  wvl_cld(i) = 4.5 + 10.*i
               end do
               do i=1,num_levs
                  call get_ext_ssa_asym_liq_lw(rel(j,num_levs+1-i),wvl_cld,ext_liq_cld(:,i),ssa_liq_cld(:,i),asym_liq_cld(:,i))
                  call get_ext_ssa_asym_ice_lw(rei_local(j,num_levs+1-i),wvl_cld,ext_ice_cld(:,i),ssa_ice_cld(:,i), &
                                               asym_ice_cld(:,i))
               end do
            else
               do i=1,num_levs
                  call get_ext_ssa_asym_liq(rel(j,num_levs+1-i),wvl_cld,ext_liq_cld(:,i),ssa_liq_cld(:,i),asym_liq_cld(:,i))
                  call get_ext_ssa_asym_ice(rei_local(j,num_levs+1-i),wvl_cld,ext_ice_cld(:,i),ssa_ice_cld(:,i),asym_ice_cld(:,i))
               end do
            endif

            !Weight the cloud optical properties by extinction
            do i=1,24 !reset vals by wavelength to 0
              dummy_ext(i) = 0.
              dummy_ssa(i) = 0.
              dummy_asym(i) = 0.
            end do
            do i=1,num_levs
               dummy_ext = dummy_ext+ext_liq_cld(:,i)*cliqwp_local(j,num_levs+1-i)
               dummy_ssa = dummy_ssa+ssa_liq_cld(:,i)*cliqwp_local(j,num_levs+1-i)
               do mm=1,24
                 dummy_asym(mm) = dummy_asym(mm)+asym_liq_cld(mm,i)*cliqwp_local(j,num_levs+1-i)*(1.+ssa_liq_cld(mm,i))
               end do
            end do

            do i=1,24
              dummy_sum = sum(cliqwp_local(j,:))
              wvl_cld_real(i) = real(wvl_cld(i))
              if (dummy_sum.gt.0.) then
                liq_ext_avg(i) = dummy_ext(i)/sum(cliqwp_local(j,:))
                liq_ext_avg_real(i) = real(liq_ext_avg(i))
                liq_ssa_avg(i) = dummy_ssa(i)/sum(cliqwp_local(j,:))
                liq_ssa_avg_real(i) = real(liq_ssa_avg(i))
                if (liq_ssa_avg_real(i)==0.) then
                   liq_ssa_avg_real(i) = -1.0e-6
                endif
                dummy_sum = 0.
                do mm=1,num_levs
                   dummy_sum = dummy_sum+cliqwp_local(j,mm)*(1.+ssa_liq_cld(i,mm))
                end do
                liq_asym_avg(i) = dummy_asym(i)/dummy_sum
                liq_asym_avg_real(i) = real(liq_asym_avg(i))
              endif
            end do

            do i=1,24 !reset vals by wavelength to 0
              dummy_ext(i) = 0.
              dummy_ssa(i) = 0.
              dummy_asym(i) = 0.
            end do
            do i=1,num_levs
               dummy_ext = dummy_ext+ext_ice_cld(:,i)*cicewp_local(j,num_levs+1-i)
               dummy_ssa = dummy_ssa+ssa_ice_cld(:,i)*cicewp_local(j,num_levs+1-i)
               do mm=1,24
                  dummy_asym(mm) = dummy_asym(mm)+asym_ice_cld(mm,i)*cicewp_local(j,num_levs+1-i)*(1.+ssa_ice_cld(mm,i))
               end do
            end do

            do i=1,24
              dummy_sum = sum(cicewp_local(j,:))
              if (dummy_sum.gt.0.) then
              ice_ext_avg(i) = dummy_ext(i)/sum(cicewp_local(j,:))
              ice_ext_avg_real(i) = real(ice_ext_avg(i))
              ice_ssa_avg(i) = dummy_ssa(i)/sum(cicewp_local(j,:))
              ice_ssa_avg_real(i) = real(ice_ssa_avg(i))
              if (ice_ssa_avg_real(i)==0.) then
                 ice_ssa_avg_real(i) = -1.0e-6
              endif
              dummy_sum = 0.
              do mm=1,num_levs
                 dummy_sum = dummy_sum+cicewp_local(j,mm)*(1.+ssa_ice_cld(i,mm))
              end do
              ice_asym_avg(i) = dummy_asym(i)/dummy_sum
              ice_asym_avg_real(i) = real(ice_asym_avg(i))
              endif
            end do

            !write(*,*) "liq_ext_avg= ",liq_ext_avg
            !write(*,*) "liq_ssa_avg= ",liq_ssa_avg
            !write(*,*) "liq_asym_avg= ",liq_asym_avg
            !!write(*,*) "ssa_liq_cld= ",ssa_liq_cld(:,1)
            !write(*,*) 'cldtau_ice = ',cldtau_ice_lw(j) 
            !write(*,*) "ice_ext_avg= ",ice_ext_avg
            !write(*,*) "ice_ssa_avg= ",ice_ssa_avg
            !write(*,*) "ice_asym_avg= ",ice_asym_avg
            !stop 'error radctl line 2000'

            ! Get BRDF parameters
            
            !modis channels
            !aer_wavelen(1) = 0.645
            !aer_wavelen(2) = 0.8585
            !aer_wavelen(3) = 0.469
            !aer_wavelen(4) = 0.555
            !aer_wavelen(5) = 1.240
            !aer_wavelen(6) = 1.640
            !aer_wavelen(7) = 2.130            
            
            a_flag = .true.

!!! Aerosol mapping of spectral bands
            num_levs = pverr
            jj = 0

            !Modtran temporary file name
            fname_index_last = 120 
            do jj=120,1,-1
              dummy_char = ipath_in(jj)
              if (dummy_char .eq. ' ') then 
                fname_index_last = fname_index_last - 1
              endif
            end do
            fname_index_first = 120
            fname_flag  = .true.
            do jj=120,1,-1
              dummy_char = ipath_in(jj)
              if (dummy_char .ne. '/' .and. fname_flag) then 
                fname_index_first = fname_index_first - 1
              else
                fname_flag = .false.
              endif
            end do
            fname_index_first = fname_index_first + 1
            fname = 'modtran_out_' 
            do jj=fname_index_first,fname_index_last-2
              dummy_char = ipath_in(jj)
              fname =  fname(1:12+jj-fname_index_first) // dummy_char 
            end do

            if (lchnk.lt.10) then
               write(fval_1,fmt='(i1)') lchnk
               fname = fname(1:12-fname_index_first+fname_index_last-1)//fval_1
               fname_length = 12 - fname_index_first + fname_index_last
            elseif (lchnk.ge.10 .and. lchnk.lt.100) then
               write(fval_2,fmt='(i2)') lchnk
               fname = fname(1:12-fname_index_first+fname_index_last-1)//fval_2
               fname_length = 12 - fname_index_first + fname_index_last+1 
            elseif (lchnk.ge.100) then
               write(fval_3,fmt='(i3)') lchnk
               fname = fname(1:12-fname_index_first+fname_index_last-1)//fval_3
               fname_length = 12 - fname_index_first + fname_index_last+2
            endif

            filler_period = '.'
            fname = fname(1:fname_length)//filler_period
            fname_length = fname_length + 1

            if (j.lt.10) then
               write(fval_1,fmt='(i1)') j
               fname = fname(1:fname_length)//fval_1
               fname_length = fname_length + 1
            elseif (j.ge.10 .and. j.lt.100) then
               write(fval_2,fmt='(i2)') j
               fname = fname(1:fname_length)//fval_2
               fname_length = fname_length + 2
            elseif (j.ge.100) then
               write(fval_3,fmt='(i3)') j
               fname = fname(1:fname_length)//fval_3
               fname_length = fname_length + 3
            endif
            !write(*,*) "fname = ",fname(1:fname_length)
            !write(*,*) "j = ",j
            !write(*,*) "lchnk = ",lchnk
            
            dummy2 = 0

            !write(*,*) 'snow_flag = ',landtype
            !stop 'error on line radctl.F90 in 1827'
            ice_flag_local = .false.
            
            !Define brdf_param2,brdf_len,brdf_wvl
            if (icefrac(j).ge.0.001) then
               ice_flag_local = .true.
            endif

            !MODIS bands: 20, 22, 23, 29, 31, 32
            spec_emissivity_len = 6
            spec_emissivity_wvl(1) = 3.75
            spec_emissivity_wvl(2) = 3.959
            spec_emissivity_wvl(3) = 4.05
            spec_emissivity_wvl(4) = 8.55
            spec_emissivity_wvl(5) = 11.03
            spec_emissivity_wvl(6) = 12.02
            if (land_flag(j)) then
               do jj=1,6
                  spec_emissivity_vals(jj) = real(emis_array(j,jj))
               end do
            elseif (icefrac(j).ge.0.001) then
               do jj=1,6
                  spec_emissivity_vals(jj) = 0.5
               end do
            else
               !ocean albedo/emissivity (from Masuda et al, Remote Sens. Env. 1999, 24(2) 313-329)
               spec_emissivity_vals(1) =  0.9745
               spec_emissivity_vals(2) =  0.97705
               spec_emissivity_vals(3) =  0.97760
               spec_emissivity_vals(4) =  0.98500
               spec_emissivity_vals(5) =  0.99252
               spec_emissivity_vals(6) =  0.98894
            endif

            do jj=1,15
              brdf_wvl(jj) = 0.
              do kk=1,3
                brdf_param2(jj,kk) = 0.
              end do
            end do
            if (land_flag(j)) then !land
               brdf_wvl(1) = 0.47
               brdf_wvl(2) = 0.56
               brdf_wvl(3) = 0.64
               brdf_wvl(4) = 0.86
               brdf_wvl(5) = 1.24
               brdf_wvl(6) = 1.64
               brdf_wvl(7) = 2.13
               do jj=1,7
                  do kk=1,3
                     brdf_param2(jj,kk) = real(brdf_param(j,jj,kk)/1000.0)
                  end do
               end do
               brdf_len = 7
               ocean_flag = .false.
            elseif (icefrac(j).ge.0.001) then
               !icy
               ice_flag_local = .true.

               brdf_len = 7

               brdf_wvl(1) = 0.47
               brdf_wvl(2) = 0.56
               brdf_wvl(3) = 0.64
               brdf_wvl(4) = 0.86
               brdf_wvl(5) = 1.24
               brdf_wvl(6) = 1.64
               brdf_wvl(7) = 2.13

               brdf_param2(1,1) = 0.76
               brdf_param2(2,1) = 0.74
               brdf_param2(3,1) = 0.72
               brdf_param2(4,1) = 0.6
               brdf_param2(5,1) = 0.25
               brdf_param2(6,1) = 0.07
               brdf_param2(7,1) = 0.06

               do ii=1,7
                 brdf_param2(ii,2) = 0. !geometric
                 brdf_param2(ii,3) = 0. !volumetric
               end do

               ocean_flag = .false.
            else
               !ocean reflectance data

               !cos_zen_real = real(coszrs(j))
               !if (cos_zen_real.lt.0.) then
               !   cos_zen_real = 0.001
               !endif
               !write(*,*) 'cos_zen = ',cos_zen_real

               ocean_flag = .true.

               !no chlorophyll, 5.0 m/s wind-speed
               !real_windspeed = real(windspeed(j))
               !call zjin_wgtlut(0.0,cos_zen_real,real_windspeed,0.2,fuspecalb,fu_wvl)

               brdf_len = 15 

               do jj=1,brdf_len
                  brdf_wvl(jj) = 0.4 + real(jj-1)/10.
                  brdf_param2(jj,1) = 0.4 + real(jj-1)/10.
                  brdf_param2(jj,2) = real_windspeed
                  brdf_param2(jj,3) = 0.
               end do
            endif

            !write(*,*) "lat index = ",latitude_index
            !write(*,*) "j = ",j
            !write(*,*) "windspeed = ",windspeed(j)
            !stop

            num_ac_wvl = 42
            num_hres = 3100
            num_lres = 910
            !write(*,*) "mcdate = ",mcdate
            mcdate_local = mcdate - 1e4*floor(mcdate/1e4)
            
            !Write output to tape5 for debugging
            latval_real = real(latval)
            lonval_real = real(lonval(j))
            co2vmr_real = real(1.0e6*co2vmr)
            gndalt_real = full_atm(1,1)
            surftemp_real = ts(j)
            gmt_real = 13.5

            !Calculate solar azimuth angle
            pi_real = 3.1415926535898
            psipo_real = real(solar_az(j))
            zen_real = acos(coszrs(j))*180./pi_real

            if (lw_flag) then
              !v1_real = 20.  !cm^-1
              !v2_real = 2200.
              !dv_input_real = 5.0
              v1_real = 200.  !cm^-1
              v2_real = 2000.
              dv_input_real = 5.0
            else
              v1_real = 300. !nm
              v2_real = 2500.
              dv_input_real = 5.0
            endif
            do jj=1,15
               phase_wvls_real(jj) = real(phase_wvls(jj))
            end do
            do jj=1,50
               phase_angles_real(jj) = real(phase_angles(jj))*180./3.1415926536
            end do
            do jj=1,num_levs-1
               z_cld_val_real(jj) = real(full_atm(jj,1))
            end do

            !write(*,*) "co2vmr = ",co2vmr_real
             !for writing out .tp5 file
             do jj=1,num_levs
                m_cloud_liq_real(jj) = real(cliqwp_local(j,pver+1-jj)/(7000.*dlog(pnm(j,pver+2-jj)/pnm(j,1+pver-jj))))
                m_cloud_ice_real(jj) = real(cicewp_local(j,pver+1-jj)/(7000.*dlog(pnm(j,pver+2-jj)/pnm(j,1+pver-jj))))
             end do
             if (.not.land_flag(j) .and. icefrac(j).lt.0.001) then
                ocean_flag = .true.
             else
                ocean_flag = .false.
             endif

            !write(*,*) "land_flag = ",land_flag(j)
            !write(*,*) "snow_flag = ",snow_flag(j)
            !write(*,*) "c_flag = ",c_flag
            !write(*,*) "o_aerosol = ",o_aerosol(:,4)
            !write(*,*) "full_atm(1,3) = ",full_atm(1,3)
            !write(*,*) "latval = ",latval_real
            !write(*,*) "lonval = ",lonval(j)
            !write(*,*) "a_flag = ",a_flag
            !write(*,*) "full_atm = ",full_atm(:,1)
            !write(*,*) "phase_wvls = ",phase_wvls
            !write(*,*) "land_flag = ",land_flag(j)
            !write(*,*) "brdf_len = ",brdf_len
            !write(*,*) "brdf_wvl = ",brdf_wvl
            !write(*,*) "brdf_param2 = ",brdf_param2
            !write(*,*) "pfs = ",pfs
            !write(*,*) "phase_wvls = ",real(phase_wvls)
            !write(*,*) "phase_angs = ",real(phase_angles)
            !write(*,*) "wvl = ",o_aerosol(:,4)

            !For all-sky calculations, perform multiple calculations and average per
            !the maximum/random cloud overlap approximation

            if (maxval(cld_local(j,:)).gt.1.0e-2) then 
               c_flag = .true.
            else
               c_flag = .false.
            endif

            !write(*,*) "my_longitude_end = ",my_longitude_end
            !   write(*,*) 'full_atm(:,1) = ',full_atm(:,1)
            !   write(*,*) 'full_atm(:,2) = ',full_atm(:,2)
            !   write(*,*) 'full_atm(:,3) = ',full_atm(:,3)
            !   write(*,*) 'full_atm(:,4) = ',full_atm(:,4)
            !   write(*,*) 'full_atm(:,5) = ',full_atm(:,5)
            !   write(*,*) 'full_atm(:,6) = ',full_atm(:,6)
            !   write(*,*) 'full_atm(:,7) = ',full_atm(:,7)
            !   write(*,*) 'full_atm(:,8) = ',full_atm(:,8)
            !   write(*,*) 'full_atm(:,9) = ',full_atm(:,9)
            !   write(*,*) 'full_atm(:,16) = ',full_atm(:,16)
            !   write(*,*) 'full_atm(:,17) = ',full_atm(:,17)
            !   write(*,*) 'cliqwp_local = ',cliqwp_local(j,:)
            !   write(*,*) 'cicewp_local = ',cicewp_local(j,:)

           !write(*,*) 'max liq cld = ',maxval(cliqwp)
           !write(*,*) 'max ice cld = ',maxval(cicewp)

            if (ice_flag_local .and. .not. lw_flag) then
               !Combine icy & ice-free conditions
               !write(*,*) 'icy combo'

               !allocate arrays
               allocate(tmp_rad_lres_clr_noice(wvlng2))
               allocate(tmp_rad_lres_clr_ice(wvlng2))
               allocate(tmp_rad_hres_clr_noice(wvlng_hres2))
               allocate(tmp_rad_hres_clr_ice(wvlng_hres2))
               allocate(tmp_rad_lres_all_noice(wvlng2))
               allocate(tmp_rad_lres_all_ice(wvlng2))
               allocate(tmp_rad_hres_all_noice(wvlng_hres2))
               allocate(tmp_rad_hres_all_ice(wvlng_hres2))

               allocate(tmp_bbupdif_clr_noice(pver))
               allocate(tmp_bbupdif_clr_ice(pver))
               allocate(tmp_bbupdif_all_noice(pver))
               allocate(tmp_bbupdif_all_ice(pver))

               allocate(tmp_bbdndif_clr_noice(pver))
               allocate(tmp_bbdndif_clr_ice(pver))
               allocate(tmp_bbdndif_all_noice(pver))
               allocate(tmp_bbdndif_all_ice(pver))

               allocate(tmp_bbdndir_clr_noice(pver))
               allocate(tmp_bbdndir_clr_ice(pver))
               allocate(tmp_bbdndir_all_noice(pver))
               allocate(tmp_bbdndir_all_ice(pver))

               allocate(tmp_flxupdif_clr_noice(pver,num_bands))
               allocate(tmp_flxupdif_clr_ice(pver,num_bands))
               allocate(tmp_flxupdif_all_noice(pver,num_bands))
               allocate(tmp_flxupdif_all_ice(pver,num_bands))

               allocate(tmp_flxdndif_clr_noice(pver,num_bands))
               allocate(tmp_flxdndif_clr_ice(pver,num_bands))
               allocate(tmp_flxdndif_all_noice(pver,num_bands))
               allocate(tmp_flxdndif_all_ice(pver,num_bands))

               allocate(tmp_flxdndir_clr_noice(pver,num_bands))
               allocate(tmp_flxdndir_clr_ice(pver,num_bands))
               allocate(tmp_flxdndir_all_noice(pver,num_bands))
               allocate(tmp_flxdndir_all_ice(pver,num_bands))

               allocate(tmp_difflux_clr_noice(wvlng2))
               allocate(tmp_difflux_clr_ice(wvlng2))
               allocate(tmp_difflux_all_noice(wvlng2))
               allocate(tmp_difflux_all_ice(wvlng2))

               !Assign brdf for sea-ice conditions
               do jj=1,15
                   brdf_wvl(jj) = 0.
                   do kk=1,3
                      brdf_param2(jj,kk) = 0.
                   end do
               end do

               brdf_len = 7
               brdf_wvl(1) = 0.47
               brdf_wvl(2) = 0.56
               brdf_wvl(3) = 0.64
               brdf_wvl(4) = 0.86
               brdf_wvl(5) = 1.24
               brdf_wvl(6) = 1.64
               brdf_wvl(7) = 2.13

               brdf_param2(1,1) = 0.76
               brdf_param2(2,1) = 0.74
               brdf_param2(3,1) = 0.72
               brdf_param2(4,1) = 0.6
               brdf_param2(5,1) = 0.25
               brdf_param2(6,1) = 0.07
               brdf_param2(7,1) = 0.06

               do ii=1,7
                 brdf_param2(ii,2) = 0. !geometric
                 brdf_param2(ii,3) = 0. !volumetric
               end do
               ocean_flag = .false.

               if (modtran_flag) then
!               if (j.eq.84 .and. lchnk.eq.118) then
		!  write(*,*) 'cloud ice is ',cicewp_local(j,:)
               if (j.lt.20 .and. lchnk.lt.20) then
!                  call run_modpcrtm
                 call run_modtran_or_pcrtm(lchnk,j,nconfig_v(j),wgtv_v(j,:),ccon_v(j,:,:),totwgt_v(j),cliqwp_local(j,:),&
                               cicewp_local(j,:),cld_local(j,:),rei(j,:),rel(j,:),fname, &
                               fname_length,surftemp_real,co2vmr_real,gndalt_real,jday_input,gmt_real,zen_real,psipo_real, &
                               a_flag,c_flag,full_atm,pnm(j,:),v1_real,v2_real,ocean_flag,brdf_len,brdf_wvl,brdf_param2, &
                               spec_emissivity_len,spec_emissivity_wvl,spec_emissivity_vals,&
                               num_levs,m_aerosol,o_aerosol,aer_wavelen,pfs,phase_wvls_real,phase_angles_real,wvl_cld, &
                               ext_liq_cld,ssa_liq_cld,asym_liq_cld,ext_ice_cld,ssa_ice_cld,asym_ice_cld,z_cld_val_real, &
                               dv_input_real,num_lres,num_hres,num_ac_wvl,debug_flag,lw_flag, & 
                               wavelength_lres,tmp_rad_lres_clr_ice,tmp_rad_lres_all_ice,solar_flux(j,:), &
                               tmp_bbupdif_clr_ice,tmp_bbdndif_clr_ice,tmp_bbdndir_clr_ice, & 
                               tmp_bbupdif_all_ice,tmp_bbdndif_all_ice,tmp_bbdndir_all_ice, & 
                               tmp_flxupdif_clr_ice,tmp_flxdndif_clr_ice,tmp_flxdndir_clr_ice, &
                               tmp_flxupdif_all_ice,tmp_flxdndif_all_ice,tmp_flxdndir_all_ice, &
                               wavelength_hres,tmp_rad_hres_clr_ice,tmp_rad_hres_all_ice,solar_zenith(j), &
                               tmp_difflux_clr_ice,tmp_difflux_all_ice)
               endif

            

               stop 'error radctl line 2361'
       endif

               !Assign brdf for ocean conditions
               do jj=1,15
                   brdf_wvl(jj) = 0.
                   do kk=1,3
                      brdf_param2(jj,kk) = 0.
                   end do
               end do

               brdf_len = 15 
               do jj=1,brdf_len
                  brdf_wvl(jj) = 0.4 + real(jj-1)/10.
                  brdf_param2(jj,1) = 0.4 + real(jj-1)/10.
                  brdf_param2(jj,2) = real_windspeed
                  brdf_param2(jj,3) = 0.
               end do

               ocean_flag = .true.
               fname =  fname//'ocean'
               fname_length = fname_length+5

               if (modtran_flag) then
               !if (j.eq.84 .and. lchnk.eq.118) then
               !if (j.eq.1 .and. lchnk.eq.1) then
               !if (j.eq.211 .and. lchnk.eq.62) then
               !if (j.eq.210 .and. lchnk.eq.132) then
		  !write(*,*) 'cloud liquid is ',cliqwp_local(j,:)
		  !write(*,*) 'cloud ice is ',cicewp_local(j,:)
               !if (j.eq.166 .and. lchnk.eq.166) then
               !if (j.lt.20 .and. lchnk.lt.20) then
!	       call run_modpcrtm
                  call run_modtran_or_pcrtm(lchnk,j,nconfig_v(j),wgtv_v(j,:),ccon_v(j,:,:),totwgt_v(j),cliqwp_local(j,:), &
                               cicewp_local(j,:),cld_local(j,:),rei(j,:),rel(j,:),fname, &
                               fname_length,surftemp_real,co2vmr_real,gndalt_real,jday_input,gmt_real,zen_real,psipo_real, &
                               a_flag,c_flag,full_atm,pnm(j,:),v1_real,v2_real,ocean_flag,brdf_len,brdf_wvl,brdf_param2, &
                               spec_emissivity_len,spec_emissivity_wvl,spec_emissivity_vals,&
                               num_levs,m_aerosol,o_aerosol,aer_wavelen,pfs,phase_wvls_real,phase_angles_real,wvl_cld, &
                               ext_liq_cld,ssa_liq_cld,asym_liq_cld,ext_ice_cld,ssa_ice_cld,asym_ice_cld,z_cld_val_real, &
                               dv_input_real,num_lres,num_hres,num_ac_wvl,debug_flag,lw_flag, & 
                               wavelength_lres,tmp_rad_lres_clr_noice,tmp_rad_lres_all_noice,solar_flux(j,:), &
                               tmp_bbupdif_clr_noice,tmp_bbdndif_clr_noice,tmp_bbdndir_clr_noice, & 
                               tmp_bbupdif_all_noice,tmp_bbdndif_all_noice,tmp_bbdndir_all_noice, & 
                               tmp_flxupdif_clr_noice,tmp_flxdndif_clr_noice,tmp_flxdndir_clr_noice, &
                               tmp_flxupdif_all_noice,tmp_flxdndif_all_noice,tmp_flxdndir_all_noice, &
                               wavelength_hres,tmp_rad_hres_clr_noice,tmp_rad_hres_all_noice,solar_zenith(j), &
                               tmp_difflux_clr_noice,tmp_difflux_all_noice)
               endif
               !stop 'error radctl line 2398'

               radiance_lres_clr(j,:) = tmp_rad_lres_clr_noice*(1.-real(icefrac(j))) + tmp_rad_lres_clr_ice*real(icefrac(j))
               radiance_lres_all(j,:) = tmp_rad_lres_all_noice*(1.-real(icefrac(j))) + tmp_rad_lres_all_ice*real(icefrac(j))
               radiance_hres_clr(j,:) = tmp_rad_hres_clr_noice*(1.-real(icefrac(j))) + tmp_rad_hres_clr_ice*real(icefrac(j))
               radiance_hres_all(j,:) = tmp_rad_hres_all_noice*(1.-real(icefrac(j))) + tmp_rad_hres_all_ice*real(icefrac(j))

               bb_updiffuse_clr(j,:) = tmp_bbupdif_clr_noice*(1.-real(icefrac(j))) + tmp_bbupdif_clr_ice*real(icefrac(j))
               bb_updiffuse_all(j,:) = tmp_bbupdif_all_noice*(1.-real(icefrac(j))) + tmp_bbupdif_all_ice*real(icefrac(j))
               bb_dndiffuse_clr(j,:) = tmp_bbdndif_clr_noice*(1.-real(icefrac(j))) + tmp_bbdndif_clr_ice*real(icefrac(j))
               bb_dndiffuse_all(j,:) = tmp_bbdndif_all_noice*(1.-real(icefrac(j))) + tmp_bbdndif_all_ice*real(icefrac(j))

               if (lchnk.eq.46 .and. j.eq.178) then
                  !write(*,*) 'bb_updiffuse_all problem = ',tmp_bbupdif_all_noice
                  !write(*,*) 'bb_dndiffuse_all problem = ',tmp_bbdndif_all_noice
                  !write(*,*) 'bb_dndirect_all problem = ',tmp_bbdndir_all_noice
               endif

               bb_dndirect_clr(j,:) = tmp_bbdndir_clr_noice*(1.-real(icefrac(j))) + tmp_bbdndir_clr_ice*real(icefrac(j))
               bb_dndirect_all(j,:) = tmp_bbdndir_all_noice*(1.-real(icefrac(j))) + tmp_bbdndir_all_ice*real(icefrac(j))

               flx_updiffuse_clr(j,:,:) = tmp_flxupdif_clr_noice*(1.-real(icefrac(j))) + tmp_flxupdif_clr_ice*real(icefrac(j))
               flx_updiffuse_all(j,:,:) = tmp_flxupdif_all_noice*(1.-real(icefrac(j))) + tmp_flxupdif_all_ice*real(icefrac(j))
               flx_dndiffuse_clr(j,:,:) = tmp_flxdndif_clr_noice*(1.-real(icefrac(j))) + tmp_flxdndif_clr_ice*real(icefrac(j))
               flx_dndiffuse_all(j,:,:) = tmp_flxdndif_all_noice*(1.-real(icefrac(j))) + tmp_flxdndif_all_ice*real(icefrac(j))

               flx_dndirect_clr(j,:,:) = tmp_flxdndir_clr_noice*(1.-real(icefrac(j))) + tmp_flxdndir_clr_ice*real(icefrac(j))
               flx_dndirect_all(j,:,:) = tmp_flxdndir_all_noice*(1.-real(icefrac(j))) + tmp_flxdndir_all_ice*real(icefrac(j))
               diffuse_flux_clr(j,:) = tmp_difflux_clr_noice*(1.-real(icefrac(j))) + tmp_difflux_clr_ice*real(icefrac(j))
               diffuse_flux_all(j,:) = tmp_difflux_all_noice*(1.-real(icefrac(j))) + tmp_difflux_all_ice*real(icefrac(j))

               !Deallocate arrays
               deallocate(tmp_rad_lres_clr_noice)
               deallocate(tmp_rad_lres_clr_ice)
               deallocate(tmp_rad_hres_clr_noice)
               deallocate(tmp_rad_hres_clr_ice)
               deallocate(tmp_rad_lres_all_noice)
               deallocate(tmp_rad_lres_all_ice)
               deallocate(tmp_rad_hres_all_noice)
               deallocate(tmp_rad_hres_all_ice)

               deallocate(tmp_bbupdif_clr_noice)
               deallocate(tmp_bbupdif_clr_ice)
               deallocate(tmp_bbupdif_all_noice)
               deallocate(tmp_bbupdif_all_ice)

               deallocate(tmp_bbdndif_clr_noice)
               deallocate(tmp_bbdndif_clr_ice)
               deallocate(tmp_bbdndif_all_noice)
               deallocate(tmp_bbdndif_all_ice)

               deallocate(tmp_bbdndir_clr_noice)
               deallocate(tmp_bbdndir_clr_ice)
               deallocate(tmp_bbdndir_all_noice)
               deallocate(tmp_bbdndir_all_ice)

               deallocate(tmp_flxupdif_clr_noice)
               deallocate(tmp_flxupdif_clr_ice)
               deallocate(tmp_flxupdif_all_noice)
               deallocate(tmp_flxupdif_all_ice)

               deallocate(tmp_flxdndif_clr_noice)
               deallocate(tmp_flxdndif_clr_ice)
               deallocate(tmp_flxdndif_all_noice)
               deallocate(tmp_flxdndif_all_ice)

               deallocate(tmp_flxdndir_clr_noice)
               deallocate(tmp_flxdndir_clr_ice)
               deallocate(tmp_flxdndir_all_noice)
               deallocate(tmp_flxdndir_all_ice)

               deallocate(tmp_difflux_clr_noice)
               deallocate(tmp_difflux_clr_ice)
               deallocate(tmp_difflux_all_noice)
               deallocate(tmp_difflux_all_ice)

               !fname_spectrum = 'modtran_spectrum_out.nc' 
               !call modtran_out_netcdf(wavelength_lres,wavelength_hres,&
               !                       radiance_lres_clr(j,:), radiance_hres_clr(j,:),&
               !                       solar_flux(j,:), solar_zenith(j), diffuse_flux_clr(j,:), &
               !                       bb_updiffuse_clr(j,:),bb_dndiffuse_clr(j,:),bb_dndirect_clr(j,:), &
               !                       flx_updiffuse_clr(j,:,:),flx_dndiffuse_clr(j,:,:),flx_dndirect_clr(j,:,:),fname_spectrum)

            else  !for ice_flag_local
               !write(*,*) 'ice-free'
               
               if (modtran_flag) then
               !if (j.eq.84 .and. lchnk.eq.118) then
               !if (j.eq.1 .and. lchnk.eq.1) then
               !if (j.eq.211 .and. lchnk.eq.62) then
               !if (j.eq.210 .and. lchnk.eq.132) then
	!	  write(*,*) 'cloud liquid is ',cliqwp_local(j,:)
	!	  write(*,*) 'cloud ice is ',cicewp_local(j,:)
               !if (j.lt.20 .and. lchnk.lt.20) then
               !if (j.eq.166 .and. lchnk.eq.166) then
!	       call run_modpcrtm
                 call run_modtran_or_pcrtm(lchnk,j,nconfig_v(j),wgtv_v(j,:),ccon_v(j,:,:),totwgt_v(j),cliqwp_local(j,:),&
                               cicewp_local(j,:),cld_local(j,:),rei(j,:),rel(j,:),fname, &
                               fname_length,surftemp_real,co2vmr_real,gndalt_real,jday_input,gmt_real,zen_real,psipo_real, &
                               a_flag,c_flag,full_atm,pnm(j,:),v1_real,v2_real,ocean_flag,brdf_len,brdf_wvl,brdf_param2, &
                               spec_emissivity_len,spec_emissivity_wvl,spec_emissivity_vals,&
                               num_levs,m_aerosol,o_aerosol,aer_wavelen,pfs,phase_wvls_real,phase_angles_real,wvl_cld, &
                               ext_liq_cld,ssa_liq_cld,asym_liq_cld,ext_ice_cld,ssa_ice_cld,asym_ice_cld,z_cld_val_real, &
                               dv_input_real,num_lres,num_hres,num_ac_wvl,debug_flag,lw_flag, & 
                               wavelength_lres,radiance_lres_clr(j,:),radiance_lres_all(j,:),solar_flux(j,:), &
                               bb_updiffuse_clr(j,:),bb_dndiffuse_clr(j,:),bb_dndirect_clr(j,:), &
                               bb_updiffuse_all(j,:),bb_dndiffuse_all(j,:),bb_dndirect_all(j,:), &
                              flx_updiffuse_clr(j,:,:),flx_dndiffuse_clr(j,:,:),flx_dndirect_clr(j,:,:), &
                               flx_updiffuse_all(j,:,:),flx_dndiffuse_all(j,:,:),flx_dndirect_all(j,:,:), &
                               wavelength_hres,radiance_hres_clr(j,:),radiance_hres_all(j,:),solar_zenith(j), &
                               diffuse_flux_clr(j,:),diffuse_flux_all(j,:))
               endif

               !lchnk is lat index, j is lon index 
               if (j.eq.101 .and. lchnk.eq.64) then  
                  !write(*,*) 'cloud liq at 64= ',cliqwp_local(j,:)
                  !write(*,*) 'cloud ice at 64= ',cicewp_local(j,:)
               endif


               !if (j.eq.43 .and. lchnk.eq.43) then
               !   write(*,*) 'bb_updiffuse_all problem = ',bb_updiffuse_all(j,:)
               !   write(*,*) 'bb_dndiffuse_all problem = ',bb_dndiffuse_all(j,:)
               !   write(*,*) 'bb_dndirect_all problem = ',bb_dndirect_all(j,:)
               !   write(*,*) 'radiance_hres_clr = ',radiance_hres_clr(j,:)
               !   write(*,*) 'radiance_hres_all = ',radiance_hres_all(j,:)
               !endif


               !write(*,*) 'diffuse_flux_clr = ',radiance_lres_clr(j,51)
               !write(*,*) 'diffuse_flux_all = ',radiance_lres_all(j,51)
               !stop 'error radctl line 2492'

            endif !end of ice_flag_local structure

         end do  ! of j loop

         !Deallocate cloud arrays
         !deallocate(out_cldfrac)
         !deallocate(out_cldliqprof)
         !deallocate(out_cldiceprof)
 
         !
         !  Output fluxes at 200 mb
         !
         call vertinterp(ncol, pcols, pverp, pint, 20000._r8, fln, fln200)
         call vertinterp(ncol, pcols, pverp, pint, 20000._r8, flnc, fln200c)
         !
         ! Convert units of longwave fields needed by rest of model from CGS to MKS
         !
         do i=1,ncol
            flnt(i)  = flnt(i)*1.e-3
            flut(i)  = flut(i)*1.e-3
            flutc(i) = flutc(i)*1.e-3
            flns(i)  = flns(i)*1.e-3
            flntc(i) = flntc(i)*1.e-3
            fln200(i)  = fln200(i)*1.e-3
            fln200c(i) = fln200c(i)*1.e-3
            flnsc(i) = flnsc(i)*1.e-3
            flwds(i) = flwds(i)*1.e-3
            lwcf(i)=flutc(i) - flut(i)
            swcf(i)=fsntoa(i) - fsntoac(i)
            do jj=1,pver
               qrl(i,jj) = 0.
            end do
         end do

         !write(*,*) 'flnt = ',flnt(1)

         !
         ! Dump longwave radiation information to history tape buffer (diagnostics)
         !
         call outfld('QRL     ',qrl(:ncol,:)/cpair,ncol,lchnk)
         call outfld('FLNT    ',flnt  ,pcols,lchnk)
         call outfld('FLUT    ',flut  ,pcols,lchnk)
         call outfld('FLUTC   ',flutc ,pcols,lchnk)
         call outfld('FLNTC   ',flntc ,pcols,lchnk)
         call outfld('FLNS    ',flns  ,pcols,lchnk)
         call outfld('FLNSC   ',flnsc ,pcols,lchnk)
         call outfld('LWCF    ',lwcf  ,pcols,lchnk)
         call outfld('SWCF    ',swcf  ,pcols,lchnk)
         call outfld('FLN200  ',fln200,pcols,lchnk)
         call outfld('FLN200C ',fln200c,pcols,lchnk)
      end if

      !   
      return
    end subroutine radctl
    
