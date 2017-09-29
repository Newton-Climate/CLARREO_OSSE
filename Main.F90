#include <misc.h>
#include <params.h>
#define PCWDETRAIN

program main


!---------------------------Code history--------------------------------
!
!
!Dynamically allocated the arrays for modtran, so we stay -CAA-6/22/09
!   below the 2GB boundary.                                  
!
!
! 
!-----------------------------------------------------------------------

  
  use shr_kind_mod,  only: r8 => shr_kind_r8    
  use shr_const_mod, only: SHR_CONST_PI
  use shr_orb_mod
  use ppgrid
  use pmgrid
  use physconst,     only: gravit, latvap, cpair, tmelt, rga
  use physconst,     only: stebol, epsilo, pstd, rh2o, latice
  use physics_types, only: physics_state
  use radae,         only: initialize_radbuffer
  use constituents,  only: pcnst, ppcnst
  use ghg_surfvals,  only: scenario_ghg, ghg_surfvals_init, &
                           ghg_surfvals_ramp, ghg_surfvals_set
  use constituents,  only: cnst_get_ind
  use Chunks
  use RadInput
  use Output
  use pkg_cldoptics, only: cldefr, cldems, cldovrlap
  use chemistry,     only: trace_gas
  use prescribed_aerosols,  only: radforce, strat_volcanic, &
                                  tauback, sulscl, carscl, &
                                  ssltscl, dustscl, volcscl, &
                                  bgscl_rf, sulscl_rf, carscl_rf, &
                                  ssltscl_rf, dustscl_rf, volcscl_rf, &
                                  scenario_carbon_scale, &
                                  scenario_prescribed_sulfur, &
                                  prescribed_sulfur, &
                                  rampyear_prescribed_sulfur, &
                                  aerosol_initialize, get_int_scales, &
                                  get_aerosol, aerint

  use aer_optics,    only: aer_optics_initialize
  use volcanicmass,  only: read_volcanic_mass
  use commap,        only: latdeg, londeg
  use filenames,     only: absems_data, aeroptics, bndtvaer, bndtvcarbonscale, &
                           bndtvghg, bndtvo, bndtvscon, bndtvvolc, ncdata,  &
                           brdf_file,ocean_refl,snow_refl,snow_brdf,phase_funcs,alb_file,alb_snow,snow_frac,emis_file
  use rgrid,         only: nlon
  use infnan,        only: bigint
  use time_manager,  only: mdbase, msbase, mbdate, mbsec, mdcur, &
                           mscur, mcdate, mcsec, calday
  use ramp_scon,     only: ramp_sconst, rampnl_scon
  use history, only: outfld !DRF
  use wv_saturation, only: aqsat !DRF
  use get_surface_properties, only: get_brdf,get_emis,emis_array,brdf_param,snow_reflectance, &
                                    ocean_reflectance,land_flag, &
                                    modis_asdir,modis_asdif,modis_aldir,modis_aldif,&
                                    snow_brdf_landtype,landtype,snow_frac,fsno_read
  use get_phase_functions, only: get_pfs,phase_functions,phase_rhs, &
                                 phase_wvls,phase_angles
  use get_previous_results, only: load_previous_results, &
                                  starting_from_scratch,plat_start

#if ( defined MODTRAN_SPMD )
  use our_wrap_mpi
#endif


  implicit none

#include <comctl.h>
#include <comsol.h>
#include <comhyb.h>
#include <comlun.h>
!
! Flags for reconstruction/input
!
  character(len=16) :: build_aermmr     ! Build AERMMR internally
                                        !   ALL  = build all aerosols from
                                        !          climatology.  Volcanics
                                        !          read from separate climo.
                                        !   NONE = read in all aerosols
                                        !   IPCC = read in sulfate
                                        !          use climatology for other
                                        !          scale carbon aerosols
                                        !          read in volcanics
  logical :: build_trace                ! Build CFCs,CH4,& N2O internally
  logical :: build_emis                 ! Build EMIS internally
  logical :: build_re                   ! Build RE internally
  logical :: build_ozone                ! Build O3 internally 
  logical :: sol_ann_mean               ! Set eccf = 1.0 for insolation
  logical :: no_o2_abs                  ! Turn off absorption by O2 in SW
  logical :: lw_flag                    ! True for LW, false for SW (DRF)

  character*16 scenario_scon
!                    ! values can be 'FIXED' or 'RAMPED'
!                    ! FIXED => scon is fixed and can either have preset or
!                    ! namelist value
!                    ! RAMPED => scon is ramped
!                    ! DEFAULT => FIXED
  integer rampYear_scon
!                    ! ramped scon fixed at this year if set to a value
!                    ! greater than zero.  Default value is 0.

!---------------------------Radiation Output----------------------------

#include <ptrrgrid.h>
  real(r8) :: qrs(plon,pver,plat)       ! Shortwave heating rate
  real(r8) :: qrl(plon,pver,plat)       ! Longwave  heating rate
  real(r8) :: flwds(plon,plat)          ! Surface longwave down flux
  real(r8) :: fsns(plon,plat)           ! Surface solar absorbed flux
  real(r8) :: fsnt(plon,plat)           ! Net column abs solar flux at model top
  real(r8) :: flns(plon,plat)           ! Srf longwave cooling (up-down) flux
  real(r8) :: flnt(plon,plat)           ! Net outgoing lw flux at model top
  real(r8) :: sols(plon,plat)           ! Direct beam solar rad. onto srf (sw)
  real(r8) :: soll(plon,plat)           ! Direct beam solar rad. onto srf (lw)
  real(r8) :: solsd(plon,plat)          ! Diffuse solar radiation onto srf (sw)
  real(r8) :: solld(plon,plat)          ! Diffuse solar radiation onto srf (lw)

  real(r8) :: FLNSC(plon,plat)          ! Net surface clear-sky longwave
  real(r8) :: FLNTC(plon,plat)          ! Net TOA clear-sky longwave
  real(r8) :: FSDS(plon,plat)           ! Downwelling surface shortwave
  real(r8) :: FSDSC(plon,plat)          ! Downwelling clear surface shortwave
  real(r8) :: FSNSC(plon,plat)          ! Net surface clear shortwave
  real(r8) :: FSNTC(plon,plat)          ! Net TOA cear shortwave
  real(r8) :: FSNIRTOA(plon,plat)       ! Solar NIR flux at TOA all-sky (DRF)
  real(r8) :: FSNIRTOAC(plon,plat)      ! Solar NIR flux at TOA clear-sky (DRF)
  real(r8) :: SOLIN(plon,plat)          ! TOA insolation
  real(r8) :: FLN200(plon,plat)         ! Net all-sky LW flux at 200 mb
  real(r8) :: FLN200C(plon,plat)        ! Net clear-sky LW flux at 200 mb
  real(r8) :: FSN200(plon,plat)         ! Net all-sky SW flux at 200 mb
  real(r8) :: FSN200C(plon,plat)        ! Net clear-sky SW flux at 200 mb
  real(r8) :: OASDIR(plon,plat)         ! Output Direct shortwave albedo (DRF)
  real(r8) :: OASDIF(plon,plat)         ! Output Diffuse shortwave albedo (DRF)
  real(r8) :: OALDIR(plon,plat)         ! Output Direct longwave albedo (DRF)
  real(r8) :: OALDIF(plon,plat)         ! Output Diffuse longwave albedo (DRF)

  real(r8) :: FLN(plon,pverp,plat)      ! Net longwave flux at each interface
  real(r8) :: FSN(plon,pverp,plat)      ! Net shortwave flux at each interface
  real(r8) :: aertau(plon,aertau_species,plat)  ! aerosol optical depths calculated by radcswmx.F90 (DRF)
  real(r8) :: o3vmr(pcols,pverr)        ! Ozone volume mixing ratio (DRF)
  real(r8) :: pbr(pcols,pverr)          ! Model mid-level pressures (dynes/cm2) (DRF)
  real(r8) :: pnm(pcols,pverrp)         ! Model interface pressures (dynes/cm2) (DRF)
  real(r8) :: eccf                      ! Earth/sun distance factor (DRF)
  real(r8) :: o3mmr(pcols,pverr)        ! Ozone mass mixing ratio (DRF)
  real(r8) :: esat(pcols,pverr)         ! saturation vapor pressure (DRF)
  real(r8) :: qsat(pcols,pverr)         ! saturation specific humidity (DRF)
  real(r8) :: rhq(pcols,pverr)          ! relative humidity from aqsat (DRF)
  real(r8) :: dummy(pcols)              !dummy variable (DRF)
  real(r8) :: cldtau_ice_sw(plon,plat)         !cloud optical depth from ice clouds from radcswmx.F90 (DRF)
  real(r8) :: cldtau_liq_sw(plon,plat)         !cloud optical depth from liquid clouds from radcswmx.F90 (DRF)
  real(r8) :: cldtau_lw(plon,plat)         !cloud optical depth from liquid clouds from radcswmx.F90 (DRF)
  
  !Allocate memory for the arrays for modtran at runtime
  real*4, allocatable :: wavelength_lres(:)     !Wavelength values from Modtran (DRF)
  real*4, allocatable :: radiance_lres_clr(:,:,:)   !Radiance values from Modtran (DRF) clear-sky
  real*4, allocatable :: radiance_lres_all(:,:,:)   !Radiance values from Modtran (DRF) all-sky
  real*4, allocatable :: wavelength_hres(:)     !Wavelength values from Modtran (DRF)
  real*4, allocatable :: radiance_hres_clr(:,:,:)   !Radiance values from Modtran (DRF) clear-sky
  real*4, allocatable :: radiance_hres_all(:,:,:)   !Radiance values from Modtran (DRF) all-sky
  real*4, allocatable :: solar_flux(:,:,:)      !TOA downwelling solar flux values from Modtran (DRF)
  real*4, allocatable :: diffuse_flux_clr(:,:,:)    !Diffuse shortwave flux from Modtran (DRF) clear-sky
  real*4, allocatable :: diffuse_flux_all(:,:,:)    !Diffuse shortwave flux from Modtran (DRF) all-sky
  real*4, allocatable :: solar_zenith(:,:)      !Solar zenith angle from Modtran (DRF)
  real*4, allocatable :: aod(:,:)         !Solar zenith angle from Modtran (DRF)
  real*4, allocatable :: tau_sulf(:,:)         !Solar zenith angle from Modtran (DRF)
  real*4, allocatable :: tau_dust(:,:)         !Solar zenith angle from Modtran (DRF)
  real*4, allocatable :: tau_soot(:,:)         !Solar zenith angle from Modtran (DRF)
  real*4, allocatable :: tau_sslt(:,:)         !Solar zenith angle from Modtran (DRF)
  real*4, allocatable :: bb_updiffuse_clr(:,:,:)   !layer boundary upward diffuse broad-band flux (DRF)
  real*4, allocatable :: bb_updiffuse_all(:,:,:)   !layer boundary upward diffuse broad-band flux (DRF)
  real*4, allocatable :: bb_dndiffuse_clr(:,:,:)   !layer boundary downward diffuse broad-band flux (DRF)
  real*4, allocatable :: bb_dndiffuse_all(:,:,:)   !layer boundary downward diffuse broad-band flux (DRF)
  real*4, allocatable :: bb_dndirect_clr(:,:,:)   !layer boundary downward direct broad-band flux (DRF)
  real*4, allocatable :: bb_dndirect_all(:,:,:)   !layer boundary downward direct broad-band flux (DRF)
  !BEGIN DRF FOR FLUXES
  real*4, allocatable :: flx_updiffuse_clr(:,:,:,:)   !layer boundary upward diffuse spectral flux (DRF)
  real*4, allocatable :: flx_updiffuse_all(:,:,:,:)   !layer boundary upward diffuse spectral flux (DRF)
  real*4, allocatable :: flx_dndiffuse_clr(:,:,:,:)   !layer boundary downward diffuse spectral flux (DRF)
  real*4, allocatable :: flx_dndiffuse_all(:,:,:,:)   !layer boundary downward diffuse spectral flux (DRF)
  real*4, allocatable :: flx_dndirect_clr(:,:,:,:)   !layer boundary downward direct spectral flux (DRF)
  real*4, allocatable :: flx_dndirect_all(:,:,:,:)   !layer boundary downward direct spectral flux (DRF)
  !END DRF FOR FLUXES

real*4, allocatable  ::  mpi_singles_block(:)
real*4, allocatable :: mpi_wvlnglres_block(:)    !Wavelength values from Modtran (DRF)
real*4, allocatable :: mpi_wvlnghres_block(:)    !Wavelength values from Modtran (DRF)
real*4, allocatable :: mpi_level_block(:)    !level values (fluxes) from Modtran (DRF)
real*4, allocatable :: mpi_flux_block(:)    !level values (fluxes) from Modtran (DRF)
!BEGIN ADD OF RADCSWMX
real(r8), allocatable :: mpi_doubles_block(:) !single values from radcswmx.F90
real(r8), allocatable :: mpi_level_block_double(:)   !level values from radcswmx.F90
!END ADD OF RADCSWMX



!
!---------------------------Local workspace-----------------------------
!

  integer  :: iargc                     ! Routine to count arguments
  integer  :: jarg                      ! Counter for arguments

  type(physics_state) :: state(plat)    ! state information
  
  real(r8) :: aermmr(pcols,pver,plat)   ! level aerosol mass mixing ratio
  real(r8) :: rh(pcols,pver,plat)       ! level relative humidity (fraction)

  real(r8) :: pmid(plon,pver,plat)      ! midpoint pressure (Pa) 
  real(r8) :: pdel(plon,pver,plat)      ! layer thickness (Pa)
  real(r8) :: lnpmid(plon,pver,plat)    ! ln(pmid)
  real(r8) :: pint(plon,pver+1,plat)    ! interface pressure (Pa)
  real(r8) :: lnpint(plon,pver+1,plat)  ! ln(pint)

  real(r8) :: pmidx10(plon,pver)        ! pressure vals * 10 for cloud_config (DRF)
  real(r8) :: pintx10(plon,pver+1)       ! pressure vals * 10 for cloud_config (DRF)
  integer :: nconfig_v(plon,plat)
  real*4 :: nconfig_real(plon,plat)
  real(r8) :: totwgt_v(plon,plat)

  integer :: ccon_v_cloudconfig(pcols,pver,15)
  real(r8) :: wgtv_v_cloudconfig(pcols,15)

  integer :: lchnk_cloudconfig
  integer :: ncol_cloudconfig
  real(r8) :: solin_cloudconfig(pcols)
  real(r8) :: qrs_cloudconfig(pcols)
  real(r8) :: fsns_cloudconfig(pcols)
  real(r8) :: fsnt_cloudconfig(pcols)
  real(r8) :: fsntc_cloudconfig(pcols)
  real(r8) :: fsntoa_cloudconfig(pcols)
  real(r8) :: fsntoac_cloudconfig(pcols)
  real(r8) :: fsds_cloudconfig(pcols)
  real(r8) :: fsdsc_cloudconfig(pcols)
  real(r8) :: fsnsc_cloudconfig(pcols)

  real(r8) :: fsnirtoa_cloudconfig(pcols)
  real(r8) :: fsnrtoac_cloudconfig(pcols)
  real(r8) :: fsnrtoaq_cloudconfig(pcols)
  real(r8) :: sols_cloudconfig(pcols)
  real(r8) :: soll_cloudconfig(pcols)
  real(r8) :: solsd_cloudconfig(pcols)
  real(r8) :: solld_cloudconfig(pcols)
  real(r8) :: aertau_cloudconfig(pcols,19,7)
  real(r8) :: aerssa_cloudconfig(pcols,19,7)
  real(r8) :: aerasm_cloudconfig(pcols,19,7)
  real(r8) :: aerfwd_cloudconfig(pcols,19,7)
  real(r8) :: fns_cloudconfig(pcols,pverp)
  real(r8) :: fcns_cloudconfig(pcols,pverp)
  real(r8) :: frcday_cloudconfig(pcols)
  real(r8) :: cldtau_ice_cloudconfig(pcols)
  real(r8) :: cldtau_liq_cloudconfig(pcols)

  real(r8) :: asdir_cloudconfig(pcols)
  real(r8) :: asdif_cloudconfig(pcols)
  real(r8) :: aldir_cloudconfig(pcols)
  real(r8) :: aldif_cloudconfig(pcols)
  logical cldovrlp

  real(r8) :: gicewp(plon,pver,plat)    ! grid-box average IWP
  real(r8) :: gliqwp(plon,pver,plat)    ! grid-box average LWP
  real(r8) :: cldice(plon,pver,plat)    ! grid-box average ice mixing ratio
  real(r8) :: cldliq(plon,pver,plat)    ! grid-box average liq. mixing ratio
  real(r8) :: ficemr(plon,pver,plat)    ! fraction of ice in condensate, by mr
  real(r8) :: cwp(plon,pver,plat)       ! total in-cloud condensate water path
  real(r8) :: lwp(plon,pver,plat)       ! liquid condensate in-cloud water path

  real(r8) :: clat(plon,plat)           ! current latitudes(radians)
  real(r8) :: clon(plon,plat)           ! current longitude (radians)

  real(r8) :: coszrs(plon,plat)         ! Cosine solar zenith angle
  integer ::  nslice                    ! number of time slices in file
  integer, pointer ::  nstep(:)         ! current time step
  !integer, pointer ::  date(:)          ! current date
  !integer, pointer ::  datesec(:)       ! current second in date
  real(r8) :: dtime                     ! length of time step (seconds)

  real(r8) :: qm1(plon,pver,ppcnst,plat)! Tracers
  real(r8) :: aerosol(plon,pver,naer_all,plat) ! aerosol mass mixing ratios
  real(r8) :: scales(naer_all)          ! scaling factors for aerosols


  real(r8) pmxrgn(plon,pverp,plat)      ! Maximum values of pressure for each
                                        !  maximally overlapped region.
                                        !  0->pmxrgn(i,1)=range of pressure for
                                        !  1st region,
                                        !  pmxrgn(i,1)->pmxrgn(i,2) for  
                                        !  2nd region, etc
  integer nmxrgn(plon,plat)             ! Number of maximally overlapped regions

  character(len=120) :: ipath1          ! Path to input data from 1st run.
  character(len=120) :: ipath2          ! Path to input data from 2nd run.
  character(len=120) :: opath           ! Path to output data.

  integer :: itime                      ! Time index
  integer :: j                          ! Latitude index
  integer :: k                          ! Index used for writing mpi flux vals
  integer :: ih2o                       ! Index for H2O
  integer :: in2o                       ! Index for N2O
  integer :: ich4                       ! Index for CH4
  integer :: if11                       ! Index for CFC11
  integer :: if12                       ! Index for CFC12
  integer :: iaer                       ! Counter for aerosol species
  integer :: kaer                       ! Index for aerosol species
  integer :: jday_min(12)               ! Julian day of the beginning of the month (DRF)
  integer :: jday_max(12)               ! Julian day of the end of the month (DRF)
  integer :: jday_mid(12)               ! Julian day of the middle of the month (DRF)
  integer :: month_selected             ! Month corresponding to calday (DRF)
  integer :: ii                         ! index for getting month_selected (DRF)
  integer :: jj                         ! index for modis albedos (DRF)

  integer :: iii			! index for zeroing out aerosols
  integer :: jjj			! index for zeroing out aerosols
  integer :: kkk			! index for zeroing out aerosols
  integer :: mmm			! index for zeroing out aerosols

  integer :: fld_option                 ! Option for field switching
                                        !   0 = no switch
                                        !   1 = swap temperatures (T and TS)
                                        !   2 = swap spec. humidity (Q)
                                        !   3 = swap clouds
                                        !   4 = swap albedos
                                        !   5 = perform swaps 1-4

  integer :: qr_option                  ! Option for outputting heating rates
                                        !   0 = no QRL/QRS on output
                                        !   1 = QRL/QRS on output


  integer :: usrsys                     ! Flag for timing

  logical :: log_print                  ! Flag for diagnostic output from orbital
  logical :: create_output_flag         ! DRF flag for whether of not output file should be created (for restart)

  real(r8) :: solar_zen                   ! Solar zenith angle (degrees) Produced from Wiscombe calculator
  real(r8) :: solar_az(plon,plat)         ! Solar azimuth angle (degrees) Produced from Wiscombe calculator
  real(r8) :: solar_hour                  ! Local solar hour used for Wiscombe calculator
  integer :: mcdate_local                 ! Year used for Wiscombe calculator
  real(r8) :: mcdate_dble                 ! Year used for call to get_brdf
  integer :: day_local                    ! Day used for Wiscombe calculator
  real(r8) :: dummy_soldia                ! Solar diameter from Wiscombe calculator (not used)
  real(r8) :: dummy_soldst                ! Solar distance from Wiscombe calculator (not used)
  real(r8) :: sun_calc_long               ! Longitude value (degrees) for Wiscombe calculator
  real(r8) :: pi_double
  real(r8) :: degrees_per_radian
  real(r8):: calday_local


!MPI-------------------------------------------------------
 
  integer :: ierr,tag,i,proc_number
  integer :: my_rank
  integer :: tot_tasks
  logical :: our_master_proc
  integer :: my_latitude_start
  integer :: my_latitude_end
  integer :: my_longitude_start
  integer :: my_longitude_end

  real(r8) ::  total_time_start, total_time_end
  real(r8) ::  r0_message_start, r0_message_end
  real(r8) ::  r1_message_start, r1_message_end
 
 
 
 
  ! code to be timed
  real(r8) :: max_cld  !DRF


!  real(r4) :: low_wavelength_block,high_wavelength_block

  integer :: total_lat,total_long  !lat=128,lon=256
  integer :: lat_chunk,long_chunk
  integer :: lat_start_offset,lon_start_offset !for calculating subsets of the whole globe

! for > 128 procs
! 1 processor has 
! only 1 latitude and a subset of the longitudes
  integer :: total_chunks
  integer :: my_chunk_number
  integer :: chunks_per_lat
  integer :: longs_per_chunk
  integer :: my_lat_number
  integer :: start_mult_longs_chunk



!NEW MPI VARIABLES--------------------------------------------------------


integer :: current_proc, proc_multiplier, procs_ten_percent   !for counting
real(r8) :: weight_total, weight_average, weight_current, weight_used 
real(r8) :: ice_weight, cloud_weight, ice_param, cloud_param
integer :: mpi_signal


!map_proc_to_lat(proc,1):= starting lat for proc
!map_proc_to_lat(proc,2):= ending lat for proc
!map_lat_to_long(lat,1):= starting long for lat
!map_lat_to_long(lat,2):= ending long for lat
integer, allocatable :: map_proc_to_lat(:,:) !every proc has the same copy
integer, allocatable :: map_lat_to_long(:,:) !this is unique for each proc

real*8, allocatable :: weight_table(:,:)






!--------------------------------------------------------------------------

  namelist /settings/ absems_data, aeroptics, bndtvaer, bndtvcarbonscale, &
                      bndtvghg, bndtvo, bndtvscon, bndtvvolc, ncdata, &
                      build_emis, build_re, scon, iyear_AD, scenario_scon, &
                      rampYear_scon, sol_ann_mean, no_o2_abs, &
                      build_trace, scenario_ghg, &
                      build_aermmr, &
                      radforce, tauback, sulscl, carscl, ssltscl, dustscl, &
                      volcscl, bgscl_rf, sulscl_rf, carscl_rf, ssltscl_rf, &
                      dustscl_rf, volcscl_rf, build_ozone, ozncyc, &
                      ipath1, ipath2, opath, fld_option, qr_option,brdf_file, &
                      ocean_refl,snow_refl,snow_brdf,phase_funcs,alb_file, &
                      alb_snow,snow_frac,emis_file

!======================================================================

! Settings:
!   1. Climate sensitivity
!      build_aermmr =  'ALL'
!      build_trace =   .TRUE.
!      build_emis =    .TRUE.
!      build_re =      .TRUE.
!      build_ozone =   .TRUE.
!      scenario_scon = 'FIXED'
!      scenario_ghg  = 'FIXED'
!      ozncyc        = .TRUE.
!      strat_volcanic= .FALSE.
!   2. Climate sensitivity
!      build_aermmr =  'IPCC'
!      build_trace =   .FALSE.
!      build_emis =    .TRUE.
!      build_re =      .TRUE.
!      build_ozone =   .TRUE.
!      scenario_scon = 'RAMPED'
!      rampYear_scon = 0 (allow constant to vary)
!      scenario_ghg  = 'FIXED'*
!      ozncyc        = .FALSE.
!      strat_volcanic= .TRUE. itime
!
!* -- since GHG loss mechanism is on, we read CH4,N2O, & CFCs from history,
!     and we input CO2 from the 1D time series now included
!     in history files
!
! Notes: 
! 1. Split build_aermmr into flags for each species?  Current
!    behavior for IPCC is complicated -- some data are read, some are
!    reconstructed, some are scaled.
! 2. build_ozone flag is 
!    misnamed -- build_ozone just interpolates using the ozone boundary
!    data rather than reading it from the history file.
! 3. GHGs -- need consistent nomenclature between "trace" gases and "ghg"
!    control variables.
! 4. Flag for solar forcing for consistency?

!
!----------------------------------------------------------------------
! INITIALIZE MPI
!----------------------------------------------------------------------
!
  mpi_signal=1
  tag=0
  my_rank = 0
  our_master_proc = .true.
  tot_tasks = 0
  my_longitude_start=0
  my_longitude_end=0

  !set these only
  total_lat=128
  total_long=256
  !total_long=64
  lat_start_offset=0
  lon_start_offset=0


  
if (my_rank == 0) then
total_time_start = our_mpi_wtime()
endif



#if ( defined MODTRAN_SPMD )


  call our_mpiinit()
  call our_mpicommrank( my_rank )
  call our_mpicommsize( tot_tasks )
!  PRINT *, "task_no is ", my_rank, " of total ", tot_tasks, " tasks"

#endif


  if (my_rank == 0) then
open( unit=60, file="rank0timing.txt", &
      position="append", action="write" )
endif

if (my_rank == 1) then
open( unit=61, file="rank1timing.txt", &
      position="append", action="write" )
endif








  if (my_rank == 0) then
    our_master_proc = .true.
  else
    our_master_proc = .false.
  endif
  
 
!  total_chunks = tot_tasks
!  my_chunk_number = my_rank + 1 

!  chunks_per_lat = total_chunks/total_lat
!  longs_per_chunk = total_long/chunks_per_lat

!  my_latitude_start = 1+lat_start_offset + ((my_chunk_number-1)/chunks_per_lat)
!  my_latitude_end = my_latitude_start

!  start_mult_longs_chunk = mod((my_chunk_number-1),chunks_per_lat)

  
!  my_longitude_start =  (start_mult_longs_chunk*longs_per_chunk)+1+lon_start_offset

!  my_longitude_end = (start_mult_longs_chunk+1)*longs_per_chunk+lon_start_offset











!PRINT *, "total tasks=",tot_tasks," lat_chunk=",lat_chunk

!PRINT *, "my rank is ", my_rank, &
!         " start lat=", my_latitude_start, &
!         " end lat=",my_latitude_end,   & 
!         " start long=", my_longitude_start,  &
!         " end long=",my_longitude_end    








!  PRINT *, "********************INIT_MPI*******************************************************************"
!  PRINT *, "my rank is ", my_rank, " my master_proc value is", our_master_proc
!  PRINT *, "my rank is ", my_rank, " my start lat is=", my_latitude_start,  &
!           " my end lat is=",my_latitude_end 
!  PRINT *, "********************INIT_MPI*******************************************************************"





!
!----------------------------------------------------------------------
! ALLOCATE MEMORY FOR MODTRAN ARRAYS
!----------------------------------------------------------------------
!
 
  allocate (wavelength_lres(wvlng))             !Wavelength values from Modtran (DRF)
  allocate (radiance_lres_clr(plon,wvlng,plat))     !Radiance values from Modtran (DRF) clear-sky
  allocate (radiance_lres_all(plon,wvlng,plat))     !Radiance values from Modtran (DRF) all-sky
  allocate (wavelength_hres(wvlng_hres))        !Wavelength values from Modtran (DRF)
  allocate (radiance_hres_clr(plon,wvlng_hres,plat))!Radiance values from Modtran (DRF) clear-sky
  allocate (radiance_hres_all(plon,wvlng_hres,plat))!Radiance values from Modtran (DRF) all-sky
  allocate (solar_flux(plon,wvlng,plat))        !TOA downwelling solar flux values from Modtran (DRF)
  allocate (diffuse_flux_clr(plon,wvlng,plat))      !Diffuse shortwave flux from Modtran (DRF) clear-sky
  allocate (diffuse_flux_all(plon,wvlng,plat))      !Diffuse shortwave flux from Modtran (DRF) all-sky
  allocate (solar_zenith(plon,plat))            !Solar zenith angle from Modtran (DRF)
  allocate (aod(plon,plat))                     !Solar zenith angle from Modtran (DRF)
  allocate (tau_sulf(plon,plat))                     !Solar zenith angle from Modtran (DRF)
  allocate (tau_dust(plon,plat))                     !Solar zenith angle from Modtran (DRF)
  allocate (tau_soot(plon,plat))                     !Solar zenith angle from Modtran (DRF)
  allocate (tau_sslt(plon,plat))                     !Solar zenith angle from Modtran (DRF)
  allocate (bb_updiffuse_clr(plon,pver,plat))   !layer boundary upward diffuse broad-band flux (DRF)
  allocate (bb_updiffuse_all(plon,pver,plat))   !layer boundary upward diffuse broad-band flux (DRF)
  allocate (bb_dndiffuse_clr(plon,pver,plat))   !layer boundary downward diffuse broad-band flux (DRF)
  allocate (bb_dndiffuse_all(plon,pver,plat))   !layer boundary downward diffuse broad-band flux (DRF)
  allocate (bb_dndirect_clr(plon,pver,plat))   !layer boundary downward direct broad-band flux (DRF)
  allocate (bb_dndirect_all(plon,pver,plat))   !layer boundary downward direct broad-band flux (DRF)
  !BEGIN DRF FOR FLUXES
  allocate (flx_updiffuse_clr(plon,pver,num_bands,plat))   !layer boundary upward diffuse broad-band flux (DRF)
  allocate (flx_updiffuse_all(plon,pver,num_bands,plat))   !layer boundary upward diffuse broad-band flux (DRF)
  allocate (flx_dndiffuse_clr(plon,pver,num_bands,plat))   !layer boundary downward diffuse broad-band flux (DRF)
  allocate (flx_dndiffuse_all(plon,pver,num_bands,plat))   !layer boundary downward diffuse broad-band flux (DRF)
  allocate (flx_dndirect_clr(plon,pver,num_bands,plat))   !layer boundary downward direct broad-band flux (DRF)
  allocate (flx_dndirect_all(plon,pver,num_bands,plat))   !layer boundary downward direct broad-band flux (DRF)
  !END DRF FOR FLUXES

  allocate (mpi_singles_block(7))
  allocate (mpi_wvlnglres_block(5*wvlng))    !Wavelength values from Modtran (DRF)
  allocate (mpi_wvlnghres_block(2*wvlng_hres))    !Wavelength values from Modtran (DRF)
  allocate (mpi_level_block(6*pver))     !flux values from Modtran (DRF)
  !BEGIN DRF FOR FLUXES
  allocate (mpi_flux_block(6*pver*2))     !flux values from Modtran (DRF) now includes just vis and nir
  !END DRF FOR FLUXES
!BEGIN ADD OF RADCSWMX
  allocate (mpi_doubles_block(30))        !flux values, aod vals from radcswmx.F90
  !write(*,*) "here in main"
  allocate (mpi_level_block_double(2*pverp+2*pver+aertau_species))     !heating rates, aod vals from radcswmx.F90
!END ADD OF RADCSWMX
!
!----------------------------------------------------------------------
! FIELD CONSTRUCTION FLAGS AND OTHER NAMELIST DEFAULTS
!----------------------------------------------------------------------
!
! DEFAULTS GIVEN BELOW ARE FOR IPCC
!

!>>> Cloud Properties <<<
  build_emis =   .TRUE.
  build_re =     .TRUE.

!>>> Solar <<<
  scon          = 1.360e6
  iyear_AD      = 1950

  scenario_scon = 'RAMPED'
  rampYear_scon = 0

  sol_ann_mean  = .FALSE.
  no_o2_abs     = .FALSE.

!>>> WMGHGs <<<
  build_trace =  .FALSE.
  scenario_ghg = 'RAMPED'    

!>>> Aerosols <<<
  build_aermmr =  'IPCC'

  radforce   = .false. 

  tauback = 0._r8
  sulscl  = 1._r8
  carscl  = 1._r8
  ssltscl = 1._r8
  dustscl = 1._r8
  volcscl = 1._r8
  
  bgscl_rf   = 0._r8
  sulscl_rf  = 0._r8
  carscl_rf  = 0._r8
  ssltscl_rf = 0._r8
  dustscl_rf = 0._r8
  volcscl_rf = 0._r8
  
!>>> Ozone <<<
  build_ozone =  .TRUE.
  ozncyc = .FALSE.

!
! Filenames
!
  ipath1 = "in1.nc"
  ipath2 = "in2.nc"
  opath  = "output.nc"

!
! Options for input selection and heating rate output
!
  fld_option = 0
  qr_option  = 0

!
!----------------------------------------------------------------------
! INPUT/OUTPUT
!----------------------------------------------------------------------

 
  lw_flag = .false.




  if (iargc() /= 3) then 
     write(*,*) "Usage: input_1.nc input_2.nc output.nc < settings.inp"
     stop
  endif


 





  jarg = 1
  call getarg(jarg, ipath1)
 
  jarg = jarg + 1
  call getarg(jarg, ipath2)

  jarg = jarg + 1
  call getarg(jarg, opath)



  OPEN (UNIT=73, FILE='/nobackupp2/nhnguye3/osse_inmcm4/clarreo_osse/trunk/NOAA_modtran_surface/set.inp')
  !read(*, nml=settings)
 


  



!if(my_rank==0) then
  read(73, nml=settings)
!endif






  if (masterproc) then
     log_print = .TRUE.
  else
     log_print = .FALSE.
  end if






!
! Input timing data (needed for solar ephemeris)
!

       

  call input_times(ipath1, nslice, nstep   ,dtime   ,mdbase  ,msbase  ,&
                   mbdate, mbsec)

  !DRF 2/22/16 call input_times(ipath1, nslice, nstep   ,dtime   ,mdbase  ,msbase  ,&
  !DRF 2/22/16                 mbdate, mbsec, date, datesec)
  !write (*,*) "Here", ipath1, nslice, nstep   ,dtime   ,mdbase  ,msbase  ,&
  !                 mbdate, mbsec, date, datesec
!
! Input magic landm * cos(lat) flag if reconstructing rel/rei
!
  if (build_re) then
     !DRF 2/22/16 call input_landmcos(ncdata)
     call input_landmcos(ipath1)
  endif

!
! Check to see if output file already exists
!     call starting_from_scratch(opath,create_output_flag)

  create_output_flag = .true.
!
! Create output file if necessary
!

if (our_master_proc) then
  if (create_output_flag) then
     call create_output(ipath1, opath, nslice,  qr_option)
     !DRF 2/20/16 call create_output(ipath1, opath, nslice, qr_option)
  endif
endif
!
!----------------------------------------------------------------------
! MODEL GRID
!----------------------------------------------------------------------
! 
! Input vertical layer information and reconstruct hybrid coefficients
!    (needed for vertical grid analysis for aerosols and top level for
!     LW calculations)

  call input_vert_grid(ipath1)
  call hycoef()

!
! Set begchunk to 1 and endchunk to plat to allocate room for 
!     absorptivity/emissivity arrays
!
! In this code, the only place where begchunk and endchunk are used
!     is in the allocation of the absorptivity/emissivity arrays
!
  begchunk = 1
  endchunk = plat
!
! Set up rgrid module for volcanic mass interpolation
!
  nlon = plon

!
!----------------------------------------------------------------------
! SOLAR FORCING
!----------------------------------------------------------------------
!
! Solar and orbital constants
!
  call shr_orb_params (iyear_AD, eccen , obliq , mvelp, obliqr, &
                       lambm0, mvelpp, log_print)

!
! Determine ramping logic
!
   if (scenario_scon == 'FIXED') then
      doRamp_scon = .false.
   else if (scenario_scon == 'RAMPED') then
      doRamp_scon = .true.
   else
      call endrun ('READ_NAMELIST: SCENARIO_SCON must be set to either FIXED or RAMPED')
   endif
!       
! Initialize namelist related scon info
!
   if (doRamp_scon) then
      call rampnl_scon( rampYear_scon )
   endif

!
!----------------------------------------------------------------------
! GHG FORCING
!----------------------------------------------------------------------
!
! Initialize trace gas indices
!  NOTE: We use trace gas provision to load either reconstructed
!  profiles or input profiles of CH4,N2O, and CFC11/12 into slots
!  in the qm1 array.  This code will not work correctly if trace_gas = false
!
  trace_gas   = .TRUE.          

  call initindx()
  call cnst_get_ind('N2O'  , in2o)
  call cnst_get_ind('CH4'  , ich4)
  call cnst_get_ind('CFC11', if11)
  call cnst_get_ind('CFC12', if12)
  call cnst_get_ind('Q'    , ih2o)

  call ghg_surfvals_init()

!
!----------------------------------------------------------------------
! OZONE FORCING
!----------------------------------------------------------------------
!

!
! Set up to read ozone if build_ozone = .TRUE.
!  

  if (build_ozone) then
     call wrap_open(bndtvo, 0, ncid_oz)
!    oznini is called once we know the date, below.
  endif

!
!----------------------------------------------------------------------
! RADIATION INITIALIZATION
!----------------------------------------------------------------------




!
! Initialize radiation
!
  call initialize_radbuffer()

  call radini  (gravit  ,cpair   ,epsilo  ,stebol  ,pstd*10.0)!

!
! Set up vapor tables needed for RH calculation for aerosol hygroscopic growth
!
  call esinti  (epsilo  ,latvap  ,latice  ,rh2o    ,cpair  , &
                tmelt   )   

!
! Various radiation settings
!
  dosw =       .TRUE.
  dolw =       .TRUE.
  indirect =   .FALSE.

!
! Establish frequency of abs/ems calculations. Logic from parse_namelist.F90
!

  iradae =     -12
  if (iradae < 0) iradae = nint((-iradae*3600.)/dtime)

!===================================================
! Radiation computations
!===================================================

  usrsys = 0
  call t_setoptionf(usrsys,1)
  call t_initializef()

  do itime = 1, nslice

     call input_data(ipath1, ipath2, fld_option, itime, &
                     build_aermmr, build_trace, build_emis, &
                     build_re, build_ozone)
     call transpose_input(build_aermmr, build_trace, build_emis, &
                          build_re, build_ozone)
     call dump_input_data()

  !write(*,*) 'itime = ',itime
  !write(*,*) 'dtime = ',dtime
  !write(*,*) 'mdbase = ',mdbase
  !write(*,*) 'msbase = ',msbase
  !write(*,*) 'mbdate = ',mbdate
  !write(*,*) 'mbsec = ',mbsec
     call calendr(nstep(itime),dtime   ,mdbase  ,msbase  ,mbdate  , &
                  mbsec   ,mdcur   ,mscur   ,mcdate  ,mcsec   , &
                  calday  )
     write(*,*) "Here cal", nstep(itime),dtime   ,mdbase  ,msbase  ,mbdate  , &
                  mbsec   ,mdcur   ,mscur   ,mcdate  ,mcsec   , &
                  calday  


!------------------------------------------------------------------------------------------------------------------------------
!Calculate the average weight for each processor


!user can set these parameters
!ice_param 
!cloud_param
!ice_weight
!cloud_weight

           mcdate_local = 2000 !hard-coded
           solar_hour = 13.5
	   !solar_hour = 23.3893 !Gaussian quadrature point 1 (of 7)
	   !solar_hour = 20.8984 !Gaussian quadrature point 2 (of 7)
	   !solar_hour = 16.8701 !Gaussian quadrature point 3 (of 7)
	   !solar_hour = 12.0000 !Gaussian quadrature point 4 (of 7)
	   !solar_hour = 7.1299 !Gaussian quadrature point 5 (of 7)
	   !solar_hour = 3.1016 !Gaussian quadrature point 6 (of 7)
	   !solar_hour = 0.6107 !Gaussian quadrature point 7 (of 7)

           jday_min(1) = 1
           jday_min(2) = 32
           jday_min(3) = 60
           jday_min(4) = 91
           jday_min(5) = 121
           jday_min(6) = 152
           jday_min(7) = 182
           jday_min(8) = 213
           jday_min(9) = 244
           jday_min(10) = 274
           jday_min(11) = 305
           jday_min(12) = 335
           jday_max(1) = 31
           jday_max(2) = 59
           jday_max(3) = 90
           jday_max(4) = 120
           jday_max(5) = 151
           jday_max(6) =181
           jday_max(7) = 212
           jday_max(8) = 243
           jday_max(9) = 273
           jday_max(10) = 304
           jday_max(11) = 334
           jday_max(12) = 365
           jday_mid(1) = 15.
           jday_mid(2) = 46.
           jday_mid(3) = 75.
           jday_mid(4) = 105.
           jday_mid(5) = 135.
           jday_mid(6) = 165.
           jday_mid(7) = 195.
           jday_mid(8) = 225.
           jday_mid(9) = 255.
           jday_mid(10) = 285.
           jday_mid(11) = 315.
           jday_mid(12) = 345.

           calday_local = calday
           if (calday_local.lt.0.) then
              calday_local = calday_local+365.
           endif
           do ii=1,12
              if (floor(calday_local).ge.jday_min(ii).and.floor(calday_local).le.jday_max(ii)) then
                 month_selected = ii
              endif
           end do
           day_local =jday_mid(month_selected)
           pi_double = 3.14159265358979
           degrees_per_radian = 180.

!DRF modify to include cloud configurations for load balancing (8/16/10)
do j = 1,plat
   call plevs0(pcols, pcols, pver, c_ps(:,j), &
                       pint(:,:,j), pmid(:,:,j), pdel(:,:,j))
   
   do ii=1,plon
        call sun_calc(mcdate_local,day_local,solar_hour,c_lat(j),sun_calc_long, &
                      solar_az(ii,j),solar_zen,dummy_soldia,dummy_soldst)
        if (solar_zen.gt.87. .and. solar_zen.lt.90.) then
             solar_zen = 87.
        endif
        !if (solar_zen.ge.90. .and. solar_zen.lt.93.) then
        if (solar_zen.ge.90.) then
             solar_zen = 93.
        endif
         coszrs(ii,j) = cos(solar_zen*pi_double/degrees_per_radian)

         if (lw_flag) then
            coszrs(ii,j) = 0.5
         endif 
   end do  
   call cldovrlap(j, pcols, pint(:,:,j), c_cloud(:,:,j), &
                nmxrgn(:,j), pmxrgn(:,:,j))

   pmidx10 = 10.*pmid(:,:,j)
   pintx10 = 10.*pint(:,:,j) 
   lchnk_cloudconfig = j
   ncol_cloudconfig = pcols

   do ii=1,plon
     asdir_cloudconfig(ii) = 0.0_r8
     asdif_cloudconfig(ii) = 0.0_r8
     aldir_cloudconfig(ii) = 0.0_r8
     aldif_cloudconfig(ii) = 0.0_r8
   end do

   if (build_re) then

       !write(*,*) 'cldefr landfrac = ',maxval(c_landfrac(:,j)),minval(c_landfrac(:,j))
       !write(*,*) 'cldefr c_t_cld = ',maxval(c_t_cld(:,:,j)),minval(c_t_cld(:,:,j))
       !write(*,*) 'cldefr c_ps = ',maxval(c_ps(:,j)),minval(c_ps(:,j))
       !write(*,*) 'cldefr pmid = ',maxval(pmid(:,:,j)),minval(pmid(:,:,j))
       !write(*,*) 'cldefr c_landmcos = ',maxval(c_landmcos(:,j)),minval(c_landmcos(:,j))
       !write(*,*) 'cldefr c_icefrac = ',maxval(c_icefrac(:,j)),minval(c_icefrac(:,j))
       !write(*,*) 'cldefr c_snowh = ',maxval(c_snowh(:,j)),minval(c_snowh(:,j))

       call cldefr(j,  pcols,    &
            c_landfrac(:,j), c_t_cld(:,:,j), c_rel(:,:,j), &
            c_rei(:,:,j), c_ps(:,j), pmid(:,:,j), c_landmcos(:,j), &
            c_icefrac(:,j), c_snowh(:,j))

       !write(*,*) 'cldefr c_rel = ',maxval(c_rel(:,:,j)),minval(c_rel(:,:,j))
       !write(*,*) 'cldefr c_rei = ',maxval(c_rei(:,:,j)),minval(c_rei(:,:,j))
   endif

   !write(*,*) 'radcswmx_inputs lchnk_cloudconfig = ',lchnk_cloudconfig
   !write(*,*) 'radcswmx_inputs ncol_cloudconfig = ',ncol_cloudconfig
   !write(*,*) 'radcswmx_inputs pintx10 = ',pintx10
   !write(*,*) 'radcswmx_inputs pmidx10 = ',pmidx10
   !write(*,*) 'radcswmx_inputs qm1 = ',maxval(qm1(:,:,:,j)),minval(qm1(:,:,:,j)),shape(qm1)
   !write(*,*) 'radcswmx_inputs q = ',maxval(c_q(:,:,j)),minval(c_q(:,:,j)),shape(c_q)
   !write(*,*) 'radcswmx_inputs c_rh = ',maxval(c_rh(:,:,j)),minval(c_rh(:,:,j)),shape(c_rh)
   !write(*,*) 'radcswmx_inputs o3mmr =',maxval(o3mmr),minval(o3mmr),shape(o3mmr)
   !write(*,*) 'radcswmx_inputs aerosol = ',maxval(aerosol(:,:,:,j)),minval(aerosol(:,:,:,j)),shape(aerosol)
   !write(*,*) 'radcswmx_inputs c_cloud = ',maxval(c_cloud(:,:,j)),minval(c_cloud(:,:,j)),shape(c_cloud)
   !write(*,*) 'pint = ',maxval(pint(:,:,j)),minval(pint(:,:,j)),shape(pint)
   !write(*,*) 'radcswmx_inputs c_icldiwp = ',maxval(c_icldiwp(:,:,j)),minval(c_icldiwp(:,:,j)),shape(c_icldiwp)
   !write(*,*) 'radcswmx_inputs c_icldlwp = ',maxval(c_icldlwp(:,:,j)),minval(c_icldlwp(:,:,j)),shape(c_icldlwp)
   !write(*,*) 'radcswmx_inputs c_rel = ',maxval(c_rel(:,:,j)),minval(c_rel(:,:,j)),shape(c_rel)
   !write(*,*) 'radcswmx_inputs c_rei = ',maxval(c_rei(:,:,j)),minval(c_rei(:,:,j)),shape(c_rei)
   !write(*,*) 'radcswmx_inputs coszrs = ',shape(coszrs)
   !write(*,*) 'radcswmx_inputs aldir = ',shape(aldir_cloudconfig)
   !write(*,*) 'radcswmx_inputs asdir = ',shape(asdir_cloudconfig)
   !write(*,*) 'radcswmx_inputs aldif = ',shape(aldif_cloudconfig)
   !write(*,*) 'radcswmx_inputs asdif = ',shape(asdif_cloudconfig)
   !write(*,*) 'nmxrgn = ',maxval(nmxrgn),minval(nmxrgn),nmxrgn

   !write(*,*) 'nmxrgn = ',nmxrgn(:,j)
   !write(*,*) 'pmxrgn = ',pmxrgn(:,j)


   cldovrlp = .true.
   call radcswmx_cldovrlp(lchnk_cloudconfig   ,ncol_cloudconfig    ,                        &
                    pintx10    ,pmidx10   ,c_q(:,:,j),c_rh(:,:,j),o3mmr   , &
                    aerosol(:,:,:,j),c_cloud(:,:,j),c_icldiwp(:,:,j),c_icldlwp(:,:,j),c_rel(:,:,j), &
                    c_rei(:,:,j),eccf    ,coszrs(:,j)  ,          &
                    scon    ,solin_cloudconfig, &
                    asdir_cloudconfig,asdif_cloudconfig,aldir_cloudconfig,aldif_cloudconfig,nmxrgn(:,j)  , &
                    pmxrgn(:,:,j)  ,qrs_cloudconfig,fsnt_cloudconfig,fsntc_cloudconfig,fsntoa_cloudconfig, &
                    fsntoac_cloudconfig ,fsnirtoa_cloudconfig,fsnrtoac_cloudconfig,fsnrtoaq_cloudconfig,fsns_cloudconfig, &
                    fsnsc_cloudconfig,fsdsc_cloudconfig,fsds_cloudconfig,sols_cloudconfig,soll_cloudconfig, &
                    solsd_cloudconfig,solld_cloudconfig,frcday_cloudconfig,            &
                    aertau_cloudconfig,aerssa_cloudconfig,aerasm_cloudconfig,aerfwd_cloudconfig,fns_cloudconfig, &
                    cldtau_ice_cloudconfig,cldtau_liq_cloudconfig, &
                    wgtv_v_cloudconfig,totwgt_v(:,j),ccon_v_cloudconfig,nconfig_v(:,j),         &
#ifndef OFFLINE
                    fcns_cloudconfig,cldovrlp    )
#else
                    fcns_cloudconfig    ,no_o2_abs,cldovrlp)
#endif
   !write(*,*) 'here at 1079'
   !stop

   do ii=1,plon
      nconfig_real(ii,j) = real(nconfig_v(ii,j))
      if (nconfig_real(ii,j)<1.) then
         nconfig_real(ii,j) = 1.
      endif
   end do
end do



!Initialize and allocate memory for weight_table 
weight_total=0.0_r8
weight_average=0.0_r8
allocate (weight_table(plon,plat))   
weight_table=0.0_r8


!Find the total weight
do j = 1, plat !run over latitudes

	do i = 1, plon  !run over longitudes
               
		        if (c_icefrac(i,j).ge.0.001) then
                      		ice_weight = 2.0_r8
                        else
				ice_weight = 1.0_r8
			end if 
                                 
                	weight_table(i,j) = 1.0_r8 * ice_weight * dble(nconfig_real(i,j))
			weight_total = weight_total + weight_table(i,j)
                                
end do
	end do


!The average weight is
weight_average = weight_total/dble(tot_tasks)

!----------------------------------------------------------------------------------

!Assign the starting and ending latitudes and longitudes for each processor
!(longitude,latitude)(i,j)
!assumes i,j both begin at 1


!Allocate memory for the maps
allocate (map_proc_to_lat(tot_tasks,2)) 
allocate (map_lat_to_long(plat,2)) 


!Initialize and set the starting values for proc 0
!current_proc runs from 0 to tot_tasks-1
current_proc = 0
weight_current = 0.0_r8
weight_used = 0.0_r8
proc_multiplier = 1
!note: floor(x) returns the largest integer < x
!this is only used when tot_tasks is GREATER THAN 10
procs_ten_percent = floor(dble(tot_tasks)/100.0_r8)


!NOTE: my_rank runs from 0 to tot_tasks-1
!while map_proc_to_lat runs from 1 to tot_tasks because array indices begin at 1
!so we need to add 1 to my_rank to map it to a proc in this array
!meaning we need map_proc_to_lat(my_rank+1,.) to access the info for my_rank

!the beginning lat point for proc 0
map_proc_to_lat(1,1)=1
!the beginning long point for the beginning lat point for proc 0
if (my_rank .eq. 0) map_lat_to_long(1,1)=1

call our_mpibarrier()
call flush(6)

outer_lat_loop: do j = 1, plat !run over latitudes

	do i = 1, plon  !run over longitudes


		!handle the situation if I am at the last processor
                if((current_proc .eq. (tot_tasks-1)) .and. (j .eq. plat)) then

                       	map_proc_to_lat(current_proc+1,2)=plat
                       	if (my_rank .eq. current_proc) map_lat_to_long(j,2)=plon

                     	exit outer_lat_loop

                end if
                        

                !add to the running weight counter
        	weight_current = weight_current + weight_table(i,j)
                weight_used = weight_used + weight_table(i,j)


                !we have exceeded the average weight, so assign the end points
                if ((weight_current .gt. weight_average) .and. (current_proc /= (tot_tasks-1))) then


                        !DEBUG INFO
                        !CHECK OF AVERAGE WEIGHT OF EACH PROCESSOR
                        if ((my_rank .eq. 0) .and. (current_proc .eq. 0)) then
 				PRINT *, "********** BEGIN LOOKING AT THE WEIGHT OF EACH PROCESSOR ***************"
                        	PRINT *, "total procs = ", tot_tasks," total weight = ", weight_total
                        	PRINT *, "my rank = ", current_proc, &
					 " my weight = ", weight_current-weight_table(i,j)
                        else if ((my_rank .eq. 0) .and. (current_proc /= tot_tasks-2)) then
				PRINT *, "my rank = ", current_proc, &
					 " my weight = ", weight_current-weight_table(i,j)
                        else if ((my_rank .eq. 0) .and. (current_proc .eq. tot_tasks-2)) then
                        	PRINT *, "my rank = ", current_proc, &
					 " my weight = ", weight_current-weight_table(i,j)
                        	PRINT *, "my rank = ", current_proc+1, &
					 " my weight = ", weight_total - (weight_used-weight_table(i,j))        
                        end if
                        call flush(6)

   
             		!assign ending latitude and longitude 
                        !don't take the current i,j
                        !careful!we may be at the first longitude
                        if(i .eq. 1) then
				map_proc_to_lat(current_proc+1,2)=j-1
        	                if (my_rank .eq. current_proc) map_lat_to_long(j,2)=plon
                        else
	                        map_proc_to_lat(current_proc+1,2)=j
        	                if (my_rank .eq. current_proc) map_lat_to_long(j,2)=i-1
                        end if


                        !bump up the current processor number and reset the running weight
                        !make sure to include the weight of the current point (i,j)
                                             
			current_proc = current_proc + 1
                        weight_current = 0.0_r8 + weight_table(i,j)
                        weight_used = weight_used
                        

                        !assign the current i,j as the starting values
                        !for the current proc
                        map_proc_to_lat(current_proc+1,1)=j
                        if (my_rank .eq. current_proc) map_lat_to_long(j,1)=i

                        !recalculate average weight
                        !for every 10 percent of the total procs
                        !do this when we have GREATER THAN 10 mpi tasks
                        if((tot_tasks>10) .and. (current_proc /= tot_tasks) .and. (proc_multiplier /= 0)) then
                       
                        	if((current_proc) .ge. proc_multiplier*procs_ten_percent) then

                                	!bump up to the next 10 percent of procs
                                        proc_multiplier = proc_multiplier + 1

                                         
                                        !find new average weight
                                        !note: current proc has been bumped up one but has not been assigned to yet
                                        !note: do not include the current point in the leftover weight
                                        weight_average = (weight_total - weight_used - weight_table(i,j))/ &
                                                         dble(tot_tasks-current_proc)

				end if
                       
                        end if !recalculate average weight
                                                                           
  		end if !assign end points and move to next proc


	end do


        !went on to the next latitude, we need to account for this
        !on the current lat strip, assign the last longitude for the end longitude
        
        !if (my_rank .eq. 0) then
        ! 	PRINT *, "proc no = ", current_proc, &
	!	" i,j = ", i," ",j
        !end if
        !call flush(6)
	
        if (my_rank .eq. current_proc) map_lat_to_long(j,2)=plon
       	!start at the next lat point and start from the beginning of the longitudes
       	!note: j has not increased yet so we need to do it for the mapping
       	if (my_rank .eq. current_proc) map_lat_to_long(j+1,1)=1


end do  outer_lat_loop


!deallocate the weight table
deallocate (weight_table)

!------------------------------------------------------------------------------------------------------------------
!debug the load balancing
call our_mpibarrier()
if (my_rank .eq. 0) PRINT *, "*******START DEBUG LOAD BALANCING*****************************************************"
call flush(6)

do current_proc = 0, tot_tasks-1

        
	if (my_rank .eq. 0) then

                PRINT *
		PRINT *, "my rank = ", current_proc
!                PRINT *, "start lat = ", map_proc_to_lat(current_proc+1,1) , &
!		           " end lat = ", map_proc_to_lat(current_proc+1,2)
		
        end if
        call flush(6)


    
        do i = map_proc_to_lat(current_proc+1,1), map_proc_to_lat(current_proc+1,2)
	
		 !Find my_longitude_start and my_longitude_end for this particular latitude
      		 !and send that information to proc 0 as well
        	 if(my_rank==current_proc .and. current_proc /= 0) then

                        call our_mpirecvint (mpi_signal, 1, 0, 1)
	      		my_longitude_start = map_lat_to_long(i,1)
      	      		my_longitude_end = map_lat_to_long(i,2)
              		call our_mpisendint (my_longitude_start, 1, 0, 2)
              		call our_mpisendint (my_longitude_end, 1, 0, 3)

                 else if (my_rank==0 .and. current_proc /= 0) then

                        call our_mpisendint (mpi_signal, 1, current_proc, 1) 
              		call our_mpirecvint (my_longitude_start, 1, current_proc, 2)
              		call our_mpirecvint (my_longitude_end, 1, current_proc, 3)

                       	PRINT *, "lat = ", i
               		PRINT *, "start long = ", my_longitude_start , &
				 " end long = ", my_longitude_end

                 else if (my_rank==0 .and. current_proc==0) then

			PRINT *, "lat = ", i
            		PRINT *, "start long = ", map_lat_to_long(i,1) , &
				 " end long = ", map_lat_to_long(i,2)

       		 end if
                 call flush(6)

         end do



         if(my_rank==0) then

	         PRINT *
        	 PRINT *,"---------------------------------------------------------------------------------------"

	 end if


        call flush(6)
    
end do


if (my_rank .eq. 0) PRINT *, "*********END DEBUG LOAD BALANCING*****************************************************"

!------------------------------------------------------------------------------------------------------------------





































     if (itime == 1) then

        plat_start = 1
!        if (.not.create_output_flag) then
!            call get_nc_varids(opath)  
!            call load_previous_results(opath,itime,FLN,FLNS,FLNSC,FLNT,FLNTC,   &
!                              FLN200,FLN200C,FSDS,FSDSC,FSN,FSNS,FSNSC,FSNT,FSNTC, &
!                              FSN200,FSN200C,QRL,QRS,SOLIN,SOLL,SOLLD,SOLS,SOLSD, &
!                              WAVELENGTH_LRES,RADIANCE_LRES_CLR,RADIANCE_LRES_ALL, &
!                              WAVELENGTH_HRES,RADIANCE_HRES_CLR,RADIANCE_HRES_ALL, &
!                              SOLAR_FLUX,DIFFUSE_FLUX_CLR,DIFFUSE_FLUX_ALL,SOLAR_ZENITH,AOD)            
!        endif

!
! Initialize aerosols.  This requires coordinates in input not
!   available until now.
!
        latdeg = c_lat
        do j = 1, plat
           londeg(1:plon,j) = c_lon(:)
        end do

        call aer_optics_initialize()
        mcdate_dble = dble(mcdate)
        call get_brdf(c_lat,c_lon,c_snowh,c_icefrac,mcdate_dble)
        call get_emis(c_lat,c_lon,c_snowh,c_icefrac,mcdate_dble)
        call get_pfs()
        !write(*,*) "pfs = ",phase_functions(:,1,1,1)

        if (build_aermmr /= 'IPCC') then
           scenario_carbon_scale   = 'FIXED'
           strat_volcanic = .false.
        else
           scenario_carbon_scale   = 'RAMPED'
           strat_volcanic = .false.
        endif

        if (build_aermmr /= 'NONE') then
           scenario_prescribed_sulfur = 'FIXED' 
           prescribed_sulfur          = 'direct' 
           rampyear_prescribed_sulfur = bigint
           call aerosol_initialize()
        endif


!
! Initialize ozone.  This requires timeing information not available
!    until now
!
        if (build_ozone) then
           call oznini()
        endif
     endif
  

!
     if (dosw .or. dolw) then
!
! doabsems logic is identical to that in advnce.F90
!
        doabsems = nstep(itime).eq.0 .or. iradae.eq.1 .or. &
                   (mod(nstep(itime)-1,iradae).eq.0 .and. nstep(itime).ne.1)
        doabsems = .true.
        
        write(*,*) nstep(itime),dtime   ,mdbase  ,msbase  ,mbdate  , &
                     mbsec   ,mdcur   ,mscur   ,mcdate  ,mcsec   , &
                     calday
!
! Read next ozone, aerosol, and volcanic mass if needed
!        
        if (build_ozone) then 
           call oznint ()
        endif

        if (build_aermmr /= 'NONE') then
           call aerint ()
           if(strat_volcanic) then
              call read_volcanic_mass
           end if
        endif

        if (doRamp_scon) call ramp_sconst
!
! Set chunk index to latitude.  This is needed to pick appropriate
!     slice of absorptivity/emissivity 3D arrays
!
! See radclwmx for location where chunking is invoked to identify abs/ems sections
!
!$OMP PARALLEL DO PRIVATE (J)





		PRINT *, " MYRANK === ", my_rank, &
                " STARTLAT === ", map_proc_to_lat(my_rank+1,1) , &
		           " ENDLAT === ", map_proc_to_lat(my_rank+1,2)

    
    
    



        !Get value for restarting calculations, first read in previous values
 
        do j= map_proc_to_lat(my_rank+1,1), map_proc_to_lat(my_rank+1,2)  !plat_start,plat   !DRF
           write(*,*) "--> itime, ilat = ",itime,j



		PRINT *, "***MYRANK === ", my_rank , " ***MYLAT === ", j, &
                " ***STARTLONG === ", map_lat_to_long(j,1), &
		" ***ENDLONG === ", map_lat_to_long(j,2)





! Geographic coordinates and zenith angle calculations
!
           clat(:,j) = c_lat(j) / (180.0 / SHR_CONST_PI)
           clon(:,j) = c_lon(:) / (180.0 / SHR_CONST_PI)

           call zenith (calday, clat(:,j), clon(:,j), coszrs(:,j), pcols)
           mcdate_local = floor(mcdate/1e4)
           solar_hour = 13.5
	   !solar_hour = 23.3893 !Gaussian quadrature point 1 (of 7)
	   !solar_hour = 20.8984 !Gaussian quadrature point 2 (of 7)
	   !solar_hour = 16.8701 !Gaussian quadrature point 3 (of 7)
	   !solar_hour = 12.0000 !Gaussian quadrature point 4 (of 7)
	   !solar_hour = 7.1299 !Gaussian quadrature point 5 (of 7)
	   !solar_hour = 3.1016 !Gaussian quadrature point 6 (of 7)
	   !solar_hour = 0.6107 !Gaussian quadrature point 7 (of 7)

           jday_min(1) = 1
           jday_min(2) = 32
           jday_min(3) = 60
           jday_min(4) = 91
           jday_min(5) = 121
           jday_min(6) = 152
           jday_min(7) = 182
           jday_min(8) = 213
           jday_min(9) = 244
           jday_min(10) = 274
           jday_min(11) = 305
           jday_min(12) = 335
           jday_max(1) = 31
           jday_max(2) = 59
           jday_max(3) = 90
           jday_max(4) = 120
           jday_max(5) = 151
           jday_max(6) =181
           jday_max(7) = 212
           jday_max(8) = 243
           jday_max(9) = 273
           jday_max(10) = 304
           jday_max(11) = 334
           jday_max(12) = 365
           jday_mid(1) = 15.
           jday_mid(2) = 46.
           jday_mid(3) = 75.
           jday_mid(4) = 105.
           jday_mid(5) = 135.
           jday_mid(6) = 165.
           jday_mid(7) = 195.
           jday_mid(8) = 225.
           jday_mid(9) = 255.
           jday_mid(10) = 285.
           jday_mid(11) = 315.
           jday_mid(12) = 345.

           calday_local = calday
           if (calday_local.lt.0.) then 
              calday_local = calday_local+365.
           endif
 
           do ii=1,12
              if (floor(calday_local).ge.jday_min(ii).and.floor(calday_local).le.jday_max(ii)) then
                 month_selected = ii
              endif
           end do

           do ii=1,plon
              do jj=1,plat
                 oasdir(ii,jj) = modis_asdir(month_selected,jj,ii)
                 oasdif(ii,jj) = modis_asdif(month_selected,jj,ii)
                 oaldir(ii,jj) = modis_aldir(month_selected,jj,ii)
                 oaldif(ii,jj) = modis_aldif(month_selected,jj,ii)
              end do
           end do

           day_local =jday_mid(month_selected)
           sun_calc_long = 0.
           pi_double = 3.14159265358979
           degrees_per_radian = 180.
           mcdate_local = 2000 !hard-coded
           do ii=1,plon
              call sun_calc(mcdate_local,day_local,solar_hour,c_lat(j),sun_calc_long, &
                            solar_az(ii,j),solar_zen,dummy_soldia,dummy_soldst)
              if (solar_zen.gt.87. .and. solar_zen.lt.90.) then
                 solar_zen = 87.
              endif
              !if (solar_zen.ge.90. .and. solar_zen.lt.93.) then
              if (solar_zen.ge.90.) then
                 solar_zen = 93.
              endif
              coszrs(ii,j) = cos(solar_zen*pi_double/degrees_per_radian)
              
              if (lw_flag) then
                  coszrs(ii,j) = 0.5
              endif
           end do

!
! Build pressure fields
!
           call plevs0(pcols, pcols, pver, c_ps(:,j), &
                       pint(:,:,j), pmid(:,:,j), pdel(:,:,j))
           lnpmid(:,:,j) = dlog(pmid(:,:,j))
           lnpint(:,:,j) = dlog(pint(:,:,j))
!
! Cloud particle size and fraction of ice
!
           
           if (build_re) then
              call cldefr(j,  pcols,    &
                   c_landfrac(:,j), c_t_cld(:,:,j), c_rel(:,:,j), &
                   c_rei(:,:,j), c_ps(:,j), pmid(:,:,j), c_landmcos(:,j), &
                   c_icefrac(:,j), c_snowh(:,j))
           endif
!
! Cloud emissivity.
!
           if (build_emis) then
!
! Reverse the code in param_cldoptics that computes icldiwp and icldlwp
!    to obtain the grid-box averaged condensate mixing ratios cldice & cldliq
!
              gicewp(:,:,j) = c_icldiwp(:,:,j) * max(0.01_r8, c_cloud(:,:,j))
              !gliqwp(:,:,j) = (c_icldlwp(:,:,j) - c_icldiwp(:,:,j)) * &
              !                max(0.01_r8, c_cloud(:,:,j)) 
              gliqwp(:,:,j) = (c_icldlwp(:,:,j)) * &
                              max(0.01_r8, c_cloud(:,:,j)) 
              cldice(:,:,j) = gicewp(:,:,j) / (pdel(:,:,j) / gravit * 1000.0) 
              cldliq(:,:,j) = gliqwp(:,:,j) / (pdel(:,:,j) / gravit * 1000.0) 



!
! Now compute the ice mixing ratio and total condensate path using expressions
!     from param_cldoptics
!
              ficemr(:,:,j) = cldice(:,:,j) / &
                   max(1.e-10_r8, (cldice(:,:,j)+cldliq(:,:,j)))
              !DRF 07/29/16 cwp(:,:,j) = c_icldlwp(:,:,j)
              cwp(:,:,j) = c_icldlwp(:,:,j)+c_icldiwp(:,:,j)

	      !if (j.eq.132) then
!		write(*,*) 'landfrac = ',c_landfrac(:,j)
!		write(*,*) 'icefrac = ',c_icefrac(:,j)
!		write(*,*) 't_cld = ',c_t_cld(210,:,j)
!		write(*,*) 'ps in main= ',c_ps(:,j)
!		write(*,*) 'c_snowh in main = ',c_snowh(:,j)
!		write(*,*) 'c_landmcos in main = ',c_landmcos(:,j)
!
!
!		write(*,*) 'cldliq size = ',size(c_icldlwp,1),size(c_icldlwp,2),size(c_icldlwp,3)
!		write(*,*) 'cldliq in Main = ',gliqwp(210,:,132)
!		write(*,*) 'gcldice in Main = ',gicewp(210,:,132)
!		write(*,*) 'cldice in Main = ',c_icldiwp(210,:,132)
!		write(*,*) 'cloud in Main = ',max(0.01_r8,c_cloud(210,:,132))
!		!write(*,*) 'pdel in Main = ',gravit,pdel(210,:,132)
!		!write(*,*) 'cwp = ',c_icldlwp(210,:,132)
!		write(*,*) 'rei = ',c_rei(210,:,132)
!		write(*,*) 'ficemr = ',ficemr(210,:,132)
!		endif

              call cldems(j, pcols, cwp(:,:,j), ficemr(:,:,j), &
                   c_rei(:,:,j), c_emis(:,:,j))
!	      if (j.eq.132) then
!		write(*,*) 'emis = ',c_emis(210,:,132)
!		endif
           endif
!
! Determine parameters for maximum/random overlap
!
           call cldovrlap(j, pcols, pint(:,:,j), c_cloud(:,:,j), &
                nmxrgn(:,j), pmxrgn(:,:,j))
!
! Reconstitute tracers
!
!
! Ramping ghg if appropriate.  Note -- just used for CO2 in IPCC runs.
!    In IPCC runs, build_trace = .false.
!
           if (ghg_surfvals_ramp()) call ghg_surfvals_set()

           if (build_trace) then
              call trcmix(j, pcols, pmid(:,:,j), &
                   qm1(:,:,in2o,j), qm1(:,:,ich4,j), &
                   qm1(:,:,if11,j), qm1(:,:,if12,j),clat(:,j))
           else
!
! Or transfer tracers
!
              qm1(:,:,in2o,j) = c_n2o(:,:,j)
              qm1(:,:,ich4,j) = c_ch4(:,:,j)
              qm1(:,:,if11,j) = c_cfc11(:,:,j)
              qm1(:,:,if12,j) = c_cfc12(:,:,j)
           endif
!
! Transfer humidity
!
           qm1(:,:,ih2o,j) = c_q(:,:,j)

!
! Aerosols: build aermmr fields if build_aermmr == 'ALL' or 'IPCC'
!
           if (build_aermmr /= 'NONE') then
              call get_int_scales(scales)

              call get_aerosol(j, pint(:,:,j), aerosol(:,:,:,j), scales)
           else
!
! Convert from mass to mixing ratio for everything except volcanics
!    Volcanic parameterizations expect mass per unit area in each layer
!
              do iaer = 1, naer_all
                 kaer = aerosol_index(iaer)
!
! The c_allaer is accumulated mass from the surface (k=pver) up to an arbitrary layer k'
!
                 aerosol(:,1:pver-1,kaer,j) = c_allaer(:,1:pver-1,j,kaer) - &
                                                c_allaer(:,2:pver,j,kaer)
                 aerosol(:,pver,kaer,j)     = c_allaer(:,pver,j,kaer)
!do j= my_latitude_start, my_latitude_end
! All aerosols besides volcanics should be expressed as mixing ratios
!
                 if (kaer /= idxVOLC) then
                    aerosol(:,:,kaer,j) = aerosol(:,:,kaer,j) / &
                                               (pdel(:,:,j) * rga)
                 endif
              end do
           endif

!           if (j .eq. 16) then
!              open(11,form='unformatted',file='aerosol_read.dat')
!              write(11) aerosol(:,:,:,j)
!              close(11)
!              stop
!           endif
           !stop
           if (build_aermmr == 'IPCC') then
!
! Determine index for sulfates
!
              kaer = minval(minloc(abs(aerosol_index - idxSUL)))
!
! The c_allaer is accumulated mass from the surface (k=pver) up to an arbitrary layer k'
!
              !DRF aerosol(:,1:pver-1,kaer,j) = c_allaer(:,1:pver-1,j,kaer) - &
              !DRF                                  c_allaer(:,2:pver,j,kaer)
              !DRF aerosol(:,pver,kaer,j)     = c_allaer(:,pver,j,kaer)


!
! All aerosols besides volcanics should be expressed as mixing ratios
!
              !DRF aerosol(:,:,kaer,j) = aerosol(:,:,kaer,j) / (pdel(:,:,j) * rga)
              aerosol(:,:,kaer,j) = c_allaer(:,:,j,kaer)
           endif
!              
! Complete radiation calculations
!
           state(j)%t           = c_t(:,:,j)
           state(j)%pmid        = pmid(:,:,j)
           state(j)%q(:,:,ih2o) = c_q(:,:,j)

!DRF !!!, LWP is positive           lwp(:,:,j) = (c_icldlwp(:,:,j) - c_icldiwp(:,:,j))
!           lwp(:,:,j) = (c_icldlwp(:,:,j) - c_icldiwp(:,:,j))
           lwp(:,:,j) = (c_icldlwp(:,:,j) )

        !if (build_aermmr == 'NONE') then
		do iii=1,plon
		   do jjj=1,pver
		      do kkk=1,naer_all
		         do mmm=1,plat
				aerosol(iii,jjj,kkk,mmm) = 0.0_r8
			 end do
		      end do
 		   end do
	         end do
        !endif

!!! BEGIN DRF CHANGES
!
!   VARIABLE MAP BETWEEN radctl.F90 and Main.F90
!   radctl.F90 variable      Main.F90 variable
!   lchnk                    j
!   pcols                     pcols
!   pmid                     pmid(:,:,j)
!   pint                     pint(:,:,j)
!   o3vmr                    not defined
!   pbr                      not defined
!   pnm                      not defined
!   eccf                     not defined
!   o3mmr                    not defined
!   calday                   calday

!
! Interpolate ozone volume mixing ratio to model levels
!
   if (build_ozone) then
      call radozn(j   ,pcols    ,pmid(:,:,j)    ,o3vmr   )


      ! Set chunk dependent radiation input
      !
      call radinp(j   ,pcols,                                &
               pmid(:,:,j)    ,pint(:,:,j)    ,o3vmr   , pbr     ,&
#ifndef OFFLINE
               pnm     ,eccf    ,o3mmr   )
#else
               pnm     ,eccf    ,o3mmr   , calday)


      do ii=1,pcols
        do jj=1,pverr
           c_o3vmr(ii,jj,j) = o3vmr(ii,jj)
        end do
      end do

      if (sol_ann_mean) then
         eccf = 1.0_r8
      endif
#endif
    else
      ! Set chunk dependent radiation input
      !
      call radinp(j   ,pcols,                                &
               pmid(:,:,j)    ,pint(:,:,j)    ,c_o3vmr(:,:,j)   , pbr     ,&
#ifndef OFFLINE
               pnm     ,eccf    ,o3mmr   )
#else
               pnm     ,eccf    ,o3mmr   , calday)

      if (sol_ann_mean) then
         eccf = 1.0_r8
      endif
#endif


   endif
   call outfld('O3VMR   ',o3vmr ,pcols, j)


!Get saturation water vapor pressure and relative humidity (DRF)
      call aqsat(state(j)%t, state(j)%pmid, esat, qsat, pcols, &
                 pcols, pver, 1, pver)

      ! calculate relative humidity
      rhq(1:pcols,1:pver) = state(j)%q(1:pcols,1:pver,1) / qsat(1:pcols,1:pver) * &
         ((1.0 - epsilo) * qsat(1:pcols,1:pver) + epsilo) / &
         ((1.0 - epsilo) * state(j)%q(1:pcols,1:pver,1) + epsilo)

      !write(*,*) 'max cicldiwp in Main = ',maxval(c_icldiwp(:,:,j))
      !write(*,*) 'max lwp in Main = ',maxval(lwp(:,:,j))
      !write(*,*) 'coszrs = ',coszrs(:,j)

      call radctl (j, pcols, in2o, ich4, if11,if12, &
                c_lwup(:,j), c_ts(:,j),c_emis(:,:,j), pmid(:,:,j), &
                pint(:,:,j), lnpmid(:,:,j), lnpint(:,:,j), c_t(:,:,j), &
                qm1(:,:,:,j), c_cloud(:,:,j), c_icldiwp(:,:,j), &
                lwp(:,:,j), coszrs(:,j), day_local,c_asdir(:,j), c_asdif(:,j), &
                c_aldir(:,j), c_aldif(:,j), pmxrgn(:,:,j), nmxrgn(:,j), &
                fsns(:,j), fsnt(:,j), fsnirtoa(:,j), flns(:,j), flnt(:,j), &
                qrs(:,:,j), qrl(:,:,j), flwds(:,j), c_rel(:,:,j), &
                c_rei(:,:,j), sols(:,j), soll(:,j), solsd(:,j), solld(:,j), &
                build_ozone, state(j), aerosol(:,:,:,j), calday, c_o3vmr(:,:,j), &
                flnsc(:,j), flntc(:,j), fsds(:,j), fsdsc(:,j), fsnsc(:,j), &
                fsntc(:,j), fsnirtoac(:,j),solin(:,j), fln200(:,j), fln200c(:,j), &
                fsn200(:,j), fsn200c(:,j), fln(:,:,j), fsn(:,:,j), aertau(:,:,j), &
                sol_ann_mean, no_o2_abs, pbr,pnm,eccf,o3mmr,esat, qsat,c_rh(:,:,j), &
                c_lat(j),c_lon,co2vmr, wavelength_lres, radiance_lres_clr(:,:,j),radiance_lres_all(:,:,j), &
                wavelength_hres, radiance_hres_clr(:,:,j),radiance_hres_all(:,:,j), solar_flux(:,:,j), &
                land_flag(:,j), c_icefrac(:,j),fsno_read(:,j),brdf_param(:,:,:,j), emis_array(:,:,j), &
                ocean_reflectance,snow_reflectance, solar_az(:,j), &
                modis_asdir(month_selected,j,:),modis_asdif(month_selected,j,:), &
                modis_aldir(month_selected,j,:),modis_aldif(month_selected,j,:), &
                phase_functions,phase_rhs,phase_wvls,phase_angles, &
                diffuse_flux_clr(:,:,j),diffuse_flux_all(:,:,j),solar_zenith(:,j),ipath1,mcdate, &
                aod(:,j), tau_sulf(:,j),tau_dust(:,j),tau_soot(:,j),tau_sslt(:,j), &
                cldtau_ice_sw(:,j),cldtau_liq_sw(:,j), cldtau_lw(:,j),&
                bb_updiffuse_clr(:,:,j),bb_updiffuse_all(:,:,j), &
                bb_dndiffuse_clr(:,:,j),bb_dndiffuse_all(:,:,j), &
                bb_dndirect_clr(:,:,j), bb_dndirect_all(:,:,j), &             
                flx_updiffuse_clr(:,:,:,j),flx_updiffuse_all(:,:,:,j), &
                flx_dndiffuse_clr(:,:,:,j),flx_dndiffuse_all(:,:,:,j), &
                flx_dndirect_clr(:,:,:,j), flx_dndirect_all(:,:,:,j), &             
                lw_flag,c_windspeed(:,j),j,map_lat_to_long(j,1),map_lat_to_long(j,2),my_rank)


         max_cld = maxval(maxval(cldtau_liq_sw,1))
         if (max_cld>0.0_r8) then
            write(*,*) "cldtau = ",max_cld
         endif



           qrl(:,:,j) = qrl(:,:,j) / cpair
           qrs(:,:,j) = qrs(:,:,j) / cpair

           fln(:,:,j) = fln(:,:,j) * 1.e-3
           fsn(:,:,j) = fsn(:,:,j) * 1.e-3
  
           !write(*,*) "starting writing of output files"
        end do ! end of loop over latitudes
     end if


!CAA no messages
!if(1<0) then


if (my_rank == 0) then
r0_message_start = our_mpi_wtime()
endif






  !loop over all the processors besides processor 0
  do proc_number=1,(tot_tasks-1)



        !All processors besides proc 0 wait here to begin until they are signaled by proc 0
	if(my_rank/=0 .and. proc_number==1) then

        	 call our_mpirecvint (mpi_signal, 1, 0, 1000)
	      	
	else if (my_rank==0) then

        	 call our_mpisendint (mpi_signal, 1, proc_number, 1000)

        end if 
              		





















    if(my_rank==proc_number .or. my_rank==0)then

if (my_rank == 1) then
r1_message_start = our_mpi_wtime()
endif



!      if(my_rank==0)then
!        !my_chunk_number = my_rank + 1 
!        my_chunk_number = proc_number + 1 
!        my_latitude_start = 1 + lat_start_offset + ((my_chunk_number-1)/chunks_per_lat)
!        my_latitude_end = my_latitude_start
!        start_mult_longs_chunk = mod((my_chunk_number-1),chunks_per_lat)
!        my_longitude_start =  (start_mult_longs_chunk*longs_per_chunk)+1+lon_start_offset
!        my_longitude_end = (start_mult_longs_chunk+1)*longs_per_chunk+lon_start_offset
!      endif

     
      
      do j=map_proc_to_lat(proc_number+1,1), map_proc_to_lat(proc_number+1,2)

        
        !Find my_longitude_start and my_longitude_end for this particular latitude
        !and send that information to proc 0 as well
        if(my_rank==proc_number) then

	      my_longitude_start = map_lat_to_long(j,1)
      	      my_longitude_end = map_lat_to_long(j,2)
              call our_mpisendint (my_longitude_start, 1, 0, 1001)
              call our_mpisendint (my_longitude_end, 1, 0, 1002)

               PRINT *, "THEPROC=== ", my_rank , " THELAT=== ", j, &
                " THESTARTLONG === ", map_lat_to_long(j,1), &
		" THEENDLONG === ", map_lat_to_long(j,2)

        else

              call our_mpirecvint (my_longitude_start, 1, proc_number, 1001)
              call our_mpirecvint (my_longitude_end, 1, proc_number, 1002)	
              PRINT *, "000PROC=== ", my_rank , " THELAT=== ", j, &
                " THESTARTLONG === ", my_longitude_start, &
		" THEENDLONG === ", my_longitude_end

        end if
      

!my_longitude_start=1
!my_longitude_end=256



        do i=my_longitude_start,my_longitude_end

          if(my_rank==proc_number)then

            !PACK AND SEND

            mpi_singles_block(1)=solar_zenith(i,j)
            mpi_singles_block(2)=aod(i,j)
            mpi_singles_block(3)=tau_sulf(i,j)
            mpi_singles_block(4)=tau_dust(i,j)
            mpi_singles_block(5)=tau_soot(i,j)
            mpi_singles_block(6)=tau_sslt(i,j)
            mpi_singles_block(7)=nconfig_real(i,j)
            tag=1

            call our_mpisendreal (mpi_singles_block, 7, 0, tag)

            mpi_wvlnghres_block(1:wvlng_hres)=radiance_hres_clr(i,:,j)
            mpi_wvlnghres_block(wvlng_hres+1:2*wvlng_hres)=radiance_hres_all(i,:,j)
            tag=2

            call our_mpisendreal (mpi_wvlnghres_block, 2*wvlng_hres, 0, tag)

            mpi_wvlnglres_block(1:wvlng)=radiance_lres_clr(i,:,j)
            mpi_wvlnglres_block(wvlng+1:2*wvlng)=radiance_lres_all(i,:,j)
            mpi_wvlnglres_block((2*wvlng)+1:3*wvlng)=diffuse_flux_clr(i,:,j)
            mpi_wvlnglres_block((3*wvlng)+1:4*wvlng)=diffuse_flux_all(i,:,j)
            mpi_wvlnglres_block((4*wvlng)+1:5*wvlng)=solar_flux(i,:,j)
            tag=3

            call our_mpisendreal (mpi_wvlnglres_block, 5*wvlng, 0, tag)

            mpi_level_block(1:pver) = bb_updiffuse_clr(i,:,j)
            mpi_level_block(1+pver:2*pver) = bb_updiffuse_all(i,:,j)
            mpi_level_block(1+2*pver:3*pver) = bb_dndiffuse_clr(i,:,j)
            mpi_level_block(1+3*pver:4*pver) = bb_dndiffuse_all(i,:,j)
            mpi_level_block(1+4*pver:5*pver) = bb_dndirect_clr(i,:,j)
            mpi_level_block(1+5*pver:6*pver) = bb_dndirect_all(i,:,j)
            tag=4
            call our_mpisendreal (mpi_level_block, 6*pver, 0, tag)


            !BEGIN ADD OF RADCSWMX
            mpi_doubles_block(1)=flns(i,j)
            mpi_doubles_block(2)=flnsc(i,j)
            mpi_doubles_block(3)=flnt(i,j)
            mpi_doubles_block(4)=flntc(i,j)
            mpi_doubles_block(5)=fln200(i,j)
            mpi_doubles_block(6)=fln200c(i,j)
            mpi_doubles_block(7)=fsds(i,j)
            mpi_doubles_block(8)=fsdsc(i,j)
            mpi_doubles_block(9)=fsns(i,j)
            mpi_doubles_block(10)=fsnsc(i,j)
            mpi_doubles_block(11)=fsnt(i,j)
            mpi_doubles_block(12)=fsntc(i,j)
            mpi_doubles_block(13)=fsn200(i,j)
            mpi_doubles_block(14)=fsn200c(i,j)
            mpi_doubles_block(15)=fsnirtoa(i,j)
            mpi_doubles_block(16)=fsnirtoac(i,j)
            mpi_doubles_block(17)=solin(i,j)
            mpi_doubles_block(18)=soll(i,j)
            mpi_doubles_block(19)=solld(i,j)
            mpi_doubles_block(20)=sols(i,j)
            mpi_doubles_block(21)=solsd(i,j)
            mpi_doubles_block(22)=cldtau_liq_sw(i,j)
            mpi_doubles_block(23)=cldtau_ice_sw(i,j)
            mpi_doubles_block(24)=cldtau_lw(i,j)
            mpi_doubles_block(25)=oasdir(i,j) 
            mpi_doubles_block(26)=oasdif(i,j) 
            mpi_doubles_block(27)=oaldir(i,j) 
            mpi_doubles_block(28)=oaldif(i,j) 
            mpi_doubles_block(29)=flwds(i,j)
       
            tag=5
            call our_mpisenddouble(mpi_doubles_block, 29, 0, tag)

            mpi_level_block_double(1:pverp)=fsn(i,:,j)
            mpi_level_block_double(1+pverp:2*pverp)=fln(i,:,j)
            mpi_level_block_double(1+2*pverp:2*pverp+pver)=qrl(i,:,j)
            mpi_level_block_double(1+2*pverp+pver:2*pverp+2*pver)=qrs(i,:,j)
            mpi_level_block_double(1+2*pverp+2*pver:aertau_species+2*pverp+2*pver)=aertau(i,:,j)
            tag=6
            call our_mpisenddouble (mpi_level_block_double, 2*pver+2*pverp+aertau_species, 0, tag)
            !END ADD OF RADCSWMX


            mpi_flux_block(1:pver) = flx_updiffuse_clr(i,:,1,j)
            mpi_flux_block(1+pver:2*pver) = flx_updiffuse_clr(i,:,2,j)
            mpi_flux_block(1+2*pver:3*pver) = flx_updiffuse_all(i,:,1,j)
            mpi_flux_block(1+3*pver:4*pver) = flx_updiffuse_all(i,:,2,j)
            mpi_flux_block(1+4*pver:5*pver) = flx_dndiffuse_clr(i,:,1,j)
            mpi_flux_block(1+5*pver:6*pver) = flx_dndiffuse_clr(i,:,2,j)
            mpi_flux_block(1+6*pver:7*pver) = flx_dndiffuse_all(i,:,1,j)
            mpi_flux_block(1+7*pver:8*pver) = flx_dndiffuse_all(i,:,2,j)
            mpi_flux_block(1+8*pver:9*pver) = flx_dndirect_clr(i,:,1,j)
            mpi_flux_block(1+9*pver:10*pver) = flx_dndirect_clr(i,:,2,j)
            mpi_flux_block(1+10*pver:11*pver) = flx_dndirect_all(i,:,1,j)
            mpi_flux_block(1+11*pver:12*pver) = flx_dndirect_all(i,:,2,j)
            tag=7
            call our_mpisendreal (mpi_flux_block, 6*2*pver, 0, tag)

          endif !rank proc_number

          if(my_rank==0)then


            !RECEIVE AND UNPACK

            tag=1  
            call our_mpirecvreal (mpi_singles_block, 7, proc_number, tag)

            solar_zenith(i,j)=mpi_singles_block(1)
            aod(i,j)=mpi_singles_block(2)
            tau_sulf(i,j)=mpi_singles_block(3)
            tau_dust(i,j)=mpi_singles_block(4)
            tau_soot(i,j)=mpi_singles_block(5)
            tau_sslt(i,j)=mpi_singles_block(6)
            nconfig_real(i,j)=mpi_singles_block(7)
            
            tag=2
            call our_mpirecvreal (mpi_wvlnghres_block, 2*wvlng_hres, proc_number, tag)

            radiance_hres_clr(i,:,j)=mpi_wvlnghres_block(1:wvlng_hres)
            radiance_hres_all(i,:,j)=mpi_wvlnghres_block(wvlng_hres+1:2*wvlng_hres)

            tag=3
            call our_mpirecvreal (mpi_wvlnglres_block, 5*wvlng, proc_number, tag)

            radiance_lres_clr(i,:,j)=mpi_wvlnglres_block(1:wvlng)
            radiance_lres_all(i,:,j)=mpi_wvlnglres_block(wvlng+1:2*wvlng)
            diffuse_flux_clr(i,:,j)=mpi_wvlnglres_block((2*wvlng)+1:3*wvlng)
            diffuse_flux_all(i,:,j)=mpi_wvlnglres_block((3*wvlng)+1:4*wvlng)
            solar_flux(i,:,j)=mpi_wvlnglres_block((4*wvlng)+1:5*wvlng)

            tag=4
            call our_mpirecvreal (mpi_level_block, 6*pver, proc_number, tag)
            bb_updiffuse_clr(i,:,j) = mpi_level_block(1:pver)
            bb_updiffuse_all(i,:,j) = mpi_level_block(1+pver:2*pver)
            bb_dndiffuse_clr(i,:,j) = mpi_level_block(1+2*pver:3*pver)
            bb_dndiffuse_all(i,:,j) = mpi_level_block(1+3*pver:4*pver)
            bb_dndirect_clr(i,:,j) = mpi_level_block(1+4*pver:5*pver)
            bb_dndirect_all(i,:,j) = mpi_level_block(1+5*pver:6*pver)

          
            !BEGIN ADD OF RADCSWMX
            tag=5  
            call our_mpirecvdouble (mpi_doubles_block, 30, proc_number, tag)

            flns(i,j) = mpi_doubles_block(1)
            flnsc(i,j) = mpi_doubles_block(2)
            flnt(i,j) = mpi_doubles_block(3)
            flntc(i,j) = mpi_doubles_block(4)
            fln200(i,j) = mpi_doubles_block(5)
            fln200c(i,j) = mpi_doubles_block(6)
            fsds(i,j) = mpi_doubles_block(7)
            fsdsc(i,j) = mpi_doubles_block(8)
            fsns(i,j) = mpi_doubles_block(9)
            fsnsc(i,j) = mpi_doubles_block(10)
            fsnt(i,j) = mpi_doubles_block(11)
            fsntc(i,j) = mpi_doubles_block(12)
            fsn200(i,j) = mpi_doubles_block(13)
            fsn200c(i,j) = mpi_doubles_block(14)
            fsnirtoa(i,j) = mpi_doubles_block(15)
            fsnirtoac(i,j) = mpi_doubles_block(16)
            solin(i,j) = mpi_doubles_block(17)
            soll(i,j) = mpi_doubles_block(18)
            solld(i,j) = mpi_doubles_block(19)
            sols(i,j) = mpi_doubles_block(20)
            solsd(i,j) = mpi_doubles_block(21)
            cldtau_liq_sw(i,j) = mpi_doubles_block(22)
            cldtau_ice_sw(i,j) = mpi_doubles_block(23)
            cldtau_lw(i,j) = mpi_doubles_block(24)
            oasdir(i,j) = mpi_doubles_block(25)
            oasdif(i,j) = mpi_doubles_block(26)
            oaldir(i,j) = mpi_doubles_block(27)
            oaldif(i,j) = mpi_doubles_block(28)
            flwds(i,j) = mpi_doubles_block(29)
           
            tag=6  
            call our_mpirecvdouble (mpi_level_block_double, 2*pver+2*pverp+aertau_species, proc_number, tag)
            fsn(i,:,j) = mpi_level_block_double(1:pverp)
            fln(i,:,j) = mpi_level_block_double(1+pverp:2*pverp)
            qrl(i,:,j) = mpi_level_block_double(1+2*pverp:2*pverp+pver)
            qrs(i,:,j) = mpi_level_block_double(1+2*pverp+pver:2*pverp+2*pver)
            aertau(i,:,j) = mpi_level_block_double(1+2*pverp+2*pver:aertau_species+2*pverp+2*pver)
            !END ADD OF RADCSWMX

            !BEGIN DRF FOR FLUXES
            tag=7
            call our_mpirecvreal (mpi_flux_block, 2*pver*6, proc_number, tag)

            flx_updiffuse_clr(i,:,1,j) = mpi_flux_block(1:pver)
            flx_updiffuse_clr(i,:,2,j) = mpi_flux_block(1+pver:2*pver) 
            flx_updiffuse_all(i,:,1,j) = mpi_flux_block(1+2*pver:3*pver) 
            flx_updiffuse_all(i,:,2,j) = mpi_flux_block(1+3*pver:4*pver) 
            flx_dndiffuse_clr(i,:,1,j) = mpi_flux_block(1+4*pver:5*pver) 
            flx_dndiffuse_clr(i,:,2,j) = mpi_flux_block(1+5*pver:6*pver) 
            flx_dndiffuse_all(i,:,1,j) = mpi_flux_block(1+6*pver:7*pver) 
            flx_dndiffuse_all(i,:,2,j) = mpi_flux_block(1+7*pver:8*pver) 
            flx_dndirect_clr(i,:,1,j) = mpi_flux_block(1+8*pver:9*pver) 
            flx_dndirect_clr(i,:,2,j) = mpi_flux_block(1+9*pver:10*pver) 
            flx_dndirect_all(i,:,1,j) = mpi_flux_block(1+10*pver:11*pver) 
            flx_dndirect_all(i,:,2,j) = mpi_flux_block(1+11*pver:12*pver) 
            !END DRF FOR FLUXES

          endif !rank 0
        enddo !i
      enddo !j

if (my_rank == 1) then
r1_message_end = our_mpi_wtime()
endif


    endif !both ranks
!    call our_mpibarrier ()
  enddo !end of loop of processor number 




if (my_rank == 0) then
r0_message_end = our_mpi_wtime()
endif


!CAA no messages
!end if


   
     !Note: this here for now since at the end of the do loop itime will be increased by 1 thus printing out
     !another slice of time data which is irrelevant to us

     if (our_master_proc) then
           !write(*,*) 'bb_upflux_clr 204_91 in Main before write = ',bb_updiffuse_clr(204,:,91)

     call write_output(opath, itime, FLN, FLNS, FLNSC, FLNT, FLNTC, &
                       FLN200, FLN200C,FLWDS, &
                       FSDS, FSDSC, FSN, FSNS, FSNSC, FSNT, FSNTC, &
                       FSN200, FSN200C, FSNIRTOA,FSNIRTOAC, &
                       QRL, QRS, SOLIN, SOLL, SOLLD, SOLS, SOLSD,c_lat,c_lon,  &
                       WAVELENGTH_LRES,RADIANCE_LRES_CLR,RADIANCE_LRES_ALL,WAVELENGTH_HRES, &
                       RADIANCE_HRES_CLR,RADIANCE_HRES_ALL,SOLAR_FLUX,DIFFUSE_FLUX_CLR,DIFFUSE_FLUX_ALL, & 
                       SOLAR_ZENITH,AOD,TAU_SULF,TAU_DUST,TAU_SOOT,TAU_SSLT,NCONFIG_REAL, &
                       CLDTAU_ICE_SW,CLDTAU_LIQ_SW,CLDTAU_LW, &
                       BB_UPDIFFUSE_CLR,BB_UPDIFFUSE_ALL, &
                       BB_DNDIFFUSE_CLR,BB_DNDIFFUSE_ALL, &
                       BB_DNDIRECT_CLR,BB_DNDIRECT_ALL, &
                       FLX_UPDIFFUSE_CLR,FLX_UPDIFFUSE_ALL, &
                       FLX_DNDIFFUSE_CLR,FLX_DNDIFFUSE_ALL, &
                       FLX_DNDIRECT_CLR,FLX_DNDIRECT_ALL,qr_option ) !DRF


     write(*,*) "finished with write_output"
 
     end if !our_master_proc
  
  end do ! end of loop over time


 

  call t_prf(0)




if (my_rank == 0) then
total_time_end = our_mpi_wtime()
endif


if (my_rank == 0) then
 write(60,*)"Total time of the program is ",total_time_end-total_time_start
 write(60,*)"Message time is ",r0_message_end-r0_message_start
endif


if (my_rank == 1) then
 write(61,*)"Message time is ",r1_message_end-r1_message_start
endif








if (my_rank == 0) then
CLOSE (UNIT=60) 
endif

if (my_rank == 1) then
CLOSE (UNIT=61) 
endif





#if ( defined MODTRAN_SPMD )
  call our_mpifinalize() 
#endif


  stop "End of program"
end program main
