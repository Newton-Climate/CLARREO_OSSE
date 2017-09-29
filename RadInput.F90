! Correct settings for operation
! Set build_aermmr = ??
!     build_emis  = .TRUE.
!     build_ozone = .TRUE. ozone data set downloaded from http://www.cesm.ucar.edu/models/atm-cam/download/
!     build_re    = .TRUE
!     build_trace = .TRUE. and ensure that the concentrations are up to date
!
!WDC -- 11/7/15 -- CMIP5 OSSE 

#include <params.h> !WDC -- 11/7/15 -- CMIP5 OSSE 

module RadInput

  use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4                                   
  use shr_const_mod, only: SHR_CONST_G 

  use, intrinsic :: IEEE_ARITHMETIC
!  use, intrinsic :: IEEE_EXCEPTIONS
!  use, intrinsic :: IEEE_FEATURES
  use ppgrid
  use pmgrid
  use error_messages
  use prescribed_aerosols, only: naer_all, &
                                 idxSUL, &
                                 idxSSLT, &
                                 idxDUSTfirst, &
                                 idxOCPHO, &
                                 idxOCPHI, &
                                 idxBCPHO, &
                                 idxBCPHI, &
                                 idxBG, &
                                 idxVOLC

  use ghg_surfvals, only:        co2vmr 

  implicit none
  save
!
! Fields input from history file
!
! 3D
!
  real(r8), allocatable :: inp_aermass(:,:,:) ! Aerosol masses per layer  
                                              !  Generic name: history vars =
                                              !    MBCPHI_V: BC hydrophilic
                                              !    MBCPHO_V: BC hydrophobic
                                              !    MBG_V   : background (SO4)
                                              !    MDUST1_V: dust size 1
                                              !    MDUST2_V: dust size 2
                                              !    MDUST3_V: dust size 3
                                              !    MDUST4_V: dust size 4
                                              !    MOCPHI_V: OC hydrophilic
                                              !    MOCPHO_V: OC hydrophobic
                                              !    MSSLT_V:  sea salt
                                              !    MSUL_V:   sulfate
                                              !    MVOLC:    volcanic
  real(r8), allocatable :: inp_allaer (:,:,:,:) ! All aerosol masses
  real(r8), allocatable :: inp_cfc11  (:,:,:) ! CFC11 mixing ratio 
  real(r8), allocatable :: inp_cfc12  (:,:,:) ! CFC12 mixing ratio 
  real(r8), allocatable :: inp_ch4    (:,:,:) ! CH4 mixing ratio 
  real(r8), allocatable :: inp_cloud  (:,:,:) ! cloud amount
  real(r8), allocatable :: inp_emis   (:,:,:) ! Cloud longwave emissivity 
  real(r8), allocatable :: inp_icldiwp(:,:,:) ! in-cloud ice water path
  real(r8), allocatable :: inp_icldlwp(:,:,:) ! in-cloud total water path
  real(r8), allocatable :: inp_fice(:,:,:)    ! cloud ice fraction
  real(r8), allocatable :: inp_tcldice(:,:,:) ! cloud ice total grid water box 
  real(r8), allocatable :: inp_tcldliq(:,:,:) ! cloud liq total grid water box
  real(r8), allocatable :: inp_n2o    (:,:,:) ! N2O mixing ratio 
  real(r8), allocatable :: inp_o3vmr  (:,:,:) ! Ozone volume mixing ratio 
  real(r8), allocatable :: inp_rel    (:,:,:) ! Liq. droplet eff. radius 
  real(r8), allocatable :: inp_rei    (:,:,:) ! Ice effective drop size  
  real(r8), allocatable :: inp_t_cld  (:,:,:) ! Atmospheric temperature for 
                                              !   cloud properties
  real(r4), allocatable :: stdlvl_q      (:,:,:)! Water vapor mixing ratio on standard-level grid
  real(r4), allocatable :: stdlvl_rh     (:,:,:)! Water vapor relative humidity on standard-level grid
  real(r4), allocatable :: stdlvl_rh_fix (:,:,:)! Water vapor relative humidity on standard-level grid
  real(r4), allocatable :: stdlvl_t      (:,:,:)! Atmospheric temperature on standard-level grid

  real(r4) :: fillvalue_q      ! Fillvalue Water vapor mixing ratio on standard-level grid
  real(r4) :: fillvalue_rh     ! Fillvalue Water vapor relative humidity on standard-level grid
  real(r4) :: fillvalue_t      ! Fillvalue Atmospheric temperature on standard-level grid

  real(r8), allocatable :: valid_values(:,:)  ! Valid values of rh, q, t

  integer :: level
  integer :: start_level
  integer :: end_level
  integer :: d_level
  integer :: integrand_first_level
  integer :: integrand_last_level
  integer :: integrand_TOA_level
  integer :: pint_first_level
  integer :: pint_last_level
  integer :: pint_TOA_level
  

  real(r8), allocatable :: inp_q       (:,:,:)! Water vapor mixing ratio on parent model grid
  real(r8), allocatable :: inp_rh      (:,:,:)! Water vapor relative humidity on parent model grid
  real(r8), allocatable :: inp_t       (:,:,:)! Atmospheric temperature on parent model grid
!
! 2D
!
  real(r8), allocatable :: inp_asdir   (:,:)  ! Albedo: shortwave, direct 
  real(r8), allocatable :: inp_asdif   (:,:)  ! Albedo: shortwave, diffuse 
  real(r8), allocatable :: inp_aldir   (:,:)  ! Albedo: longwave, direct 
  real(r8), allocatable :: inp_aldif   (:,:)  ! Albedo: longwave, diffuse 
  real(r8), allocatable :: inp_icefrac (:,:)  ! Fractional sea-ice amount 
  real(r8), allocatable :: inp_landfrac(:,:)  ! Fractional land amount 
  real(r8), allocatable :: inp_landmcos(:,:)  ! Landm coslat field used for rel/rei 
  real(r8), allocatable :: inp_ps      (:,:)  ! Surface pressure 
  real(r8), allocatable :: inp_snowh   (:,:)  ! Snow height (equivalent water depth) 
  real(r8), allocatable :: inp_lwup    (:,:)  ! Upwelling longwave flux 
  real(r8), allocatable :: inp_ts      (:,:)  ! Radiative surface temperature (DRF) 
  real(r8), allocatable :: inp_windspeed  (:,:)  ! surface wind speed (DRF) 
!
! 1D
!
  real(r8), allocatable :: inp_lat     (:)    ! latitudes 
  real(r8), allocatable :: inp_lon     (:)    ! longitudes 

  integer, parameter :: stdlvl_pver  = SLEV        ! WDC -- 11/7/15 -- CMIP5 OSSE number of standard vertical levels
  integer, parameter :: stdlvl_pverp = (stdlvl_pver + 1 ) ! WDC -- 11/7/15 -- CMIP5 OSSE number of standard vertical levels + 1

!  WDC -- 11/7/15 -- CMIP5 OSSE
! Standard level coordinates  WDC -- 11/7/15 -- CMIP5 OSSE
!  WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8) :: stdlvl_hyam(stdlvl_pver) ! vertical coordinate formula term: a(k) WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8) :: stdlvl_hyai(stdlvl_pverp) ! vertical coordinate formula term: a(k) WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8) :: stdlvl_hybm(stdlvl_pver) ! vertical coordinate formula term: b(k) WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8) :: stdlvl_hybi(stdlvl_pverp) ! vertical coordinate formula term: b(k) WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8) :: stdlvl_ps0 ! vertical coordinate formula term: reference pressure WDC -- 11/7/15 -- CMIP5 OSSE 
  real(r8), allocatable :: stdlvl_ps(:,:)  ! surface_air_pressure WDC -- 11/7/15 -- CMIP5 OSSE

  real(r8), allocatable :: stdlvl_pint(:,:,:) ! standard-level layer interface pressures  WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: stdlvl_pint_plus_toa(:,:,:) ! standard-level layer interface pressures  WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: stdlvl_pdel(:,:,:) ! standard-level pressure increment         WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable, private :: pint(:,:,:)  ! interface pressure (Pa)          WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable, private :: pdel(:,:,:)  ! pressure increment (Pa)          WDC -- 11/7/15 -- CMIP5 OSSE

!  real(r8), allocatable :: P2S(:,:,:,:) ! Matrix for linear interpolation from standard to cloud coordinates
  real(r8), allocatable :: S2P(:,:,:,:) ! Matrix for linear interpolation from cloud to standard coordinates

  real(r8), allocatable :: stdlvl_sat_spec_hum(:,:,:) ! saturated spec humidity      WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: stdlvl_sat_vapor_mass(:,:,:) ! saturated precipitable water      WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: stdlvl_int_sat_vapor_mass(:,:,:) ! cumulative saturated precipitable water      WDC -- 11/7/15 -- CMIP5 OSSE

  real(r8), allocatable :: sat_spec_hum(:,:,:) ! saturated spec humidity, on standard-level grid       WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: sat_vapor_mass(:,:,:) ! sat. precip water, on standard-level grid      WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: int_sat_vapor_mass(:,:,:) ! cum. sat. precip water, on standard-level grid      WDC -- 11/7/15 -- CMIP5 OSSE

  real(r8), allocatable :: stdlvl_spec_hum(:,:,:) ! spec humidity      WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: stdlvl_vapor_mass(:,:,:) ! precipitable water      WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: stdlvl_int_vapor_mass(:,:,:) ! cumulative precipitable water      WDC -- 11/7/15 -- CMIP5 OSSE

  real(r8), allocatable :: spec_hum(:,:,:) ! spec humidity, on standard-level grid       WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: vapor_mass(:,:,:) ! precip water, on standard-level grid      WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: int_vapor_mass(:,:,:) ! cum. precip water, on standard-level grid      WDC -- 11/7/15 -- CMIP5 OSSE

  real(r8), allocatable :: stdlvl_temp(:,:,:) ! temperature      WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: stdlvl_internal_energy(:,:,:) ! internal energy      WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: stdlvl_int_internal_energy(:,:,:) ! cumulative internal energy      WDC -- 11/7/15 -- CMIP5 OSSE

  real(r8), allocatable :: temp(:,:,:) ! temperature, on standard-level grid       WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: internal_energy(:,:,:) ! internal energy, on standard-level grid      WDC -- 11/7/15 -- CMIP5 OSSE
  real(r8), allocatable :: int_internal_energy(:,:,:) ! cum. internal energy, on standard-level grid      WDC -- 11/7/15 -- CMIP5 OSSE

  integer, allocatable, target     ::  date(:)     ! current date
  double precision, allocatable, target     ::  day_of_yr(:)   ! current calendar day of yr

!  real, allocatable :: inp_co2vmr      (:)    ! co2vmr DRF

! Indices for aerosol mmmr fields expected by radcswmx
  integer, public, parameter :: aerosol_index(naer_all) = &
      (/ idxSUL, &
         idxSSLT, &
         idxDUSTfirst, &
         idxDUSTfirst+1, &
         idxDUSTfirst+2, &
         idxDUSTfirst+3, &
         idxOCPHO, &
         idxBCPHO, &
         idxOCPHI, &
         idxBCPHI, &
         idxBG, &
         idxVOLC /)

!
! Vertical cloud coordinate has to be flipped
!
  logical :: flip_z

  integer :: z0
  integer :: z1
  integer :: z0p
  integer :: z1p
  integer :: dz

  real(r8), allocatable :: inp_temp  (:,:,:) ! buffer for flipping z

CONTAINS 

  subroutine input_data(path1, path2, fld_option, itime, &
             build_aermmr, build_trace, build_emis, build_re, build_ozone)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Input instantaneous history fields output each radiation time step
!    required to regenerate instantaneous radiation fluxes
!
! Note: Timing data required to reconstitute day_of_yr is read using
!    input_times
! 
! Method: 
! Use NetCDF wrapper routines to input a single time slice.
! Data is tranposed onto internal arrays for parallelization in Chunks
! module.
!
! Author: W. Collins
! 
!-----------------------------------------------------------------------

    implicit none
#include <netcdf.inc>
#include <comhyb.h>
!
! Input arguments
!
    character(len=*), intent(in) :: path1        ! path to data
    character(len=*), intent(in) :: path2        ! path to data
    integer,      intent(in) :: fld_option   ! Option for field switching
                                             !   0 = no switch
                                             !   1 = swap temperatures
                                             !   2 = swap spec. humidity (Q)
                                             !   3 = swap clouds
                                             !   4 = swap albedos
                                             !   5 = perform swaps 1-4

    integer, intent(in) :: itime             ! time slice
    character(len=16), intent(in) :: build_aermmr ! Build AERMMR internally
    logical, intent(in) :: build_trace       ! Build CFCs,CH4,& N2O internally
    logical, intent(in) :: build_emis        ! Build EMIS internally
    logical, intent(in) :: build_re          ! Build RE internally
    logical, intent(in) :: build_ozone       ! Build O3 internally 
!
! Local variables
!
    integer, parameter :: SWAP_TEMP  = 1     ! Flag for temp. swap
    integer, parameter :: SWAP_HUMID = 2     ! Flag for humid. swap
    integer, parameter :: SWAP_CLDS  = 3     ! Flag for temp. swap
    integer, parameter :: SWAP_ALBS  = 4     ! Flag for temp. swap
    integer, parameter :: SWAP_ALL   = 5     ! Flag for swaps 1-4

    integer :: nfid1                         ! NetCDF id for path1
    integer :: nfid2                         ! NetCDF id for path2
    integer :: NFIDX                         ! nfid1 or nfid2
    integer :: start(4)                      ! start indices
    integer :: kount(4)                      ! kount indices
    integer :: istat                         ! allocate status
    integer :: varid                         ! variable id
    character(len = 10) :: routine = "input_data"

    integer whacked                          ! Number of bad cloud points
    integer iaer                             ! Aerosol index
    integer ii,jj,kk                         ! DRF indices for recreating IWP
    real(r8) dummy_var                       ! DRF for assigning IWP, LWP values
! Names of aerosol mass fields on history file
    character(len=8), parameter :: aerosol_name(naer_all) =  &
     (/"MSUL_V  "&
      ,"MSSLT_V "&
      ,"MDUST1_V"&
      ,"MDUST2_V"&
      ,"MDUST3_V"&
      ,"MDUST4_V"&
      ,"MOCPHO_V"&
      ,"MBCPHO_V"&
      ,"MOCPHI_V"&
      ,"MBCPHI_V"&
      ,"MBG_V   "&
      ,"MVOLC   "/)

    real(r8), parameter :: ieee_var_type = 1.0

!----------------------------------------------------------------------

    call wrap_open(path1, nf_nowrite, nfid1)
    call wrap_open(path2, nf_nowrite, nfid2)

!
! Setup temp input buffer for flipping vertical coordinate
!
    allocate(inp_temp(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_temp",(plon*plat*pver) )    
!
!WDC 11/23/15 -- null out everything
!    
! 3D
    if (build_aermmr == 'NONE' .or. build_aermmr == 'IPCC') then
       allocate(inp_aermass(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_aermass",(plon*plat*pver) )
    endif
       
    if (build_aermmr == 'NONE') then
       allocate(inp_allaer(plon,plat,pver,naer_all), stat = istat) 
       call alloc_err(istat, routine, "inp_allaer",(plon*plat*pver*naer_all) )
       
       inp_aermass = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
       do iaer = 1, naer_all
          inp_allaer(:,:,:,aerosol_index(iaer)) = inp_aermass
       end do
    endif

    if (build_aermmr == 'IPCC') then
       inp_aermass = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
    endif
    
    if (.not. build_trace) then 
       allocate(inp_cfc11(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_cfc11",(plon*plat*pver) )
       inp_cfc11 = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
       
       allocate(inp_cfc12(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_cfc12",(plon*plat*pver) )
       inp_cfc12 = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
    
       allocate(inp_ch4(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_ch4",(plon*plat*pver) )
       inp_ch4 = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
       write(*,*) "Need to fix this for our CH4 calculation"

       allocate(inp_n2o(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_n2o",(plon*plat*pver) )
       inp_n2o = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
    endif

    if (.not. build_emis) then
       allocate(inp_emis(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_emis",(plon*plat*pver) )
       inp_emis = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
    endif
    
    if (.not. build_ozone) then
       allocate(inp_o3vmr(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_o3vmr",(plon*plat*pver) )
       inp_o3vmr = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
    endif
        
    if (.not. build_re) then 
       allocate(inp_rel(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_rel",(plon*plat*pver) )
       inp_rel = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
       
       allocate(inp_rei(plon,plat,pver), stat = istat) 
       call alloc_err(istat, routine, "inp_rei",(plon*plat*pver) )
       inp_rei = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
    endif

!2D

    allocate(inp_asdir(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_asdir",(plon*plat) )
    inp_asdir = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
    
    allocate(inp_asdif(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_asdif",(plon*plat) )
    inp_asdif = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
    
    allocate(inp_aldir(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_aldir",(plon*plat) )
    inp_aldir = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE
    
    allocate(inp_aldif(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_aldif",(plon*plat) )
    inp_aldif = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN) !WDC -- 11/7/15 -- CMIP5 OSSE

!0D

    !co2vmr = IEEE_VALUE(ieee_var_type, IEEE_QUIET_NAN)
    co2vmr = 0.0003949625 !DRF 03/08/16

!
!WDC 11/23/15 -- end of zeroing out
!   
    allocate(inp_cloud(plon,plat,pver), stat = istat) !WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "inp_cloud",(plon*plat*pver) ) !WDC -- 11/7/15 -- CMIP5 OSSE
    start = (/1,   1,   1,   itime/) 
    kount = (/plon,plat,pver,1/) 
    if (fld_option == SWAP_CLDS .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    call wrap_inq_varid(NFIDX, "CLOUD", varid)
    call wrap_get_vara_realx(NFIDX, varid, start, kount, inp_temp) !WDC -- 11/7/15 -- CMIP5 OSSE
    inp_cloud = inp_temp(:,:,z0:z1:dz)
    
    inp_cloud = inp_cloud / 100.0 !DRF 08/03/16
!
! Read in Q, RH, and T on original grid, ther interpolate
!    
    allocate(stdlvl_q(plon,plat,stdlvl_pver), stat = istat) 
    call alloc_err(istat, routine, "stdlvl_q",(plon*plat*stdlvl_pver) )
    start = (/1,   1,   1,   itime/) 
    kount = (/plon,plat,stdlvl_pver,1/) 
    if (fld_option == SWAP_HUMID .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    call wrap_inq_varid(NFIDX, "Q", varid)
    call wrap_get_vara_real4(NFIDX, varid, start, kount, stdlvl_q)
    call wrap_get_att_real(NFIDX, varid, "_FillValue", fillvalue_q)

    allocate(stdlvl_rh(plon,plat,stdlvl_pver), stat = istat) 
    call alloc_err(istat, routine, "stdlvl_rh",(plon*plat*stdlvl_pver) )
    allocate(stdlvl_rh_fix(plon,plat,stdlvl_pver), stat = istat) 
    call alloc_err(istat, routine, "stdlvl_rh_fix",(plon*plat*stdlvl_pver) )

    start = (/1,   1,   1,   itime/) 
    kount = (/plon,plat,stdlvl_pver,1/) 
    if (fld_option == SWAP_HUMID .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    call wrap_inq_varid(NFIDX, "RH", varid) !WDC -- 11/7/15 -- CMIP5 OSSE
    call wrap_get_vara_real4(NFIDX, varid, start, kount, stdlvl_rh)
    call wrap_get_att_real(NFIDX, varid, "_FillValue", fillvalue_rh)

    allocate(stdlvl_t(plon,plat,stdlvl_pver), stat = istat) 
    call alloc_err(istat, routine, "stdlvl_t",(plon*plat*stdlvl_pver) )
    start = (/1,   1,   1,   itime/) 
    kount = (/plon,plat,stdlvl_pver,1/) 
    if (fld_option == SWAP_TEMP .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    call wrap_inq_varid(NFIDX, "T", varid)
    call wrap_get_vara_real4(NFIDX, varid, start, kount, stdlvl_t)
    call wrap_get_att_real(NFIDX, varid, "_FillValue", fillvalue_t)
    
    allocate(inp_q(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_q",(plon*plat*pver) )

    allocate(inp_rh(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_rh",(plon*plat*pver) )

    allocate(inp_t(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_t",(plon*plat*pver) )
    
!
! End of reading in Q, T, and RH
!
    allocate(inp_t_cld(plon,plat,pver), stat = istat) 
    call alloc_err(istat, routine, "inp_t_cld", (plon*plat*pver) )

    allocate(inp_tcldice(plon,plat,pver), stat = istat) !WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "inp_tcldice",(plon*plat*pver) ) !WDC -- 11/7/15 -- CMIP5 OSSE
    start = (/1,   1,   1,   itime/) !WDC -- 11/7/15 -- CMIP5 OSSE
    kount = (/plon,plat,pver,1/)     !WDC -- 11/7/15 -- CMIP5 OSSE
    call wrap_inq_varid(NFIDX, "CLDIWP", varid) !WDC -- 11/7/15 -- CMIP5 OSSE
    call wrap_get_vara_realx(NFIDX, varid, start, kount, inp_temp) !WDC -- 11/7/15 -- CMIP5 OSSE
    inp_tcldice = inp_temp(:,:,z0:z1:dz)
    
    allocate(inp_tcldliq(plon,plat,pver), stat = istat) !WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "inp_tcldliq",(plon*plat*pver) ) !WDC -- 11/7/15 -- CMIP5 OSSE
    start = (/1,   1,   1,   itime/) !WDC -- 11/7/15 -- CMIP5 OSSE
    kount = (/plon,plat,pver,1/)     !WDC -- 11/7/15 -- CMIP5 OSSE
    call wrap_inq_varid(NFIDX, "CLDLWP", varid) !WDC -- 11/7/15 -- CMIP5 OSSE
    call wrap_get_vara_realx(NFIDX, varid, start, kount, inp_temp) !WDC -- 11/7/15 -- CMIP5 OSSE
    inp_tcldliq = inp_temp(:,:,z0:z1:dz)

! 2D 

    allocate(stdlvl_ps(plon,plat), stat = istat)                    !WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_ps",(plon*plat) )        !WDC -- 11/7/15 -- CMIP5 OSSE
    start = (/1,   1,   itime, -1/)                           !WDC -- 11/7/15 -- CMIP5 OSSE
    kount = (/plon,plat,1,     -1/)                           !WDC -- 11/7/15 -- CMIP5 OSSE
    call wrap_inq_varid(nfid1, "PS", varid)                  !WDC -- 11/7/15 -- CMIP5 OSSE
    call wrap_get_vara_realx(nfid1, varid, start, kount, stdlvl_ps) !WDC -- 11/7/15 -- CMIP5 OSSE

    if (build_re) then
       allocate(inp_icefrac(plon,plat), stat = istat) 
       call alloc_err(istat, routine, "inp_icefrac",(plon*plat) )
       start = (/1,   1,   itime, -1/) 
       kount = (/plon,plat,1,     -1/) 
       call wrap_inq_varid(nfid1, "sic", varid) !WDC -- 11/7/15 -- CMIP5 OSSE
       call wrap_get_vara_realx(nfid1, varid, start, kount, inp_icefrac)
!
! Correct sea ice fraction to 0 to 1 range, if in percent
!
       inp_icefrac = inp_icefrac / 100.0 !DRF (07/26/16)
       !inp_icefrac = inp_icefrac * 10**(-nint(maxval(log10(inp_icefrac))))

       allocate(inp_landfrac(plon,plat), stat = istat) 
       call alloc_err(istat, routine, "inp_landfrac",(plon*plat) )
       start = (/1,   1,   itime, -1/) 
       kount = (/plon,plat,1,     -1/) 
       call wrap_inq_varid(nfid1, "LANDFRAC", varid)
       call wrap_get_vara_realx(nfid1, varid, start, kount, inp_landfrac)
!
! Correct landfrac to 0 to 1 range, if in percent
!
       inp_landfrac = inp_landfrac / 100.0 !DRF (07/26/16)
       !DRF (07/26/16) inp_landfrac = inp_landfrac * 10**(-nint(maxval(log10(inp_landfrac))))

       !write(*,*) 'inp_landfrac ',maxval(inp_landfrac)
 
       allocate(inp_snowh(plon,plat), stat = istat) 
       call alloc_err(istat, routine, "inp_snowh",(plon*plat) )
       start = (/1,   1,   itime, -1/) 
       kount = (/plon,plat,1,     -1/) 
       call wrap_inq_varid(nfid1, "snd", varid) !WDC -- 11/7/15 -- CMIP5 OSSE
       call wrap_get_vara_realx(nfid1, varid, start, kount, inp_snowh)
    endif

    allocate(inp_lwup(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_lwup",(plon*plat) )
    start = (/1,   1,   itime, -1/)
    kount = (/plon,plat,1    , -1/)
    if (fld_option == SWAP_TEMP .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    call wrap_inq_varid(NFIDX, "LWUP_r", varid) !WDC -- 11/7/15 -- CMIP5 OSSE
    call wrap_get_vara_realx(NFIDX, varid, start, kount, inp_lwup) !WDC -- 11/7/15 -- CMIP5 OSSE

    allocate(inp_ps(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_ps",(plon*plat) )
    start = (/1,   1,   itime, -1/) 
    kount = (/plon,plat,1,     -1/)
!
! WDC -- We are now using cloud coordinates as the primary vertical coordinates
! 
    call wrap_inq_varid(nfid1, "cps", varid)
    call wrap_get_vara_realx(nfid1, varid, start, kount, inp_ps)

    allocate(inp_ts(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_ts",(plon*plat) )
    start = (/1,   1,   itime, -1/)
    kount = (/plon,plat,1    , -1/)
    if (fld_option == SWAP_TEMP .or. fld_option == SWAP_ALL) then
       NFIDX = nfid2
    else
       NFIDX = nfid1
    endif
    call wrap_inq_varid(NFIDX, "ts", varid)
    call wrap_get_vara_realx(NFIDX, varid, start, kount, inp_ts)

    allocate(inp_windspeed(plon,plat), stat = istat)
    call alloc_err(istat, routine, "inp_windspeed",(plon*plat) )
    start = (/1,   1,   itime, -1/) !WDC -- 11/7/15 -- CMIP5 OSSE
    kount = (/plon,plat,1,     -1/) !WDC -- 11/7/15 -- CMIP5 OSSE
    call wrap_inq_varid(nfid1, "sfcWind", varid) !WDC -- 11/7/15 -- CMIP5 OSSE
    call wrap_get_vara_realx(nfid1, varid, start, kount, inp_windspeed) !WDC -- 11/7/15 -- CMIP5 OSSE

! !D
    allocate(inp_lat(plat), stat = istat) 
    call alloc_err(istat, routine, "inp_lat",(plat) )
    call wrap_inq_varid(nfid1, "lat", varid)
    call wrap_get_var_realx(nfid1, varid, inp_lat)

    allocate(inp_lon(plon), stat = istat) 
    call alloc_err(istat, routine, "inp_lon",(plon) )
    call wrap_inq_varid(nfid1, "lon", varid)
    call wrap_get_var_realx(nfid1, varid, inp_lon)

!WDC -- 11/7/15 -- CMIP5 OSSE

!
! First, allocate all the intermediate result storage required
!
    allocate(stdlvl_pdel(plon,plat,stdlvl_pver), stat = istat) ! cloud mid-point pressures        WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_pdel", (plon*stdlvl_pver*plat)) !WDC -- 11/7/15 -- CMIP5 OSSE

    allocate(stdlvl_pint(plon,plat,stdlvl_pverp), stat = istat) ! cloud layer interface pressures  WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_pint", (plon*stdlvl_pverp*plat)) !WDC -- 11/7/15 -- CMIP5 OSSE

    allocate(stdlvl_pint_plus_toa(plon,plat,(stdlvl_pverp+1)), stat = istat) ! cloud layer interface pressures  WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_pint_plus_toa", (plon*(stdlvl_pverp+1)*plat)) !WDC -- 11/7/15 -- CMIP5 OSSE

    allocate(pdel(plon,plat,pver), stat = istat)  ! midpoint pressure (Pa)           WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "pdel", (plon*pver*plat)) !WDC -- 11/7/15 -- CMIP5 OSSE

    allocate(pint(plon,plat,pverp), stat = istat)  ! interface pressure (Pa)          WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "pint", (plon*pverp*plat)) !WDC -- 11/7/15 -- CMIP5 OSSE

    allocate(stdlvl_sat_spec_hum(plon,plat,stdlvl_pver)) ! saturated spec humidity, on cloud grid       WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_sat_spec_hum", (plon*plat*stdlvl_pver))
    allocate(stdlvl_sat_vapor_mass(plon,plat,stdlvl_pverp)) ! sat. precip water, on cloud grid      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_sat_vapor_mass", (plon*plat*stdlvl_pverp)) ! sat. precip water, on cloud grid      WDC -- 11/7/15 -- CMIP5 OSSE
    allocate(stdlvl_int_sat_vapor_mass(plon,plat,stdlvl_pverp+1)) ! cum. sat. precip water, on cloud grid      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_int_sat_vapor_mass", (plon*plat*(stdlvl_pverp+1)))

    allocate(sat_spec_hum(plon,plat,pver)) ! saturated spec humidity      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "sat_spec_hum", (plon*plat*pver))
    allocate(sat_vapor_mass(plon,plat,pver)) ! saturated precipitable water      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "sat_vapor_mass", (plon*plat*pver))
    allocate(int_sat_vapor_mass(plon,plat,pverp)) ! cumulative saturated precipitable water      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "xint_sat_vapor_mass", (plon*plat*pverp))

    allocate(stdlvl_spec_hum(plon,plat,stdlvl_pver)) ! spec humidity, on cloud grid       WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_spec_hum", (plon*plat*stdlvl_pver))
    allocate(stdlvl_vapor_mass(plon,plat,stdlvl_pverp)) ! precip water, on cloud grid      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_vapor_mass", (plon*plat*stdlvl_pverp)) ! sat. precip water, on cloud grid      WDC -- 11/7/15 -- CMIP5 OSSE
    allocate(stdlvl_int_vapor_mass(plon,plat,stdlvl_pverp+1)) ! cum. precip water, on cloud grid      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_int_vapor_mass", (plon*plat*(stdlvl_pverp+1)))

    allocate(spec_hum(plon,plat,pver)) ! spec humidity      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "spec_hum", (plon*plat*pver))
    allocate(vapor_mass(plon,plat,pver)) ! precipitable water      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "vapor_mass", (plon*plat*pver))
    allocate(int_vapor_mass(plon,plat,pverp)) ! cumulative precipitable water      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "xint_vapor_mass", (plon*plat*pverp))

    allocate(stdlvl_temp(plon,plat,stdlvl_pver)) ! temperature, on cloud grid       WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_temp", (plon*plat*stdlvl_pver))
    allocate(stdlvl_internal_energy(plon,plat,stdlvl_pverp)) ! internal energy, on cloud grid      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_internal_energy", (plon*plat*stdlvl_pverp)) !  internal energy, on cloud grid      WDC -- 11/7/15 -- CMIP5 OSSE
    allocate(stdlvl_int_internal_energy(plon,plat,stdlvl_pverp+1)) ! cum. internal energy, on cloud grid      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "stdlvl_int_internal_energy", (plon*plat*(stdlvl_pverp+1)))

    allocate(temp(plon,plat,pver)) ! temperature      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "temp", (plon*plat*pver))
    allocate(internal_energy(plon,plat,pver)) ! internal energy      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "internal_energy", (plon*plat*pver))
    allocate(int_internal_energy(plon,plat,pverp)) ! cumulative internal energy      WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "int_internal_energy", (plon*plat*pverp))

    allocate(inp_icldlwp(plon,plat,pver), stat = istat) !WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "inp_icldlwp",(plon*plat*pver) ) !WDC -- 11/7/15 -- CMIP5 OSSE
    allocate(inp_icldiwp(plon,plat,pver), stat = istat) !WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "inp_icldiwp",(plon*plat*pver) )!WDC -- 11/7/15 -- CMIP5 OSSE
    allocate(inp_fice(plon,plat,pver), stat = istat) !WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "inp_fice",(plon*plat*pver) )!WDC -- 11/7/15 -- CMIP5 OSSE

!    allocate(P2S(plon,plat,pverp,stdlvl_pverp), stat = istat) !WDC -- 11/7/15 -- CMIP5 OSSE
!    call alloc_err(istat, routine, "P2S", (plon*plat*pverp*stdlvl_pverp))
    allocate(S2P(plon,plat,stdlvl_pverp+1,pverp), stat = istat) !WDC -- 11/7/15 -- CMIP5 OSSE
    call alloc_err(istat, routine, "S2P", (plon*plat*(stdlvl_pverp+1) *pverp))

!
! 
!
    allocate(valid_values(plon,plat), stat = istat)
    call alloc_err(istat, routine, "valid_values", (plon * plat))

! Construct the required pressure fields
! WDC CMIP5 OSSE -- units of stdlvl_ps0 and ps0 may not be consistent
    stdlvl_ps0 = stdlvl_ps0 * 10.0**(nint(log10(ps0)) - nint(log10(stdlvl_ps0)))
    stdlvl_ps = stdlvl_ps * 10.0**(nint(log10(maxval(inp_ps))) - nint(log10(maxval(stdlvl_ps))))

    stdlvl_pint(:,:,:) = &
         stdlvl_ps0 * spread(spread(stdlvl_hyai, 1, plat), 1, plon) + &                 !WDC -- 11/7/15 -- CMIP5 OSSE
         spread(stdlvl_ps, 3, stdlvl_pverp) * spread(spread(stdlvl_hybi, 1, plat), 1, plon)   !WDC -- 11/7/15 -- CMIP5 OSSE

    stdlvl_pdel(:,:,:) = stdlvl_pint(:,:,2:stdlvl_pverp) - stdlvl_pint(:,:,1:stdlvl_pver) !WDC -- 11/7/15 -- CMIP5 OSSE

    pint(:,:,:)  = &
         ps0 * spread(spread(hyai, 1, plat), 1, plon) + &                 !WDC -- 11/7/15 -- CMIP5 OSSE 
         spread(inp_ps, 3, pverp) * spread(spread(hybi, 1, plat), 1, plon) !WDC -- 11/7/15 -- CMIP5 OSSE 

    pint(:,:,1) = 100.0 !DRF 08/07/16 (no zeros in interface pressure)

    pdel(:,:,:) = pint(:,:,2:pverp) - pint(:,:,1:pver) !WDC -- 11/7/15 -- CMIP5 OSSE

!
! Determine sign of pressure increase
!
    if (all(stdlvl_pint(:,:,1) .le. stdlvl_pint(:,:,stdlvl_pverp))) then
       start_level = 1
       end_level = stdlvl_pver
       d_level = 1
       integrand_first_level = 2
       integrand_last_level = stdlvl_pverp
       integrand_TOA_level = 1
       pint_first_level = 2
       pint_last_level = stdlvl_pverp+1
       pint_TOA_level = 1
    else
       start_level = stdlvl_pver
       end_level = 1
       d_level = -1
       integrand_first_level = 1
       integrand_last_level = stdlvl_pver
       integrand_TOA_level = stdlvl_pverp
       pint_first_level = 1
       pint_last_level = stdlvl_pverp
       pint_TOA_level = stdlvl_pverp+1
    endif

!
! Construct the interpolation matrices
!
!    call build_interpolation_matrix(P2S, pint, pverp, stdlvl_pint, stdlvl_pverp, .false.)

    stdlvl_pint_plus_toa(:,:,pint_first_level:pint_last_level) = stdlvl_pint
    stdlvl_pint_plus_toa(:,:,pint_TOA_level) = 0.0
    call build_interpolation_matrix(S2P, stdlvl_pint_plus_toa, stdlvl_pverp+1, pint, pverp, .false.)

! 
! Interpolate specific humidity
!
! Fix q
!
    stdlvl_spec_hum = stdlvl_q
    valid_values = stdlvl_q(:,:,start_level)
    if (any(valid_values < 0.0 .or. valid_values > 1.0 .or. &
            valid_values == fillvalue_q)) then
       write(*,*) "q panic attack", minval(valid_values),maxval(valid_values)
       stop
    endif

    do level = start_level + d_level, end_level, d_level
       where (stdlvl_spec_hum(:,:,level) .eq. fillvalue_q)
          stdlvl_spec_hum(:,:,level) = valid_values
       elsewhere
          valid_values = stdlvl_spec_hum(:,:,level)
       endwhere
    end do

    stdlvl_vapor_mass(:,:,integrand_first_level:integrand_last_level) = stdlvl_spec_hum * stdlvl_pdel
    stdlvl_vapor_mass(:,:,integrand_TOA_level) = stdlvl_spec_hum(:,:,start_level) * stdlvl_pint(:,:,start_level)
    
    call cumsum(stdlvl_vapor_mass, stdlvl_int_vapor_mass, stdlvl_pverp)

    call interpolate(stdlvl_int_vapor_mass, int_vapor_mass, S2P, stdlvl_pverp+1, pverp)
        
    call difference( int_vapor_mass,  vapor_mass, pverp)

    spec_hum = vapor_mass / pdel
! 
! Interpolate saturated specific humidity
!
    stdlvl_rh_fix = stdlvl_rh 
    valid_values = stdlvl_rh(:,:,start_level)
    if (any(valid_values < 0.0 .or. valid_values > 250.0 .or. &
            valid_values == fillvalue_rh)) then
       write(*,*) "rh panic attack", minval(valid_values),maxval(valid_values)
       !stop
    endif

    do level = start_level + d_level, end_level, d_level
       where (stdlvl_rh_fix(:,:,level) .eq. fillvalue_rh)
          stdlvl_rh_fix(:,:,level) = valid_values
       elsewhere
          valid_values = stdlvl_rh_fix(:,:,level)
       endwhere
       where (stdlvl_rh_fix(:,:,level) < 0.0)
	  stdlvl_rh_fix(:,:,level) = 0.0
       endwhere
       where (stdlvl_rh_fix(:,:,level) > 100.0)
	  stdlvl_rh_fix(:,:,level) = 100.0
       endwhere
    end do

    stdlvl_sat_spec_hum = stdlvl_spec_hum / stdlvl_rh_fix * 100.0

    stdlvl_sat_vapor_mass(:,:,integrand_first_level:integrand_last_level) = stdlvl_sat_spec_hum * stdlvl_pdel
    stdlvl_sat_vapor_mass(:,:,integrand_TOA_level) = stdlvl_sat_spec_hum(:,:,start_level) * stdlvl_pint(:,:,start_level)
    
    call cumsum(stdlvl_sat_vapor_mass, stdlvl_int_sat_vapor_mass, stdlvl_pverp)

    call interpolate(stdlvl_int_sat_vapor_mass, int_sat_vapor_mass, S2P, stdlvl_pverp+1, pverp)
        
    call difference( int_sat_vapor_mass,  sat_vapor_mass, pverp)

    sat_spec_hum = sat_vapor_mass / pdel

!
! Construct inp_q, inp_rh
!
    !DRF fixing erroneous q values 08/22/16
    do ii = 1,plon 
	do jj = 1,plat
		do kk = 1,pver
			if (spec_hum(ii,jj,kk) < 0.0) then
				spec_hum(ii,jj,kk) = 0.0
!				write(*,*) 'changing q at ',ii,jj,kk
			end if
		end do
	end do
    end do
    !DRF end fix of erroneous q values

    inp_q = spec_hum
    inp_rh = 100.0 * spec_hum / sat_spec_hum

    do ii = 1,plon 
	do jj = 1,plat
		do kk = 1,pver
			if (inp_rh(ii,jj,kk) < 0.0) then
				inp_rh(ii,jj,kk) = 0.0
				write(*,*) 'changing rh at ',ii,jj,kk,spec_hum(ii,jj,kk),sat_spec_hum(ii,jj,kk)
			end if
		end do
	end do
    end do

    !whacked = count(inp_rh < 0  .or. inp_rh > 170.0 .or. inp_q < 0.0 .or. inp_q > 1.0)
    whacked = count(inp_rh < 0  .or. spec_hum < 0.0)
    if (whacked > 0) then
       write(*,*)  whacked, " instances of bogus q/rh out of  ", plon * plat * pver
       write(*,*) 'rh',minval(inp_rh),maxval(inp_rh), minval(stdlvl_rh_fix), maxval(stdlvl_rh_fix)
       write(*,*) 'q',minval(inp_q),maxval(inp_q), minval(stdlvl_spec_hum), maxval(stdlvl_spec_hum)
       stop
    end if

! 
! Interpolate temperature
!

    stdlvl_temp = stdlvl_t
    valid_values = stdlvl_t(:,:,start_level)
    if (any(valid_values < 155.0 .or. valid_values > 373.15 .or. &
            valid_values == fillvalue_t)) then
       write(*,*) "t panic attack", minval(valid_values),maxval(valid_values)
       stop
    endif

    do level = start_level + d_level, end_level, d_level
       where (stdlvl_temp(:,:,level) .eq. fillvalue_t)
          stdlvl_temp(:,:,level) = valid_values
       elsewhere
          valid_values = stdlvl_temp(:,:,level)
       endwhere
    end do

    stdlvl_internal_energy(:,:,integrand_first_level:integrand_last_level) = stdlvl_temp * stdlvl_pdel
    stdlvl_internal_energy(:,:,integrand_TOA_level) = stdlvl_temp(:,:,start_level) * stdlvl_pint(:,:,start_level)
    
    call cumsum(stdlvl_internal_energy, stdlvl_int_internal_energy, stdlvl_pverp)

    call interpolate(stdlvl_int_internal_energy, int_internal_energy, S2P, stdlvl_pverp+1, pverp)
        
    call difference( int_internal_energy,  internal_energy, pverp)

    temp = internal_energy / pdel

    inp_t = temp

    whacked = count(inp_t < 155.0 .or. inp_t > 373.15)
    if (whacked > 0) then
       write(*,*)  whacked, " instances of bogus t out of  ", plon * plat * pver
       write(*,*) minval(inp_t), maxval(inp_t), minval(stdlvl_t), maxval(stdlvl_t)
       stop
    end if

!
! Construct cloud condensate amounts
!
    where (inp_cloud > 0.0 & 
           .and. inp_tcldliq > 0 .and. inp_tcldice > 0)  !WDC 2/29/16 
!
! tcldliq/tcldice are in kg/kg: CESM expects icldlwp/iwp to be in gm/m^2
!
       !DRF 07/28/16 inp_icldlwp = inp_tcldliq / inp_cloud * (1000.0 / SHR_CONST_G) * &
       !DRF 07/28/16              abs(pdel) * 10000.0
       !DRF 07/28/16 inp_icldiwp = inp_tcldice / inp_cloud * (1000.0 / SHR_CONST_G) * &
       !DRF 07/28/16              abs(pdel) * 10000.0
       inp_icldlwp = inp_tcldliq / inp_cloud * (1000.0 / SHR_CONST_G) * &
                     abs(pdel)
       inp_icldiwp = inp_tcldice / inp_cloud * (1000.0 / SHR_CONST_G) * &
                     abs(pdel)


    elsewhere
       inp_icldlwp = 0.0
       inp_icldiwp = 0.0
       inp_tcldliq = 0.0
       inp_tcldice = 0.0
       inp_cloud   = 0.0
    end where

	!write(*,*) 'inp_tcldliq = ',inp_tcldliq(210,132,:)
	!write(*,*) 'inp_tcldice = ',inp_tcldice(210,132,:)
	!write(*,*) 'inp_cloud = ',inp_cloud(210,132,:)
	!write(*,*) 'const = ',SHR_CONST_G
	!write(*,*) 'pdel = ',abs(pdel(210,132,:))

    where (inp_tcldliq > 0 .or. inp_tcldice > 0) 
       inp_fice = inp_tcldice / (inp_tcldliq + inp_tcldice)
    elsewhere
       inp_fice = 0.0
    end where

    whacked = count(inp_cloud .lt. 0 .or. &
                    inp_tcldice .lt. 0 .or. inp_tcldliq .lt. 0 .or. &
                    inp_icldlwp .lt. 0 .or. inp_icldiwp .lt. 0 .or. &
                    inp_fice .lt. 0)
    if (whacked > 0) then
       write(*,*)  whacked, " instances of < 0 cloud amounts out of  ", plon * plat * pver
    end if

!
! Setup t_cld
!
    !inp_t_cld = stdlvl_t
    inp_t_cld = inp_t !DRF

! Close the files
    call wrap_close(nfid1)
    call wrap_close(nfid2)

    return

  end subroutine input_data

  subroutine input_times(path, nslice, nstep   ,dtime   ,mdbase  ,msbase  ,&
                         mbdate, mbsec)
!-----------------------------------------------------------------------
!
! Purpose:
! Input timing data required to reconstitute day_of_yr
!
! Method:
! Use NetCDF wrapper routines to input date and calendar day
! Output corresponds to following fields in history files:

! Author: W. Collins
!
!-----------------------------------------------------------------------
    implicit none
#include <netcdf.inc>

!
! Input arguments
!
    character(len=*), intent(in) :: path     ! path to data
!
! Output arguments
!
    integer, intent(out) ::  nslice      ! number of time slices in file
    integer, pointer     ::  nstep(:)    ! current time step
    integer, pointer     ::  date(:)     ! current date
    integer, pointer     ::  datesec(:)  ! current second in date
    real(r8), intent(out) :: dtime       ! length of time step (seconds)
    integer, intent(out) ::  mdbase      ! base day of run (e.g., 0)
    integer, intent(out) ::  msbase      ! base seconds of base day (e.g., 0)
    integer, intent(out) ::  mbdate      ! base date (yyyymmdd format) of run
    integer, intent(out) ::  mbsec       ! base seconds of base date (e.g., 0)
!
! Local variables
!
    integer :: nfid                      ! NetCDF id
    integer :: istat                     ! allocate status
    integer :: varid                     ! variable id
    integer :: dimid                     ! dimension id
    character (len = (nf_max_name)) :: dimname   ! dimension name
    character (len = 11) :: routine = "input_times"

!----------------------------------------------------------------------

    call wrap_open(path, nf_nowrite, nfid)

!
! Get number of time slices
!
    call wrap_inq_dimid(nfid, "time", dimid)
    call wrap_inq_dim(nfid, dimid, dimname, nslice)
!
! Allocate space for time steps, date, day_of_yr
!   

    allocate(nstep(nslice), stat = istat)
    call alloc_err(istat, routine, "nstep", nslice)

    allocate(date(nslice), stat = istat)
    call alloc_err(istat, routine, "date", nslice)

    allocate(datesec(nslice), stat = istat)
    call alloc_err(istat, routine, "datesec", nslice)

    allocate(day_of_yr(nslice), stat = istat)
    call alloc_err(istat, routine, "day_of_yr", nslice)

!
! Get the data
!
    call wrap_inq_varid(nfid, "date", varid)
    call wrap_get_var_int(nfid, varid, date)

    call wrap_inq_varid(nfid, "calday", varid)
    call wrap_get_var_realx(nfid, varid, day_of_yr)

    nstep(1) = 0;
    mdbase = 0
    msbase = 0
    mbdate = date(1)
    mbsec = 0
    dtime = 600

!
    call wrap_close(nfid)

    return

  end subroutine input_times


  subroutine input_vert_grid(path)
!----------------------------------------------------------------------- 
! 
! Purpose
! Reads in vertical grid data
! 
! Method: 
! Acquires following fields
!    hyai
!    hybi
!    hyam
!    hybm
!    p0
!WDC -- 11/7/15 -- CMIP5 OSSE a, a_bnds, b, b_bnds
! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------
    implicit none
#include <netcdf.inc>
#include <comhyb.h>
!
! Input arguments
!
    character(len=*), intent(in) :: path     ! path to data
!
! Local variables
!
    integer :: nfid                      ! NetCDF id
    integer :: varid                     ! variable id
    character (len = 15) :: routine = "input_vert_grid"
    real(r8) :: hybrid_a(pver) ! vertical coordinate formula term: a(k+1/2)  WDC -- 11/7/15 -- CMIP5 OSSE
    real(r8) :: hybrid_b(pver) ! vertical coordinate formula term: b(k+1/2)  WDC -- 11/7/15 -- CMIP5 OSSE
    real(r8) :: hybrid_a_bnds(2, pver) ! vertical coordinate formula term: a(k+1/2)  WDC -- 11/7/15 -- CMIP5 OSSE
    real(r8) :: hybrid_b_bnds(2, pver) ! vertical coordinate formula term: b(k+1/2)  WDC -- 11/7/15 -- CMIP5 OSSE
    real(r8) :: hy_temp(pverp)

!
! Get vertical coordinate info in comhyb common block for constructing
!     pressure fields
!
    call wrap_open(path, nf_nowrite, nfid)
   
    call wrap_inq_varid(nfid, "b", varid)
    call wrap_get_var_realx(nfid, varid, hybrid_b)

    flip_z = (hybrid_b(1) > hybrid_b(pver))
    if (flip_z) then
       z0 = pver
       z0p = pverp
       z1 = 1
       z1p = 1
       dz = -1
    else
       z0 = 1
       z0p = z0
       z1 = pver
       z1p = pverp
       dz = 1
    end if

    hybm = hybrid_b(z0:z1:dz)

    call wrap_inq_varid(nfid, "a", varid)
    call wrap_get_var_realx(nfid, varid, hybrid_a)

    hyam = hybrid_a(z0:z1:dz)

    call wrap_inq_varid(nfid, "cp0", varid)   !WDC -- 11/7/15 -- CMIP5 OSSE 
    call wrap_get_var_realx(nfid, varid, ps0) !WDC -- 11/7/15 -- CMIP5 OSSE 

    call wrap_inq_varid(nfid, "a_bnds", varid)   !WDC -- 11/7/15 -- CMIP5 OSSE 
    call wrap_get_var_realx(nfid, varid, hybrid_a_bnds) !WDC -- 11/7/15 -- CMIP5 OSSE 

    call wrap_inq_varid(nfid, "b_bnds", varid)   !WDC -- 11/7/15 -- CMIP5 OSSE 
    call wrap_get_var_realx(nfid, varid, hybrid_b_bnds) !WDC -- 11/7/15 -- CMIP5 OSSE 

    hy_temp(1:pver) = hybrid_a_bnds(1, 1:pver) !WDC -- 11/7/15 -- CMIP5 OSSE 
    hy_temp(pverp)  = hybrid_a_bnds(2, pver)   !WDC -- 11/7/15 -- CMIP5 OSSE 

    hyai = hy_temp(z0p:z1p:dz)

    hy_temp(1:pver) = hybrid_b_bnds(1, 1:pver) !WDC -- 11/7/15 -- CMIP5 OSSE 
    hy_temp(pverp)  = hybrid_b_bnds(2, pver)   !WDC -- 11/7/15 -- CMIP5 OSSE 

    hybi = hy_temp(z0p:z1p:dz)
    
    call wrap_inq_varid(nfid, "hyam", varid)   !WDC -- 11/7/15 -- CMIP5 OSSE 
    call wrap_get_var_realx(nfid, varid, stdlvl_hyam) !WDC -- 11/7/15 -- CMIP5 OSSE 

    call wrap_inq_varid(nfid, "hyai", varid)   !WDC -- 11/7/15 -- CMIP5 OSSE 
    call wrap_get_var_realx(nfid, varid, stdlvl_hyai) !WDC -- 11/7/15 -- CMIP5 OSSE 

    call wrap_inq_varid(nfid, "hybm", varid)   !WDC -- 11/7/15 -- CMIP5 OSSE 
    call wrap_get_var_realx(nfid, varid, stdlvl_hybm) !WDC -- 11/7/15 -- CMIP5 OSSE 

    call wrap_inq_varid(nfid, "hybi", varid)   !WDC -- 11/7/15 -- CMIP5 OSSE 
    call wrap_get_var_realx(nfid, varid, stdlvl_hybi) !WDC -- 11/7/15 -- CMIP5 OSSE 

    call wrap_inq_varid(nfid, "P0", varid)
    call wrap_get_var_realx(nfid, varid, stdlvl_ps0)

    call wrap_close(nfid)

    return

 end subroutine input_vert_grid

  subroutine dump_input_data
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Dump allocated arrays for input data 
! 
! Method: 
! Standard F90 deallocate calls
! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------

!
! Local variables
!
    integer :: istat           ! status variable
    character (len = 15) :: routine = "dump_input_data"

    if (allocated(inp_aermass))   then 
       deallocate(inp_aermass, stat = istat)
       call dealloc_err(istat, routine, "inp_aermass")
    endif
    if (allocated(inp_allaer))   then 
       deallocate(inp_allaer, stat = istat)
       call dealloc_err(istat, routine, "inp_allaer")
    endif
    if (allocated(inp_asdir))    then 
       deallocate(inp_asdir, stat = istat)
       call dealloc_err(istat, routine, "inp_asdir")
    endif
    if (allocated(inp_asdif))    then 
       deallocate(inp_asdif, stat = istat)
       call dealloc_err(istat, routine, "inp_asdif")
    endif
    if (allocated(inp_aldir))    then 
       deallocate(inp_aldir, stat = istat)
       call dealloc_err(istat, routine, "inp_aldir")
    endif
    if (allocated(inp_aldif))    then 
       deallocate(inp_aldif, stat = istat)
       call dealloc_err(istat, routine, "inp_aldif")
    endif
    if (allocated(inp_cfc11))    then 
       deallocate(inp_cfc11, stat = istat)
       call dealloc_err(istat, routine, "inp_cfc11")
    endif
    if (allocated(inp_cfc12))    then 
       deallocate(inp_cfc12, stat = istat)
       call dealloc_err(istat, routine, "inp_cfc12")
    endif
    if (allocated(inp_ch4))      then 
       deallocate(inp_ch4, stat = istat)
       call dealloc_err(istat, routine, "inp_ch4")
    endif
    if (allocated(inp_temp))    then 
       deallocate(inp_temp, stat = istat)
       call dealloc_err(istat, routine, "inp_temp")
    endif
    if (allocated(inp_cloud))    then 
       deallocate(inp_cloud, stat = istat)
       call dealloc_err(istat, routine, "inp_cloud")
    endif
    if (allocated(inp_emis))     then 
       deallocate(inp_emis, stat = istat)
       call dealloc_err(istat, routine, "inp_emis")
    endif
    if (allocated(inp_icldiwp))  then 
       deallocate(inp_icldiwp, stat = istat)
       call dealloc_err(istat, routine, "inp_icldiwp")
    endif
    if (allocated(inp_icldlwp))  then 
       deallocate(inp_icldlwp, stat = istat)
       call dealloc_err(istat, routine, "inp_icldlwp")
    endif
    if (allocated(inp_fice))  then 
       deallocate(inp_fice, stat = istat)
       call dealloc_err(istat, routine, "inp_fice")
    endif
    if (allocated(inp_tcldice))  then 
       deallocate(inp_tcldice, stat = istat)
       call dealloc_err(istat, routine, "inp_tcldice")
    endif
    if (allocated(inp_tcldliq))  then 
       deallocate(inp_tcldliq, stat = istat)
       call dealloc_err(istat, routine, "inp_tcldliq")
    endif
    if (allocated(inp_landfrac)) then 
       deallocate(inp_landfrac, stat = istat)
       call dealloc_err(istat, routine, "inp_landfrac")
    endif
    if (allocated(inp_landmcos)) then 
       deallocate(inp_landmcos, stat = istat)
       call dealloc_err(istat, routine, "inp_landmcos")
    endif
    if (allocated(inp_icefrac)) then 
       deallocate(inp_icefrac, stat = istat)
       call dealloc_err(istat, routine, "inp_icefrac")
    endif
    if (allocated(inp_snowh)) then 
       deallocate(inp_snowh, stat = istat)
       call dealloc_err(istat, routine, "inp_snowh")
    endif
    if (allocated(inp_n2o))      then 
       deallocate(inp_n2o, stat = istat)
       call dealloc_err(istat, routine, "inp_n2o")
    endif
    if (allocated(inp_o3vmr))    then 
       deallocate(inp_o3vmr, stat = istat)
       call dealloc_err(istat, routine, "inp_o3vmr")
    endif
    if (allocated(inp_ps))       then 
       deallocate(inp_ps, stat = istat)
       call dealloc_err(istat, routine, "inp_ps")
    endif
    if (allocated(stdlvl_q))        then 
       deallocate(stdlvl_q, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_q")
    endif
    if (allocated(stdlvl_rh))        then 
       deallocate(stdlvl_rh, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_rh")
    endif
    if (allocated(stdlvl_rh_fix))        then 
       deallocate(stdlvl_rh_fix, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_rh_fix")
    endif
    if (allocated(inp_rel))      then 
       deallocate(inp_rel, stat = istat)
       call dealloc_err(istat, routine, "inp_rel")
    endif
    if (allocated(inp_rei))      then 
       deallocate(inp_rei, stat = istat)
       call dealloc_err(istat, routine, "inp_rei")
    endif
    if (allocated(stdlvl_t))        then 
       deallocate(stdlvl_t, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_t")
    endif
    if (allocated(inp_q))        then 
       deallocate(inp_q, stat = istat)
       call dealloc_err(istat, routine, "inp_q")
    endif
    if (allocated(inp_rh))        then 
       deallocate(inp_rh, stat = istat)
       call dealloc_err(istat, routine, "inp_rh")
    endif
    if (allocated(inp_t))        then 
       deallocate(inp_t, stat = istat)
       call dealloc_err(istat, routine, "inp_t")
    endif
    if (allocated(inp_t_cld))        then 
       deallocate(inp_t_cld, stat = istat)
       call dealloc_err(istat, routine, "inp_t_cld")
    endif
    if (allocated(inp_lwup))       then 
       deallocate(inp_lwup, stat = istat)
       call dealloc_err(istat, routine, "inp_lwup")
    endif
    if (allocated(inp_ts))       then 
       deallocate(inp_ts, stat = istat)
       call dealloc_err(istat, routine, "inp_ts")
    endif
    if (allocated(inp_windspeed))       then 
       deallocate(inp_windspeed, stat = istat)
       call dealloc_err(istat, routine, "inp_windspeed")
    endif
    if (allocated(inp_lat))       then 
       deallocate(inp_lat, stat = istat)
       call dealloc_err(istat, routine, "inp_lat")
    endif
    if (allocated(inp_lon))       then 
       deallocate(inp_lon, stat = istat)
       call dealloc_err(istat, routine, "inp_lon")
    endif
    if (allocated(pdel))           then
       deallocate(pdel, stat = istat)
       call dealloc_err(istat, routine, "pdel")
    endif
    if (allocated(pint))           then
       deallocate(pint, stat = istat)
       call dealloc_err(istat, routine, "pint")
    endif
    if (allocated(stdlvl_ps))           then
       deallocate(stdlvl_ps, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_ps")
    endif
    if (allocated(stdlvl_pint))           then
       deallocate(stdlvl_pint, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_pint")
    endif
    if (allocated(stdlvl_pint_plus_toa))           then
       deallocate(stdlvl_pint_plus_toa, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_pint_plus_toa")
    endif
    if (allocated(stdlvl_pdel))           then
       deallocate(stdlvl_pdel, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_pdel")
    endif

    if (allocated(stdlvl_sat_spec_hum)) then
       deallocate(stdlvl_sat_spec_hum, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_sat_spec_hum")
    endif
    if (allocated(stdlvl_sat_vapor_mass)) then
       deallocate(stdlvl_sat_vapor_mass, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_sat_vapor_mass")
    endif
    if (allocated(stdlvl_int_sat_vapor_mass)) then
       deallocate(stdlvl_int_sat_vapor_mass, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_int_sat_vapor_mass")
    endif

    if (allocated(sat_spec_hum)) then
       deallocate(sat_spec_hum, stat = istat)
       call dealloc_err(istat, routine, "sat_spec_hum")
    endif
    if (allocated(sat_vapor_mass)) then
       deallocate(sat_vapor_mass, stat = istat)
       call dealloc_err(istat, routine, "sat_vapor_mass")
    endif
    if (allocated(int_sat_vapor_mass)) then
       deallocate(int_sat_vapor_mass, stat = istat)
       call dealloc_err(istat, routine, "int_sat_vapor_mass")
    endif

    if (allocated(stdlvl_spec_hum)) then
       deallocate(stdlvl_spec_hum, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_spec_hum")
    endif
    if (allocated(stdlvl_vapor_mass)) then
       deallocate(stdlvl_vapor_mass, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_vapor_mass")
    endif
    if (allocated(stdlvl_int_vapor_mass)) then
       deallocate(stdlvl_int_vapor_mass, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_int_vapor_mass")
    endif

    if (allocated(spec_hum)) then
       deallocate(spec_hum, stat = istat)
       call dealloc_err(istat, routine, "spec_hum")
    endif
    if (allocated(vapor_mass)) then
       deallocate(vapor_mass, stat = istat)
       call dealloc_err(istat, routine, "vapor_mass")
    endif
    if (allocated(int_vapor_mass)) then
       deallocate(int_vapor_mass, stat = istat)
       call dealloc_err(istat, routine, "int_vapor_mass")
    endif

!HERE
    if (allocated(stdlvl_temp)) then
       deallocate(stdlvl_temp, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_temp")
    endif
    if (allocated(stdlvl_internal_energy)) then
       deallocate(stdlvl_internal_energy, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_internal_energy")
    endif
    if (allocated(stdlvl_int_internal_energy)) then
       deallocate(stdlvl_int_internal_energy, stat = istat)
       call dealloc_err(istat, routine, "stdlvl_int_internal_energy")
    endif

    if (allocated(temp)) then
       deallocate(temp, stat = istat)
       call dealloc_err(istat, routine, "temp")
    endif
    if (allocated(internal_energy)) then
       deallocate(internal_energy, stat = istat)
       call dealloc_err(istat, routine, "internal_energy")
    endif
    if (allocated(int_internal_energy)) then
       deallocate(int_internal_energy, stat = istat)
       call dealloc_err(istat, routine, "int_internal_energy")
    endif

    if (allocated(valid_values)) then
       deallocate(valid_values, stat = istat)
       call dealloc_err(istat, routine, "valid_values")
    endif

    if (allocated(S2P)) then
       deallocate(S2P, stat = istat)
       call dealloc_err(istat, routine, "S2P")
    endif
!    if (allocated(P2S)) then
!       deallocate(P2S, stat = istat)
!       call dealloc_err(istat, routine, "P2S")
!    endif
    if (allocated(date)) then
       deallocate(date, stat = istat)
       call dealloc_err(istat, routine, "date")
    endif
    if (allocated(day_of_yr)) then
       deallocate(day_of_yr, stat = istat)
       call dealloc_err(istat, routine, "day_of_yr")
    endif

    return

  end subroutine dump_input_data

  subroutine cumsum(array_a, array_b, n)
    use shr_kind_mod, only: r8 => shr_kind_r8 
    use pmgrid

    implicit none

    integer  :: n
    real(r8), intent(in)    :: array_a(plon, plat, n)
    real(r8), intent(inout) :: array_b(plon, plat, n+1)
    
    integer  :: i

    array_b(:,:,1) = 0.0
    do i = 1, n
       array_b(:,:,i+1) = array_b(:,:,i) + array_a(:,:,i)
    end do

  end subroutine cumsum

  subroutine difference(array_a, array_b, np1)
    use shr_kind_mod, only: r8 => shr_kind_r8 
    use pmgrid

    implicit none

    integer  :: np1
    real(r8), intent(in)    :: array_a(plon, plat, np1)
    real(r8), intent(inout) :: array_b(plon, plat, np1-1)
    
    integer  :: i

    array_b(:,:,1:np1-1) = array_a(:,:,2:np1) - array_a(:,:,1:np1-1)

  end subroutine difference

  subroutine interpolate(array_a, array_b, wgts, np1_a, np1_b)
    use shr_kind_mod, only: r8 => shr_kind_r8 
    use pmgrid

    implicit none

    integer, intent(in) :: np1_a
    integer, intent(in) :: np1_b
    real(r8), intent(in) :: array_a(plon, plat, np1_a)
    real(r8), intent(inout) :: array_b(plon, plat, np1_b)
    real(r8), intent(in) :: wgts(plon, plat, np1_a, np1_b)

    integer :: i

    do i = 1, np1_b
       array_b(:,:,i) = sum(array_a * wgts(:,:,:,i),3)
    end do

  end subroutine interpolate

  subroutine build_interpolation_matrix(wgts, pint_from, np1_from, pint_to, np1_to,print_diag)
    use pmgrid
    use shr_kind_mod, only: r8 => shr_kind_r8                                   

    implicit none

    real(r8), intent(in)    :: pint_from(plon, plat, np1_from)
    real(r8), intent(in)    :: pint_to(plon, plat, np1_to)
    integer, intent(in)     :: np1_from
    integer, intent(in)     :: np1_to
    real(r8), intent(inout) :: wgts(plon, plat, np1_from, np1_to)
    logical, intent(in)     :: print_diag

    integer :: ilon
    integer :: ilat
    integer :: ifrom
    integer :: ito
    integer :: min_from
    integer :: max_from
    logical :: intersection
    logical :: from_to_same_signed
    logical :: p_maximum
    logical :: p_minimum
    real(r8) :: product
    real(r8) :: d
    real(r8) :: d1

    wgts = 0.0

    do ilon = 1, plon
       do ilat = 1, plat
          min_from = 1
          max_from = np1_from - 1
          from_to_same_signed = &
               ((pint_from(ilon,ilat,np1_from) - pint_from(ilon,ilat,1)) * &
                (pint_to(ilon,ilat,np1_to)     - pint_to(ilon,ilat,1))) .gt. 0.0
          do ito = 1, np1_to
             intersection = .false.
             do ifrom = min_from, max_from
                d =  (pint_to(ilon, ilat, ito)       - pint_from(ilon, ilat, ifrom))
                d1 = (pint_from(ilon, ilat, ifrom+1) - pint_to(ilon, ilat, ito))
                product = d * d1
                if (product .ge. 0.0) then
                   if (print_diag) write(*,*) pint_to(ilon, ilat, ito),pint_from(ilon, ilat, ifrom),pint_from(ilon, ilat, ifrom+1)
                   intersection = .true.
                   wgts(ilon, ilat, ifrom, ito)   = d1 / (d1 + d)
                   wgts(ilon, ilat, ifrom+1, ito) = d  / (d1 + d)
                   if (from_to_same_signed) then
                      min_from = ifrom
                   else
                      max_from = ifrom
                   endif
                   exit
                endif
             end do
!logic for edge cases
             if (intersection .eq. .false.) then
                p_minimum = (pint_to(ilon, ilat, ito) .lt. minval(pint_from(ilon, ilat, :)))
                p_maximum = (pint_to(ilon, ilat, ito) .gt. maxval(pint_from(ilon, ilat, :)))
                if (p_minimum) then 
                   wgts(ilon, ilat, minloc(pint_from(ilon, ilat, :)), ito) = 1.0
                else if (p_maximum) then 
                   wgts(ilon, ilat, maxloc(pint_from(ilon, ilat, :)), ito) = 1.0   
                else
                   write(*,*) " Error in interpolation array construction"
                   stop
                endif
             endif
          end do
       end do
    end do

  end subroutine build_interpolation_matrix

  subroutine input_landmcos(ncdata)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Input landm_coslat field from initial conditions files
! 
! Method: 
! Use NetCDF wrapper routines to input landm_coslat
! Output corresponds to following fields in initial files:
!     landm_coslat
!
! Author: W. Collins
! 
!-----------------------------------------------------------------------
    implicit none
#include <netcdf.inc>

!
! Input arguments
!
    character(len=*), intent(in) :: ncdata     ! path to data
!
! Output arguments
!
!
! Local variables
!
    integer :: nfid                      ! NetCDF id
    integer :: istat                     ! allocate status
    integer :: varid                     ! variable id

    character (len = 14) :: routine = "input_landmcos"

    double precision, allocatable :: lat(:) ! Latitude
    integer :: latvarid                  ! Latitude variable ID
    double precision, parameter :: DEG2RAD = 3.141592653589793 / 180.0
    
!----------------------------------------------------------------------

    call wrap_open(ncdata, nf_nowrite, nfid)

!
! Allocate space for landm_coslat
!
    allocate(inp_landmcos(plon,plat), stat = istat) 
    call alloc_err(istat, routine, "inp_landmcos",(plon*plat) )
!
! Get the landm_coslat data
!
    call wrap_inq_varid(nfid, "LANDFRAC", varid)
    call wrap_get_var_realx(nfid, varid, inp_landmcos)
!
! Allocate space for latitude
!
    allocate(lat(plat), stat = istat)
    call alloc_err(istat, routine, "lat", (plat))
!
! Get the data
!
    call wrap_inq_varid(nfid, "lat", latvarid)
    call wrap_get_var_realx(nfid, latvarid, lat)
!
! Weight the landfraction
!
    inp_landmcos = inp_landmcos * spread(cos(lat * DEG2RAD), 1, plon)
!
! set to fraction between 0 and 1
!
    inp_landmcos = inp_landmcos * 10**(-nint(maxval(log10(inp_landmcos))))
!
! Deallocate space for latitutde
!
    deallocate(lat, stat = istat)
    call alloc_err(istat, routine, "lat",(plon*plat) )

    call wrap_close(nfid)
    
    return

  end subroutine input_landmcos

end module RadInput

