#include <misc.h>
#include <params.h>

module get_surface_properties

  !Purpose:
  !   get surface properties including land/ocean and brdf
  !   Auther: D. Feldman

  !  brdf properties from MODIS MCD43C1 product

  use ppgrid
  use Chunks, only: c_landfrac
  use pmgrid
  use filenames, only: brdf_file,ocean_refl,snow_refl,alb_file,snow_brdf,alb_snow,snow_frac,emis_file
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  !quantities of interest
  logical, public :: land_flag(plon,plat)
  real(r8), public :: fsno_read(plon,plat)
  real(r8), public :: brdf_param(plon,7,3,plat) !Ross-Li BRDF parameters from Modis
  real(r8), public :: emis_array(plon,6,plat) !Ross-Li BRDF parameters from Modis
  real(r8), public :: ocean_reflectance(24,2) !column 1 = wavelength,2=reflectance
  real(r8), public :: snow_reflectance(969,4) 
  real(r8), public :: modis_asdir(12,plat,plon)
  real(r8), public :: modis_asdif(12,plat,plon)
  real(r8), public :: modis_aldir(12,plat,plon)
  real(r8), public :: modis_aldif(12,plat,plon)
  real(r8), public :: snow_brdf_landtype(7,3,16)
  integer, public :: landtype(plon,plat)
  !quantities of interest

  save

contains
  subroutine get_emis(c_lat,c_lon,c_snowh,c_icefrac,mcdate_input)

    use ioFileMod, only: getfil
    implicit none

    include 'netcdf.inc'

    !Input arguments !DRF
    real(r8), intent(in) :: c_lat(plat)          ! Longwave up flux at surface
    real(r8), intent(in) :: c_lon(plon)          ! Cloud emissivity
    real(r8), intent(in) :: c_snowh(plon,plat)   ! snow height
    real(r8), intent(in) :: c_icefrac(plon,plat) ! sea-ice fraction
    real(r8), intent(in) ::  mcdate_input         ! current date (yyyymmdd format) (e.g., 021105)

    !local variables
    integer :: i,j,k,m,n
    integer :: lat_id, lon_id
    integer :: index1,index2

    !variables from emis_file
    integer :: emis_id,mon_id,band_id,count_id
    integer :: num_lat,num_lon,num_mon,num_band
        
    character(len=256) :: locfn   !local file
    integer :: nc_id
    real(r8) :: coszrs  !cosine solar zenith angle
    real(r8), allocatable :: latitude_emis(:)
    real(r8), allocatable :: rh_opac(:) 
    real(r8), allocatable :: longitude(:)
    real(r8), allocatable :: dummy_array(:)

    integer :: jday_mid(12)
    integer :: fsno_id,fsno_index
    real(r8) :: fsno_index1,fsno_index2
    integer :: month_selected
    integer :: start2(4)
    integer :: count2(4)

    !!!!!!!!!!!!!!!!!!!!!
    write (6, '(2x, a)') '_______________________________________________________'
    write (6, '(2x, a)') '_______ initializing land emissivity properties _______'
    write (6, '(2x, a)') '_______________________________________________________'

    call getfil(emis_file, locfn, 0)
    call wrap_open(locfn, 0, nc_id)
    call wrap_inq_dimid(nc_id, 'latitude', lat_id)
    call wrap_inq_dimid(nc_id, 'longitude', lon_id)
    call wrap_inq_dimid(nc_id, 'month',mon_id)
    call wrap_inq_dimid(nc_id, 'wavelength',band_id)

    call wrap_inq_dimlen(nc_id, lat_id, num_lat)
    call wrap_inq_dimlen(nc_id, lon_id, num_lon)
    call wrap_inq_dimlen(nc_id, mon_id, num_mon)
    call wrap_inq_dimlen(nc_id, band_id, num_band)
    
    allocate(latitude_emis(num_lat))
    allocate(longitude(num_lon))
    allocate(dummy_array(num_lon))

    call wrap_inq_varid(nc_id, "emis", emis_id)
    
    call wrap_get_var_realx(nc_id,lat_id,latitude_emis)
    call wrap_get_var_realx(nc_id,lon_id,longitude)
   
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

    !fsno_index1 = floor((mcdate_input-20000000.)/10000.)*12.
    !fsno_index2 = floor((mcdate_input-20000000.-fsno_index1*10000./12.)/100.)
    fsno_index1 = mod(mcdate_input,10000.)
    fsno_index2 = (fsno_index1 - mod(fsno_index1,100.))/100.

    month_selected = int(fsno_index2)
    if (month_selected.eq.0) then
       month_selected = 12
    endif

    write(*,*) 'mcdate_input = ',mcdate_input
    write(*,*) 'fsno_index1 = ',fsno_index1
    write(*,*) 'fsno_index2 = ',fsno_index2
    write(*,*) 'month_selected = ',month_selected
    write(*,*) 'num_lon = ',num_lon

    do j=1,num_band
      do m=1,num_lat
         start2 = (/1,m,j,month_selected/)
         count2 = (/num_lon,1,1,1/)
         call wrap_get_vara_realx(nc_id, emis_id, start2, count2, dummy_array)             
         emis_array(:,j,m) = dummy_array
      end do
    end do
    deallocate(dummy_array)


  end subroutine get_emis

  subroutine get_brdf(c_lat,c_lon,c_snowh,c_icefrac,mcdate_input)
    
    use ioFileMod, only: getfil
    implicit none

    include 'netcdf.inc'

    !Input arguments !DRF
    real(r8), intent(in) :: c_lat(plat)          ! Longwave up flux at surface
    real(r8), intent(in) :: c_lon(plon)          ! Cloud emissivity
    real(r8), intent(in) :: c_snowh(plon,plat)   ! snow height
    real(r8), intent(in) :: c_icefrac(plon,plat) ! sea-ice fraction
    real(r8), intent(in) ::  mcdate_input         ! current date (yyyymmdd format) (e.g., 021105)

    !local variables
    integer :: i,j,k,m,n
    integer :: lat_id, lon_id
    integer :: index1,index2
        
    character(len=256) :: locfn   !local file
    integer :: nc_id
    real(r8) :: coszrs  !cosine solar zenith angle
    real(r8), allocatable :: latitude_brdf(:)
    real(r8), allocatable :: rh_opac(:) 
    real(r8), allocatable :: longitude(:)
    real(r8), allocatable :: brdf_array(:,:,:,:)
    real(r8), allocatable :: snow_brdf_array(:,:,:,:)
    real(r8), allocatable :: dummy_array(:)
    real(r8), allocatable :: dummy_array3(:)
    real(r8), allocatable :: dummy_array4(:)
    real(r8), allocatable :: dummy_snowbrdf(:)
    real(r8), allocatable :: local_landfrac(:,:)
    integer, allocatable :: dummy_landtype(:)
    real(r8), allocatable :: snow_asdir(:,:)
    real(r8), allocatable :: snow_asdif(:,:)
    real(r8), allocatable :: snow_aldir(:,:)
    real(r8), allocatable :: snow_aldif(:,:)
    !real(r8), allocatable :: fsno_read(:,:)

    !variables from brdf_file
    integer :: brdf_id,param_id,mon_id,band_id,count_id,landtype_id
    integer :: num_lat,num_lon,num_mon,num_band,num_param

    !variables from ocean_refl file
    integer :: wvl_id,refl_id
    integer :: asdir_id,asdif_id,aldir_id,aldif_id
    integer :: fsno_id,fsno_index
    real(r8) :: fsno_index1,fsno_index2
    integer :: month_selected
    integer :: num_wvl,num_refl
    integer :: start3(2)
    integer :: count3(2)

    integer :: start(5)
    integer :: count(5)

    integer :: start2(3)
    integer :: count2(3)

    integer :: start4(4)
    integer :: count4(4)

    integer :: jday_mid(12)
    integer :: ii, jj

    real(r8) :: solar_az                    ! Solar azimuth angle (degrees) Produced from Wiscombe calculator
    real(r8) :: solar_zen                   ! Solar zenith angle (degrees) Produced from Wiscombe calculator
    real(r8) :: solar_zen_rad               ! Solar zenith angle (radians) Produced from Wiscombe calculator
    real(r8) :: half_pi                     ! constant (pi/2)
    real(r8) :: solar_hour                  ! Local solar hour used for Wiscombe calculator
    real(r8) :: mcdate_local                 ! Year used for Wiscombe calculator
    integer :: mcdate_local2                ! Integer year used for Wiscombe calculator
    real(r8) :: dummy_soldia                ! Solar diameter from Wiscombe calculator (not used)
    real(r8) :: dummy_soldst                ! Solar distance from Wiscombe calculator (not used)
    real(r8) :: sun_calc_long               ! Longitude value (degrees) for Wiscombe calculator

    ! variables for 6S brdf integration
    real*4 :: pws   !wind-speed (m/s)
    real*4 :: paw   !azimuth argument over which to be integrated
    real*4 :: xsal  !salt concentration (ppt)
    real*4 :: pwl   !wavelength (um)
    real*4 :: pcl   !chlorophyll concentration (mg/m^3)
    real*4 :: rfoam !unused output of 6S brdf program
    real*4 :: rwat  !unused output of 6S brdf program
    real*4 :: rglit !unused output of 6S brdf program
    integer :: mu    !number of zenith angles being input
    integer :: np    !number of aziumth angles being input
    real*4 :: rp     ! ??
    real*4 :: rm(-1:1) !angles over which to be integrated
    real*4 :: brdfint(-1:1) !output of oceabrdf routine

    real(r8) :: gws_iso  ! multiplication factor for white-sky albedo isotropic factor (DRF)
    real(r8) :: gws_vol  ! multiplication factor for white-sky albedo volumetric factor (DRF)
    real(r8) :: gws_geo  ! multiplication factor for white-sky albedo geometric factor (DRF)
    real(r8) :: gbs_iso(3)  ! multiplication factor for black-sky albedo isotropic factor (DRF)
    real(r8) :: gbs_vol(3)  ! multiplication factor for black-sky albedo volumetric factor (DRF)
    real(r8) :: gbs_geo(3)  ! multiplication factor for black-sky albedo geometric factor (DRF)
    real(r8) :: brdf_param2(7,3) !brdf parameters for albedo integration (both black- and white-sky) (DRF)
    real(r8) :: brdf_wvl(7) !brdf wavelength grid (DRF)
    real(r8) :: brdf_mult_vis(7) !multiplying factor for wavelength-integration over visible (DRF)
    real(r8) :: brdf_mult_nir(7) !multiplying factor for wavelength-integration over near-ir (DRF)
    real(r8) :: spec_alb_ws(7)  !white-sky spectral albedo (DRF)
    real(r8) :: spec_alb_bs(7)  !black-sky spectral albedo (DRF)
    
    logical :: dummy
    real(r8) :: dummy2

    write (6, '(2x, a)') '_______________________________________________________'
    write (6, '(2x, a)') '_______ initializing land surface properties __________'
    write (6, '(2x, a)') '_______________________________________________________'

    call getfil(brdf_file, locfn, 0)
    call wrap_open(locfn, 0, nc_id)
    call wrap_inq_dimid(nc_id, 'latitude', lat_id)
    call wrap_inq_dimid(nc_id, 'longitude', lon_id)
    call wrap_inq_dimid(nc_id, 'parameter_number',param_id)
    call wrap_inq_dimid(nc_id, 'month',mon_id)
    call wrap_inq_dimid(nc_id, 'wavelength',band_id)

    call wrap_inq_dimlen(nc_id, lat_id, num_lat)
    call wrap_inq_dimlen(nc_id, lon_id, num_lon)
    call wrap_inq_dimlen(nc_id, param_id, num_param)
    call wrap_inq_dimlen(nc_id, mon_id, num_mon)
    call wrap_inq_dimlen(nc_id, band_id, num_band)
    
    allocate(latitude_brdf(num_lat))
    allocate(longitude(num_lon))
    allocate(brdf_array(num_band,num_param,num_lat,num_lon))
    allocate(dummy_array(num_lon))
    allocate(local_landfrac(num_lon,num_lat))

    call wrap_inq_varid(nc_id, "brdf", brdf_id)
    
    call wrap_get_var_realx(nc_id,lat_id,latitude_brdf)
    call wrap_get_var_realx(nc_id,lon_id,longitude)
   
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

    !Get month of selection
    write(*,*) 'mcdate_input = ',mcdate_input
    fsno_index1 = mod(mcdate_input,10000.)
    fsno_index2 = (fsno_index1 - mod(fsno_index1,100.))/100.

    !fsno_index1 = floor((mcdate_input-20000000.)/10000.)*12.
    !fsno_index2 = floor((mcdate_input-20000000.-fsno_index1*10000./12.)/100.)

    month_selected = int(fsno_index2)
    if (month_selected.eq.0) then
       month_selected = 12
    endif

    write(*,*) 'mcdate_input = ',mcdate_input
    write(*,*) 'fsno_index1 = ',fsno_index1
    write(*,*) 'fsno_index2 = ',fsno_index2
    write(*,*) 'month_selected = ',month_selected
    write(*,*) 'num_lon = ',num_lon

    do j=1,num_band
      do k=1,num_param
         do m=1,num_lat
            start = (/1,m,k,j,month_selected/)
            count = (/num_lon,1,1,1,1/)
            call wrap_get_vara_realx(nc_id, brdf_id, start, count, dummy_array)             
            brdf_array(j,k,m,:) = dummy_array

	    !if (m.eq.166) then
!		write(*,*) 'brdf_array read = ',brdf_array(j,k,m,166)
!	    endif

         end do
      end do
    end do
 
    !Reassign to format of brdf_param array
    do i=1,plon
       index1 = minloc(abs(c_lon(i)-longitude),1)
       do j=1,plat
          !Get lat/lon
          index2 = minloc(abs(c_lat(j)-latitude_brdf),1)
          land_flag(i,j) = .false.

          !use cam land fraction

          local_landfrac(i,j) = c_landfrac(i,j) !for debugging purposes
          if (c_landfrac(i,j)<0.75) then
             land_flag(i,j) = .false.
          else
             land_flag(i,j) = .true.
          end if

          !Read in brdf information
          do m=1,num_band
             !Switch 2nd and 3rd Ross-Li parameters
             brdf_param(i,m,1,j) = brdf_array(m,1,index2,index1)    !iso
             brdf_param(i,m,2,j) = brdf_array(m,3,index2,index1) !flipped !geo
             brdf_param(i,m,3,j) = brdf_array(m,2,index2,index1)  !vol
          end do
       end do
    end do

    !write(*,*) 'brdf_param1 read',brdf_param(166,:,1,166)
    !write(*,*) 'brdf_param2 read',brdf_param(166,:,2,166)
    !write(*,*) 'brdf_param3 read',brdf_param(166,:,3,166)

    !free up resources
    deallocate(latitude_brdf)
    deallocate(longitude)
    deallocate(brdf_array)
    deallocate(dummy_array)
    deallocate(local_landfrac)

    call wrap_close(nc_id)

    !!!!!!!!!!!!!!!!!!!!!
    write (6, '(2x, a)') '_______________________________________________________'
    write (6, '(2x, a)') '_______ initializing ocean surface properties _________'
    write (6, '(2x, a)') '_______________________________________________________'

    !Now get data for ocean refl
    call getfil(ocean_refl, locfn, 0)
    call wrap_open(locfn, 0, nc_id)
    call wrap_inq_dimid(nc_id, 'wavelength', wvl_id)
    call wrap_inq_varid(nc_id, 'reflectance', refl_id)

    call wrap_inq_dimlen(nc_id, wvl_id, num_wvl)
    allocate(dummy_array3(num_wvl))

    call wrap_get_var_realx(nc_id,wvl_id,dummy_array3)
    ocean_reflectance(:,1) = dummy_array3

    start3 = (/1,1/)
    count3 = (/1,24/)
    call wrap_get_vara_realx(nc_id,refl_id,start3,count3,dummy_array3)
    ocean_reflectance(:,2) = dummy_array3
    deallocate(dummy_array3)
    call wrap_close(nc_id)

    !!!!!!!!!!!!!!!!!!!!!
    write (6, '(2x, a)') '_______________________________________________________'
    write (6, '(2x, a)') '_______ initializing ice/snow surface properties ______'
    write (6, '(2x, a)') '_______________________________________________________'

    call getfil(snow_refl, locfn, 0)
    call wrap_open(locfn, 0, nc_id)
    call wrap_inq_dimid(nc_id, 'wavelength', wvl_id)
    call wrap_inq_varid(nc_id, 'reflectance', refl_id)

    call wrap_inq_dimlen(nc_id, wvl_id, num_wvl)
    allocate(dummy_array4(num_wvl))

    call wrap_get_var_realx(nc_id,wvl_id,dummy_array4)
    snow_reflectance(:,1) = dummy_array4

    !fine,medium,coarse grain size
    start3 = (/1,1/)
    count3 = (/1,969/)
    call wrap_get_vara_realx(nc_id,refl_id,start3,count3,dummy_array4)
    snow_reflectance(:,2) = dummy_array4

    start3 = (/2,1/)
    count3 = (/1,969/)
    call wrap_get_vara_realx(nc_id,refl_id,start3,count3,dummy_array4)
    snow_reflectance(:,3) = dummy_array4
  
    start3 = (/3,1/)
    count3 = (/1,969/)
    call wrap_get_vara_realx(nc_id,refl_id,start3,count3,dummy_array4)
    snow_reflectance(:,4) = dummy_array4
    deallocate(dummy_array4)
    call wrap_close(nc_id)


    !!!!!!!!!!!!!!!!!!!!!
    write (6, '(2x, a)') '_______________________________________________________'
    write (6, '(2x, a)') '_______ read in snow-fraction--------------------------'
    write (6, '(2x, a)') '_______________________________________________________'

    call getfil(snow_frac, locfn, 0)
    call wrap_open(locfn, 0, nc_id)

    call wrap_inq_varid(nc_id, 'snc', fsno_id)
    !allocate(fsno_read(plon,plat))

    fsno_index1 = mod(mcdate_input,10000.)
    fsno_index2 = (fsno_index1 - mod(fsno_index1,100.))/100.
    fsno_index = int(fsno_index2)

    !fsno_index1 = floor((mcdate_input-20000000.)/10000.)*12.
    !fsno_index2 = floor((mcdate_input-20000000.-fsno_index1*10000./12.)/100.)
    !fsno_index = int(fsno_index1+fsno_index2)

    write(*,*) 'fsno_index = ',fsno_index

    !Now recombine brdf values
    start2 = (/1,   1,   fsno_index/) 
    count2 = (/plon,plat,1/) 
    call wrap_get_vara_realx(nc_id, fsno_id, start2, count2, fsno_read)
    call wrap_close(nc_id)

    fsno_read = fsno_read / 100.0  !DRF (07/26/16)

    !!!!!!!!!!!!!!!!!!!!!
    write (6, '(2x, a)') '_______________________________________________________'
    write (6, '(2x, a)') '_______ initializing ice/snow brdf land type info _____'
    write (6, '(2x, a)') '_______________________________________________________'

    allocate(snow_brdf_array(num_band,num_param,num_lat,num_lon))
    allocate(dummy_array(num_lon))
    allocate(latitude_brdf(num_lat))
    allocate(longitude(num_lon))

    !write(*,*) 'snow_brdf = ',snow_brdf
    !write(*,*) 'brdf_file = ',brdf_file
    !write(*,*) 'ocean_refl = ',ocean_refl
    !write(*,*) 'snow_refl = ',snow_refl
    !write(*,*) 'alb_file = ',alb_file
    !write(*,*) 'alb_snow = ',alb_snow
    !write(*,*) 'snow_frac = ',snow_frac

    call getfil(snow_brdf, locfn, 0)
    call wrap_open(locfn, 0, nc_id)
    call wrap_inq_varid(nc_id, 'brdf', brdf_id)
    call wrap_inq_dimid(nc_id, 'latitude', lat_id)
    call wrap_inq_dimid(nc_id, 'longitude', lon_id)
    call wrap_inq_dimlen(nc_id, lat_id, num_lat)
    call wrap_inq_dimlen(nc_id, lon_id, num_lon)

    call wrap_get_var_realx(nc_id,lat_id,latitude_brdf)
    call wrap_get_var_realx(nc_id,lon_id,longitude)

    do j=1,num_band
       do k=1,num_param
          do m=1,num_lat
             start = (/1,m,k,j,month_selected/)
             count = (/num_lon,1,1,1,1/)
             call wrap_get_vara_realx(nc_id, brdf_id, start, count, dummy_array) 


             snow_brdf_array(j,k,m,:) = dummy_array
          end do
       end do
    end do
    call wrap_close(nc_id)

    !Reassign to format of brdf_param array
    do i=1,plon
       index1 = minloc(abs(c_lon(i)-longitude),1)
       do j=1,plat
          !Get lat/lon
          index2 = minloc(abs(c_lat(j)-latitude_brdf),1)

          !Read in brdf information
          do m=1,num_band
              !Switch 2nd and 3rd Ross-Li parameters

              brdf_param(i,m,1,j) = (1.-fsno_read(i,j))*brdf_param(i,m,1,j)+fsno_read(i,j)*snow_brdf_array(m,1,index2,index1) !iso
              brdf_param(i,m,2,j) = (1.-fsno_read(i,j))*brdf_param(i,m,2,j)+fsno_read(i,j)*snow_brdf_array(m,3,index2,index1) !flipped geo
              brdf_param(i,m,3,j) = (1.-fsno_read(i,j))*brdf_param(i,m,3,j)+fsno_read(i,j)*snow_brdf_array(m,2,index2,index1) !vol

           end do
       end do
    end do

    !Deallocate data
    deallocate(snow_brdf_array)
    deallocate(dummy_array)
    deallocate(latitude_brdf)
    deallocate(longitude)

    !!!!!!!!!!!!!!!!!!!!!
    write (6, '(2x, a)') '_______________________________________________________'
    write (6, '(2x, a)') '_______ initializing MODIS Snow Albedo For Landtype----'
    write (6, '(2x, a)') '_______________________________________________________'

    call getfil(alb_snow, locfn, 0)
    call wrap_open(locfn, 0, nc_id)
    call wrap_inq_varid(nc_id, 'Albedo_BSA_vis', asdir_id)
    call wrap_inq_varid(nc_id, 'Albedo_WSA_vis', asdif_id)
    call wrap_inq_varid(nc_id, 'Albedo_BSA_nir', aldir_id)
    call wrap_inq_varid(nc_id, 'Albedo_WSA_nir', aldif_id)
    allocate(dummy_array4(16))

    allocate(snow_asdir(12,16))
    allocate(snow_asdif(12,16))
    allocate(snow_aldir(12,16))
    allocate(snow_aldif(12,16))
  
    do i=1,12
      start3 = (/1,i/)
      count3 = (/16,1/)
      call wrap_get_vara_realx(nc_id,asdir_id,start3,count3,dummy_array4)
      snow_asdir(i,:) = dummy_array4
      
      call wrap_get_vara_realx(nc_id,asdif_id,start3,count3,dummy_array4)
      snow_asdif(i,:) = dummy_array4

      call wrap_get_vara_realx(nc_id,aldif_id,start3,count3,dummy_array4)
      snow_aldif(i,:) = dummy_array4
      
      call wrap_get_vara_realx(nc_id,aldir_id,start3,count3,dummy_array4)
      snow_aldir(i,:) = dummy_array4
    end do
    deallocate(dummy_array4)
    call wrap_close(nc_id)

    !!!!!!!!!!!!!!!!!!!!!
    write (6, '(2x, a)') '_______________________________________________________'
    write (6, '(2x, a)') '_______ initializing MODIS Snow-Free Albedo properties_'
    write (6, '(2x, a)') '_______________________________________________________'

    call getfil(alb_file, locfn, 0)
    call wrap_open(locfn, 0, nc_id)
    call wrap_inq_varid(nc_id, 'ASDIR_r', asdir_id)
    call wrap_inq_varid(nc_id, 'ASDIF_r', asdif_id)
    call wrap_inq_varid(nc_id, 'ALDIR_r', aldir_id)
    call wrap_inq_varid(nc_id, 'ALDIF_r', aldif_id)
    allocate(dummy_array4(plon))

    do i=1,12
      do j=1,plat
      start2 = (/1,j,i/)
      count2 = (/plon,1,1/)
      call wrap_get_vara_realx(nc_id,asdir_id,start2,count2,dummy_array4)
      modis_asdir(i,j,:) = dummy_array4
      call wrap_get_vara_realx(nc_id,asdif_id,start2,count2,dummy_array4)
      modis_asdif(i,j,:) = dummy_array4
      call wrap_get_vara_realx(nc_id,aldir_id,start2,count2,dummy_array4)
      modis_aldir(i,j,:) = dummy_array4
      call wrap_get_vara_realx(nc_id,aldif_id,start2,count2,dummy_array4)
      modis_aldif(i,j,:) = dummy_array4

      !replace modis data with angular dependence for ocean albedo 
      !also replace snow/ice albedo information

      do m=1,plon
         mcdate_local = 2000.
         mcdate_local2 = 2000
         solar_hour = 13.5
         sun_calc_long = 0.
 
         call sun_calc(mcdate_local2,jday_mid(i),solar_hour,c_lat(j),sun_calc_long, &
                         solar_az,solar_zen,dummy_soldia,dummy_soldst)

            !if (i.eq.55 .and. j.eq.87) then
            !write(*,*) 'landtype = ',landtype(i,j)
            !write(*,*) 'brdf1 = ',brdf_param2(:,1)
            !write(*,*) 'brdf2 = ',brdf_param2(:,2)
            !write(*,*) 'brdf3 = ',brdf_param2(:,3)
            !stop 'error get_surf.F90 line 423'
            !endif

         coszrs = cos(solar_zen*3.1415926535898/180.)
         solar_zen_rad = solar_zen*3.1415926535898/180.
         half_pi = 3.1415926535898/2.

         if (.not.land_flag(m,j)) then
            !write(*,*) "changing modis albedo"  !!REMOVE DRF!!

            call albocean(coszrs,modis_asdir(i,j,m),modis_aldir(i,j,m), &
                       modis_asdif(i,j,m),modis_aldif(i,j,m))

            !if (i.eq.1 .and. m.eq.1 .and. j.eq.119) then
            !   write(*,*) 'zen 119 = ',solar_zen
            !   write(*,*) 'c_lat = ',c_lat(j)
            !   write(*,*) 'modis_asdir = ',modis_asdir(i,j,m)
            !   write(*,*) 'modis_asdif = ',modis_asdif(i,j,m)
            !   stop 'error get_surface_properties.F90 line 380'
            !endif
            !if (i.eq.4 .and. m.eq.137 .and. j.eq.30) then
            !   write(*,*) 'zen 22 = ',coszrs 
            !   write(*,*) 'modis_asdir = ',modis_asdir(i,j,m)
            !   write(*,*) 'modis_asdif = ',modis_asdif(i,j,m)
            !endif

            !write(*,*) "coszrs = ",coszrs
            !write(*,*) "modis_asdir = ",modis_asdir(i,j,m)
            !write(*,*) "modis_aldir = ",modis_aldir(i,j,m)
            !write(*,*) "modis_asdif = ",modis_asdif(i,j,m)
            !write(*,*) "modis_aldif = ",modis_aldif(i,j,m)
            !stop

            !write(*,*) "changing modis albedo"  !!REMOVE DRF!!

            !pws = 5.
            !xsal = 35.
            !pcl = 0.
            !pwl = 0.55
            !mu = 1
            !np = 1

            !paw = 30. !degrees
            !rm(-1) = paw*3.1415927/180.
            !rm(0) = 0.5
            !rm(1) = 0.5
            !rp = rm(-1)

            !call oceabrdf(pws,paw,xsal,pcl,pwl,rfoam,rwat,rglit,mu,np,rm,rp,brdfint)
            !write(*,*) 'brdfint = ',brdfint(1)

         endif

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

         brdf_wvl(1) = 0.47
         brdf_wvl(2) = 0.56
         brdf_wvl(3) = 0.64
         brdf_wvl(4) = 0.86
         brdf_wvl(5) = 1.24
         brdf_wvl(6) = 1.64
         brdf_wvl(7) = 2.13

         brdf_mult_vis(1) = (1.-fsno_read(m,j))*0.4364+fsno_read(m,j)*0.5096
         brdf_mult_vis(2) = (1.-fsno_read(m,j))*0.2366+fsno_read(m,j)*0.2437
         brdf_mult_vis(3) = (1.-fsno_read(m,j))*0.3265+fsno_read(m,j)*0.2465
        
         !from literature
         !brdf_mult_vis(1) = 0.4364
         !brdf_mult_vis(2) = 0.2366
         !brdf_mult_vis(3) = 0.3265 

         !from calculations
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
         brdf_mult_nir(4) =  0.5318
         brdf_mult_nir(5) =  0.2790
         brdf_mult_nir(6) =  0.1539
         brdf_mult_nir(7) =  0.1446

         !brdf_mult_nir(4) =  0.5271
         !brdf_mult_nir(5) =  0.1795
         !brdf_mult_nir(6) =  0.0000
         !brdf_mult_nir(7) =  0.2755

         !brdf_mult_nir(4) = (1.-fsno_read(m,j))*0.5271+fsno_read(m,j)*0.5318
         !brdf_mult_nir(5) = (1.-fsno_read(m,j))*0.1795+fsno_read(m,j)*0.2790
         !brdf_mult_nir(6) = (1.-fsno_read(m,j))*0.0 + fsno_read(m,j)*0.1539
         !brdf_mult_nir(7) = (1.-fsno_read(m,j))*0.2755+fsno_read(m,j)*0.1446

!         if (ice_flag(m,j)) then
!
!            !integrate up ice-brdf
!            brdf_param2(1,1) = 0.76
!            brdf_param2(2,1) = 0.74
!            brdf_param2(3,1) = 0.72
!            brdf_param2(4,1) = 0.6
!            brdf_param2(5,1) = 0.25
!            brdf_param2(6,1) = 0.07
!            brdf_param2(7,1) = 0.06
!
!            do ii=1,7
!               brdf_param2(ii,2) = 0. !volumetric
!               brdf_param2(ii,3) = 0. !geometric
!            end do

!            do ii=1,7
!               spec_alb_ws(ii) = 0.
!               spec_alb_ws(ii) = spec_alb_ws(ii)+gws_iso*brdf_param2(ii,1)
!               spec_alb_ws(ii) = spec_alb_ws(ii)+gws_vol*brdf_param2(ii,2)
!               spec_alb_ws(ii) = spec_alb_ws(ii)+gws_geo*brdf_param2(ii,3)
!
!               spec_alb_bs(ii) = 0.
!               spec_alb_bs(ii) = spec_alb_bs(ii)+ brdf_param2(ii,1)* &
!                    (gbs_iso(1) + gbs_iso(2)*solar_zen_rad**2 + gbs_iso(3)*solar_zen_rad**3)
!               spec_alb_bs(ii) = spec_alb_bs(ii)+ brdf_param2(ii,2)* &
!                    (gbs_vol(1) + gbs_vol(2)*solar_zen_rad**2 + gbs_vol(3)*solar_zen_rad**3)
!               spec_alb_bs(ii) = spec_alb_bs(ii)+ brdf_param2(ii,3)* &
!                    (gbs_geo(1) + gbs_geo(2)*solar_zen_rad**2 + gbs_geo(3)*solar_zen_rad**3)
!            end do
!            modis_asdir(i,j,m) = 0.
!            modis_asdif(i,j,m) = 0.
!            modis_aldir(i,j,m) = 0.
!            modis_aldif(i,j,m) = 0.
!
!            do ii=1,7
!               modis_asdir(i,j,m) = modis_asdir(i,j,m) + spec_alb_bs(ii)*brdf_mult_vis(ii) 
!               modis_asdif(i,j,m) = modis_asdif(i,j,m) + spec_alb_ws(ii)*brdf_mult_vis(ii) 
!               modis_aldir(i,j,m) = modis_aldir(i,j,m) + spec_alb_bs(ii)*brdf_mult_nir(ii) 
!               modis_aldif(i,j,m) = modis_aldif(i,j,m) + spec_alb_ws(ii)*brdf_mult_nir(ii) 
!            end do
!            !write(*,*) 'asdir = ',modis_asdir(i,j,m)
!            !write(*,*) 'asdif = ',modis_asdif(i,j,m)
!            !write(*,*) 'aldir = ',modis_aldir(i,j,m)
!            !write(*,*) 'aldif = ',modis_aldif(i,j,m)
!            !stop
!         endif
      end do 
    end do 
    end do

    !write(*,*) 'max asdir = ',maxval(maxval(modis_asdir(4,:,:),1))
    !write(*,*) 'min asdir = ',minval(minval(modis_asdir(4,:,:),1))
    !write(*,*) 'max asdif = ',maxval(maxval(modis_asdif(4,:,:),1))
    !write(*,*) 'min asdif = ',minval(minval(modis_asdif(4,:,:),1))
    !stop 'error get_surface_properties.F90 line 599'

    deallocate(dummy_array4)
    deallocate(snow_asdir)
    deallocate(snow_asdif)
    deallocate(snow_aldir)
    deallocate(snow_aldif)

    call wrap_close(nc_id)

    ! broadcast brdf info to all nodes

  end subroutine get_brdf

end module get_surface_properties


