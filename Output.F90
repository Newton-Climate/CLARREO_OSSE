module Output

  use shr_kind_mod, only: r8 => shr_kind_r8                                   

  use ppgrid
  use pmgrid
  use error_messages

  implicit none
  save

!
! File ID
!
  integer :: nfid                          ! NetCDF id

!
! Variable IDs
!
  integer :: vid_FLN                       ! Variable id for FLN
  integer :: vid_FLNS                      ! Variable id for FLNS
  integer :: vid_FLNSC                     ! Variable id for FLNSC
  integer :: vid_FLNT                      ! Variable id for FLNT
  integer :: vid_FLNTC                     ! Variable id for FLNTC
  integer :: vid_FLN200                    ! Variable id for FLN200
  integer :: vid_FLN200C                   ! Variable id for FLN200C
  integer :: vid_FLWDS                     ! Variable id for FLWDS (DRF)
  integer :: vid_FSDS                      ! Variable id for FSDS
  integer :: vid_FSDSC                     ! Variable id for FSDSC
  integer :: vid_FSN                       ! Variable id for FSN
  integer :: vid_FSNS                      ! Variable id for FSNS
  integer :: vid_FSNSC                     ! Variable id for FSNSC
  integer :: vid_FSNT                      ! Variable id for FSNT
  integer :: vid_FSNTC                     ! Variable id for FSNTC
  integer :: vid_FSN200                    ! Variable id for FSN200
  integer :: vid_FSN200C                   ! Variable id for FSN200C
  integer :: vid_FSNIRTOA                  ! Variable id for FSNIRTOA
  integer :: vid_FSNIRTOAC                 ! Variable id for FSNIRTOAC
  integer :: vid_QRL                       ! Variable id for QRL
  integer :: vid_QRS                       ! Variable id for QRS
  integer :: vid_SOLIN                     ! Variable id for SOLIN
  integer :: vid_SOLL                      ! Variable id for SOLL
  integer :: vid_SOLLD                     ! Variable id for SOLLD
  integer :: vid_SOLS                      ! Variable id for SOLS
  integer :: vid_SOLSD                     ! Variable id for SOLSD
  integer :: vid_WAVELENGTH_LRES           ! Variable id for WAVELENGTH_LRES
  integer :: vid_RADIANCE_LRES_CLR         ! Variable id for RADIANCE_LRES_CLR
  integer :: vid_RADIANCE_LRES_ALL         ! Variable id for RADIANCE_LRES_ALL
  integer :: vid_WAVELENGTH_HRES           ! Variable id for WAVELENGTH_HRES
  integer :: vid_RADIANCE_HRES_CLR         ! Variable id for RADIANCE_LRES_CLR
  integer :: vid_RADIANCE_HRES_ALL         ! Variable id for RADIANCE_LRES_ALL
  integer :: vid_SOLAR_FLUX                ! Variable id for SOLAR_FLUX
  integer :: vid_DIFFUSE_FLUX_CLR          ! Variable id for DIFFUSE_FLUX_CLR
  integer :: vid_DIFFUSE_FLUX_ALL          ! Variable id for DIFFUSE_FLUX_ALL
  integer :: vid_SOLAR_ZENITH              ! Variable id for SOLAR_ZENITH
  integer :: vid_AOD                       ! Variable id for AOD
  !integer :: vid_AERTAU                       ! Variable id for AERTAU (radcswmx.F90)
  integer :: vid_TAU_SULF                       ! Variable id for TAU_SULF 
  integer :: vid_TAU_DUST                       ! Variable id for TAU_SULF 
  integer :: vid_TAU_SOOT                       ! Variable id for TAU_SULF 
  integer :: vid_TAU_SSLT                       ! Variable id for TAU_SULF 
  integer :: vid_NCONFIG_REAL
  integer :: vid_BB_UPDIFFUSE_CLR                   ! Variable id for BB_UPDIFFUSE_CLR
  integer :: vid_BB_UPDIFFUSE_ALL                   ! Variable id for BB_UPDIFFUSE_ALL
  integer :: vid_BB_DNDIFFUSE_CLR                   ! Variable id for BB_DNDIFFUSE_CLR
  integer :: vid_BB_DNDIFFUSE_ALL                   ! Variable id for BB_DNDIFFUSE_ALL
  integer :: vid_BB_DNDIRECT_CLR                   ! Variable id for BB_DNDIRECT_CLR
  integer :: vid_BB_DNDIRECT_ALL                   ! Variable id for BB_DNDIRECT_ALL
  integer :: vid_CLDTAU_ICE_SW                     ! Variable id for CLDICE_TAU sw
  integer :: vid_CLDTAU_LIQ_SW                     ! Variable id for CLDLIQ_TAU sw
  integer :: vid_CLDTAU_LW                     ! Variable id for CLD_TAU lw
  integer :: vid_VIS_UPDIFFUSE_CLR                   ! Variable id for VIS_UPDIFFUSE_CLR
  integer :: vid_VIS_UPDIFFUSE_ALL                   ! Variable id for VIS_UPDIFFUSE_ALL
  integer :: vid_VIS_DNDIFFUSE_CLR                   ! Variable id for VIS_DNDIFFUSE_CLR
  integer :: vid_VIS_DNDIFFUSE_ALL                   ! Variable id for VIS_DNDIFFUSE_ALL
  integer :: vid_VIS_DNDIRECT_CLR                   ! Variable id for VIS_DNDIRECT_CLR
  integer :: vid_VIS_DNDIRECT_ALL                   ! Variable id for VIS_DNDIRECT_ALL
  integer :: vid_LAT                                ! Variable id for latitude values
  integer :: vid_LON                                ! Variable id for longitude values

  !BEGIN DRF FOR FLUXES
  integer :: vid_NIR_UPDIFFUSE_CLR                   ! Variable id for NIR_UPDIFFUSE_CLR
  integer :: vid_NIR_UPDIFFUSE_ALL                   ! Variable id for NIR_UPDIFFUSE_ALL
  integer :: vid_NIR_DNDIFFUSE_CLR                   ! Variable id for NIR_DNDIFFUSE_CLR
  integer :: vid_NIR_DNDIFFUSE_ALL                   ! Variable id for NIR_DNDIFFUSE_ALL
  integer :: vid_NIR_DNDIRECT_CLR                   ! Variable id for NIR_DNDIRECT_CLR
  integer :: vid_NIR_DNDIRECT_ALL                   ! Variable id for NIR_DNDIRECT_ALL

  !integer :: vid_MODIS_ASDIR                       ! Variable id for MODIS_ASDIR
  !integer :: vid_MODIS_ASDIF                       ! Variable id for MODIS_ASDIF
  !integer :: vid_MODIS_ALDIR                       ! Variable id for MODIS_ALDIR
  !integer :: vid_MODIS_ALDIF                       ! Variable id for MODIS_ALDIF

  !END DRF FOR FLUXES


!
! Dimension ids
!
  integer :: did_lat                       ! Latitude dimension ID
  integer :: did_lon                       ! Longitude dimension ID
  integer :: did_time                      ! Time dimension ID
  integer :: did_lev                       ! Level dimension ID
  integer :: did_ilev                      ! Interface dimension ID
  !DRF
  integer :: did_wvl                       ! Wavelength dimension ID
  integer :: did_wvl_hres                  ! Wavelength dimension ID

CONTAINS

  !subroutine create_output(ipath, opath, nslice, date, datesec, qr_option)
  subroutine create_output(ipath, opath, nslice, qr_option)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Create output file for computed radiation fluxes
!
! Method: 
! Use NetCDF wrapper routines to create file.  
!
! Output dimensions:
!    lat
!    lon
!    time
!    pver (for QRS and QRL)
!
! Output fields:
!    FLN
!    FLNS
!    FLNSC
!    FLNT
!    FLNTC
!    FLN200
!    FLN200C
!    FLWDS
!    FSDS
!    FSDSC
!    FSN
!    FSNS
!    FSNSC
!    FSNT
!    FSNTC
!    FSN200
!    FSN200C
!    FSNIRTOA
!    FSNIRTOAC
!    QRL
!    QRS
!    SOLIN
!    SOLL
!    SOLLD
!    SOLS
!    SOLSD
!
! Output global attributes
!    date
!    datesec
!    input path file
!
! Author: W. Collins
! 
!-----------------------------------------------------------------------

    implicit none
#include <netcdf.inc>
!
! Input arguments
!
    character(len=*), intent(in) :: ipath        ! path to input data
    character(len=*), intent(in) :: opath        ! path to output data
    integer, intent(in)      :: nslice       ! number of time slices 
    !integer, pointer         :: date(:)      ! Date YYYYMMDD
    !integer, pointer         :: datesec(:)   ! Seconds in date
    integer, intent(in)      :: qr_option    ! Option for output of QRL/QRS
                                             !   0 = no QRL/QRS on output
                                             !   1 = QRL/QRS on output

!
! Local variables
!
    integer :: vdims(4)                      ! dimension ids
    integer :: vdims2(4) !DRF
    integer :: vdims3(4) !DRF
    integer :: vdims_single
    integer :: nvdims                        ! number of dimensions
    integer :: nvdims2                        ! number of dimensions for DRF variables
    integer :: ret                           ! return code

    character (len = 13) :: routine = "create_output"

!----------------------------------------------------------------------

    call wrap_create(opath, nf_clobber, nfid)
!
! Define dimensions
!    
    call wrap_def_dim(nfid, "lat", plat, did_lat)
    call wrap_def_dim(nfid, "lon", plon, did_lon)
    call wrap_def_dim(nfid, "time", nf_unlimited, did_time)
    if (qr_option == 1) then
       call wrap_def_dim(nfid, "lev",  plev,  did_lev)
       call wrap_def_dim(nfid, "ilev", plevp, did_ilev)
    endif
    call wrap_def_dim(nfid,"wavelength",wvlng,did_wvl)
    call wrap_def_dim(nfid,"wavelength_hres",wvlng_hres,did_wvl_hres)
!
! Define variables
!
    nvdims = 3
    vdims = (/did_lon, did_lat, did_time, -1/)
    vdims2 = (/did_lon, did_lat, did_wvl, did_time/)
    vdims3 = (/did_lon, did_lat, did_wvl_hres, did_time/)

    call wrap_def_var(nfid, "FLNS",    nf_double, nvdims, vdims, vid_FLNS)
    call wrap_def_var(nfid, "FLNSC",   nf_double, nvdims, vdims, vid_FLNSC) 
    call wrap_def_var(nfid, "FLNT",    nf_double, nvdims, vdims, vid_FLNT) 
    call wrap_def_var(nfid, "FLNTC",   nf_double, nvdims, vdims, vid_FLNTC) 
    call wrap_def_var(nfid, "FLN200",  nf_double, nvdims, vdims, vid_FLN200) 
    call wrap_def_var(nfid, "FLN200C", nf_double, nvdims, vdims, vid_FLN200C) 
    call wrap_def_var(nfid, "FLWDS",   nf_double, nvdims, vdims, vid_FLWDS) !DRF
    call wrap_def_var(nfid, "FSDS",    nf_double, nvdims, vdims, vid_FSDS) 
    call wrap_def_var(nfid, "FSDSC",   nf_double, nvdims, vdims, vid_FSDSC) 
    call wrap_def_var(nfid, "FSNS",    nf_double, nvdims, vdims, vid_FSNS) 
    call wrap_def_var(nfid, "FSNSC",   nf_double, nvdims, vdims, vid_FSNSC) 
    call wrap_def_var(nfid, "FSNT",    nf_double, nvdims, vdims, vid_FSNT) 
    call wrap_def_var(nfid, "FSNTC",   nf_double, nvdims, vdims, vid_FSNTC) 
    call wrap_def_var(nfid, "FSN200",  nf_double, nvdims, vdims, vid_FSN200) 
    call wrap_def_var(nfid, "FSN200C", nf_double, nvdims, vdims, vid_FSN200C) 
    call wrap_def_var(nfid, "FSNIRTOA", nf_double, nvdims, vdims, vid_FSNIRTOA) 
    call wrap_def_var(nfid, "FSNIRTOAC", nf_double, nvdims, vdims, vid_FSNIRTOAC) 
    call wrap_def_var(nfid, "SOLIN",   nf_double, nvdims, vdims, vid_SOLIN) 
    call wrap_def_var(nfid, "SOLL",    nf_double, nvdims, vdims, vid_SOLL) 
    call wrap_def_var(nfid, "SOLLD",   nf_double, nvdims, vdims, vid_SOLLD) 
    call wrap_def_var(nfid, "SOLS",    nf_double, nvdims, vdims, vid_SOLS) 
    call wrap_def_var(nfid, "SOLSD",   nf_double, nvdims, vdims, vid_SOLSD) 
    !call wrap_def_var(nfid, "MODIS_ASDIR",   nf_double, nvdims, vdims, vid_MODIS_ASDIR) 
    !call wrap_def_var(nfid, "MODIS_ASDIF",   nf_double, nvdims, vdims, vid_MODIS_ASDIF) 
    !call wrap_def_var(nfid, "MODIS_ALDIR",   nf_double, nvdims, vdims, vid_MODIS_ALDIR) 
    !call wrap_def_var(nfid, "MODIS_ALDIF",   nf_double, nvdims, vdims, vid_MODIS_ALDIF) 


    !DRF
    nvdims2 = 1
    vdims_single = did_wvl
    call wrap_def_var(nfid, "WAVELENGTH_LRES",nf_real, nvdims2, vdims_single, vid_WAVELENGTH_LRES)
    vdims_single = did_wvl_hres
    call wrap_def_var(nfid, "WAVELENGTH_HRES",nf_real, nvdims2, vdims_single, vid_WAVELENGTH_HRES)

    vdims_single = did_lat
    call wrap_def_var(nfid, "lat",nf_double, nvdims2, vdims_single, vid_LAT)
    vdims_single = did_lon
    call wrap_def_var(nfid, "lon",nf_double, nvdims2, vdims_single, vid_LON)

    nvdims2 = 4
    call wrap_def_var(nfid, "RADIANCE_LRES_CLR",nf_real, nvdims2, vdims2, vid_RADIANCE_LRES_CLR)
    call wrap_def_var(nfid, "RADIANCE_LRES_ALL",nf_real, nvdims2, vdims2, vid_RADIANCE_LRES_ALL)
    call wrap_def_var(nfid, "RADIANCE_HRES_CLR",nf_real, nvdims2, vdims3, vid_RADIANCE_HRES_CLR)
    call wrap_def_var(nfid, "RADIANCE_HRES_ALL",nf_real, nvdims2, vdims3, vid_RADIANCE_HRES_ALL)
    call wrap_def_var(nfid, "SOLAR_FLUX",nf_real, nvdims2, vdims2, vid_SOLAR_FLUX)
    call wrap_def_var(nfid, "DIFFUSE_FLUX_CLR",nf_real, nvdims2, vdims2, vid_DIFFUSE_FLUX_CLR)
    call wrap_def_var(nfid, "DIFFUSE_FLUX_ALL",nf_real, nvdims2, vdims2, vid_DIFFUSE_FLUX_ALL)
    call wrap_def_var(nfid, "SOLAR_ZENITH",nf_real, nvdims, vdims, vid_SOLAR_ZENITH)
    call wrap_def_var(nfid, "AOD",nf_real, nvdims, vdims, vid_AOD)
    !call wrap_def_var(nfid, "AERTAU",nf_double, nvdims, vdims, vid_AERTAU)
    call wrap_def_var(nfid, "TAU_SULF",nf_real, nvdims, vdims, vid_TAU_SULF)
    call wrap_def_var(nfid, "TAU_DUST",nf_real, nvdims, vdims, vid_TAU_DUST)
    call wrap_def_var(nfid, "TAU_SOOT",nf_real, nvdims, vdims, vid_TAU_SOOT)
    call wrap_def_var(nfid, "TAU_SSLT",nf_real, nvdims, vdims, vid_TAU_SSLT)
    call wrap_def_var(nfid, "NCONFIG_REAL",nf_real, nvdims, vdims, vid_NCONFIG_REAL)
    call wrap_def_var(nfid, "CLDTAU_ICE_SW",nf_double, nvdims, vdims, vid_CLDTAU_ICE_SW)
    call wrap_def_var(nfid, "CLDTAU_LIQ_SW",nf_double, nvdims, vdims, vid_CLDTAU_LIQ_SW)
    call wrap_def_var(nfid, "CLDTAU_LW",nf_double, nvdims, vdims, vid_CLDTAU_LW)
   
    nvdims = 4
    vdims = (/did_lon, did_lat, did_lev, did_time/)
    call wrap_def_var(nfid, "VIS_UPDIFFUSE_CLR",nf_real, nvdims, vdims, vid_VIS_UPDIFFUSE_CLR)
    call wrap_def_var(nfid, "VIS_UPDIFFUSE_ALL",nf_real, nvdims, vdims, vid_VIS_UPDIFFUSE_ALL)
    call wrap_def_var(nfid, "VIS_DNDIFFUSE_CLR",nf_real, nvdims, vdims, vid_VIS_DNDIFFUSE_CLR)
    call wrap_def_var(nfid, "VIS_DNDIFFUSE_ALL",nf_real, nvdims, vdims, vid_VIS_DNDIFFUSE_ALL)
    call wrap_def_var(nfid, "VIS_DNDIRECT_CLR",nf_real, nvdims, vdims, vid_VIS_DNDIRECT_CLR)
    call wrap_def_var(nfid, "VIS_DNDIRECT_ALL",nf_real, nvdims, vdims, vid_VIS_DNDIRECT_ALL)

    call wrap_def_var(nfid, "NIR_UPDIFFUSE_CLR",nf_real, nvdims, vdims, vid_NIR_UPDIFFUSE_CLR)
    call wrap_def_var(nfid, "NIR_UPDIFFUSE_ALL",nf_real, nvdims, vdims, vid_NIR_UPDIFFUSE_ALL)
    call wrap_def_var(nfid, "NIR_DNDIFFUSE_CLR",nf_real, nvdims, vdims, vid_NIR_DNDIFFUSE_CLR)
    call wrap_def_var(nfid, "NIR_DNDIFFUSE_ALL",nf_real, nvdims, vdims, vid_NIR_DNDIFFUSE_ALL)
    call wrap_def_var(nfid, "NIR_DNDIRECT_CLR",nf_real, nvdims, vdims, vid_NIR_DNDIRECT_CLR)
    call wrap_def_var(nfid, "NIR_DNDIRECT_ALL",nf_real, nvdims, vdims, vid_NIR_DNDIRECT_ALL)

    call wrap_def_var(nfid, "BB_UPDIFFUSE_CLR",nf_real, nvdims, vdims, vid_BB_UPDIFFUSE_CLR)
    call wrap_def_var(nfid, "BB_UPDIFFUSE_ALL",nf_real, nvdims, vdims, vid_BB_UPDIFFUSE_ALL)
    call wrap_def_var(nfid, "BB_DNDIFFUSE_CLR",nf_real, nvdims, vdims, vid_BB_DNDIFFUSE_CLR)
    call wrap_def_var(nfid, "BB_DNDIFFUSE_ALL",nf_real, nvdims, vdims, vid_BB_DNDIFFUSE_ALL)
    call wrap_def_var(nfid, "BB_DNDIRECT_CLR",nf_real, nvdims, vdims, vid_BB_DNDIRECT_CLR)
    call wrap_def_var(nfid, "BB_DNDIRECT_ALL",nf_real, nvdims, vdims, vid_BB_DNDIRECT_ALL)

    !DRF 
    if (qr_option == 1) then
       nvdims = 4
       vdims = (/did_lon, did_lat, did_lev, did_time/)
       
       call wrap_def_var(nfid, "QRL", nf_double, nvdims, vdims, vid_QRL) 
       call wrap_def_var(nfid, "QRS", nf_double, nvdims, vdims, vid_QRS) 

       vdims = (/did_lon, did_lat, did_ilev, did_time/)
       
       call wrap_def_var(nfid, "FLN", nf_double, nvdims, vdims, vid_FLN) 
       call wrap_def_var(nfid, "FSN", nf_double, nvdims, vdims, vid_FSN)   
    endif
!
! Add attributes
!
    call wrap_put_att_text(nfid, nf_global, "input_history_file", ipath)

    !ret = nf_put_att_int (nfid, nf_global, "date", nf_int, nslice, date)
    !if (ret/=NF_NOERR) call handle_error (ret)

    !ret = nf_put_att_int (nfid, nf_global, "datesec", nf_int, nslice, datesec)
    !if (ret/=NF_NOERR) call handle_error (ret)

    call wrap_close(nfid)

    return

  end subroutine create_output

  subroutine write_output(opath, itime, FLN, FLNS, FLNSC, FLNT, FLNTC, &
                          FLN200, FLN200C,FLWDS, &
                          FSDS, FSDSC, FSN, FSNS, FSNSC, FSNT, FSNTC, &
                          FSN200, FSN200C, FSNIRTOA,FSNIRTOAC, &
                          QRL, QRS, SOLIN, SOLL, SOLLD, SOLS, SOLSD, IN_LATITUDE,IN_LONGITUDE,&
                          WAVELENGTH_LRES, RADIANCE_LRES_CLR, RADIANCE_LRES_ALL,WAVELENGTH_HRES, &
                          RADIANCE_HRES_CLR,RADIANCE_HRES_ALL,SOLAR_FLUX,DIFFUSE_FLUX_CLR,DIFFUSE_FLUX_ALL, &
                          SOLAR_ZENITH,AOD,TAU_SULF,TAU_DUST,TAU_SOOT,TAU_SSLT,NCONFIG_REAL, &
                          CLDTAU_ICE_SW,CLDTAU_LIQ_SW,CLDTAU_LW, &
                          BB_UPDIFFUSE_CLR,BB_UPDIFFUSE_ALL, &
                          BB_DNDIFFUSE_CLR,BB_DNDIFFUSE_ALL, &
                          BB_DNDIRECT_CLR,BB_DNDIRECT_ALL, &
                          FLX_UPDIFFUSE_CLR,FLX_UPDIFFUSE_ALL, &
                          FLX_DNDIFFUSE_CLR,FLX_DNDIFFUSE_ALL, &
                          FLX_DNDIRECT_CLR,FLX_DNDIRECT_ALL,qr_option)


!----------------------------------------------------------------------- 
! 
! Purpose: 
! Write output file for computed radiation fluxes
!
! Method: 
! Use NetCDF wrapper routines to write file.  
!
! Output fields:
!    FLN
!    FLNS
!    FLNSC
!    FLNT
!    FLNTC
!    FLN200
!    FLN200C
!    FLWDS
!    FSDS
!    FSDSC
!    FSN
!    FSNS
!    FSNSC
!    FSNT
!    FSNTC
!    FSN200
!    FSN200C
!    FSNIRTOA
!    FSNIRTOAC
!    QRL
!    QRS
!    SOLIN
!    SOLL
!    SOLLD
!    SOLS
!    SOLSD
!    MODIS_ASDIR
!    MODIS_ASDIF
!    MODIS_ALDIR
!    MODIS_ALDIF
!    RADIANCE
!
! Author: W. Collins
!         modified by D. Feldman 
!-----------------------------------------------------------------------

    implicit none
#include <netcdf.inc>
!
! Input arguments
!
    character(len=*), intent(in) :: opath        ! path to output data
    integer, intent(in)      :: itime        ! time slice

    real(r8), intent(in) :: fln(plon,pverp,plat)! Net longwave flux at each interface
    real(r8), intent(in) :: flns(plon,plat)    ! Srf longwave cooling (up-down) flux
    real(r8), intent(in) :: flnsc(plon,plat)   ! Net surface clear-sky longwave
    real(r8), intent(in) :: flnt(plon,plat)    ! Net outgoing lw flux at model top
    real(r8), intent(in) :: flntc(plon,plat)   ! Net TOA clear-sky longwave
    real(r8), intent(in) :: fln200(plon,plat)  ! Net all-sky longwave at 200mb
    real(r8), intent(in) :: fln200c(plon,plat) ! Net clear-sky longwave at 200mb
    real(r8), intent(in) :: flwds(plon,plat)   ! Downwelling all-sky longwave flux at surface (DRF)
    real(r8), intent(in) :: fsds(plon,plat)    ! Downwelling surface shortwave
    real(r8), intent(in) :: fsdsc(plon,plat)   ! Downwelling clear surface shortwave
    real(r8), intent(in) :: fsn(plon,pverp,plat)! Net shortwave flux at each interface
    real(r8), intent(in) :: fsns(plon,plat)    ! Surface solar absorbed flux
    real(r8), intent(in) :: fsnsc(plon,plat)   ! Net surface clear shortwave
    real(r8), intent(in) :: fsnt(plon,plat)    ! Net column abs solar flux at model top
    real(r8), intent(in) :: fsntc(plon,plat)   ! Net TOA cear shortwave
    real(r8), intent(in) :: fsn200(plon,plat)  ! Net all-sky shortwave at 200mb
    real(r8), intent(in) :: fsn200c(plon,plat) ! Net clear-sky shortwave at 200mb
    real(r8), intent(in) :: fsnirtoa(plon,plat) ! Net all-sky nir flux at toa
    real(r8), intent(in) :: fsnirtoac(plon,plat) ! Net clear-sky nir flux at toa
    real(r8), intent(in) :: qrl(plon,pver,plat)! Longwave  heating rate
    real(r8), intent(in) :: qrs(plon,pver,plat)! Shortwave heating rate
    real(r8), intent(in) :: solin(plon,plat)   ! TOA insolation
    
    real(r8), intent(in) :: sols(plon,plat)    ! Direct beam solar rad. onto srf (sw)
    real(r8), intent(in) :: soll(plon,plat)    ! Direct beam solar rad. onto srf (lw)
    real(r8), intent(in) :: solsd(plon,plat)   ! Diffuse solar radiation onto srf (sw)
    real(r8), intent(in) :: solld(plon,plat)   ! Diffuse solar radiation onto srf (lw)
    
    !real(r8), intent(in) :: modis_asdir(plon,plat) ! modis direct albedo data being fed into cam (sw)
    !real(r8), intent(in) :: modis_asdif(plon,plat) ! modis diffuse albedo data being fed into cam (sw)
    !real(r8), intent(in) :: modis_aldir(plon,plat) ! modis direct albedo data being fed into cam (lw)
    !real(r8), intent(in) :: modis_aldif(plon,plat) ! modis diffuse albedo data being fed into cam (lw)


    !real(r8), intent(in) :: aertau(plon,aertau_species,plat)  ! aerosol optical depths calculated by radcswmx.F90 (DRF)
    real*4, intent(in) :: wavelength_lres(wvlng) ! DRF sw wavelength values from modtran
    real*4, intent(in) :: radiance_lres_clr(plon,wvlng,plat) ! DRF sw radiance values from modtran clear-sky
    real*4, intent(in) :: radiance_lres_all(plon,wvlng,plat) ! DRF sw radiance values from modtran all-sky
    real*4, intent(in) :: wavelength_hres(wvlng_hres) ! DRF sw wavelength values from modtran
    real*4, intent(in) :: radiance_hres_clr(plon,wvlng_hres,plat) ! DRF sw radiance values from modtran clear-sky
    real*4, intent(in) :: radiance_hres_all(plon,wvlng_hres,plat) ! DRF sw radiance values from modtran all-sky
    real*4, intent(in) :: solar_flux(plon,wvlng,plat) ! DRF sw TOA flux values from modtran
    real*4, intent(in) :: diffuse_flux_clr(plon,wvlng,plat) ! DRF sw radiance values from modtran clear-sky
    real*4, intent(in) :: diffuse_flux_all(plon,wvlng,plat) ! DRF sw radiance values from modtran all-sky
    real*4, intent(in) :: solar_zenith(plon,plat) ! DRF solar zenith angle from modtran 
    real*4, intent(in) :: aod(plon,plat) ! DRF solar zenith angle from modtran 
    real*4, intent(in) :: tau_sulf(plon,plat) ! DRF solar zenith angle from modtran 
    real*4, intent(in) :: tau_dust(plon,plat) ! DRF solar zenith angle from modtran 
    real*4, intent(in) :: tau_soot(plon,plat) ! DRF solar zenith angle from modtran 
    real*4, intent(in) :: tau_sslt(plon,plat) ! DRF solar zenith angle from modtran 
    real*4, intent(in) :: nconfig_real(plon,plat) ! # cloud configs
    real(r8), intent(in) :: cldtau_ice_sw(plon,plat) ! DRF solar zenith angle from modtran 
    real(r8), intent(in) :: cldtau_liq_sw(plon,plat) ! DRF solar zenith angle from modtran 
    real(r8), intent(in) :: cldtau_lw(plon,plat) ! DRF solar zenith angle from modtran 
    real(r8), intent(in) :: in_latitude(plat)     ! DRF latitude coordinate values
    real(r8), intent(in) :: in_longitude(plon)     ! DRF longitude coordinate values
    real*4, intent(in) :: bb_updiffuse_clr(plon,pver,plat) ! DRF upwelling diffuse broad-band flux from MODTRAN
    real*4, intent(in) :: bb_updiffuse_all(plon,pver,plat) ! DRF upwelling diffuse broad-band flux from MODTRAN
    real*4, intent(in) :: bb_dndiffuse_clr(plon,pver,plat) ! DRF downwelling diffuse broad-band flux from MODTRAN
    real*4, intent(in) :: bb_dndiffuse_all(plon,pver,plat) ! DRF downwelling diffuse broad-band flux from MODTRAN
    real*4, intent(in) :: bb_dndirect_clr(plon,pver,plat) ! DRF downwelling direct broad-band flux from MODTRAN
    real*4, intent(in) :: bb_dndirect_all(plon,pver,plat) ! DRF downwelling direct broad-band flux from MODTRAN

    real*4, intent(in) :: flx_updiffuse_clr(plon,pver,num_bands,plat) ! DRF upwelling diffuse broad-band flux from MODTRAN
    real*4, intent(in) :: flx_updiffuse_all(plon,pver,num_bands,plat) ! DRF upwelling diffuse broad-band flux from MODTRAN
    real*4, intent(in) :: flx_dndiffuse_clr(plon,pver,num_bands,plat) ! DRF downwelling diffuse broad-band flux from MODTRAN
    real*4, intent(in) :: flx_dndiffuse_all(plon,pver,num_bands,plat) ! DRF downwelling diffuse broad-band flux from MODTRAN
    real*4, intent(in) :: flx_dndirect_clr(plon,pver,num_bands,plat) ! DRF downwelling direct broad-band flux from MODTRAN
    real*4, intent(in) :: flx_dndirect_all(plon,pver,num_bands,plat) ! DRF downwelling direct broad-band flux from MODTRAN
    
    integer, intent(in)      :: qr_option    ! Option for output of QRL/QRS
                                             !   0 = no QRL/QRS on output
                                             !   1 = QRL/QRS on output

!
! Local variables
!
    integer :: start(4)                 ! starting point on edges
    integer :: count(4)                 ! number of elements to write

    !DRF
    integer :: start2(4)
    integer :: count2(4)

    integer :: start_single
    integer :: count_single

    real(r8) :: buffer(plon,plat,pver)  ! Buffer for transposing QRL/QRS output
    real(r8) :: bufferp(plon,plat,pverp)! Buffer for transposing QRL/QRS output
    real*4, allocatable :: bufferr_modtran(:,:,:) !Buffer for handling all of the Modtran output

    !real*4 :: bufferr(plon,plat,wvlng)! Buffer for transposing RADIANCE output
    !real*4 :: bufferr_hres(plon,plat,wvlng_hres)! Buffer for transposing RADIANCE_HRES output
    !real*4 :: bufferr_lev(plon,plat,pver) ! Buffer for transposing fluxes
    integer :: ilon                     ! Longitude index
    character (len = 12) :: routine = "write_output"

!----------------------------------------------------------------------

    call wrap_open(opath, nf_write, nfid)
!
! Write variables
!
    start = (/1,    1,    itime, -1/)
    count = (/plon, plat, 1,     -1/)
    call wrap_put_vara_realx(nfid, vid_FLNS,    start, count, FLNS)
    call wrap_put_vara_realx(nfid, vid_FLNSC,   start, count, FLNSC) 
    call wrap_put_vara_realx(nfid, vid_FLNT,    start, count, FLNT) 
    call wrap_put_vara_realx(nfid, vid_FLNTC,   start, count, FLNTC) 
    call wrap_put_vara_realx(nfid, vid_FLN200,  start, count, FLN200) 
    call wrap_put_vara_realx(nfid, vid_FLN200C, start, count, FLN200C) 
    call wrap_put_vara_realx(nfid, vid_FLWDS,   start, count, FLWDS)
    call wrap_put_vara_realx(nfid, vid_FSDS,    start, count, FSDS) 
    call wrap_put_vara_realx(nfid, vid_FSDSC,   start, count, FSDSC) 
    call wrap_put_vara_realx(nfid, vid_FSNS,    start, count, FSNS) 
    call wrap_put_vara_realx(nfid, vid_FSNSC,   start, count, FSNSC) 
    call wrap_put_vara_realx(nfid, vid_FSNT,    start, count, FSNT) 
    call wrap_put_vara_realx(nfid, vid_FSNTC,   start, count, FSNTC) 
    call wrap_put_vara_realx(nfid, vid_FSN200,  start, count, FSN200) 
    call wrap_put_vara_realx(nfid, vid_FSN200C, start, count, FSN200C) 
    call wrap_put_vara_realx(nfid, vid_FSNIRTOA, start, count, FSNIRTOA) 
    call wrap_put_vara_realx(nfid, vid_FSNIRTOAC, start, count, FSNIRTOAC) 
    call wrap_put_vara_realx(nfid, vid_SOLIN,   start, count, SOLIN) 
    call wrap_put_vara_realx(nfid, vid_SOLL,    start, count, SOLL) 
    call wrap_put_vara_realx(nfid, vid_SOLLD,   start, count, SOLLD) 
    call wrap_put_vara_realx(nfid, vid_SOLS,    start, count, SOLS) 
    call wrap_put_vara_realx(nfid, vid_SOLSD,   start, count, SOLSD) 
    !call wrap_put_vara_realx(nfid, vid_MODIS_ASDIR,   start, count, MODIS_ASDIR) 
    !call wrap_put_vara_realx(nfid, vid_MODIS_ASDIF,   start, count, MODIS_ASDIF) 
    !call wrap_put_vara_realx(nfid, vid_MODIS_ALDIR,   start, count, MODIS_ALDIR) 
    !call wrap_put_vara_realx(nfid, vid_MODIS_ALDIF,   start, count, MODIS_ALDIF) 

    call wrap_put_vara_real(nfid, vid_SOLAR_ZENITH,   start, count, SOLAR_ZENITH) 
    call wrap_put_vara_real(nfid, vid_AOD,   start, count, AOD) 
    call wrap_put_vara_real(nfid, vid_TAU_SULF,   start, count, TAU_SULF) 
    call wrap_put_vara_real(nfid, vid_TAU_DUST,   start, count, TAU_DUST) 
    call wrap_put_vara_real(nfid, vid_TAU_SOOT,   start, count, TAU_SOOT) 
    call wrap_put_vara_real(nfid, vid_TAU_SSLT,   start, count, TAU_SSLT) 
    call wrap_put_vara_real(nfid, vid_NCONFIG_REAL,   start, count, NCONFIG_REAL)
    call wrap_put_vara_realx(nfid,vid_CLDTAU_ICE_SW, start, count, CLDTAU_ICE_SW) 
    call wrap_put_vara_realx(nfid,vid_CLDTAU_LIQ_SW, start, count, CLDTAU_LIQ_SW) 
    call wrap_put_vara_realx(nfid,vid_CLDTAU_LW, start, count, CLDTAU_LW) 
    
    start2 = (/1,1, 1, itime/)
    count2 = (/plon, plat, wvlng,1/)
    
    allocate(bufferr_modtran(plon,plat,wvlng))
    do ilon = 1, plon
      bufferr_modtran(ilon,:,:) = transpose(radiance_lres_clr(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_RADIANCE_LRES_CLR,start2, count2, bufferr_modtran)  !DRF
    deallocate(bufferr_modtran)
 
    start_single = 1
    count_single = wvlng
    call wrap_put_vara_real(nfid, vid_WAVELENGTH_LRES,start_single, count_single, WAVELENGTH_LRES) !DRF
    count_single = wvlng_hres
    call wrap_put_vara_real(nfid, vid_WAVELENGTH_HRES,start_single, count_single, WAVELENGTH_HRES) !DRF

    !write(*,*) 'in_latitude = ',in_latitude
    !write(*,*) 'in_longitude= ',in_longitude

    count_single = plat
    call wrap_put_vara_realx(nfid, vid_LAT,start_single, count_single, IN_LATITUDE) !DRF
    count_single = plon 
    call wrap_put_vara_realx(nfid, vid_LON,start_single, count_single, IN_LONGITUDE) !DRF
    
    allocate(bufferr_modtran(plon,plat,wvlng))
    do ilon = 1, plon
      bufferr_modtran(ilon,:,:) = transpose(radiance_lres_all(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_RADIANCE_LRES_ALL,start2, count2, bufferr_modtran)  !DRF
    deallocate(bufferr_modtran)

    allocate(bufferr_modtran(plon,plat,wvlng))
    do ilon = 1, plon
      bufferr_modtran(ilon,:,:) = transpose(diffuse_flux_clr(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_DIFFUSE_FLUX_CLR,start2, count2, bufferr_modtran)  !DRF
    deallocate(bufferr_modtran)

    allocate(bufferr_modtran(plon,plat,wvlng))
    do ilon = 1, plon
      bufferr_modtran(ilon,:,:) = transpose(diffuse_flux_all(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_DIFFUSE_FLUX_ALL,start2, count2, bufferr_modtran)  !DRF
    deallocate(bufferr_modtran)

    !DRF
    allocate(bufferr_modtran(plon,plat,wvlng))
    do ilon = 1, plon
      bufferr_modtran(ilon,:,:) = transpose(solar_flux(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_SOLAR_FLUX,start2, count2, bufferr_modtran)
    deallocate(bufferr_modtran)
    !DRF

    allocate(bufferr_modtran(plon,plat,wvlng_hres))
    count2 = (/plon, plat, wvlng_hres,1/)
    do ilon = 1, plon
      bufferr_modtran(ilon,:,:) = transpose(radiance_hres_clr(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_RADIANCE_HRES_CLR,start2, count2, bufferr_modtran)
    deallocate(bufferr_modtran)
     
    allocate(bufferr_modtran(plon,plat,wvlng_hres))
    do ilon = 1, plon
      bufferr_modtran(ilon,:,:) = transpose(radiance_hres_all(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_RADIANCE_HRES_ALL,start2, count2, bufferr_modtran)
    deallocate(bufferr_modtran)

    !write(*,*) 'updiffuse_flux = ',flx_updiffuse_clr(ilon,26,1,50)
    allocate(bufferr_modtran(plon,plat,pver))
    start = (/1,    1,    1,   itime/)
    count = (/plon, plat, pver ,1    /)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(flx_updiffuse_clr(ilon,:,1,:))
    end do
    call wrap_put_vara_real(nfid, vid_VIS_UPDIFFUSE_CLR,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(flx_updiffuse_all(ilon,:,1,:))
    end do
    call wrap_put_vara_real(nfid, vid_VIS_UPDIFFUSE_ALL,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(flx_dndiffuse_clr(ilon,:,1,:))
    end do
    call wrap_put_vara_real(nfid, vid_VIS_DNDIFFUSE_CLR,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(flx_dndiffuse_all(ilon,:,1,:))
    end do
    call wrap_put_vara_real(nfid, vid_VIS_DNDIFFUSE_ALL,start, count, bufferr_modtran)
    !! ok to here
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(flx_dndirect_clr(ilon,:,1,:))
    end do
    call wrap_put_vara_real(nfid, vid_VIS_DNDIRECT_CLR,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(flx_dndirect_all(ilon,:,1,:))
    end do
    call wrap_put_vara_real(nfid, vid_VIS_DNDIRECT_ALL,start, count, bufferr_modtran)

    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(flx_updiffuse_clr(ilon,:,2,:))
    end do
    call wrap_put_vara_real(nfid, vid_NIR_UPDIFFUSE_CLR,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(flx_updiffuse_all(ilon,:,2,:))
    end do
    call wrap_put_vara_real(nfid, vid_NIR_UPDIFFUSE_ALL,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(flx_dndiffuse_clr(ilon,:,2,:))
    end do
    call wrap_put_vara_real(nfid, vid_NIR_DNDIFFUSE_CLR,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(flx_dndiffuse_all(ilon,:,2,:))
    end do
    call wrap_put_vara_real(nfid, vid_NIR_DNDIFFUSE_ALL,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(flx_dndirect_clr(ilon,:,2,:))
    end do
    call wrap_put_vara_real(nfid, vid_NIR_DNDIRECT_CLR,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(flx_dndirect_all(ilon,:,2,:))
    end do
    call wrap_put_vara_real(nfid, vid_NIR_DNDIRECT_ALL,start, count, bufferr_modtran)

    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(bb_updiffuse_clr(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_BB_UPDIFFUSE_CLR,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(bb_updiffuse_all(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_BB_UPDIFFUSE_ALL,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(bb_dndiffuse_clr(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_BB_DNDIFFUSE_CLR,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(bb_dndiffuse_all(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_BB_DNDIFFUSE_ALL,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(bb_dndirect_clr(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_BB_DNDIRECT_CLR,start, count, bufferr_modtran)
    do ilon =1, plon
      bufferr_modtran(ilon,:,:) = transpose(bb_dndirect_all(ilon,:,:))
    end do
    call wrap_put_vara_real(nfid, vid_BB_DNDIRECT_ALL,start, count, bufferr_modtran)
    deallocate(bufferr_modtran)

    if (qr_option == 1) then
       start = (/1,    1,    1,   itime/)
       count = (/plon, plat, plev,1    /)
       
       do ilon = 1, plon
          buffer(ilon,:,:) = transpose(qrl(ilon,:,:))
       end do
       call wrap_put_vara_realx(nfid, vid_QRL, start, count, buffer) 
       do ilon = 1, plon
          buffer(ilon,:,:) = transpose(qrs(ilon,:,:))
       end do
       call wrap_put_vara_realx(nfid, vid_QRS, start, count, buffer) 

       start = (/1,    1,    1,   itime/)
       count = (/plon, plat, plevp,1    /)
       
       do ilon = 1, plon
          bufferp(ilon,:,:) = transpose(fln(ilon,:,:))
       end do
       call wrap_put_vara_realx(nfid, vid_FLN, start, count, bufferp) 
       do ilon = 1, plon
          bufferp(ilon,:,:) = transpose(fsn(ilon,:,:))
       end do
       call wrap_put_vara_realx(nfid, vid_FSN, start, count, bufferp) 


    endif

    call wrap_close(nfid)

    return

  end subroutine write_output

  subroutine get_nc_varids(opath)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! If output file already exists, put in 
!
! Method: 
! Use NetCDF wrapper routines to get var ids 
!
!

     use ioFileMod, only: getfil
    implicit none
#include <netcdf.inc>


!
! Input arguments
!
    character(len=*), intent(in) :: opath        ! path to output data

!Local variables

      integer :: nc_id
      character(len=256) :: locfn !local file

!
! Variable IDs
!

         call getfil(opath, locfn, 0)
         call wrap_open(locfn, 0, nc_id)

         !Get ID for the regular fields
         call wrap_inq_varid(nc_id, "FLN", vid_FLN)
         call wrap_inq_varid(nc_id, "FLNS", vid_FLNS)
         call wrap_inq_varid(nc_id, "FLNSC", vid_FLNSC)
         call wrap_inq_varid(nc_id, "FLNT", vid_FLNT)
         call wrap_inq_varid(nc_id, "FLNTC", vid_FLNTC)
         call wrap_inq_varid(nc_id, "FLN200", vid_FLN200)
         call wrap_inq_varid(nc_id, "FLN200C", vid_FLN200C)
         call wrap_inq_varid(nc_id, "FLWDS", vid_FLWDS)
         call wrap_inq_varid(nc_id, "FSDS", vid_FSDS)
         call wrap_inq_varid(nc_id, "FSDSC", vid_FSDSC)
         call wrap_inq_varid(nc_id, "FSN", vid_FSN)
         call wrap_inq_varid(nc_id, "FSNS", vid_FSNS)
         call wrap_inq_varid(nc_id, "FSNSC", vid_FSNSC)
         call wrap_inq_varid(nc_id, "FSNT", vid_FSNT)
         call wrap_inq_varid(nc_id, "FSNTC", vid_FSNTC)
         call wrap_inq_varid(nc_id, "FSN200", vid_FSN200)
         call wrap_inq_varid(nc_id, "FSN200C", vid_FSN200C)
         call wrap_inq_varid(nc_id, "FSNIRTOA", vid_FSNIRTOA)
         call wrap_inq_varid(nc_id, "FSNIRTOAC", vid_FSNIRTOAC)
         call wrap_inq_varid(nc_id, "QRS", vid_QRS)
         call wrap_inq_varid(nc_id, "QRL", vid_QRL)
         call wrap_inq_varid(nc_id, "SOLIN", vid_SOLIN)
         call wrap_inq_varid(nc_id, "SOLL", vid_SOLL)
         call wrap_inq_varid(nc_id, "SOLS", vid_SOLS)
         call wrap_inq_varid(nc_id, "SOLSD", vid_SOLSD)
         call wrap_inq_varid(nc_id, "SOLLD", vid_SOLLD)
 
         !call wrap_inq_varid(nc_id, "MODIS_ASDIR", vid_MODIS_ASDIR)
         !call wrap_inq_varid(nc_id, "MODIS_ASDIF", vid_MODIS_ASDIF)
         !call wrap_inq_varid(nc_id, "MODIS_ALDIR", vid_MODIS_ALDIR)
         !call wrap_inq_varid(nc_id, "MODIS_ALDIF", vid_MODIS_ALDIF)

         !Get ID numbers of the Modtran added data fields
         call wrap_inq_varid(nc_id, "WAVELENGTH_LRES", vid_WAVELENGTH_LRES)
         call wrap_inq_varid(nc_id, "RADIANCE_LRES_CLR", vid_RADIANCE_LRES_CLR)
         call wrap_inq_varid(nc_id, "RADIANCE_LRES_ALL", vid_RADIANCE_LRES_ALL)
         call wrap_inq_varid(nc_id, "WAVELENGTH_HRES",vid_WAVELENGTH_HRES)
         call wrap_inq_varid(nc_id, "RADIANCE_HRES_CLR",vid_RADIANCE_HRES_CLR)
         call wrap_inq_varid(nc_id, "RADIANCE_HRES_ALL",vid_RADIANCE_HRES_ALL)
         call wrap_inq_varid(nc_id, "SOLAR_FLUX",vid_SOLAR_FLUX)
         call wrap_inq_varid(nc_id, "DIFFUSE_FLUX_CLR",vid_DIFFUSE_FLUX_CLR)
         call wrap_inq_varid(nc_id, "DIFFUSE_FLUX_ALL",vid_DIFFUSE_FLUX_ALL)
         call wrap_inq_varid(nc_id, "SOLAR_ZENITH",vid_SOLAR_ZENITH)
         call wrap_inq_varid(nc_id, "AOD", vid_AOD)
         call wrap_inq_varid(nc_id, "TAU_SULF", vid_TAU_SULF)
         call wrap_inq_varid(nc_id, "TAU_DUST", vid_TAU_DUST)
         call wrap_inq_varid(nc_id, "TAU_SOOT", vid_TAU_SOOT)
         call wrap_inq_varid(nc_id, "TAU_SSLT", vid_TAU_SSLT)
         call wrap_inq_varid(nc_id, "NCONFIG_REAL", vid_NCONFIG_REAL)
         call wrap_inq_varid(nc_id, "CLDTAU_ICE_SW", vid_CLDTAU_ICE_SW)
         call wrap_inq_varid(nc_id, "CLDTAU_LIQ_SW", vid_CLDTAU_LIQ_SW)
         call wrap_inq_varid(nc_id, "CLDTAU_LW", vid_CLDTAU_LW)
         call wrap_inq_varid(nc_id, "VIS_UPDIFFUSE_CLR", vid_VIS_UPDIFFUSE_CLR)
         call wrap_inq_varid(nc_id, "VIS_UPDIFFUSE_ALL", vid_VIS_UPDIFFUSE_ALL)
         call wrap_inq_varid(nc_id, "VIS_DNDIFFUSE_CLR", vid_VIS_DNDIFFUSE_CLR)
         call wrap_inq_varid(nc_id, "VIS_DNDIFFUSE_ALL", vid_VIS_DNDIFFUSE_ALL)
         call wrap_inq_varid(nc_id, "VIS_DNDIRECT_CLR", vid_VIS_DNDIRECT_CLR)
         call wrap_inq_varid(nc_id, "VIS_DNDIRECT_ALL", vid_VIS_DNDIRECT_ALL)
         call wrap_close(nc_id)

  end subroutine get_nc_varids


end module Output
