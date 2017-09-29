#include <misc.h>
#include <params.h>
      
module get_previous_results
      
                                !Purpose:
                                !   get previous results from calculations for restart
                                !   Auther: D. Feldman
      
                                !  
      
      use ppgrid
      use pmgrid
      use shr_kind_mod, only: r8 => shr_kind_r8
      
      implicit none
      
                                !quantities of interest
      integer, public :: plat_start !latitude index to restart calculations
                                !quantities of interest
      
      save
      
contains

      subroutine load_previous_results(opath,itime,FLN,FLNS,FLNSC,FLNT,FLNTC,   &
        FLN200,FLN200C,FSDS,FSDSC,FSN,FSNS,FSNSC,FSNT,FSNTC, &
        FSN200,FSN200C,QRL,QRS,SOLIN,SOLL,SOLLD,SOLS,SOLSD, &
        WAVELENGTH_LRES,RADIANCE_LRES_CLR,RADIANCE_LRES_ALL,WAVELENGTH_HRES, &
        RADIANCE_HRES_CLR,RADIANCE_HRES_ALL,SOLAR_FLUX,DIFFUSE_FLUX_CLR,DIFFUSE_FLUX_ALL, &
        SOLAR_ZENITH,AOD)
!     Read already-written results from netcdf file
      
      use ioFileMod, only: getfil
      implicit none
      
      include 'netcdf.inc'
      
      !Input arguments !DRF
      character(len=120), intent(in) :: opath ! Path to output data.
      integer, intent(in) :: itime             ! time slice
      real(r8), intent(out) :: QRS(plon,pver,plat) ! Shortwave heating rate
      real(r8), intent(out) :: QRL(plon,pver,plat) ! Longwave  heating rate
      real(r8), intent(out) :: FSNS(plon,plat) ! Surface solar absorbed flux
      real(r8), intent(out) :: FSNT(plon,plat) ! Net column abs solar flux at model top
      real(r8), intent(out) :: FLNS(plon,plat) ! Srf longwave cooling (up-down) flux
      real(r8), intent(out) :: FLNT(plon,plat) ! Net outgoing lw flux at model top
      real(r8), intent(out) :: SOLS(plon,plat) ! Direct beam solar rad. onto srf (sw)
      real(r8), intent(out) :: SOLL(plon,plat) ! Direct beam solar rad. onto srf (lw)
      real(r8), intent(out) :: SOLSD(plon,plat) ! Diffuse solar radiation onto srf (sw)
      real(r8), intent(out) :: SOLLD(plon,plat) ! Diffuse solar radiation onto srf (lw)
      real(r8), intent(out) :: FLNSC(plon,plat) ! Net surface clear-sky longwave
      real(r8), intent(out) :: FLNTC(plon,plat) ! Net TOA clear-sky longwave
      real(r8), intent(out) :: FSDS(plon,plat) ! Downwelling surface shortwave
      real(r8), intent(out) :: FSDSC(plon,plat) ! Downwelling clear surface shortwave
      real(r8), intent(out) :: FSNSC(plon,plat) ! Net surface clear shortwave
      real(r8), intent(out) :: FSNTC(plon,plat) ! Net TOA cear shortwave
      real(r8), intent(out) :: SOLIN(plon,plat) ! TOA insolation
      real(r8), intent(out) :: FLN200(plon,plat) ! Net all-sky LW flux at 200 mb
      real(r8), intent(out) :: FLN200C(plon,plat) ! Net clear-sky LW flux at 200 mb
      real(r8), intent(out) :: FSN200(plon,plat) ! Net all-sky SW flux at 200 mb
      real(r8), intent(out) :: FSN200C(plon,plat) ! Net clear-sky SW flux at 200 mb
      real(r8), intent(out) :: FLN(plon,pverp,plat) ! Net longwave flux at each interface
      real(r8), intent(out) :: FSN(plon,pverp,plat) ! Net shortwave flux at each interface
      real*4, intent(out) :: wavelength_lres(wvlng) !Wavelength values from Modtran (DRF)
      real*4, intent(out) :: radiance_lres_clr(plon,wvlng,plat) !Radiance values from Modtran (DRF) clear-sky
      real*4, intent(out) :: radiance_lres_all(plon,wvlng,plat) !Radiance values from Modtran (DRF) all-sky
      real*4, intent(out) :: wavelength_hres(wvlng_hres) !Wavelength values from Modtran (DRF)
      real*4, intent(out) :: radiance_hres_clr(plon,wvlng_hres,plat) !Radiance values from Modtran (DRF) clear-sky
      real*4, intent(out) :: radiance_hres_all(plon,wvlng_hres,plat) !Radiance values from Modtran (DRF) all-sky
      real*4, intent(out) :: solar_flux(plon,wvlng,plat) ! TOA downwelling solar flux values from Modtran (DRF)
      real*4, intent(out) :: diffuse_flux_clr(plon,wvlng,plat) !Diffuse shortwave flux from Modtran (DRF) clear-sky
      real*4, intent(out) :: diffuse_flux_all(plon,wvlng,plat) !Diffuse shortwave flux from Modtran (DRF) all-sky
      real*4, intent(out) :: solar_zenith(plon,plat) !Solar zenith angle from Modtran (DRF)
      real*4, intent(out) :: aod(plon,plat) !Solar zenith angle from Modtran (DRF)

      !local variables
      integer :: i,j,k,m,n
      integer :: lat_id, lon_id,lev_id,ilev_id,wvl_lres_id,wvl_hres_id
      integer :: index1,index2
      integer :: rad_lres_clr_id,rad_lres_all_id
      integer :: rad_hres_clr_id,rad_hres_all_id
      integer :: diffuse_flux_clr_id,diffuse_flux_all_id
      integer :: solar_flux_id,aod_id,solar_zenith_id
     
      !Regular cam_output fields
      integer :: qrs_id,qrl_id,fsns_id,fsnt_id,flns_id,flnt_id,sols_id,soll_id,solsd_id
      integer :: solld_id,flnsc_id,flntc_id,fsds_id,fsdsc_id,fsnsc_id,fsntc_id,solin_id
      integer :: fln200_id,fln200c_id,fsn200_id,fsn200c_id,fln_id,fsn_id

      character(len=256) :: locfn !local file
      integer :: nc_id
      
      integer :: num_lat,num_lon,num_lev,num_ilev,num_lres,num_hres
      logical :: lexist
      integer :: start_4d(4)
      integer :: count_4d(4)
      character(len = 10) :: routine = "input_data"
      integer :: istat                         ! allocate status

      !Readin variables
      real(r8), allocatable :: wvl_lres_readin(:)
      real(r8), allocatable :: wvl_hres_readin(:)
      real(r8), allocatable :: fln_readin(:,:,:)
      real(r8), allocatable :: fsn_readin(:,:,:)
      real(r8), allocatable :: qrs_readin(:,:,:)
      real(r8), allocatable :: qrl_readin(:,:,:)
      real(r8), allocatable :: diffuse_flux_clr_readin(:,:,:)
      real(r8), allocatable :: diffuse_flux_all_readin(:,:,:)
      real(r8), allocatable :: solar_flux_readin(:,:,:)
      real(r8), allocatable :: radiance_lres_clr_readin(:,:,:)
      real(r8), allocatable :: radiance_lres_all_readin(:,:,:)
      real(r8), allocatable :: radiance_hres_clr_readin(:,:,:)
      real(r8), allocatable :: radiance_hres_all_readin(:,:,:)
      real(r8), allocatable :: solar_zenith_readin(:,:)
      real(r8), allocatable :: aod_readin(:,:)

      !Dummy variable
      real(r8), allocatable :: dummy_array(:)
      
      write (6, '(2x, a)') '_______________________________________________________'
      write (6, '(2x, a)') '_______ reading results from previous runs   __________'
      write (6, '(2x, a)') '_______________________________________________________'
      
      plat_start = 1
      
         !Check to see if output file (opath) already exists
         call getfil(opath, locfn, 0)
         call wrap_open(locfn, 0, nc_id)

         !Get data dimensions
         call wrap_inq_dimid(nc_id, "lat", lat_id)
         call wrap_inq_dimid(nc_id, "lon", lon_id)
         call wrap_inq_dimid(nc_id, "lev",lev_id)
         call wrap_inq_dimid(nc_id, "ilev",ilev_id)
         call wrap_inq_dimid(nc_id, "wavelength",wvl_lres_id)
         call wrap_inq_dimid(nc_id, "wavelength_hres",wvl_hres_id)
         
         call wrap_inq_dimlen(nc_id, lat_id, num_lat)
         call wrap_inq_dimlen(nc_id, lon_id, num_lon)
         call wrap_inq_dimlen(nc_id, lev_id, num_lev)
         call wrap_inq_dimlen(nc_id, ilev_id, num_ilev)
         call wrap_inq_dimlen(nc_id, wvl_lres_id, num_lres)
         call wrap_inq_dimlen(nc_id, wvl_hres_id, num_hres)
         
         !Get ID for the regular fields
         call wrap_inq_varid(nc_id, "QRS", qrs_id)
         call wrap_inq_varid(nc_id, "QRL", qrl_id)
         call wrap_inq_varid(nc_id, "FSNS", fsns_id)
         call wrap_inq_varid(nc_id, "FSNT", fsnt_id)
         call wrap_inq_varid(nc_id, "FLNS", flns_id)
         call wrap_inq_varid(nc_id, "FLNT", flnt_id)
         call wrap_inq_varid(nc_id, "SOLS", sols_id)
         call wrap_inq_varid(nc_id, "SOLL", soll_id)
         call wrap_inq_varid(nc_id, "SOLSD", solsd_id)
         call wrap_inq_varid(nc_id, "SOLLD", solld_id)
         call wrap_inq_varid(nc_id, "FLNSC", flnsc_id)
         call wrap_inq_varid(nc_id, "FLNTC", flntc_id)
         call wrap_inq_varid(nc_id, "FSDS", fsds_id)
         call wrap_inq_varid(nc_id, "FSDSC", fsdsc_id)
         call wrap_inq_varid(nc_id, "FSNSC", fsnsc_id)
         call wrap_inq_varid(nc_id, "FSNTC", fsntc_id)
         call wrap_inq_varid(nc_id, "SOLIN", solin_id)
         call wrap_inq_varid(nc_id, "FLN200", fln200_id)
         call wrap_inq_varid(nc_id, "FLN200C", fln200c_id)
         call wrap_inq_varid(nc_id, "FSN200", fsn200_id)
         call wrap_inq_varid(nc_id, "FSN200C", fsn200c_id)
         call wrap_inq_varid(nc_id, "FLN", fln_id)
         call wrap_inq_varid(nc_id, "FSN", fsn_id)

         !Get ID numbers of the Modtran added data fields
         call wrap_inq_varid(nc_id, "WAVELENGTH_LRES", wvl_lres_id)
         call wrap_inq_varid(nc_id, "RADIANCE_LRES_CLR", rad_lres_clr_id)
         call wrap_inq_varid(nc_id, "RADIANCE_LRES_ALL", rad_lres_all_id)
         call wrap_inq_varid(nc_id, "WAVELENGTH_HRES",wvl_hres_id)
         call wrap_inq_varid(nc_id, "RADIANCE_HRES_CLR",rad_hres_clr_id)
         call wrap_inq_varid(nc_id, "RADIANCE_HRES_ALL",rad_hres_all_id)
         call wrap_inq_varid(nc_id, "SOLAR_FLUX",solar_flux_id)
         call wrap_inq_varid(nc_id, "DIFFUSE_FLUX_CLR",diffuse_flux_clr_id)
         call wrap_inq_varid(nc_id, "DIFFUSE_FLUX_ALL",diffuse_flux_all_id)
         call wrap_inq_varid(nc_id, "SOLAR_ZENITH",solar_zenith_id)
         call wrap_inq_varid(nc_id, "AOD", aod_id)

         start_4d = (/1,   1,   itime, -1/) 
         count_4d = (/plon,plat,1,     -1/) 
 
         !Read in 2-d variables (assuming 1 time-step) 
         call wrap_get_vara_realx(nc_id, fsns_id, start_4d, count_4d,FSNS )         
         call wrap_get_vara_realx(nc_id, fsnt_id, start_4d, count_4d,FSNT )         
         call wrap_get_vara_realx(nc_id, flns_id, start_4d, count_4d,FLNS )         
         call wrap_get_vara_realx(nc_id, flnt_id, start_4d, count_4d,FLNT )         
         call wrap_get_vara_realx(nc_id, sols_id, start_4d, count_4d,SOLS )         
         call wrap_get_vara_realx(nc_id, soll_id, start_4d, count_4d,SOLL )         
         call wrap_get_vara_realx(nc_id, solsd_id, start_4d, count_4d,SOLSD )         
         call wrap_get_vara_realx(nc_id, solld_id, start_4d, count_4d,SOLLD )         
         call wrap_get_vara_realx(nc_id, flnsc_id, start_4d, count_4d,FLNSC )         
         call wrap_get_vara_realx(nc_id, flntc_id, start_4d, count_4d,FLNTC )         
         call wrap_get_vara_realx(nc_id, fsds_id, start_4d, count_4d,FSDS )         
         call wrap_get_vara_realx(nc_id, fsdsc_id, start_4d, count_4d,FSDSC )         
         call wrap_get_vara_realx(nc_id, fsnsc_id, start_4d, count_4d,FSNSC )         
         call wrap_get_vara_realx(nc_id, fsntc_id, start_4d, count_4d,FSNTC )         
         call wrap_get_vara_realx(nc_id, solin_id, start_4d, count_4d,SOLIN )         
         call wrap_get_vara_realx(nc_id, fln200_id, start_4d, count_4d,FLN200 )         
         call wrap_get_vara_realx(nc_id, fln200c_id, start_4d, count_4d,FLN200C )         
         call wrap_get_vara_realx(nc_id, fsn200_id, start_4d, count_4d,FSN200 )         
         call wrap_get_vara_realx(nc_id, fsn200c_id, start_4d, count_4d,FSN200C )         

 
         allocate(solar_zenith_readin(plon,plat),stat=istat)
         call alloc_err(istat, routine, "solar_zenith_readin",(plon*plat) )
         call wrap_get_vara_realx(nc_id, solar_zenith_id, start_4d, count_4d, solar_zenith_readin)         

         allocate(aod_readin(plon,plat),stat=istat)
         call alloc_err(istat, routine, "aod_readin",(plon*plat) )
         call wrap_get_vara_realx(nc_id, aod_id, start_4d, count_4d, aod_readin)        
 
         plat_start = 0
         do i=1,plat
            if (solar_zenith_readin(1,i).gt.0.0) then
               plat_start = i
            endif
         end do
          
         !Read in 3-d variables, normal data fields
         start_4d = (/1,   1,   1,     itime/) 
         count_4d = (/plon,plat,pverp, 1/)
         allocate(fln_readin(plon,plat,pverp),stat=istat)
         call alloc_err(istat, routine, "fln_readin",(plon*plat*pverp) )
         call wrap_get_vara_realx(nc_id, fln_id, start_4d, count_4d, fln_readin)         
         allocate(fsn_readin(plon,plat,pverp),stat=istat)
         call alloc_err(istat, routine, "fsn_readin",(plon*plat*pverp) )
         call wrap_get_vara_realx(nc_id, fsn_id, start_4d, count_4d, fsn_readin)         
         
         start_4d = (/1,   1,   1,     itime/) 
         count_4d = (/plon,plat,pver, 1/)
         allocate(qrs_readin(plon,plat,pver),stat=istat)
         call alloc_err(istat, routine, "qrs_readin",(plon*plat*pver) )
         call wrap_get_vara_realx(nc_id, qrs_id, start_4d, count_4d, qrs_readin)         
         allocate(qrl_readin(plon,plat,pver),stat=istat)
         call alloc_err(istat, routine, "qrl_readin",(plon*plat*pver) )
         call wrap_get_vara_realx(nc_id, qrl_id, start_4d, count_4d, qrl_readin)         


         !Read in 3-d variables (assuming 1 time-stp), low-res
         allocate(diffuse_flux_clr_readin(plon,plat,wvlng),stat=istat)
         call alloc_err(istat, routine, "diffuse_flux_clr_readin",(plon*plat*wvlng) )
         allocate(diffuse_flux_all_readin(plon,plat,wvlng),stat=istat)
         call alloc_err(istat, routine, "diffuse_flux_all_readin",(plon*plat*wvlng) )
         allocate(solar_flux_readin(plon,plat,wvlng),stat=istat)
         call alloc_err(istat, routine, "solar_flux_readin",(plon*plat*wvlng) )
         allocate(radiance_lres_clr_readin(plon,plat,wvlng),stat=istat)
         call alloc_err(istat, routine, "radiance_lres_clr_readin",(plon*plat*wvlng) )
         allocate(radiance_lres_all_readin(plon,plat,wvlng),stat=istat)
         call alloc_err(istat, routine, "radiance_lres_all_readin",(plon*plat*wvlng) )

         start_4d = (/1,   1,   1,     itime/) 
         count_4d = (/plon,plat,wvlng, 1/) 
         call wrap_get_vara_realx(nc_id, diffuse_flux_clr_id, start_4d, count_4d, diffuse_flux_clr_readin)         
         call wrap_get_vara_realx(nc_id, diffuse_flux_all_id, start_4d, count_4d, diffuse_flux_all_readin)         
         call wrap_get_vara_realx(nc_id, solar_flux_id, start_4d, count_4d, solar_flux_readin)         
         call wrap_get_vara_realx(nc_id, rad_lres_clr_id, start_4d, count_4d, radiance_lres_clr_readin)         
         call wrap_get_vara_realx(nc_id, rad_lres_all_id, start_4d, count_4d, radiance_lres_all_readin)         

         !Read in 3-d variables (assuming 1 time-stp), high-res
         allocate(radiance_hres_clr_readin(plon,plat,wvlng_hres),stat=istat)
         call alloc_err(istat, routine, "radiance_hres_clr_readin",(plon*plat*wvlng_hres) )
         allocate(radiance_hres_all_readin(plon,plat,wvlng_hres),stat=istat)
         call alloc_err(istat, routine, "radiance_hres_all_readin",(plon*plat*wvlng_hres) )
         start_4d = (/1,   1,   1,     itime/) 
         count_4d = (/plon,plat,wvlng_hres, 1/) 
         call wrap_get_vara_realx(nc_id, rad_hres_clr_id, start_4d, count_4d, radiance_hres_clr_readin)         
         call wrap_get_vara_realx(nc_id, rad_hres_all_id, start_4d, count_4d, radiance_hres_all_readin)         

         allocate(wvl_lres_readin(wvlng),stat=istat)
         call alloc_err(istat, routine, "wvl_lres_readin",(wvlng) )
         start_4d = (/1,   1,   1,     itime/) 
         allocate(wvl_hres_readin(wvlng_hres),stat=istat)
         call alloc_err(istat, routine, "wvl_hres_readin",(wvlng_hres) )
         call wrap_get_var_realx(nc_id,wvl_lres_id,wvl_lres_readin)
         call wrap_get_var_realx(nc_id,wvl_hres_id,wvl_hres_readin)

        !Convert readin values to actual outputs
         do i=1,wvlng
            wavelength_lres(i) = real(wvl_lres_readin(i))
         end do
         do i=1,wvlng_hres
            wavelength_hres(i) = real(wvl_hres_readin(i))
         end do

          do j=1,plat_start
             do i=1,plon
                FLN(i,:,j) = fln_readin(i,j,:)
                FSN(i,:,j) = fsn_readin(i,j,:)
                QRS(i,:,j) = qrs_readin(i,j,:)
                QRL(i,:,j) = qrl_readin(i,j,:)
                solar_zenith(i,j) = real(solar_zenith_readin(i,j))                
                aod(i,j) = real(aod_readin(i,j))                
                diffuse_flux_clr(i,:,j) = real(diffuse_flux_clr_readin(i,j,:))                
                diffuse_flux_all(i,:,j) = real(diffuse_flux_all_readin(i,j,:))                
                solar_flux(i,:,j) = real(solar_flux_readin(i,j,:))                
                radiance_lres_clr(i,:,j) = real(radiance_lres_clr_readin(i,j,:))                
                radiance_lres_all(i,:,j) = real(radiance_lres_all_readin(i,j,:))                
                radiance_hres_clr(i,:,j) = real(radiance_hres_clr_readin(i,j,:))                
                radiance_hres_all(i,:,j) = real(radiance_hres_all_readin(i,j,:))                
             end do
         end do
         plat_start = plat_start + 1

         write(*,*) 'read successfully,restarting at index ',plat_start

    !free up resources
         deallocate(wvl_lres_readin)
         deallocate(wvl_hres_readin)
         deallocate(fln_readin)
         deallocate(fsn_readin)
         deallocate(qrs_readin)
         deallocate(qrl_readin)
         deallocate(solar_zenith_readin)
         deallocate(aod_readin)
         deallocate(diffuse_flux_clr_readin)
         deallocate(diffuse_flux_all_readin)
         deallocate(solar_flux_readin)
         deallocate(radiance_lres_clr_readin)
         deallocate(radiance_lres_all_readin)
         deallocate(radiance_hres_clr_readin)
         deallocate(radiance_hres_all_readin)
         call wrap_close(nc_id)
  end subroutine load_previous_results

  subroutine starting_from_scratch(opath,create_output_flag)
! Purpose:
!   To check to see whether the output file alread exists
!
      use ioFileMod, only: getfil
      implicit none
 
      character(len=120), intent(in) :: opath ! Path to output data.
      logical, intent(out) :: create_output_flag

      !local variable
      logical :: lexist
      inquire (file=opath,exist=lexist)

      if (lexist) then
        create_output_flag = .false.
      else
        create_output_flag = .true.
      endif

  end subroutine starting_from_scratch

end module get_previous_results


