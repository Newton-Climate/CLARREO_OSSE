#include <misc.h>
#include <params.h>

module write_modtran_tp5

   use ppgrid
   use pmgrid
   implicit none
   save

contains

subroutine modtran_out_netcdf(wavelength_lres,wavelength_hres,radiance_lres, radiance_hres, &
                              solar_flux,solar_zenith, diffuse_flux,&
                              bb_updiffuse,bb_dndiffuse,bb_dndirect, &
                              fname)

!---------------
!
! Purpose: to write a single spectrum output from Modtran to a netcdf file
!          differs from the regular modtran_write_tp5 in that this is output from the internal
!          call to Modtran
!
!
!  Author: D. Feldman, U.C. Berkekely, 2010
!

     use ioFileMod, only: getfil
     implicit none

     include 'netcdf.inc'

     real*4, intent(in) :: wavelength_lres(wvlng2)  ! Wavelength from Modtran (DRF)
     real*4, intent(in) :: wavelength_hres(wvlng_hres2)  ! Wavelength from Modtran (DRF)
     real*4, intent(in) :: radiance_lres(wvlng2)   ! Radiance from Modtran (DRF)
     real*4, intent(in) :: radiance_hres(wvlng_hres2)   ! Radiance from Modtran (DRF)
     real*4, intent(in) :: solar_flux(wvlng2)   ! TOA solar flux from Modtran (DRF)
     real*4, intent(in) :: solar_zenith  !Solar zenith angle from Modtran (DRF)
     real*4, intent(in) :: diffuse_flux(wvlng2)   ! TOA upwelling reflected flux from Modtran (DRF)
     character(len=*), intent(in) :: fname ! path to write to netcdf file (DRF)
     real*4, intent(in) :: bb_updiffuse(pver)
     real*4, intent(in) :: bb_dndiffuse(pver)
     real*4, intent(in) :: bb_dndirect(pver)


     !local variables
     integer :: vid_WAVELENGTH_LRES           ! Variable id for WAVELENGTH_LRES
     integer :: vid_RADIANCE_LRES             ! Variable id for RADIANCE_LRES
     integer :: vid_WAVELENGTH_HRES           ! Variable id for WAVELENGTH_LRES
     integer :: vid_RADIANCE_HRES             ! Variable id for RADIANCE_LRES
     integer :: vid_SOLAR_FLUX                ! Variable id for SOLAR_FLUX
     integer :: vid_DIFFUSE_FLUX              ! Variable id for DIFFUSE_FLUX
     integer :: vid_BB_UPDIFFUSE            ! Variable id for BB_UPDIFFUSE
     integer :: vid_BB_DNDIFFUSE            ! Variable id for BB_DNDIFFUSE
     integer :: vid_BB_DNDIRECT             ! Variable id for BB_DNDIRECT

     integer :: vid_SOLAR_ZENITH              ! Variable id for SOLAR_ZENITH
     real*4  ::  bufferr(2,pver)

!
! Dimension ids
!
     integer :: did_lat                       ! Latitude dimension ID
     integer :: did_lon                       ! Longitude dimension ID
!DRF
     integer :: did_wvl                       ! Wavelength dimension ID
     integer :: did_wvl_hres                  ! Wavelength dimension ID
     integer :: did_plev                      ! level dimension ID
     integer :: did_angle                     ! angle dimension ID
     integer :: did_band                     ! band dimension ID

     integer :: plat_local
     integer :: plon_local
     integer :: plev_local
     integer :: angle_local
     integer :: band_local

     integer :: vdims                      ! dimension ids
     integer :: vdims2 !DRF
     integer :: vdims3 !DRF
     integer :: vdims4(2) !DRF
     integer :: nvdims                        ! number of dimensions
     integer :: ret                           ! return code

!
! File ID
!
     integer :: nfid                          ! NetCDF id

!Dimension info
    integer :: start                 ! starting point on edges
    integer :: count                 ! number of elements to write
    integer :: start2(2)                 ! starting point on edges
    integer :: count2(2)                 ! number of elements to write

!!!!!!!!!!!!
! Begin code 
     call wrap_create(fname, nf_clobber, nfid)

!
! Define dimensions
!    
    plat_local = 1
    plon_local = 1
    plev_local = 26
    band_local = 2

    call wrap_def_dim(nfid, "lat", plat_local, did_lat)
    call wrap_def_dim(nfid, "lon", plon_local, did_lon)
    call wrap_def_dim(nfid,"wavelength",wvlng,did_wvl)
    call wrap_def_dim(nfid,"wavelength_hres",wvlng_hres,did_wvl_hres)
    call wrap_def_dim(nfid,"levs",plev_local,did_plev)
    call wrap_def_dim(nfid,"angle",angle_local,did_angle)
    call wrap_def_dim(nfid,"band",band_local,did_band)

!
! Define variables
!
    nvdims = 1
    vdims = did_wvl
    vdims2 = did_wvl_hres
    vdims3 = did_angle

    !DRF
    call wrap_def_var(nfid, "WAVELENGTH_LRES",nf_real, nvdims, vdims, vid_WAVELENGTH_LRES)
    call wrap_def_var(nfid, "WAVELENGTH_HRES",nf_real, nvdims, vdims2, vid_WAVELENGTH_HRES)
    call wrap_def_var(nfid, "RADIANCE_LRES",nf_real, nvdims, vdims, vid_RADIANCE_LRES)
    call wrap_def_var(nfid, "RADIANCE_HRES",nf_real, nvdims, vdims2, vid_RADIANCE_HRES)
    call wrap_def_var(nfid, "SOLAR_FLUX",nf_real, nvdims, vdims, vid_SOLAR_FLUX)
    call wrap_def_var(nfid, "DIFFUSE_FLUX",nf_real, nvdims, vdims, vid_DIFFUSE_FLUX)
    call wrap_def_var(nfid, "SOLAR_ZENITH",nf_real, nvdims, vdims3, vid_SOLAR_ZENITH)

    nvdims = 1
    vdims3 = did_plev
    call wrap_def_var(nfid, "BB_UPDIFFUSE",nf_real, nvdims,vdims3, vid_BB_UPDIFFUSE)
    call wrap_def_var(nfid, "BB_DNDIFFUSE",nf_real, nvdims,vdims3, vid_BB_DNDIFFUSE)
    call wrap_def_var(nfid, "BB_DNDIRECT",nf_real, nvdims,vdims3, vid_BB_DNDIRECT)

    !Write attributes to file
    !ret = nf_put_att_int (nfid, nf_global, "date", nf_int, nslice, date)
    call wrap_close(nfid)

    !!!!!!!!!!!!!!!!!!!!!!!!!
    !Write actual data to file
    call wrap_open(fname, nf_write, nfid)
    start = 1
    count = wvlng
    
    call wrap_put_vara_real(nfid, vid_WAVELENGTH_LRES,start, count, wavelength_lres) !DRF
    call wrap_put_vara_real(nfid, vid_RADIANCE_LRES,start, count, radiance_lres)  !DRF
    call wrap_put_vara_real(nfid, vid_SOLAR_FLUX,start, count, solar_flux)  !DRF
    call wrap_put_vara_real(nfid, vid_DIFFUSE_FLUX,start, count, diffuse_flux)  !DRF

    count = wvlng_hres
    call wrap_put_vara_real(nfid, vid_WAVELENGTH_HRES,start,count, wavelength_hres) !DRF
    call wrap_put_vara_real(nfid, vid_RADIANCE_HRES,start,count, radiance_hres) !DRF

    count = 1
    call wrap_put_vara_real(nfid, vid_SOLAR_ZENITH,start,count, solar_zenith) !DRF

    start = 1
    count = pver
    call wrap_put_vara_real(nfid, vid_BB_UPDIFFUSE,start, count,bb_updiffuse)  !DRF
    call wrap_put_vara_real(nfid, vid_BB_DNDIFFUSE,start, count,bb_dndiffuse)  !DRF
    call wrap_put_vara_real(nfid, vid_BB_DNDIRECT,start, count,bb_dndirect)  !DRF

    call wrap_close(nfid)

    return

end subroutine modtran_out_netcdf


subroutine modtran4_write_tp5(surf_temp,co2vmr,surf_alt,jday,gmt,zen_real,az_real,a_flag,  &
                             full_atm,v1,v2,ocean_flag,brdf_len,brdf_wvl,brdf_param,    &
                             spec_albedo_len,spec_albedo_wvl,spec_emissivity_vals, &
                             c_flag,num_levs,m_aerosol,o_aerosol,m_cloud_liq,&
                             m_cloud_ice,wvl_cld,liq_ext_cld,liq_ssa_cld,liq_asym_cld, &
                             ice_ext_cld,ice_ssa_cld,ice_asym_cld,z_cld_levs,          &
                             aer_wavelen,fname,fname_length,phase_funcs,phase_wvls,    &
                             phase_angles,dv,tp5_file,debug_flag,lw_flag)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! write a tape5 file for analysis that can be used in Modtran5.3.0.0 externally
! 
!
! Author: D. Feldman U.C. Berkeley
! 
!-----------------------------------------------------------------------
!
! Input arguments
!

    real*4 :: surf_temp
    real*4 :: co2vmr
    real*4 :: surf_alt
    integer jday,num_levs
    real*4 :: gmt
    logical a_flag
    real*4 :: full_atm(num_levs,28)
    real*4 :: v1
    real*4 :: v2
    logical ocean_flag
    integer brdf_len
    real*4 :: brdf_wvl(brdf_len)
    real*4 :: brdf_param(brdf_len,3)
    integer spec_albedo_len
    real*4 :: spec_albedo_wvl(6)
    real*4 :: spec_emissivity_vals(6)
    logical c_flag
    real*4 :: m_aerosol(num_levs,4)
    real*4 :: o_aerosol(42,12)
    real*4 :: m_cloud_liq(num_levs)
    real*4 :: m_cloud_ice(num_levs)
    real*4 :: wvl_cld(24)
    real*4 :: liq_ext_cld(24)
    real*4 :: liq_ssa_cld(24)
    real*4 :: liq_asym_cld(24)
    real*4 :: ice_ext_cld(24)
    real*4 :: ice_ssa_cld(24)
    real*4 :: ice_asym_cld(24)
    real*4 :: z_cld_levs(num_levs-1)
    real*4 :: aer_wavelen(42)
    character*120 fname
    integer fname_length
    real*4 :: phase_funcs(50,15,4)
    real*4 :: phase_wvls(15)
    real*4 :: phase_angles(50)
    real*4:: dv
    real*4 :: zen_real
    real*4 :: az_real
    logical debug_flag
    logical lw_flag
    real*4  emis_broadband, alb_broadband

!   output argument
    character, intent(out):: tp5_file(1000)*110

!
!---------------------------Local variables-----------------------------
!
   integer i, j, k, jj, kk              ! index
   real*4 ii
   logical c_flag_local
   logical a_flag_local
   integer tp5_linenumber

!-----------------------------------------------------------------------
!Modtran tape5 variables
!

!Card 1
   character MODTRN,SPEED
   integer MODEL,ITYPE,IEMSCT,IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT
   real*4 TPTEMP
   character*7 SURREF

!Card 1A
   character DIS,DISAZM
   integer NSTR
   integer ISUN
   character LSUN
   real*4 CO2MX
   character*10 H2OSTR,O3STR
   character LSUNFL,LBMNAM,LFLTNM,H2OAER,LDATDR
   real*4 SOLCON

!Card 1A2
   character*17 BMNAME

!Card 2
   character*2 APLUS
   integer IHAZE
   character CNOVAM
   integer ISEASN
   character*3 ARUSS
   integer IVULCN,ICSTL,ICLD,IVSA
   real*4 VIS,WSS,WHH,RAINRT,GNDALT

!Optional Card 2A
   real*4 CTHIK,CALT,CEXT
   integer NCRALT,NCRSPC
   real*4 CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP

!Card 2C
   integer ML,IRD1,IRD2
   character*20 HMODEL
   real*4 REE
   integer NMOLYC

!Cards 2C1,2C2,2C2X
   real*4 ZM,P,T,WMOL(12),WMOLX(13)
   character*14 JCHAR
   character JCHARX

!Card 2C3
   real*4 AHAZE(4),RRATZ

!Card 2D,2D1,2D2
   integer IREG(4)
   real*4 AWCCON
   character*18 TITLE
   real*4 VARSPC(42),EXTC(42),ABSC(42),ASYM(42) 

!Card 2E1,2E2
   real*4 ZCLD,CLD,CLDICE,RR
   real*4 WAVLEN,EXTC_LIQ,ABSC_LIQ,ASYM_LIQ
   real*4 EXTC_ICE,ABSC_ICE,ASYM_ICE

!Card 3
   real*4 H1,H2,ANGLE,mRANGE,BETA,RO
   integer LENN
   real*4 PHI

!Cards 3A1,3A2
   integer IPARM,IPH,IDAY,ISOURC
   real*4 PARM1,PARM2,PARM3,PARM4,TIME,PSIPO,ANGLEM,G

! Card 3B1
   integer NANGLS,NWLF

!Card 3C1
   real*4 ANGF(50)

!Card 3C2
   real*4 WLF(15)

!Card 3C3-6
   real*4 F(15)

!Card 4
  real*4 IV1,IV2,IDV,FWHM
  character YFLAG,XFLAG
  character*8 DLIMIT
  character*7 FLAGS
  integer MLFLX

!Card 4A
   integer NSURF
   real*4 AATEMP
   real*4 DH2O
   character MLTRFL

!Card 4B1
   character*10 CBRDF

!Card 4B2
   integer NWVSRF
   real*4 SURFZN
   real*4 SURFAZ     

!Card 4B3
   real*4 WVSURF
   real*4 PARAMS(5)

!DRF CHANGE
!Card 4L1
   integer NWVALB
   
!Card 4L2
   real*4 WVSURFALB
   real*4 SURFALBVAL
!DRF CHANGE
   
 
!Card 5
   integer IRPT

   logical debug_flag_local

!-----------------------------------------------------------------------
!Open file

   debug_flag_local = debug_flag
   if (debug_flag_local) then
     open(unit=50,file=fname(1:fname_length)//'mod4.tp5',status='new')
   endif
   c_flag_local = c_flag
   a_flag_local = a_flag
   !c_flag_local = .false.
   !a_flag_local = .false.

!Assign Card 1
   MODTRN = 'C'
   SPEED = 'M'
   MODEL = 7
   ITYPE = 3
   if (lw_flag) then
     IEMSCT = 1
     IMULT = 1
   else  
     IEMSCT = 2
     IMULT = -1
   endif
   M1 = 0
   M2 = 0
   M3 = 0
   M4 = 0
   M5 = 0
   M6 = 0
   MDEF = 2
   IM = 1
   NOPRNT = 0
   TPTEMP = surf_temp
   if (lw_flag) then
     !convert emissivity from spectral to broadband
     !this follows the publication Bo-Hui Tang et al, Optics Express 19(1) 185-192 (2011)
     !"Estimation of broadband surface emissivity from narrowband emissivities"
     !emis_broadband = 0.0127+0.7852*spec_emissivity_vals(4) &
     !                 -0.0151*spec_emissivity_vals(5)+0.139*spec_emissivity_vals(6) 
     emis_broadband = 0.999
     !write(*,*) 'emis_broadband = ',emis_broadband
     alb_broadband = 1.-emis_broadband
     !write(*,*) 'alb_broadband = ',alb_broadband
   else
     SURREF = 'BRDF'
   endif
!Write Card 1
   if (debug_flag_local) then
       if (.not.lw_flag) then
         write(50,'(2A1,I3,12I5,F8.3,A7)') MODTRN,SPEED,MODEL,ITYPE,IEMSCT, &
          IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT,TPTEMP,SURREF
         write(tp5_file(1),'(2A1,I3,12I5,F8.3,A7)') MODTRN,SPEED,MODEL,ITYPE,IEMSCT, &
          IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT,TPTEMP,SURREF
       else
         write(50,'(2A1,I3,12I5,F8.3,F7.5)') MODTRN,SPEED,MODEL,ITYPE,IEMSCT, &
          IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT,TPTEMP,alb_broadband
         write(tp5_file(1),'(2A1,I3,12I5,F8.3,F7.5)') MODTRN,SPEED,MODEL,ITYPE,IEMSCT, &
          IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT,TPTEMP,alb_broadband
       endif
   endif

!Assign Card 1A
   if (lw_flag) then 
     DIS = 'F'
   else
     DIS = 'T'
   endif
   DISAZM = 'F'
   NSTR = 8
   LSUN = 'F'
   ISUN = 1
   CO2MX = co2vmr
   H2OSTR = '       1.0'
   O3STR = '       1.0'
   LSUNFL = 'F'
   LBMNAM = 'T'
   LFLTNM = 'F'
   H2OAER = 'F'
   LDATDR = 'F'
   SOLCON = 1360.0

!Write Card 1A
   if (debug_flag_local) then
      write(50,'(2A1,I3,A1,I4,F10.5,2A10,4(1X,A1),(1X,A1),F10.3)') &
        DIS,DISAZM,NSTR,LSUN,ISUN,CO2MX,H2OSTR,O3STR,LSUNFL,LBMNAM, &
        LFLTNM,H2OAER,LDATDR,SOLCON
   endif

!Assign Card 1A2
   if (lw_flag) then
     BMNAME = 'DATA/B2001_01.BIN'
   else
     BMNAME = 'DATA/B2001_15.BIN'
   endif
!Write Card 1A2
   if (debug_flag_local) then
      write(50,'(A17)') BMNAME
   endif

!Assign Card 2
   APLUS = ' '
   CNOVAM = ' '
   ISEASN = 0
   if (a_flag_local .and. .not.lw_flag) then
     ARUSS = 'USS'
     IHAZE = 7
   else
     ARUSS = '   '
     IHAZE = 0
   endif
   IVULCN = 0
   ICSTL = 3
   if (c_flag_local) then
      ICLD = 1
   else
      ICLD = 0
   endif
   IVSA = 0
   VIS = 0.
   WSS = 0.
   WHH = 0.
   RAINRT = 0.
   GNDALT = surf_alt

!Write Card 2
   if (debug_flag_local) then
     write(50,'(A2,I3,A1,I4,A3,I2,3I5,5F10.4)')                    &
       APLUS,IHAZE,CNOVAM,ISEASN,ARUSS,                            &
       IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT
   endif

!Optional Card 2A
   if (c_flag_local) then
     CTHIK = -9.
     CALT = -9.
     CEXT = -9.
     NCRALT = num_levs
     NCRSPC = 24
     CWAVLN = 0.55
     CCOLWD = -9.
     CCOLIP = -9.
     CHUMID = 110.
     ASYMWD = -9.
     ASYMIP = -9.

     if (debug_flag_local) then
        write(50,'(3F8.3,2I4,6F8.3)')  &
          CTHIK,CALT,CEXT,NCRALT,NCRSPC,CWAVLN,CCOLWD,CCOLIP,CHUMID, &
          ASYMWD,ASYMIP
     endif
   endif

!Assign Card 2C
   ML = num_levs
   IRD1 = 1

   if (a_flag_local .or. c_flag_local) then
     IRD2 = 2
   else
     IRD2 = 0
   endif
   HMODEL = '   USERDEFINED'
   REE = 0.0

!Write Card 2C
   if (debug_flag_local) then
      write(50,'(3I5,A20,F10.0)')                               &
         ML,IRD1,IRD2,HMODEL,REE
   endif

!-------------

   if (lw_flag) then 
     JCHAR = 'AAAAAAAAAAAAAA'
   else
     JCHAR = 'AAAAAAAA1AAAAA'
   endif
   JCHARX = 'A'
   do i=1,num_levs

     !Assign Cards 2C1,2C2,2C2X
     ZM = full_atm(i,1)
     P = full_atm(i,2)
     T = full_atm(i,3)
     do j=1,12
       WMOL(j) = full_atm(i,j+3)
     end do
     do j=1,13
       WMOLX(j) = full_atm(i,j+15)
     end do
     AHAZE(1) = m_aerosol(i,1)
     RRATZ = 0.
     AHAZE(2) = m_aerosol(i,2)
     AHAZE(3) = m_aerosol(i,3)
     AHAZE(4) = m_aerosol(i,4)

     !Write Cards 2C1
     if (debug_flag_local) then
       !Write Card 2C1
       write(50,'(F10.4,5F10.4,1A14,1X,A1)') &
         ZM,P,T,WMOL(1),WMOL(2),WMOL(3),JCHAR,JCHARX

        !Write Card 2C2
        write(50,'(8F10.4,/F10.4)') &
          WMOL(4),WMOL(5),WMOL(6),WMOL(7),WMOL(8),WMOL(9), &
          WMOL(10),WMOL(11),WMOL(12)

     !Write Card 2C2X
       write(50,'(8F10.4,/5F10.4)') &
         WMOLX(1),WMOLX(2),WMOLX(3),WMOLX(4),WMOLX(5),WMOLX(6),WMOLX(7),WMOLX(8),WMOLX(9), &
         WMOLX(10),WMOLX(11),WMOLX(12),WMOLX(13)
     endif

     if (IRD2.eq.2) then
       !Write Card 2C3
       write(50,'(10X,F10.7,10X,4F10.7)') &
          AHAZE(1),RRATZ,AHAZE(2),AHAZE(3),AHAZE(4)
     endif
   end do

!-------------

!!! Apparently Card 2E1/2E2 read before 2D1/2D2

if (c_flag_local) then
  !Card 2E1, 2E2

  do j=1,26
    ZCLD = full_atm(j,1)
    CLD = m_cloud_liq(j)
    CLDICE = m_cloud_ice(j)
    RR = 0.
    
    write(50,'(4(F10.5))') &
      ZCLD,CLD,CLDICE,RR
  end do
  
  if (debug_flag_local) then
    do j=1,NCRSPC
      write(50,'(7(F10.5))') &
        wvl_cld(j),liq_ext_cld(j),liq_ssa_cld(j),liq_asym_cld(j), &
        ice_ext_cld(j),ice_ssa_cld(j),ice_asym_cld(j)
    end do
  endif
endif

!-------------
! Assign Card 2D
   IREG(1) = 42
   IREG(2) = 42
   IREG(3) = 42
   IREG(4) = 42

  if (a_flag_local) then
! Write Card 2D
   if (debug_flag_local) then
     write(50,'(4I5)') &
       IREG(1),IREG(2),IREG(3),IREG(4)
   endif

   !Cards 2D1 & 2D2
   do i=1,4

     !Card 2D1
     AWCCON = 0.1
     TITLE = '   AER-USERDEFINED'

     !Assign optical properties from o_aerosol
     do j=1,42
       VARSPC(j) = aer_wavelen(j)
       EXTC(j) = o_aerosol(j,1+(i-1)*3)
       ABSC(j) = o_aerosol(j,2+(i-1)*3)
       ASYM(j) = o_aerosol(j,3+(i-1)*3)
     end do
 
     if (debug_flag_local) then
     write(50,'(E10.3,A18)') &
       AWCCON,TITLE
     
     do j=1,42/3
         !Card 2D2
         write(50,'(3(F6.2,2F7.5,F6.4))') &
           VARSPC(1+(j-1)*3),EXTC(1+(j-1)*3),ABSC(1+(j-1)*3),ASYM(1+(j-1)*3), &
           VARSPC(2+(j-1)*3),EXTC(2+(j-1)*3),ABSC(2+(j-1)*3),ASYM(2+(j-1)*3), &
           VARSPC(3+(j-1)*3),EXTC(3+(j-1)*3),ABSC(3+(j-1)*3),ASYM(3+(j-1)*3)
       end do 
     endif
   end do
endif  
       
!Assign Card 3
   H1 = 700.
   H2 = surf_alt
   ANGLE = 180.
   mRANGE = 0.
   BETA = 0.
   RO = 0.
   LENN = 0
   PHI = 0.

!Write Card 3
   if (debug_flag_local) then
     write(50,'(6F10.4,I5,5X,F10.3)') &
       H1,H2,ANGLE,mRANGE,BETA,RO,LENN,PHI
   endif


!Assign Card 3A1
   IPARM = 2

   if (a_flag_local) then
     IPH = 1
   else
     IPH = 2
   endif
   IDAY = 93 !jday
   ISOURC = 0

  if (.not.lw_flag) then
!Write Card 3A1
   if (debug_flag_local) then
     write(50,'(4I5)') &
       IPARM,IPH,IDAY,ISOURC
   endif

!Assign Card 3A2
   PARM1 = az_real
   PARM2 = zen_real
   PARM3 = 0.
   PARM4 = 0.
   TIME = gmt
   PSIPO = 0.
   ANGLEM = 0.
   G = 0.

!Write Card 3A2
   if (debug_flag_local) then
     write(50,'(8F10.3)') &
       PARM1,PARM2,PARM3,PARM4,TIME,PSIPO,ANGLEM,G
   endif

!Assign Card 3B1
   NANGLS = 50
   NWLF =  15

if (a_flag_local) then

!Write Card 3B1
   if (debug_flag_local) then
     write(50,'(2I5)') &
       NANGLS,NWLF
   endif

!Assign Card 3C1
   do i=1,50
     ANGF(i) = phase_angles(i)
   end do

!Write Card 3C1
   if (debug_flag_local) then
   do i=1,6
     write(50,'(8(1X,F9.2))') &
       ANGF(1+(i-1)*8),ANGF(2+(i-1)*8),ANGF(3+(i-1)*8),ANGF(4+(i-1)*8), &
       ANGF(5+(i-1)*8),ANGF(6+(i-1)*8),ANGF(7+(i-1)*8),ANGF(8+(i-1)*8)
   end do
   write(50,'(2(1X,F9.2))') &
       ANGF(49),ANGF(50)
   endif

!Assign Card 3C2
   do i=1,15
     WLF(i) = phase_wvls(i)
   end do

!Write Card 3C2
   if (debug_flag_local) then
     write(50,'(8(1X,F9.3))') &
       WLF(1),WLF(2),WLF(3),WLF(4), &
       WLF(5),WLF(6),WLF(7),WLF(8)
     write(50,'(7(1X,F9.2))') &
       WLF(9),WLF(10),WLF(11),WLF(12), &
       WLF(13),WLF(14),WLF(15)
   endif

!Card 3C3-3C6
   do i=1,4
     do j=1,50
       do k=1,15
         F(k) = phase_funcs(j,k,i)
       end do
       if (debug_flag_local) then
       write(50,'(8(1X,E9.3))') &
         F(1),F(2),F(3),F(4), &
         F(5),F(6),F(7),F(8)
       write(50,'(7(1X,E9.3))') &
         F(9),F(10),F(11),F(12), &
         F(13),F(14),F(15)
       endif
     end do
   end do
endif !end of a_flag
endif !end of .not.lw_flag

!Assign Card 4
  IV1 = v1
  IV2 = v2
  IDV = dv
  !FWHM = dv*2.
  FWHM = dv
  YFLAG = 'R'
  DLIMIT = '  OUTPUT'
  if (lw_flag) then
    FLAGS = 'W6AA  T'
    XFLAG = 'W'
  else
    FLAGS = 'N6AA  T'
    XFLAG = 'N'
  endif
  MLFLX = num_levs

!Write Card 4
  if (debug_flag_local) then
    write(50,'(4F10.3,2A1,A8,A7,I3)') &
      IV1,IV2,IDV,FWHM,YFLAG,XFLAG,DLIMIT,FLAGS,MLFLX
  endif

!Assign Card 4A
  NSURF = 1.0
  AATEMP = -1.0

!Write Card 4A
  if (.not.lw_flag) then
    if (debug_flag_local) then
       write(50,'(I1,F9.0)') &
         NSURF,AATEMP
    endif
  endif 

  if (lw_flag) then
!    !Assing Card 4L1, 4L2, 4L3
!     NWVALB = spec_albedo_len
!     if (debug_flag_local) then
!      write(50,'(I2)') &
!        NWVALB
!     endif
!
!     do i=1,6
!       WVSURFALB = spec_albedo_wvl(i)
!       SURFALBVAL = 1.-spec_emissivity_vals(i)
!       if (debug_flag_local) then
!           write(50,'(2F10.6)') &
!             WVSURFALB,SURFALBVAL
!       endif
!     end do
  else
  !Assign Card 4B1
     CBRDF = '   Ross-Li'
  !Write Card 4B1
    if (debug_flag_local) then
      write(50,'(A10)') &
        CBRDF
    endif

  !Assign Card 4B2
     NWVSRF = brdf_len
     SURFZN = 0.
     SURFAZ = 0.

    !Write Card 4B2
    if (debug_flag_local) then
      write(50,'(I5,2F10.3)') &
        NWVSRF,SURFZN,SURFAZ
    endif

    do i=1,brdf_len
      !Assign Card 4B3
       WVSURF = brdf_wvl(i)
       PARAMS(1) = brdf_param(i,1)
       PARAMS(2) = brdf_param(i,2)
       PARAMS(3) = brdf_param(i,3)
       PARAMS(4) = 2.0
       PARAMS(5) = 1.0
       !Write Card 4B3
       if (debug_flag_local) then
         write(50,'(6F10.6)') &
           WVSURF,PARAMS(1),PARAMS(2),PARAMS(3),PARAMS(4),PARAMS(5)
       endif
     end do
  endif
!Assign Card 5
  IRPT = 0

!Write Card 5
  if (debug_flag_local) then
    write(50,'(I5)') &
      IRPT

!Close file
    close(unit=50)
   endif



    end subroutine modtran4_write_tp5

subroutine modtran5_write_tp5(surf_temp,co2vmr,surf_alt,jday,gmt,zen_real,az_real,a_flag,  &
                             full_atm,v1,v2,ocean_flag,brdf_len,brdf_wvl,brdf_param,    &
                             spec_albedo_len,spec_albedo_wvl,spec_emissivity_vals, &
                             c_flag,num_levs,m_aerosol,o_aerosol,m_cloud_liq,&
                             m_cloud_ice,wvl_cld,liq_ext_cld,liq_ssa_cld,liq_asym_cld, &
                             ice_ext_cld,ice_ssa_cld,ice_asym_cld,z_cld_levs,          &
                             aer_wavelen,fname,fname_length,phase_funcs,phase_wvls,    &
                             phase_angles,dv,tp5_file,debug_flag,lw_flag)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! write a tape5 file for analysis that can be used in Modtran5.3.0.0 externally
! 
!
! Author: D. Feldman U.C. Berkeley
! 
!-----------------------------------------------------------------------
!
! Input arguments
!

    real*4 :: surf_temp
    real*4 :: co2vmr
    real*4 :: surf_alt
    integer jday,num_levs
    real*4 :: gmt
    logical a_flag
    real*4 :: full_atm(num_levs,28)
    real*4 :: v1
    real*4 :: v2
    logical ocean_flag
    integer brdf_len
    real*4 :: brdf_wvl(15)
    real*4 :: brdf_param(15,3)
    integer spec_albedo_len
    real*4 :: spec_albedo_wvl(6)
    real*4 :: spec_emissivity_vals(6)
    logical c_flag
    real*4 :: m_aerosol(num_levs,4)
    real*4 :: o_aerosol(42,12)
    real*4 :: m_cloud_liq(num_levs)
    real*4 :: m_cloud_ice(num_levs)
    real*4 :: wvl_cld(24)
    real*4 :: liq_ext_cld(24)
    real*4 :: liq_ssa_cld(24)
    real*4 :: liq_asym_cld(24)
    real*4 :: ice_ext_cld(24)
    real*4 :: ice_ssa_cld(24)
    real*4 :: ice_asym_cld(24)
    real*4 :: z_cld_levs(num_levs-1)
    real*4 :: aer_wavelen(42)
    character*120 fname
    integer fname_length
    real*4 :: phase_funcs(50,15,4)
    real*4 :: phase_wvls(15)
    real*4 :: phase_angles(50)
    real*4:: dv
    real*4 :: zen_real
    real*4 :: az_real
    logical debug_flag
    logical lw_flag
    real*4  emis_broadband, alb_broadband
    
!   output argument

    character, intent(out):: tp5_file(1000)*110


!
!---------------------------Local variables-----------------------------
!
   integer i, j, k, jj, kk              ! index
   real*4 ii
   logical c_flag_local
   logical a_flag_local
   integer tp5_linenumber

!-----------------------------------------------------------------------
!Modtran tape5 variables
!

!Card 1
   character MODTRN,SPEED,BINARY,LYMOLC,CKPRNT
   integer MODEL,ITYPE,IEMSCT,IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT,I_RD2C
   real*4 TPTEMP
   character*7 SURREF

!Card 1A
   character DIS,DISAZM,DISALB
   integer NSTR
   integer ISUN
   real*4 CO2MX,SFWHM
   character*10 H2OSTR,O3STR
   character LSUNFL,LBMNAM,LFLTNM,H2OAER,LDATDR,C_PROF
   real*4 SOLCON
   character CDASTM,CDTDIR
   real*4 ASTMC,ASTMX,ASTMO,AERRH
   integer NSSALB

!Card 1A2
   character*7 BMNAME

!Card 2
   character*2 APLUS
   integer IHAZE
   character CNOVAM
   integer ISEASN
   character*3 ARUSS
   integer IVULCN,ICSTL,ICLD,IVSA
   real*4 VIS,WSS,WHH,RAINRT,GNDALT

!Optional Card 2A
   real*4 CTHIK,CALT,CEXT
   integer NCRALT,NCRSPC
   real*4 CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP

!Card 2C
   integer ML,IRD1,IRD2
   character*20 HMODEL
   real*4 REE
   integer NMOLYC

!Cards 2C1,2C2,2C2X
   real*4 ZM,P,T,WMOL(12),WMOLX(13)
   character*14 JCHAR
   character JCHARX,JCHARY  

!Card 2C3
   real*4 AHAZE(4),RRATZ

!Card 2D,2D1,2D2
   integer IREG(4)
   real*4 AWCCON
   character*18 TITLE
   real*4 VARSPC(42),EXTC(42),ABSC(42),ASYM(42) 

!Card 2E1,2E2
   real*4 ZCLD,CLD,CLDICE,RR
   real*4 WAVLEN,EXTC_LIQ,ABSC_LIQ,ASYM_LIQ
   real*4 EXTC_ICE,ABSC_ICE,ASYM_ICE

!Card 3
   real*4 H1ALT,H2ALT,OBSZEN,HRANGE,BETA,RAD_E
   integer LENN
   real*4 BCKZEN,CKRANG

!Cards 3A1,3A2
   integer IPARM,IPH,IDAY,ISOURC
   real*4 PARM1,PARM2,PARM3,PARM4,TIME,TRUEAZ,ANGLEM,G

! Card 3B1
   integer NANGLS,NWLF

!Card 3C1
   real*4 ANGF(50)

!Card 3C2
   real*4 WLF(15)

!Card 3C3-6
   real*4 F(15)

!Card 4
  real*4 IV1,IV2,IDV,FWHM
  character YFLAG,XFLAG
  character*8 DLIMIT
  character*7 FLAGS
  integer MLFLX

!Card 4A
   integer NSURF
   real*4 AATEMP
   real*4 DH2O
   character MLTRFL

!Card 4B1
   character*10 CBRDF

!Card 4B2
   integer NWVSRF
   real*4 SURFZN
   real*4 SURFAZ     

!Card 4B3
   real*4 WVSURF
   real*4 PARAMS(5)
 
   character*63 SEABRDF_LINE1 !for ross-sea brdf hard-coding
   character*63 SEABRDF_LINE2
   character*63 SEABRDF_LINE3

!DRF CHANGE
!Card 4L1
   character*20 SALBFL
   integer NWVALB
   
!Card 4L2
   character*20 CSALB
   real*4 WVSURFALB
   real*4 SURFALBVAL
!DRF CHANGE

!Card 5
   integer IRPT

   logical debug_flag_local

!-----------------------------------------------------------------------
!Open file

   debug_flag_local = debug_flag
   if (debug_flag_local) then
     write(*,*) 'fname debug write = ',fname(1:fname_length)
     open(unit=50,file=fname(1:fname_length)//'.tp5',status='new')
   endif
   c_flag_local = c_flag
   a_flag_local = a_flag

   if (lw_flag) then
      a_flag_local = .false.
   endif

   !c_flag_local = .false.
   !a_flag_local = .false.

!Assign Card 1
   if (lw_flag) then
     MODTRN = 'C'
   else
     MODTRN = 'M'
   endif
   SPEED = 'M'
   BINARY = 'F'
   LYMOLC = ' '
   MODEL = 7
   if (lw_flag) then
     ITYPE = 2
     IEMSCT = 1
     IMULT = -1
   else
     ITYPE = 2
     IEMSCT = 2
     IMULT = -1
   endif
   M1 = 0
   M2 = 0
   M3 = 0
   M4 = 0
   M5 = 0
   M6 = 0
   MDEF = 2
   I_RD2C = 1
   CKPRNT = 'F'
   IM = 1
   NOPRNT = -1
   TPTEMP = surf_temp
   if (lw_flag) then
     !convert emissivity from spectral to broadband
     !this follows the publication Bo-Hui Tang et al, Optics Express 19(1) 185-192 (2011)
     !"Estimation of broadband surface emissivity from narrowband emissivities"
     !emis_broadband = 0.0127+0.7852*spec_emissivity_vals(4) &
     !                 -0.0151*spec_emissivity_vals(5)+0.139*spec_emissivity_vals(6) 
     emis_broadband = 0.9999
     alb_broadband = 1.-emis_broadband
   else
     SURREF = 'BRDF'
   endif

!Write Card 1

   if (debug_flag_local) then
      if (lw_flag) then
        write(50,'(4A1,I1,11I5,A1,I4,F8.3,F7.5)') MODTRN,SPEED,BINARY,LYMOLC,MODEL, &
           ITYPE,IEMSCT,IMULT,M1,M2,M3,M4,M5,M6,MDEF,I_RD2C,CKPRNT,NOPRNT,TPTEMP,alb_broadband
      else
        write(50,'(4A1,I1,11I5,A1,I4,F8.3,A7)') MODTRN,SPEED,BINARY,LYMOLC,MODEL, &
           ITYPE,IEMSCT,IMULT,M1,M2,M3,M4,M5,M6,MDEF,I_RD2C,CKPRNT,NOPRNT,TPTEMP,SURREF
      endif
   endif
   tp5_linenumber=1
   if (lw_flag) then
     write(tp5_file(tp5_linenumber),'(4A1,I1,11I5,A1,I4,F8.3,F7.5)') MODTRN,SPEED,BINARY,LYMOLC,MODEL, &
       ITYPE,IEMSCT,IMULT,M1,M2,M3,M4,M5,M6,MDEF,I_RD2C,CKPRNT,NOPRNT,TPTEMP,alb_broadband
   else
     write(tp5_file(tp5_linenumber),'(4A1,I1,11I5,A1,I4,F8.3,A7)') MODTRN,SPEED,BINARY,LYMOLC,MODEL, &
       ITYPE,IEMSCT,IMULT,M1,M2,M3,M4,M5,M6,MDEF,I_RD2C,CKPRNT,NOPRNT,TPTEMP,SURREF
   endif
   tp5_linenumber=tp5_linenumber+1

!Assign Card 1A
   if (lw_flag) then 
     DIS = 'F'
     DISALB = 'F'
   else
     DIS = 'T'
     DISALB = 'T'
   endif
   DISAZM = 'F'
   NSTR = 8
   SFWHM = 1.
   CO2MX = co2vmr
   H2OSTR = '          '
   O3STR = '          '
   C_PROF = ' '
   LSUNFL = 'F'
   LBMNAM = 'T'
   LFLTNM = 'F'
   H2OAER = 'F'
   CDTDIR = 'F'
   if (lw_flag) then
      SOLCON = 0.0
   else
      SOLCON = 1360.0
   endif 
   CDASTM = ' '
   ASTMC = 0.
   ASTMX = 0.
   ASTMO = 0.
   AERRH = 0.
   NSSALB = 0
   
!Write Card 1A
   
   if (debug_flag_local) then
      write(50,'(3A1,I3,F4.2,F10.5,2A10,2A1,4(1X,A1),F10.3,A1,F9.3,3F10.2,I10)') &
        DIS,DISAZM,DISALB,NSTR,SFWHM,CO2MX,H2OSTR,O3STR,LSUNFL,C_PROF,LBMNAM, &
        LFLTNM,H2OAER,CDTDIR,SOLCON,CDASTM,ASTMC,ASTMX,ASTMO,AERRH,NSSALB
   endif
   write(tp5_file(tp5_linenumber),'(3A1,I3,F4.2,F10.5,2A10,2A1,4(1X,A1),F10.3,A1,F9.3,3F10.2,I10)') &
     DIS,DISAZM,DISALB,NSTR,SFWHM,CO2MX,H2OSTR,O3STR,LSUNFL,C_PROF,LBMNAM, &
     LFLTNM,H2OAER,CDTDIR,SOLCON,CDASTM,ASTMC,ASTMX,ASTMO,AERRH,NSSALB
   tp5_linenumber=tp5_linenumber+1

!Assign Card 1A2
   BMNAME = '15_2008'
   if (lw_flag) then
     BMNAME = '01_2008'
     !BMNAME = '15_2008'
   else
     BMNAME = '15_2008'
   endif

!Write Card 1A2
   if (debug_flag_local) then
      write(50,'(A7)') BMNAME
   endif
   write(tp5_file(tp5_linenumber),'(A7)') BMNAME
   tp5_linenumber=tp5_linenumber+1

!Assign Card 2
   APLUS = ' '
   CNOVAM = ' '
   ISEASN = 0
   if (a_flag_local) then
     ARUSS = 'USS'
     IHAZE = 7 
   else
     ARUSS = '   '
     IHAZE = 0
   endif
   IVULCN = 0
   !if (lw_flag) then
   !   ICSTL = 0
   !else
      ICSTL = 3
   !endif
   if (c_flag_local) then
      ICLD = 1
      IHAZE = -1
   else
      ICLD = 0
   endif
   IVSA = 0
   VIS = 0.
   WSS = 0.
   WHH = 0.
   RAINRT = 0. 
   GNDALT = surf_alt

!Write Card 2
   if (debug_flag_local) then
      write(50,'(A2,I3,A1,I4,A3,I2,3I5,5F10.4)')            &
        APLUS,IHAZE,CNOVAM,ISEASN,ARUSS,                            &
        IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT
   endif
   write(tp5_file(tp5_linenumber),'(A2,I3,A1,I4,A3,I2,3I5,5F10.4)')                    &
     APLUS,IHAZE,CNOVAM,ISEASN,ARUSS,                            &
     IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT
   tp5_linenumber=tp5_linenumber+1

!Optional Card 2A
   if (c_flag_local .and. ICLD.gt.0) then
     CTHIK = -9.
     CALT = -9.
     CEXT = -9.
     NCRALT = num_levs
     NCRSPC = 24
     CWAVLN = 0.55
     CCOLWD = -9.
     CCOLIP = -9.
     CHUMID = 110.
     ASYMWD = -9.
     ASYMIP = -9.

     if (debug_flag_local) then
        write(50,'(3F8.3,2I4,6F8.3)')  &
          CTHIK,CALT,CEXT,NCRALT,NCRSPC,CWAVLN,CCOLWD,CCOLIP,CHUMID, &
          ASYMWD,ASYMIP
     endif
     write(tp5_file(tp5_linenumber),'(3F8.3,2I4,6F8.3)')  &
       CTHIK,CALT,CEXT,NCRALT,NCRSPC,CWAVLN,CCOLWD,CCOLIP,CHUMID, &
       ASYMWD,ASYMIP
     tp5_linenumber=tp5_linenumber+1
   endif

!Assign Card 2C
   ML = num_levs
   IRD1 = 1

   if (a_flag_local .or. c_flag_local) then
     IRD2 = 2   
   else
     IRD2 = 0
   endif
   HMODEL = '   USERDEFINED'
   REE = 0.0
   NMOLYC = 0

!Write Card 2C
   if (debug_flag_local) then
      write(50,'(3I5,A20,F10.0,I5)')                               &
         ML,IRD1,IRD2,HMODEL,REE,NMOLYC
   endif
   write(tp5_file(tp5_linenumber),'(3I5,A20,F10.0,I5)')                               &
      ML,IRD1,IRD2,HMODEL,REE,NMOLYC
   tp5_linenumber=tp5_linenumber+1

!-------------

   !Assign Card 2C1
   if (lw_flag) then 
      JCHAR = 'AAAAAAAAAAAAAA'
   else
      JCHAR = 'AAAAAAAA1AAAAA'
   endif
   JCHARX = 'A'
   JCHARY = 'A'
   do i=1,num_levs

     !Assign Cards 2C1,2C2,2C2X
     ZM = full_atm(i,1)
     P = full_atm(i,2) 
     T = full_atm(i,3) 
     do j=1,12
       WMOL(j) = full_atm(i,j+3)
     end do
     do j=1,13
       WMOLX(j) = full_atm(i,j+15)
     end do
     AHAZE(1) = m_aerosol(i,1)
     RRATZ = 0.
     AHAZE(2) = m_aerosol(i,2)
     AHAZE(3) = m_aerosol(i,3)
     AHAZE(4) = m_aerosol(i,4)

     !Write Cards 2C1
     if (debug_flag_local) then
        write(50,'(F10.4,5E10.4,1A14,1X,2A1)') &
          ZM,P,T,WMOL(1),WMOL(2),WMOL(3),JCHAR,JCHARX,JCHARY
     endif
     write(tp5_file(tp5_linenumber),'(F10.4,5E10.4,1A14,1X,2A1)') &
       ZM,P,T,WMOL(1),WMOL(2),WMOL(3),JCHAR,JCHARX,JCHARY
     tp5_linenumber=tp5_linenumber+1

     !Write Card 2C2
     if (debug_flag_local) then
        write(50,'(8E10.4)') &
          WMOL(4),WMOL(5),WMOL(6),WMOL(7),WMOL(8),WMOL(9), &
          WMOL(10),WMOL(11)
        write(50,'(1F10.4)') WMOL(12)
     endif
     write(tp5_file(tp5_linenumber),'(8E10.4)') &
       WMOL(4),WMOL(5),WMOL(6),WMOL(7),WMOL(8),WMOL(9), &
       WMOL(10),WMOL(11)
     tp5_linenumber=tp5_linenumber+1
     write(tp5_file(tp5_linenumber),'(1F10.4)') WMOL(12)
     tp5_linenumber=tp5_linenumber+1

     !Write Card 2C2X
     if (debug_flag_local) then
        write(50,'(8E10.4)') &
          WMOLX(1),WMOLX(2),WMOLX(3),WMOLX(4),WMOLX(5),WMOLX(6),WMOLX(7),WMOLX(8)
        write(50,'(5E10.4)') &
          WMOLX(9), WMOLX(10),WMOLX(11),WMOLX(12),WMOLX(13)
     endif
     write(tp5_file(tp5_linenumber),'(8E10.4)') &
       WMOLX(1),WMOLX(2),WMOLX(3),WMOLX(4),WMOLX(5),WMOLX(6),WMOLX(7),WMOLX(8)
     tp5_linenumber=tp5_linenumber+1
     write(tp5_file(tp5_linenumber),'(5E10.4)') &
       WMOLX(9), WMOLX(10),WMOLX(11),WMOLX(12),WMOLX(13)
     tp5_linenumber=tp5_linenumber+1
     
     if (IRD2.eq.2) then
       !Write Card 2C3
       if (debug_flag_local) then
          write(50,'(10X,F10.8,10X,4F10.8)') & 
             AHAZE(1),RRATZ,AHAZE(2),AHAZE(3),AHAZE(4)
       endif
       write(tp5_file(tp5_linenumber),'(10X,F10.7,10X,4F10.7)') & 
          AHAZE(1),RRATZ,AHAZE(2),AHAZE(3),AHAZE(4)
       tp5_linenumber=tp5_linenumber+1
     endif
   end do

!-------------

!!! Apparently Card 2E1/2E2 read before 2D1/2D2

if (c_flag_local) then
  !Card 2E1, 2E2

  do j=1,num_levs
    ZCLD = full_atm(j,1)
    CLD = m_cloud_liq(j)
    CLDICE = m_cloud_ice(j)
    RR = 0.

    if (debug_flag_local) then
      write(50,'(4(F10.5))') &
        ZCLD,CLD,CLDICE,RR
    endif
    write(tp5_file(tp5_linenumber),'(4(F10.5))') &
      ZCLD,CLD,CLDICE,RR
    tp5_linenumber=tp5_linenumber+1
  end do
 
  do j=1,NCRSPC
    if (debug_flag_local) then
       write(50,'(7(F10.5))') &
         wvl_cld(j),liq_ext_cld(j),liq_ssa_cld(j),liq_asym_cld(j), &
         ice_ext_cld(j),ice_ssa_cld(j),ice_asym_cld(j)
    endif
    write(tp5_file(tp5_linenumber),'(7(F10.5))') &
      wvl_cld(j),liq_ext_cld(j),liq_ssa_cld(j),liq_asym_cld(j), &
      ice_ext_cld(j),ice_ssa_cld(j),ice_asym_cld(j)
    tp5_linenumber=tp5_linenumber+1
  end do 
endif

! Assign Card 2D
   IREG(1) = 42
   IREG(2) = 42
   IREG(3) = 42
   IREG(4) = 42

  if (a_flag_local) then
! Write Card 2D
   if (debug_flag_local) then
     write(50,'(4I5)') &
       IREG(1),IREG(2),IREG(3),IREG(4)
   endif
   write(tp5_file(tp5_linenumber),'(4I5)') &
     IREG(1),IREG(2),IREG(3),IREG(4)
   tp5_linenumber=tp5_linenumber+1

   !Cards 2D1 & 2D2
   do i=1,4

     !Card 2D1
     AWCCON = 0.1
     TITLE = '   AER-USERDEFINED'

     !Assign optical properties from o_aerosol
     do j=1,42
       VARSPC(j) = aer_wavelen(j)
       EXTC(j) = o_aerosol(j,1+(i-1)*3)
       ABSC(j) = o_aerosol(j,2+(i-1)*3)
       ASYM(j) = o_aerosol(j,3+(i-1)*3)
     end do
     
     if (debug_flag_local) then
        write(50,'(E10.3,A18)') &
           AWCCON,TITLE
     endif
     write(tp5_file(tp5_linenumber),'(E10.3,A18)') &
       AWCCON,TITLE
     tp5_linenumber=tp5_linenumber+1
 
     do j=1,42/3

       !Card 2D2
       if (debug_flag_local) then
          write(50,'(3(F6.2,2F7.5,F6.4))') &
            VARSPC(1+(j-1)*3),EXTC(1+(j-1)*3),ABSC(1+(j-1)*3),ASYM(1+(j-1)*3), &
            VARSPC(2+(j-1)*3),EXTC(2+(j-1)*3),ABSC(2+(j-1)*3),ASYM(2+(j-1)*3), &
            VARSPC(3+(j-1)*3),EXTC(3+(j-1)*3),ABSC(3+(j-1)*3),ASYM(3+(j-1)*3) 
       endif
       write(tp5_file(tp5_linenumber),'(3(F6.2,2F7.5,F6.4))') &
         VARSPC(1+(j-1)*3),EXTC(1+(j-1)*3),ABSC(1+(j-1)*3),ASYM(1+(j-1)*3), &
         VARSPC(2+(j-1)*3),EXTC(2+(j-1)*3),ABSC(2+(j-1)*3),ASYM(2+(j-1)*3), &
         VARSPC(3+(j-1)*3),EXTC(3+(j-1)*3),ABSC(3+(j-1)*3),ASYM(3+(j-1)*3) 
       tp5_linenumber=tp5_linenumber+1
     end do
   end do
endif


!Assign Card 3
   H1ALT = 700.
   if (lw_flag) then 
     H2ALT = surf_alt
   else
     H2ALT = surf_alt
   endif
   !H2ALT = 0.
   OBSZEN = 180.
   HRANGE = 0.
   BETA = 0.
   RAD_E = 0.
   LENN = 0
   BCKZEN = 0.
   CKRANG = 0.

!Write Card 3
   if (debug_flag_local) then
      write(50,'(6F10.4,I5,5X,2F10.3)') & 
        H1ALT,H2ALT,OBSZEN,HRANGE,BETA,RAD_E,LENN,BCKZEN,CKRANG
   endif
   write(tp5_file(tp5_linenumber),'(6F10.4,I5,5X,2F10.3)') & 
     H1ALT,H2ALT,OBSZEN,HRANGE,BETA,RAD_E,LENN,BCKZEN,CKRANG
   tp5_linenumber=tp5_linenumber+1

!Assign Card 3A1
   IPARM = 2 

   if (a_flag_local) then
     IPH = 1
   else
     IPH = 2
   endif
   IDAY = 93 !jday
   ISOURC = 0

!Write Card 3A1
   if (.not.lw_flag) then

     if (debug_flag_local) then
       write(50,'(4I5)') &
         IPARM,IPH,IDAY,ISOURC
     endif
     write(tp5_file(tp5_linenumber),'(4I5)') &
       IPARM,IPH,IDAY,ISOURC
     tp5_linenumber=tp5_linenumber+1

  !Assign Card 3A2
     PARM1 = az_real
     PARM2 = zen_real
     PARM3 = 0.
     PARM4 = 0.
     TIME = gmt
     TRUEAZ = 0.
     ANGLEM = 0.
     G = 0.

  !Write Card 3A2
     if (debug_flag_local) then
        write(50,'(8F10.3)') &
          PARM1,PARM2,PARM3,PARM4,TIME,TRUEAZ,ANGLEM,G
     endif
     write(tp5_file(tp5_linenumber),'(8F10.3)') &
       PARM1,PARM2,PARM3,PARM4,TIME,TRUEAZ,ANGLEM,G
     tp5_linenumber=tp5_linenumber+1
  !Assign Card 3B1
     NANGLS = 50
     NWLF =  15

  if (a_flag_local) then

    !Write Card 3B1
       if (debug_flag_local) then
          write(50,'(2I5)') &
             NANGLS,NWLF
       endif
       write(tp5_file(tp5_linenumber),'(2I5)') &
         NANGLS,NWLF
       tp5_linenumber=tp5_linenumber+1

    !Assign Card 3C1
       do i=1,50
         ANGF(i) = phase_angles(i) !in degrees
       end do

    !Write Card 3C1
       do i=1,6
         if (debug_flag_local) then
            write(50,'(8(1X,F9.2))') &
              ANGF(1+(i-1)*8),ANGF(2+(i-1)*8),ANGF(3+(i-1)*8),ANGF(4+(i-1)*8), &
              ANGF(5+(i-1)*8),ANGF(6+(i-1)*8),ANGF(7+(i-1)*8),ANGF(8+(i-1)*8) 
         endif
         write(tp5_file(tp5_linenumber),'(8(1X,F9.2))') &
           ANGF(1+(i-1)*8),ANGF(2+(i-1)*8),ANGF(3+(i-1)*8),ANGF(4+(i-1)*8), &
           ANGF(5+(i-1)*8),ANGF(6+(i-1)*8),ANGF(7+(i-1)*8),ANGF(8+(i-1)*8) 
         tp5_linenumber=tp5_linenumber+1
       end do
       if (debug_flag_local) then
          write(50,'(2(1X,F9.2))') &
             ANGF(49),ANGF(50)
       endif
       write(tp5_file(tp5_linenumber),'(2(1X,F9.2))') &
           ANGF(49),ANGF(50)
       tp5_linenumber=tp5_linenumber+1

  !Assign Card 3C2
     do i=1,15
       WLF(i) = phase_wvls(i)
     end do

  !Write Card 3C2
     if (debug_flag_local) then
        write(50,'(8(1X,F9.3))') &
          WLF(1),WLF(2),WLF(3),WLF(4), &
          WLF(5),WLF(6),WLF(7),WLF(8) 
        write(50,'(7(1X,F9.2))') &
          WLF(9),WLF(10),WLF(11),WLF(12), &
          WLF(13),WLF(14),WLF(15)
     endif
     write(tp5_file(tp5_linenumber),'(8(1X,F9.3))') &
       WLF(1),WLF(2),WLF(3),WLF(4), &
       WLF(5),WLF(6),WLF(7),WLF(8) 
     tp5_linenumber=tp5_linenumber+1
     write(tp5_file(tp5_linenumber),'(7(1X,F9.2))') &
       WLF(9),WLF(10),WLF(11),WLF(12), &
       WLF(13),WLF(14),WLF(15)
     tp5_linenumber=tp5_linenumber+1

  !Card 3C3-3C6
     do i=1,4
       do j=1,50
         do k=1,15
           F(k) = phase_funcs(j,k,i)
         end do
         if (debug_flag_local) then
            write(50,'(8(1X,E9.3))') &
              F(1),F(2),F(3),F(4), &
              F(5),F(6),F(7),F(8) 
            write(50,'(7(1X,E9.3))') &
              F(9),F(10),F(11),F(12), &
              F(13),F(14),F(15)
         endif
         write(tp5_file(tp5_linenumber),'(8(1X,E9.3))') &
           F(1),F(2),F(3),F(4), &
           F(5),F(6),F(7),F(8) 
         tp5_linenumber=tp5_linenumber+1
         write(tp5_file(tp5_linenumber),'(7(1X,E9.3))') &
           F(9),F(10),F(11),F(12), &
           F(13),F(14),F(15)
         tp5_linenumber=tp5_linenumber+1
       end do    
     end do
  endif !end of a_flag
endif !end of .not.lw_flag

!Assign Card 4
  IV1 = v1
  IV2 = v2
  IDV = dv
  FWHM = dv*2.
  YFLAG = 'R'
  DLIMIT = '        '
  if (lw_flag) then
    XFLAG = 'W'
    !FLAGS = ' 2    T'
    FLAGS = 'WR    T'
  else
    FLAGS = 'N6AA  T'
    XFLAG = 'N'
  endif
  MLFLX = num_levs-1

!Write Card 4
  if (debug_flag_local) then
    write(50,'(4F10.3,2A1,A8,A7,I3)') &
       IV1,IV2,IDV,FWHM,YFLAG,XFLAG,DLIMIT,FLAGS,MLFLX
  endif
  write(tp5_file(tp5_linenumber),'(4F10.3,2A1,A8,A7,I3)') &
    IV1,IV2,IDV,FWHM,YFLAG,XFLAG,DLIMIT,FLAGS,MLFLX
  tp5_linenumber=tp5_linenumber+1

!Assign Card 4A
  NSURF = 1
  AATEMP = -1.0
  DH2O = 0.
  MLTRFL = 'F'

!Write Card 4A
   if (.not.lw_flag) then
     if (debug_flag_local) then
       write(50,'(I1,2F9.2,A1)') &
          NSURF,AATEMP,DH2O,MLTRFL
     endif
     write(tp5_file(tp5_linenumber),'(I1,2F9.2,A1)') &
       NSURF,AATEMP,DH2O,MLTRFL
     tp5_linenumber=tp5_linenumber+1
   endif

   if (lw_flag) then
!    !Assing Card 4L1, 4L2, 4L3
!     SALBFL = 'DATA/spec_alb.dat'
!     if (debug_flag_local) then
!       write(50,'(A20)') &
!          SALBFL
!     endif
!     write(tp5_file(tp5_linenumber),'(A20)') &
!        SALBFL
!     tp5_linenumber=tp5_linenumber+1
!     CSALB = '    CCM3 sea ice'
!     if (debug_flag_local) then
!       write(50,'(A20)') &
!          CSALB
!     endif
!     write(tp5_file(tp5_linenumber),'(A20)') &
!        CSALB
!     tp5_linenumber=tp5_linenumber+1
!
!     do i=1,6
!       WVSURFALB = spec_albedo_wvl(i)
!       SURFALBVAL = 1.-spec_emissivity_vals(i)
!       if (debug_flag_local) then
!           write(50,'(2F10.6)') &
!             WVSURFALB,SURFALBVAL
!       endif
!       write(tp5_file(tp5_linenumber),'(2F10.6)') &
!             WVSURFALB,SURFALBVAL
!       tp5_linenumber=tp5_linenumber+1
!     end do
   else  !for if lw_flag
   !Assign Card 4B1
     if (ocean_flag) then
        CBRDF = 'ROSS-SEA  '
        !CBRDF = 'COX-MUNK  '
     else
        CBRDF = '   Ross-Li'  
     endif

     !Write Card 4B1
     if (debug_flag_local) then
        write(50,'(A10)') &
          CBRDF
     endif 
     write(tp5_file(tp5_linenumber),'(A10)') &
        CBRDF
     tp5_linenumber=tp5_linenumber+1

     !Assign Card 4B2
     if (ocean_flag) then
       NWVSRF = 3
       SURFZN = 0.
       SURFAZ = 0.

       !NWVSRF = 15
       !SURFZN = 0.
       !SURFAZ = 0.
     else
       NWVSRF = brdf_len
       SURFZN = 0.
       SURFAZ = 0.
     endif

     !Write Card 4B2
     if (debug_flag_local) then
        write(50,'(I5,2F10.3)') &
          NWVSRF,SURFZN,SURFAZ
     endif
     write(tp5_file(tp5_linenumber),'(I5,2F10.3)') &
       NWVSRF,SURFZN,SURFAZ
     tp5_linenumber=tp5_linenumber+1

     if (ocean_flag) then

       !do i=1,brdf_len
       !  !Assign Card 4B3
       !     WVSURF = brdf_wvl(i)
       !     PARAMS(1) = brdf_param(i,1)
       !     PARAMS(2) = brdf_param(i,2)
     !Write Card 4B3
       !  if (debug_flag_local) then
       !     write(50,'(6F10.6)') &
       !       WVSURF,PARAMS(1),PARAMS(2)
       !  endif
       !  write(tp5_file(tp5_linenumber),'(6F10.6)') &
       !    WVSURF,PARAMS(1),PARAMS(2)
       !  tp5_linenumber=tp5_linenumber+1
       !end do

       SEABRDF_LINE1  = '  0.30480 1.36984 0.000e+000 0.12570 0.11225    1.13446 0.0 0.0'
       SEABRDF_LINE2  = '  0.46030 1.34286 0.000e+000 0.12570 0.11225    1.13446 0.0 0.0'
       SEABRDF_LINE3  = '  0.92040 1.32322 1.000e-006 0.12570 0.11225    1.13446 0.0 0.0'
       if (debug_flag_local) then
          write(50,'(A63)') &
             SEABRDF_LINE1
          write(50,'(A63)') &
             SEABRDF_LINE2
          write(50,'(A63)') &
             SEABRDF_LINE3
       endif
       write(tp5_file(tp5_linenumber),'(A63)') &
          SEABRDF_LINE1
       tp5_linenumber=tp5_linenumber+1
       write(tp5_file(tp5_linenumber),'(A63)') &
          SEABRDF_LINE2
       tp5_linenumber=tp5_linenumber+1
       write(tp5_file(tp5_linenumber),'(A63)') &
          SEABRDF_LINE3
       tp5_linenumber=tp5_linenumber+1
     else
       do i=1,brdf_len
       !Assign Card 4B3
         WVSURF = brdf_wvl(i)
         PARAMS(1) = brdf_param(i,1)
         PARAMS(2) = brdf_param(i,2)
         PARAMS(3) = brdf_param(i,3)
         PARAMS(4) = 2.0
         PARAMS(5) = 1.0
      !Write Card 4B3
         if (debug_flag_local) then
            write(50,'(6F10.6)') &
              WVSURF,PARAMS(1),PARAMS(2),PARAMS(3),PARAMS(4),PARAMS(5)
         endif
         write(tp5_file(tp5_linenumber),'(6F10.6)') &
           WVSURF,PARAMS(1),PARAMS(2),PARAMS(3),PARAMS(4),PARAMS(5)
         tp5_linenumber=tp5_linenumber+1
       end do
     endif
  endif !of if lw_flag
!Assign Card 5
  IRPT = 0

!Write Card 5
   if (debug_flag_local) then
      write(50,'(I5)') IRPT
   endif
   write(tp5_file(tp5_linenumber),'(I5)') &
     IRPT
   tp5_linenumber=tp5_linenumber+1

!Close file
   if (debug_flag_local) then
      close(unit=50)
   endif
      return
    end subroutine modtran5_write_tp5


subroutine modtran_tape5_reader(surf_temp,co2vmr,surf_alt,gmt,zen_real,az_real,  &
                             full_atm,v1,v2,brdf_wvl,brdf_param,    &
                             spec_albedo_len,spec_albedo_wvl,spec_emissivity_vals, &
                             m_aerosol,o_aerosol,m_cloud_liq,&
                             m_cloud_ice,wvl_cld,liq_ext_cld,liq_ssa_cld,liq_asym_cld, &
                             ice_ext_cld,ice_ssa_cld,ice_asym_cld,z_cld_levs,          &
                             aer_wavelen,fname,fname_length,phase_funcs,phase_wvls,    &
                             phase_angles,dv,brdf_len,num_levs,a_flag,c_flag)

!---------
!Purpose: this routine will read in the .tp5 file that modtran_write_tp5 wrote

    character(len=120), intent(in) :: fname ! path to write to netcdf file (DRF)
    integer,intent(in) :: fname_length
    integer, intent(in) :: brdf_len
    integer, intent(in) :: num_levs
    real*4, intent(out) :: surf_temp
    real*4, intent(out) :: co2vmr
    real*4, intent(out) :: surf_alt
    real*4 , intent(out) :: gmt
    real*4 , intent(out) :: full_atm(26,28)
    real*4 , intent(out) :: v1
    real*4 , intent(out) :: v2
    real*4 , intent(out) :: brdf_wvl(15)
    real*4 , intent(out) :: brdf_param(15,3)
    integer, intent(out) :: spec_albedo_len
    real*4 , intent(out) :: spec_albedo_wvl(15)
    real*4 , intent(out) :: spec_emissivity_vals(15)
    real*4 , intent(out) :: m_aerosol(num_levs,4)
    real*4 , intent(out) :: o_aerosol(42,12)
    real*4 , intent(out) :: m_cloud_liq(num_levs-1)
    real*4 , intent(out) :: m_cloud_ice(num_levs-1)
    real*4 , intent(out) :: wvl_cld(24)
    real*4 , intent(out) :: liq_ext_cld(24)
    real*4 , intent(out) :: liq_ssa_cld(24)
    real*4 , intent(out) :: liq_asym_cld(24)
    real*4 , intent(out) :: ice_ext_cld(24)
    real*4 , intent(out) :: ice_ssa_cld(24)
    real*4 , intent(out) :: ice_asym_cld(24)
    real*4 , intent(out) :: z_cld_levs(num_levs-1)
    real*4 , intent(out) :: aer_wavelen(42)
    real*4 , intent(out) :: phase_funcs(50,15,4)
    real*4 , intent(out) :: phase_wvls(15)
    real*4 , intent(out) :: phase_angles(50)
    real*4 , intent(out) :: dv
    real*4 , intent(out) :: zen_real
    real*4 , intent(out) :: az_real
    logical, intent(in) :: a_flag
    logical, intent(in) :: c_flag
    
!
!---------------------------Local variables-----------------------------
!
   integer i, j, k, jj, kk              ! index
   real*4 ii
   logical c_flag_local
   logical a_flag_local


!Card 1
   character MODTRN,SPEED
   integer MODEL,ITYPE,IEMSCT,IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT
   real*4 TPTEMP
   character*7 SURREF

!Card 1A
   character DIS,DISAZM
   integer NSTR
   character LSUN
   integer ISUN
   real*4 CO2MX
   character*10 H2OSTR,O3STR
   character LSUNFL,LBMNAM,LFLTNM,H2OAER,LDATDR
   real*4 SOLCON

!Card 1A2
   character*17 BMNAME

!Card 2
   character*2 APLUS
   integer IHAZE
   character CNOVAM
   integer ISEASN
   character*3 ARUSS
   integer IVULCN,ICSTL,ICLD,IVSA
   real*4 VIS,WSS,WHH,RAINRT,GNDALT

!Optional Card 2A
   real*4 CTHIK,CALT,CEXT
   integer NCRALT,NCRSPC
   real*4 CWAVLN,CCOLWD,CCOLIP,CHUMID,ASYMWD,ASYMIP

!Card 2C
   integer ML,IRD1,IRD2
   character*20 HMODEL
   real*4 REE

!Cards 2C1,2C2,2C2X
   real*4 ZM,P,T,WMOL(12),WMOLX(13)
   character*14 JCHAR
   character JCHARX  

!Card 2C3
   real*4 AHAZE(4),RRATZ

!Card 2D,2D1,2D2
   integer IREG(4)
   real*4 AWCCON
   character*18 TITLE
   real*4 VARSPC(42),EXTC(42),ABSC(42),ASYM(42) 

!Card 2E1,2E2
   real*4 ZCLD,CLD,CLDICE,RR
   real*4 WAVLEN,EXTC_LIQ,ABSC_LIQ,ASYM_LIQ
   real*4 EXTC_ICE,ABSC_ICE,ASYM_ICE

!Card 3
   real*4 H1,H2,ANGLE,mRANGE,BETA,RO
   integer LENN
   real*4 PHI

!Cards 3A1,3A2
   integer IPARM,IPH,IDAY,ISOURC
   real*4 PARM1,PARM2,PARM3,PARM4,TIME,PSIPO,ANGLEM,G

! Card 3B1
   integer NANGLS,NWLF

!Card 3C1
   real*4 ANGF(50)

!Card 3C2
   real*4 WLF(15)

!Card 3C3-6
   real*4 F(15)

!Card 4
  real*4 IV1,IV2,IDV,FWHM
  character YFLAG,XFLAG
  character*8 DLIMIT
  character*7 FLAGS
  integer MLFLX

!Card 4A
   integer NSURF
   real*4 AATEMP

!DRF MODIFIED
!Card 4L1
   integer NWVALB
   
!Card 4L2
   real*4 WVSURFALB
   real*4 SURFALBVAL
!DRF MODIFIED


!Card 4B1
   character*10 CBRDF

!Card 4B2
   integer NWVSRF
   real*4 SURFZN
   real*4 SURFAZ     

!Card 4B3
   real*4 WVSURF
   real*4 PARAMS(5)

!Card 5
   integer IRPT

!-----------------------------------------------------------------------
!Open file
!-----------------------------------------------------------------------

    open(unit=50,file=fname(1:fname_length)//'.mod5.tp5',status='old')


!Read Card 1
    read(50,'(2A1,I3,12I5,F8.3,A7)') MODTRN,SPEED,MODEL,ITYPE,IEMSCT, &
        IMULT,M1,M2,M3,M4,M5,M6,MDEF,IM,NOPRNT,TPTEMP,SURREF
!Assign vals for Card 1  
    surf_temp = TPTEMP

!Read Card 1A
   read(50,'(2A1,I3,A1,I4,F10.5,2A10,4(1X,A1),(1X,A1),F10.3)') &
     DIS,DISAZM,NSTR,LSUN,ISUN,CO2MX,H2OSTR,O3STR,LSUNFL,LBMNAM, &
     LFLTNM,H2OAER,LDATDR,SOLCON
!Assign vals for Card 1A  
    co2vmr = CO2MX

!Read Card 1A2
   read(50,'(A17)') BMNAME

!Read Card 2
   read(50,'(A2,I3,A1,I4,A3,I2,3I5,5F10.4)')                    &
     APLUS,IHAZE,CNOVAM,ISEASN,ARUSS,                            &
     IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,RAINRT,GNDALT
   surf_alt = GNDALT

!Read Card 2C
   read(50,'(3I5,A20,F10.0)')                               &
      ML,IRD1,IRD2,HMODEL,REE

   do i=1,num_levs
     !Read Cards 2C1
     read(50,'(F10.4,5F10.4,1A14,1X,A1)') &
       ZM,P,T,WMOL(1),WMOL(2),WMOL(3),JCHAR,JCHARX     

     !Write Card 2C2
     read(50,'(8F10.4,/F10.4)') &
       WMOL(4),WMOL(5),WMOL(6),WMOL(7),WMOL(8),WMOL(9), &
       WMOL(10),WMOL(11),WMOL(12)

     !Write Card 2C2X
     read(50,'(8F10.4,/5F10.4)') &
       WMOLX(1),WMOLX(2),WMOLX(3),WMOLX(4),WMOLX(5),WMOLX(6),WMOLX(7),WMOLX(8),WMOLX(9), &
       WMOLX(10),WMOLX(11),WMOLX(12),WMOLX(13)

     !Assign values
     full_atm(i,1) = ZM
     full_atm(i,2) = P
     full_atm(i,3) = T
     do j=1,12
       full_atm(i,j+3) = WMOL(j)
     end do
     do j=1,13
       full_atm(i,j+15) = WMOLX(j)
     end do
     
     if (IRD2.eq.2) then
       !Read Card 2C3
       read(50,'(10X,F10.8,10X,4F10.8)') & 
          AHAZE(1),RRATZ,AHAZE(2),AHAZE(3),AHAZE(4)
       m_aerosol(i,1) = AHAZE(1)
       m_aerosol(i,2) = AHAZE(2)
       m_aerosol(i,3) = AHAZE(3)
       m_aerosol(i,4) = AHAZE(4)
     endif
   end do 

  if (a_flag) then
    ! Read Card 2D
       read(50,'(4I5)') &
          IREG(1),IREG(2),IREG(3),IREG(4)

   !Cards 2D1 & 2D2
   do i=1,4
     !Card 2D1
     read(50,'(E10.3,A18)') &
       AWCCON,TITLE

     do j=1,42/3
       !Card 2D2
       read(50,'(3(F6.2,2F7.5,F6.4))') &
         VARSPC(1+(j-1)*3),EXTC(1+(j-1)*3),ABSC(1+(j-1)*3),ASYM(1+(j-1)*3), &
         VARSPC(2+(j-1)*3),EXTC(2+(j-1)*3),ABSC(2+(j-1)*3),ASYM(2+(j-1)*3), &
         VARSPC(3+(j-1)*3),EXTC(3+(j-1)*3),ABSC(3+(j-1)*3),ASYM(3+(j-1)*3) 
     end do

     !Assign optical properties from o_aerosol
     do j=1,42
       aer_wavelen(j) = VARSPC(j)
       o_aerosol(j,1+(i-1)*3) = EXTC(j)
       o_aerosol(j,2+(i-1)*3) = ABSC(j)
       o_aerosol(j,3+(i-1)*3) = ASYM(j)
     end do
   end do
  endif

  !read Card 2E1,2E2
  do j=1,num_levs-1
    m_cloud_liq(j) = 0.
    m_cloud_ice(j) = 0.
    z_cld_levs(j) = 0.
  end do
  do j=1,24
    wvl_cld(j) = 0.
    liq_ext_cld(j) = 0.
    liq_ssa_cld(j) = 0.
    liq_asym_cld(j) = 0.
    ice_ext_cld(j) = 0.
    ice_ssa_cld(j) = 0.
    ice_asym_cld(j) = 0.
  end do

!Read Card 3
   read(50,'(6F10.4,I5,5X,F10.3)') & 
     H1,H2,ANGLE,mRANGE,BETA,RO,LENN,PHI

!Read Card 3A1
   read(50,'(4I5)') &
     IPARM,IPH,IDAY,ISOURC

!Read Card 3A2
   read(50,'(8F10.3)') &
     PARM1,PARM2,PARM3,PARM4,TIME,PSIPO,ANGLEM,G
   gmt = TIME
   az_real = PARM1
   zen_real = PARM2

if (a_flag) then

!Read Card 3B1
   read(50,'(2I5)') &
     NANGLS,NWLF


!Read Card 3C1
   do i=1,6
     read(50,'(8(1X,F9.2))') &
       ANGF(1+(i-1)*8),ANGF(2+(i-1)*8),ANGF(3+(i-1)*8),ANGF(4+(i-1)*8), &
       ANGF(5+(i-1)*8),ANGF(6+(i-1)*8),ANGF(7+(i-1)*8),ANGF(8+(i-1)*8) 
   end do
   read(50,'(2(1X,F9.2))') &
       ANGF(49),ANGF(50)

!Assign Card 3C1
   do i=1,50
     phase_angles(i) = ANGF(i)
   end do

!Read Card 3C2
   read(50,'(8(1X,F9.3))') &
     WLF(1),WLF(2),WLF(3),WLF(4), &
     WLF(5),WLF(6),WLF(7),WLF(8) 
   read(50,'(7(1X,F9.2))') &
     WLF(9),WLF(10),WLF(11),WLF(12), &
     WLF(13),WLF(14),WLF(15)

!Assign Card 3C2
   do i=1,15
     phase_wvls(i) = WLF(i)
   end do

!Card 3C3-3C6
   do i=1,4
     do j=1,50
       read(50,'(8(1X,E9.3))') &
         F(1),F(2),F(3),F(4), &
         F(5),F(6),F(7),F(8) 
       read(50,'(7(1X,E9.3))') &
         F(9),F(10),F(11),F(12), &
         F(13),F(14),F(15)
       do k=1,15
         phase_funcs(j,k,i) = F(k)
       end do
     end do    
   end do
endif !end of a_flag


!Read Card 4
  read(50,'(4F10.3,2A1,A8,A7,I3)') &
    IV1,IV2,IDV,FWHM,YFLAG,XFLAG,DLIMIT,FLAGS,MLFLX
  v1 = IV1
  v2 = IV2
  dv = IDV

!Read Card 4A
   read(50,'(I1,F9.0)') &
     NSURF,AATEMP

!Read Card 4L1 & 4L2 if appropriate
   if (SURREF.eq.'LAMBER') then
      read(50,'(I2)') &
          NWVALB
      spec_albedo_len = NWVALB
      do i=1,NWVALB
         read(50,'(2F10.6)') &
            WVSURFALB,SURFALBVAL
         spec_albedo_wvl(i) = WVSURFALB
         spec_emissivity_vals(i) = 1.-SURFALBVAL
      end do
   else
     !Read Card 4B1
      read(50,'(A10)') &
         CBRDF

     !Write Card 4B2
      read(50,'(I5,2F10.3)') &
         NWVSRF,SURFZN,SURFAZ

      do i=1,brdf_len
        !Read Card 4B3
        read(50,'(6F10.6)') &
          WVSURF,PARAMS(1),PARAMS(2),PARAMS(3),PARAMS(4),PARAMS(5)
      !Assign Card 4B3
        brdf_wvl(i) = WVSURF
        brdf_param(i,1) = PARAMS(1)
        brdf_param(i,2) = PARAMS(2)
        brdf_param(i,3) = PARAMS(3)
      end do
   endif


!Read Card 5
   read(50,'(I5)') &
     IRPT

    !close(unit=50)

      return

    end subroutine modtran_tape5_reader




subroutine modis_albedo_out(d_modis_asdir,d_modis_asdif,d_modis_aldir,d_modis_aldif, &
                            c_modis_asdir,c_modis_asdif,c_modis_aldir,c_modis_aldif)  

!---------------
!
! Purpose: to write a single spectrum output from Modtran to a netcdf file
!          differs from the regular modtran_write_tp5 in that this is output from the internal
!          call to Modtran
!
!
!  Author: D. Feldman, U.C. Berkekely, 2010
!

     use shr_kind_mod, only: r8 => shr_kind_r8
     use ioFileMod, only: getfil

     include 'netcdf.inc'

     real(r8), intent(in) :: d_modis_asdir(plat,plon)  !modis asdir from daac
     real(r8), intent(in) :: d_modis_asdif(plat,plon)  !modis asdif from daac
     real(r8), intent(in) :: d_modis_aldir(plat,plon)  !modis aldir from daac
     real(r8), intent(in) :: d_modis_aldif(plat,plon)  !modis aldif from daac

     real(r8), intent(in) :: c_modis_asdir(plat,plon)  !modis asdir calculated 
     real(r8), intent(in) :: c_modis_asdif(plat,plon)  !modis asdif calculated
     real(r8), intent(in) :: c_modis_aldir(plat,plon)  !modis aldir calculated
     real(r8), intent(in) :: c_modis_aldif(plat,plon)  !modis aldif calculated


     !local variables
     real(r8) bufferr(plon,plat)  !for writing to file
     integer :: ii,jj

     integer :: vid_DMODIS_ASDIR          ! Variable id for dmodis_asdir
     integer :: vid_DMODIS_ASDIF          ! Variable id for dmodis_asdif
     integer :: vid_DMODIS_ALDIR          ! Variable id for dmodis_aldir
     integer :: vid_DMODIS_ALDIF          ! Variable id for dmodis_aldif
     integer :: vid_CMODIS_ASDIR          ! Variable id for cmodis_asdir
     integer :: vid_CMODIS_ASDIF          ! Variable id for cmodis_asdif
     integer :: vid_CMODIS_ALDIR          ! Variable id for cmodis_aldir
     integer :: vid_CMODIS_ALDIF          ! Variable id for cmodis_aldif

!
! Dimension ids
!
     integer :: did_lat                       ! Latitude dimension ID
     integer :: did_lon                       ! Longitude dimension ID
     integer :: did_time                      ! Time dimension ID
!DRF

     integer :: vdims(2)                      ! dimension ids
     !integer :: vdims                      ! dimension ids
     integer :: nvdims                        ! number of dimensions
     integer :: ret                           ! return code
     integer :: itime                         ! integer for time dimension

!
! File ID
!
     integer :: nfid                          ! NetCDF id

!Dimension info
    integer :: start(2)                 ! starting point on edges
    integer :: count(2)                 ! number of elements to write
    !integer :: start                 ! starting point on edges
    !integer :: count                 ! number of elements to write
    integer :: plat_local
    integer :: plon_local
    character(len=21) :: fname ! path to write to netcdf file (DRF)

    fname = 'modis_alb_out_test.nc'

! create file
    itime = 1
    plat_local = 1
    plon_local = 1

    !write(*,*) 'plon = ',plon
    !write(*,*) 'plat = ',plat

    call wrap_create(fname, nf_clobber, nfid)
    call wrap_def_dim(nfid, "lat", plat, did_lat)
    call wrap_def_dim(nfid, "lon", plon, did_lon)
    call wrap_def_dim(nfid, "time", nf_unlimited, did_time)

    nvdims = 2
    vdims = (/did_lon, did_lat/)

    call wrap_def_var(nfid, "DMODIS_ASDIR",nf_double, nvdims, vdims, vid_DMODIS_ASDIR)
    call wrap_def_var(nfid, "DMODIS_ASDIF",nf_double, nvdims, vdims, vid_DMODIS_ASDIF)
    call wrap_def_var(nfid, "DMODIS_ALDIR",nf_double, nvdims, vdims, vid_DMODIS_ALDIR)
    call wrap_def_var(nfid, "DMODIS_ALDIF",nf_double, nvdims, vdims, vid_DMODIS_ALDIF)

    call wrap_def_var(nfid, "CMODIS_ASDIR",nf_double, nvdims, vdims, vid_CMODIS_ASDIR)
    call wrap_def_var(nfid, "CMODIS_ASDIF",nf_double, nvdims, vdims, vid_CMODIS_ASDIF)
    call wrap_def_var(nfid, "CMODIS_ALDIR",nf_double, nvdims, vdims, vid_CMODIS_ALDIR)
    call wrap_def_var(nfid, "CMODIS_ALDIF",nf_double, nvdims, vdims, vid_CMODIS_ALDIF)

    call wrap_close(nfid)

    call wrap_open(fname, nf_write, nfid)
    
    start = (/1,    1/)
    count = (/plon, plat/)

    do ii=1,plon
       do jj=1,plat
          bufferr(ii,jj) = d_modis_asdir(jj,ii)
       end do
    end do
    call wrap_put_vara_realx(nfid, vid_DMODIS_ASDIR, start, count, bufferr)

    do ii=1,plon
       do jj=1,plat
          bufferr(ii,jj) = d_modis_asdif(jj,ii)
       end do
    end do
    call wrap_put_vara_realx(nfid, vid_DMODIS_ASDIF, start, count, bufferr)

    do ii=1,plon
       do jj=1,plat
          bufferr(ii,jj) = d_modis_aldif(jj,ii)
       end do
    end do
    call wrap_put_vara_realx(nfid, vid_DMODIS_ALDIF, start, count, bufferr)

    do ii=1,plon
       do jj=1,plat
          bufferr(ii,jj) = d_modis_aldir(jj,ii)
       end do
    end do
    call wrap_put_vara_realx(nfid, vid_DMODIS_ALDIR, start, count, bufferr)

    do ii=1,plon
       do jj=1,plat
          bufferr(ii,jj) = c_modis_asdir(jj,ii)
       end do
    end do
    call wrap_put_vara_realx(nfid, vid_CMODIS_ASDIR, start, count, bufferr)

    do ii=1,plon
       do jj=1,plat
          bufferr(ii,jj) = c_modis_asdif(jj,ii)
       end do
    end do
    call wrap_put_vara_realx(nfid, vid_CMODIS_ASDIF, start, count, bufferr)

    do ii=1,plon
       do jj=1,plat
          bufferr(ii,jj) = c_modis_aldif(jj,ii)
       end do
    end do
    call wrap_put_vara_realx(nfid, vid_CMODIS_ALDIF, start, count, bufferr)

    do ii=1,plon
       do jj=1,plat
          bufferr(ii,jj) = c_modis_aldir(jj,ii)
       end do
    end do
    call wrap_put_vara_realx(nfid, vid_CMODIS_ALDIR, start, count, bufferr)

    call wrap_close(nfid)

    return

    end subroutine modis_albedo_out

subroutine write_rrtm_input(surf_temp,co2vmr,surf_alt,jday,gmt,zen_real,a_flag,  &
                             full_atm,spec_albedo_len,spec_albedo_wvl,spec_emissivity_vals, &
                             c_flag,num_levs,m_aerosol,o_aerosol,m_cloud_liq,&
                             m_cloud_ice,cld_frac,r_ice,r_liq, &
                             aer_wavelen,fname,fname_length,    &
                             debug_flag,lw_flag)




!----------------------------------------------------------------------- 
! 
! Purpose: 
! write a tape5 file for analysis that can be used in RRTM (LW or SW)
! 
!
! Author: D. Feldman Lawrence Berkeley National Laboratory
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
    use shr_kind_mod, only: r8 => shr_kind_r8
    use cloud_optics
    implicit none

    real*4,intent(in) :: surf_temp
    real*4 :: co2vmr
    real*4 :: surf_alt
    integer jday,num_levs
    real*4 :: gmt
    logical a_flag
    real*4 :: full_atm(num_levs,28)
    integer spec_albedo_len
    real*4 :: spec_albedo_wvl(6)
    real*4 :: spec_emissivity_vals(6)
    logical c_flag
    real*4 :: m_aerosol(num_levs,4)
    real*4 :: o_aerosol(42,12)
    real*4 :: m_cloud_liq(num_levs)
    real*4 :: m_cloud_ice(num_levs)
    real*4 :: cld_frac(num_levs)
    real(r8) :: r_ice(num_levs)
    real(r8) :: r_liq(num_levs)
    real*4 :: wvl_cld(24)
    real*4 :: liq_ext_cld(24)
    real*4 :: liq_ssa_cld(24)
    real*4 :: liq_asym_cld(24)
    real*4 :: ice_ext_cld(24)
    real*4 :: ice_ssa_cld(24)
    real*4 :: ice_asym_cld(24)
    real*4 :: aer_wavelen(42)
    character*120 fname
    integer fname_length
    real*4 :: phase_funcs(50,15,4)
    real*4 :: phase_wvls(15)
    real*4 :: phase_angles(50)
    real*4:: dv
    real*4 :: zen_real
    logical debug_flag
    logical lw_flag
    real*4  emis_broadband, alb_broadband
    
!   Local variables

    integer i,j,index
    real*4 diff_val

    ! Record 1.1
    character *80 CXID

    ! Record 1.2
    integer IATM,IXSECT,NUMANGS,IOUT,ICLD

    ! Record 1.4
    real*4 :: TBOUND
    integer IEMIS,IREFLECT
    real*4 :: SEMISS(16)
    real*4 :: wvl_emis(16)

    ! Record 3.1
    integer MODEL,IBMAX,NOPRNT,NMOL,IPUNCH,MUNITS,RE,CO2MX,REF_LAT

    ! Record 3.2
    real *4 :: HBOUND
    real *4 :: HTOA 

    ! Record 3.3b
    real *4 :: PBND(num_levs)
 
    ! Record 3.4
    integer IMMAX
    character *13 HMOD

    ! Record 3.5
    real *4 :: ZM(num_levs)
    real *4 :: PM(num_levs)
    real *4 :: TM(num_levs)
    character *1 JCHARP
    character *1 JCHART
    character *32 JCHAR
    character *1 JLONG

    ! Record 3.7
    character *3 ENDCHAR

    ! Record C1.1
    integer INFLAG,ICEFLAG,LIQFLAG
    character *1 TESTCHAR

    integer :: c_level_flag(num_levs)
    integer :: c_layer_number(num_levs)
    real*4  :: c_cwp(num_levs)
    real*4  :: c_fracice(num_levs)
    real*4  :: c_radius_ice(num_levs)
    real*4  :: c_radius_liq(num_levs)

    real(r8):: LAY(num_levs)
    real*4  :: CLDFRAC(num_levs)
    real(r8):: CWP(num_levs)
    real(r8):: FRACICE(num_levs)
    real(r8):: EFFRADICE(num_levs)
    real(r8):: EFFRADLIQ(num_levs)

!   Record 1.1
    CXID = 'RRTM input created by Dan Feldman'

!   Record 1.2
    IATM = 1
    IXSECT = 0
    NUMANGS = 4
    IOUT = 99

    if (c_flag) then
      ICLD = 2  !maximum-random overlap assumption
    else
      ICLD = 0
    endif

!   Record 1.4
    TBOUND  = surf_temp
    IEMIS   = 2
    IREFLECT = 0

    wvl_emis(1) = 514.29
    wvl_emis(2) = 24.29
    wvl_emis(3) = 17.94
    wvl_emis(4) = 15.08
    wvl_emis(5) = 13.24
    wvl_emis(6) = 11.20
    wvl_emis(7) = 9.73
    wvl_emis(8) = 8.87
    wvl_emis(9) = 7.83
    wvl_emis(10) = 6.98
    wvl_emis(11) = 6.16
    wvl_emis(12) = 5.18
    wvl_emis(13) = 4.62
    wvl_emis(14) = 4.32
    wvl_emis(15) = 4.02
    wvl_emis(16) = 3.46

    do i=1,16
       diff_val = 1000.
       index = 1
       do j=1,6
          if (abs(spec_albedo_wvl(j)-wvl_emis(i)).lt.diff_val) then
             diff_val = abs(spec_albedo_wvl(j)-wvl_emis(i))
             index = j
          endif
       end do
       SEMISS(i) = spec_emissivity_vals(index)
    end do  

!   Record 3.1
    MODEL = 0
    IBMAX = -1*num_levs
    NOPRNT= 0
    NMOL  = 32
    IPUNCH= 1
    MUNITS= 1
    RE    = 0
    CO2MX = 0
    REF_LAT = 0

!   Record 3.3b
    do i=1,num_levs
       PBND(i) = full_atm(i,2)
    end do
    HBOUND = full_atm(1,2)
    HTOA = full_atm(num_levs,2)

!   Record 3.4
    IMMAX = -1*num_levs
    HMOD  = ' user-defined'

!   Record 3.5
    do i=1,num_levs
       ZM(i) = full_atm(i,1)
       PM(i) = full_atm(i,2)
       TM(i) = full_atm(i,3)
       JCHARP = 'A'
       JCHART = 'A'
       JLONG = ' '
       JCHAR = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    end do

    ! Record 3.7
    ENDCHAR = '%%%'

    ! Record 3.1.1
    if (c_flag) then
       INFLAG = 2
       ICEFLAG = 3
       LIQFLAG = 1
       TESTCHAR = ' '
    endif

    do i=1,num_levs
       c_level_flag(i) = 0
       c_layer_number(i) = i
       if (m_cloud_liq(i).gt.0. .or. m_cloud_ice(i).gt.0.) then
          c_level_flag(i) = 1
          c_cwp(i) = m_cloud_liq(i) + m_cloud_ice(i)
          c_fracice(i) = m_cloud_ice(i)/(m_cloud_ice(i)+m_cloud_liq(i))
          c_radius_ice(i) = r_ice(i)
          c_radius_liq(i) = r_liq(i)
       endif
    end do

!   WRITE TAPE5_TS
    open(unit=50,file=fname(1:fname_length)//'.rrtm',status='new')
    open(unit=51,file=fname(1:fname_length)//'.rrtm_cld',status='new')
   
    ! Record 1.1 
    write(50,'(A80)') &
        CXID
  
    ! Record 1.2
    write(50,'(49X,I1,19X,I1,13X,I2,2X,I3,4X,I1)') &
       IATM,IXSECT,NUMANGS,IOUT,ICLD

    ! Record 1.4
    write(50,'(E10.3,1X,I1,2X,I1,16E5.3)') &
       TBOUND,IEMIS,IREFLECT,SEMISS(1),SEMISS(2),SEMISS(3),SEMISS(4), &
       SEMISS(5),SEMISS(6),SEMISS(7),SEMISS(8),SEMISS(9),SEMISS(10), &
       SEMISS(11),SEMISS(12),SEMISS(13),SEMISS(14),SEMISS(15),SEMISS(16)

    ! Record 3.1
    write(50,'(I5,5X,I5,5X,I5,I5,I5,3X,I2,F10.3,20X,F10.3)') &
       MODEL,IBMAX,NOPRNT,NMOL,IPUNCH,MUNITS,RE,CO2MX

    ! Record 3.2
    write(50,'(F10.3,F10.3)') &
       HBOUND,HTOA

    ! Record 3.3b
    write(50,'(8F10.3)') &
       PBND(1),PBND(2),PBND(3),PBND(4),PBND(5),PBND(7),PBND(8)

    write(50,'(8F10.3)') &
       PBND(9),PBND(10),PBND(11),PBND(12),PBND(13),PBND(14),PBND(15)

    write(50,'(8F10.3)') &
       PBND(16),PBND(17),PBND(18),PBND(19),PBND(20),PBND(21),PBND(22)

    write(50,'(8F10.3)') &
       PBND(22),PBND(23),PBND(24),PBND(25),PBND(26)

    ! Record 3.4
    write(50,'(I5,A24)') &
       IMMAX,HMOD

    ! Record 3.5 & 3.6.1 to 3.6.N
    do i=1,num_levs
       ! Record 3.5
       write(50,'(E10.3,E10.3,E10.3,5X,A1,A1,3X,A28)') &
          ZM(i),PM(i),TM(i),JCHARP,JCHART,JCHAR

       ! Record 3.6.1 to 3.6.N
       write(50,'(8E10.4)') &
          full_atm(i,4),full_atm(i,5),full_atm(i,6),full_atm(i,7),full_atm(i,8), &
          full_atm(i,9),full_atm(i,10),full_atm(i,11),full_atm(i,12)

       ! zero out all of the other radiatively active species
       do j=1,3
         write(50,'(8E10.4)') &
            full_atm(i,13),full_atm(i,13),full_atm(i,13),full_atm(i,13),full_atm(i,13), &
            full_atm(i,13),full_atm(i,13),full_atm(i,13),full_atm(i,13)
       end do
    end do

    ! Record 3.7
    write(50,'(A3)') &
       ENDCHAR

    if (c_flag) then
      ! Record C1.1
      write(51,'(4X,I1,4X,I1,4X,I1)') &
         INFLAG,ICEFLAG,LIQFLAG

      ! Record C1.2
      do i=1,num_levs
         if (c_level_flag(i).eq.1) then
           LAY = c_level_flag(i)
           CLDFRAC = cld_frac(i)
           CWP = c_cwp(i)
           FRACICE = c_fracice(i)
           EFFRADICE = r_ice(i)
           EFFRADLIQ = r_liq(i)

           write(51,'(A1,1X,I3,E10.5,E10.5,E10.5,E10.5,E10.5)') &
              TESTCHAR,LAY,CLDFRAC,CWP,FRACICE,EFFRADICE,EFFRADLIQ
         endif
      end do
      close(unit=51)
    endif 
    close(unit=50)

    end subroutine write_rrtm_input

subroutine run_modtran_or_pcrtm(lat_index,lon_index,nconfig_v,wgtv_v,ccon_v,totwgt_v,cliqwp,cicewp,cld,r_ice,r_liq, &
                        fname,fname_length,surftemp_real, co2vmr_real,gndalt_real, &
                        jday_input,gmt_real,zen_real,psipo_real,a_flag,c_flag_input,full_atm,pnm,v1_real,v2_real,ocean_flag, &
                        brdf_len,brdf_wvl,brdf_param2,spec_albedo_len,spec_albedo_wvl,spec_emissivity_vals, &
                        num_levs,m_aerosol,o_aerosol, aer_wavelen,pfs,phase_wvls_real,phase_angles_real, &
                        wvl_cld,ext_liq_cld,ssa_liq_cld,asym_liq_cld,ext_ice_cld,ssa_ice_cld,asym_ice_cld,z_cld_val_real, &
                        dv_input_real,num_lres,num_hres,num_ac_wvl,debug_flag,lw_flag, &           
                        wavelength_lres,radiance_lres_clr,radiance_lres_all,solar_flux, &
                        bb_updiffuse_clr,bb_dndiffuse_clr,bb_dndirect_clr,bb_updiffuse_all,bb_dndiffuse_all,bb_dndirect_all, &
                        flx_updiffuse_clr,flx_dndiffuse_clr,flx_dndirect_clr,flx_updiffuse_all,flx_dndiffuse_all, &
                        flx_dndirect_all,wavelength_hres,radiance_hres_clr,radiance_hres_all,solar_zenith, &
                        diffuse_flux_clr,diffuse_flux_all)
                        
!---------------
!
! Purpose: to manage the execution of Modtran from radctl.F90
!
!  Author: D. Feldman, U.C. Berkekely, 2010
!
                       
    use shr_kind_mod, only: r8 => shr_kind_r8
    use cloud_optics
    implicit none
#include <ptrrgrid.h>

! Input arguments
    integer, intent(in) :: lat_index        ! DRF
    integer, intent(in) :: lon_index        ! DRF
    integer, intent(in) :: nconfig_v        ! number of cloud configurations
    real(r8), intent(in) :: wgtv_v(15)      ! Fractional area for cloud configurations
    integer, intent(in) :: ccon_v(pver,15)  !Binary cloud configurations (0=no,1=yes)
    real(r8), intent(in) :: totwgt_v        ! Sum of wgtv_v = total fractional area
    real(r8), intent(in) :: cicewp(pver)    ! in-cloud cloud ice water path
    real(r8), intent(in) :: cliqwp(pver)    ! in-cloud cloud liquid water path
    real(r8), intent(in) :: cld(pver)       ! cloud-fraction
    real(r8), intent(in) :: r_ice(pver)     ! ice cloud effective radius
    real(r8), intent(in) :: r_liq(pver)     ! liquid cloud effective radius
    character(len=120), intent(in) :: fname !to be sent to modtran
    real(r8), intent (in) :: pnm(pverrp)   ! Model interface pressures (dynes/cm2) (DRF)

    integer, intent(in) :: fname_length !DRF
    real*4, intent(in)  :: surftemp_real                ! surface temp (for submission to driver) DRF
    real*4, intent(in)  :: co2vmr_real                ! input co2vmr (DRF)
    real*4, intent(in)  :: gndalt_real                ! ground altitude (for submission to driver) DRF
    integer, intent(in) :: jday_input            ! Julian day input
    real*4, intent(in)  :: gmt_real                ! gmt time (for submission to driver) DRF
    real*4, intent(in)  :: zen_real                ! zenith_angle (for submission to driver) DRF
    real*4, intent(in)  :: psipo_real                  ! solar azimuth angle that modtran requires (DRF)
    logical, intent(in) :: a_flag
    logical, intent(in) :: c_flag_input
    logical, intent(in) :: ocean_flag
    logical, intent(in) :: debug_flag
    logical, intent(in) :: lw_flag               !false if doing shortwave, true if doing longwave
    real*4, intent(in)  :: full_atm(pver,28)             !full-atmosphere of Z,P,T,H2O,CO2, etc ...
    real*4, intent(in)  :: v1_real                ! v1 in wavelength (for submission to driver) DRF
    real*4, intent(in)  :: v2_real                ! v2 in wavelengt (for submission to driver) DRF
    integer, intent(in) :: brdf_len !DRF
    real*4, intent(in)  :: brdf_wvl(15) !DRF
    real*4, intent(in)  :: brdf_param2(15,3) !For entry into driver.f !DRF
    integer, intent(in) :: spec_albedo_len !DRF
    real*4, intent(in)  :: spec_albedo_wvl(6)    !for LW emissivity DRF
    real*4, intent(in)  :: spec_emissivity_vals(6)
    integer, intent(in) :: num_levs                 !DRF
    real*4, intent(in)  :: m_aerosol(pver,4)             !profile of total aerosol mass (in km^-1 per layer)
    real*4, intent(in)  :: o_aerosol(42,12)           !absorption, extinction, ssa spectral data for aerosols
    real*4, intent(in)  :: aer_wavelen(42)              !DRF
    real*4, intent(in)  :: pfs(50,15,4) !DRF phase functions at appropriate relative humidity
    real*4, intent(in)  :: phase_wvls_real(15)        ! phase function wavelengths (for submission to driver) DRF
    real*4, intent(in)  :: phase_angles_real(50)        ! phase function angles (degrees) (for submission to driver) DRF
    real(r8),intent(in) :: wvl_cld(24)
    real(r8),intent(in) :: ext_liq_cld(24,pver)
    real(r8),intent(in) :: ssa_liq_cld(24,pver)
    real(r8),intent(in) :: asym_liq_cld(24,pver)
    real(r8),intent(in) :: ext_ice_cld(24,pver)
    real(r8),intent(in) :: ssa_ice_cld(24,pver)
    real(r8),intent(in) :: asym_ice_cld(24,pver)
    real*4, intent(in)  :: z_cld_val_real(pver-1)        ! cloud altitude values (for submission to driver) DRF
    real*4, intent(in)  :: dv_input_real               ! output spectral resolution (for submission to driver DRF
    integer, intent(in) :: num_hres
    integer, intent(in) :: num_lres        !DRF
    integer, intent(in) :: num_ac_wvl !DRF

! Output arguments
    real*4, intent(out) :: wavelength_lres(wvlng2)  ! Wavelength from Modtran (DRF)
    real*4, intent(out) :: wavelength_hres(wvlng_hres2)  ! Wavelength from Modtran (DRF)
    real*4, intent(out) :: radiance_lres_clr(wvlng2)   ! Radiance from Modtran (DRF)
    real*4, intent(out) :: radiance_lres_all(wvlng2)   ! Radiance from Modtran (DRF)
    real*4, intent(out) :: radiance_hres_clr(wvlng_hres2)   ! Radiance from Modtran (DRF)
    real*4, intent(out) :: radiance_hres_all(wvlng_hres2)   ! Radiance from Modtran (DRF)
    real*4, intent(out) :: solar_flux(wvlng2)   ! TOA solar flux from Modtran (DRF)
    real*4, intent(out) :: diffuse_flux_clr(wvlng2)   !SW diffuse flux from Modtran (DRF)
    real*4, intent(out) :: diffuse_flux_all(wvlng2)   !SW diffuse flux from Modtran (DRF)
    real*4, intent(out) :: bb_updiffuse_clr(pver)   !broadband upwelling diffuse flux from MODTRAN (DRF)
    real*4, intent(out) :: bb_updiffuse_all(pver)   !broadband upwelling diffuse flux from MODTRAN (DRF)
    real*4, intent(out) :: bb_dndiffuse_clr(pver)   !broadband downwelling diffuse flux from MODTRAN (DRF)
    real*4, intent(out) :: bb_dndiffuse_all(pver)   !broadband downwelling diffuse flux from MODTRAN (DRF)
    real*4, intent(out) :: bb_dndirect_clr(pver)   !broadband downwelling direct flux from MODTRAN (DRF)
    real*4, intent(out) :: bb_dndirect_all(pver)   !broadband downwelling direct flux from MODTRAN (DRF)
    real*4, intent(out) :: flx_updiffuse_clr(pver,num_bands)   !spectral upwelling diffuse flux from MODTRAN (DRF)
    real*4, intent(out) :: flx_updiffuse_all(pver,num_bands)   !spectral upwelling diffuse flux from MODTRAN (DRF)
    real*4, intent(out) :: flx_dndiffuse_clr(pver,num_bands)   !spectral downwelling diffuse flux from MODTRAN (DRF)
    real*4, intent(out) :: flx_dndiffuse_all(pver,num_bands)   !spectral downwelling diffuse flux from MODTRAN (DRF)
    real*4, intent(out) :: flx_dndirect_clr(pver,num_bands)   !spectral downwelling direct flux from MODTRAN (DRF)
    real*4, intent(out) :: flx_dndirect_all(pver,num_bands)   !spectral downwelling direct flux from MODTRAN (DRF)
    real*4, intent(out) :: solar_zenith   !Solar zenith angle from Modtran (DRF)

! Local variables
    real, allocatable :: tmp_radiance_lres(:,:)  !Temporary arrays for cloud-overlap approximation averaging
    real, allocatable :: tmp_radiance_hres(:,:)
    real, allocatable :: tmp_diffuse_flux(:,:)
    real, allocatable :: tmp_bb_updiffuse(:,:)
    real, allocatable :: tmp_bb_dndiffuse(:,:)
    real, allocatable :: tmp_bb_dndirect(:,:)
    real, allocatable :: tmp_flx_updiffuse(:,:,:)
    real, allocatable :: tmp_flx_dndiffuse(:,:,:)
    real, allocatable :: tmp_flx_dndirect(:,:,:)

    !For a single calculation of modtran
    real, allocatable :: rad_lres(:)
    real, allocatable :: rad_hres(:)
    real, allocatable :: diff_flux(:)
    real, allocatable :: bb_updiff(:)
    real, allocatable :: bb_dndiff(:)
    real, allocatable :: bb_dndir(:)
    real, allocatable :: flx_updiff(:,:)
    real, allocatable :: flx_dndiff(:,:)
    real, allocatable :: flx_dndir(:,:)

    integer num_profs !DRF Number of times to call cloud overlap maximum-random overlap profile generator
    integer i,jj,k,mm,kk
    real*4 :: m_cloud_liq_real(pver)              !profile of liquid water cloud DRF
    real*4 :: m_cloud_ice_real(pver)              !profile of liquid water cloud DRF
    real*4 :: cld_frac_real(pver)                 !cloud fraction
    real*4 min_cld_val                !min val for cloud determination
    real(r8) dummy_ext(24),dummy_ssa(24),dummy_asym(24),dummy_sum
    integer :: max_cld_index             !highest vertical layer with clouds present
    integer :: num_cld_layer             !number of cloud layers (for adding an optically thin cloud)
    logical c_flag

    character* (1) fval_1  !DRF for modtran temp file naming
    character* (2) fval_2  !ditto
    character* (3) fval_3  !ditto
    character* (1) filler_period !ditto
    character(len=120) :: fname_cld !to be sent to modtran
    integer :: fname_length_cld !DRF
    real*4 cext_input    !DRF 
    character tp5_file(1000)*110
    logical tp5_flag
    real*4 wvl_cld_real(24)
    real*4 liq_ext_avg_real(24),liq_ssa_avg_real(24),liq_asym_avg_real(24)
    real*4 ice_ext_avg_real(24),ice_ssa_avg_real(24),ice_asym_avg_real(24)
    real(r8) liq_ext_avg(24),liq_ssa_avg(24),liq_asym_avg(24)
    real(r8) ice_ext_avg(24),ice_ssa_avg(24),ice_asym_avg(24)
    character(len=23) :: fname_spectrum ! path to write to netcdf file (DRF)
    integer vis_nir_index  !separation index for vis/nir (700 nm)
    real*4 wvl_vis_nir     !separation value for vis/nir (700 nm)
    integer end_nir_index  !max index greater than 0 
    real*4 dummy_trapval   !variable used to calculate trapezoidal rule
    real*4 dummy_mcloud
    integer clr_flx_index !DRF
    !for surface reentrant problem in Modtran
    logical PASS1
    logical debug_flag_local
    COMMON /GTBRDFC/ PASS1

! begin regular code
    num_profs = 16
    min_cld_val = 1.e-5    !DRF modification
    clr_flx_index = 1
    debug_flag_local = debug_flag

    !if (lat_index.eq.56 .and. lon_index.eq.86) then
    !   debug_flag_local = .true.
    !endif

    do jj=1,wvlng2
       radiance_lres_all(jj) = 0.
       radiance_lres_clr(jj) = 0.
       diffuse_flux_all(jj) = 0.
       diffuse_flux_clr(jj) = 0.
    end do
    do jj=1,wvlng_hres2
       radiance_hres_all(jj) = 0.
       radiance_hres_clr(jj) = 0.
    end do
    do jj=1,pver
       bb_updiffuse_all(jj) = 0.
       bb_dndiffuse_all(jj) = 0.
       bb_dndirect_all(jj) = 0.
       bb_updiffuse_clr(jj) = 0.
       bb_dndiffuse_clr(jj) = 0.
       bb_dndirect_clr(jj) = 0.
    end do
    do jj=1,pver
       do i=1,num_bands
          flx_updiffuse_all(jj,i) = 0.
          flx_dndiffuse_all(jj,i) = 0.
          flx_dndirect_all(jj,i) = 0.
          flx_updiffuse_clr(jj,i) = 0.
          flx_dndiffuse_clr(jj,i) = 0.
          flx_dndirect_clr(jj,i) = 0.
       end do
    end do

    if (c_flag_input) then
        c_flag = .true.
       
        !allocate output variables
        allocate(tmp_radiance_lres(num_profs,wvlng2))
        allocate(tmp_diffuse_flux(num_profs,wvlng2))
        allocate(tmp_radiance_hres(num_profs,wvlng_hres2))
        allocate(tmp_bb_updiffuse(num_profs,pver))
        allocate(tmp_bb_dndiffuse(num_profs,pver))
        allocate(tmp_bb_dndirect(num_profs,pver))
        allocate(tmp_flx_updiffuse(pver,wvlng2,num_profs))
        allocate(tmp_flx_dndiffuse(pver,wvlng2,num_profs))
        allocate(tmp_flx_dndirect(pver,wvlng2,num_profs))

        allocate(rad_lres(wvlng2))
        allocate(rad_hres(wvlng_hres2))
        allocate(diff_flux(wvlng2))
        allocate(bb_updiff(pver))
        allocate(bb_dndiff(pver))
        allocate(bb_dndir(pver))
        allocate(flx_updiff(pver,wvlng2))
        allocate(flx_dndiff(pver,wvlng2))
        allocate(flx_dndir(pver,wvlng2))

        do jj=1,num_profs
            do i=1,wvlng2
               tmp_radiance_lres(jj,i) = 0.
               tmp_diffuse_flux(jj,i) = 0.
            end do
        end do
             
        do jj=1,num_profs
            do i=1,wvlng_hres2
                tmp_radiance_hres(jj,i) = 0.
            end do
        end do
        do jj=1,num_profs
            do i=1,pver
                tmp_bb_updiffuse(jj,i) = 0.
                tmp_bb_dndiffuse(jj,i) = 0.
                tmp_bb_dndirect(jj,i) = 0.
                do k=1,wvlng2
                   tmp_flx_updiffuse(i,k,jj) = 0.
                   tmp_flx_dndiffuse(i,k,jj) = 0.
                   tmp_flx_dndirect(i,k,jj) = 0.
                end do
            end do 
        end do
        
        do jj=1,nconfig_v
            do i=1,num_levs
                if (ccon_v(pver+1-i,jj).eq.1) then   !Implement cloud configuration
                    m_cloud_liq_real(i) = real(cliqwp(pver+1-i)/(7000.*dlog(pnm(pver+2-i)/pnm(1+pver-i))))
                    m_cloud_ice_real(i) = real(cicewp(pver+1-i)/(7000.*dlog(pnm(pver+2-i)/pnm(1+pver-i))))
                    !m_cloud_ice_real(i) = m_cloud_ice_real(i) - m_cloud_liq_real(i) !to deal with unusual definition of ICLDLWP in history file (ICLDLWP = LIQ + ICE)

                else
                    m_cloud_liq_real(i) = 0.
                    m_cloud_ice_real(i) = 0.
                endif
                cld_frac_real(i) = real(cld(pver+1-i))
            end do
            !write(*,*) 'm_cloud_liq in run_modtran = ', m_cloud_liq_real

            !Get cloud optical properties to deliver to modtran
            !!!
            !Weight the cloud optical properties by extinction
            do i=1,24 !reset vals by wavelength to 0
              dummy_ext(i) = 0.
              dummy_ssa(i) = 0.
              dummy_asym(i) = 0.
            end do
            do i=1,num_levs
               dummy_ext = dummy_ext+ext_liq_cld(:,i)*cliqwp(num_levs+1-i)
               dummy_ssa = dummy_ssa+ssa_liq_cld(:,i)*cliqwp(num_levs+1-i)
               do mm=1,24
                 dummy_asym(mm) = dummy_asym(mm)+asym_liq_cld(mm,i)*cliqwp(num_levs+1-i)*(1.+ssa_liq_cld(mm,i))
               end do
            end do

            do i=1,24
              dummy_sum = sum(cliqwp)
              wvl_cld_real(i) = real(wvl_cld(i))
              if (dummy_sum.gt.0.) then
                 !write(*,*) 'nonzero liq_ext'


                 liq_ext_avg(i) = dummy_ext(i)/sum(cliqwp)
                 liq_ext_avg_real(i) = real(liq_ext_avg(i))
                 liq_ssa_avg(i) = dummy_ssa(i)/sum(cliqwp)
                 liq_ssa_avg_real(i) = real(liq_ssa_avg(i))
                 if (liq_ssa_avg_real(i)==0.) then
                    liq_ssa_avg_real(i) = -1.0e-6
                 endif
                 dummy_sum = 0.
                 do mm=1,num_levs
                    dummy_sum = dummy_sum+cliqwp(mm)*(1.+ssa_liq_cld(i,mm))
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
               dummy_ext = dummy_ext+ext_ice_cld(:,i)*cicewp(num_levs+1-i)
               dummy_ssa = dummy_ssa+ssa_ice_cld(:,i)*cicewp(num_levs+1-i)
               do mm=1,24
                  dummy_asym(mm) = dummy_asym(mm)+asym_ice_cld(mm,i)*cicewp(num_levs+1-i)*(1.+ssa_ice_cld(mm,i))
               end do
            end do

            do i=1,24
              dummy_sum = sum(cicewp)
              if (dummy_sum.gt.0.) then
                ice_ext_avg(i) = dummy_ext(i)/sum(cicewp)
                ice_ext_avg_real(i) = real(ice_ext_avg(i))
                ice_ssa_avg(i) = dummy_ssa(i)/sum(cicewp)
                ice_ssa_avg_real(i) = real(ice_ssa_avg(i))
                if (ice_ssa_avg_real(i)==0.) then
                   ice_ssa_avg_real(i) = -1.0e-6
                endif
                dummy_sum = 0.
                do mm=1,num_levs
                   dummy_sum = dummy_sum+cicewp(mm)*(1.+ssa_ice_cld(i,mm))
                end do
                ice_asym_avg(i) = dummy_asym(i)/dummy_sum
                ice_asym_avg_real(i) = real(ice_asym_avg(i))
              endif
            end do

            !Set top two cloud layer to 0
             m_cloud_liq_real(pver) = 0.0
             m_cloud_ice_real(pver) = 0.0
             m_cloud_liq_real(pver-1) = 0.0
             m_cloud_ice_real(pver-1) = 0.0

             !Make sure that a single-layer cloud has a thin layer above to avoid modtran fatal error with single cloud layer
             if (sum(m_cloud_liq_real).gt.min_cld_val) then
                 !get max index of clouds
                 max_cld_index = 0.
                 num_cld_layer = 0.
                 do k=1,pver-2  !top model level does not have clouds
                     if (m_cloud_liq_real(k).gt.min_cld_val) then
                         max_cld_index = k
                         num_cld_layer = num_cld_layer + 1
                     endif
                 end do

                 !add a small layer if necessary
                 if (num_cld_layer.eq.1) then
                     !write(*,*) "adding thin layer"
                     if (max_cld_index.eq.pver-2) then
                         m_cloud_liq_real(pver-3)=min_cld_val
                     else
                         m_cloud_liq_real(max_cld_index+1)=min_cld_val
                     endif
                 endif
                 if (num_cld_layer.eq.0) then  !if all layers are too thin, then set to 0
                     do k=1,pver
                         m_cloud_liq_real(k) = 0.
                     end do
                 endif
             endif

             if (sum(m_cloud_ice_real).gt.min_cld_val) then
                 !get max index of clouds
                 max_cld_index = 0.
                 num_cld_layer = 0.
                 do k=1,pver-2  !top model level does not have clouds
                    if (m_cloud_ice_real(k).gt.min_cld_val) then
                        max_cld_index = k
                        num_cld_layer = num_cld_layer + 1
                    endif
                 end do
           
                 !add a small layer if necessary
                 if (num_cld_layer.eq.1) then
                     !write(*,*) "adding thin layer"
                     if (max_cld_index.eq.pver-2) then
                          m_cloud_ice_real(pver-3)=min_cld_val
                     else
                          m_cloud_ice_real(max_cld_index+1)=min_cld_val
                     endif
                 endif
                 if (num_cld_layer.eq.0) then  !if all layers are too thin, then set to 0
                     do k=1,pver
                        m_cloud_ice_real(k) = 0.
                     end do
                 endif
             endif

             !Change fname length
             filler_period = '.'
             if (jj.lt.10) then
                 write(fval_1,fmt='(i1)') jj
                 fname_cld = fname(1:fname_length)//filler_period//fval_1
                 fname_length_cld = fname_length + 2
             elseif (jj.ge.10 .and. jj.lt.100) then
                 write(fval_2,fmt='(i2)') jj
                 fname_cld = fname(1:fname_length)//filler_period //fval_2
                 fname_length_cld = fname_length + 3
             endif 

             !final determination on cloud values
             c_flag = .false.
             cext_input = sum(m_cloud_liq_real) 
             if (cext_input.gt.min_cld_val) then
                 c_flag = .true.
             endif

             !write(*,*) 'ccon_v = ',ccon_v(:,jj)
             !write(*,*) 'cliqwp = ',cliqwp
             !write(*,*) 'cicewp = ',cicewp
             !write(*,*) 'cext write_modtran line 2809 = ',cext_input

             cext_input = sum(m_cloud_ice_real) 
             if (cext_input.gt.min_cld_val) then
                 c_flag = .true.
             endif

             dummy_mcloud = sum(m_cloud_liq_real)

             ! call RRTM
             if (debug_flag_local) then
                 call write_rrtm_input(surftemp_real,co2vmr_real,gndalt_real,jday_input,gmt_real,zen_real,a_flag,  &
                             full_atm,spec_albedo_len,spec_albedo_wvl,spec_emissivity_vals, &
                             c_flag,num_levs,m_aerosol,o_aerosol,m_cloud_liq_real,&
                             m_cloud_ice_real,cld_frac_real,r_ice,r_liq, &
                             aer_wavelen,fname,fname_length,debug_flag,lw_flag)
             endif
             
             call modtran5_write_tp5(surftemp_real,co2vmr_real,gndalt_real,jday_input,gmt_real, &
                             zen_real,psipo_real,a_flag,  &
                             full_atm,v1_real,v2_real,ocean_flag,brdf_len,brdf_wvl,brdf_param2,    &
                             spec_albedo_len,spec_albedo_wvl,spec_emissivity_vals, &
                             c_flag,num_levs,m_aerosol,o_aerosol,m_cloud_liq_real, &
                             m_cloud_ice_real,wvl_cld_real,liq_ext_avg_real,liq_ssa_avg_real,liq_asym_avg_real, &
                             ice_ext_avg_real,ice_ssa_avg_real,ice_asym_avg_real,z_cld_val_real,aer_wavelen, &
                             fname_cld,fname_length_cld,pfs,phase_wvls_real,phase_angles_real,dv_input_real,tp5_file, &
                             debug_flag_local,lw_flag)
             tp5_flag = .true.

                PASS1 = .TRUE.
!                call init_modpcrtm(fname_cld,fname_length_cld,tp5_flag,tp5_file,num_levs,wavelength_lres,tmp_radiance_lres(jj,:), &
!                   wavelength_hres,tmp_radiance_hres(jj,:),num_lres,num_hres,num_ac_wvl, &
!                   solar_flux,tmp_diffuse_flux(jj,:), &
!                   tmp_bb_updiffuse(jj,:), tmp_bb_dndiffuse(jj,:),tmp_bb_dndirect(jj,:), &
!                   tmp_flx_updiffuse(:,:,jj),tmp_flx_dndiffuse(:,:,jj),tmp_flx_dndirect(:,:,jj))

!             call driver(fname_cld,fname_length_cld,tp5_flag,tp5_file,num_levs,wavelength_lres,tmp_radiance_lres(jj,:), &
!                         wavelength_hres,tmp_radiance_hres(jj,:),num_lres,num_hres,num_ac_wvl, &
!                         solar_flux,tmp_diffuse_flux(jj,:), &
!                         tmp_bb_updiffuse(jj,:), tmp_bb_dndiffuse(jj,:),tmp_bb_dndirect(jj,:), &
!                         tmp_flx_updiffuse(:,:,jj),tmp_flx_dndiffuse(:,:,jj),tmp_flx_dndirect(:,:,jj))

             !if (debug_flag_local) then
             !     write(*,*) 'radiance_hres run_modtran = ',tmp_radiance_hres(jj,:)
             !endif

    !if (lat_index.eq.46 .and. lon_index.eq.178) then
    !   write(*,*) 'tmp_bb_updiffuse = ',tmp_bb_updiffuse(jj,:)
    !   write(*,*) 'tmp_bb_dndiffuse = ',tmp_bb_dndiffuse(jj,:)
    !   write(*,*) 'tmp_bb_dndirect = ',tmp_bb_dndirect(jj,:)
    !endif

!             !call driver(fname_cld,fname_length_cld,tp5_flag,tp5_file,num_levs,wavelength_lres,rad_lres, &
!             !            wavelength_hres,rad_hres,num_lres,num_hres,num_ac_wvl, &
!             !            solar_flux,diff_flux,bb_updiff,bb_dndiff,bb_dndir, &
!             !            flx_updiff,flx_dndiff,flx_dndir)


             ! stop 'error write_modtran line 2851'


             !do kk=1,wvlng2
             !   tmp_radiance_lres(jj,kk) = rad_lres(kk)
             !end do
             !do kk=1,wvlng_hres2
             !   tmp_radiance_hres(jj,kk) = rad_hres(kk)
             !end do
             !do kk=1,pver
             !   tmp_bb_updiffuse(jj,kk) = bb_updiff(kk) 
             !   tmp_bb_dndiffuse(jj,kk) = bb_dndiff(kk) 
             !   tmp_bb_dndirect(jj,kk) = bb_dndir(kk) 
             !end do
 
             !do kk=1,pver
             !  do mm=1,wvlng2
             !     tmp_flx_updiffuse(kk,mm,jj) = flx_updiff(kk,mm)
             !     tmp_flx_dndiffuse(kk,mm,jj) = flx_dndiff(kk,mm)
             !     tmp_flx_dndirect(kk,mm,jj) = flx_dndir(kk,mm)
             !  end do
             !end do



             !if (debug_flag_local .and. c_flag) then
             !    call modtran4_write_tp5(surftemp_real,co2vmr_real,gndalt_real,jday_input,gmt_real, &
             !                   zen_real,psipo_real,a_flag,  &
             !                   full_atm,v1_real,v2_real,ocean_flag,brdf_len,brdf_wvl,brdf_param2,    &
             !                   spec_albedo_len,spec_albedo_wvl,spec_emissivity_vals, &
             !                   c_flag,num_levs,m_aerosol,o_aerosol,m_cloud_liq_real, &
             !                   m_cloud_ice_real,wvl_cld_real,liq_ext_avg_real,liq_ssa_avg_real,liq_asym_avg_real, &
             !                   ice_ext_avg_real,ice_ssa_avg_real,ice_asym_avg_real,z_cld_val_real,aer_wavelen, &
             !                   fname_cld,fname_length_cld,pfs,phase_wvls_real,phase_angles_real,dv_input_real,tp5_file, &
             !                   debug_flag_local,lw_flag)

             !    solar_zenith = zen_real
             !    fname_spectrum = 'modtran_spectrum_out.nc' 
             !    call modtran_out_netcdf(wavelength_lres,wavelength_hres,&
             !                            tmp_radiance_lres(jj,:), tmp_radiance_hres(jj,:),&
             !                            solar_flux, solar_zenith, tmp_diffuse_flux(jj,:),&
             !                            tmp_bb_updiffuse(jj,:),tmp_bb_dndiffuse(jj,:),tmp_bb_dndirect(jj,:), &
             !                            fname_spectrum)
             !    stop 'error write_modtran_tp5 line 2574'
             !endif

             !for trapezoidal rule to separate visible and near-ir
             wvl_vis_nir = 700.
             vis_nir_index = minloc(abs(wavelength_lres-wvl_vis_nir),1)    
             do kk=1,wvlng2
                if(wavelength_lres(kk).gt.10.0) then 
                    end_nir_index = kk
                endif
             end do
        end do !end of jj=1,nconfig_v

        do jj=1,nconfig_v
             radiance_lres_all = radiance_lres_all + tmp_radiance_lres(jj,:)*real(wgtv_v(jj)/totwgt_v)
             diffuse_flux_all = diffuse_flux_all + tmp_diffuse_flux(jj,:)*real(wgtv_v(jj)/totwgt_v)
             radiance_hres_all = radiance_hres_all+ tmp_radiance_hres(jj,:)*real(wgtv_v(jj)/totwgt_v)
             bb_updiffuse_all = bb_updiffuse_all+ tmp_bb_updiffuse(jj,:)*real(wgtv_v(jj)/totwgt_v)
             bb_dndiffuse_all = bb_dndiffuse_all+ tmp_bb_dndiffuse(jj,:)*real(wgtv_v(jj)/totwgt_v)
             bb_dndirect_all = bb_dndirect_all+ tmp_bb_dndirect(jj,:)*real(wgtv_v(jj)/totwgt_v)

             !Convert to just vis and nir using the trapezoidal rule
             wvl_vis_nir = 700.
             vis_nir_index = minloc(abs(wavelength_lres-wvl_vis_nir),1)    
             do kk=1,wvlng2
                if(wavelength_lres(kk).gt.10.0) then 
                   end_nir_index = kk
                endif
             end do

             do kk=1,num_levs
                dummy_trapval = tmp_flx_updiffuse(kk,1,jj)
                do k=2,vis_nir_index-1
                   dummy_trapval = dummy_trapval+2.*tmp_flx_updiffuse(kk,k,jj)
                end do
                dummy_trapval = dummy_trapval+tmp_flx_updiffuse(kk,vis_nir_index,jj)
                dummy_trapval = dummy_trapval*(wavelength_lres(vis_nir_index)-wavelength_lres(1))/(2.*vis_nir_index)
                flx_updiffuse_all(kk,1) = flx_updiffuse_all(kk,1)+dummy_trapval*real(wgtv_v(jj)/totwgt_v)

                dummy_trapval = tmp_flx_dndiffuse(kk,1,jj)
                do k=2,vis_nir_index-1
                   dummy_trapval = dummy_trapval+2.*tmp_flx_dndiffuse(kk,k,jj)
                end do
                dummy_trapval = dummy_trapval+tmp_flx_dndiffuse(kk,vis_nir_index,jj)
                dummy_trapval = dummy_trapval*(wavelength_lres(vis_nir_index)-wavelength_lres(1))/(2.*vis_nir_index) 
                flx_dndiffuse_all(kk,1) = flx_dndiffuse_all(kk,1)+dummy_trapval*real(wgtv_v(jj)/totwgt_v)

                dummy_trapval = tmp_flx_dndirect(kk,1,jj)
                do k=2,vis_nir_index-1
                   dummy_trapval = dummy_trapval+2.*tmp_flx_dndirect(kk,k,jj)
                end do
                dummy_trapval = dummy_trapval+tmp_flx_dndirect(kk,vis_nir_index,jj)
                dummy_trapval = dummy_trapval*(wavelength_lres(vis_nir_index)-wavelength_lres(1))/(2.*vis_nir_index) 
                flx_dndirect_all(kk,1) = flx_dndirect_all(kk,1)+dummy_trapval*real(wgtv_v(jj)/totwgt_v)

                !Now do near-infrared
                dummy_trapval = tmp_flx_updiffuse(kk,vis_nir_index,jj)
                do k=vis_nir_index+1,end_nir_index-1
                   dummy_trapval = dummy_trapval+2.*tmp_flx_updiffuse(kk,k,jj)
                end do
                dummy_trapval = dummy_trapval+tmp_flx_updiffuse(kk,end_nir_index,jj)
                dummy_trapval = dummy_trapval*  &
                         (wavelength_lres(end_nir_index)-wavelength_lres(vis_nir_index+1))/(2.*(end_nir_index-vis_nir_index+1))
                flx_updiffuse_all(kk,2) = flx_updiffuse_all(kk,2)+dummy_trapval*real(wgtv_v(jj)/totwgt_v)

                dummy_trapval = tmp_flx_dndiffuse(kk,vis_nir_index,jj)
                do k=vis_nir_index+1,end_nir_index-1
                   dummy_trapval = dummy_trapval+2.*tmp_flx_dndiffuse(kk,k,jj)
                end do
                dummy_trapval = dummy_trapval+tmp_flx_dndiffuse(kk,end_nir_index,jj)
                dummy_trapval = dummy_trapval*  &
                      (wavelength_lres(end_nir_index)-wavelength_lres(vis_nir_index+1))/(2.*(end_nir_index-vis_nir_index+1))
                flx_dndiffuse_all(kk,2) = flx_dndiffuse_all(kk,2)+dummy_trapval*real(wgtv_v(jj)/totwgt_v)

                dummy_trapval = tmp_flx_dndirect(kk,vis_nir_index,jj)
                do k=vis_nir_index+1,end_nir_index-1
                        dummy_trapval = dummy_trapval+2.*tmp_flx_dndirect(kk,k,jj)
                end do
                dummy_trapval = dummy_trapval+tmp_flx_dndirect(kk,end_nir_index,jj)
                dummy_trapval = dummy_trapval*  &
                    (wavelength_lres(end_nir_index)-wavelength_lres(vis_nir_index+1))/(2.*(end_nir_index-vis_nir_index+1))
                flx_dndirect_all(kk,2) = flx_dndirect_all(kk,2)+dummy_trapval*real(wgtv_v(jj)/totwgt_v)
             end do 
        end do

        !Also do clear-sky calculations
        c_flag = .false.
        call modtran5_write_tp5(surftemp_real,co2vmr_real,gndalt_real,jday_input,gmt_real, &
                             zen_real,psipo_real,a_flag,  &
                             full_atm,v1_real,v2_real,ocean_flag,brdf_len,brdf_wvl,brdf_param2,    &
                             spec_albedo_len,spec_albedo_wvl,spec_emissivity_vals, &
                             c_flag,num_levs,m_aerosol,o_aerosol,m_cloud_liq_real, &
                             m_cloud_ice_real,wvl_cld_real,liq_ext_avg_real,liq_ssa_avg_real,liq_asym_avg_real, &
                             ice_ext_avg_real,ice_ssa_avg_real,ice_asym_avg_real,z_cld_val_real,aer_wavelen, &
                             fname,fname_length,pfs,phase_wvls_real,phase_angles_real,dv_input_real,tp5_file,&
                             debug_flag_local,lw_flag)
        tp5_flag = .true.
        solar_zenith = zen_real
                PASS1 = .TRUE.

!        call init_modpcrtm(fname,fname_length,tp5_flag,tp5_file,num_levs,wavelength_lres,radiance_lres_clr, &
!                    wavelength_hres,radiance_hres_clr,num_lres,num_hres,num_ac_wvl, &
!                    solar_flux,diffuse_flux_clr, &
!                    bb_updiffuse_clr, bb_dndiffuse_clr,bb_dndirect_clr, &
!                    tmp_flx_updiffuse(:,:,clr_flx_index),tmp_flx_dndiffuse(:,:,clr_flx_index), &
!                    tmp_flx_dndirect(:,:,clr_flx_index))

!        call driver(fname,fname_length,tp5_flag,tp5_file,num_levs,wavelength_lres,radiance_lres_clr, &
!                    wavelength_hres,radiance_hres_clr,num_lres,num_hres,num_ac_wvl, &
!                    solar_flux,diffuse_flux_clr, &
!                    bb_updiffuse_clr, bb_dndiffuse_clr,bb_dndirect_clr, &
!                    tmp_flx_updiffuse(:,:,clr_flx_index),tmp_flx_dndiffuse(:,:,clr_flx_index), &
!                    tmp_flx_dndirect(:,:,clr_flx_index))

        !if (debug_flag_local) then
        !     call modtran4_write_tp5(surftemp_real,co2vmr_real,gndalt_real,jday_input,gmt_real, &
        !                     zen_real,psipo_real,a_flag,  &
        !                     full_atm,v1_real,v2_real,ocean_flag,brdf_len,brdf_wvl,brdf_param2,    &
        !                     spec_albedo_len,spec_albedo_wvl,spec_emissivity_vals, &
        !                     c_flag,num_levs,m_aerosol,o_aerosol,m_cloud_liq_real, &
        !                     m_cloud_ice_real,wvl_cld_real,liq_ext_avg_real,liq_ssa_avg_real,liq_asym_avg_real, &
        !                     ice_ext_avg_real,ice_ssa_avg_real,ice_asym_avg_real,z_cld_val_real,aer_wavelen, &
        !                     fname,fname_length,pfs,phase_wvls_real,phase_angles_real,dv_input_real,tp5_file,debug_flag_local,lw_flag)

        !     fname_spectrum = 'modtran_spectrum_out.nc' 
        !     call modtran_out_netcdf(wavelength_lres,wavelength_hres,&
        !                             radiance_lres_clr, radiance_hres_clr,&
        !                             solar_flux, solar_zenith, diffuse_flux_clr, &
        !                             bb_updiffuse_clr,bb_dndiffuse_clr,bb_dndirect_clr, &
        !                             fname_spectrum)

        !     stop 'error write_modtran_tp5 line 2694'
        !endif

        !for trapezoidal rule to separate visible and near-ir
        wvl_vis_nir = 700.
        vis_nir_index = minloc(abs(wavelength_lres-wvl_vis_nir),1)    
        do kk=1,wvlng2
             if(wavelength_lres(kk).gt.10.0) then 
                 end_nir_index = kk
             endif
        end do

        !write(*,*) "done with modtran"
        do kk=1,num_levs
             dummy_trapval = tmp_flx_updiffuse(kk,1,1)
             do k=2,vis_nir_index-1
                 dummy_trapval = dummy_trapval+2.*tmp_flx_updiffuse(kk,k,1)
             end do
             dummy_trapval = dummy_trapval+tmp_flx_updiffuse(kk,vis_nir_index,1)
             dummy_trapval = dummy_trapval*(wavelength_lres(vis_nir_index)-wavelength_lres(1))/(2.*vis_nir_index) 
             flx_updiffuse_clr(kk,1) = dummy_trapval

             dummy_trapval = tmp_flx_dndiffuse(kk,1,1)
             do k=2,vis_nir_index-1
                 dummy_trapval = dummy_trapval+2.*tmp_flx_dndiffuse(kk,k,1)
             end do
             dummy_trapval = dummy_trapval+tmp_flx_dndiffuse(kk,vis_nir_index,1)
             dummy_trapval = dummy_trapval*(wavelength_lres(vis_nir_index)-wavelength_lres(1))/(2.*vis_nir_index) 
             flx_dndiffuse_clr(kk,1) = dummy_trapval

             dummy_trapval = tmp_flx_dndirect(kk,1,1)
             do k=2,vis_nir_index-1
                 dummy_trapval = dummy_trapval+2.*tmp_flx_dndirect(kk,k,1)
             end do
             dummy_trapval = dummy_trapval+tmp_flx_dndirect(kk,vis_nir_index,1)
             dummy_trapval = dummy_trapval*(wavelength_lres(vis_nir_index)-wavelength_lres(1))/(2.*vis_nir_index) 
             flx_dndirect_clr(kk,1) = dummy_trapval

             dummy_trapval = tmp_flx_updiffuse(kk,vis_nir_index,1)
             do k=vis_nir_index+1,end_nir_index-1
                 dummy_trapval = dummy_trapval+2.*tmp_flx_updiffuse(kk,k,1)
             end do
             dummy_trapval = dummy_trapval+tmp_flx_updiffuse(kk,end_nir_index,1)
             dummy_trapval = dummy_trapval*  &
                      (wavelength_lres(end_nir_index)-wavelength_lres(vis_nir_index+1))/(2.*(end_nir_index-vis_nir_index+1))
             flx_updiffuse_clr(kk,2) = dummy_trapval

             dummy_trapval = tmp_flx_dndiffuse(kk,vis_nir_index,1)
             do k=vis_nir_index+1,end_nir_index-1
                 dummy_trapval = dummy_trapval+2.*tmp_flx_dndiffuse(kk,k,1)
             end do
             dummy_trapval = dummy_trapval+tmp_flx_dndiffuse(kk,end_nir_index,1)
             dummy_trapval = dummy_trapval*  &
                      (wavelength_lres(end_nir_index)-wavelength_lres(vis_nir_index+1))/(2.*(end_nir_index-vis_nir_index+1))
             flx_dndiffuse_clr(kk,2) = dummy_trapval

             dummy_trapval = tmp_flx_dndirect(kk,vis_nir_index,1)
             do k=vis_nir_index+1,end_nir_index-1
                 dummy_trapval = dummy_trapval+2.*tmp_flx_dndirect(kk,k,1)
             end do
             dummy_trapval = dummy_trapval+tmp_flx_dndirect(kk,end_nir_index,1)
             dummy_trapval = dummy_trapval*  &
                      (wavelength_lres(end_nir_index)-wavelength_lres(vis_nir_index+1))/(2.*(end_nir_index-vis_nir_index+1))
             flx_dndirect_clr(kk,2) = dummy_trapval
        end do

        !Deallocate arrays
        deallocate(tmp_radiance_lres)
        deallocate(tmp_radiance_hres)
        deallocate(tmp_diffuse_flux)
        deallocate(tmp_bb_updiffuse)
        deallocate(tmp_bb_dndiffuse)
        deallocate(tmp_bb_dndirect)
        deallocate(tmp_flx_updiffuse)
        deallocate(tmp_flx_dndiffuse)
        deallocate(tmp_flx_dndirect)


        deallocate(rad_lres)
        deallocate(rad_hres)
        deallocate(diff_flux)
        deallocate(bb_updiff)
        deallocate(bb_dndiff)
        deallocate(bb_dndir)
        deallocate(flx_updiff)
        deallocate(flx_dndiff)
        deallocate(flx_dndir)


    else
        !cloud-free conditions
        allocate(tmp_flx_updiffuse(pver,wvlng2,num_profs))
        allocate(tmp_flx_dndiffuse(pver,wvlng2,num_profs))
        allocate(tmp_flx_dndirect(pver,wvlng2,num_profs))

        do jj=1,num_profs
           do i=1,pver
              do k=1,wvlng2
                 tmp_flx_updiffuse(i,k,jj) = 0.
                 tmp_flx_dndiffuse(i,k,jj) = 0.
                 tmp_flx_dndirect(i,k,jj) = 0.
              end do
              m_cloud_liq_real(i) = 0.
              m_cloud_ice_real(i) = 0.
           end do 
        end do

        c_flag = .false.
        call modtran5_write_tp5(surftemp_real,co2vmr_real,gndalt_real,jday_input,gmt_real, &
                             zen_real,psipo_real,a_flag,  &
                             full_atm,v1_real,v2_real,ocean_flag,brdf_len,brdf_wvl,brdf_param2,    &
                             spec_albedo_len,spec_albedo_wvl,spec_emissivity_vals, &
                             c_flag,num_levs,m_aerosol,o_aerosol,m_cloud_liq_real, &
                             m_cloud_ice_real,wvl_cld_real,liq_ext_avg_real,liq_ssa_avg_real,liq_asym_avg_real, &
                             ice_ext_avg_real,ice_ssa_avg_real,ice_asym_avg_real,z_cld_val_real,aer_wavelen, &
                             fname,fname_length,pfs,phase_wvls_real,phase_angles_real,dv_input_real,tp5_file,&
                             debug_flag_local,lw_flag)
        tp5_flag = .true.
        solar_zenith = zen_real
                PASS1 = .TRUE.

!        call init_modpcrtm(fname,fname_length,tp5_flag,tp5_file,num_levs,wavelength_lres,radiance_lres_clr, &
!                    wavelength_hres,radiance_hres_clr,num_lres,num_hres,num_ac_wvl, &
!                    solar_flux,diffuse_flux_clr, &
!                    bb_updiffuse_clr, bb_dndiffuse_clr,bb_dndirect_clr, &
!                    tmp_flx_updiffuse(:,:,clr_flx_index),tmp_flx_dndiffuse(:,:,clr_flx_index), &
!                    tmp_flx_dndirect(:,:,clr_flx_index))

!        call driver(fname,fname_length,tp5_flag,tp5_file,num_levs,wavelength_lres,radiance_lres_clr, &
!                    wavelength_hres,radiance_hres_clr,num_lres,num_hres,num_ac_wvl, &
!                    solar_flux,diffuse_flux_clr, &
!                    bb_updiffuse_clr, bb_dndiffuse_clr,bb_dndirect_clr, &
!                    tmp_flx_updiffuse(:,:,clr_flx_index),tmp_flx_dndiffuse(:,:,clr_flx_index), &
!                    tmp_flx_dndirect(:,:,clr_flx_index))

        !if (debug_flag_local) then
        !     call modtran4_write_tp5(surftemp_real,co2vmr_real,gndalt_real,jday_input,gmt_real, &
        !                     zen_real,psipo_real,a_flag,  &
        !                     full_atm,v1_real,v2_real,ocean_flag,brdf_len,brdf_wvl,brdf_param2,    &
        !                     spec_albedo_len,spec_albedo_wvl,spec_emissivity_vals, &
        !                     c_flag,num_levs,m_aerosol,o_aerosol,m_cloud_liq_real, &
        !                     m_cloud_ice_real,wvl_cld_real,liq_ext_avg_real,liq_ssa_avg_real,liq_asym_avg_real, &
        !                     ice_ext_avg_real,ice_ssa_avg_real,ice_asym_avg_real,z_cld_val_real,aer_wavelen, &
        !                     fname,fname_length,pfs,phase_wvls_real,phase_angles_real,dv_input_real,tp5_file,debug_flag_local,lw_flag)

        !     fname_spectrum = 'modtran_spectrum_out.nc' 
        !     call modtran_out_netcdf(wavelength_lres,wavelength_hres,&
        !                             radiance_lres_clr, radiance_hres_clr,&
        !                             solar_flux, solar_zenith, diffuse_flux_clr, &
        !                             bb_updiffuse_clr,bb_dndiffuse_clr,bb_dndirect_clr, &
        !                             fname_spectrum)

        !     stop 'error write_modtran_tp5 line 2822'
        !endif

        !for trapezoidal rule to separate visible and near-ir
        wvl_vis_nir = 700.
        vis_nir_index = minloc(abs(wavelength_lres-wvl_vis_nir),1)    
        do kk=1,wvlng2
             if(wavelength_lres(kk).gt.10.0) then 
                 end_nir_index = kk
             endif
        end do

        !write(*,*) "done with modtran"
        do kk=1,num_levs
             dummy_trapval = tmp_flx_updiffuse(kk,1,1)
             do k=2,vis_nir_index-1
                 dummy_trapval = dummy_trapval+2.*tmp_flx_updiffuse(kk,k,1)
             end do
             dummy_trapval = dummy_trapval+tmp_flx_updiffuse(kk,vis_nir_index,1)
             dummy_trapval = dummy_trapval*(wavelength_lres(vis_nir_index)-wavelength_lres(1))/(2.*vis_nir_index) 
             flx_updiffuse_clr(kk,1) = dummy_trapval

             dummy_trapval = tmp_flx_dndiffuse(kk,1,1)
             do k=2,vis_nir_index-1
                 dummy_trapval = dummy_trapval+2.*tmp_flx_dndiffuse(kk,k,1)
             end do
             dummy_trapval = dummy_trapval+tmp_flx_dndiffuse(kk,vis_nir_index,1)
             dummy_trapval = dummy_trapval*(wavelength_lres(vis_nir_index)-wavelength_lres(1))/(2.*vis_nir_index) 
             flx_dndiffuse_clr(kk,1) = dummy_trapval

             dummy_trapval = tmp_flx_dndirect(kk,1,1)
             do k=2,vis_nir_index-1
                 dummy_trapval = dummy_trapval+2.*tmp_flx_dndirect(kk,k,1)
             end do
             dummy_trapval = dummy_trapval+tmp_flx_dndirect(kk,vis_nir_index,1)
             dummy_trapval = dummy_trapval*(wavelength_lres(vis_nir_index)-wavelength_lres(1))/(2.*vis_nir_index) 
             flx_dndirect_clr(kk,1) = dummy_trapval

             dummy_trapval = tmp_flx_updiffuse(kk,vis_nir_index,1)
             do k=vis_nir_index+1,end_nir_index-1
                 dummy_trapval = dummy_trapval+2.*tmp_flx_updiffuse(kk,k,1)
             end do
             dummy_trapval = dummy_trapval+tmp_flx_updiffuse(kk,end_nir_index,1)
             dummy_trapval = dummy_trapval*  &
                      (wavelength_lres(end_nir_index)-wavelength_lres(vis_nir_index+1))/(2.*(end_nir_index-vis_nir_index+1))
             flx_updiffuse_clr(kk,2) = dummy_trapval

             dummy_trapval = tmp_flx_dndiffuse(kk,vis_nir_index,1)
             do k=vis_nir_index+1,end_nir_index-1
                 dummy_trapval = dummy_trapval+2.*tmp_flx_dndiffuse(kk,k,1)
             end do
             dummy_trapval = dummy_trapval+tmp_flx_dndiffuse(kk,end_nir_index,1)
             dummy_trapval = dummy_trapval*  &
                      (wavelength_lres(end_nir_index)-wavelength_lres(vis_nir_index+1))/(2.*(end_nir_index-vis_nir_index+1))
             flx_dndiffuse_clr(kk,2) = dummy_trapval

             dummy_trapval = tmp_flx_dndirect(kk,vis_nir_index,1)
             do k=vis_nir_index+1,end_nir_index-1
                 dummy_trapval = dummy_trapval+2.*tmp_flx_dndirect(kk,k,1)
             end do
             dummy_trapval = dummy_trapval+tmp_flx_dndirect(kk,end_nir_index,1)
             dummy_trapval = dummy_trapval*  &
                      (wavelength_lres(end_nir_index)-wavelength_lres(vis_nir_index+1))/(2.*(end_nir_index-vis_nir_index+1))
             flx_dndirect_clr(kk,2) = dummy_trapval
        end do

        !Assign clear-sky values to all-sky values
        radiance_lres_all = radiance_lres_clr
        radiance_hres_all = radiance_hres_clr
        diffuse_flux_all = diffuse_flux_clr
        bb_updiffuse_all = bb_updiffuse_clr
        bb_dndiffuse_all = bb_dndiffuse_clr
        bb_dndirect_all = bb_dndirect_clr
        flx_updiffuse_all(:,:) = flx_updiffuse_clr(:,:)
        flx_dndiffuse_all(:,:) = flx_dndiffuse_clr(:,:)
        flx_dndirect_all(:,:) = flx_dndirect_clr(:,:)
    endif ! end of if c_flag


    end subroutine run_modtran_or_pcrtm
   
end module write_modtran_tp5
