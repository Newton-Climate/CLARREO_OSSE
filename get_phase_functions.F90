#include <misc.h>
#include <params.h>

module get_phase_functions

  !Purpose:
  !   get aerosol phase functions from Modtran aerosol properties
  !   Auther: D. Feldman

  !  brdf properties from MODIS MCD43C1 product

  use filenames, only: phase_funcs
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  !quantities of interest
  real(r8), public :: phase_rhs(8) ! in % rh
  real(r8), public :: phase_wvls(15) ! in um
  real(r8), public :: phase_angles(50) ! in degrees
  real(r8), public :: phase_functions(50,15,8,4) !50 angles, 15 wavelength,8 relative humidities,4 species
  !quantities of interest

  save

contains

  subroutine get_pfs()

    use ioFileMod, only: getfil
    implicit none

    include 'netcdf.inc'

    !local variables for loops
    integer :: i,j,k,m,n
        
    character(len=256) :: locfn   !local file
    integer :: nc_id

    real(r8) :: dummy_array(4)

    !identifier from phase function file
    integer :: pf_id,rh_id,wvl_id,ang_id

    integer :: start(4)
    integer :: count(4)

    write (6, '(2x, a)') '_______________________________________________________'
    write (6, '(2x, a)') '_______ initializing phase function properties ________'
    write (6, '(2x, a)') '_______________________________________________________'

    call getfil(phase_funcs, locfn, 0)
    call wrap_open(locfn, 0, nc_id)
    
    call wrap_inq_dimid(nc_id, 'rh',rh_id)
    call wrap_inq_dimid(nc_id, 'wavelength',wvl_id)
    call wrap_inq_dimid(nc_id, 'angle',ang_id)
    call wrap_inq_varid(nc_id, 'phase_function', pf_id)

    call wrap_get_var_realx(nc_id,rh_id,phase_rhs)
    call wrap_get_var_realx(nc_id,wvl_id,phase_wvls)
    call wrap_get_var_realx(nc_id,ang_id,phase_angles)
    
    do i=1,50
       do j=1,15
          do k=1,8
             start = (/1,k,j,i/)
             count = (/4,1,1,1/)
             call wrap_get_vara_realx(nc_id,pf_id,start,count,dummy_array)
             do m=1,4
                phase_functions(i,j,k,m) = dummy_array(m)
             end do
          end do
       end do
    end do

    call wrap_close(nc_id)

  end subroutine get_pfs

end module get_phase_functions
