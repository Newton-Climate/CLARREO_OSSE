#include <misc.h>
module phys_grid
!----------------------------------------------------------------------- 
! 
! Purpose: Definition of physics computational horizontal grid.
!
! Method: Variables are private; interface routines used to extract
!         information for use in user code.
! 
! Entry points:
!      get_ncols_p         get number of columns for a given chunk
!      scatter_field_to_chunk >DUMMY<
!                          distribute longitude/latitude field
!                          to decomposed chunk data structure
!
! Author: Patrick Worley and John Drake
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
   use ppgrid, only: pcols, pver, begchunk, endchunk
   use pmgrid, only: plon, plat, beglat, endlat
   use abortutils, only: endrun
#if ( defined SPMD )
   use spmd_dyn, only: proc, npes, nsmps, proc_smp_map
   use mpishorthand
#endif

   implicit none

   save

contains
!
!========================================================================
!
   integer function get_ncols_p(lchunkid)
!----------------------------------------------------------------------- 
! 
! Purpose: Return number of columns in chunk given the local chunk id.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id

   get_ncols_p = pcols

   return
   end function get_ncols_p

   subroutine scatter_field_to_chunk(fdim,mdim,ldim, &
                                     nlond,globalfield,localchunks)
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute longitude/latitude field
!          to decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam, masterproc
!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: nlond     ! declared number of longitudes
   real(r8), intent(in) :: globalfield(fdim,nlond,mdim,plat,ldim) 
                                    ! global field

   real(r8), intent(out):: localchunks(fdim,pcols,mdim, &
                                       begchunk:endchunk,ldim) 
                                    ! local chunks


   localchunks(:,1:min(pcols,nlond),:,:,:) = &
        globalfield(:,1:min(pcols,nlond),:,:,:)

   return
   end subroutine scatter_field_to_chunk

end module phys_grid
