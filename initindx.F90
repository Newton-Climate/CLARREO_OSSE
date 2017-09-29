subroutine initindx
!----------------------------------------------------------------------- 
! 
! Purpose: Register constituents and physics buffer fields.
! 
! Author:    CSM Contact: M. Vertenstein, Aug. 1997
!            B.A. Boville, Oct 2001
! 
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst, ppcnst, cnst_add, advected, cnst_chk_dim, cnst_name
#ifndef OFFLINE
  use phys_buffer,  only: pbuf_init
#endif
  use chemistry,    only: trace_gas, chem_register
#ifndef OFFLINE
  use cldcond,      only: cldcond_register
#endif
  use physconst,    only: mwdry, cpair, mwh2o, cph2o
#ifndef OFFLINE
  use tracers, only: tracers_register
#endif
!  use constituents, only: dcconnam, sflxnam, hadvnam, vadvnam, fixcnam, tendnam, tottnam
  use constituents, only: dcconnam, sflxnam
#ifndef OFFLINE
  use check_energy, only: check_energy_register
  use aerosol_intr, only: aerosol_register_cnst
#endif

  implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!---------------------------Local variables-----------------------------
!
  integer m            ! loop index
  integer mm           ! constituent index 
!-----------------------------------------------------------------------

! Initialize physics buffer
#ifndef OFFLINE
  call pbuf_init()
#endif

! Register water vapor.
! ***** N.B. ***** This must be the first call to cnst_add so that
!                  water vapor is constituent 1.
  call cnst_add('Q', advected, mwh2o, cph2o, 1.E-12_r8, mm, &
                longname='Specific humidity', readiv=.true.)
!
! Register cloud water
#ifndef OFFLINE
  call cldcond_register()
#endif

!
! Register chemical constituents
  if (trace_gas) then
     call chem_register()
  endif
!
! register aerosols
#ifndef OFFLINE
  call aerosol_register_cnst()
#endif

! Register advected test tracers and determine starting index
#ifndef OFFLINE
  call tracers_register()
#endif

!
! All tracers registered, check that the dimensions are correct
  call cnst_chk_dim()
!
! Set default names for non-water advected and non-advected tracers
! Set names of advected and non-advected tracer diagnostics
!
  do m=1,ppcnst
     dcconnam(m) = 'DC'//cnst_name(m)
     sflxnam(m)  = 'SF'//cnst_name(m)
  end do
!  do m=1,pcnst
!     hadvnam(m)  = 'HA'//cnst_name(m)
!     vadvnam(m)  = 'VA'//cnst_name(m)
!     fixcnam(m)  = 'DF'//cnst_name(m)
!     tendnam(m)  = 'TE'//cnst_name(m)
!     tottnam(m)  = 'TA'//cnst_name(m)
!  end do

#ifndef OFFLINE
  call check_energy_register()
#endif

end subroutine initindx
