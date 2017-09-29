#include <misc.h>

!
! Local dumbed down version by W. Collins
!

module time_manager

   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils,   only: endrun
   use infnan,       only: bigint

   implicit none
   save


   integer ::  mdbase                    ! base day of run (e.g., 0)
   integer ::  msbase                    ! base seconds of base day (e.g., 0)
   integer ::  mbdate                    ! base date (yyyymmdd format) of run
   integer ::  mbsec                     ! base seconds of base date (e.g., 0)
   integer ::  mdcur                     ! current day (0, 1, ...)
   integer ::  mscur                     ! current seconds of current day (0, ..., 86400)
   integer ::  mcdate                    ! current date (yyyymmdd format) (e.g., 021105)
   integer ::  mcsec                     ! current seconds of current date (0, ..., 86400)
   real(r8) :: calday                    ! current calendar day

   
! Public methods

   public ::&
      get_curr_date,            &! return date components at end of current timestep
      get_perp_date,            &! return components of the perpetual date, and current time of day
      get_curr_calday,          &! return calendar day at end of current timestep
      is_perpetual,             &! return true if perpetual calendar is in use
      timemgr_datediff,         &! computes differences between dates
      is_end_curr_day            ! flag for end of current day

 
!=========================================================================================
contains
!=========================================================================================

subroutine get_curr_date(yr, mon, day, tod, offset)

! Return date components valid at end of current timestep with an optional
! offset (positive or negative) in seconds.

   implicit none
   
! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.
   
   if (present(offset)) then
      if (offset .ne. 0) then 
         write (*,*) "get_curr_date: offset != 0"
         stop
      endif
   endif

   yr  = mcdate / 10000
   mon = mod(mcdate / 100, 100)
   day = mod(mcdate, 100)
   tod = mcsec

end subroutine get_curr_date
!=========================================================================================

subroutine get_perp_date(yr, mon, day, tod, offset)

! Return time of day valid at end of current timestep and the components
! of the perpetual date (with an optional offset (positive or negative) in seconds.

   implicit none
   
! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.

   write(*,*) "get_perp_date called -- should never occur"
   stop

end subroutine get_perp_date

!=========================================================================================

function get_curr_calday(offset)

! Return calendar day at end of current timestep with optional offset.
! Calendar day 1.0 = 0Z on Jan 1.

   implicit none
   
! Arguments
   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.
! Return value
   real(r8) :: get_curr_calday

   if (present(offset)) then
      if (offset .ne. 0) then 
         write (*,*) "get_curr_date: offset != 0"
         stop
      endif
   endif

   get_curr_calday = calday
   return

end function get_curr_calday
!=========================================================================================

function is_perpetual()

! Return true on last timestep.

   implicit none
   
! Return value
   logical :: is_perpetual
!-----------------------------------------------------------------------------------------

   is_perpetual = .FALSE.
   return

end function is_perpetual

!=========================================================================================

subroutine timemgr_datediff(ymd1, tod1, ymd2, tod2, days)

! Calculate the difference (ymd2,tod2) - (ymd1,tod1) and return the result in days.

   implicit none

! Arguments
   integer, intent(in) ::&
      ymd1,    &! date1 in yyyymmdd format
      tod1,    &! time of day relative to date1 (seconds past 0Z)
      ymd2,    &! date2 in yyyymmdd format
      tod2      ! time of day relative to date2 (seconds past 0Z)

   real(r8) :: days ! (ymd2,tod2)-(ymd1,tod1) in days

!-----------------------------------------------------------------------------------------

end subroutine timemgr_datediff

!=========================================================================================

function is_end_curr_day()

! Return true if current timestep is last timestep in current day.

   implicit none
   
! Return value
   logical :: is_end_curr_day

   is_end_curr_day = .false.

end function is_end_curr_day

end module time_manager
