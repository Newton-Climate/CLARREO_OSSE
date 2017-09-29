module error_messages

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! General purpose routines for issuing error messages.
   ! 
   ! Author: B. Eaton
   ! 
   !----------------------------------------------------------------------- 
   use abortutils, only: endrun

   implicit none
   save
   private
   public :: &
      alloc_err,      &! Issue error message after non-zero return from an allocate statement.
#ifdef OFFLINE
      dealloc_err,    &! Issue error message after non-zero return from a deallocate statement.
#endif
      handle_err,     &! Issue error message after non-zero return from anything
      handle_ncerr     ! Handle error returns from netCDF library procedures.

!##############################################################################
contains
!##############################################################################

   subroutine alloc_err( istat, routine, name, nelem )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Issue error message after non-zero return from an allocate statement.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      integer, intent(in) ::&
         istat           ! status from allocate statement
      character(len=*), intent(in) ::&
         routine,       &! routine that called allocate
         name            ! name of array
      integer, intent(in) ::&
         nelem           ! number of elements attempted to allocate
      !-----------------------------------------------------------------------

      if ( istat .ne. 0 ) then
         write(*,*)'ERROR trying to allocate memory in routine: ' &
                   //trim(routine)
         write(*,*)'  Variable name: '//trim(name)
         write(*,*)'  Number of elements: ',nelem
         call endrun ('ALLOC_ERR')
      end if

      return

   end subroutine alloc_err

#ifdef OFFLINE
!##############################################################################

   subroutine dealloc_err( istat, routine, name )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Issue error message after non-zero return from an deallocate statement.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      implicit none

      integer, intent(in) ::&
         istat           ! status from allocate statement
      character(len=*), intent(in) ::&
         routine,       &! routine that called allocate
         name            ! name of array
      !-----------------------------------------------------------------------

      if ( istat .ne. 0 ) then
         write(*,*)'ERROR trying to deallocate memory in routine: ' &
                   //trim(routine)
         write(*,*)'  Variable name: '//trim(name)
         call endrun
      end if

   end subroutine dealloc_err

!##############################################################################
#endif

!##############################################################################

   subroutine handle_err( istat, msg )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Issue error message after non-zero return from anything.
      !
      ! Author: T. Henderson
      !----------------------------------------------------------------------- 

      integer,          intent(in) :: istat  ! status, zero = "no error"
      character(len=*), intent(in) :: msg    ! error message to print
      !-----------------------------------------------------------------------

      if ( istat .ne. 0 ) then
         call endrun (trim(msg))
      end if

      return

   end subroutine handle_err

!##############################################################################

   subroutine handle_ncerr( ret, mes )
      
      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Check netCDF library function return code.  If error detected 
      ! issue error message then abort.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      integer, intent(in) ::&
         ret                 ! return code from netCDF library routine
      character(len=*), intent(in) ::&
         mes                 ! message to be printed if error detected
      !-----------------------------------------------------------------------

      if ( ret .ne. NF_NOERR ) then
         write(6,*) mes
         write(6,*) nf_strerror( ret )
         call endrun ('HANDLE_NCERR')
      endif

      return

   end subroutine handle_ncerr

!##############################################################################

end module error_messages
