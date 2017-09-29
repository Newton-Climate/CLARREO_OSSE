c##############################################################################

      subroutine alloc_err( istat, routine, name, nelem )

c Issue error message after non-zero return from an allocate statement.

      implicit none

c Input arguments:
      integer :: istat   ! status from allocate statement
      character(len=*) ::
     $  routine          ! routine that called allocate
     $, name             ! name of array
      integer :: nelem   ! number of elements attempted to allocate
c-----------------------------------------------------------------------

      if ( istat .ne. 0 ) then
         write(*,*)'utils: ERROR trying to allocate memory in routine:'
     $          //trim(routine)
         write(*,*)'  Variable name: '//trim(name) 
         write(*,*)'  Number of elements: ',nelem
         call endrun
      end if
      return
      end

c##############################################################################

      subroutine endrun

c Abort the model.

      implicit none
c ----------------------------------------------------------------------

      call abort
      end

c##############################################################################

