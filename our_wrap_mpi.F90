#include <misc.h>
#include <params.h>

module our_wrap_mpi


#if ( defined MODTRAN_SPMD )
   use mpi
#endif


   implicit none


!   private

! PUBLIC: Public interfaces

!   public :: shr_mpi_chkerr
!   public :: shr_mpi_send
!   public :: shr_mpi_recv
!   public :: shr_mpi_bcast
!   public :: shr_mpi_sum
!   public :: shr_mpi_min
!   public :: shr_mpi_max
!   public :: shr_mpi_commsize
!   public :: shr_mpi_commrank
!   public :: shr_mpi_initialized
!   public :: shr_mpi_abort
!   public :: shr_mpi_barrier
!   public :: shr_mpi_init
!   public :: shr_mpi_finalize


    integer :: wvlng_lres_block, wvlng_hres_block
        














!---------------------------------------------------------------------------
!
! Purpose:
!
! 	Wrapper routines for the MPI (Message Passing) library for the
!	distributed memory (SPMD) version of the code. Also data with
!	"shorthand" names for the MPI data types.
!
! Author: Many
!
!---------------------------------------------------------------------------
!
! Compile these routines only when SPMD is defined
!


contains


#if (defined MODTRAN_SPMD)


!****************************************************************
 
   subroutine our_mpiinit
!
! End of all MPI communication
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils, only: endrun
   implicit none
 
   integer :: ierr   !MP error code
 
   call mpi_init( ierr )
   if (ierr/=mpi_success) then
      write(6,*)'mpi_init failed ierr=',ierr
      call endrun
   end if
 

   return
   end subroutine our_mpiinit
 
!****************************************************************
 
   subroutine our_mpitypevector()
!
! End of all MPI communication
!
   use abortutils, only: endrun
   implicit none
 
   integer :: ierr   !MP error code
 
   !call mpi_init( ierr )
   !if (ierr/=mpi_success) then
   !   write(6,*)'mpi_init failed ierr=',ierr
   !   call endrun
   !end if
 

   call MPI_TYPE_vector(WVLNGTH_LRES,1,PLON,MPI_REAL, wvlng_lres_block, ierr)
   call MPI_TYPE_COMMIT(wvlng_lres_block, ierr)

   call MPI_TYPE_vector(WVLNGTH_HRES,1,PLON,MPI_REAL, wvlng_lres_block, ierr)
   call MPI_TYPE_COMMIT(wvlng_lres_block, ierr)

   return
   end subroutine our_mpitypevector
 
!****************************************************************



























 
   subroutine our_mpifinalize
!
! End of all MPI communication
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils, only: endrun
   implicit none
 
   integer :: ierr   !MP error code
 
   call mpi_finalize (ierr)
   if (ierr/=mpi_success) then
      write(6,*)'mpi_finalize failed ierr=',ierr
      call endrun
   end if
 
   return
   end subroutine our_mpifinalize
 
!****************************************************************
 
!===============================================================================
!===============================================================================

SUBROUTINE our_mpicommsize(size)
   

   IMPLICIT none

   !----- arguments ---
   integer,intent(out)      :: size
   
   !----- local ---
   integer                  :: ierr  !MPI error code

!-------------------------------------------------------------------------------
! PURPOSE: MPI commsize
!-------------------------------------------------------------------------------

   call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
   if (ierr/=mpi_success) then
      write(6,*)'mpi_comm_size failed ierr=',ierr
      call endrun
   end if
  

END SUBROUTINE our_mpicommsize

!*******************************************************************************

SUBROUTINE our_mpicommrank(rank)


   IMPLICIT none

   !----- arguments ---
   integer,intent(out)                :: rank

   !----- local ---
   integer :: ierr   !MPI error code
   

!-------------------------------------------------------------------------------
! PURPOSE: MPI commrank
!-------------------------------------------------------------------------------

   call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
   if (ierr/=mpi_success) then
      write(6,*)'mpi_comm_rank failed ierr=',ierr
      call endrun
   end if

END SUBROUTINE our_mpicommrank

!*******************************************************************************


!****************************************************************
   ! 1 real
   subroutine mpibcastreal (buffer,count,root)
!
! Broadcasts a message from masterproc to all threads
!
   use shr_kind_mod, only: r4 => shr_kind_r4
   
   use abortutils, only: endrun
   implicit none
 
   real (r4), intent(inout):: buffer
   integer, intent(in):: root
   integer, intent(in):: count
   
 
   integer ierr   !MP error code
 
   call t_startf ('mpi_bcast')
   call mpi_bcast (buffer,count,MPI_REAL,root,MPI_COMM_WORLD, ierr)
   if (ierr/=mpi_success) then
      write(6,*)'mpi_bcast failed ier=',ierr
      call endrun
   end if
   call t_stopf ('mpi_bcast')
 
   return
   end subroutine mpibcastreal
!****************************************************************
 
   subroutine mpibcasthresreal (buffer,root)
!
! Broadcasts a message from masterproc to all threads
!
   use shr_kind_mod, only: r4 => shr_kind_r4
   
   use abortutils, only: endrun
   implicit none
 
   real (r4), intent(inout):: buffer(*)
   integer, intent(in):: root
  
   
 
   integer ierr   !MP error code
 
   call t_startf ('mpi_bcast')
   call mpi_bcast (buffer,wvlng_hres_block,MPI_REAL,root,MPI_COMM_WORLD, ierr)
   if (ierr/=mpi_success) then
      write(6,*)'mpi_bcast failed ier=',ierr
      call endrun
   end if
   call t_stopf ('mpi_bcast')
 
   return
   end subroutine mpibcasthresreal
!****************************************************************

!****************************************************************

   subroutine our_mpibarrier ()

   use abortutils, only: endrun

   implicit none
!
! MPI barrier, have threads wait until all threads have reached this point
!
    
   integer ierr   !MP error code
 
   call mpi_barrier (MPI_COMM_WORLD, ierr)
   if (ierr.ne.mpi_success) then
      write(6,*)'mpi_barrier failed ierr=',ierr
      call endrun
   end if
 
   return
   end subroutine our_mpibarrier
 
!****************************************************************
   !1 int
   subroutine our_mpisendint (buf, count, dest, tag)
!
! Does a blocking send
!
   use abortutils, only: endrun
   implicit none
 
   integer, intent(in):: buf
   integer, intent(in):: count
   integer, intent(in):: dest
   integer, intent(in):: tag
 
   integer ierr   !MP error code
 
!   call t_startf ('mpi_send')
   call mpi_send (buf, count, MPI_INTEGER, dest, tag, MPI_COMM_WORLD, ierr)
   if (ierr/=mpi_success) then
      write(6,*)'mpi_send failed ierr=',ierr
      call endrun
   end if
!   nsend = nsend + 1
!   nwsend = nwsend + count
!   call t_stopf ('mpi_send')
 
   return
   end subroutine our_mpisendint
 
!****************************************************************
   !1 int
   subroutine our_mpirecvint (buf, count, source, tag) 
!
! Does a blocking receive
!
   use abortutils, only: endrun
   implicit none
 
   integer, intent(out):: buf
   integer, intent(in):: count
   integer, intent(in):: source
   integer, intent(in):: tag
 
   integer status (MPI_STATUS_SIZE) ! Status of message
   integer ierr   !MP error code
 
!   call t_startf ('mpi_recv')
   call mpi_recv (buf, count, MPI_INTEGER, source, tag, MPI_COMM_WORLD, status, ierr)
   if (ierr/=mpi_success) then
      write(6,*)'mpi_recv failed ierr=',ierr
      call endrun
   end if
!   nrecv = nrecv + 1
!   nwrecv = nwrecv + count
!   call t_stopf ('mpi_recv')
 
   return
   end subroutine our_mpirecvint

!****************************************************************
 
   subroutine our_mpisendreal (buf, count, dest, tag)
!
! Does a blocking send
!
   use shr_kind_mod, only: r4 => shr_kind_r4
   use abortutils, only: endrun
   implicit none
 
   real (r4), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: dest
   integer, intent(in):: tag
 
   integer ierr   !MP error code
 
!   call t_startf ('mpi_send')
   call mpi_send (buf, count, MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr)
   if (ierr/=mpi_success) then
      write(6,*)'mpi_send failed ierr=',ierr
      call endrun
   end if
!   nsend = nsend + 1
!   nwsend = nwsend + count
!   call t_stopf ('mpi_send')
 
   return
   end subroutine our_mpisendreal
 
!****************************************************************
   subroutine our_mpirecvreal (buf, count, source, tag)
!
! Does a blocking receive
!
   use shr_kind_mod, only: r4 => shr_kind_r4
   use abortutils, only: endrun
   implicit none
 
   real (r4), intent(out):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: source
   integer, intent(in):: tag
 
   integer status (MPI_STATUS_SIZE) ! Status of message
   integer ierr   !MP error code
 
!   call t_startf ('mpi_recv')
   call mpi_recv (buf, count, MPI_REAL, source, tag, MPI_COMM_WORLD, status, ierr)
   if (ierr/=mpi_success) then
      write(6,*)'mpi_recv failed ierr=',ierr
      call endrun
   end if
!   nrecv = nrecv + 1
!   nwrecv = nwrecv + count
!   call t_stopf ('mpi_recv')
 
   return
   end subroutine our_mpirecvreal
 
!****************************************************************
!****************************************************************

   subroutine our_mpisenddouble (buf, count, dest, tag)
!
! Does a blocking send
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils, only: endrun
   implicit none
 
   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: dest
   integer, intent(in):: tag
 
   integer ierr   !MP error code
 
!   call t_startf ('mpi_send')
   call mpi_send (buf, count, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
   if (ierr/=mpi_success) then
      write(6,*)'mpi_send failed ierr=',ierr
      call endrun
   end if
!   nsend = nsend + 1
!   nwsend = nwsend + count
!   call t_stopf ('mpi_send')
 
   return
   end subroutine our_mpisenddouble
 
!****************************************************************
   subroutine our_mpirecvdouble (buf, count, source, tag)
!
! Does a blocking receive
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils, only: endrun
   implicit none
 
   real (r8), intent(out):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: source
   integer, intent(in):: tag
 
   integer status (MPI_STATUS_SIZE) ! Status of message
   integer ierr   !MP error code
 
!   call t_startf ('mpi_recv')
   call mpi_recv (buf, count, MPI_DOUBLE_PRECISION, source, tag, MPI_COMM_WORLD, status, ierr)
   if (ierr/=mpi_success) then
      write(6,*)'mpi_recv failed ierr=',ierr
      call endrun
   end if
!   nrecv = nrecv + 1
!   nwrecv = nwrecv + count
!   call t_stopf ('mpi_recv')
 
   return
   end subroutine our_mpirecvdouble


!****************************************************************
!****************************************************************
 





   double precision function our_MPI_Wtime()
 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils, only: endrun

   implicit none
!
! MPI barrier, have threads wait until all threads have reached this point
!
    
 
   
 
   real(r8) the_time 

   the_time = MPI_Wtime()
   
   our_MPI_Wtime = the_time

   end function our_MPI_Wtime
 
!****************************************************************































!
! If SPMD is not turned on
!
#else
   subroutine our_wrap_mpi
   use abortutils, only: endrun
   implicit none
!
! A unused stub routine to make the compiler happy when SPMD is
! turned off (which means you don't need anything in this file).
!
   call endrun ('(OUR_WRAP_MPI): This should not be called at all')
   end subroutine our_wrap_mpi

#endif


end module our_wrap_mpi
