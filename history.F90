module history
  use shr_kind_mod, only: r8 => shr_kind_r8                                   
contains
  subroutine outfld (fname, field, idim, c)
!----------------------------------------------------------------------- 
! 
! Purpose: Accumulate (or take min, max, etc. as appropriate) input field
!          into its history buffer for appropriate tapes
! 
! Method: Search for fname among fields on history tapes.  If found, do the
!         accumulation.  If not found, return silently.
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Arguments
!
    character(len=*), intent(in) :: fname ! Field name--should be 8 chars long

    integer, intent(in) :: idim           ! Longitude dimension of field array
    integer, intent(in) :: c              ! chunk (physics) or latitude (dynamics) index

    real(r8), intent(in) :: field(idim,*) ! Array containing field values
    return
  end subroutine outfld

end module history
