!=======================================================================
! Fiscm Advection-Diffusion Routines
! Copyright:    2008(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Comments:     FISCM Global Type and Associated Functions
!=======================================================================

Module mod_AD
use gparms
use mod_igroup

implicit none
contains

!----------------------------------------------------
! advection-diffusion driver
!----------------------------------------------------
subroutine adv_diff(ng,g,time)
  integer     , intent(in   ) :: ng
  type(igroup), intent(inout), dimension(ng) :: g
  real(sp)    , intent(in   ) :: time
  integer :: n
  do n=1,ng
	
	select case(g(n)%problem_dimension)
	case(0)
	  cycle
	case(1)
	  write(*,*)'1D simulation not setup'
	  stop
	case(2)
	  call advect2D( g(n) )
	  if(g(n)%hdiff_type == HDIFF_CONSTANT) call rw_hdiff_const( g(n) )
	  if(g(n)%hdiff_type == HDIFF_VISSER)   call rw_hdiff_visser( g(n) )
	case(3)
	  call advect3D( g(n) )
	  if(g(n)%hdiff_type == HDIFF_CONSTANT) call rw_hdiff_const( g(n) )
	  if(g(n)%hdiff_type == HDIFF_VISSER)   call rw_hdiff_visser( g(n) )
!   	  if(g%vdiff_type == VDIFF_CONSTANT) call vdiff_const( g(n) )
!	  if(g%vdiff_type == VDIFF_VISSER)   call vdiff_visser( g(n) )
	case default
	  write(*,*)'problem_dimension must be [0,2,3]'
	  stop
	end select
	
  end do

end subroutine adv_diff
	
!----------------------------------------------------
! Random-Walk horizontal diffusion with constant
!    turbulent eddy diffusivity
!----------------------------------------------------
subroutine rw_hdiff_const(g)
  type(igroup), intent(inout) :: g

end subroutine rw_hdiff_const


!----------------------------------------------------
! Random-Walk horizontal diffusion with spatially
!   variable turbulent eddy diffusivity
!----------------------------------------------------
subroutine rw_hdiff_visser(g)
  type(igroup), intent(inout) :: g

end subroutine rw_hdiff_visser

!----------------------------------------------------
! 2-D Advection (use vertically-averaged velocity)
!----------------------------------------------------
subroutine advect2D(g)
  type(igroup), intent(inout) :: g

end subroutine advect2D

!----------------------------------------------------
! 3-D Advection 
!----------------------------------------------------
subroutine advect3D(g)
  type(igroup), intent(inout) :: g

end subroutine advect3D

End Module mod_AD
	
	
	
