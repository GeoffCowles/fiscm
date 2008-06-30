!=======================================================================
! Fiscm Biology Example
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

Module Bio

Use gparms, only : sp,fstr
Use mod_igroup

Implicit none
integer,  parameter :: nstages = 4 
real(sp), parameter :: fcomp_settle = 0.5

contains

subroutine init_bio(g,Nind)
  type(igroup), intent(inout) :: g
  integer,      intent(in)    :: Nind
  real(sp),     pointer :: x(:)
  real(sp),     pointer :: y(:)
  
  x => get_from_state('x',g)
  write(*,*) x(1:5)
  y => get_from_state('y',g)
  write(*,*) y(1:5)
  x = x*x

end subroutine init_bio

subroutine advance_bio(g)
  type(igroup), intent(inout) :: g

real(sp),     pointer :: x(:)

x => get_from_state('x',g)
write(*,*) x(1:5)

end subroutine advance_bio


End Module Bio