!=======================================================================
! Fiscm NetCDF Output Routines
! Copyright:    2008(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Comments:  Requires fortran90 NetCDF 3x libraries
!=======================================================================
module output_routines
use gparms
use mod_igroup
implicit none

contains
!-------------------------------------------
! output driver
!-------------------------------------------
subroutine output(ng,g,time)
  integer, intent(in) :: ng
  type(igroup),intent(inout), dimension(ng) :: g
  real(sp), intent(in) :: time
  integer :: n

   do n=1,ng
     if( (time-g(n)%tlast_out) < g(n)%intvl_out .or. time < g(n)%start_out) cycle
     if(g(n)%frame_out == 0) call write_header(g(n),time)
     !call output_group(g(n),time)
   end do
end subroutine output

!-------------------------------------------
! write netcdf header
!-------------------------------------------
subroutine write_header(g,time)
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: time
  write(*,*)'here',time

  g%frame_out = 1
  

end subroutine write_header
  

end module output_routines