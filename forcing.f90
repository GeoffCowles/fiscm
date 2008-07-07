!=======================================================================
! Fiscm NetCDF Forcing Routines
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
module forcing
use gparms
use mod_igroup
use netcdf
implicit none


integer, private :: ffile_id
character(len=fstr), private :: force_file


contains

!-------------------------------------------
! check forcing:
!   1.) see if forcing file is necessary
!   2.) if yes, make sure it exists
!   3.) if yes, make sure variables are there
!   4.) make sure time ranges make sense i.e.
!     for forwards sim:  
!       a.) first spawning after begin time
!       b.) last spawning before end time
!       c.) forcing file comprises begin/end time
!     for backwards sim:
!       a.) analogous   
!-------------------------------------------
subroutine check_forcing(g,ffile,tbeg,tend,simdir)
  type(igroup),intent(inout) :: g
  logical :: need_ext = .false.
  character(len=fstr) :: ffile
  real(sp) :: tbeg,tend
  integer  :: simdir
  real(sp) :: tspawn_min,tspawn_max
  real(sp), pointer :: tspawn(:)
  integer :: N

  !set dim
  N = g%Nind

  !check spawning time vs. specified model begin/end time
  call get_state('tspawn',g,tspawn)
  tspawn_min = minval(tspawn(1:N))
  tspawn_max = maxval(tspawn(1:N))
  if(simdir > 0)then
	if(tspawn_min < tbeg)then
	  write(*,*)'fatal error: first spawn of group: ',g%id,' occurs before simulation begins'
	  write(*,*)'first spawn occurs at: ',tspawn_min
	  write(*,*)'simulation begins at: ',tbeg
	  stop
	else if(tspawn_max > tend)then
      write(*,*)'fatal error: last spawn of group: ',g%id,' occurs after simulation ends'
	  write(*,*)'last spawn occurs at: ',tspawn_max
	  write(*,*)'simulation ends at: ',tend
	  stop
	endif
  else if(simdir < 0)then
	if(tspawn_min < tend)then
      write(*,*)'fatal error: last spawn of group: ',g%id,' occurs after simulation ends'
	  write(*,*)'last spawn occurs at: ',tspawn_max
	  write(*,*)'simulation ends at: ',tend
	  stop
	else if(tspawn_max > tbeg)then
	  write(*,*)'fatal error: last spawn of group: ',g%id,' occurs before simulation begins'
	  write(*,*)'first spawn occurs at: ',tspawn_min
      write(*,*)'simulation begins at: ',tbeg
      stop
	endif
  endif

  !check if any state state variables require external vars

  !if no ext vars required and sim is 0-Dimensional, don't need forcing file
  if(.not. need_ext .and. g%space_dim < 2)return

  !if need forcing file, check if it has been specified and if it exists
  if(ffile == "NONE")then
	write(*,*)'Fatal error: group ',g%id, 'requires forcing but NONE specified'
	write(*,*)'for forcing_file in fiscm.nml'
	stop
  endif

  force_file = ffile

  !get time range from forcing file

  !make sure time range is suitable given problem end and begin and spawn time

  !make sure necessary external vars are present in the file

end subroutine check_forcing
  

End Module forcing