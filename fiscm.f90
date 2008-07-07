!=======================================================================
! Fiscm Main 
! Copyright:    2008(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Major Todos:
!    1.) migrate to D.S. netcdf libraries
!    2.) migrate to D.S. time type
!    3.) openmp doloops to multithread heavier sections 
!        parallel random number generator?
!=======================================================================
Module fiscm_data
  use gparms
  use mod_igroup
  implicit none
  type(igroup), allocatable :: igroups(:)
  real(sp) :: beg_time,end_time
  real(sp) :: deltaT
  integer  :: sim_direction
  integer  :: ngroups

  Namelist /NML_FISCM/ &
     & beg_time,       & 
     & end_time,       &
     & deltaT,  &
     & Ngroups 
  
End Module fiscm_data

Program fiscm
  use gparms 
  use fiscm_data
  use mod_igroup
  use mod_AD
  use bio
  use output_routines
  integer :: i,n

  !---------------------------------------------------
  ! read primary simulation data and setup groups
  !---------------------------------------------------
  call setup

  !---------------------------------------------------
  ! initialize groups and data
  !---------------------------------------------------
  call init_bio(igroups(1),igroups(1)%Tnind)

  !---------------------------------------------------
  ! setup and simulation summary
  !---------------------------------------------------
  do n=1,ngroups
    call print_group_summary(igroups(n))
  end do

  !---------------------------------------------------
  ! main loop over time 
  !---------------------------------------------------
  t = beg_time
  do while (t <= end_time)
	
	call adv_diff(ngroups,igroups,t)
	
	do n=1,ngroups
      if(igroups(n)%biology)call advance_bio(igroups(n),t)
    end do

    call output(ngroups,igroups,t)

    t = t + deltaT
  end do

  !---------------------------------------------------
  ! cleanup data
  !---------------------------------------------------
  deallocate(igroups)


End Program fiscm

Subroutine setup
  use gparms
  use fiscm_data
  use mod_igroup
  implicit none
  logical :: fexist
  integer, parameter :: iunit = 33
  integer :: n,ios

  !--------------------------------------------------------
  ! check for existence of primary control file
  !--------------------------------------------------------
  inquire(file='fiscm.nml',exist=fexist)
  if(.not.fexist)then
    write(*,*)'fatal error: namelist file: fiscm.nml does not exist, stopping...'
    stop 
  endif
  open(unit=iunit,file='fiscm.nml',form='formatted')
  read(unit=iunit,nml=nml_fiscm,iostat=ios)
  if(ios /= 0)then
    write(*,*)'fatal error: could not read fiscm namelist from fiscm.nml'
    stop
  endif

  !set values, check validity and report
  if(ngroups < 1)then
    write(*,*)'Fatal error: number of groups < 1' ; stop 
  endif

  if(end_time /= beg_time)then
    sim_direction = sign(1., end_time-beg_time)
  else
    write(*,*)'fatal error: begin and end time are identical'
    stop
  endif

  write(*,*)'begin time:  ',beg_time
  write(*,*)'end time  :  ',end_time
  if(sim_direction ==  1)  write(*,*)'direction:      forward'
  if(sim_direction == -1)  write(*,*)'direction:      backward'
  write(*,*)'time step:   ',deltaT
  write(*,*)'num groups:  ',ngroups

  !read and allocate individual groups
  allocate(igroups(ngroups))
  do n=1,ngroups
    igroups(n) = group_(iunit,n)
  end do

  !open output files and assign netcdf ids to each group
  call cdf_out(ngroups,igroups,0.0,NCDO_HEADER)
  

End Subroutine setup
