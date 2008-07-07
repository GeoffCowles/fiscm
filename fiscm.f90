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
  character(len=fstr) :: forcing_file

  Namelist /NML_FISCM/ &
      beg_time,       & 
      end_time,       &
      deltaT,         &
      Ngroups,        &
      forcing_file
  
End Module fiscm_data

Program fiscm
  use gparms 
  use fiscm_data
  use mod_igroup
  use mod_AD
  use bio
  use output_routines
  use forcing
  implicit none
  integer :: i,n
  real(sp):: t = 0.0

  !---------------------------------------------------
  ! read primary simulation data and setup groups
  !---------------------------------------------------
  call setup

  !---------------------------------------------------
  ! initialize groups and data
  !---------------------------------------------------
  do n=1,ngroups
    if(igroups(n)%biology)call init_bio(igroups(n),igroups(n)%Tnind)
  end do

  !---------------------------------------------------
  ! define variables for output
  !---------------------------------------------------
  call cdf_out(ngroups,igroups,t,NCDO_ADD_STATES)

  !---------------------------------------------------
  ! setup and simulation summary
  !---------------------------------------------------
  do n=1,ngroups
    call print_group_summary(igroups(n))
  end do

  !---------------------------------------------------
  ! open input file
  ! ensure input file contains necessary data
  !---------------------------------------------------
  do n=1,ngroups
    call check_forcing(igroups(n),forcing_file,beg_time,end_time,sim_direction)
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

    call cdf_out(ngroups,igroups,t,NCDO_OUTPUT)

    if(.not. checkstatus(ngroups,igroups,t))then
      write(*,*)'no active particles left in the simulation'
      write(*,*)'shutting down prematurely'
      exit
    endif

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
  use output_routines
  implicit none
  logical :: fexist
  integer, parameter :: iunit = 33
  integer :: n,ios
  real(sp) :: uno = 1.0
  real(sp) :: t   = 0.0

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
    sim_direction = sign(uno, end_time-beg_time)
  else
    write(*,*)'fatal error: begin and end time are identical'
    stop
  endif

  !convert time step to days
  deltaT = deltaT*sec_2_day

  write(*,*)'begin time:  ',beg_time
  write(*,*)'end time  :  ',end_time
  if(sim_direction ==  1)  write(*,*)'direction:      forward'
  if(sim_direction == -1)  write(*,*)'direction:      backward'
  write(*,*)'time step(s) ',deltaT*day_2_sec
  write(*,*)'num groups:  ',ngroups

  !read and allocate individual groups
  allocate(igroups(ngroups))
  do n=1,ngroups
    igroups(n) = group_(iunit,n)
  end do

  !open output files and assign netcdf ids to each group
  call cdf_out(ngroups,igroups,t,NCDO_HEADER)
  
End Subroutine setup
