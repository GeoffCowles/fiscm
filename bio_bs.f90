!=======================================================================
! Fiscm Biology Example 
!
! Description
!  - i-state model for a simple creature 
!  - uses temperature-dependent stage duration 
!
! Comments:
!  - contains two routines:
!     init_bio: read optional namelist to set params
!               set initial positions (x,y,s)
!               set spawning times
!               add global information to NetCDF output
!  advance_bio: advance the biological state in time 
!               e.g. grow, die, settle, reproduce
!
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!
!=======================================================================

Module Bio 

Use gparms
Use mod_igroup

Implicit none

!parameters
integer,parameter  :: nstages = 4
real(sp) :: days_at_liberty = 1.0
real(sp) :: some_other_var  = 1.0


!namelist can be read from input to overide parameters
  Namelist /NML_BS/ &
     & days_at_liberty,     &
     & some_other_var

contains 

!------------------------------------------------------------------
! initialize a group
!   - read the control namelist (optional)
!   - add state variables and initialize them
!   - modify state variables (position, spawn time, initial weight, etc.)
!------------------------------------------------------------------
subroutine init_bio(g,Nind_start)
  use output_routines
  use utilities
  implicit none
  type(igroup), intent(inout) :: g
  integer,      intent(in)    :: Nind_start
  integer :: iunit,ios
  logical :: fexist
  real(sp), pointer :: tspawn(:)
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  real(sp), pointer :: s(:)
  real(sp) :: ff
  integer :: i,ii,nhead

  !set current problem size
  if(Nind_start < 0)then
    write(*,*)'must specify more than 0 individuals in init_bio'
    stop
  else if(Nind_start > g%Tnind)then
    write(*,*)'number of individuals in init_bio (',Nind_start,') exceeds problem size',g%Tnind
    stop
  endif
  g%Nind = Nind_start

  
  !read the namelist (if any) 
  if(g%paramfile /= "NONE")then
  inquire(file=trim(g%paramfile),exist=fexist)
  if(.not.fexist)then
    write(*,*)'fatal error: namelist file: ',trim(g%paramfile),' does not exist, stopping...'
    stop 
  endif
  iunit = 33
  open(unit=iunit,file=trim(g%paramfile),form='formatted')
  read(unit=iunit,nml=nml_bs,iostat=ios)
  if(ios /= 0)then
    write(*,*)'fatal error: could not read bs namelist from: ',trim(g%paramfile)
    stop
  endif
  close(iunit)
  endif ! file /= NONE

  !set the spawning time
  call get_state('tspawn',g,tspawn)
  tspawn = 0.0
  
  !set the spawning location
  !-----------------------------------
  ! 2D,3D -> initialize x,y 
  !-----------------------------------
  if(g%space_dim == 2)then
    call get_state('x',g,x)
    call get_state('y',g,y)
    !gom
    !call random_square(g%Nind,8.751e5_sp,9.01e5_sp,-7.55e4_sp,-4.567e4_sp,x,y) 
    !do i=1,g%Nind
    !  x(i) = 270000 + float(i-1)*1000. 
    !  y(i) = 154000. 
    !end do
  
    !fake_forcing
    !do i=1,g%Nind
    !  x(i) = 0.0
    !  y(i) = 0.0 
    !end do

    !initial condition file
    write(*,*)'reading ic file 2d'
    open(unit=33,file=trim(g%init_pos_file),form='formatted')
    read(33,*)ii,nhead
    do i=1,nhead
      read(33,*)
    end do
    do i=1,g%Nind
      read(33,*)x(i),y(i),ff,tspawn(i)
    end do
    close(33)
    
    
  endif
  !-----------------------------------
  ! 3D -> initialize s-coordinate
  !-----------------------------------
  if(g%space_dim == 3)then
    call get_state('s',g,s)
    call get_state('x',g,x)
    call get_state('y',g,y)
    s = 0;
    !do i=1,g%Nind
    !  s(i) = -float(i-1)/float(g%Nind-1)
    !end do
   !initial condition file
    write(*,*)'reading ic file'
    open(unit=33,file=trim(g%init_pos_file),form='formatted')
    read(33,*)ii,nhead
    write(*,*)nhead
    do i=1,nhead
      read(33,*)
    end do
    do i=1,g%Nind
      read(33,*)x(i),y(i),s(i),tspawn(i)
    end do
    close(33)

  endif
  write(*,*)'done with ic'

!  skagit
!  do i=1,g%Nind
!    x(i) = 540006 + float(i)*100.
!    y(i) = 5.3526256e6 + float(i)*100.
!  end do
!  test data
!  do i=1,g%Nind
!    x(i) = float(i-1)*100. 
!    y(i) = 0.0 
!  end do

  !add parameters to netcdf file header 
  call add_cdfglobs(g,"info","some kind of info")
  call add_cdfglobs(g,"nstages",4)
  call add_cdfglobs(g,"days_at_liberty",days_at_liberty)
  call add_cdfglobs(g,"some_other_var",some_other_var)

  
end subroutine init_bio 

!------------------------------------------------------------------
! advance the biology (this routine is called from the main loop at
!                      a time interval of DT_bio )
! use temperature-stage dependence from 
!
!     Incze, L.S., Naimie, C., 2000. Modeling the transport of lobster
!      (homarus americanus) larvae and postlarvae in the Gulf of
!      Maine. Fish. Oceanogr. 9, 99-113. 
!------------------------------------------------------------------
subroutine advance_bio(g,mtime)  
  type(igroup), intent(inout) :: g
  real(sp),     intent(in   ) :: mtime
  real(sp),     pointer :: tspawn(:)
  real(sp),     pointer :: T(:)
  integer ,     pointer :: stage(:)
  integer ,     pointer :: istatus(:)
  integer               :: i,N
  real(sp)              :: deltaT,D

  !construct pointers to access and modify state variables for the group
  call get_state('status',g,istatus)
  call get_state('stage',g,stage)
  call get_state('tspawn',g,tspawn)
  
  !set problem size
  N = g%nind
  deltaT = g%DT_bio*sec_2_day !time step in days 

!  following interface is cleaner but compiler cannot construct
!  generic interface if arguments are the same even if the
!  return value is of different type
!  PASD   => get_state('PASD',g)
!  istatus => get_state('status',g)
!  stage  => get_state('stage',g)

  do i=1,N !main loop

    if(istatus(i) /= ACTIVE)cycle
    
       if( (mtime - tspawn(i)) > days_at_liberty)then
         istatus(i) = SETTLED 
       endif;
 
  end do !end main loop


end subroutine advance_bio 


End Module Bio 

