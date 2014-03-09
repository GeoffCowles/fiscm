!=======================================================================
! Fiscm Vertical Diffusion Example 
!
! Description
!
! Comments:
!  - contains two routines:
!     init_bio: read optional namelist to set params
!               set initial positions (x,y,s)
!               add global information to NetCDF output
!  advance_bio: advance the biological state in time 
!
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!
!=======================================================================

Module Bio 

Use gparms
Use mod_igroup

Implicit none




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
    write(*,*)'here'
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
    read(33,*)ii
    do i=1,g%Nind
      read(33,*)x(i),y(i),s(i),ff
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

  call get_state('status',g,istatus)
  N = g%nind
  do i=1,N !main loop
    istatus(i) = ACTIVE
  end do



end subroutine advance_bio 


End Module Bio 

