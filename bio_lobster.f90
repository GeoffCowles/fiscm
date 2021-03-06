!=======================================================================
! Fiscm Biology Example 
!
! Description
!  - i-state model for Homerus americanus (from L. Incze et al.) 
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
  integer , pointer :: pid(:)
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  real(sp), pointer :: s(:)
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

  

  !set the spawning time
  call get_state('tspawn',g,tspawn)
  call get_state('pid',g,pid)
  tspawn = 0.0
  pid = 0.0
  
  !set the spawning location
  !-----------------------------------
  ! 2D,3D -> initialize x,y 
  !-----------------------------------
  if(g%space_dim > 1)then
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
    
  endif
  !-----------------------------------
  ! 3D -> initialize s-coordinate
  !-----------------------------------
  if(g%space_dim > 2)then

    call get_state('s',g,s)
    call get_state('x',g,x)
    call get_state('y',g,y)
    allocate(zpini(g%Nind))
    s = 0;

   !initial condition file
    open(unit=33,file=trim(g%init_pos_file),form='formatted')
    read(33,*)ii,nhead
    write(*,*)nhead
    do i=1,nhead
      read(33,*)
    end do
    do i=1,g%Nind
      read(33,*)x(i),y(i),zpini(i),tspawn(i),pid(i)
    end do
    if(sz_cor == 0)then
      s = zpini
    elseif(sz_cor == 1)then
      write(*,*)'error initiating vertical positions, not setup for z-coord'
      stop
    endif 
    close(33)

  endif



  !add parameters to netcdf file header 
  call add_cdfglobs(g,"info","some kind of info")
  call add_cdfglobs(g,"nstages",4)

  
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
  real(sp),     pointer :: PASD(:)
  real(sp),     pointer :: T(:)
  real(sp),     pointer :: s(:)
  integer ,     pointer :: stage(:)
  integer ,     pointer :: istatus(:)
  integer               :: i,N
  real(sp)              :: deltaT,D

  !construct pointers to access and modify state variables for the group
  call get_state('PASD',g,PASD)
  call get_state('T',g,T)
  call get_state('status',g,istatus)
  call get_state('stage',g,stage)
  call get_state('s',g,s)
  
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
    !update PASD using stage-based Duration
    select case (stage(i))
    case(1)
      D = 851.*(T(i)-0.84)**(-1.91)
    case(2)
      D = 200*(T(i)-4.88)**(-1.47)
    case(3)
      D = 252*(T(i)-5.30)**(-1.45)
    case(4)
      D = 703.5*(T(i))**(-1.26)  !!!.358833*T(i)**2 - 14.316*T(i) + 156.895
    case(5)
      D = 0.0
    case default
      write(*,*)'stage: ',stage(i),' not a valid stage for a Homerus'
      stop
    end select

    PASD(i) = PASD(i) + deltaT/(D)

    !settle the post-larvae 
    if(PASD(i) > 1.0)then
      stage(i) = stage(i) + 1
      PASD(i)  = 0.0

    endif
  
    !settle at stage 5
    if(stage(i) == 5) istatus(i) = SETTLED
   
    !set layer 
    if(stage(i) < 4)then
      s(i) = -1. !bottom layer
    else
      s(i) = 0.  !neustonic
    endif
 
  end do !end main loop


end subroutine advance_bio 


End Module Bio 
