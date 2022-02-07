!=======================================================================
! Fiscm Constant Bearing Example
!
! Description
!
! Comments:
!  - contains three routines:
!     init_bio: read optional namelist to set params
!               set initial positions (x,y,s)
!               add global information to NetCDF output
!  advance_bio: advance the biological state in time
!  advance_bio_movement: called in advect, implements horizontal
!                        movement due to behaviour
!
! !REVISION HISTORY:
!  Original author(s): G. Cowles
!  Modified by J. Ounsley
!
!=======================================================================

Module Bio

Use gparms
Use mod_igroup
Use ocean_model

Implicit none

! Parameters, to be configured
! The swimming behaviour type
! 0 = passive - no swimming
! 1 = constant bearing - 1 body length per second
! 2 = negative rheotaxis - swim with current
! 3 = positive rheotaxis - swim against current
! 4 = constant bearing range - range of bearings around const bearing
! 5 = temperature and salinity follwoing - follow currrents but stick to preferred
!                                          temp and salinity values Mork et al. 2012
! 6 = thermotaxis, stick to temperature range 4 - 8 C, with random movement in this
!                                           range, 1.6 bl per second or 0.2 m per second
!                                           Booker et al. 2008
! 7 = salinity preference only following mork
! 8 = salinity preference with constant bearing
! 9 = salinity preference with constant bearing - expontential weighting between the two
! 10 = Follow depth gradient
! 11 = Follow depth gradient when less than 100m weighted against constant bearing
! 12 = Booker et al. 2008 using temperature preference range of 8 - 11
integer  :: behaviour_type = 0
! The biology type
! 0 = no biology routine
! 1 = Mork et al. 2012, exponential growth
integer  :: bio_type = 0
real(sp) :: swim_speed = 1.0_sp ! Body lengths per second
real(sp) :: const_bearing = 0.0_sp
real(sp) :: bearing_range = 0.0_sp
integer  :: collision_avoidance = 0 ! Do we want to follow currents on collision
real(sp) :: transition_time = 1.0_sp ! Transition to second behaviour in vattenfall


! Namelist to configure parameters
namelist /NML_SALMON/     &
     & behaviour_type,    &
     & bio_type,          &
     & swim_speed,        &
     & const_bearing,     &
     & bearing_range,	  &
     & collision_avoidance, &
     & transition_time


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
  real(sp), pointer :: bearing(:)
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
     read(unit=iunit,nml=nml_salmon,iostat=ios)
     if(ios /= 0)then
        write(*,*)'fatal error: could not read salmon namelist from: ',trim(g%paramfile)
        stop
     endif
     close(iunit)
  endif ! file /= NONE
  write(*,*)'behaviour_type: ',behaviour_type
  write(*,*)'bio_type: ',bio_type
  write(*,*)'swim_speed: ',swim_speed
  write(*,*)'const_bearing: ',const_bearing
  write(*,*)'collision_avoidance: ',collision_avoidance
  if (behaviour_type == 4)  write(*,*)'bearing_range: ',bearing_range


  !-----------------------------------
  ! 3D -> initialize s-coordinate
  !-----------------------------------
  if(g%space_dim == 3)then
    ! JO - TODO Initialise particle size attribute here
    !      Could generate these in the ic file or create them
    !      here?


    ! JO - TODO This should probabily happen somewhere else
    allocate(zpini(g%Nind))

    call get_state('s',g,s)
    call get_state('x',g,x)
    call get_state('y',g,y)
    call get_state('tspawn',g,tspawn)
    s = 0;
    write(*,*)'reading ic file'
    open(unit=33,file=trim(g%init_pos_file),form='formatted')
    read(33,*)ii
    do i=1,g%Nind
      read(33,*)x(i),y(i),s(i),tspawn(i)
    end do
    close(33)

    ! Set initial values
    zpini = s

    ! Randomly generate bearing given
    ! the bearing range
    if (behaviour_type == 4 .or. behaviour_type == 9 .or.  &
		behaviour_type == 8 .or. behaviour_type == 11 .or. &
		behaviour_type == 13)then
       call get_state('bearing',g,bearing)
       do i=1,g%Nind
          bearing(i) = const_bearing + ran_from_range(-bearing_range/2,bearing_range/2)
       end do

    end if

  endif
  write(*,*)'done with ic'


  !add parameters to netcdf file header
  call add_cdfglobs(g,"behaviour_type",behaviour_type)
  call add_cdfglobs(g,"bio_type",bio_type)
  call add_cdfglobs(g,"swim_speed",swim_speed)
  call add_cdfglobs(g,"const_bearing",const_bearing)
  call add_cdfglobs(g,"bearing_range",bearing_range)
  call add_cdfglobs(g,"collision_avoidance",collision_avoidance)
  call add_cdfglobs(g,"transition_time",transition_time)


end subroutine init_bio

!------------------------------------------------------------------
! advance the biology (this routine is called from the main loop at
!  a time interval of DT_bio )
!
!------------------------------------------------------------------
subroutine advance_bio(g,mtime)
  type(igroup), intent(inout) :: g
  real(sp),     intent(in   ) :: mtime
  real(sp)              :: deltaT,D

  ! Determine the biological routine
  if(bio_type == 0) then
     ! Do nothing
  else if (bio_type == 1) then
     call bio_mork(g)
  else
     write (*,*) 'WARNING::::::: Unknown biology type ', bio_type, '. Biology ignored.'
  endif

end subroutine advance_bio

!----------------------------------------------------------------------
! Added by J. Ounsley
! Routine for biology replicating that of
! Mork et al. 2012
! Exponential growth given by rule
! L = 7.34e^0.0059*doy,
! or increase of 1.006 per day
!----------------------------------------------------------------------
subroutine bio_mork(g)
  type(igroup), intent(inout) :: g
  real(sp),     pointer       :: length(:)
  integer,   pointer       :: istatus(:)
  real(sp)                    :: deltaT

  call get_state('length',g,length)
  call get_state('status',g,istatus) ! Need status to check if we have spawned


  ! Update the length to the new length
  ! given time difference (in days) and
  ! previous length
  deltaT = g%DT_bio*sec_2_day
  where(istatus == ACTIVE)
     length = length*exp(0.0059 * deltaT)
  end where

  nullify(length)
  nullify(istatus)

end subroutine bio_mork

!----------------------------------------------------------------------
! Added by J. Ounsley
! Perform any changes to particle position due to biological movement
!
! Called in do_adv_diff which is called in main loop at time interval
! dT
! Implement swimming behaviour with constant bearing
!----------------------------------------------------------------------
subroutine advance_bio_movement(g,dT,time)
    type(igroup), intent(inout) :: g
    real(sp),     intent(in   ) :: dT,time

    integer :: i,np
    real(sp),    pointer :: x(:) !pointers allow for direct modification
    real(sp),    pointer :: y(:)
    integer,     pointer :: istatus(:)
    integer,     pointer :: cell(:)

    ! Temp variables for change
    real(sp), allocatable :: pdxt(:), pdyt(:)
    real(sp), allocatable :: pu(:), pv(:)
    real(sp)              :: speed

    !construct pointers to access and modify state variables for the group
    call get_state('x',g,x)
    call get_state('y',g,y)
    call get_state('status',g,istatus) ! Need status to check if we go out of water
    call get_state('cell',g,cell)

    np = g%nind

    allocate(pdxt(np))  ;
    allocate(pdyt(np))  ;
    allocate(pu(np))  ;
    allocate(pv(np))  ;

    pu = 0
    pv = 0
    !print *,x(1)

    ! Calculate velocity due to behaviour
    if(behaviour_type == 0) then
       call passive(g,pu,pv)
    else if (behaviour_type == 1) then
       call constant_bearing(g,pu,pv)
    else if (behaviour_type == 2) then
       call negative_rheotaxis(g,pu,pv)
    else if (behaviour_type == 3) then
       call positive_rheotaxis(g,pu,pv)
    else if (behaviour_type == 4) then
       call const_bearing_range(g,pu,pv)
    else if (behaviour_type == 5) then
       call mork(g,x,y,cell,istatus,pu,pv)
    else if (behaviour_type == 6) then
       call booker(g,x,y,cell,istatus,pu,pv)
    else if (behaviour_type == 7) then
       call salinity_gradient(g,x,y,cell,istatus,pu,pv)
    else if (behaviour_type == 8) then
       call salinity_const_bearing(g,x,y,cell,istatus,pu,pv)
    else if (behaviour_type == 9) then
       call salinity_const_bearing_2(g,x,y,cell,istatus,pu,pv)
    else if (behaviour_type == 10) then
       call depth_gradient(g,x,y,cell,istatus,pu,pv)
    else if (behaviour_type == 11) then
       call depth_const_bearing(g,x,y,cell,istatus,pu,pv)
    else if (behaviour_type == 12) then
       call booker_2(g,x,y,cell,istatus,pu,pv)
    else if (behaviour_type == 13) then
       call vattenfall(g,cell,pu,pv,time)
    else
       write (*,*) 'WARNING::::::: Unknown behaviour type ', behaviour_type, '. Behaviour ignored.'
    endif

    ! Update particle position with calculated velocity
    do i=1,np
      if(spherical==0)then
        pdxt(i) = x(i) + pu(i)*dT
        pdyt(i) = y(i) + pv(i)*dT
      elseif (spherical==1)then
          pdyt(i) = y(i)  + pv(i)*dt/(radius_earth*d2r)
          pdxt(i) = x(i)  + pu(i)*dt/(radius_earth*d2r*COS(pdyt(i)*d2r) + 1.0E-6)
      end if

    end do

    ! Check to see if we go out of the area
    call find_element(np,pdxt,pdyt,cell,istatus)
    !print *,x(1)

    ! If a particle exited the domain, try following the
    ! current instead, for non neg rheotaxis behaviour
    !print *,behaviour_type
    !print *,collision_avoidance
    if (behaviour_type /= 2 .and. collision_avoidance == 1) then
    	!print *,'Collision, performing rheotaxis'
   		call negative_rheotaxis(g,pu,pv)
    	do i=1,np
    		if (istatus(i)==EXITED) then
    			!print *,'Collision, performing rheotaxis'
			    if(spherical==0)then
			        pdxt(i) = x(i) + pu(i)*dT
			        pdyt(i) = y(i) + pv(i)*dT
			    elseif (spherical==1)then
			        pdyt(i) = y(i)  + pv(i)*dt/(radius_earth*d2r)
			        pdxt(i) = x(i)  + pu(i)*dt/(radius_earth*d2r*COS(pdyt(i)*d2r) + 1.0E-6)
			    end if
			end if
    	end do

    	! Check position again,
    	! hopefully still active
	    where(istatus==EXITED)
	      istatus=ACTIVE
	    end where

	    call find_element(np,pdxt,pdyt,cell,istatus)

    end if


    where (istatus==ACTIVE)
      x = pdxt
      y = pdyt
    end where

    !!--reset position of particles which are lost from domain to last known position
    where(istatus==EXITED)
      istatus=ACTIVE
    end where

    nullify(x)
    nullify(y)
    nullify(istatus)
    nullify(cell)

end subroutine advance_bio_movement


!----------------------------------------------------------------------
! Added by J. Ounsley
!
! No active movement, simply return 0 change in
! position
!----------------------------------------------------------------------
subroutine passive(g,pu,pv)
  type(igroup), intent(in) :: g
  real(sp), intent(out)    :: pu(:),pv(:) ! Swimming velocity

  pu(:) = 0
  pv(:) = 0

end subroutine! passive

!----------------------------------------------------------------------
! Added by J. Ounsley
!
! Calculate velocity due to swimming at constant bearing,
! at 1 body length per second, with given angle.!
!----------------------------------------------------------------------
subroutine constant_bearing(g,pu,pv)
	type(igroup), intent(in) :: g
	real(sp), intent(out)    :: pu(:),pv(:) ! Swimming velocity

	real(sp), pointer  :: u(:),v(:) ! Local forcing
	real(sp), pointer  :: length(:) ! Length needed for speed
	integer, pointer   :: cell(:) ! Length needed for speed
	!real(sp)           :: theta(g%nind)   ! For direction of current

	integer			 :: i

	call get_state('u',g,u)
	call get_state('v',g,v)
	call get_state('length',g,length)
	call get_state('cell',g,cell)

	!theta = atan2(u,(v + 1.0E-6))

	! If in boundary cell, then move with currents
	do i=1,g%nind
	 	if(isbce(cell(i))==1 .and. collision_avoidance==1)then
	 		! Compare the angle to the const bearing,
	 		! If they differ by more than 90 then go
	 		! opposite to current

	 		!print *,theta
	 		!if((theta(i)/d2r-const_bearing) >90 .or. (theta(i)/d2r-const_bearing) <-90)then
		 	!	theta(i) = theta(i) + pi
			!end if
			pu(i) = swim_speed * length(i) * u(i) / (sqrt(u(i)**2+v(i)**2) + 1.0E-6)
	  		pv(i) = swim_speed * length(i) * v(i) / (sqrt(u(i)**2+v(i)**2) + 1.0E-6)
			!pu(i) = swim_speed * length(i) * sin(theta(i))
	  		!pv(i) = swim_speed * length(i) * cos(theta(i))
		else
			pu(i) = swim_speed * length(i) * sin(const_bearing*d2r)
	  		pv(i) = swim_speed * length(i) * cos(const_bearing*d2r)
	  	end if
	end do

	! Swim at swim_speed bl per second,
	! at const_bearing degrees clockwise from north
	! pu(:) = swim_speed * length * sin(const_bearing*d2r)
	! pv(:) = swim_speed * length * cos(const_bearing*d2r)

	nullify(length)
	nullify(cell)
  	nullify(u)
  	nullify(v)

end subroutine! constant_bearing

!----------------------------------------------------------------------
! Added by J. Ounsley
!
! Swim with the local current at 1 body length per second
!----------------------------------------------------------------------
subroutine negative_rheotaxis(g,pu,pv)
  type(igroup), intent(in) :: g
  real(sp), intent(out)   :: pu(:),pv(:) ! Swimming velocity

  real(sp), pointer  :: u(:),v(:) ! Local forcing
  real(sp), pointer  :: length(:) ! Length needed for speed

  call get_state('u',g,u)
  call get_state('v',g,v)
  call get_state('length',g,length)



  ! swim_speed body lengths per second
  pu(:) = swim_speed * length * u / (sqrt(u**2+v**2) + 1.0E-6)
  pv(:) = swim_speed * length * v / (sqrt(u**2+v**2) + 1.0E-6)

  nullify(length)
  nullify(v)
  nullify(u)


end subroutine! negative_rheotaxis

!----------------------------------------------------------------------
! Added by J. Ounsley
!
! Swim against the local current at 1 body length per second
!----------------------------------------------------------------------
subroutine positive_rheotaxis(g,pu,pv)
  type(igroup), intent(in) :: g
  real(sp), intent(out)   :: pu(:),pv(:) !  Swimming velocity

  call negative_rheotaxis(g,pu,pv)
  pu = - pu
  pv = - pv

end subroutine! positive_rheotaxis

subroutine const_bearing_range(g,pu,pv)
	type(igroup), intent(in) :: g
	real(sp), intent(out)    :: pu(:),pv(:) ! Swimming velocity

	real(sp), pointer  :: u(:),v(:) ! Local forcing
	real(sp), pointer  :: length(:) ! Length needed for speed
	real(sp), pointer  :: bearing(:) ! Bearing relative to constant bearing
	integer, pointer   :: cell(:)

	integer			 :: i

	call get_state('length',g,length)
	call get_state('bearing',g,bearing)
	call get_state('cell',g,cell)
  	call get_state('u',g,u)
	call get_state('v',g,v)

  ! If in boundary cell, then move with currents
	do i=1,g%nind
	 	if(isbce(cell(i))==1 .and. collision_avoidance==1)then
			pu(i) = swim_speed * length(i) * u(i) / (sqrt(u(i)**2+v(i)**2) + 1.0E-6)
	  		pv(i) = swim_speed * length(i) * v(i) / (sqrt(u(i)**2+v(i)**2) + 1.0E-6)
		else
			pu(i) = swim_speed * length(i) * sin(bearing(i)*d2r)
	  		pv(i) = swim_speed * length(i) * cos(bearing(i)*d2r)
	  	end if
	end do

  	nullify(length)
  	nullify(bearing)
  	nullify(cell)
  	nullify(u)
  	nullify(v)

end subroutine! const_bearing_range

!--------------------------------------------
! Based on Mork et al. 2012
! Particles follow the current at normal speed but move to within
! preferred range of temperatures following gradient of local
! cell with additional linear increase in speed of movement
!
!--------------------------------------------
subroutine mork(g,x,y,cell,istatus,pu,pv)
  type(igroup), intent(in) :: g
  real(sp), intent(in)     :: x(:),y(:)
  integer, intent(in)      :: cell(:),istatus(:)
  real(sp), intent(out)    :: pu(:),pv(:)
  real(sp)                 :: dx_t(g%nind),dy_t(g%nind),dx_s(g%nind),dy_s(g%nind) ! Swimming direction
  real(sp)                 :: dv_t(g%nind),dv_s(g%nind) ! Gradient size

  real(sp), pointer  :: length(:) ! Length needed for speed
  real(sp), pointer  :: temp(:) ! Temperature at current location
  real(sp), pointer  :: salinity(:) ! Salinity at current location
  real(sp), pointer  :: s(:) ! Depth of particles
  real(sp), pointer  :: u(:),v(:) ! Local forcing
  integer			 :: i



  call get_state('length',g,length)
  call get_state('salinity',g,salinity)
  call get_state('T',g,temp)

  call get_state('s',g,s)

  !call get_state('cell',g,cell)
  call get_state('u',g,u)
  call get_state('v',g,v)


  ! Swim with the currents
  call negative_rheotaxis(g,pu,pv)


  !call gradient(g%Nind,x,y,s,cell,istatus,'temp',dx,dy,dv,4)
  call gradient(g%Nind,x,y,s,cell,istatus,'temp',dx_t,dy_t,dv_t,4)
  call gradient(g%Nind,x,y,s,cell,istatus,'salinity',dx_s,dy_s,dv_s,4)


  ! Move in direction of gradient, with additional speed
  ! linearly increasing with magnitude of temperature
  ! difference from preference up to max 1.
  where(temp < 8)
     pu = pu + .5 * swim_speed * length * dx_t * min(8-temp,1.0)
     pv = pv + .5 * swim_speed * length * dy_t * min(8-temp,1.0)
  elsewhere(temp > 14)
     pu = pu -.5 * swim_speed * length * dx_t * min(temp-14,1.0)
     pv = pv -.5 * swim_speed * length * dy_t * min(temp-14,1.0)
  elsewhere(salinity < 35)
     pu = pu + .5 * swim_speed * length * dx_s * min(35-salinity,1.0)
     pv = pv + .5 * swim_speed * length * dy_s * min(35-salinity,1.0)
  end where

  ! If in boundary cell, then move with currents
  do i=1,g%nind
     if(isbce(cell(i))==1 .and. collision_avoidance==1)then
	pu(i) = swim_speed * length(i) * u(i) / (sqrt(u(i)**2+v(i)**2) + 1.0E-6)
	pv(i) = swim_speed * length(i) * v(i) / (sqrt(u(i)**2+v(i)**2) + 1.0E-6)
     end if
  end do

  nullify(length)
  nullify(temp)
  nullify(s)
  nullify(salinity)

end subroutine mork

!--------------------------------------------
! Based on Booker et al. 2008
! Particles move to within preferred range
! of temperatures following gradient of local
! cell.
!
!--------------------------------------------
subroutine booker(g,x,y,cell,istatus,pu,pv)
  use utilities, only : ran_from_range
  type(igroup), intent(in) :: g
  real(sp), intent(in)     :: x(:),y(:)
  integer, intent(in)      :: cell(:),istatus(:)
  real(sp), intent(out)    :: pu(:),pv(:)
  real(sp)                 :: dx_t(g%nind),dy_t(g%nind) ! Swimming direction
  real(sp)                 :: dv_t(g%nind) ! Gradient size

  real(sp), pointer  :: length(:) ! Length needed for speed
  real(sp), pointer  :: temp(:) ! Temperature at current location
  real(sp), pointer  :: s(:) ! Depth of particles
  real(sp), pointer  :: u(:),v(:) ! Local forcing
  integer			 :: i

  call get_state('length',g,length)
  call get_state('T',g,temp)
  call get_state('s',g,s)

  !call get_state('cell',g,cell)
  call get_state('u',g,u)
  call get_state('v',g,v)


  !call gradient(g%Nind,x,y,s,cell,istatus,'temp',dx,dy,dv,4)
  call gradient(g%Nind,x,y,s,cell,istatus,'temp',dx_t,dy_t,dv_t,4)

  ! Move in direction of gradient, with speed
  ! linearly increasing with magnitude of temperature
  ! difference from preference up to max 1.
  where(temp < 4)
     pu = swim_speed * length  * dx_t
     pv = swim_speed * length  * dy_t
  elsewhere(temp > 8)
     pu = - swim_speed * length  * dx_t
     pv = - swim_speed * length  * dy_t
  elsewhere ! Random direction
     pu = ran_from_range(0.0_sp,360.0_sp)
     pv = swim_speed * length  * sin(pu*d2r)
     pu = swim_speed * length  * cos(pu*d2r)
  end where

  ! If in boundary cell, then move with currents
  do i=1,g%nind
     if(isbce(cell(i))==1 .and. collision_avoidance==1)then
	pu(i) = swim_speed * length(i) * u(i) / (sqrt(u(i)**2+v(i)**2) + 1.0E-6)
	pv(i) = swim_speed * length(i) * v(i) / (sqrt(u(i)**2+v(i)**2) + 1.0E-6)
     end if
  end do

  nullify(length)
  nullify(s)
  nullify(temp)

end subroutine booker

!--------------------------------------------
! Based on Booker et al. 2008
! Particles move to within preferred range
! of temperatures following gradient of local
! cell.
! Temperature range is modified from original
! to aim for 8-11 range
!--------------------------------------------
subroutine booker_2(g,x,y,cell,istatus,pu,pv)
  use utilities, only : ran_from_range
  type(igroup), intent(in) :: g
  real(sp), intent(in)     :: x(:),y(:)
  integer, intent(in)      :: cell(:),istatus(:)
  real(sp), intent(out)    :: pu(:),pv(:)
  real(sp)                 :: dx_t(g%nind),dy_t(g%nind) ! Swimming direction
  real(sp)                 :: dv_t(g%nind) ! Gradient size

  real(sp), pointer  :: length(:) ! Length needed for speed
  real(sp), pointer  :: temp(:) ! Temperature at current location
  real(sp), pointer  :: s(:) ! Depth of particles

  call get_state('length',g,length)
  call get_state('T',g,temp)
  call get_state('s',g,s)


  !call gradient(g%Nind,x,y,s,cell,istatus,'temp',dx,dy,dv,4)
  call gradient(g%Nind,x,y,s,cell,istatus,'temp',dx_t,dy_t,dv_t,4)

  ! Move in direction of gradient, with speed
  ! linearly increasing with magnitude of temperature
  ! difference from preference up to max 1.
  where(temp < 8)
     pu = swim_speed * length  * dx_t * max(1,int(8-temp))
     pv = swim_speed * length  * dy_t * max(1,int(8-temp))
  elsewhere(temp > 11)
     pu = - swim_speed * length  * dx_t * max(1,int(temp-11))
     pv = - swim_speed * length  * dy_t * max(1,int(temp-11))
  elsewhere ! Random direction
     pu = ran_from_range(0.0_sp,360.0_sp)
     pv = swim_speed * length  * sin(pu*d2r)
     pu = swim_speed * length  * cos(pu*d2r)
  end where

  nullify(length)
  nullify(s)
  nullify(temp)

end subroutine booker_2

subroutine salinity_gradient(g,x,y,cell,istatus,pu,pv)
  type(igroup), intent(in) :: g
  real(sp), intent(in)     :: x(:),y(:)
  integer, intent(in)      :: cell(:),istatus(:)
  real(sp), intent(out)    :: pu(:),pv(:)
  real(sp)                 :: dx_s(g%nind),dy_s(g%nind) ! Swimming direction
  real(sp)                 :: dv_s(g%nind) ! Gradient size

  real(sp), pointer  :: u(:) ! Local velocity
  real(sp), pointer  :: v(:) !
  real(sp), pointer  :: length(:) ! Length needed for speed
  real(sp), pointer  :: salinity(:) ! Salinity at current location
  real(sp), pointer  :: s(:) ! Depth of particles

  call get_state('u',g,u)
  call get_state('v',g,v)
  call get_state('length',g,length)
  call get_state('salinity',g,salinity)

  call get_state('s',g,s)

  call gradient(g%Nind,x,y,s,cell,istatus,'salinity',dx_s,dy_s,dv_s,4)


  ! Move in direction of gradient, with constant speed
  ! Otherwise follow currents
  where(salinity < 35)
     pu = swim_speed * length * dx_s
     pv = swim_speed * length * dy_s
  elsewhere
 ! swim_speed body lengths per second
     pu = swim_speed * length * u / (sqrt(u**2+v**2) + 1.0E-6)
     pv = swim_speed * length * v / (sqrt(u**2+v**2) + 1.0E-6)
  end where


  nullify(length)
  nullify(s)
  nullify(salinity)
  nullify(u)
  nullify(v)

end subroutine salinity_gradient


!
! Move in a constant bearing but follow salinity gradients
! when salinity is less than 34
! Additionally allow for collision avoidance, moving with current
! when in boundary cell
!
subroutine salinity_const_bearing(g,x,y,cell,istatus,pu,pv)
  type(igroup), intent(in) :: g
  real(sp), intent(in)     :: x(:),y(:)
  integer, intent(in)      :: istatus(:)
  real(sp), intent(out)    :: pu(:),pv(:)
  real(sp)                 :: dx_s(g%nind),dy_s(g%nind) ! Swimming direction
  real(sp)                 :: dv_s(g%nind) ! Gradient size

  real(sp), pointer  :: u(:) ! Local velocity
  real(sp), pointer  :: v(:) !
  real(sp), pointer  :: bearing(:) !
  real(sp), pointer  :: length(:) ! Length needed for speed
  real(sp), pointer  :: salinity(:) ! Salinity at current location
  real(sp), pointer  :: s(:) ! Depth of particles
  integer, pointer   :: cell(:) !

  call get_state('u',g,u)
  call get_state('v',g,v)
  call get_state('length',g,length)
  call get_state('salinity',g,salinity)
  call get_state('bearing',g,bearing)
  call get_state('s',g,s)


  call gradient(g%Nind,x,y,s,cell,istatus,'salinity',dx_s,dy_s,dv_s,4)


  ! Move in direction of salinity gradient, with constant speed
  ! Otherwise, if on boundary follow current otherwise swim
  ! in given direction
  where(salinity < 34.0)
     pu = swim_speed * length * dx_s
     pv = swim_speed * length * dy_s
  elsewhere(isbce(cell)==1 .and. collision_avoidance==1)
     ! swim_speed body lengths per second
     pu = swim_speed * length * u / (sqrt(u**2+v**2) + 1.0E-6)
     pv = swim_speed * length * v / (sqrt(u**2+v**2) + 1.0E-6)
  elsewhere
  	 pu = swim_speed * length * sin(bearing*d2r)
	 pv = swim_speed * length * cos(bearing*d2r)
  end where


  nullify(length)
  nullify(s)
  nullify(salinity)
  nullify(bearing)
  nullify(u)
  nullify(v)

end subroutine salinity_const_bearing

!
! Move in weighted combination  of  constant bearing and
! salinity gradient with exponential weighting with threshold of
! salinity influence at 34
!
subroutine salinity_const_bearing_2(g,x,y,cell,istatus,pu,pv)
  type(igroup), intent(in) :: g
  real(sp), intent(in)     :: x(:),y(:)
  integer, intent(in)      :: istatus(:)
  real(sp), intent(out)    :: pu(:),pv(:)
  real(sp)                 :: dx_s(g%nind),dy_s(g%nind) ! Swimming direction
  real(sp)                 :: dv_s(g%nind) ! Gradient size
  real(sp)                 :: weight(g%nind) ! Weighting on salinity

  real(sp), pointer  :: u(:) ! Local velocity
  real(sp), pointer  :: v(:) !
  real(sp), pointer  :: bearing(:) !
  real(sp), pointer  :: length(:) ! Length needed for speed
  real(sp), pointer  :: salinity(:) ! Salinity at current location
  real(sp), pointer  :: s(:) ! Depth of particles
  integer, pointer   :: cell(:) !
  integer			 :: i

  call get_state('u',g,u)
  call get_state('v',g,v)
  call get_state('length',g,length)
  call get_state('salinity',g,salinity)
  call get_state('bearing',g,bearing)
  call get_state('s',g,s)


  call gradient(g%Nind,x,y,s,cell,istatus,'salinity',dx_s,dy_s,dv_s,4)

	do i=1,g%nind

		weight(i) = min(0.1*exp(34-salinity(i)),1.0)
		!weight(i) = 0
		pu(i) = swim_speed * length(i) * (sin(bearing(i)*d2r) * (1-weight(i)) + weight(i)*dx_s(i))
  		pv(i) = swim_speed * length(i) * (cos(bearing(i)*d2r) * (1-weight(i)) + weight(i)*dy_s(i))
	end do

	!print * ,weight(1)
	!print * ,weight(2)
  	! Move in direction of salinity gradient, and bearing
  	! weighted by salinity less than 34
	!pu = swim_speed*length*( sin(bearing*d2r) * (1-weight) + weight*dx_s)
	!pv = swim_speed*length*( cos(bearing*d2r) * (1-weight) + weight*dy_s)


  nullify(length)
  nullify(s)
  nullify(salinity)
  nullify(bearing)
  nullify(u)
  nullify(v)

end subroutine salinity_const_bearing_2


!
! Move along depth gradient
!
subroutine depth_gradient(g,x,y,cell,istatus,pu,pv)
  type(igroup), intent(in) :: g
  real(sp), intent(in)     :: x(:),y(:)
  integer, intent(in)      :: istatus(:)
  real(sp), intent(out)    :: pu(:),pv(:)
  real(sp)                 :: dx_s(g%nind),dy_s(g%nind) ! Swimming direction
  real(sp)                 :: dv_s(g%nind) ! Gradient size

  real(sp), pointer  :: u(:) ! Local velocity
  real(sp), pointer  :: v(:) !
  real(sp), pointer  :: bearing(:) !
  real(sp), pointer  :: length(:) ! Length needed for speed
  real(sp), pointer  :: h(:) ! Bathymetry at current location positive down (m)
  real(sp), pointer  :: s(:) ! Depth of particles
  integer, pointer   :: cell(:) !

  call get_state('u',g,u)
  call get_state('v',g,v)
  call get_state('length',g,length)
  call get_state('h',g,h)
  call get_state('bearing',g,bearing)
  call get_state('s',g,s)


  call gradient(g%Nind,x,y,cell,istatus,'h',dx_s,dy_s,dv_s,4)


  pu = swim_speed * length * dx_s
  pv = swim_speed * length * dy_s

  nullify(length)
  nullify(s)
  nullify(h)
  nullify(bearing)
  nullify(u)
  nullify(v)

end subroutine depth_gradient


!
! Move in weighted linear combination of constant bearing and
! depth gradient. Depth gradient only influential upto 100 m
!
subroutine depth_const_bearing(g,x,y,cell,istatus,pu,pv)
  type(igroup), intent(in) :: g
  real(sp), intent(in)     :: x(:),y(:)
  integer, intent(in)      :: istatus(:)
  real(sp), intent(out)    :: pu(:),pv(:)
  real(sp)                 :: dx_s(g%nind),dy_s(g%nind) ! Swimming direction
  real(sp)                 :: dv_s(g%nind) ! Gradient size
  real(sp)                 :: weight(g%nind) ! Weighting on salinity

  real(sp), pointer  :: u(:) ! Local velocity
  real(sp), pointer  :: v(:) !
  real(sp), pointer  :: bearing(:) !
  real(sp), pointer  :: length(:) ! Length needed for speed
  real(sp), pointer  :: h(:) ! Bathymetry at current location positive down (m)
  real(sp), pointer  :: s(:) ! Depth of particles
  integer, pointer   :: cell(:) !
  integer			 :: i

  call get_state('u',g,u)
  call get_state('v',g,v)
  call get_state('length',g,length)
  call get_state('h',g,h)
  call get_state('bearing',g,bearing)
  call get_state('s',g,s)


  call gradient(g%Nind,x,y,cell,istatus,'h',dx_s,dy_s,dv_s,4)


  weight = max((100-h)/100,0.0)

  pu = swim_speed * length * (sin(bearing*d2r) * (1-weight) + weight*dx_s)
  pv = swim_speed * length * (cos(bearing*d2r) * (1-weight) + weight*dy_s)

  ! When at depths less than 5 meters, move down gradient
  !where(h < 100.0)
  !   pu = swim_speed * length * dx_s
  !   pv = swim_speed * length * dy_s
  !elsewhere
  !   ! Otherwise move in given bearing
  !   pu = swim_speed * length * sin(bearing*d2r)
  !   pv = swim_speed * length * cos(bearing*d2r)
  !end where

  nullify(length)
  nullify(s)
  nullify(h)
  nullify(bearing)
  nullify(u)
  nullify(v)

end subroutine depth_const_bearing

! Vattenfall simulations
! Look at behavioural transition from a
! SE to N after a given time period
subroutine vattenfall(g,cell,pu,pv,time)
  use utilities, only : ran_from_range
	type(igroup), intent(in) :: g
	real(sp), intent(out)    :: pu(:),pv(:) ! Swimming velocity
	real(sp), intent(in)     :: time

	real(sp), pointer  :: u(:),v(:) ! Local forcing
	real(sp), pointer  :: length(:) ! Length needed for speed
	real(sp), pointer  :: bearing(:) ! Bearing relative to constant bearing
	real(sp), pointer  :: tspawn(:) ! spawn time to change direction
        integer, pointer   :: stage(:)
	integer, pointer   :: cell(:) !

	integer			 :: i

	call get_state('length',g,length)
	call get_state('bearing',g,bearing)
  	call get_state('u',g,u)
	call get_state('v',g,v)
        call get_state('stage',g,stage)
	call get_state('tspawn',g,tspawn)


	! Check if bearing should change to N
	do i=1,g%nind
		!print *,tspawn(i)
		!print *,stage(i)
		if (stage(i) == 0 .and. (time - tspawn(i)) > transition_time)then
			stage(i) = 1
			! Change bearing to N
			! plus variability
        		bearing(i) = 0 + ran_from_range(-bearing_range/2,bearing_range/2)
		end if
	end do

	! If in boundary cell, then move with currents
	do i=1,g%nind
	 	if(isbce(cell(i))==1 .and. collision_avoidance==1)then
			pu(i) = swim_speed * length(i) * u(i) / (sqrt(u(i)**2+v(i)**2) + 1.0E-6)
	  		pv(i) = swim_speed * length(i) * v(i) / (sqrt(u(i)**2+v(i)**2) + 1.0E-6)
		else
			pu(i) = swim_speed * length(i) * sin(bearing(i)*d2r)
	  		pv(i) = swim_speed * length(i) * cos(bearing(i)*d2r)
	  	end if
	end do

  	nullify(length)
  	nullify(bearing)
  	nullify(u)
  	nullify(v)
	nullify(stage)
	nullify(tspawn)

end subroutine! vattenfall



End Module Bio
