!=======================================================================
! Fiscm Main 
!
! Description
!   Main routine for FISCM - Setup and Time Looper 
!    
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!=======================================================================
Module fiscm_data
  use gparms
  use mod_igroup
  use forcing
  implicit none
  type(igroup), allocatable :: igroups(:)
  !character(len=fstr) :: runcontrol 
  real(sp) :: beg_time,end_time
  real(sp) :: beg_time_days,end_time_days
  real(sp) :: fiscm_time,fiscm_time_days,tshift
  real(sp) :: fbeg,fend
  real(sp) :: tsmin,tsmax
  real(sp) :: deltaT
  integer  :: Ireport 
  integer  :: sim_direction
  integer  :: ngroups
  integer  :: nfiles_in
  character(len=fstr) :: forcing_file(max_nf)
  character(len=fstr) :: output_prefix
  
  Namelist /NML_FISCM/ &
      beg_time_days,   & 
      end_time_days,   &
      deltaT,          &
      ireport,         &
      Ngroups,         &
      Nfiles_in,       &     
      forcing_file,    &
      output_prefix,   &
      spherical,       &
      sz_cor,          &
      fix_dep,         &
      dvm_bio,         &
     wind_type,        &
     dvmh_up,          &
     dvmh_dn  
End Module fiscm_data

Program fiscm
  use gparms 
  use fiscm_data
  use mod_igroup
  use mod_driver
  use bio
  use output_routines
  use forcing
  use utilities
  implicit none
  integer :: i,n,its,nits
  real(sp):: t = 0
  integer :: nvars ,uniqvars
  character(len=fstr) :: needvars(max_state_vars)

  !===================================================================
  !initialize the model (read namelists, allocate data, set ic's
  !===================================================================

  !----------------------------------------------------------
  ! read primary simulation data and setup groups
  !----------------------------------------------------------
  fbeg = -1.
  fend = -1.
  call setup

  !----------------------------------------------------------
  ! initialize biology 
  !  - set spawning times
  !  - set initial locations
  !  - set initial conditions for other states (e.g. length)
  !----------------------------------------------------------
  do n=1,ngroups
    if(igroups(n)%biology)call init_bio(igroups(n),igroups(n)%Tnind)
  end do

  !----------------------------------------------------------
  ! define variables for output
  !----------------------------------------------------------
  call cdf_out(ngroups,igroups,0,t,NCDO_ADD_STATES)

  !----------------------------------------------------------
  ! setup connection with forcing 
  !   determine which forcing vars needed from states
  !   determine which forcing vars needed for adv/diffusion 
  !   read in the mesh [ocean model dependent]
  !   setup frames to store data read from netcdf
  !----------------------------------------------------------
  nvars = 0
  call get_ext_varnames(ngroups,igroups,nvars,needvars)
  call ocean_model_init(ngroups,igroups,nvars,needvars)
  call get_unique_strings(nvars,needvars,uniqvars);nvars = uniqvars

  if(nvars > 0 .and. maxval(igroups%space_dim) > 1)call setup_forcing(nvars,needvars)

  !----------------------------------------------------------
  ! find elements containing particles 
  !----------------------------------------------------------
  call update_element(ngroups,igroups) 
  if(.not. checkstatus(ngroups,igroups,t))then
    write(*,*)'all dead'
    stop
  endif

  !----------------------------------------------------------
  ! summarize all states to the screen 
  !----------------------------------------------------------
  do n=1,ngroups
    call print_group_summary(igroups(n))
  end do

  !----------------------------------------------------------
  ! ensure times (spawning, simulation, forcing)
  ! are compatible
  !----------------------------------------------------------
   tsmin=beg_time
   tsmax=end_time
  call get_spawn_range(ngroups,igroups,tsmin,tsmax)
   beg_time=FLOOR(tsmin)
  call check_times(beg_time,end_time,tsmin,tsmax,fbeg,fend)

  !===================================================================
  ! => begin main loop over time 
  !===================================================================
  call drawline("-") ; write(*,*)'Beginning Sim' ; call drawline("-")

  t = beg_time
  nits = day_2_sec*(end_time-beg_time)/deltaT
    !---------------------------------------------------------
   if(maxval(igroups%space_dim) > 1) then
    call update_forcing(t,3)
    call sz_ini(ngroups,igroups)
    call activate(ngroups,igroups,t,sim_direction)
    call interp_forcing(ngroups,igroups,3) 
    call cdf_out(ngroups,igroups,its,t,NCDO_OUTPUT)
    call do_bio(ngroups,igroups,t,its)
    !---------------------------------------------------------
    !call exchange_forcing
    endif
  do its=1,nits
    t = t + deltaT*sec_2_day

    !  write(*,'(I10,F10.2)')its,t
    !---------------------------------------------------------
    ! update forcing to time t    
    !---------------------------------------------------------
    if(maxval(igroups%space_dim) > 1) call update_forcing(t,4) 

    !---------------------------------------------------------
    ! check for activation through spawning 
    !---------------------------------------------------------
    call activate(ngroups,igroups,t,sim_direction) 
    !---------------------------------------------------------
    ! advect and diffuse particles 
    !---------------------------------------------------------
    call do_adv_diff(ngroups,igroups,deltaT,t)
    !---------------------------------------------------------
    ! interpolate external vars at new particle positions 
    !---------------------------------------------------------
    call interp_forcing(ngroups,igroups,4) !4 is new current field

    !---------------------------------------------------------
    ! output the particle states to a NetCDF file 
    !---------------------------------------------------------
    call cdf_out(ngroups,igroups,its,t,NCDO_OUTPUT)

    !---------------------------------------------------------
    ! report progress to screen 
    !---------------------------------------------------------
    if(mod(its-1,Ireport)==0)then
      if(.not. checkstatus(ngroups,igroups,t))exit
    endif


    !---------------------------------------------------------
    ! advance the biology in time 
    !---------------------------------------------------------
    call do_bio(ngroups,igroups,t,its)
    !---------------------------------------------------------
    ! update the time step 
    !---------------------------------------------------------
    !t = t + deltaT
    call exchange_forcing

  end do
  !===================================================================
  ! <=  end main loop over time 
  !===================================================================

  !---------------------------------------------------------
  ! cleanup data
  !---------------------------------------------------------
  call drawline("-") ; write(*,*)'Finalizing Sim' ; call drawline("-")
  deallocate(igroups)

End Program fiscm

Subroutine setup
  use gparms
  use fiscm_data
  use mod_igroup
  use output_routines
  use utilities
  implicit none
  logical :: fexist
  integer, parameter :: iunit = 33
  integer :: n,ios
  real(sp) :: uno = 1.0
  real(sp) :: t   = 0.0
  character(len=fstr) :: buffer

  !--------------------------------------------------------
  ! get primary control file name from the command line 
  !--------------------------------------------------------
  n = iargc()
  if(n /= 1)then
    write(*,*)'incorrect usage: '
    write(*,*)'Usage:'
    write(*,*)'   fiscm file.nml'
    write(*,*)'where file.nml is the main parameter file'
    write(*,*)'stopping'
    stop
  endif
  call getarg(1,buffer)
  read(buffer,*) runcontrol 

  !--------------------------------------------------------
  ! check for existence of primary control file
  !--------------------------------------------------------
  call drawline("-")
  write(*,*)'reading main parameters from: ',trim(runcontrol)
  call drawline("-")
  inquire(file=trim(runcontrol),exist=fexist)
  if(.not.fexist)then
    write(*,*)'fatal error: namelist file: ',trim(runcontrol),' does not exist, stopping...'
    stop 
  endif

  !--------------------------------------------------------
  ! open and read global namelist:  nml_fiscm 
  !--------------------------------------------------------
  output_prefix = 'fiscm'
  open(unit=iunit,file=trim(runcontrol),form='formatted')
  read(unit=iunit,nml=nml_fiscm,iostat=ios)
  if(ios /= 0)then
    write(*,*)'fatal error: could not read fiscm namelist from',trim(runcontrol)
    stop
  endif

  !set begin/end time 
  beg_time = beg_time_days
  end_time = end_time_days

  !sanity on number of gruops 
  if(ngroups < 1)then
    write(*,*)'Fatal error: number of groups < 1' ; stop 
  endif

  !set simulation direction based on order of begin/end times
  if(end_time /= beg_time)then
    sim_direction = sign(uno, end_time-beg_time)
  else
    write(*,*)'fatal error: begin and end time are identical'
    stop
  endif

  !check for existence and open forcing file (if needed) 
  !if(trim(forcing_file) /= 'NONE' .and. trim(forcing_file) /= 'noe')then
   if(nfiles_in > 0) then
    fbeg = -1.0
    fend = -1.0
!    call open_forcing_file(trim(forcing_file),fbeg,fend) 
     call open_forcing_file(forcing_file,nfiles_in,fbeg,fend)
  endif

  !report
  write(*,*)'begin time:  ',beg_time
  write(*,*)'end time  :  ',end_time
  if(sim_direction ==  1)  write(*,*)'direction:      forward'
  if(sim_direction == -1)  write(*,*)'direction:      backward'
  write(*,*)'time step(s) ',deltaT
  write(*,*)'num groups:  ',ngroups
  write(*,*)'Ireport   :  ',Ireport
  call drawline("-")

  !---------------------------------------------------
  !read and allocate individual groups
  !---------------------------------------------------
  allocate(igroups(ngroups))
  do n=1,ngroups
    igroups(n) = group_(iunit,n,deltaT,output_prefix)
  end do
  close(iunit)

  !---------------------------------------------------
  !add state variables
  !---------------------------------------------------
  do n=1,ngroups
    call group_addstates(igroups(n))  
  end do

  !open output files and assign netcdf ids to each group
  call cdf_out(ngroups,igroups,0,t,NCDO_HEADER)
  call drawline("-")
  write(*,*)'setup complete'
  call drawline("-")
  
End Subroutine setup


!------------------------------------------------------------------
! check times 
!------------------------------------------------------------------
Subroutine check_times(mbeg,mend,smin,smax,fbeg,fend) 
  use gparms
  use utilities
  implicit none
  real(sp), intent(in) :: mbeg
  real(sp), intent(in) :: mend 
  real(sp), intent(in) :: smin
  real(sp), intent(in) :: smax
  real(sp), intent(in) :: fbeg
  real(sp), intent(in) :: fend 
  !-----------------------------
  logical istop
  real(sp) :: tmin,tmax

  call drawline('-')
  write(*,*)'checking coherency between: '
  write(*,*)'   spawning times '
  write(*,*)'   simulation begin/end time'
  write(*,*)'   model forcing begin/end time (if applicable)'
  call drawline('.')
  write(*,*)'time of latest   spawning: ',smax 
  write(*,*)'time of earliest spawning: ',smin
  if(fbeg >= 0.0 .and. fend >= 0.0)then
    write(*,*)'netcdf forcing start time: ',fbeg
    write(*,*)'netcdf forcing end time  : ',fend
  endif
  write(*,*)'model simulation beg time: ',mbeg
  write(*,*)'model simulation end time: ',mend

  istop = .false.
  tmin = min(mend,mbeg)
  tmax = max(mend,mbeg)

  !make sure spawning range is contained inside the simulation range
  if(smin < tmin .or. smin > tmax .or. smax < tmin .or. smax > tmax)then
    write(*,*)'fatal error: spawning not contained within simulation bounds' 
    istop = .true. 
  endif
  !make sure simulation range is contained within forcing data time range 
  !(if forcing active)
  if(fbeg >= 0.0 .and. fend >= 0.0)then
    if(tmin < fbeg .or. tmin > fend .or. tmax < fbeg .or. tmax > fend)then
      write(*,*)'fatal error: simulation not contained within forcing bounds' 
     istop = .true. 
    endif
  endif
  if(istop) stop

End Subroutine check_times

