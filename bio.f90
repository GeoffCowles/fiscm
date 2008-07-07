!=======================================================================
! Fiscm Biology Example (Homerus americanus from Incze et al.)
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

Use gparms
Use mod_igroup

Implicit none

!parameters
integer,parameter  :: nstages = 4
real(sp) :: fcomp_settle = 0.5
real(sp) :: SCF = 1.0         


!namelist can be read from input to overide parameters
  Namelist /NML_LOBSTER/ &
     & fcomp_settle,     &
     & SCF

contains 

!------------------------------------------------------------------
! initialize a group
!   - read the control namelist (optional)
!   - add state variables and initialize them
!   - modify state variables (position, spawn time, initial weight, etc.)
!------------------------------------------------------------------
subroutine init_bio(g,Nind_start)
  use output_routines
  type(igroup), intent(inout) :: g
  integer,      intent(in)    :: Nind_start
  integer :: iunit,ios
  logical :: fexist

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
    write(*,*)'fatal error: namelist file: fiscm.nml does not exist, stopping...'
    stop 
  endif
  write(*,*)'initializing namelist using paramfile',g%paramfile
  open(unit=iunit,file=trim(g%paramfile),form='formatted')
  read(unit=iunit,nml=nml_lobster,iostat=ios)
  if(ios /= 0)then
    write(*,*)'fatal error: could not read lobster namelist from: ',trim(g%paramfile)
    stop
  endif
  endif ! file /= NONE
  

  !add parameters to netcdf file header 
  call add_cdfglobs(g,"info","some kind of info")
  call add_cdfglobs(g,"nstages",4)
  call add_cdfglobs(g,"fcomp_settle",fcomp_settle)
  call add_cdfglobs(g,"SCF",SCF)

end subroutine init_bio 

!------------------------------------------------------------------
! advance the biology (this routine is called from the main loop at
!                      a time interval of DT_bio )
!------------------------------------------------------------------
subroutine advance_bio(g,mtime)
  type(igroup), intent(inout) :: g
  real(sp),     intent(in   ) :: mtime
  real(sp),     pointer :: PASD(:)
  real(sp),     pointer :: T(:)
  integer ,     pointer :: stage(:)
  integer ,     pointer :: istatus(:)
  integer               :: i,N
  real(sp)              :: deltaT,D

  !construct pointers to access and modify state variables for the group
  call get_state('PASD',g,PASD)
  call get_state('T',g,T)
  call get_state('status',g,istatus)
  call get_state('stage',g,stage)
  
  !set problem size
  N = g%nind
  deltaT = g%DT_bio

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
      write(*,*)'stage 1 duration: ',D,deltaT,PASD(i)
    case(2)
      D = 200*(T(i)-4.88)**(-1.47)
      write(*,*)'stage 2 duration: ',D,deltaT
    case(3)
      D = 252*(T(i)-5.30)**(1.45)
    case(4)
      D = .3583*T(i)**2 - 14.316*T(i) + 156.895
    case(5)
      D = 0.0
    case default
      write(*,*)'stage: ',stage(i),' not a valid stage for a Homerus'
      stop
    end select

    PASD(i) = PASD(i) + SCF*deltaT/D

    !settle the post-larvae 
    if(stage(i)==4 .and. PASD(i) > fcomp_settle)then
       istatus(i) = SETTLED 
       stage(i) = 5
    else
      if(PASD(i) > 1.0)then
        stage(i) = stage(i) + 1
        PASD(i)  = 0.0
      endif
    endif  
 
  end do !end main loop


end subroutine advance_bio 


End Module Bio 
