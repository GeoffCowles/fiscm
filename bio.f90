!=======================================================================
! Fiscm Biology Example
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

!lobster model data - could put in a namelist
integer,  parameter :: nstages = 4 
real(sp), parameter :: fcomp_settle = 0.5
real(sp), parameter :: SCF = 1.0   !stage duration correction factor

contains

subroutine init_bio(g,Nind)
  type(igroup), intent(inout) :: g
  integer,      intent(in)    :: Nind

  !set current problem size
  if(Nind < 0)then
	write(*,*)'must specify more than 0 individuals in init_bio'
	stop
  else if(Nind > g%Tnind)then
	write(*,*)'number of individuals in init_bio (',Nind,') exceeds problem size',g%Tnind
	stop
  endif
  g%Nind = Nind

  !add state vectors
  g%Nind = Nind
  call add_state(g,'stage','lobster_stage','-',NETCDF_YES,1)
  call add_state(g,'PASD','currency','-',NETCDF_YES,0.0)
  call add_state(g,'T','temperature','C',NETCDF_YES,15.0)
  
end subroutine init_bio

subroutine advance_bio(g,mtime)
  type(igroup), intent(inout) :: g
  real(sp),     intent(in   ) :: mtime
  real(sp),     pointer :: PASD(:)
  real(sp),     pointer :: T(:)
  integer ,     pointer :: stage(:)
  integer ,     pointer :: status(:)
  integer               :: i,N
  real(sp)              :: deltaT,D

  call get_state('PASD',g,PASD)
  call get_state('T',g,T)
  call get_state('status',g,status)
  call get_state('stage',g,stage)
  
  !set problem size
  N = g%nind
  deltaT = g%DT_bio

!  following interface is better but compiler cannot construct
!  generic interface if arguments are the same even if the
!  return value is of different type
!  PASD   => get_state('PASD',g)
!  status => get_state('status',g)
!  stage  => get_state('stage',g)

  do i=1,N !main loop
	
    if(status(i) /= ACTIVE)cycle
    
    !update PASD using stage-based Duration
    select case (stage(i))
    case(1)
      D = 851.*(T(i)-0.84)**-1.91
    case(2)
      D = 200*(T(i)-4.88)**-1.47
    case(3)
      D = 252*(T(i)-5.30)**1.45
    case(4)
      D = .3583*T(i)**2 - 14.316*T(i) + 156.895
    case(5)
      D = 0.0
    case default
      write(*,*)'stage: ',stage(i),' not a valid stage for a lobster'
      stop
    end select

    PASD(i) = PASD(i) + SCF*(deltaT/86400)*D

    !settle the post-larvae 
    if(stage(i)==4 .and. PASD(i) > fcomp_settle)then
	  status(i) = SETTLED 
	  stage(i) = 5
	else
	  if(PASD(i) > 1.0)then
	    stage(i) = stage(i) + 1
	    PASD(i)  = 0.0
	  endif
	endif  
 
  end do !end main loop

  !debug
  !write(*,*)mtime/(3600*24),stage(1)

end subroutine advance_bio


End Module Bio