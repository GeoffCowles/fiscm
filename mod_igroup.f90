!=======================================================================
! Fiscm IGroup Type
! Copyright:    2008(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Comments:     FISCM - Groups of Like Individuals 
! 
!=======================================================================

Module mod_igroup

Use gparms
Use mod_pvar

Implicit None

!Group type
type igroup 
  integer             :: id
  integer             :: nstate
  integer             :: problem_dimension
  character(len=fstr) :: group_name
  integer             :: Tnind
  integer             :: nind
  real(sp)            :: hdiff_const_val
  real(sp)            :: vdiff_const_val
  integer             :: hdiff_type
  integer             :: vdiff_type
  real(sp)            :: DT_bio       !time step for biology
  logical             :: biology
  real(sp)            :: intvl_out
  real(sp)            :: tlast_out
  real(sp)            :: start_out
  integer             :: frame_out
  integer             :: fid_out
  character(len=fstr) :: fname_out
  type(pvar_list)     :: state
end type igroup 

interface add_state
  module procedure add_state_fvec
  module procedure add_state_ivec
end interface
  
interface get_state
  module procedure get_state_fvec
  module procedure get_state_ivec
end interface
contains
!---------------------------------------------------------
!Create the initial group type (passive particle)
!---------------------------------------------------------
function group_(i,id) result (g)
  type(igroup) :: g
  integer, intent(in) :: i,id

  if(i < 1)then
	write(*,*) 'error creating group: size must be greater than 0'
	stop 
  endif

  !initial conditions for scalar components of group type, should initialize within type definition, but
  !this is not allowed on some older compilers
  g%id    = id
  g%Tnind = i
  g%nind  = i      !gwc?
  g%DT_bio = 3600. !gwc?
  g%state = pvar_list_()
  g%nstate = 0
  g%group_name = ""
  g%problem_dimension = 0
  g%hdiff_type = HDIFF_NONE
  g%vdiff_type = HDIFF_NONE
  g%frame_out  = 0
  g%start_out  = 0.0
  g%intvl_out  = 7200.
  g%tlast_out  = g%start_out - g%intvl_out

  !(x) - x location of particle in the domain
  call add_state(g,'x','x location of particle','m',NETCDF_YES,1.0)

  !(y) - y location of particle in the domain
  call add_state(g,'y','y location of particle','m',NETCDF_YES,10.0_sp)
  
  !(status) - particle status
  call add_state(g,'status','particle status','-',NETCDF_YES,1)

  !(pathlength) - integrated path length
  call add_state(g,'pathlength','integrated trajectory','m',NETCDF_YES,0.0)

end function group_

subroutine add_state_fvec(g,varname,longname,units,output,init_val)
   type(igroup)        :: g 
   character(len=*)    :: varname
   character(len=*)    :: longname
   character(len=*)    :: units
   integer             :: output
   real(sp)            :: init_val
   type(pvar)          :: new
   integer             :: i

   !need to make sure we arent duplicating state names!
   allocate(new%fvar(g%Tnind)) 
   do i=1,g%Tnind
	 new%fvar(i) = init_val
   end do
    
   !check output flag
   if(output < 0 .or. output > 1)then
	 write(*,*)'error adding state variable: ',varname
	 write(*,*)'output argument must be 0 (no output), or 1 (output)'
     stop
   endif
 
   new%output   = output
   new%istype   = ftype
   new%idim     = g%Tnind
   new%varname  = trim(varname)
   new%longname = trim(longname)
   new%units    = trim(units)
   new%output   = output

   call add_pvar_node(g%state,new)
   g%nstate = g%nstate + 1;
end subroutine add_state_fvec

subroutine add_state_ivec(g,varname,longname,units,output,init_val)
   type(igroup)        :: g 
   character(len=*)    :: varname
   character(len=*)    :: longname
   character(len=*)    :: units
   integer             :: output
   integer             :: init_val
   type(pvar)          :: new
   integer             :: i

   !need to make sure we arent duplicating state names!
 
   allocate(new%ivar(g%Tnind)) 
   
   do i=1,g%Tnind
	 new%ivar(i) = init_val
   end do
  
   !check output flag
   if(output < 0 .or. output > 1)then
	 write(*,*)'error adding state variable: ',varname
	 write(*,*)'output argument must be 0 (no output), or 1 (output)'
     stop
   endif

   new%output   = output 
   new%idim     = g%Tnind
   new%varname  = trim(varname)
   new%longname = trim(longname)
   new%units    = trim(units)
   new%istype   = itype
   call add_pvar_node(g%state,new)
   g%nstate = g%nstate + 1;
end subroutine add_state_ivec

!function get_state_fvec(varname,g) result(v)
subroutine get_state_fvec(varname,g,v)
   character(len=*)    :: varname
   type(igroup)        :: g
   real(sp), pointer   :: v(:)
   type(pvar), pointer :: p
   
   p => get_pvar(g%state,varname)
   v => p%fvar
   
end subroutine get_state_fvec
!end function get_state_fvec

!function get_state_ivec(varname,g) result(v)
subroutine get_state_ivec(varname,g,v)
   character(len=*)    :: varname
   type(igroup)        :: g
   integer, pointer    :: v(:)
   type(pvar), pointer :: p
   
   p => get_pvar(g%state,varname)
   v => p%ivar
   
end subroutine get_state_ivec
!end function get_state_ivec


subroutine print_group_summary(g) 
   type(igroup)        :: g

   !group parameters
   write(*,*)
   write(*,*)'==========================Group Summary=========================='
   write(*,'(A25,I10)')'total individuals:       ',g%Tnind
   write(*,'(A25,I10)')'allocated individuals:   ',g%nind
   call print_state_vars(g%state)

   
end subroutine print_group_summary

End Module mod_igroup
