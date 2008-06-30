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
!=======================================================================

Module mod_igroup

Use gparms
Use mod_pvar

Implicit None

!Group type
type igroup 
  integer             :: nstate
  character(len=fstr) :: species
  integer             :: Tnind
  integer             :: nind
  real(sp)            :: hdiff
  logical             :: biology
  type(pvar_list)     :: state
end type igroup 
  
contains
!---------------------------------------------------------
!Create the initial group type (passive particle)
!---------------------------------------------------------
function group_(i) result (g)
  type(igroup) :: g
  integer, intent(in) :: i

  if(i < 1)then
	write(*,*) 'error creating group: size must be greater than 0'
	stop 
  endif

  g%Tnind = i
  g%nind  = 0
  g%state = pvar_list_()
  g%nstate = 0

  !(x) - x location of particle in the domain
  call add_state(g,'x','x location of particle','m',NETCDF_YES,1.0)
  call add_state(g,'y','y location of particle','m',NETCDF_YES,10.0_sp)
  write(*,*)'state added?'
  

end function group_

subroutine add_state(g,varname,longname,units,output,init_val)
   type(igroup)        :: g 
   character(len=*)    :: varname
   character(len=*)    :: longname
   character(len=*)    :: units
   integer             :: output
   real(sp)            :: init_val
   type(pvar)          :: new
   integer             :: i

   write(*,*)'adding the state: ',trim(varname),g%Tnind
   !need to make sure we arent duplicating state names!
   allocate(new%fvar(g%Tnind)) 
   do i=1,g%Tnind
	 new%fvar(i) = float(i)*init_val
   end do
   write(*,*)'allocated'
   new%idim     = g%Tnind
   new%varname  = trim(varname)
   new%longname = trim(longname)
   new%units    = trim(units)
   call add_pvar_node(g%state,new)
   g%nstate = g%nstate + 1;
end subroutine add_state

function get_from_state(varname,g) result(v)
   character(len=*)    :: varname
   type(igroup)        :: g
   real(sp), pointer   :: v(:)
   type(pvar), pointer :: p
   
   p => get_pvar(g%state,varname)
   v => p%fvar
   
end function get_from_state

End Module mod_igroup
