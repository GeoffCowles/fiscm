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

!group - for i/o
integer  :: Tnind,nind,space_dim,nstate
integer  :: hdiff_type,vdiff_type
real(sp) :: DT_bio
real(sp) :: tlast_out,start_out,intvl_out
real(sp) :: group_deltat
real(sp) :: hdiff_const_val,vdiff_const_val
logical  :: biology
character(len=fstr) :: group_name,fname_out
character(len=fstr) :: statefile,paramfile

Namelist /NML_GROUP/     &
   & Tnind,              &
   & space_dim,          &
   & group_name,         &
   & hdiff_type,         &
   & hdiff_const_val,    &
   & vdiff_type,         &
   & vdiff_const_val,    &
   & DT_bio,             &
   & biology,            &
   & intvl_out,          &
   & start_out,          &
   & nstate,             &
   & statefile,          &
   & paramfile           

character(len=fstr) :: state_varname
character(len=fstr) :: state_longname
character(len=fstr) :: state_units
integer             :: state_netcdf_out
integer             :: state_vartype
integer             :: state_initval_int
real(sp)            :: state_initval_flt
character(len=fstr) :: state_from_ext_var

Namelist /NML_STATEVAR/   &
   & state_varname,             &
   & state_longname,            &
   & state_units,               &
   & state_netcdf_out,          &
   & state_vartype,             &
   & state_initval_int,         &
   & state_initval_flt,         &
   & state_from_ext_var

!Group type
type igroup 
  integer             :: id
  integer             :: Tnind
  integer             :: nind
  integer             :: space_dim
  character(len=fstr) :: name
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
  character(len=fstr) :: statefile
  character(len=fstr) :: paramfile
  type(pvar_list)     :: state
  integer             :: nstate
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
!function group_(i,id,bio) result (g)
function group_(fid,id) result(g)
  type(igroup) :: g
  integer, intent(in) :: fid
  integer, intent(in) :: id
  integer :: ios,i,ns
  logical :: fexist
  real(sp) :: fval
  integer  :: ival

  !read the group namelist from unit fid
  read(unit=fid,nml=nml_group,iostat=ios)
  if(ios /= 0)then
    write(*,*)'fatal error: could not read group namelist from fiscm.nml'
    stop
  endif

  !check problem size
  if(Tnind < 1)then
    write(*,*) 'error creating group: size must be greater than 0'
    stop 
  endif
  !check problem dimension
  if(space_dim < 0 .or. space_dim > 3)then
    write(*,*)'fatal error: space dim must be 0,2,3'
    stop
  endif

  !convert time step to days
  DT_bio = DT_bio*sec_2_day
  intvl_out = intvl_out*sec_2_day

  g%id     = id           
  g%Tnind  = Tnind   
  g%nind   = Tnind    
  g%space_dim = space_dim
  g%name   = group_name   
  g%hdiff_type = hdiff_type
  g%hdiff_const_val = hdiff_const_val
  g%vdiff_type = vdiff_type
  g%vdiff_const_val = vdiff_const_val
  g%DT_bio = DT_bio 
  g%biology = biology
  g%intvl_out = intvl_out
  g%tlast_out = 0.0
  g%start_out = start_out
  g%frame_out = 0
  g%fid_out   = 0
  g%fname_out = ""
  g%nstate = 0
  g%statefile = statefile
  g%paramfile = paramfile
  g%state = pvar_list_() 

  !(x,y,z) - particle position
  if(g%space_dim > 1)then
    fval=1.;call add_state(g,'x','x location of particle','m',NETCDF_YES,fval)
    fval=0.;call add_state(g,'y','y location of particle','m',NETCDF_YES,fval)
    fval=0.;call add_state(g,'pathlength','integrated trajectory','m',NETCDF_YES,fval)
  endif

  if(g%space_dim == 3)then
    fval=0.;call add_state(g,'z','z location of particle','m',NETCDF_YES,fval)
  endif

  !(status) - particle status
  call add_state(g,'status','particle status','-',NETCDF_YES,1)

  !(status) - spawn time
  fval=0.;call add_state(g,'tspawn','time of spawning','sec',NETCDF_NO,fval)

  !--------------------------------------------------------------------
  !user-defined state variables
  !--------------------------------------------------------------------
  if(statefile=="NONE" .and. nstate > 0)then
    write(*,*)'fatal error: group: ',trim(group_name) 
    write(*,*)'number of state variables from namelist: ',nstate
    write(*,*)'no state variable file specified'
    stop
  endif

  if(nstate > 1)then
    inquire(file=statefile,exist=fexist)
    if(.not.fexist)then
      write(*,*)'fatal error: statefile ',trim(statefile),' does not exist'
      stop
    endif
    open(unit=fid+1,file=trim(statefile),form='formatted')
    do ns = 1,nstate

      !read the group namelist from unit fid
      read(unit=fid+1,nml=nml_statevar,iostat=ios)
      if(ios /= 0)then
        write(*,*)'fatal error: could not read statevar namelist from: ',trim(statefile)
        stop
      endif
      if(state_vartype == itype)then
        call add_state(g,trim(state_varname),trim(state_longname),trim(state_units), &
                     state_netcdf_out,state_initval_int,state_from_ext_var)
      else if(state_vartype == ftype)then
          call add_state(g,trim(state_varname),trim(state_longname),trim(state_units), &
                     state_netcdf_out,state_initval_flt,state_from_ext_var)
      else
        write(*,*)'fatal error: not setup for variable type ',state_vartype,' in group_'
        stop
      endif
    end do
    close(34)
  end if !g%nstate > 1


end function group_

subroutine add_state_fvec(g,varname,longname,units,output,init_val,ext_name)
   type(igroup)        :: g 
   character(len=*)    :: varname
   character(len=*)    :: longname
   character(len=*)    :: units
   integer             :: output
   real(sp)            :: init_val
   character(len=*), optional :: ext_name

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
  if(present(ext_name))then
     new%from_ext_var = trim(ext_name)
   else
     new%from_ext_var = "NONE"
   endif

   call add_pvar_node(g%state,new)
   g%nstate = g%nstate + 1;
end subroutine add_state_fvec

subroutine add_state_ivec(g,varname,longname,units,output,init_val,ext_name)
   type(igroup)        :: g 
   character(len=*)    :: varname
   character(len=*)    :: longname
   character(len=*)    :: units
   integer             :: output
   integer             :: init_val
   character(len=*), optional :: ext_name
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
   if(present(ext_name))then
     new%from_ext_var = trim(ext_name)
   else
     new%from_ext_var = "NONE"
   endif
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

   write(*,*)
   write(*,*)'==========================Group Summary=========================='
   write(*,'(A21,I10)')'group id number   :: ',g%id
   write(*,'(A21,A20)')'group name        :: ',g%name
   write(*,'(A21,I10)')'total individuals :: ',g%Tnind
   write(*,'(A21,I10)')'assigned slots    :: ',g%nind
   write(*,'(A21,I10)')'spatial dimension :: ',g%space_dim
  
   select case(g%hdiff_type)
     case(HDIFF_NONE) 
       write(*,'(A21,A20)')'hdiff type        :: ','HDIFF_NONE'
     case(HDIFF_CONSTANT) 
       write(*,'(A21,A20)')'hdiff type        :: ','HDIFF_CONSTANT'
       write(*,'(A21,F10.6)')'hdiff constant    :: ',g%hdiff_const_val
     case(HDIFF_VISSER)   
       write(*,'(A21,A20)')'hdiff type        :: ','HDIFF_VISSER'
   end select
   select case(g%vdiff_type)
     case(VDIFF_NONE) 
       write(*,'(A21,A20)')'vdiff type        :: ','VDIFF_NONE'
     case(VDIFF_CONSTANT) 
       write(*,'(A21,A20)')'vdiff type        :: ','VDIFF_CONSTANT'
       write(*,'(A21,A10)')'biology           :: ','ACTIVE'
     case(VDIFF_VISSER)   
       write(*,'(A21,A20)')'vdiff type        :: ','VDIFF_VISSER'
   end select
   if(g%biology)then
     write(*,'(A21,A10)')'biology           :: ','ACTIVE'
     write(*,'(A21,F10.2)')'bio time step(s)  :: ',g%DT_bio*day_2_sec
   else
     write(*,'(A21,A10)')'biology           :: ','INACTIVE'
   endif
   write(*,'(A21,F10.2)')'output interval(s):: ',g%intvl_out*day_2_sec
   write(*,'(A21,F10.2)')'output starts at  :: ',g%start_out
   write(*,'(A21,I10)')'num. state vars   :: ',g%nstate

   call print_state_vars(g%state)

end subroutine print_group_summary

!---------------------------------------------------------
! check status of particles, report, return false if
! no active particles left
!---------------------------------------------------------
function checkstatus(ng,g,time) result(validsim)
   implicit none
   logical :: validsim 
   integer, intent(in)     :: ng
   type(igroup),intent(in) :: g(ng)
   real(sp), intent(in)    :: time
   integer :: nACTIVE  
   integer :: nSETTLED 
   integer :: nDEAD    
   integer :: nUNKNOWN 
   integer :: nEXITED  
   integer :: nTOTAL  
   integer, pointer :: istatus(:)
   integer :: n,NI
   integer, allocatable :: mark(:)

   !assume we still have some active particles
   validsim = .true.

   !initialize counters
   nACTIVE  = 0
   nSETTLED = 0
   nDEAD    = 0
   nUNKNOWN = 0
   nEXITED  = 0
   nTOTAL   = 0

   do n=1,ng
     NI = g(n)%Nind
     call get_state('status',g(n),istatus)
     allocate(mark(NI)) ; mark = 0
     where(istatus(1:NI) == SETTLED)mark = 1
       nSETTLED = nSETTLED + sum(mark)
       mark = 0
     where(istatus(1:NI) == DEAD)mark = 1
       nDEAD = nDEAD + sum(mark)
       mark = 0
     where(istatus(1:NI) == ACTIVE)mark = 1
       nACTIVE = nACTIVE + sum(mark)
       mark = 0
     where(istatus(1:NI) == EXITED)mark = 1
       nEXITED = nEXITED + sum(mark)
       mark = 0
     where(istatus(1:NI) == UNKNOWN)mark = 1
       nUNKNOWN = nUNKNOWN + sum(mark)
       mark = 0
     deallocate(mark)
     nTOTAL = nTOTAL + NI 
   end do
!   write(*,*)time,nTOTAL,nSETTLED,nDEAD,nACTIVE,nEXITED,nUNKNOWN

   if(nACTIVE == 0) validsim = .false.

end function checkstatus


End Module mod_igroup
