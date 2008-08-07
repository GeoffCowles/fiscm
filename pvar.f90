!=======================================================================
! Fiscm Polymorphic Type 
!
! Description
!   Defines a Polymorphic Type to store user and fiscm-defined state  
!     variables
!   Enable access to data stored in linked-lists through pointers
!    
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!
!=======================================================================

Module mod_pvar


Use gparms
Use netcdf
Implicit None

Integer, parameter :: itype = 1
Integer, parameter :: ftype = 2
Integer, parameter :: ltype = 3
Integer, parameter :: ctype = 4


!-------------------------------------------------------
!Polymorphic Type
!-------------------------------------------------------
type pvar
  integer :: idim
  integer :: istype
  logical :: isintern
  character(len=fstr) :: varname
  character(len=fstr) :: longname
  character(len=fstr) :: units
  character(len=fstr) :: from_ext_var
  integer             :: output
  integer             :: nc_vid
!  switch to allocatable type components when gfortran 4.2+ 
!  in general use
!  integer,  allocatable, dimension(:) :: ivar
!  real(sp), allocatable, dimension(:) :: fvar
  integer,  pointer, dimension(:) :: ivar
  real(sp), pointer, dimension(:) :: fvar
end type pvar
  
!-------------------------------------------------------
!Polymorphic Type - Node (for linked-list)
!-------------------------------------------------------
type pvar_node
  type(pvar) :: v
  type(pvar_node), pointer :: next
end type pvar_node

!-------------------------------------------------------
!Linked List of polymorphic types
!-------------------------------------------------------
type pvar_list
  type(pvar_node), pointer :: begin
end type pvar_list
  
contains

!-------------------------------------------------------
! initialize a pvar_list
!-------------------------------------------------------
function pvar_list_() result(plist)
  type(pvar_list) :: plist
  integer :: astat
  allocate(plist%begin,stat=astat)
  if(astat /= 0)then
    write(*,*)'error creating linked list of state variables'
    stop
  endif
  nullify(plist%begin%next)
end function pvar_list_

!-------------------------------------------------------
! add a new node to the list
!-------------------------------------------------------
subroutine add_pvar_node(plist,p)
  type(pvar_list) :: plist
  type(pvar)      :: p
  type(pvar_node), pointer :: plst,pnew

  plst => plist%begin
  pnew => plst%next
  do 
   if(.not.associated(pnew))then
     allocate(plst%next)
     if(p%istype == ftype)then
       allocate(plst%next%v%fvar(p%idim))
     else if(p%istype == itype)then
       allocate(plst%next%v%ivar(p%idim))
     endif
     plst%next%v = p
     nullify(plst%next%next)
     exit
   endif

   plst => pnew
   pnew  => pnew%next
 end do

end subroutine add_pvar_node
  
 
!-------------------------------------------------------
!Search linked-list for a particular variable by name
! and pass pointer to data
!-------------------------------------------------------
function get_pvar(plist,vname) result(p)
  type(pvar), pointer :: p
  type(pvar_list) :: plist
  character(len=*), intent(in) :: vname
  logical :: found
  type(pvar_node), pointer :: plst,pnew
  found = .false.
  
  plst => plist%begin
  pnew => plst%next
  if(.not.associated(plst%next))then
    write(*,*)'plist has no nodes'
    stop
  endif
  
  
  do 
    if(.not.associated(plst%next)) exit

    if(pnew%v%varname == vname)then
      p => plst%next%v
      found = .true.
      exit
    else
      plst => pnew
      pnew => pnew%next
    endif

  end do

  if(.not.found)then
    write(*,*)'error in get_pvar'
    write(*,*)'variable: ',trim(vname),' is not a state variable for group'
    stop
  endif
  

end function get_pvar
    
subroutine print_state_vars(plist)
  use utilities
  type(pvar_list) :: plist
  type(pvar_node), pointer :: plst,pnew
  
  plst => plist%begin
  pnew => plst%next

  if(.not.associated(plst%next))then
    write(*,*)'plist has no nodes'
    stop
  endif
  call drawline("-")
  write(*,*)'state variable|          long name           | units    | external var'
  call drawline("-")

  do 
    if(.not.associated(plst%next)) exit

    write(*,'(A15,A1,A30,A1,A10,A1,A30)')pnew%v%varname,'|',pnew%v%longname,'|',pnew%v%units,'|',pnew%v%from_ext_var
    plst => pnew
    pnew => pnew%next
  end do

end subroutine print_state_vars

subroutine write_cdf_header_vars(plist,fid,dims)
  type(pvar_list) :: plist
  integer, intent(in) :: fid
  integer, intent(in) :: dims(2)
  type(pvar_node), pointer :: plst,pnew
  type(pvar), pointer :: dum
  integer :: ierr
  
  plst => plist%begin
  pnew => plst%next

  if(.not.associated(plst%next))then
    write(*,*)'plist has no nodes'
    stop
  endif
  do 
    if(.not.associated(plst%next)) exit
      !dump the header if user selected to dump
         if(pnew%v%output == NETCDF_YES)then
        dum => pnew%v
        select case(dum%istype)
        case(ftype)
          ierr = nf90_def_var(fid,dum%varname,nf90_float,dims,dum%nc_vid)
          if(ierr /= nf90_noerr) then
            write(*,*)'error defining netcdf header for variable: ',dum%varname
            write(*,*)trim(nf90_strerror(ierr)) ; stop 
          end if
        case(itype)
          ierr = nf90_def_var(fid,dum%varname,nf90_int  ,dims,dum%nc_vid)
          if(ierr /= nf90_noerr) then
            write(*,*)'error defining netcdf header for variable: ',dum%varname
            write(*,*)trim(nf90_strerror(ierr)) ; stop 
          end if
        case default
          write(*,*)'not setup for netcdf dumping state variable of type:',dum%istype
          stop
        end select
        ierr = nf90_put_att(fid,dum%nc_vid,"long_name",dum%longname)
        if(ierr /= nf90_noerr) then
          write(*,*)'error defining netcdf header for variable: ',dum%varname
          write(*,*)trim(nf90_strerror(ierr)) ; stop 
        end if

        ierr = nf90_put_att(fid,dum%nc_vid,"units",dum%units)
        if(ierr /= nf90_noerr) then
           write(*,*)'error defining netcdf header for variable: ',dum%varname
            write(*,*)trim(nf90_strerror(ierr)) ; stop 
        end if
      endif

    plst => pnew
    pnew => pnew%next
  end do

end subroutine write_cdf_header_vars

subroutine write_cdf_data(plist,fid,frame)
  type(pvar_list) :: plist
  integer, intent(in) :: fid
  integer, intent(in) :: frame
  type(pvar_node), pointer :: plst,pnew
  type(pvar), pointer :: dum
  integer :: ierr
  integer :: dims(2)

  dims(1) = 1
  dims(2) = frame
  
  plst => plist%begin
  pnew => plst%next

  if(.not.associated(plst%next))then
    write(*,*)'plist has no nodes'
    stop
  endif
  do 
    if(.not.associated(plst%next)) exit
      !dump the variable if user selected to dump
      if(pnew%v%output == NETCDF_YES)then
        dum => pnew%v
        select case(dum%istype)
        case(ftype)
          ierr = nf90_put_var(fid, dum%nc_vid, dum%fvar,START=dims)
          if(ierr /= nf90_noerr) then
            write(*,*)'error dumping frame for variable: ',dum%varname
            write(*,*)trim(nf90_strerror(ierr)) ; stop 
          end if
        case(itype)
          ierr = nf90_put_var(fid, dum%nc_vid, dum%ivar,START=dims)
          if(ierr /= nf90_noerr) then
             write(*,*)'error dumping frame for variable: ',dum%varname
             write(*,*)trim(nf90_strerror(ierr)) ; stop 
          end if
        case default
          write(*,*)'not setup for netcdf dumping state variable of type:',dum%istype
          stop
        end select

      endif
    plst => pnew
    pnew => pnew%next
  end do
end subroutine write_cdf_data
  
End Module mod_pvar









