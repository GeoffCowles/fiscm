!=======================================================================
! Fiscm Drivers 
!
! Description
!  - Main driver for calling ocean-model-specific advect/diffuse/find routines 
!
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!
!=======================================================================

Module mod_driver
use gparms
use mod_igroup
use ocean_model
use bio

implicit none
contains
!  s from z initialization
subroutine sz_ini(ng,g)
  integer     , intent(in   ) :: ng
  type(igroup), intent(inout), dimension(ng) :: g
  integer   :: n,np

   do n=1,ng
    np = g(n)%nind
    if(g(n)%space_dim == 3)then
    call sz_trans(np,g(n))
    endif
   end do
end subroutine  

!----------------------------------------------------
! advection-diffusion driver
!   call the appropriate advection and diffusion
!   routines, group-dependent
!----------------------------------------------------
subroutine do_adv_diff(ng,g,dT,time)
  integer     , intent(in   ) :: ng
  type(igroup), intent(inout), dimension(ng) :: g
  real(sp), intent(in) :: dT,time
  integer :: n,np
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)

  do n=1,ng

    !set group dimension
    np = g(n)%Tnind

    select case(g(n)%space_dim)
    case(0)
      cycle
    case(1)
      write(*,*)'error in driver'
      write(*,*)'1D simulation not setup'
      stop
    !--------------------------------------------------------------------
    case(2) ! 2D Problem
    !--------------------------------------------------------------------
      call advect2D( g(n) , dT,g(n)%nind)
      if(g(n)%hdiff_type == HDIFF_CONSTANT) call rw_hdiff_constant( g(n) ,dT) 
      if(g(n)%hdiff_type == HDIFF_VARIABLE) call rw_hdiff_variable( g(n),dT)
    !--------------------------------------------------------------------
    case(3) ! 3D Problem
    !--------------------------------------------------------------------
      call advect3D( g(n) ,dT ,g(n)%nind,time)
      if(g(n)%hdiff_type == HDIFF_CONSTANT) call rw_hdiff_constant( g(n) ,dT) 
      if(g(n)%hdiff_type == HDIFF_VARIABLE) call rw_hdiff_variable( g(n),dT)
      if(g(n)%vdiff_type == VDIFF_VARIABLE) call rw_vdiff(g(n), dT, g(n)%vdiff_substeps)  
      if(g(n)%vdiff_type == VDIFF_SPLINED ) call rw_vdiff_splined(g(n), dT, g(n)%vdiff_substeps)  
      if(g(n)%vdiff_type == VDIFF_BINNED  ) call rw_vdiff_binned(g(n), dT, g(n)%vdiff_substeps)  
    case default
      write(*,*)'space_dim must be [0,2,3]'
      stop
    end select

  end do

  !update elements
  call update_element(ng,g)

end subroutine do_adv_diff

!---------------------------------------------------------
! driver to advance the biology in time 
!---------------------------------------------------------
subroutine do_bio(ng,g,t,its)
  integer     , intent(in   ) :: ng
  type(igroup), intent(inout), dimension(ng) :: g
  real(sp), intent(in) :: t
  integer,  intent(in) :: its
  integer :: n
  do n=1,ng
    if(g(n)%biology .and. mod(its,g(n)%intvl_bio) == 0)then 
       call advance_bio(g(n),t)
    endif
  end do
end subroutine do_bio

!---------------------------------------------------------

!---------------------------------------------------------
! driver to interpolate forcing onto particle positions 
!---------------------------------------------------------
subroutine interp_forcing(ng,g,iframe)
  integer     , intent(in   ) :: ng,iframe
  type(igroup), intent(inout), dimension(ng) :: g
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  integer,  pointer :: cell(:)
  integer,  pointer :: istatus(:)
  real(sp), pointer :: s(:)
  real(sp), pointer :: f(:)
  real(sp), pointer :: i(:)

  integer :: n,vtype,v
  character(len=fstr) :: evar,svar

  
  do n=1,ng

    !skip 0-d groups 
    if(g(n)%space_dim < 2)cycle

    !get the horizontal particle positions for the group
    call get_state('x',g(n),x)
    call get_state('y',g(n),y)
    call get_state('cell',g(n),cell)
    call get_state('status',g(n),istatus)

    !get the vertical position if sim is 3-d
    if(g(n)%space_dim == 3)then
      call get_state('s',g(n),s)
    endif
    !loop over external variables, interpolate onto points
    do v=1,g(n)%next
  
      !get varname and type
      svar  = g(n)%ext_var(v,1)
      evar  = g(n)%ext_var(v,2)
      vtype = 1 !istype(g(n),svar) gwc debug

      !get pointer to var
      if(vtype == flt_type)then
        call get_state(svar,g(n),f) 
      elseif(vtype == int_type)then
        call get_state(evar,g(n),i) 
      else
        write(*,*)'error in interp_forcing'
        write(*,*)'cant deal with state variable of type: ',vtype
        stop
      endif
   
      !interpolate - 2D
      !interp in fvcom_drivers.f90
      if(g(n)%space_dim == 2)then
        if(vtype == flt_type)then
          call interp(g(n)%nind,x,y,cell,istatus,evar,f,iframe)
        elseif(vtype == int_type)then
        !  call interp(g(n)%nind,x,y,cell,istatus,evar,i)
        endif
      else !3D
        if(vtype == flt_type)then
          call interp(g(n)%nind,x,y,s,cell,istatus,evar,f,iframe) 
        elseif(vtype == int_type)then
        !  call interp(g(n)%nind,x,y,s,cell,istatus,evar,i)
        endif
      endif

    end do !loop over external vars

  end do !group loop
nullify(x)
nullify(y)
nullify(istatus)
nullify(cell)  
nullify(s)

end subroutine interp_forcing

!----------------------------------------------------
! determine in which element each particle in the  
! group resides.
!----------------------------------------------------
subroutine update_element(ng,g)
  integer, intent(in) :: ng
  type(igroup), intent(inout) :: g(ng)
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  integer , pointer :: cell(:)
  integer , pointer :: istatus(:)
  integer :: np,n

  do n=1,ng

    if(g(n)%space_dim < 2)cycle

    !set dimensinos
    np = g(n)%nind

    !set pointers to states 
    call get_state('x',g(n),x)
    call get_state('y',g(n),y)
    call get_state('cell',g(n),cell)
    call get_state('status',g(n),istatus)
    !update element containing cell 
    call find_element(np,x(1:np),y(1:np),cell(1:np),istatus(1:np))
  end do
nullify(x)
nullify(y)
nullify(cell)
nullify(istatus)

end subroutine update_element

  
!---------------------------------------------------------
! activate particles after spawning time reached 
!   for each particle:
!     if the status is still initial (status = 0) and
!     we have exceeded the spawning time (or < for backwards)
!     set status to active 
!---------------------------------------------------------
subroutine activate(ng,g,t,direction)
  integer     , intent(in   ) :: ng
  type(igroup), intent(inout), dimension(ng) :: g
  real(sp), intent(in) :: t
  integer,  intent(in) :: direction
  !---------------------------------
  real(sp), pointer :: tspawn(:)
  integer , pointer :: istatus(:)
  integer :: n,p,np

  !forward model
  if(direction > 0)then
    do n=1,ng
      call get_state('tspawn',g(n),tspawn)
      call get_state('status',g(n),istatus)
      np = g(n)%nind
      do p=1,np
        if(istatus(p) == UNKNOWN .and. t >= tspawn(p)) istatus(p) = ACTIVE
        !gwc if(t <= zptini(p)) istatus(p) = UNKNOWN
      end do
    end do
  !reverse model
  else
    do n=1,ng
      call get_state('tspawn',g(n),tspawn)
      call get_state('status',g(n),istatus)
      np = g(n)%nind
      do p=1,np
        if(istatus(p) == 0 .and. t <= tspawn(p)) istatus(p) = 1
        !gwc if(t >= zptini(p)) istatus(p) = UNKNOWN

      end do
    end do
  endif
       
end subroutine activate

  

End Module mod_driver
