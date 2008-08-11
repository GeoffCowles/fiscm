!=======================================================================
! Fiscm FVCOM Driver 
!
! Description
!    - Read and setup mesh
!    - Advect
!    - Diffuse
!    - Locate Particles
!    
! Comments: 
!    - Advection currently uses 1st Order Euler Step, update to a stiffer 
!      solver (for example 4-stage RK used by C. Chen)
!    - Vertical diffusion uses Vissers modified random walk
!    - Vertical diffusion improvement possibilities:
!         a.) splines (previously used by J. Pringle, R. Ji, and M. Huret) 
!         b.) tensioned splines to avoid needing to smooth diffusivity 
!             before splining
!         c.) binned random walk (Thygesen and Adlandsvik, MEPS 347,2007)
!    - In theory fiscm could be driven by an alternate ocean model if 
!      these routines (initialize,advect,diffuse,find) are recoded
!
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!=======================================================================
Module Ocean_Model 
use gparms
use mod_igroup
use forcing
implicit none

!mesh params
integer, parameter :: level_based = 1
integer, parameter :: layer_based = 2

!Runge-Kutta integration coefficients
integer, parameter :: nstage = 1
real(sp),parameter :: alpha(1) = 1.0 

!integer, parameter :: nstage = 4
!real(sp),parameter :: alpha(4) = (/0.0_sp,0.5_sp,0.5_sp,1.0_sp/)
!real(sp), parameter:: beta(4)  = (/1.0_sp/6.0_sp,1.0_sp/3.0_sp, 1.0_sp/3.0_sp,1.0_sp/6.0_sp/)
!real(sp), parameter:: delta(4) = (/0.0_sp,0.5_sp,0.5_sp,1.0_sp/)

!dimensions
integer :: N_lev
integer :: N_lay
integer :: N_verts
integer :: N_elems 
integer :: Max_Elems

!mesh 
logical :: mesh_setup = .false.
real(sp), pointer :: xm(:)
real(sp), pointer :: ym(:)
real(sp), pointer :: xc(:)
real(sp), pointer :: yc(:)
real(sp), pointer :: hm(:)
real(sp), pointer :: aw0(:,:)
real(sp), pointer :: awx(:,:)
real(sp), pointer :: awy(:,:)
integer, pointer  :: tri(:,:)
integer, pointer  :: ntve(:)
integer, pointer  :: nbve(:,:)
real(sp), pointer :: siglay(:,:)
real(sp), pointer :: siglev(:,:)
real(sp), pointer :: esiglay(:,:)
real(sp), pointer :: esiglev(:,:)
logical :: grid_metrics

interface interp
  module procedure interp_float2D
  module procedure interp_float3D
end interface

interface interp_from_nodes
  module procedure interp_flt_from_nodes 
end interface

contains

!----------------------------------------------------
! Read the mesh and interpolation coefficients 
! Add mesh to output files for viz
!----------------------------------------------------
subroutine ocean_model_init(ng,g,lsize,varlist)
  use utilities, only : drawline,cfcheck
  integer, intent(in)      :: ng
  type(igroup), intent(in) :: g(ng)
  integer, intent(inout)   :: lsize
  character(len=*)         :: varlist(max_state_vars)
  !----------------------------------------------------
  character(len=mstr) :: msg
  integer :: dimid,varid,fid
  character(len=fstr) :: dname
  integer :: subset(3)
  integer :: i,n,ierr,ofid,k
  integer :: x_vid,y_vid,h_vid,nv_vid,nele_did,node_did,three_did,zeta_vid

  !return if group spatial dimensions are all 0-d or 1-d
  if(maxval(g%space_dim) < 2)return

  !get the forcing file netcdf id
  fid = get_ncfid()

  !add required time dependent variables to the list
  do n=1,ng
    if(g(n)%space_dim == 2)then
      lsize = lsize + 1 ; varlist(lsize) = 'ua'
      lsize = lsize + 1 ; varlist(lsize) = 'va'
      lsize = lsize + 1 ; varlist(lsize) = 'zeta'
      lsize = lsize + 1 ; varlist(lsize) = 'h'
    elseif(g(n)%space_dim ==3)then
      lsize = lsize + 1 ; varlist(lsize) = 'u'
      lsize = lsize + 1 ; varlist(lsize) = 'v'
      lsize = lsize + 1 ; varlist(lsize) = 'zeta'
      lsize = lsize + 1 ; varlist(lsize) = 'omega'
      lsize = lsize + 1 ; varlist(lsize) = 'h'
      lsize = lsize + 1 ; varlist(lsize) = 'kh'
    endif
  end do

  !determine number of elements
  msg = "dimension 'nele' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, 'nele', dimid ),msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, N_elems ))

  !determine number of nodes 
  msg = "dimension 'node' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, 'node', dimid ),msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, N_verts ))

  !determine number of layers
  msg = "dimension 'siglay' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, 'siglay', dimid ),msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, N_lay ))
  N_lev = N_lay + 1

  !allocate dataspace
  allocate(xm(N_verts))
  allocate(ym(N_verts))
  allocate(xc(N_elems))
  allocate(yc(N_elems))
  allocate(hm(N_verts))
  allocate(aw0(N_elems,3)) ; aw0 = a3rd
  allocate(awx(N_elems,3)) ; awx = zero
  allocate(awy(N_elems,3)) ; awy = zero
  allocate(tri(N_elems,3))
  allocate(siglay(N_verts,N_lay))
  allocate(siglev(N_verts,N_lev))
  allocate(esiglay(N_elems,N_lay))
  allocate(esiglev(N_elems,N_lev))

  !read in mesh
  msg = "error reading x coordinate"
  call ncdchk( nf90_inq_varid(fid,'x',varid),msg )
  call ncdchk(nf90_get_var(fid, varid, xm),msg)
  msg = "error reading y coordinate"
  call ncdchk( nf90_inq_varid(fid,'y',varid),msg )
  call ncdchk(nf90_get_var(fid, varid, ym),msg)
  msg = "error reading h coordinate"
  call ncdchk( nf90_inq_varid(fid,'h',varid),msg )
  call ncdchk(nf90_get_var(fid, varid, hm),msg)
  msg = "error reading nv coordinate"
  call ncdchk( nf90_inq_varid(fid,'nv',varid),msg )
  call ncdchk(nf90_get_var(fid, varid, tri),msg)
  msg = "error reading aw0"
  !read aw0 if they exist, otherwise use 1st order interpolation
  if(ncdscan( nf90_inq_varid(fid,'aw0',varid),msg ) )then
    call ncdchk(nf90_get_var(fid, varid, aw0),msg)
    msg = "error reading awx"
    call ncdchk( nf90_inq_varid(fid,'awx',varid),msg )
    call ncdchk(nf90_get_var(fid, varid, awx),msg)
    msg = "error reading awy"
    call ncdchk( nf90_inq_varid(fid,'awy',varid),msg )
    call ncdchk(nf90_get_var(fid, varid, awy),msg)
  else 
    write(*,*)'WARNING:::::: AW0 AWX AWY Do NOT exist in forcing file'
    write(*,*)'Proceeding with: 1st Order interpolation' 
    write(*,*)'AW0 = 1/3; AWX=AWY = 0'
    write(*,*)'In the future, select [grid metrics] in your NetCDF namelist'
  endif
  msg = "error reading siglay"
  call ncdchk( nf90_inq_varid(fid,'siglay',varid),msg )
  call ncdchk(nf90_get_var(fid, varid, siglay),msg)
  msg = "error reading siglev"
  call ncdchk( nf90_inq_varid(fid,'siglev',varid),msg )
  call ncdchk(nf90_get_var(fid, varid, siglev),msg)

  !read secondary connectivity (nbve/ntve) 
  grid_metrics = .false.
  msg = "dimension 'maxelem' not in the netcdf dataset"
  if(ncdscan( nf90_inq_dimid(fid,'maxelem',dimid),msg ) )then
    call ncdchk(nf90_inquire_dimension(fid, dimid, dname, Max_Elems))
    allocate(ntve(N_verts))
    allocate(nbve(N_verts,Max_Elems))
    msg = "error reading ntve"
    call ncdchk( nf90_inq_varid(fid,'ntve',varid),msg )
    call ncdchk(nf90_get_var(fid, varid, ntve),msg)
    msg = "error reading nbve"
    call ncdchk( nf90_inq_varid(fid,'nbve',varid),msg )
    call ncdchk(nf90_get_var(fid, varid, nbve),msg)
    grid_metrics = .true.
  else 
    write(*,*)'WARNING:::::: NTVE/NBVE Do NOT exist in forcing file'
    write(*,*)'This will slow down the element search procedure'
    write(*,*)'Proceeding with: 1st Order interpolation' 
    write(*,*)'In the future, select [grid metrics] in your NetCDF namelist'
  endif

  !calculate cell center coordinates
  do i=1,N_elems
    subset = tri(i,1:3)
    xc(i)  = a3rd*(sum(xm(subset)))
    yc(i)  = a3rd*(sum(ym(subset)))
  end do

  !calculate cell-center siglay/siglev
  do i=1,N_elems
    subset      = tri(i,1:3)
    do k=1,N_lay
      esiglay(i,k)  = a3rd*(sum(siglay(subset,k)))
      esiglev(i,k)  = a3rd*(sum(siglev(subset,k)))
    end do
    esiglev(i,N_lev)  = a3rd*(sum(siglev(subset,N_lev)))
  end do

  !determine nbve/ntve - secondary connectivity, used
  !for searching element containing point

  !mark boundary elements


  !calculate node-based interpolation coefficients  

  call drawline('-')
  write(*,*)'FVCOM mesh stats '
  call drawline('-')
  write(*,*)'Number of elements:: ',N_elems
  write(*,*)'Number of nodes   :: ',N_verts
  write(*,*)'Number of sig levs:: ',N_lev  
  write(*,*)'xmin              :: ',minval(xm)
  write(*,*)'xmax              :: ',maxval(xm)
  write(*,*)'ymin              :: ',minval(ym)
  write(*,*)'ymax              :: ',maxval(ym)

  !flag that mesh is setup
  mesh_setup = .true.

  !------------------------------------------------------
  !dump mesh to mesh.nc for viz
  !------------------------------------------------------
  call cfcheck( nf90_create("mesh.nc",nf90_clobber,ofid) ) 

  !dimensions
  call cfcheck(nf90_def_dim(ofid,"nele",N_elems, nele_did) )
  call cfcheck(nf90_def_dim(ofid,"node",N_verts, node_did) )
  call cfcheck(nf90_def_dim(ofid,"three",3, three_did) )

  !x
  call cfcheck( nf90_def_var(ofid,"x",nf90_float,(/node_did/), x_vid) )
  call cfcheck( nf90_put_att(ofid, x_vid,"long_name","nodal x-coordinate") )
  call cfcheck( nf90_put_att(ofid, x_vid,"units","meters") )
  call cfcheck( nf90_put_att(ofid, x_vid,"grid","TWOD_MESH") )

  !y
  call cfcheck( nf90_def_var(ofid,"y",nf90_float,(/node_did/), y_vid) )
  call cfcheck( nf90_put_att(ofid, y_vid,"long_name","nodal y-coordinate") )
  call cfcheck( nf90_put_att(ofid, y_vid,"units","meters") )
  call cfcheck( nf90_put_att(ofid, y_vid,"grid","TWOD_MESH") )

  !h
  call cfcheck( nf90_def_var(ofid,"h",nf90_float,(/node_did/), h_vid) )
  call cfcheck( nf90_put_att(ofid, h_vid,"long_name","Bathymetry") )
  call cfcheck( nf90_put_att(ofid, h_vid,"units","meters") )
  call cfcheck( nf90_put_att(ofid, h_vid,"positive","down") )
  call cfcheck( nf90_put_att(ofid, h_vid,"standard_name","depth") )
  call cfcheck( nf90_put_att(ofid, h_vid,"grid","fvcom_grid") )

  !el
  call cfcheck( nf90_def_var(ofid,"zeta",nf90_float,(/node_did/), zeta_vid) )
  call cfcheck( nf90_put_att(ofid, zeta_vid,"long_name","Water Surface Elevation") )
  call cfcheck( nf90_put_att(ofid, zeta_vid,"units","meters") )
  call cfcheck( nf90_put_att(ofid, zeta_vid,"positive","up") )
  call cfcheck( nf90_put_att(ofid, zeta_vid,"standard_name","sea_surface_elevation") )
  call cfcheck( nf90_put_att(ofid, zeta_vid,"type","data") )

  !nv
  call cfcheck( nf90_def_var(ofid,"nv",nf90_int,(/nele_did,three_did/), nv_vid) )
  call cfcheck( nf90_put_att(ofid, nv_vid,"long_name","nodes surrounding element") )

  !globals
  call cfcheck( nf90_put_att(ofid,nf90_global,"source"     ,"FVCOM") )
  call cfcheck( nf90_put_att(ofid,nf90_global,"Conventions"     ,"CF-1.0") )

  call cfcheck( nf90_enddef(ofid) )

  call cfcheck( nf90_put_var(ofid,x_vid,xm)) 
  call cfcheck( nf90_put_var(ofid,y_vid,ym)) 
  call cfcheck( nf90_put_var(ofid,h_vid,hm)) 
  call cfcheck( nf90_put_var(ofid,zeta_vid,hm*0.0)) 
  call cfcheck( nf90_put_var(ofid,nv_vid,tri,START=(/1,1/))) 

  !close the file
  call cfcheck( nf90_close(ofid) )


end subroutine ocean_model_init

!----------------------------------------------------
! Random-Walk horizontal diffusion with constant
!    turbulent eddy diffusivity
!
! m. huret use a Gauss dist. , this probably wont
! converge to the correct continuous form of diffusion
! here we will use a uniform random walk
!----------------------------------------------------
subroutine rw_hdiff_constant(g, dT)
  use utilities, only : normal,unitrand
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: dT
  
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  integer,  pointer :: istatus(:)
  integer  :: i,np
  real(sp) :: tscale
  
  !set pointers to x,y particle positions
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('status',g,istatus)

  !set dimensions for loops and time step
  np = g%nind

  !set diffusive time scale
  tscale = sqrt(2.*dT*g%hdiff_const_val)

  !horizontal random walk
  do i=1,np
    if(istatus(i) < 1)cycle
!    x(i) = x(i) + normal()*tscale
!    y(i) = y(i) + normal()*tscale
    x(i) = x(i) + unitrand()*tscale
    y(i) = y(i) + unitrand()*tscale
  end do

  !nullify pointers
  nullify(x)
  nullify(y)
  nullify(istatus)

end subroutine rw_hdiff_constant


!----------------------------------------------------
! Random-Walk horizontal diffusion with spatially
!   variable turbulent eddy diffusivity
!
! Use Visser's modified random walk to compute step
!----------------------------------------------------
subroutine rw_hdiff_variable(g, dT)
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: dT

  write(*,*)'error in rw_hdiff_variable'
  write(*,*)'hdiff visser is not yet setup'
  stop

end subroutine rw_hdiff_variable

!----------------------------------------------------
! Random-Walk vertical diffusion 
!
!   - use eddy diffusivity from the model (kh)
!   - use Vissers modified random walk to compute jump
!----------------------------------------------------
subroutine rw_vdiff(g, dT, nstep)
  use utilities, only : unitrand
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: dT
  integer,  intent(in) :: nstep
  !----------------------------
  integer,  pointer :: istatus(:)
  integer,  pointer :: cell(:)
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  real(sp), pointer :: s(:)
  real(sp), pointer :: h(:)
  real(sp), allocatable :: kh(:)
  real(sp), allocatable :: kh2(:)
  real(sp), allocatable :: ds(:)
  real(sp), allocatable :: dkh_ds(:)
  real(sp), allocatable :: zeta(:)
  real(sp), allocatable :: s_shift(:)
  real(sp), parameter :: delta_s = 0.05
  real(sp) :: deltaT,fac,randy,dz,depth,dkh_dz
  integer  :: n,p,np

  !set problem size and time step
  np = g%nind
  deltaT = dT/float(nstep)

  !set pointers to particle positions and status
  call get_state('status',g,istatus)
  call get_state('cell',g,cell)
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('s',g,s)
  call get_state('h',g,h)

  !allocate local data
  allocate(s_shift(np)); s_shift = zero 
  allocate(zeta(np)); zeta = zero 
  allocate(kh(np))  ; kh   = zero 
  allocate(kh2(np)) ; kh2  = zero 
  allocate(ds(np))  ; ds   = zero 
  allocate(dkh_ds(np)) ; dkh_ds  = zero 

  !set constants
  fac = (2./rvar)*deltaT  ![ 2*r^-1*deltaT], r = variance of uniform rw, set in gparms

  call interp(np,x,y,cell,istatus,'h',h)
  call interp(np,x,y,cell,istatus,'zeta',zeta)

  ! ==> loop over substeps
  do n=1,nstep
   
    !--------------------------------------------------
    ! calculate d(kh)/d(s) - brute force
    !--------------------------------------------------

    !set derivative step (backward near free surface)
    ds = delta_s
    where(s+delta_s > 0)ds = -delta_s

    !evaluate kh at both locations
    call interp(np,x,y,s,cell,istatus,'kh',kh)
    call interp(np,x,y,s+ds,cell,istatus,'kh',kh2)

    !form the derivative d(kh)/d(s)
    dkh_ds = (kh2-kh)/ds

    !function evaluation at [z + 0.5*dkh/dz*deltat] - Visser
    s_shift = s + dkh_ds*deltaT/((h+zeta)**2)
    call interp(np,x,y,s_shift,cell,istatus,'kh',kh)

    ! => main loop over particles
    do p=1,np
      if(istatus(p) < 1)cycle

      !update particle position using Visser modified random walk 
      depth  = h(p)+zeta(p)
      dkh_dz = dkh_ds(p)/depth
      dz     = dkh_dz*deltaT + unitrand()*sqrt(fac*kh(p)) !Visser-modified
      s(p)   = s(p) + dz/depth
!      dz     = unitrand()*sqrt(2.*kh(p))                 !naive

      !set boundary conditions at free surface and bottom
      s(p) = max(s(p),-(2.0+s(p))) !reflect off bottom
      s(p) = min(s(p),0.0)         !don't pierce free surface

    end do
    ! <= end particle loop

  end do
  ! <== end loop over substeps

  !deallocate workspace and nullify pointers
  deallocate(ds)
  deallocate(kh)
  deallocate(kh2)
  deallocate(dkh_ds)
  deallocate(zeta)
  deallocate(s_shift)
  nullify(x)
  nullify(y)
  nullify(s)
  nullify(h)
  nullify(istatus)

end subroutine rw_vdiff

!----------------------------------------------------
! Random-Walk vertical diffusion using splines
!   - spline vertical diffusivity to smooth profile
!   - adjust at ends
!   - use vissers modified random walk to compute jump
!
! Follow M. Huret's Implementation for the spline
!----------------------------------------------------
subroutine rw_vdiff_splined(g, dT, nstep)
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: dT
  integer,  intent(in) :: nstep
  real(sp) :: kh(N_lev)
  integer, pointer :: istatus(:)
  integer :: n,p,np

  call get_state('status',g,istatus)
  np = g%nind

  ! => main loop over particles
  do p=1,np
    if(istatus(p) < 1)cycle

    !interpolate eddy diffusivity at model control points in column

    !spline eddy diffusivity - assume constant during step

    !loop over substeps
    do n=1,nstep
      
      !evaluate eddy diffusivity value and gradient at x,y,s

      !calculate random step using Visser formulation 
       
      !use a random mixed layer, 1m from top and bottom 

    end do 

  end do 
  ! <= end main particle loop
  

end subroutine rw_vdiff_splined 

!----------------------------------------------------
! Random-Walk vertical diffusion using bins
!
!  - Thygesen and Adlandsvik, MEPS, v347, 2007
!----------------------------------------------------
subroutine rw_vdiff_binned(g, dT, nstep)
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: dT
  integer,  intent(in) :: nstep
end subroutine rw_vdiff_binned

!----------------------------------------------------
! 2-D Advection (use vertically-averaged velocity)
!----------------------------------------------------
subroutine advect2D(g,deltaT)
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: deltaT
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  integer , pointer :: cell(:)
  integer , pointer :: istatus(:)
  real(sp), allocatable :: ua(:)
  real(sp), allocatable :: va(:)
  integer  :: k,i,np

  !set problem size 
  np = g%nind

  !set pointers to states 
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('cell',g,cell)
  call get_state('status',g,istatus)

  !allocate local data space
  allocate(ua(np)) ; ua = 0.0
  allocate(va(np)) ; va = 0.0

  !interpolate velocity field to particle position
  call interp(np,x,y,cell,istatus,'ua',ua)
  call interp(np,x,y,cell,istatus,'va',va)

  !advance in time with Euler step
  do i=1,np 
    if(istatus(i) < 1)cycle
    x(i) = x(i) + alpha(1)*deltaT*ua(i)
    y(i) = y(i) + alpha(1)*deltaT*va(i)
  end do

  !update element containing particle => cell 
  call find_element(np,x,y,cell,istatus)

  !apply boundary conditions (reflect off solid)

  !update interpolated velocity 
  call interp(np,x,y,cell,istatus,'ua',ua)
  call interp(np,x,y,cell,istatus,'va',va)
  
  !nullify pointers
  nullify(x)
  nullify(y)
  deallocate(ua)
  deallocate(va)
  nullify(cell)
  nullify(istatus)

end subroutine advect2D



!----------------------------------------------------
! 3-D Advection 
!----------------------------------------------------
subroutine advect3D(g,deltaT)
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: deltaT
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  real(sp), pointer :: s(:)
  real(sp), pointer :: z(:)
  real(sp), pointer :: h(:)
  integer , pointer :: cell(:)
  integer , pointer :: istatus(:)
  real(sp), allocatable :: u(:)
  real(sp), allocatable :: v(:)
  real(sp), allocatable :: omega(:)
  real(sp), allocatable :: zeta(:)
  integer  :: k,i,np

  !set dimensions for loops and time step
  np = g%nind

  !set pointers to states 
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('s',g,s)
  call get_state('z',g,z)
  call get_state('h',g,h)
  call get_state('cell',g,cell)
  call get_state('status',g,istatus)

  !allocate local data space
  allocate(zeta(np) ) ; zeta  = 0.0
  allocate(u(np)    ) ; u     = 0.0
  allocate(v(np)    ) ; v     = 0.0
  allocate(omega(np)) ; omega = 0.0

  !interpolate velocity field,depth,free surface to particle position
  call interp(np,x,y,s,cell,istatus,'u',u)
  call interp(np,x,y,s,cell,istatus,'v',v)
  call interp(np,x,y,s,cell,istatus,'omega',omega)
  call interp(np,x,y,cell,istatus,'h',h)
  call interp(np,x,y,cell,istatus,'zeta',zeta)

  !advance
  do i=1,np 
    if(istatus(i) < 1)cycle
    x(i) = x(i) + alpha(1)*deltaT*u(i)
    y(i) = y(i) + alpha(1)*deltaT*v(i)
    s(i) = s(i) + alpha(1)*deltaT*omega(i)/(h(i)+zeta(i)) 
  end do

  !adjust s-coordinate value at surface and bottom boundaries
  do i=1,np
    s(i) = max(s(i),-(2.0+s(i))) !mirror bottom
    s(i) = min(s(i),0.0)         !stay below fs
  end do

  !calculate z value of particle (for visualization)
  do i=1,np
    z(i) = s(i)*(h(i)+zeta(i)) + zeta(i)  
  end do

  !update element containing particle => cell
  call find_element(np,x,y,cell,istatus)

  !apply boundary conditions => reflect off solid

  !disassociate pointers
  nullify(x)
  nullify(y)
  nullify(s)
  nullify(z)
  nullify(h)
  nullify(cell)
  nullify(istatus)

end subroutine advect3D

!----------------------------------------------------
! interpolation routine 
!   interpolate scalar field [vname] to vector [v] 
!   at particle position [x,y] in cell [cell]
!   if [istatus] < 0 (inactive) do not interpolate
!----------------------------------------------------
subroutine interp_float2D(np,x,y,cell,istatus,vname,v)
  integer, intent(in)    :: np
  real(sp),intent(in)    :: x(np)
  real(sp),intent(in)    :: y(np)
  integer, intent(in)    :: cell(np)
  integer, intent(in)    :: istatus(np)
  character(len=*)       :: vname
  real(sp),intent(inout) :: v(np)
  !-------------------------------
  real(sp), pointer :: field(:)  !pointer to FVCOM field 
  integer :: i,vts(3),icell,d1

  !determine dimension of forcing variable

  !point to forcing variable
  call get_forcing(trim(vname),field) 

  !determine the dimensions
  d1 = size(field,1)

  !interpolate to element-based quantities
  !gwc debug = 0th order
  if(d1 == N_elems)then  
    do i=1,np
      icell = cell(i)
      if(istatus(i) < 1 .or. icell == 0)cycle 
      v(i) = field(icell)
    end do
  elseif(d1 == N_verts)then !vertex-based 2-D array
    do i=1,np
      icell = cell(i)
      if(istatus(i) < 1 .or. icell == 0)cycle 
      vts  = tri(icell,1:3)
      v(i) = interp_from_nodes(icell,x(i),y(i),field(vts))
    end do
  else
    write(*,*)'field has horizontal dimensions that is not nodes or elements'
    write(*,*)'do not know how to interpolate'
    stop
  endif


end subroutine interp_float2D  

!----------------------------------------------------
! interpolation routine 
!   interp a 3D scalar to particle positions in 3-space
!----------------------------------------------------
subroutine interp_float3D(np,x,y,s,cell,istatus,vname,v)
  integer, intent(in)    :: np
  real(sp),intent(in)    :: x(np)
  real(sp),intent(in)    :: y(np)
  real(sp),intent(in)    :: s(np)
  integer, intent(in)    :: cell(np)
  integer, intent(in)    :: istatus(np)
  character(len=*)       :: vname
  real(sp),intent(inout) :: v(np)
  !-------------------------------

  real(sp), pointer :: field(:,:)  !pointer to FVCOM field 
  integer  :: i,d1,d2,icell,k1,k2,ktype
  real(sp) :: f1,f2,v1,v2
  integer  :: vts(3)


  !point to forcing variable
  call get_forcing(trim(vname),field) 

  !determine the dimensions
  d1 = size(field,1)
  d2 = size(field,2)

  !set vertical type (layers or levels)
  if(d2 == N_lev)then
    ktype = level_based 
  else if(d2 == N_lay)then
    ktype = layer_based 
  else
    write(*,*)'error in interp_float3D'
    write(*,*)'second dimension of field is not layers or levels'
    write(*,*)'==: ',d2
    stop
  endif

  !------------------------------------------------
  !interpolate to element-based quantities
  !------------------------------------------------
  if(d1 == N_elems)then  
    do i=1,np
      !set cell containing particle
      icell = cell(i)

      !loop if particle dead or out of domain
      if(istatus(i) < 1 .or. icell == 0)cycle 

      !determine the vertical layer brackets and interp coeffs
      call get_vert_interpcoefs(icell,s(i),k1,k2,f1,f2,ktype)
      
      !interpolate
      v(i) = f1*field(icell,k1) + f2*field(icell,k2)

    end do
  !------------------------------------------------
  !interpolate to node-based quantities
  !------------------------------------------------
  elseif(d1 == N_verts)then
    do i=1,np
      !set cell containing particle
      icell = cell(i)

      !loop if particle dead or out of domain
      if(istatus(i) < 1 .or. icell == 0)cycle 
        
 
      !determine the vertical layer brackets and interp coeffs
      call get_vert_interpcoefs(icell,s(i),k1,k2,f1,f2,ktype)

      !set vertices surrounding triangle
      vts  = tri(icell,1:3)
      v1 = interp_from_nodes(icell,x(i),y(i),field(vts,k1))
      v2 = interp_from_nodes(icell,x(i),y(i),field(vts,k2))
      v(i) = f1*v1 + f2*v2

      if(trim(vname)=='khh')then
        write(*,*)'interpping kh'
        write(*,*)i,s(i),k1,k2,v1,v2,f1,f2
      endif
    end do
  else
    write(*,*)'error in interp_float3D'
    write(*,*)'field has horizontal dimensions that is not nodes or elements'
    write(*,*)'do not know how to interpolate'
    stop
  endif


end subroutine interp_float3D

!----------------------------------------------------
! element search routine 
!----------------------------------------------------
subroutine find_element(np,x,y,cell,istatus,option)
  integer, intent(in)    :: np
  real(sp),intent(in)    :: x(np)
  real(sp),intent(in)    :: y(np)
  integer, intent(inout) :: cell(np)
  integer, intent(inout) :: istatus(np)
  character(len=*), optional :: option
  integer :: p

  do p=1,np
    if(istatus(p) < 0)cycle

    !try a quick find
    cell(p) = find_element_lazy(x(p),y(p),cell(p))
    if(cell(p) /= 0)cycle

    !failed, try a robust find
    cell(p) = find_element_robust(x(p),y(p)) 

    !failed, update status to lost
    if(cell(p) /= 0)cycle
    istatus(p) = EXITED
  end do

end subroutine find_element

!----------------------------------------------------
! find the element in which a point resides: lazy 
!    - look in last known element 
!    - look in last known elements neighbors
!    - should be reasonably robust
!----------------------------------------------------
function find_element_lazy(xp,yp,last) result(elem)
  use utilities
  implicit none
  real(sp), intent(in) :: xp
  real(sp), intent(in) :: yp 
  integer,  intent(in) :: last
  integer  :: elem,j,k,test,iney
  !---------------------------
  integer  :: pts(3)
  real(sp) :: xv(3),yv(3)

  !make sure we have a mesh  
  if(.not.mesh_setup)then
    write(*,*)'error in find_element_robust'
    write(*,*)'mesh is not setup yet'
    stop
  endif
  
  !initialize
  elem = 0 ; if(last == 0)return 

  !load vertices 
  pts(1:3) = tri(last,1:3)
  xv = xm(pts)
  yv = ym(pts)

  if(isintriangle(last,xp,yp,xv,yv))then     
    elem = last
  else
    if(grid_metrics)then
    outer: do j=1,3
      test = tri(last,j)
           do k=1,ntve(test)
              iney = nbve(test,k)
              pts(1:3) = tri(iney,1:3)
              xv = xm(pts)
              yv = ym(pts)
              if(isintriangle(iney,xp,yp,xv,yv))then
                 elem  = iney
                 exit outer
              end if
           end do
        end do outer
     end if !grid_metrics on
  endif

  return
end function find_element_lazy


!----------------------------------------------------
! find the element in which a point resides: robust
!    - search outward in increasing radius
!    - search only the closest max_check points
!----------------------------------------------------
function find_element_robust(xp,yp) result(elem)
  use utilities
  implicit none
  real(sp), intent(in) :: xp
  real(sp), intent(in) :: yp
  integer  :: elem
  !---------------------------
  real(sp) :: radlist(N_elems,1)
  real(sp) :: radlast
  real(sp) :: xv(3),yv(3)
  integer  :: pts(3),min_loc,locij(2),cnt
  integer, parameter :: max_check = 15
  
  !make sure we have a mesh  
  if(.not.mesh_setup)then
    write(*,*)'error in find_element_robust'
    write(*,*)'mesh is not setup yet'
    stop
  endif
  
  !initialize
  elem = 0
  cnt = 0
  radlist(1:N_elems,1) = sqrt((xc(1:N_elems)-xp)**2 + (yc(1:N_elems)-yp)**2)
  radlast = -1.0_sp
  l1:  do while(cnt < max_check)
    cnt = cnt+1
    locij   = minloc(radlist,radlist>radlast)
    min_loc = locij(1)
    if(min_loc == 0) then
      exit l1
    end if

    pts(1:3) = tri(min_loc,1:3) ; xv = xm(pts) ; yv = ym(pts)
    if(isintriangle(min_loc,xp,yp,xv,yv))then
      elem = min_loc
      exit l1
    end if
    radlast = radlist(min_loc,1)
  end do l1

  return
 end function find_element_robust

 function interp_flt_from_nodes(cell,x,y,fin) result(fout)
   integer,  intent(in) :: cell
   real(sp), intent(in) :: x
   real(sp), intent(in) :: y
   real(sp), intent(in) :: fin(3) 
   real(sp) :: fout
   real(sp) :: f0,fx,fy,deltaX,deltaY

   f0 = aw0(cell,1)*fin(1) + aw0(cell,2)*fin(2) + aw0(cell,3)*fin(3) 
   fx = awx(cell,1)*fin(1) + awx(cell,2)*fin(2) + awx(cell,3)*fin(3) 
   fy = awy(cell,1)*fin(1) + awy(cell,2)*fin(2) + awy(cell,3)*fin(3) 
   deltaX = x-xc(cell)
   deltaY = y-yc(cell)
   fout = f0 + fx*deltaX + fy*deltaY

 end function interp_flt_from_nodes

 subroutine get_vert_interpcoefs(cell,sloc,k1,k2,f1,f2,ktype) 
   integer,  intent(in)  :: cell
   real(sp), intent(in)  :: sloc
   integer,  intent(out) :: k1
   integer,  intent(out) :: k2
   real(sp), intent(out) :: f1
   real(sp), intent(out) :: f2
   integer,  intent(in ) :: ktype
   !---------------------------
   real(sp) :: ds,my_sloc
   integer  :: i,k

   my_sloc = max(sloc,-one)
   my_sloc = min(sloc,-tinynum)

   !level-based data 
   if(ktype == level_based)then
     do k=1,N_lev 
       if(my_sloc >= esiglev(cell,k))exit
     end do
     k1 = k - 1
     k2 = k 
     ds = esiglev(cell,k1)-esiglev(cell,k2) 
     f1 = (my_sloc-esiglev(cell,k2))/ds 
     f2 = 1.-f1 
   !layer-based data
   elseif(ktype == layer_based)then
     do k=1,N_lay 
       if(my_sloc > esiglay(cell,k))exit
     end do
     if(k==1)then
       k1 = 1
       k2 = 1
       f1 = one 
       f2 = zero
     elseif(k >= N_lay)then
       k1 = N_lay
       k2 = N_lay
       f1 = one 
       f2 = zero
     else
       k1 = k - 1
       k2 = k 
       ds = esiglay(cell,k1)-esiglay(cell,k2) 
       f1 = (my_sloc-esiglev(cell,k2))/ds 
       f2 = 1.-f1 
     endif
   !error            
   else
     write(*,*)'error in get_vert_interpcoefs'
     write(*,*)'ktype must be :',level_based,' or ',layer_based
     write(*,*)'but is: ',ktype
     stop
   endif
 end subroutine get_vert_interpcoefs

End Module Ocean_Model 
