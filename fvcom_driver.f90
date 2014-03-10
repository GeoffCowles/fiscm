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
integer, parameter :: mstage = 4

real(sp),parameter :: A_RK(4) = (/0.0_sp,0.5_sp,0.5_sp,1.0_sp/)

real(sp), parameter:: B_RK(4)  = (/1.0_sp/6.0_sp,1.0_sp/3.0_sp, 1.0_sp/3.0_sp,1.0_sp/6.0_sp/)

real(sp), parameter:: C_RK(4) = (/0.0_sp,0.5_sp,0.5_sp,1.0_sp/)
!Added by Xinyou Lin
  integer, allocatable :: nbe(:,:)            !!INDICES OF ELMNT NEIGHBORS^M      
  integer, allocatable :: isonb(:)            !!NODE MARKER = 0,1,2   ^M
  integer, allocatable :: isbce(:)     
  integer, allocatable :: nbvt(:,:)
!dimensions
integer :: N_lev
integer :: N_lay
integer :: N_verts
integer :: N_elems 
integer :: Max_Elems

!mesh 
integer :: iunit,ios
logical :: mesh_setup = .false.
real(sp), pointer :: xm(:),xm0(:)
real(sp), pointer :: ym(:),ym0(:)
real(sp), pointer :: xc(:),xc0(:)
real(sp), pointer :: yc(:),yc0(:)
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

!added by Xinyou
real(sp), pointer :: a1u(:,:)
real(sp), pointer :: a2u(:,:)
character(len=10)  :: x_char,y_char,h_char,u_char,v_char, &
  kh_char,viscofm_char,ua_char,va_char,nv_char,nbe_char,aw0_char,awx_char,awy_char,a1u_char, &
  a2u_char,nele_char,node_char,zeta_char,omega_char,         &
  siglay_char,siglev_char,wu_char,wv_char
  Namelist /NML_NCVAR/  &
           x_char,   &
           y_char,   &
           h_char,   &
          nv_char,   &
        nele_char,   &
        node_char,   &
         nbe_char,   &
         aw0_char,   &
         awx_char,   &
         awy_char,   &
         a1u_char,   &
         a2u_char,   &
      siglay_char,   &
      siglev_char,   &
          ua_char,   &
          va_char,   &
        zeta_char,   &
           h_char,   &
          wu_char,   &
          wv_char,   &
           u_char,   &
           v_char,   &
       omega_char,   &
          kh_char,   &
     viscofm_char 
 
!

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
  !--------------------------------------------------------
  ! open and read  time dependent variables namelist:  nml_ncvar
  !--------------------------------------------------------
  open(unit=iunit,file=trim(runcontrol),form='formatted')
  read(unit=iunit,nml=nml_ncvar,iostat=ios)
  if(ios /= 0)then
    write(*,*)'fvcom:fatal error: could not read fiscm namelist from',trim(runcontrol)
    stop
  endif

  do n=1,ng
    if(g(n)%space_dim == 2)then
      lsize = lsize + 1 ; varlist(lsize) = ua_char
      lsize = lsize + 1 ; varlist(lsize) = va_char
      lsize = lsize + 1 ; varlist(lsize) = zeta_char
      lsize = lsize + 1 ; varlist(lsize) = h_char
      if (wind_type == 1)then
      lsize = lsize + 1 ; varlist(lsize) = wu_char
      lsize = lsize + 1 ; varlist(lsize) = wv_char
      endif
    elseif(g(n)%space_dim ==3)then
      lsize = lsize + 1 ; varlist(lsize) = u_char
      lsize = lsize + 1 ; varlist(lsize) = v_char
      lsize = lsize + 1 ; varlist(lsize) = zeta_char
      lsize = lsize + 1 ; varlist(lsize) = omega_char
      lsize = lsize + 1 ; varlist(lsize) = h_char
      if (wind_type == 1)then
      lsize = lsize + 1 ; varlist(lsize) = wu_char
      lsize = lsize + 1 ; varlist(lsize) = wv_char
      endif
      lsize = lsize + 1 ; varlist(lsize) = kh_char
      if (g(n)%hdiff_type ==2)then
        lsize = lsize + 1 ; varlist(lsize) = viscofm_char
      endif
    endif
  end do

  !determine number of elements
  msg = "dimension 'nele' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, nele_char, dimid ),msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, N_elems ))

  !determine number of nodes 
  msg = "dimension 'node' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, node_char, dimid ),msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, N_verts ))

  !determine number of layers
  msg = "dimension 'siglay' not in the netcdf dataset"
  call ncdchk(nf90_inq_dimid(fid, siglay_char, dimid ),msg)
  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, N_lay ))
  N_lev = N_lay + 1

  !allocate dataspace
  allocate(xm(N_verts))
  allocate(ym(N_verts))
  allocate(xm0(N_verts))
  allocate(ym0(N_verts))

  allocate(xc0(N_elems))
  allocate(yc0(N_elems))
  allocate(xc(N_elems))
  allocate(yc(N_elems))
  allocate(hm(N_verts))
  allocate(aw0(N_elems,3)) ; aw0 = a3rd
  allocate(awx(N_elems,3)) ; awx = zero
  allocate(awy(N_elems,3)) ; awy = zero
  allocate(tri(N_elems,3))
  allocate(nbe(N_elems,3))
  allocate(siglay(N_verts,N_lay))
  allocate(siglev(N_verts,N_lev))
  allocate(esiglay(N_elems,N_lay))
  allocate(esiglev(N_elems,N_lev))

  allocate(a2u(N_elems,4)) ;a2u   = zero
  allocate(a1u(N_elems,4)) ;a1u   = zero

  !----------------Node, Boundary Condition, and Control Volume-----------------------!

 
  !ALLOCATE(NBE(0:N_elems,3))          ;NBE      = 0  !!INDICES OF ELEMENT NEIGHBORS
  ALLOCATE(NTVE(0:N_verts))           ;NTVE     = 0
  ALLOCATE(ISONB(0:N_verts))          ;ISONB    = 0  !!NODE MARKER = 0,1,2
  ALLOCATE(ISBCE(0:N_elems))          ;ISBCE    = 0




  !read in mesh
 

  msg = "error reading x coordinate"
  call ncdchk( nf90_inq_varid(fid,x_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, xm),msg)
  msg = "error reading y coordinate"
  call ncdchk( nf90_inq_varid(fid,y_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, ym),msg)

  msg = "error reading h coordinate"
  call ncdchk( nf90_inq_varid(fid,h_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, hm),msg)
  msg = "error reading nv coordinate"
  call ncdchk( nf90_inq_varid(fid,nv_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, tri),msg)
  if(ncdscan( nf90_inq_varid(fid,nbe_char,varid),msg ) )then
    msg = "error reading nbe"
    call ncdchk( nf90_inq_varid(fid,nbe_char,varid),msg )
    call ncdchk(nf90_get_var(fid, varid, nbe),msg)
  else
    write(*,*)'WARNING:::::: NBE is not in the forcing file' 
    write(*,*)'will try to compute internally'
  endif
  msg = "error reading aw0"
  !read aw0 if they exist, otherwise use 1st order interpolation
  if(ncdscan( nf90_inq_varid(fid,aw0_char,varid),msg ) )then
    call ncdchk(nf90_get_var(fid, varid, aw0),msg)
    msg = "error reading awx"
    call ncdchk( nf90_inq_varid(fid,awx_char,varid),msg )
    call ncdchk(nf90_get_var(fid, varid, awx),msg)
    msg = "error reading awy"
    call ncdchk( nf90_inq_varid(fid,awy_char,varid),msg )
    call ncdchk(nf90_get_var(fid, varid, awy),msg)
    msg = "error reading a1u"
    call ncdchk( nf90_inq_varid(fid,a1u_char,varid),msg )
    call ncdchk(nf90_get_var(fid, varid, a1u),msg)
    msg = "error reading a2u"
    call ncdchk( nf90_inq_varid(fid,a2u_char,varid),msg )
    call ncdchk(nf90_get_var(fid, varid, a2u),msg)

  else 
    write(*,*)'WARNING:::::: AW0 AWX AWY Do NOT exist in forcing file'
    write(*,*)'Proceeding with: 1st Order interpolation' 
    write(*,*)'AW0 = 1/3; AWX=AWY = 0'
    write(*,*)'In the future, select [grid metrics] in your NetCDF namelist'
  endif
  msg = "error reading siglay"
  call ncdchk( nf90_inq_varid(fid,siglay_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, siglay),msg)
  msg = "error reading siglev"
  call ncdchk( nf90_inq_varid(fid,siglev_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, siglev),msg)

  !read secondary connectivity (nbve/ntve) 
  !grid_metrics = .false.
  !msg = "dimension 'maxelem' not in the netcdf dataset"
  !if(ncdscan( nf90_inq_dimid(fid,'maxelem',dimid),msg ) )then
  !  call ncdchk(nf90_inquire_dimension(fid, dimid, dname, Max_Elems))
  !  allocate(ntve(N_verts))
  !  allocate(nbve(N_verts,Max_Elems))
  !  msg = "error reading ntve"
  !  call ncdchk( nf90_inq_varid(fid,'ntve',varid),msg )
  !  call ncdchk(nf90_get_var(fid, varid, ntve),msg)
  !  msg = "error reading nbve"
  !  call ncdchk( nf90_inq_varid(fid,'nbve',varid),msg )
  !  call ncdchk(nf90_get_var(fid, varid, nbve),msg)
  !  grid_metrics = .true.
  !else 
! Revised by Xinyou Lin  in Jan ,2009
!    write(*,*)'WARNING:::::: NTVE/NBVE Do NOT exist in forcing file'
!    write(*,*)'This will slow down the element search procedure'
!    write(*,*)'Proceeding with: 1st Order interpolation' 
!    write(*,*)'In the future, select [grid metrics] in your NetCDF namelist'

!  endif


  !calculate cell center coordinates
     if(spherical == 1)then
       do i=1,N_verts
       
       xm0(i) = xm(i)
       ym0(i) = ym(i)
       if(xm0(i) >= 0.0_sp .and. xm0(i) <=180.0_sp)then
          xm0(i) = xm0(i) + 180.0_sp
       elseif( xm0(i) > 180.0_sp .and. xm0(i) <=360.0_sp)  then
         xm0(i) = xm0(i) - 180.0_sp
       endif
       enddo
     endif


  do i=1,N_elems
    subset = tri(i,1:3)
    xc(i)  = a3rd*(sum(xm(subset)))
    yc(i)  = a3rd*(sum(ym(subset)))
    if(spherical == 1)then
    xc0(i)  = a3rd*(sum(xm0(subset)))
    yc0(i)  = a3rd*(sum(ym0(subset)))
    endif
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
  ! calculate secondary connectivity (nbve/ntve)
  !determine nbve/ntve - secondary connectivity, used
  !for searching element containing point

  !mark boundary elements
  write(*,*)'tri grid edge'
  CALL TRIANGLE_GRID_EDGE
  write(*,*)'done tge' 
  grid_metrics = .true.



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
  integer,  pointer :: cell(:)
  integer  :: i,np
  real(sp) :: tscale
  real(sp), allocatable :: pdxt(:), pdyt(:)
  
  !set pointers to x,y particle positions
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('status',g,istatus)
  call get_state('cell',g,cell)

  !set dimensions for loops and time step
  np = g%nind

  allocate(pdxt(np))  ;
  allocate(pdyt(np))  ;
  !set diffusive time scale
  tscale = sqrt(2.*dT*g%hdiff_const_val)

  !horizontal random walk
     if(spherical == 0 )then
      where(istatus == ACTIVE)
       pdxt = x + normal()*tscale
       pdyt = y + normal()*tscale
      end where
     elseif (spherical == 1)then
      where(istatus == ACTIVE)
       pdxt = x  + normal()*tscale/(tpi*COS(y) + 1.0E-6)
       pdyt = y  + normal()*tscale/tpi
      end where

      where( pdxt < 0.0_SP)
      pdxt = pdxt + 360.0_SP
      end where
      where( pdxt > 360.0_SP)
      pdxt = pdxt - 360.0_SP
      end where

      where( pdyt > 90.0_SP)
      pdyt = 180.0_SP - pdyt
      end where
      where( pdyt < -90.0_SP)
      pdyt =  - 180.0_SP - pdyt
      end where
     endif


     

  call find_element(np,pdxt,pdyt,cell,istatus)
  !!--Update Only Particle Still in Water
  
    where(istatus==ACTIVE)
      x  = pdxt
      y  = pdyt
    end where
  !!--reset position of particles which are lost from domain to last known position
    where(istatus==EXITED)
      istatus=ACTIVE
    end where

  !nullify pointers
  nullify(x)
  nullify(y)
  nullify(istatus)
  nullify(cell)
  deallocate(pdxt)
  deallocate(pdyt)

end subroutine rw_hdiff_constant


!----------------------------------------------------
! Random-Walk horizontal diffusion with spatially
!   variable turbulent eddy diffusivity
!
! Use eddy diffusivity from the model (viscofm)
! Use Visser's naive random walk to compute step
!----------------------------------------------------
subroutine rw_hdiff_variable(g, dT)
  use utilities, only : normal
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: dT
  !----------------------------
  integer,  pointer :: istatus(:)
  integer,  pointer :: cell(:)
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  real(sp), pointer :: s(:)
  real(sp), allocatable :: viscofm(:), pdxt(:), pdyt(:)
  real(sp), allocatable :: tscale(:)
  integer  :: i,np

  !set problem size and time step
  np = g%nind

  !set pointers to particle positions and status
  call get_state('status',g,istatus)
  call get_state('cell',g,cell)
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('s',g,s)
  !allocate local data
  allocate(viscofm(np))  ; viscofm   = zero 
  allocate(tscale(np))  ; tscale   = zero 
  allocate(pdxt(np))  ;
  allocate(pdyt(np))  ;

    !evaluate kh at both locations
    call interp(np,x,y,s,cell,istatus,viscofm_char,viscofm,3)

    !update particle position using Visser modified random walk 
     tscale = sqrt(2.*dT*viscofm)
  !horizontal random walk
     if(spherical == 0 )then
      where(istatus == ACTIVE)
       pdxt = x + normal()*tscale
       pdyt = y + normal()*tscale
      end where
     elseif (spherical == 1)then
      where(istatus == ACTIVE)
       pdxt = x  + normal()*tscale/(tpi*COS(y) + 1.0E-6)
       pdyt = y  + normal()*tscale/tpi
      end where

      where( pdxt < 0.0_SP)
      pdxt = pdxt + 360.0_SP
      end where
      where( pdxt > 360.0_SP)
      pdxt = pdxt - 360.0_SP
      end where

      where( pdyt > 90.0_SP)
      pdyt = 180.0_SP - pdyt
      end where
      where( pdyt < -90.0_SP)
      pdyt =  - 180.0_SP - pdyt
      end where
     endif


     

  call find_element(np,pdxt,pdyt,cell,istatus)
  !!--Update Only Particle Still in Water
  
    where(istatus==ACTIVE)
      x  = pdxt
      y  = pdyt
    end where
  !!--reset position of particles which are lost from domain to last known position
    where(istatus==EXITED)
      istatus=ACTIVE
    end where



  !deallocate workspace and nullify pointers
  deallocate(viscofm)
  deallocate(tscale)
  deallocate(pdxt)
  deallocate(pdyt)
  nullify(x)
  nullify(y)
  nullify(s)
  nullify(istatus)
  nullify(cell)
  





end subroutine rw_hdiff_variable

!----------------------------------------------------
! Random-Walk vertical diffusion 
!
!   - use eddy diffusivity from the model (kh)
!   - use Vissers modified random walk to compute jump
!----------------------------------------------------
subroutine rw_vdiff(g, dT, nstep)
  use utilities, only : normal,unitrand,ran1
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: dT
  integer,  intent(in) :: nstep
  !----------------------------
  integer,  pointer :: istatus(:)
  integer,  pointer :: cell(:)
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  real(sp), pointer :: s(:)
  real(sp), pointer :: z(:)
  real(sp), pointer :: h(:)
  real(sp), allocatable :: kh(:)
  real(sp), allocatable :: kh2(:)
  real(sp), allocatable :: ds(:)
  real(sp), allocatable :: dkh_ds(:)
  real(sp), allocatable :: zeta(:)
  real(sp), allocatable :: s_shift(:)
  real(sp), parameter :: delta_s = 0.001
  real(sp) :: deltaT,fac,randy,dz,depth,dkh_dz,vsink
  integer  :: n,p,np

  !set problem size and time step
  vsink = g%vsink
  np = g%nind
  deltaT = dT/float(nstep)

  !set pointers to particle positions and status
  call get_state('status',g,istatus)
  call get_state('cell',g,cell)
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('s',g,s)
  call get_state('z',g,z)
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

  call interp(np,x,y,cell,istatus,h_char,h,3)
  call interp(np,x,y,cell,istatus,zeta_char,zeta,3)

  ! ==> loop over substeps
  do n=1,nstep
   
    !--------------------------------------------------
    ! calculate d(kh)/d(s) - brute force
    !--------------------------------------------------

    !set derivative step (backward near free surface)
    ds = delta_s
    where(s+delta_s > 0)ds = -delta_s

    !evaluate kh at both locations
    call interp(np,x,y,s,cell,istatus,kh_char,kh,3)
    call interp(np,x,y,s+ds,cell,istatus,kh_char,kh2,3)

    !form the derivative d(kh)/d(s)
    dkh_ds = (kh2-kh)/ds

    !function evaluation at [z + 0.5*dkh/dz*deltat] - Visser
    s_shift = s + ahalf*dkh_ds*deltaT/((h+zeta)**2)
    call interp(np,x,y,s_shift,cell,istatus,kh_char,kh,3)

    ! => main loop over particles
!gwc    do p=1,1000
!gwc      write(33,*)normal(),ran1()
!gwc    end do
!gwc    stop
    do p=1,np
      if(istatus(p) < 1)cycle

      !update particle position using Visser modified random walk 
      depth  = h(p)+zeta(p)
      dkh_dz = dkh_ds(p)/depth
      dz     = dkh_dz*deltaT + normal(-deltaT*vsink*depth,1.0d0)*sqrt(fac*kh(p)) !Visser-modified
!      dz     = dkh_dz*deltaT + unitrand()*sqrt(fac*kh(p)) !Visser-modified
!      dz     = 2*unitrand()*sqrt(2.*3*kh(p)*deltaT)  !naive unitrand
!      dz     = normal(-deltaT*.001*depth,1.d0)*sqrt(2.*kh(p)*deltaT)                 !naive normal
      s(p)   = s(p) + dz/depth

      !set boundary conditions at free surface and bottom
      s(p) = max(s(p),-(2.0+s(p))) !reflect off bottom
      s(p) = min(s(p),0.0)         !don't pierce free surface

    end do
    ! <= end particle loop

  end do
  ! <== end loop over substeps

  !--Calculate Particle Location in Cartesian Vertical Coordinate----------------!
  z = s*(h+zeta)  + zeta

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
  nullify(z)
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

subroutine sz_trans(np,g)
  integer, intent(in) :: np
  type(igroup), intent(inout) :: g
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  real(sp), pointer :: s(:)
  real(sp), pointer :: z(:)
  real(sp), pointer :: h(:)
  real(sp), pointer :: zeta(:)
  integer , pointer :: cell(:)
  integer , pointer :: istatus(:)
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('s',g,s)
  call get_state('z',g,z)
  call get_state('h',g,h)
  call get_state('zeta',g,zeta)
  call get_state('status',g,istatus)
  call get_state('cell',g,cell)

  call interp(np,x,y,cell,istatus,zeta_char,zeta,3)
  call interp(np,x,y,cell,istatus,h_char,h,3)
  if(sz_cor == 1)then
!    z  = -z + zeta 
     s  = (z - zeta)/(h + zeta)
  elseif(sz_cor == 0)then
     z  = s*(h + zeta) + zeta
  endif
  nullify(x)
  nullify(y)
  nullify(s)
  nullify(z)
  nullify(h)
  nullify(cell)
  nullify(zeta)
  nullify(istatus)



end subroutine sz_trans

!---------------------------------------------------
! 2-D Advection 
!----------------------------------------------------
subroutine advect2D(g,deltaT,np)
  integer, intent(in) :: np
  integer  :: k,i,ns
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: deltaT
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  real(sp), pointer :: h(:)
  integer , pointer :: cell(:)
  integer , pointer :: istatus(:)
  
  real(sp), dimension(np) :: u,u1,u2,v,v1,v2
  real(sp), dimension(np) :: pdx,pdy
  real(sp), dimension(np) :: pdxt,pdyt

  real(sp), dimension(np,0:mstage) :: chix,chiy  
  real(sp), parameter              :: eps  = 1.0E-5

!!!!!!!

  !set dimensions for loops and time step
  !np = g%nind
  !set pointers to states
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('h',g,h)
  call get_state('cell',g,cell)
  call get_state('status',g,istatus)
  !--Initialize Stage Functional Evaluations
  chix = 0.0_sp
  chiy = 0.0_sp
  !--Loop over RK Stages
  do ns=1,mstage

     !!Particle Position at Stage N (x,y,sigma)
     pdx(:) = x(:)  + a_rk(ns)*deltaT*chix(:,ns-1)
     pdy(:) = y(:)  + a_rk(ns)*deltaT*chiy(:,ns-1)
     !!Calculate Velocity Field for Stage N Using C_RK Coefficients
     !interpolate velocity field to particle position
  call interp(np,pdx,pdy,cell,istatus,ua_char,u1,3)
  call interp(np,pdx,pdy,cell,istatus,va_char,v1,3)
  call interp(np,pdx,pdy,cell,istatus,ua_char,u2,4)
  call interp(np,pdx,pdy,cell,istatus,va_char,v2,4)
   u  = (1.0_sp-c_rk(ns))*u1 + c_rk(ns)*u2
   v  = (1.0_sp-c_rk(ns))*v1 + c_rk(ns)*v2
   chix(:,ns)  = u(:)
   chiy(:,ns)  = v(:)

end do

  !--Sum Stage Contributions to get Updated Particle Positions-------------------!
  pdxt(:)  = x(:)
  pdyt(:)  = y(:)
  do ns=1,mstage
     pdxt(:) = pdxt(:) + deltaT*chix(:,ns)*b_rk(ns)*FLOAT(istatus(:))
     pdyt(:) = pdyt(:) + deltaT*chiy(:,ns)*b_rk(ns)*FLOAT(istatus(:))
  end do
  call find_element(np,pdxt,pdyt,cell,istatus)
  !!--Update Only Particle Still in Water
  
    where(istatus==ACTIVE)
      x  = pdxt
      y  = pdyt
    end where
  !!--reset position of particles which are lost from domain to last known position
    where(istatus==EXITED)
      istatus=ACTIVE
    end where

 ! call find_element(np,x,y,cell,istatus)
  !disassociate pointers
  nullify(x)
  nullify(y)
  nullify(h)
  nullify(cell)
  nullify(istatus)

end subroutine advect2D
!---------------------------------------------------
! 3-D Advection 
!----------------------------------------------------
subroutine advect3D(g,deltaT,np,time)
  integer, intent(in) :: np
  integer  :: k,i,ns,ni
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: deltaT,time
  real(sp), pointer :: x(:)
  real(sp), pointer :: y(:)
  real(sp), pointer :: s(:)
  real(sp), pointer :: z(:)
  real(sp), pointer :: h(:)
  real(sp), pointer :: zeta(:)
  integer , pointer :: cell(:)
  integer , pointer :: istatus(:)
  
  real(sp), dimension(np) :: u,u1,u2,v,v1,v2,w,w1,w2,wm
  real(sp), dimension(np) :: zeta1,zeta2,pdx,pdy,pdz
  real(sp), dimension(np) :: wu,wu1,wu2,wv,wv1,wv2
  real(sp), dimension(np) :: pdxt,pdyt,pdzt
  real(sp), dimension(np) :: zn

  real(sp), dimension(np,0:mstage) :: chix,chiy,chiz  
  real(sp), parameter              :: eps  = 1.0E-5
  real(sp)  :: diel
!!!!!!!

  diel=time/3600.0-int(time/3600.0/24.0)*24.0
  !set dimensions for loops and time step
  !np = g%nind
  !set pointers to states
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('s',g,s)
  call get_state('z',g,z)  
  call get_state('h',g,h)
  call get_state('zeta',g,zeta)
  call get_state('cell',g,cell)
  call get_state('status',g,istatus)

  !disassociate pointers
  !--Initialize Stage Functional Evaluations
  chix = 0.0_sp
  chiy = 0.0_sp
  chiz = 0.0_sp
  pdx  =x
  pdy  =y
  pdz  =s
  !--Loop over RK Stages
  do ns=1,mstage

     !!Particle Position at Stage N (x,y,sigma)
     
     if (spherical == 0)then 
     pdx(:) = x(:)  + a_rk(ns)*deltaT*chix(:,ns-1)
     pdy(:) = y(:)  + a_rk(ns)*deltaT*chiy(:,ns-1)
     pdz(:) = s(:)  + a_rk(ns)*deltaT*chiz(:,ns-1)
!!!!!
     elseif (spherical == 1)then 
     pdx(:) = x(:)  + a_rk(ns)*deltaT*chix(:,ns-1)/(tpi*COS(pdy(:)) + 1.0E-6)
     pdy(:) = y(:)  + a_rk(ns)*deltaT*chiy(:,ns-1)/tpi
     pdz(:) = s(:)  + a_rk(ns)*deltaT*chiz(:,ns-1)
          
     where( pdx < 0.0_SP)
     pdx = pdx + 360.0_SP
     end where
     where( pdx > 360.0_SP)
     pdx = pdx - 360.0_SP
     end where

     where( pdy > 90.0_SP)
     pdy = 180.0_SP - pdy 
     end where
     where( pdy < -90.0_SP)
     pdy =  - 180.0_SP - pdy
     end where

     endif
!!!!!
     !!Adjust Sigma Position to Reflect Off Bottom (Mirroring)
     pdz = max(pdz,-(2.0+pdz))

     !!Adjust Sigma Position to Remain Below Free Surface
     pdz = min(pdz,0.0_sp)

     !!Calculate Velocity Field for Stage N Using C_RK Coefficients
     !interpolate velocity field to particle position

  !write(*,*)np,pdx,pdy,pdz,cell,istatus,u_char,u1
  call interp(np,pdx,pdy,pdz,cell,istatus,u_char,u1,3)
  call interp(np,pdx,pdy,pdz,cell,istatus,v_char,v1,3)
  call interp(np,pdx,pdy,pdz,cell,istatus,omega_char,w1,3) !wts means omega
  call interp(np,pdx,pdy,pdz,cell,istatus,u_char,u2,4)
  call interp(np,pdx,pdy,pdz,cell,istatus,v_char,v2,4)
  call interp(np,pdx,pdy,pdz,cell,istatus,omega_char,w2,4)
  call interp(np,pdx,pdy,cell,istatus,h_char,h,3)
  call interp(np,pdx,pdy,cell,istatus,zeta_char,zeta1,3)
  call interp(np,pdx,pdy,cell,istatus,zeta_char,zeta2,4)
!!!!!88!!!!!
     if (wind_type == 1)then
  call interp(np,pdx,pdy,cell,istatus,wu_char,wu1,3)
  call interp(np,pdx,pdy,cell,istatus,wv_char,wv1,3)
  call interp(np,pdx,pdy,cell,istatus,wu_char,wu2,4)
  call interp(np,pdx,pdy,cell,istatus,wv_char,wv2,4)
   wu  = (1.0_sp-c_rk(ns))*wu1 + c_rk(ns)*wu2
   wv  = (1.0_sp-c_rk(ns))*wv1 + c_rk(ns)*wv2

!!!!!!!!!!
   u  = (1.0_sp-c_rk(ns))*u1 + c_rk(ns)*u2  
   v  = (1.0_sp-c_rk(ns))*v1 + c_rk(ns)*v2  
   u  = u + wu*0.02
   v  = v + wv*0.02
     elseif (wind_type == 0) then
   u  = (1.0_sp-c_rk(ns))*u1 + c_rk(ns)*u2
   v  = (1.0_sp-c_rk(ns))*v1 + c_rk(ns)*v2

     endif
   w  = (1.0_sp-c_rk(ns))*w1 + c_rk(ns)*w2
   zeta  = (1.0_sp-c_rk(ns))*zeta1  + c_rk(ns)*zeta2

!Added by Xinyou Lin for DVM modelling:WP=WP+WM
   if(dvm_bio == 1.0)then
      do ni = 1, np
       if(6.0 <= diel .and. diel <18.0)then
        pdzt(ni) = (1 + pdz(ni))*(h(ni)+zeta(ni))  
!        w(ni) = w(ni) - 0.0075*tanh((pdzt(ni) -  2)*3.14159)
!         w(ni) = w(ni) - 0.0075*tanh((pdzt(ni) -  10))
         w(ni) = w(ni) - 0.0075*tanh((pdzt(ni) -  dvmh_dn)) !from bottom
       else
        pdzt(ni) = -pdz(ni)*(h(ni)+zeta(ni))
        w(ni) = w(ni) + 0.0075*tanh((pdzt(ni) - dvmh_up))   !from surface
       endif
       enddo
   endif
!///////////////////////

   chix(:,ns)  = u(:)
   chiy(:,ns)  = v(:)
   chiz(:,ns)  = w(:)/(h(:)+zeta(:))  !delta_sigma/deltaT = ww/D
!   chiz(:,ns) = chiz(:,ns) - .005 !!gwc try sinking - this is CFL-constrained
     !!Limit vertical motion in very shallow water
     where( (h + zeta) < eps)
        chiz(:,ns) = 0.0_sp
     end where

end do

  !--Sum Stage Contributions to get Updated Particle Positions-------------------!
  pdxt(:)  = x(:)
  pdyt(:)  = y(:)
  pdzt(:)  = s(:)
  do ns=1,mstage
     if (spherical == 0)then 
     pdxt(:) = pdxt(:) + deltaT*chix(:,ns)*b_rk(ns)*FLOAT(istatus(:))
     pdyt(:) = pdyt(:) + deltaT*chiy(:,ns)*b_rk(ns)*FLOAT(istatus(:))
     pdzt(:) = pdzt(:) + deltaT*chiz(:,ns)*b_rk(ns)*FLOAT(istatus(:))
!!!!!
     elseif (spherical == 1)then 
     pdxt(:) = pdxt(:)  + a_rk(ns)*deltaT*chix(:,ns-1)/(tpi*COS(pdy(:)) + 1.0E-6)
     pdyt(:) = pdyt(:)  + a_rk(ns)*deltaT*chiy(:,ns-1)/tpi
     pdzt(:) = pdzt(:)  + a_rk(ns)*deltaT*chiz(:,ns-1)
          
     where( pdxt < 0.0_SP)
     pdxt = pdxt + 360.0_SP
     end where
     where( pdxt > 360.0_SP)
     pdxt = pdxt - 360.0_SP
     end where

     where( pdyt > 90.0_SP)
     pdyt = 180.0_SP - pdyt 
     end where
     where( pdyt < -90.0_SP)
     pdyt =  - 180.0_SP - pdyt
     end where

     endif
!!!!!
  end do
  call find_element(np,pdxt,pdyt,cell,istatus)
  !!--Update Only Particle Still in Water
  
    where(istatus==ACTIVE)
      x  = pdxt
      y  = pdyt
    end where
  !!--reset position of particles which are lost from domain to last known position
    where(istatus==EXITED)
      istatus=ACTIVE
    end where
  zn = s*(h + zeta) + zeta  
  s(:)  = pdzt(:) 
!  s(:)  = s(:) - .05*deltaT/(h + zeta) !try constant dz sink
  !--Adjust Depth of Updated Particle Positions----------------------------------!
!  s = max(s,-(2.0+s))                 !Reflect off Bottom
  s = max(s,-1.0)                       !Do Not Pierce Bottom  
!  s = -1.0 !gwc  
  s = min(s,0.0_sp)                      !Don t Pierce Free Surface

  !--Evaluate Bathymetry and Free Surface Height at Updated Particle Position----!
  call interp(np,x,y,cell,istatus,h_char,h,3)
  call interp(np,x,y,cell,istatus,zeta_char,zeta,4)

  !--Sigma adjustment if fixed depth tracking------------------------------------!
  if(fix_dep == 1)then
      if(sz_cor == 1)then   ! when initial depth is specified in z
     ! s = (-zpini)/(h+zeta)  !  THIS IS REALLY PDZN = ((-LAG%ZPIN+EP) - EP)/(HP+EP)
     !                        !  WHERE ZPINI IS THE SPECIFIED FIXED DEPTH RELATIVE TO THE SS
      s  = zpini/(h + zeta)
      elseif(sz_cor == 0)then ! when initial depth is specified in s
      s  = zpini              !  WHERE ZPINI IS THE SPECIFIED FIXED SIGMA
      endif

      s = max(s,-1.0_SP)     ! Depth can change though if particle goes into shallower areas

  endif


  !--Calculate Particle Location in Cartesian Vertical Coordinate----------------!
! z = h !s*(h+zeta)  + zeta
 z = s*(h+zeta)  + zeta

  
  !disassociate pointers
  nullify(x)
  nullify(y)
  nullify(s)
  nullify(z)
  nullify(h)
  nullify(zeta)
  nullify(cell)
  nullify(istatus)

end subroutine advect3D

!----------------------------------------------------
! interpolation routine 
!   interpolate scalar field [vname] to vector [v] 
!   at particle position [x,y] in cell [cell]
!   if [istatus] < 0 (inactive) do not interpolate
!----------------------------------------------------
subroutine interp_float2D(np,x,y,cell,istatus,vname,v,iframe)
  integer, intent(in)    :: np,iframe
  real(sp),intent(in)    :: x(np)
  real(sp),intent(in)    :: y(np)
  integer, intent(in)    :: cell(np)
  integer, intent(in)    :: istatus(np)
  character(len=*)       :: vname
  real(sp),intent(inout) :: v(np)
  !-------------------------------
  real(sp), pointer :: field(:)  !pointer to FVCOM field 
  integer :: i,vts(3),icell,d1
  !added by xinyou
  integer  :: e1,e2,e3,n1,n2,n3
  real(sp) :: xoc,yoc,dvdx,dvdy
  real(sp) :: e0,ex,ey,tmpx
  !determine dimension of forcing variable

  !point to forcing variable
  call get_forcing(trim(vname),field,iframe) 

  !initialize v
  v = 0.0

  !determine the dimensions
  d1 = size(field,1)

  !interpolate to element-based quantities
  !gwc debug = 0th order
  if(d1 == N_elems)then  
    do i=1,np
      icell = cell(i)
      if(istatus(i) < 1 .or. icell == 0)cycle 
      e1  = nbe(icell,1)
      e2  = nbe(icell,2)
      e3  = nbe(icell,3)
        if (spherical == 0)then
      xoc = x(i) - xc(icell)
      yoc = y(i) - yc(icell)
        elseif(spherical == 1)then
            
          if(x(i) >=90.0 .and. x(i) <=270.0)then
          xoc = x(i) - xc(icell)
          yoc = y(i) - yc(icell)
          else
              if(x(i) >= 0.0_sp .and. x(i) <=180.0_sp)then
                tmpx = x(i) + 180.0_sp
              elseif( x(i) > 180.0_sp .and. x(i) <=360.0_sp)  then
                tmpx = x(i) - 180.0_sp
              endif
          xoc = tmpx - xc0(icell)
          yoc = y(i) - yc0(icell)
          endif
        endif
        
     dvdx = a1u(icell,1)*field(icell)+a1u(icell,2)*field(e1)     &
          + a1u(icell,3)*field(e2)   +a1u(icell,4)*field(e3)
     dvdy = a2u(icell,1)*field(icell)+a2u(icell,2)*field(e1)     &
          + a2u(icell,3)*field(e2)   +a2u(icell,4)*field(e3)
     v(i) = field(icell) + dvdx*xoc + dvdy*yoc

!      v(i) = field(icell)
!      v(i) =interp_from_elems(icell,x(i),y(i))     
    end do
  elseif(d1 == N_verts)then !vertex-based 2-D array
    do i=1,np
      icell = cell(i)
      if(istatus(i) < 1 .or. icell == 0)cycle 
     n1  = tri(icell,1)
     n2  = tri(icell,2)
     n3  = tri(icell,3)
     
        if (spherical == 0)then
      xoc = x(i) - xc(icell)
      yoc = y(i) - yc(icell)
        elseif(spherical == 1)then

          if(x(i) >=90.0 .and. x(i) <=270.0)then
          xoc = x(i) - xc(icell)
          yoc = y(i) - yc(icell)
          else
              if(x(i) >= 0.0_sp .and. x(i) <=180.0_sp)then
                tmpx = x(i) + 180.0_sp
              elseif( x(i) > 180.0_sp .and. x(i) <=360.0_sp)  then
                tmpx = x(i) - 180.0_sp
              endif
          xoc = tmpx - xc0(icell)
          yoc = y(i) - yc0(icell)
          endif
        endif

     !----Linear Interpolation of Free Surface Height---------------------------------!
    
     e0 = aw0(icell,1)*field(n1)+aw0(icell,2)*field(n2)+aw0(icell,3)*field(n3)
     ex = awx(icell,1)*field(n1)+awx(icell,2)*field(n2)+awx(icell,3)*field(n3)
     ey = awy(icell,1)*field(n1)+awy(icell,2)*field(n2)+awy(icell,3)*field(n3)
   v(i) = e0 + ex*xoc + ey*yoc

!      vts  = tri(icell,1:3)
!      v(i) = interp_from_nodes(icell,x(i),y(i),field(vts))
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
subroutine interp_float3D(np,x,y,s,cell,istatus,vname,v,iframe)
  integer, intent(in)    :: np,iframe
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
  real(sp) :: f1,f2,v1,v2,tmpx
  integer  :: vts(3)
  !added by xinyou
  integer  :: e1,e2,e3,k,n1,n2,n3
  real(sp) :: xoc,yoc,dvdx,dvdy,ve01,ve02
  real(sp) :: e0,ex,ey,ep01,ep02
  !point to forcing variable
  !get_forcing is in forcing.f90
  
  call get_forcing(trim(vname),field,iframe) 
  
  !initialize v
  v = 0.0

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
      !xinyou
      e1  = nbe(icell,1)
      e2  = nbe(icell,2)
      e3  = nbe(icell,3)
        if (spherical == 0)then
      xoc = x(i) - xc(icell)
      yoc = y(i) - yc(icell)
        elseif(spherical == 1)then

          if(x(i) >=90.0 .and. x(i) <=270.0)then
          xoc = x(i) - xc(icell)
          yoc = y(i) - yc(icell)
          else
              if(x(i) >= 0.0_sp .and. x(i) <=180.0_sp)then
                tmpx = x(i) + 180.0_sp
              elseif( x(i) > 180.0_sp .and. x(i) <=360.0_sp)  then
                tmpx = x(i) - 180.0_sp
              endif
          xoc = tmpx - xc0(icell)
          yoc = y(i) - yc0(icell)
          endif
        endif


        k = k1
     dvdx = a1u(icell,1)*field(icell,k)+a1u(icell,2)*field(e1,k)     &
          + a1u(icell,3)*field(e2,k)   +a1u(icell,4)*field(e3,k)
     dvdy = a2u(icell,1)*field(icell,k)+a2u(icell,2)*field(e1,k)     &
          + a2u(icell,3)*field(e2,k)   +a2u(icell,4)*field(e3,k)
     ve01 = field(icell,k) + dvdx*xoc + dvdy*yoc
        k = k2
     dvdx = a1u(icell,1)*field(icell,k)+a1u(icell,2)*field(e1,k)     &
          + a1u(icell,3)*field(e2,k)   +a1u(icell,4)*field(e3,k)
     dvdy = a2u(icell,1)*field(icell,k)+a2u(icell,2)*field(e1,k)     &
          + a2u(icell,3)*field(e2,k)   +a2u(icell,4)*field(e3,k)
     ve02 = field(icell,k) + dvdx*xoc + dvdy*yoc
 

       v(i) = f1*ve01 + f2*ve02
!      v1 = interp_from_elems(icell,x(i),y(i),field(vts,k1))
!      v2 = interp_from_elems(icell,x(i),y(i),field(vts,k2))
      
      !interpolate
!      v(i) = f1*field(icell,k1) + f2*field(icell,k2)
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
     n1  = tri(icell,1)
     n2  = tri(icell,2)
     n3  = tri(icell,3)
        if (spherical == 0)then
      xoc = x(i) - xc(icell)
      yoc = y(i) - yc(icell)
        elseif(spherical == 1)then

          if(x(i) >=90.0 .and. x(i) <=270.0)then
          xoc = x(i) - xc(icell)
          yoc = y(i) - yc(icell)
          else
              if(x(i) >= 0.0_sp .and. x(i) <=180.0_sp)then
                tmpx = x(i) + 180.0_sp
              elseif( x(i) > 180.0_sp .and. x(i) <=360.0_sp)  then
                tmpx = x(i) - 180.0_sp
              endif
          xoc = tmpx - xc0(icell)
          yoc = y(i) - yc0(icell)
          endif
        endif

     !----Linear Interpolation of Free Surface Height---------------------------------!
      k = k1
     e0 = aw0(icell,1)*field(n1,k)+aw0(icell,2)*field(n2,k)+aw0(icell,3)*field(n3,k)
     ex = awx(icell,1)*field(n1,k)+awx(icell,2)*field(n2,k)+awx(icell,3)*field(n3,k)
     ey = awy(icell,1)*field(n1,k)+awy(icell,2)*field(n2,k)+awy(icell,3)*field(n3,k)
   ep01 = e0 + ex*xoc + ey*yoc
      k = k2
     e0 = aw0(icell,1)*field(n1,k)+aw0(icell,2)*field(n2,k)+aw0(icell,3)*field(n3,k)
     ex = awx(icell,1)*field(n1,k)+awx(icell,2)*field(n2,k)+awx(icell,3)*field(n3,k)
     ey = awy(icell,1)*field(n1,k)+awy(icell,2)*field(n2,k)+awy(icell,3)*field(n3,k)
   ep02 = e0 + ex*xoc + ey*yoc
   v(i) = f1*ep01 + f2*ep02
      !set vertices surrounding triangle
      !vts  = tri(icell,1:3)
      !v1 = interp_from_nodes(icell,x(i),y(i),field(vts,k1))
      !v2 = interp_from_nodes(icell,x(i),y(i),field(vts,k2))
      !v(i) = f1*v1 + f2*v2
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
    if(cell(p) /= 0)cycle

    !failed, update status to lost
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
   ! N_lev = N_lay + 1
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
       f1 = (my_sloc-esiglay(cell,k2))/ds 
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

 subroutine   TRIANGLE_GRID_EDGE
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: TEMP,NB_TMP,CELLS,NBET
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: CELLCNT
  INTEGER                              :: I,J,II,JJ,NTMP,NCNT,NFLAG,JJB
  INTEGER                              :: N1,N2,N3,J1,J2,J3,MX_NBR_ELEM 

  !----------------------------INITIALIZE----------------------------------------!

  ISBCE = 0
  ISONB = 0
  NBE   = 0
  !----DETERMINE NBE(i=1:n,j=1:3): INDEX OF 1 to 3 NEIGHBORING ELEMENTS----------!
  
  ALLOCATE(NBET(N_elems,3)) ; NBET = 0
  ALLOCATE(CELLS(N_verts,50)) ; CELLS = 0
  ALLOCATE(CELLCNT(N_verts))  ; CELLCNT = 0
  DO I=1,N_elems
     N1 = tri(I,1) ; CELLCNT(N1) = CELLCNT(N1)+1
     N2 = tri(I,2) ; CELLCNT(N2) = CELLCNT(N2)+1
     N3 = tri(I,3) ; CELLCNT(N3) = CELLCNT(N3)+1
     CELLS(tri(I,1),CELLCNT(N1)) = I
     CELLS(tri(I,2),CELLCNT(N2)) = I
     CELLS(tri(I,3),CELLCNT(N3)) = I
  END DO
  DO I=1,N_elems
     N1 = tri(I,1)
     N2 = tri(I,2)
     N3 = tri(I,3)
     DO J1 = 1,CELLCNT(N1)
        DO J2 = 1,CELLCNT(N2)
           IF((CELLS(N1,J1) == CELLS(N2,J2)).AND. CELLS(N1,J1) /= I)NBE(I,3) = CELLS(N1,J1)
        END DO
     END DO
     DO J2 = 1,CELLCNT(N2)
        DO J3 = 1,CELLCNT(N3)
           IF((CELLS(N2,J2) == CELLS(N3,J3)).AND. CELLS(N2,J2) /= I)NBE(I,1) = CELLS(N2,J2)
        END DO
     END DO
     DO J1 = 1,CELLCNT(N1)
        DO J3 = 1,CELLCNT(N3)
           IF((CELLS(N1,J1) == CELLS(N3,J3)).AND. CELLS(N1,J1) /= I)NBE(I,2) = CELLS(N3,J3)
        END DO
     END DO
  END DO
 DEALLOCATE(CELLS,CELLCNT)
  !   IF(MSR)WRITE(IPT,*)  '!  NEIGHBOR FINDING      :    COMPLETE'
  !
  !--ENSURE ALL ELEMENTS HAVE AT LEAST ONE NEIGHBOR------------------------------!
  !
  NFLAG = 0
  DO I=1,N_elems
     IF(SUM(NBE(I,1:3))==0)THEN
        NFLAG = 1
        WRITE(*,*)'ELEMENT ',I,' AT ',XC(I),YC(I),' HAS NO NEIGHBORS'
        STOP
     END IF
  END DO
  IF(NFLAG == 1) STOP
  !
  !----IF ELEMENT ON BOUNDARY SET ISBCE(I)=1 AND ISONB(J)=1 FOR BOUNDARY NODES J-!
  !

  DO I=1,N_elems
     IF(MIN(NBE(I,1),NBE(I,2),NBE(I,3))==0)THEN    !!ELEMENT ON BOUNDARY
        ISBCE(I) = 1
        IF(NBE(I,1) == 0)THEN
           ISONB(tri(I,2)) = 1 ; ISONB(tri(I,3)) = 1
        END IF
        IF(NBE(I,2) ==0) THEN
           ISONB(tri(I,1)) = 1 ; ISONB(tri(I,3)) = 1
        END IF
        IF(NBE(I,3) ==0) THEN
           ISONB(tri(I,1)) = 1 ; ISONB(tri(I,2)) = 1
        END IF
     END IF
  END DO

  !==============================================================================|
  !             DEFINE NTVE, NBVE, NBVT                                          !
  !                                                                              !
  ! ntve(1:m):           total number of the surrounding triangles               !
  !                      connected to the given node                             !
  ! nbve(1:m, 1:ntve+1): the identification number of surrounding                !
  !                      triangles with a common node (counted clockwise)        !
  ! nbvt(1:m,ntve(1:m)): the idenfication number of a given node over            !
  !                      each individual surrounding triangle(counted            !
  !                      clockwise)                                              !
  !==============================================================================|

  !
  !----DETERMINE MAX NUMBER OF SURROUNDING ELEMENTS------------------------------!
  !
  MX_NBR_ELEM = 0
  DO I=1,N_verts
     NCNT = 0
     DO J=1,N_elems
        IF( FLOAT(tri(J,1)-I)*FLOAT(tri(J,2)-I)*FLOAT(tri(J,3)-I) == 0.0_SP) &
             NCNT = NCNT + 1
     END DO
     MX_NBR_ELEM = MAX(MX_NBR_ELEM,NCNT)
  END DO

  !
  !----ALLOCATE ARRAYS BASED ON MX_NBR_ELEM--------------------------------------!
  !
  ALLOCATE(NBVE(N_verts,MX_NBR_ELEM+1))
  ALLOCATE(NBVT(N_verts,MX_NBR_ELEM+1))
  !
  !--DETERMINE NUMBER OF SURROUNDING ELEMENTS FOR NODE I = NTVE(I)---------------!
  !--DETERMINE NBVE - INDICES OF NEIGHBORING ELEMENTS OF NODE I------------------!
  !--DETERMINE NBVT - INDEX (1,2, or 3) OF NODE I IN NEIGHBORING ELEMENT---------!
  !
  DO I=1,N_verts
     NCNT=0
     DO J=1,N_elems
        IF (FLOAT(tri(J,1)-I)*FLOAT(tri(J,2)-I)*FLOAT(tri(J,3)-I) == 0.0_SP)THEN
           NCNT = NCNT+1
           NBVE(I,NCNT)=J
           IF((tri(J,1)-I) == 0) NBVT(I,NCNT)=1
           IF((tri(J,2)-I) == 0) NBVT(I,NCNT)=2
           IF((tri(J,3)-I) == 0) NBVT(I,NCNT)=3
        END IF
     ENDDO
     NTVE(I)=NCNT
  ENDDO

  ALLOCATE(NB_TMP(N_verts,MX_NBR_ELEM+1))
  DO I=1,N_verts
     IF(ISONB(I) == 0) THEN
        NB_TMP(1,1)=NBVE(I,1)
        NB_TMP(1,2)=NBVT(I,1)
        DO J=2,NTVE(I)+1
           II=NB_TMP(J-1,1)
           JJ=NB_TMP(J-1,2)
           NB_TMP(J,1)=NBE(II,JJ+1-INT((JJ+1)/4)*3)
           JJ=NB_TMP(J,1)
           IF((tri(JJ,1)-I) == 0) NB_TMP(J,2)=1
           IF((tri(JJ,2)-I) == 0) NB_TMP(J,2)=2
           IF((tri(JJ,3)-I) == 0) NB_TMP(J,2)=3
        ENDDO

        DO J=2,NTVE(I)+1
           NBVE(I,J)=NB_TMP(J,1)
        ENDDO

        DO J=2,NTVE(I)+1
           NBVT(I,J)=NB_TMP(J,2)
        ENDDO

        NTMP=NTVE(I)+1
        IF(NBVE(I,1) /= NBVE(I,NTMP)) THEN
           PRINT*, I,'NBVE(I) NOT CORRECT!!'
           STOP
        ENDIF
        IF(NBVT(I,1) /= NBVT(I,NTMP)) THEN
           PRINT*, I,'NBVT(I) NOT CORRECT!!'
           STOP
        END IF
     ELSE
        JJB=0

        DO J=1,NTVE(I)
           JJ=NBVT(I,J)
           IF(NBE(NBVE(I,J),JJ+2-INT((JJ+2)/4)*3) == 0) THEN
              JJB=JJB+1
              NB_TMP(JJB,1)=NBVE(I,J)
              NB_TMP(JJB,2)=NBVT(I,J)
           END IF
        ENDDO

        IF(JJB /= 1) THEN
           PRINT*, 'ERROR IN ISONB !,I,J', I,J
           !PAUSE
        END IF

        DO J=2,NTVE(I)
           II=NB_TMP(J-1,1)
           JJ=NB_TMP(J-1,2)
           NB_TMP(J,1)=NBE(II,JJ+1-INT((JJ+1)/4)*3)
           JJ=NB_TMP(J,1)
           IF((tri(JJ,1)-I) == 0) NB_TMP(J,2)=1
           IF((tri(JJ,2)-I) == 0) NB_TMP(J,2)=2
           IF((tri(JJ,3)-I) == 0) NB_TMP(J,2)=3
        ENDDO

        DO J=1,NTVE(I)
           NBVE(I,J)=NB_TMP(J,1)
           NBVT(I,J)=NB_TMP(J,2)
        ENDDO
        NBVE(I,NTVE(I)+1)=0

     END IF
  END DO
  DEALLOCATE(NB_TMP)

  RETURN


 end subroutine TRIANGLE_GRID_EDGE

End Module Ocean_Model 
