!=======================================================================
! Fiscm FVCOM Driver 
!
! Description
!    - Read and setup mesh
!    - Store calculated grid metrics for future use
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
!  18/11/2016: J. Ounsley
!    - Including grid metrics cache subroutines to avoid recalculation
!      on successive runs. Adapted from Pierre Cazenave's FVCOM 
!      implementation. The  file name and whether to force creation
!      of the metrics file is configurable in the NML_NCVAR namelist
!    - Fixed bug in conversion of local cartesian velocities to 
!      spherical coordinates, and updated find_element methods to handle
!      0/360 degree boundary
!  1/4/2017: J. Ounsley
!    - Added method for horizontal diffusion with repeating
!      attempts upon collision with boundary
!  7/4/2017: J. Ounsley
!    - Added method for calculating the gradient of node based 
!      field within a cell
!  30/5/2017: J. Ounsley
!    - Fixed bug with failure to properly initialise z coordinates
!=======================================================================
Module Ocean_Model 
use gparms
use mod_igroup
use forcing
use utilities
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
  integer, allocatable :: isonb(:)            !!NODE MARKER = 0,1,2  JO - not used, but useful output
  integer, allocatable :: isbce(:)            !! JO - not used, but useful output
  integer, allocatable :: nbvt(:,:)           !! JO - not used, but useful output
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
  a2u_char,art_char,nele_char,node_char,zeta_char,omega_char,         &
  siglay_char,siglev_char,wu_char,wv_char  

  ! JO - Allow configuration of existing metrics file
  character(len=fstr) :: metrics_file
  logical :: force_metrics
  
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
         art_char,   &
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
     viscofm_char,   &
       ! JO - Allow configuration of metrics file
       !      This file stores calculated grid metrics to save time on
       !      future runs. Is alsoo useful for post. proc.
       !      TODO consider wether this should be on this namelist
       metrics_file, & ! The name of the metrics file to be loaded/saved
       force_metrics   ! Flag to indicate a forced load of the metrics file

logical :: grid_metrics

interface interp
  module procedure interp_float2D
  module procedure interp_float3D
end interface interp

! JO - This is no longer in use, interp is called
!      instead
interface interp_from_nodes
  module procedure interp_flt_from_nodes 
end interface interp_from_nodes

interface gradient
   module procedure gradient_float2D
   module procedure gradient_float3D
end interface gradient

contains

!----------------------------------------------------
! Read the mesh and interpolation coefficients 
! Add mesh to output files for viz
! If exist, load pre calculated metrics
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
  logical :: FEXIST
  ! JO - now configurable throught the namelist 
  !character(len=fstr) :: METIN, METOUT

  ! JO - Initialise the metric file name and whether
  !      or not to overwrite current metrics file
  metrics_file = "./metrics.nc"
  force_metrics = .false.

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
      ! JO - kh_char only for active vdiff (not present in SSM model)
      if (g(n)%vdiff_type > 0) then
        lsize = lsize + 1 ; varlist(lsize) = kh_char
      endif
      if (g(n)%hdiff_type ==HDIFF_VARIABLE)then
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
  allocate(siglay(N_verts,N_lay))
  allocate(siglev(N_verts,N_lev))
  allocate(esiglay(N_elems,N_lay))
  allocate(esiglev(N_elems,N_lev))

  allocate(a2u(N_elems,4)) ;a2u   = zero
  allocate(a1u(N_elems,4)) ;a1u   = zero

  !----------------Node, Boundary Condition, and Control Volume-----------------------!
  allocate(nbe(N_elems,3))          ;nbe      = 0  !! Indices of element neighbours
  ALLOCATE(NTVE(0:N_verts))           ;NTVE     = 0
  ALLOCATE(ISONB(N_verts))          ;ISONB    = 0  !!NODE MARKER = 0,1,2
  ALLOCATE(ISBCE(N_elems))          ;ISBCE    = 0

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
  msg = "error reading siglay"
  call ncdchk( nf90_inq_varid(fid,siglay_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, siglay),msg)
  msg = "error reading siglev"
  call ncdchk( nf90_inq_varid(fid,siglev_char,varid),msg )
  call ncdchk(nf90_get_var(fid, varid, siglev),msg)

  !read secondary connectivity (nbve/ntve) 
  ! JO Assume that if maxelem is not in the 
  ! file then the following also are not
  ! nbe, ntve, nbve, ntsn, nbsn
  ! awx, awy, aw0, a1u, a2u
  ! NOTE isonb and isbce are useful metrics for
  ! plotting boundaries that are not stored in .nc files
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
     msg = "error reading nbe"
     call ncdchk( nf90_inq_varid(fid,nbe_char,varid),msg )
     call ncdchk(nf90_get_var(fid, varid, nbe),msg)
   
     msg = "error reading aw0"
     call ncdchk(nf90_inq_varid(fid,aw0_char,varid),msg )
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
     grid_metrics = .true.
    else 
       write(*,*)'WARNING:::::: grid metrics Do NOT exist in forcing file'
       write(*,*)'Will try to compute internally' 
       write(*,*)'Proceeding with: 1st Order interpolation' 
       write(*,*)'AW0 = 1/3; AWX=AWY = 0'
       write(*,*)'In the future, select [grid metrics] in your NetCDF namelist'
    endif

  !calculate cell center coordinates
  ! JO - The following translates local coordinates
  ! away from 360 / 0 degree boundary to avoid issues with 
  ! interpolation across this boundary, xc0 and yc0 are used in
  ! interp function
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
    
    ! Revised by J. Ounsley 18/11/16
    ! If we did not find grid metrics or we are forcing overwrite
    ! then calculate and save, alternatvively load if file already exists
    if(.not.grid_metrics .or. force_metrics)then
       INQUIRE(FILE=metrics_file, EXIST=FEXIST)
       IF (FEXIST .and. .not.force_metrics) THEN
          !--READ SAVED GRID METRICS FROM NETCDF FILE
          WRITE (*,*) 'Attempting to load saved grid metrics'
          WRITE (*,*) 'To recalculate delete metrics.nc' !or
          !WRITE (*,*) 'set force_metrics to true in NML_NCVAR namelist'
          CALL NCD_READ_METRICS(metrics_file,FEXIST)
          ! Correct the positions using offsets VXMIN and VYMIN from netCDF metrics file !
          ! This process would otherwise have taken place in TRIANGLE_GRID_EDGE so we need
          ! to replicate that behaviour here					       !
          !--------------SHIFT GRID TO UPPER RIGHT CARTESIAN-----------------------------!
          !VX = VX - VXMIN
          !VY = VY - VYMIN
          !--------------CALCULATE GLOBAL ELEMENT CENTER GRID COORDINATES----------------!
          !DO I=1,N   
          !   XC(I)  = (VX(NV(I,1)) + VX(NV(I,2)) + VX(NV(I,3)))/3.0_SP
          !   YC(I)  = (VY(NV(I,1)) + VY(NV(I,2)) + VY(NV(I,3)))/3.0_SP
          !END DO
          !XC(0) = 0.0_SP ; YC(0) = 0.0_SP    
       END IF
       ! Check that read metrics succeeded
       IF (FEXIST) THEN
          WRITE (*,*) 'Grid metrics loaded successfully'          
       ELSE 
          IF(force_metrics)THEN
             WRITE (*,*) 'Forced generation of grid metrics.'          
          ELSE
             WRITE (*,*) 'Grid metrics cannot be loaded.'          
          END IF
          WRITE (*,*) 'Calling triangle_grid_edge...'
          CALL TRIANGLE_GRID_EDGE
          WRITE(*,*)  'Finished'
          !--SAVE THE GRID METRICS TO NETCDF FILE
          WRITE (*,*) 'Saving grid metrics'
          CALL NCD_WRITE_METRICS(metrics_file)
          WRITE (*,*) 'Grid metrics saved'
          WRITE (*,*)
       END IF
    end if ! grid_metrics

    ! At this point grid metrics were either in forcing file,
    ! loaded from metrics file or calculated from scratch
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
  integer  :: i,p,np
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
       do p=1,np
         if(istatus(p)==ACTIVE)then
           pdxt(p) = x(p) + normal()*tscale
           pdyt(p) = y(p) + normal()*tscale
         else
           pdxt(p) = x(p)
           pdyt(p) = y(p)
         endif
       end do
     elseif (spherical == 1)then
       do p=1,np
         if(istatus(p)==ACTIVE)then
             ! JO Fixed bug here where degree to radians conversion was missing
             pdyt(p) = y(p)  + normal()*tscale/(radius_earth*d2r)
             pdxt(p) = x(p)  + normal()*tscale/(radius_earth*d2r*COS(pdyt(p)*d2r) + 1.0E-6)
          else
           pdxt(p) = x(p)
           pdyt(p) = y(p)
         endif
       end do

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
     endif !spherical


     

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
! Author JO
!
! Adapted version of the horizontal random-walk
! that moves to current following advection upon
! collision
!
! Random-Walk horizontal diffusion with constant
!    turbulent eddy diffusivity
!
! m. huret use a Gauss dist. , this probably wont
! converge to the correct continuous form of diffusion
! here we will use a uniform random walk
!----------------------------------------------------
subroutine rw_hdiff_constant_rheotaxis(g, dT)
  use utilities, only : normal,unitrand
    type(igroup), intent(inout) :: g
    real(sp), intent(in) :: dT

    real(sp), pointer :: x(:)
    real(sp), pointer :: y(:)
    real(sp), pointer :: u(:)
    real(sp), pointer :: v(:)
    integer,  pointer :: istatus(:)
    integer,  pointer :: cell(:)
    integer  :: i,np,p,attempt
    real(sp) :: tscale, dist
    real(sp), allocatable :: pdxt(:), pdyt(:)

    !set pointers to x,y particle positions
    call get_state('x',g,x)
    call get_state('y',g,y)
    call get_state('u',g,u)
    call get_state('v',g,v)
    call get_state('status',g,istatus)
    call get_state('cell',g,cell)

    !set dimensions for loops and time step
    np = g%nind

    allocate(pdxt(np))  ;
    allocate(pdyt(np))  ;
    !set diffusive time scale
    tscale = sqrt(2.*dT*g%hdiff_const_val)

    !horizontal random walk
    !$OMP PARALLEL DO SCHEDULE(STATIC) IF(multithread) 
    do p=1,np
       ! Is this particle currently active?
       if(istatus(p)==ACTIVE)then

         

           if(spherical==0)then
                pdxt(p) = x(p) + normal()*tscale
                pdyt(p) = y(p) + normal()*tscale

           elseif (spherical == 1)then
                pdyt(p) = y(p)  + normal()*tscale/(radius_earth*d2r)
                pdxt(p) = x(p)  + normal()*tscale/(radius_earth*d2r*COS(pdyt(p)*d2r) + 1.0E-6)

                if( pdxt(p) < 0.0_SP)  pdxt(p) = pdxt(p) + 360.0_SP
                if( pdxt(p) > 360.0_SP) pdxt(p) = pdxt(p) - 360.0_SP

                if( pdyt(p) > 90.0_SP)  pdyt(p) = 180.0_SP - pdyt(p)
                if( pdyt(p) < -90.0_SP) pdyt(p) =  - 180.0_SP - pdyt(p)

           endif !if spherical

             ! Use some of the elements from the find_element 
             ! routine to find this particle specifically
             ! If the particle is not in a valid cell, keep looking

       cell(p) = find_element_lazy(pdxt(p),pdyt(p),cell(p))
       if(cell(p) == 0) then ! Did not find valid cell, keep looking
           cell(p) = find_element_robust(pdxt(p),pdyt(p)) 
       end if        
         
       ! Did we fail to find a valid cell?
       ! If so, move with current
       if(cell(p) == 0) then
	   !print *,'Collision on diffusion' 
	   dist = normal()
	   if(spherical==0)then
               pdxt(p) = x(p) + dist*tscale*u(p) / (sqrt(u(p)**2+v(p)**2) + 1.0E-6)
               pdyt(p) = y(p) + dist*tscale*v(p) / (sqrt(u(p)**2+v(p)**2) + 1.0E-6)
           elseif(spherical==1)then
               pdyt(p) = y(p) + dist*tscale*v(p) / (radius_earth*d2r*sqrt(u(p)**2+v(p)**2) + 1.0E-6)
               pdxt(p) = x(p) + dist*tscale*u(p) / (radius_earth*d2r*COS(pdyt(p)*d2r)*sqrt(u(p)**2+v(p)**2) + 1.0E-6)


               if( pdxt(p) < 0.0_SP)  pdxt(p) = pdxt(p) + 360.0_SP
               if( pdxt(p) > 360.0_SP) pdxt(p) = pdxt(p) - 360.0_SP

               if( pdyt(p) > 90.0_SP)  pdyt(p) = 180.0_SP - pdyt(p)
               if( pdyt(p) < -90.0_SP) pdyt(p) =  - 180.0_SP - pdyt(p)
           end if
       end if

       cell(p) = find_element_lazy(pdxt(p),pdyt(p),cell(p)) 
       if(cell(p) == 0) then ! Did not find valid cell, keep looking
           cell(p) = find_element_robust(pdxt(p),pdyt(p)) 
       end if        
       
       ! Did we fail to find a valid cell?
       ! If so mark as exited
       if(cell(p) == 0) then 
            istatus(p) = EXITED       
            !print *,'Failed to avoid collision after rheotaxis'
       end if

       else ! Not active, so don't change position
          pdxt(p) = x(p)
          pdyt(p) = y(p)
       end if ! if active 

    end do ! do p=1,np 
    !$OMP END PARALLEL DO


    ! Update positions of non exited particles
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

end subroutine rw_hdiff_constant_rheotaxis

  !----------------------------------------------------
  ! Author JO
  !
  ! Adapted version of the horizontal random-walk
  ! that repeats random sampling upon loosing the
  ! particle, in an attempt to avoid getting stuck
  ! upon collisions
  !
  ! Random-Walk horizontal diffusion with constant
  !    turbulent eddy diffusivity
  !
  ! m. huret use a Gauss dist. , this probably wont
  ! converge to the correct continuous form of diffusion
  ! here we will use a uniform random walk
  !----------------------------------------------------
  subroutine rw_hdiff_constant_collision(g, dT)
    use utilities, only : normal,unitrand
    type(igroup), intent(inout) :: g
    real(sp), intent(in) :: dT

    real(sp), pointer :: x(:)
    real(sp), pointer :: y(:)
    integer,  pointer :: istatus(:)
    integer,  pointer :: cell(:)
    integer  :: i,np,p,attempt
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
    !$OMP PARALLEL DO SCHEDULE(STATIC) IF(multithread) 
    do p=1,np
       ! Is this particle currently active?
       if(istatus(p)==ACTIVE)then

          do attempt=1,g%hdiff_collision_repeats ! try to sample valid point multiple times

             if(spherical==0)then
                pdxt(p) = x(p) + normal()*tscale
                pdyt(p) = y(p) + normal()*tscale

             elseif (spherical == 1)then
                pdyt(p) = y(p)  + normal()*tscale/(radius_earth*d2r)
                pdxt(p) = x(p)  + normal()*tscale/(radius_earth*d2r*COS(pdyt(p)*d2r) + 1.0E-6)

                if( pdxt(p) < 0.0_SP)  pdxt(p) = pdxt(p) + 360.0_SP
                if( pdxt(p) > 360.0_SP) pdxt(p) = pdxt(p) - 360.0_SP

                if( pdyt(p) > 90.0_SP)  pdyt(p) = 180.0_SP - pdyt(p)
                if( pdyt(p) < -90.0_SP) pdyt(p) =  - 180.0_SP - pdyt(p)

             endif !if spherical

             ! Use some of the elements from the find_element 
             ! routine to find this particle specifically
             ! If the particle is not in a valid cell, keep looking

             cell(p) = find_element_lazy(pdxt(p),pdyt(p),cell(p))
             if(cell(p) /= 0) exit ! Found valid cell, so stop looking

             cell(p) = find_element_robust(pdxt(p),pdyt(p)) 
             if(cell(p) /= 0) exit ! Found valid cell, so stop looking             
        
          end do

          ! Did we fail to find a valid cell?
          ! If so mark as exited
          if(cell(p) == 0) then 
             istatus(p) = EXITED       
             !print *,'Failed to avoid collision after ', g%hdiff_collision_repeats, ' repeats of hdiff'
          end if

       else ! Not active, so don't change position
          pdxt(p) = x(p)
          pdyt(p) = y(p)
       endif ! if active 

    end do ! do p=1,np 
    !$OMP END PARALLEL DO


    ! Update positions of non exited particles
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

  end subroutine rw_hdiff_constant_collision


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
  integer :: counter

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

    !JO - Fixed bug, need to initialise pdxt and pdyt to initial
    !     positions, such that inactive particles are not changed
    pdxt = x
    pdyt = y

  !horizontal random walk
     if(spherical == 0 )then
      where(istatus == ACTIVE)
       pdxt = x + normal()*tscale
       pdyt = y + normal()*tscale
      end where
     elseif (spherical == 1)then
      where(istatus == ACTIVE)
          ! JO Fixed bug here where degree to radians conversion
          pdyt = y  + normal()*tscale/(radius_earth*d2r)
          pdxt = x  + normal()*tscale/(radius_earth*d2r*COS(pdyt*d2r) + 1.0E-6)
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
  ! JO Want all particles to have h and zeta
  ! interpolated, so set a temp status to 1
  ! for this purpose
  !integer , pointer :: istatus(:)
  integer           :: temp_status(np)
  temp_status = 1

  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('s',g,s)
  call get_state('z',g,z)
  call get_state('h',g,h)
  call get_state('zeta',g,zeta)
  call get_state('cell',g,cell)


  call interp(np,x,y,cell,temp_status,zeta_char,zeta,3)
  call interp(np,x,y,cell,temp_status,h_char,h,3)
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
  real(sp), pointer :: u_particle(:) ! JO - Used to track particle velocity (output)
  real(sp), pointer :: v_particle(:) ! JO - Used to track particle velocity (output)
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
    ! JO - This subroutine is the most likely candidate for
    !      parallelisation as this, and called subroutines, is
    !      where the code spends most of its time. 
    !      It should also be possible to parallelise advection on each
    !      particle independantly

  diel=time/3600.0-int(time/3600.0/24.0)*24.0
  !set dimensions for loops and time step
  !np = g%nind
  !set pointers to states
  call get_state('x',g,x)
  call get_state('y',g,y)
  call get_state('s',g,s)
  call get_state('z',g,z)  
  call get_state('h',g,h)
  call get_state('u',g,u_particle) 
  call get_state('v',g,v_particle) 
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
     
     ! PARALLEL candidate
     if (spherical == 0)then 
     pdx(:) = x(:)  + a_rk(ns)*deltaT*chix(:,ns-1)
     pdy(:) = y(:)  + a_rk(ns)*deltaT*chiy(:,ns-1)
     pdz(:) = s(:)  + a_rk(ns)*deltaT*chiz(:,ns-1)
!!!!!
     elseif (spherical == 1)then 
     ! JO - Corrected bug in transformation here, degree to radian transform
     pdy(:) = y(:)  + a_rk(ns)*deltaT*chiy(:,ns-1)/(radius_earth*d2r)
     pdx(:) = x(:)  + a_rk(ns)*deltaT*chix(:,ns-1)/(radius_earth*d2r*COS(pdy(:)*d2r) + 1.0E-6) 
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
          ! PARALLEL candidate
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

end do !--Loop over RK Stages
  !--Sum Stage Contributions to get Updated Particle Positions-------------------!
  pdxt(:)  = x(:)
  pdyt(:)  = y(:)
  pdzt(:)  = s(:)
  u_particle(:) = 0
  v_particle(:) = 0
  do ns=1,mstage
     if (spherical == 0)then 
          ! JO - Removing multiplication with istatus
          ! as this is handled below
     pdxt(:) = pdxt(:) + deltaT*chix(:,ns)*b_rk(ns)!*FLOAT(istatus(:))
     pdyt(:) = pdyt(:) + deltaT*chiy(:,ns)*b_rk(ns)!*FLOAT(istatus(:))
     pdzt(:) = pdzt(:) + deltaT*chiz(:,ns)*b_rk(ns)!*FLOAT(istatus(:))    



    ! JO - Store the interpolated particle velocity
    u_particle(:) = u_particle(:) + chix(:,ns)*b_rk(ns)!*FLOAT(istatus(:))
    v_particle(:) = v_particle(:) + chiy(:,ns)*b_rk(ns)!*FLOAT(istatus(:))

!!!!!
     elseif (spherical == 1)then 
     ! JO - Corrected bug in transformation here with degree to radian transform
     !      Additionally, incorrect stage indices was in use, possible copy paste
     !      error
     pdyt(:) = pdyt(:)  + b_rk(ns)*deltaT*chiy(:,ns)/(radius_earth*d2r)
     pdxt(:) = pdxt(:)  + b_rk(ns)*deltaT*chix(:,ns)/(radius_earth*d2r*COS(pdyt(:)*d2r) + 1.0E-6)
     pdzt(:) = pdzt(:)  + b_rk(ns)*deltaT*chiz(:,ns)
          
     ! JO - Store the interpolated particle velocity
     u_particle(:) = u_particle(:) + b_rk(ns)*chix(:,ns)/(radius_earth*d2r*COS(pdyt(:)*d2r) + 1.0E-6)
     v_particle(:) = v_particle(:) + b_rk(ns)*chiy(:,ns)/(radius_earth*d2r)

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
  end do !--Sum Stage Contributions to get Updated Particle Positions-------------------!
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
  nullify(u_particle)
  nullify(v_particle)
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
             ! JO - The following translates local coordinates
             ! away from 360 / 0 degree boundary to avoid issues with 
             ! interpolation across this boundary, xc0 and yc0 are similarly
             ! translated
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

       ! JO
       ! PARALLEL candidate - Test with 100 & 5000
       ! particles show no gain in efficiency
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
       ! JO - The following translates local coordinates
       ! away from 360 / 0 degree boundary to avoid issues with 
       ! interpolation across this boundary, xc0 and yc0 are similarly
       ! translated
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

  !------------------------------------------------
  ! JO
  ! Calculate the gradient vector of a field
  ! through taking the normal of the plane described
  ! by the field at the current cell
  ! Sets dx,dy to the horizontal components of the
  ! gradient with unit magnitude in horizontal, pointing up the slope.
  ! dv is the change in the variable over unit distance
  !
  ! TODO: 
  ! 1) This does not take into account any interpolation
  ! coefficients
  ! 2) Currently only implemented for node based fields
  !-------------------------------------------------
  subroutine gradient_float2D(np,x,y,cell,istatus,vname,dx,dy,dv,iframe)
    integer, intent(in)    :: np,iframe
    real(sp),intent(in)    :: x(np)
    real(sp),intent(in)    :: y(np)
    integer, intent(in)    :: cell(np)
    integer, intent(in)    :: istatus(np)
    character(len=*)       :: vname
    real(sp),intent(inout) :: dx(np),dy(np),dv(np)
    !-------------------------------
    real(sp), pointer :: field(:)  !pointer to FVCOM field 
    integer :: i,icell,d1
    integer  :: n1,n2,n3
    real(sp) :: p1(3), p2(3), p3(3), normal(3)
 
    !point to forcing variable
    call get_forcing(trim(vname),field,iframe) 

    !initialize pu,pv
    dx = 0.0
    dy = 0.0
    dv = 0.0

    !determine the dimensions
    d1 = size(field,1)

    if(d1 == N_verts)then !vertex-based 2-D array
       do i=1,np
          icell = cell(i)
          if(istatus(i) < 1 .or. icell == 0)cycle 
          ! These nodes should be in a clockwise orientation
          n1  = tri(icell,1)
          n2  = tri(icell,2)
          n3  = tri(icell,3)
          
          ! The positions of the nodes, and their value
          p1 = (/ xm(n1),ym(n1),field(n1) /)
          p2 = (/ xm(n2),ym(n2),field(n2) /)
          p3 = (/ xm(n3),ym(n3),field(n3) /)
          
          ! Returns the normal vector, 
          ! right hand rule, 
          ! p2-p1 X p3-p1 into page.
          ! p2
          ! ^
          ! |
          ! p1-->p3
          ! The x and y components will point up
          ! the gradient        
          call cross_product_3(p2-p1,p3-p1,normal)
          ! Scale the normal vector to have magnitude
          ! of 1 in the horizontal. The 3rd component
          ! is then the reciprocal of size of the vertical gradient,
          ! or 1 / (change in variable over unit distance)
          ! we need to take the reciprocal for the actual
          ! magnitude of the gradient
          normal = normal/sqrt(normal(1)**2+normal(2)**2)
          ! Assign horizontal direction
          dx(i) = normal(1)
          dy(i) = normal(2)
          dv(i) = 1.0/normal(3) ! Magnitude of gradient
       end do
    else
       write(*,*)'field has horizontal dimensions that is not nodes'
       write(*,*)'do not know how to interpolate'
       stop
    endif

  end subroutine gradient_float2D! gradient
 



  !------------------------------------------------
  ! JO
  ! Calculate the gradient vector of a field
  ! through taking the normal of the plane described
  ! by the field at the current cell
  ! Sets dx,dy to the horizontal components of the
  ! gradient with unit magnitude in horizontal, pointing up the slope.
  ! dv is the change in the variable over unit distance
  !
  ! TODO: 
  ! 1) This does not take into account any interpolation
  ! coefficients aw0, awx, awy
  ! 2) Currently only implemented for node based fields
  !-------------------------------------------------
  subroutine gradient_float3D(np,x,y,s,cell,istatus,vname,dx,dy,dv,iframe)
    integer, intent(in)    :: np,iframe
    real(sp),intent(in)    :: x(np)
    real(sp),intent(in)    :: y(np)
    real(sp),intent(in)    :: s(np)
    integer, intent(in)    :: cell(np)
    integer, intent(in)    :: istatus(np)
    character(len=*)       :: vname
    real(sp),intent(inout) :: dx(np),dy(np),dv(np)
    !-------------------------------
    real(sp), pointer :: field(:,:)  !pointer to FVCOM field 
    real(sp) :: f1,f2, dx1, dx2, dy1, dy2, dv1, dv2
    integer :: i,icell,k1,k2,ktype,d1,d2
    integer  :: n1,n2,n3,k
    real(sp) :: p1(3), p2(3), p3(3), normal(3)
 
    !point to forcing variable
    k = 1
    call get_forcing(trim(vname),field,iframe) 

    !initialize pu,pv
    dx = 0.0
    dy = 0.0
    dv = 0.0

    !determine the dimensions
    d1 = size(field,1)
    d2 = size(field,2)

   !set vertical type (layers or levels)
    if(d2 == N_lev)then
       ktype = level_based 
    else if(d2 == N_lay)then
       ktype = layer_based 
    else
       write(*,*)'error in gradient_float3D'
       write(*,*)'second dimension of field is not layers or levels'
       write(*,*)'==: ',d2
       stop
    endif


    if(d1 == N_verts)then !vertex-based 2-D array
       do i=1,np
          icell = cell(i)
          if(istatus(i) < 1 .or. icell == 0)cycle 

          ! Determine the vertical layer brackets and interp coeffs
          call get_vert_interpcoefs(icell,s(i),k1,k2,f1,f2,ktype)

          !print *, i,s(i),k1,k2,f1,f2

          ! These nodes should be in a clockwise orientation
          n1  = tri(icell,1)
          n2  = tri(icell,2)
          n3  = tri(icell,3)
          
          ! Layer 1
          k = k1
          ! The positions of the nodes, and the value of
          ! the field at that level/layer
          p1 = (/ xm(n1),ym(n1),field(n1,k) /)
          p2 = (/ xm(n2),ym(n2),field(n2,k) /)
          p3 = (/ xm(n3),ym(n3),field(n3,k) /)
          
          ! Returns the normal vector, 
          ! right hand rule, 
          ! p2-p1 X p3-p1 into page.
          ! p2
          ! ^
          ! |
          ! p1-->p3
          ! The x and y components will point up
          ! the gradient        
          call cross_product_3(p2-p1,p3-p1,normal)
          ! Scale the normal vector to have magnitude
          ! of 1 in the horizontal. The 3rd component
          ! is then the inverse of the vertical gradient,
          ! or 1 / (change in variable over unit distance)
          ! we need to take the reciprocal for the actual
          ! magnitude of the gradient
          if(normal(1) /= 0 .or. normal(2) /= 0)  normal = normal/sqrt(normal(1)**2.+normal(2)**2.)
          ! Assign values
          dx1 = normal(1)
          dy1 = normal(2)
          dv1 = 1.0/normal(3)
        
          ! Layer 2
          k = k2
          ! The positions of the nodes, and the value of
          ! the field at that level/layer
          p1 = (/ xm(n1),ym(n1),field(n1,k) /)
          p2 = (/ xm(n2),ym(n2),field(n2,k) /)
          p3 = (/ xm(n3),ym(n3),field(n3,k) /)
          
          ! Returns the normal vector, 
          ! right hand rule, 
          ! p2-p1 X p3-p1 into page.
          ! p2
          ! ^
          ! |
          ! p1-->p3
          ! The x and y components will point up
          ! the gradient        
          call cross_product_3(p2-p1,p3-p1,normal)
          ! Scale the normal vector to have magnitude
          ! of 1 in the horizontal. The 3rd component
          ! is then the inverse of the vertical gradient,
          ! or 1 / (change in variable over unit distance)
          ! we need to take the reciprocal for the actual
          ! magnitude of the gradient
          if(normal(1) /= 0 .or. normal(2) /= 0)  normal = normal/sqrt(normal(1)**2.+normal(2)**2.)
          ! Assign values
          dx2 = normal(1)
          dy2 = normal(2)
          dv2 = 1.0/normal(3)
             
          ! Interpolate between layers
          dx(i) = f1*dx1 + f2*dx2
          dy(i) = f1*dy1 + f2*dy2
          dv(i) = f1*dv1 + f2*dv2

       end do
    else
       write(*,*)'field has horizontal dimensions that is not nodes'
       write(*,*)'do not know how to interpolate'
       stop
    endif

  end subroutine gradient_float3D! gradient
 

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

  ! JO
  ! The program spends most of its time here
  ! particularly with large advection per timestep
  ! which induces robust searches.
  ! Most obvious candidate for parallelisation    
  !$OMP PARALLEL DO SCHEDULE(STATIC) IF(multithread)
  do p=1,np
    if(istatus(p) < 0)cycle

    !try a quick find
    cell(p) = find_element_lazy(x(p),y(p),cell(p))
    if(cell(p) /= 0)cycle

    !failed, try a robust find
    !print *,'Failed lazy element find for p:',p
      cell(p) = find_element_robust(x(p),y(p)) 
      if(cell(p) /= 0)cycle
       !print *,'Failed robust element find for p:',p

    !failed, update status to lost
    istatus(p) = EXITED
  end do
  !$OMP END PARALLEL DO

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
  integer  :: elem,i,j,k,test,iney
  !---------------------------
  integer  :: pts(3)
  real(sp) :: xv(3),yv(3),xp_test

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

    ! JO - Use the correct vertex, spherical
    !      are translated away from 0/360
    if(spherical==0)then
       xv = xm(pts)
       yv = ym(pts)
    else if(spherical==1)then
       xv = xm0(pts)
       yv = ym0(pts)
    end if

    ! JO - Need to take into account elements that
    !      cross the 0 degree boundary
    if(spherical==0)then
       xp_test = xp
    elseif(spherical==1)then
       if(xp >= 0.0_sp .and. xp <=180.0_sp)then
          xp_test = xp + 180.0_sp
       elseif( xp > 180.0_sp .and. xp <=360.0_sp)  then
          xp_test = xp - 180.0_sp
       endif
    endif


  if(isintriangle(last,xp_test,yp,xv,yv))then     
    elem = last
  else
    if(grid_metrics)then
    outer: do j=1,3
      test = tri(last,j)
           do k=1,ntve(test)
              iney = nbve(test,k)
              pts(1:3) = tri(iney,1:3)

                ! JO - Use the correct vertex, spherical
                !      are translated away from 0/360
                if(spherical==0)then
                   xv = xm(pts)
                   yv = ym(pts)
                else if(spherical==1)then
                   xv = xm0(pts)
                   yv = ym0(pts)
                end if

                ! JO - Need to take into account elements that
                !      cross the 0 degree boundary 
                if(spherical==0)then
                   xp_test = xp
                elseif(spherical==1)then
                   if(xp >= 0.0_sp .and. xp <=180.0_sp)then
                      xp_test = xp + 180.0_sp
                   elseif( xp > 180.0_sp .and. xp <=360.0_sp)  then
                      xp_test = xp - 180.0_sp
                   endif
                endif

              if(isintriangle(iney,xp_test,yp,xv,yv))then
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

  ! JO - A temp variable for spherical adjustment
  real(sp) :: xp_test

  !make sure we have a mesh  
  if(.not.mesh_setup)then
    write(*,*)'error in find_element_robust'
    write(*,*)'mesh is not setup yet'
    stop
  endif
  
  !initialize
  elem = 0
  cnt = 0

    ! JO - Need to take into account elements that
    !      cross the 0 degree boundary       
    if(spherical==0)then
       xp_test = xp
    else if(spherical==1)then       
       if(xp >= 0.0_sp .and. xp <=180.0_sp)then
          xp_test = xp + 180.0_sp
       elseif( xp > 180.0_sp .and. xp <=360.0_sp)  then
          xp_test = xp - 180.0_sp
       endif
    endif


    if(spherical==0)then
       radlist(1:N_elems,1) = sqrt((xc(1:N_elems)-xp_test)**2 + (yc(1:N_elems)-yp)**2)
    else if(spherical==1) then
       ! JO - Use the spherical centres (xc0,yc0), which are translated
       !      away from the 0/360 boundary
       radlist(1:N_elems,1) = sqrt((xc0(1:N_elems)-xp_test)**2 + (yc0(1:N_elems)-yp)**2)
    end if

    !print *,xp,yp  
  radlast = -1.0_sp
  l1:  do while(cnt < max_check)
    cnt = cnt+1
    locij   = minloc(radlist,radlist>radlast)
    min_loc = locij(1)
    if(min_loc == 0) then
        print *,'min_loc=0'
      exit l1
    end if

    pts(1:3) = tri(min_loc,1:3)

       ! JO - Use the correct vertices, spherical (xm0,ym0)
       !      are translated away from 0/360
       if(spherical==0)then
          xv = xm(pts)
          yv = ym(pts)
       else if(spherical==1)then
          xv = xm0(pts)
          yv = ym0(pts)
       end if

       !print *,'pts',pts
       !print *,'xv',xv
       !print *,'yv',yv
    if(isintriangle(min_loc,xp_test,yp,xv,yv))then
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

  !----------------------------------------------------------------------
  ! Revised by J. Ounsley 18/11/16
  ! Adapted to allow for saving and loading of derived metrics
  ! for successive use.
  !----------------------------------------------------------------------
 subroutine   TRIANGLE_GRID_EDGE
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: TEMP,NB_TMP,CELLS,NBET
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: CELLCNT
  INTEGER                              :: I,J,II,JJ,NTMP,NCNT,NFLAG,JJB
  INTEGER                              :: N1,N2,N3,J1,J2,J3

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
  Max_Elems = 0
  DO I=1,N_verts
     NCNT = 0
     DO J=1,N_elems
        IF( FLOAT(tri(J,1)-I)*FLOAT(tri(J,2)-I)*FLOAT(tri(J,3)-I) == 0.0_SP) &
             NCNT = NCNT + 1
     END DO
     Max_Elems = MAX(Max_Elems,NCNT)
  END DO

  !
  !----ALLOCATE ARRAYS BASED ON MX_NBR_ELEM--------------------------------------!
  !
  ALLOCATE(NBVE(N_verts,Max_Elems+1))
  ALLOCATE(NBVT(N_verts,Max_Elems+1))
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

  ALLOCATE(NB_TMP(N_verts,Max_Elems+1))
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

  ! New subroutine to create a metrics netCDF file (to be read in on
  ! successive launches of the offline Lagrangian code).
  ! Pierre Cazenave, Plymouth Marine Laboratory
  !
  ! Modified by James Ounsley 18/11/2016
  ! Adapted for use in FSICM
  SUBROUTINE NCD_WRITE_METRICS(INFILE)
    !--------------------------------------------------------------------------!
    ! WRITE METRICS CALCULATED BY TRIANGLE_GRID_EDGE TO NETCDF FOR LATER REUSE.
    !
    ! INPUTS:
    ! INFILE - netCDF output file string [max=100] [CHARACTER]
    !
    ! Copied from FVCOM offlag Pierre Cazenave (07/01/2015) Plymouth 
    ! Marine Laboratory.
    ! Adapted by J. Ounsley 18/11/2016 for FISCM
    !--------------------------------------------------------------------------!

    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    CHARACTER(LEN=100), INTENT(IN)        :: INFILE
    !--------------------------------------------------------------------------!
    CHARACTER(LEN=100)                  :: TSTRING,NETCDF_CONVENTION
    INTEGER                             :: IERR
    INTEGER                             :: THREE_DID,MAXELEM_DID,NELE_DID,&
         NODE_DID
    INTEGER                             :: NBE_VID,NTVE_VID,NBVE_VID,NBVT_VID,&
         ISONB_VID,ISBCE_VID
    INTEGER                             :: NC_FID
    ! Add the coordinate offsets (VXMIN and VYMIN). Thanks to Judith Wolf
    ! (National Oceanography Centre, Liverpool, UK) for this bug report.
    INTEGER                             :: VXMIN_VID, VYMIN_VID, MX_NBR_ELEM_VID
    INTEGER                             :: STAT1DM(1),STAT1DN(1)
    INTEGER                             :: STAT2D3(2),STAT2DM(2)
    INTEGER, ALLOCATABLE                :: START(:), CNT(:)
    !--------------------------------------------------------------------------!

    IERR = NF90_CREATE(TRIM(INFILE),NF90_CLOBBER,NC_FID)
    IF(IERR /= NF90_NOERR) THEN
       WRITE(*,*)'Error creating',TRIM(INFILE)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    ! netCDF Convention string
    NETCDF_CONVENTION = 'CF-1.0'

    ! Global Attributes
    TSTRING = "FVCOM Offline Lagrangian Grid Metrics File"
    IERR = NF90_PUT_ATT(NC_FID,NF90_GLOBAL,"title",TSTRING)
    IERR = NF90_PUT_ATT(NC_FID,NF90_GLOBAL,"institution","SMAST, PML")
    IERR = NF90_PUT_ATT(NC_FID,NF90_GLOBAL,"source","OFFLINE_FVCOM")
    IERR = NF90_PUT_ATT(NC_FID,NF90_GLOBAL,"Conventions",TRIM(NETCDF_CONVENTION))

    ! Dimensioning
    IERR = NF90_DEF_DIM(NC_FID,"node",N_verts,NODE_DID)
    IERR = NF90_DEF_DIM(NC_FID,"nele",N_elems,NELE_DID)
    IERR = NF90_DEF_DIM(NC_FID,"maxelem",Max_Elems+1,MAXELEM_DID)
    IERR = NF90_DEF_DIM(NC_FID,"three",3,THREE_DID)

    STAT2DM      = (/NODE_DID,MAXELEM_DID/)     !!Static 2d var [nodes,9]
    STAT2D3      = (/NELE_DID,THREE_DID/)       !!Static 2d var [elements,3]
    STAT1DM      = (/NODE_DID/)                 !!Static 1d var [nodes]
    STAT1DN      = (/NELE_DID/)                 !!Static 1d var [elements]

    ! Variable Definitions

    !!===== NBE ===============================================================!
    IERR = NF90_DEF_VAR(NC_FID,"nbe",NF90_INT,STAT2D3,NBE_VID)
    IERR = NF90_PUT_ATT(NC_FID,NBE_VID,"long_name","Element neighbour indices")
    IERR = NF90_PUT_ATT(NC_FID,NBE_VID,"short_name","NBE")

    !!===== NTVE ==============================================================!
    IERR = NF90_DEF_VAR(NC_FID,"ntve",NF90_INT,STAT1DM,NTVE_VID)
    IERR = NF90_PUT_ATT(NC_FID,NTVE_VID,"long_name",&
         "Number of surrounding elements")
    IERR = NF90_PUT_ATT(NC_FID,NTVE_VID,"short_name","NTVE")

    !!===== NBVE ==============================================================!
    IERR = NF90_DEF_VAR(NC_FID,"nbve",NF90_INT,STAT2DM,NBVE_VID)
    IERR = NF90_PUT_ATT(NC_FID,NBVE_VID,"long_name",&
         "Indices of a nodes neighbouring elements")
    IERR = NF90_PUT_ATT(NC_FID,NBVE_VID,"short_name","NBVE")

    !!===== NBVT ==============================================================!
    IERR = NF90_DEF_VAR(NC_FID,"nbvt",NF90_INT,STAT2DM,NBVT_VID)
    IERR = NF90_PUT_ATT(NC_FID,NBVT_VID,"long_name", &
         "Index (1, 2 or 3) of a node in a neighbouring element")
    IERR = NF90_PUT_ATT(NC_FID,NBVT_VID,"short_name","NBVT")

    !!===== ISBCE =============================================================!
    IERR = NF90_DEF_VAR(NC_FID,"isbce",NF90_INT,STAT1DN,ISBCE_VID)
    IERR = NF90_PUT_ATT(NC_FID,ISBCE_VID,"long_name",&
         "Boundary element (1) or not (0)")
    IERR = NF90_PUT_ATT(NC_FID,ISBCE_VID,"short_name","ISBCE")

    !!===== ISONB =============================================================!
    IERR = NF90_DEF_VAR(NC_FID,"isonb",NF90_INT,STAT1DM,ISONB_VID)
    IERR = NF90_PUT_ATT(NC_FID,ISONB_VID,"long_name",&
         "Boundary node marker (0, 1 or 2)")
    IERR = NF90_PUT_ATT(NC_FID,ISONB_VID,"short_name","ISONB")

    !JW
    !!===== VXMIN, VYMIN ======================================================!
    !IERR = NF90_DEF_VAR(NC_FID,"vxmin",NF90_FLOAT,VXMIN_VID)
    !IERR = NF90_PUT_ATT(NC_FID,VXMIN_VID,"long_name","Minimum x-value in mesh")
    !IERR = NF90_PUT_ATT(NC_FID,VXMIN_VID,"short_name","VXMIN")

    !IERR = NF90_DEF_VAR(NC_FID,"vymin",NF90_FLOAT,VYMIN_VID)
    !IERR = NF90_PUT_ATT(NC_FID,VYMIN_VID,"long_name","Minimum y-value in mesh")
    !IERR = NF90_PUT_ATT(NC_FID,VXMIN_VID,"short_name","VYMIN")
    !JW

    !!===== Max_Elems =======================================================!
    IERR = NF90_DEF_VAR(NC_FID,"mx_nbr_elem",NF90_INT,MX_NBR_ELEM_VID)
    IERR = NF90_PUT_ATT(NC_FID,MX_NBR_ELEM_VID,"long_name", &
         "Maximum number of surrounding elements")
    IERR = NF90_PUT_ATT(NC_FID,MX_NBR_ELEM_VID,"short_name","Max_Elems")

    ! End definition section
    IERR = NF90_ENDDEF(NC_FID)

    IF(IERR /= NF90_NOERR) THEN
       WRITE(*,*)'ERROR OPENING ',TRIM(INFILE)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    ! Write NBE [N,3]
    ALLOCATE(START(2))
    ALLOCATE(CNT(2))
    START = 1
    CNT(1) = N_elems
    CNT(2) = 3
    IERR = NF90_PUT_VAR(NC_FID,NBE_VID,NBE,START,CNT)
    IF(IERR /= NF90_NOERR)THEN
       WRITE(*,*)'Error putting variable: nbe'
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF
    DEALLOCATE(START)
    DEALLOCATE(CNT)

    ! Write NTVE [M]
    ALLOCATE(START(1))
    ALLOCATE(CNT(1))
    START = 1
    CNT = N_verts
    IERR = NF90_PUT_VAR(NC_FID,NTVE_VID,NTVE,START,CNT)
    IF(IERR /= NF90_NOERR)THEN
       WRITE(*,*)'Error putting variable: ntve'
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    ! Write ISONB [M]
    IERR = NF90_PUT_VAR(NC_FID,ISONB_VID,ISONB,START,CNT)
    IF(IERR /= NF90_NOERR)THEN
       WRITE(*,*)'Error putting variable: isonb'
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF
    DEALLOCATE(START)
    DEALLOCATE(CNT)

    ! Write NBVE [M,9]
    ALLOCATE(START(2))
    ALLOCATE(CNT(2))
    START = 1
    CNT(1) = N_verts
    CNT(2) = Max_Elems
    IERR = NF90_PUT_VAR(NC_FID,NBVE_VID,NBVE,START,CNT)
    IF(IERR /= NF90_NOERR)THEN
       WRITE(*,*)'Error putting variable: nbve'
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    ! Write NBVT [M,9]
    IERR = NF90_PUT_VAR(NC_FID,NBVT_VID,NBVT,START,CNT)
    IF(IERR /= NF90_NOERR)THEN
       WRITE(*,*)'Error putting variable: nbvt'
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF
    DEALLOCATE(START)
    DEALLOCATE(CNT)

    ! Write ISBCE [N]
    ALLOCATE(START(1))
    ALLOCATE(CNT(1))
    START = 1
    CNT = N_elems
    IERR = NF90_PUT_VAR(NC_FID,ISBCE_VID,ISBCE,START,CNT)
    IF(IERR /= NF90_NOERR)THEN
       WRITE(*,*)'Error putting variable: isbce'
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF
    DEALLOCATE(START)
    DEALLOCATE(CNT)

    !JW
    ! Write VXMIN, VYMIN
    !IERR = NF90_PUT_VAR(NC_FID,VXMIN_VID,VXMIN)
    !IF(IERR /= NF90_NOERR)THEN
    !   WRITE(*,*)'Error putting variable: vxmin'
    !   WRITE(*,*)TRIM(NF90_STRERROR(IERR))
    !   STOP
    !END IF
    !IERR = NF90_PUT_VAR(NC_FID,VYMIN_VID,VYMIN)
    !IF(IERR /= NF90_NOERR)THEN
    !   WRITE(*,*)'Error putting variable: vymin'
    !   WRITE(*,*)TRIM(NF90_STRERROR(IERR))
    !   STOP
    !END IF
    !JW

    ! Write MX_NBR_ELEM
    IERR = NF90_PUT_VAR(NC_FID,MX_NBR_ELEM_VID,Max_Elems)
    IF(IERR /= NF90_NOERR)THEN
       WRITE(*,*)'Error putting variable: mx_nbr_elem'
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    ! Close file
    IERR = NF90_CLOSE(NC_FID)

  END SUBROUTINE NCD_WRITE_METRICS

  SUBROUTINE NCD_READ_METRICS(INFILE,SUCCESS)
    !--------------------------------------------------------------------------!
    ! READ DATA FROM METRICS NETCDF FILE
    !
    ! INPUTS:
    ! INFILE - netCDF input file string [max=100] [CHARACTER]
    !
    ! Copied from Pierre Cazenave's implementation in FVCOM offlag
    ! (08/01/2015) Plymouth Marine Laboratory.
    ! Adapted for FISCM by J. Ounsley 18/11/2016
    !--------------------------------------------------------------------------!
    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    CHARACTER(LEN=100), INTENT(IN)              :: INFILE
    LOGICAL, INTENT(OUT)                        :: SUCCESS
    !--------------------------------------------------------------------------!
    INTEGER            :: IERR    
    INTEGER            :: VXMIN_VID, VYMIN_VID, MX_NBR_ELEM_VID
    !--------------------------------------------------------------------------!
    INTEGER             :: NC_FID,varid
    CHARACTER(len=mstr) :: msg
    CHARACTER(len=fstr) :: dname
    INTEGER             :: temp_nele, temp_node

    SUCCESS = .FALSE.

    ! Open netCDF Datafile
    IERR = NF90_OPEN(TRIM(INFILE),NF90_NOWRITE,NC_FID)
    IF(IERR /= NF90_NOERR)THEN
       WRITE(*,*)'ERROR READING ',TRIM(INFILE)
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
    END IF


    ! JO Check that the mesh of this metrics file has the
    ! same number of nodes and elements as the forcing file
    ! determine number of elements
    msg = "dimension 'nele' not in the netcdf dataset"
    call ncdchk(nf90_inq_dimid(NC_FID, "nele", varid ),msg)
    call ncdchk(nf90_inquire_dimension(NC_FID, varid, dname, temp_nele ))

    ! determine number of nodes 
    msg = "dimension 'node' not in the netcdf dataset"
    call ncdchk(nf90_inq_dimid(NC_FID, "node", varid ),msg)
    call ncdchk(nf90_inquire_dimension(NC_FID, varid, dname, temp_node ))

    if (temp_nele /= N_elems .and. temp_node /= N_verts) then
       write(*,*) "WARNING::: variables 'node' and 'nele' in metrics file"
       write(*,*) "do not match values in forcing file."
       ! Do not read the metrics file
       return
    end if


    !!===== Max_Elems =======================================================!
    IERR = NF90_INQ_VARID(NC_FID,'mx_nbr_elem',MX_NBR_ELEM_VID)
    IF(IERR /=NF90_NOERR)THEN
       WRITE(*,*)'error getting variable id: mx_nbr_elem'
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF
    IERR = NF90_GET_VAR(NC_FID,MX_NBR_ELEM_VID,Max_Elems)
    IF(IERR /= NF90_NOERR)THEN
       WRITE(*,*)'Error getting variable: mx_nbr_elem'
       WRITE(*,*)TRIM(NF90_STRERROR(IERR))
       STOP
    END IF

    !!===== NBE ===============================================================!
    msg = "error reading nbe"
    call ncdchk(nf90_inq_varid(NC_FID,'nbe',varid),msg)
    call ncdchk(nf90_get_var(NC_FID,varid,nbe),msg)

    !!===== NTVE ==============================================================!
    msg = "error reading ntve"
    call ncdchk(nf90_inq_varid(NC_FID,'ntve',varid),msg)
    call ncdchk(nf90_get_var(NC_FID,varid,ntve),msg)

    !!===== NBVE ==============================================================!
    ! pica:
    ! We reallocate these here because if were reading the metrics from netCDF,
    ! the call to ALLOC_VARS happens before weve read in the data, so we end up
    ! with incorrectly shaped arrays.
    ALLOCATE(NBVE(N_verts,Max_Elems))
    msg = "error reading nbve"
    call ncdchk(nf90_inq_varid(NC_FID,'nbve',varid),msg)
    call ncdchk(nf90_get_var(NC_FID,varid,nbve),msg)

    !!===== NBVT ==============================================================!
    ! pica:
    ! We reallocate these here because if were reading the metrics from netCDF,
    ! the call to ALLOC_VARS happens before weve read in the data, so we end up
    ! with incorrectly shaped arrays.
    ALLOCATE(NBVT(N_verts,Max_Elems))
    msg = "error reading nbvt"
    call ncdchk(nf90_inq_varid(NC_FID,'nbvt',varid),msg)
    call ncdchk(nf90_get_var(NC_FID,varid,nbvt),msg)

    !!===== ISBCE =============================================================!
    msg = "error reading isbce"
    call ncdchk(nf90_inq_varid(NC_FID,'isbce',varid),msg)
    call ncdchk(nf90_get_var(NC_FID,varid,isbce),msg)

    !!===== ISONB =============================================================!
    msg = "error reading isonb"
    call ncdchk(nf90_inq_varid(NC_FID,'isonb',varid),msg)
    call ncdchk(nf90_get_var(NC_FID,varid,isonb),msg)

    !JW
    ! Read VXMIN, VYMIN
    !   CALL GETSVAR(NC_FID,LEN_TRIM('vxmin'),'vxmin',1,1,VXMIN)
    !   CALL GETSVAR(NC_FID,LEN_TRIM('vymin'),'vymin',1,1,VYMIN)
    ! pica: Because the VXMIN and VYMIN values are single values, we
    ! can''t use GETSVAR to write them to netCDF, so do so manually here
    ! instead.
    !IERR = NF90_INQ_VARID(NC_FID,'vxmin',VXMIN_VID)
    !IF(IERR /=NF90_NOERR)THEN
    !   WRITE(*,*)'error getting variable id: vxmin'
    !   WRITE(*,*)TRIM(NF90_STRERROR(IERR))
    !   STOP
    !END IF
    !IERR = NF90_GET_VAR(NC_FID,VXMIN_VID,VXMIN)
    !IF(IERR /= NF90_NOERR)THEN
    !   WRITE(*,*)'Error getting variable: vxmin'
    !   WRITE(*,*)TRIM(NF90_STRERROR(IERR))
    !   STOP
    !END IF
    !IERR = NF90_INQ_VARID(NC_FID,'vymin',VYMIN_VID)
    !IF(IERR /=NF90_NOERR)THEN
    !   WRITE(*,*)'error getting variable id: vymin'
    !   WRITE(*,*)TRIM(NF90_STRERROR(IERR))
    !   STOP
    !END IF
    !IERR = NF90_GET_VAR(NC_FID,VYMIN_VID,VYMIN)
    !IF(IERR /= NF90_NOERR)THEN
    !   WRITE(*,*)'Error getting variable: vymin'
    !   WRITE(*,*)TRIM(NF90_STRERROR(IERR))
    !   STOP
    !END IF
    !JW

    ! Close file
    IERR = NF90_CLOSE(NC_FID)

    SUCCESS = .TRUE.

    RETURN
  END SUBROUTINE NCD_READ_METRICS

End Module Ocean_Model 
