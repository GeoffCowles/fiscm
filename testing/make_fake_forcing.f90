!=======================================================================
! Junk program to make fake forcing to test advection, biology, and  
!   vertical diffusivity in FISCM
!
! Comments: 
!    Originally from the OSCAR shallow water equation solver
!
! Compile:
! gfortran -c make_fake_forcing.f90 -I/usr/local/netcdf/gfortran/include ; /bin/sh /usr/local/netcdf/gfortran/bin/libtool  --mode=link gfortran make_fake_forcing.o -o make_fake_forcing -L/usr/local/netcdf/gfortran/lib -lnetcdf
!
! gfortran -c make_fake_forcing.f90 -I/opt/netcdf3/gfortran/include ; gfortran  make_fake_forcing.o -o make_fake_forcing -L/opt/netcdf3/gfortran/lib  -lnetcdf
! /usr/local/bin/gfortran -c make_fake_forcing.f90 -I/opt/netcdf3/gnu32/include; /usr/local/bin/gfortran -o make_fake_forcing make_fake_forcing.o -L/opt/netcdf3/gnu32/lib -lnetcdf
!
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!
!=======================================================================

Program make_fake_forcing
  use netcdf
  implicit none
  integer,  parameter :: sp = selected_real_kind(15,307)

  !mesh dimensions and range (use 51 for case_type ==1)
  integer,  parameter :: il = 51
  integer,  parameter :: jl = 51
  real(sp), parameter :: xt = 10000.
  real(sp), parameter :: yt = 10000.
  integer :: N_levels
  integer :: N_layers 

  !time range
  integer, parameter :: nframes = 11
 
  !variables  
  real(sp), allocatable :: x(:)
  real(sp), allocatable :: y(:)
  real(sp), allocatable :: xc(:)
  real(sp), allocatable :: yc(:)
  real(sp), allocatable :: h(:)
  real(sp), allocatable :: zeta(:)
  integer , allocatable :: tri(:,:)

  real(sp), allocatable :: u(:,:)
  real(sp), allocatable :: v(:,:)
  real(sp), allocatable :: omega(:,:)
  real(sp), allocatable :: t(:,:)
  real(sp), allocatable :: s(:,:)
  real(sp), allocatable :: kh(:,:)
  real(sp), allocatable :: siglay(:,:)
  real(sp), allocatable :: siglev(:,:)
  real(sp), allocatable :: zobs(:)
  real(sp), allocatable :: khobs(:)
  
  integer  :: nx,ny,nnode,ncell,i,j,k,n,nn,n2,maxp,case_type
  real(sp) :: dx,dy,xmin,xmax,ymin,ymax,theta,pi,bathy
  integer  :: x_vid,y_vid,h_vid,tri_vid,u_vid,v_vid,t_vid,time_vid,siglay_vid,sigshift_vid,siglev_vid
  integer  :: omega_vid,kh_vid
  integer  :: nele_did,node_did,three_did,zeta_vid,siglay_did,siglev_did,time_did
  integer  :: ofid

  integer  :: N_verts,N_elems
  !------------------------------------------------------------------------------
  ! set case type 
  !  1.) oscillating solid body rotation with increasing temp
  !  2.) still flow with steady, reasonable kh profile, non-constant sigma
  !  3.) still flow with steady, reasonable kh profile, constant sigma spacing
  !  4.) still flow with constant in times/space kh profile, constant sigma 
  !------------------------------------------------------------------------------
  case_type = 3

  if(case_type == 1)then
    N_levels = 5
    N_layers = N_levels-1
    bathy    = 10.
  elseif(case_type ==2)then
    N_levels = 48
    N_layers = N_levels-1
    bathy    = 40.
    allocate(khobs(N_levels))
    allocate(zobs(N_levels))
    open(unit=33,file='kh_data.dat',form='formatted')
    do i=1,N_levels
      read(33,*)zobs(i),khobs(i)
    end do
    close(33)
  elseif(case_type ==3)then
    N_levels = 81
    N_layers = N_levels-1
    bathy    = 40.
    allocate(khobs(N_levels))
    allocate(zobs(N_levels))
    open(unit=33,file='kh_data2.dat',form='formatted')
    do i=1,N_levels
      read(33,*)zobs(i),khobs(i)
    end do
    close(33)
  elseif(case_type ==4)then
    N_levels = 81
    N_layers = N_levels-1
    bathy    = 40.
    allocate(khobs(N_levels))
    allocate(zobs(N_levels))
    open(unit=33,file='kh_data3.dat',form='formatted')
    do i=1,N_levels
      read(33,*)zobs(i),khobs(i)
    end do
    close(33)
  endif
  

  !-------------------------------------------------------
  ! generate a mesh
  !-------------------------------------------------------
  ny = jl-1
  nx = il-1
  dx = xt/float(nx)
  dy = yt/float(ny)
  nn = 0
  nnode = il*jl
  ncell = nx*ny*2
  allocate(x(il*jl))
  allocate(y(il*jl))
  allocate(h(il*jl))
  allocate(zeta(il*jl))
  do i=1,il
    do j=1,jl
      nn = nn + 1
        x(nn) = float(i-1)*dx
        y(nn) = float(j-1)*dy
    end do
  end do
  xmin = minval(x(:))
  xmax = maxval(x(:))
  ymin = minval(y(:))
  ymax = maxval(y(:))
  x = x - .5*(xmax-xmin)
  y = y - .5*(ymax-ymin)

  !bathymetry
  N_verts = il*jl
  do i=1,N_verts
    h(i) = bathy
    zeta(i) = 0.0
  end do

  allocate(tri(nx*ny*2,3))
  nn = 1
  do i=1,nx
    do j=1,ny
      n2 = nn + 1
      tri(nn,1) = j + (i-1)*jl
      tri(nn,2) = j + i*jl
      tri(nn,3) = j+1 + (i-1)*jl
      tri(n2,1) = j + i*jl
      tri(n2,2) = j+1 + i*jl
      tri(n2,3) = j+1 + (i-1)*jl
      nn = nn + 2
    end do
  end do
  maxp = max( maxval(tri(:,1)),maxval(tri(:,2)),maxval(tri(:,3)))

  N_elems = nx*ny*2
  allocate(xc(N_elems)) ; xc = 0.0_sp
  allocate(yc(N_elems)) ; yc = 0.0_sp
  do i=1,N_elems
    xc(i) = (x(tri(i,1)) + x(tri(i,2)) + x(tri(i,3)) ) /3.
    yc(i) = (y(tri(i,1)) + y(tri(i,2)) + y(tri(i,3)) ) /3.
  end do

  !-----------------------------------------------------------
  ! allocate data 
  !-----------------------------------------------------------
  allocate(u(N_elems,N_layers)) ; u = 0.0_sp
  allocate(v(N_elems,N_layers)) ; v = 0.0_sp
  allocate(omega(N_verts,N_levels)) ; omega = 0.0_sp
  allocate(kh(N_verts,N_levels)) ; kh = 0.0_sp
  allocate(t(N_verts,N_layers)) ; t = 0.0_sp
  allocate(siglay(N_verts,N_layers)) ; siglay = 0.0_sp
  allocate(siglev(N_verts,N_layers+1)) ; siglev = 0.0_sp

  do i=1,N_verts
    do k=1,N_layers+1
      if(case_type == 1)then
        siglev(i,k) = -float(k-1)/float(N_layers-1) 
      elseif(case_type ==2 .or. case_type==3 .or. case_type==4)then
        siglev(i,k) = -zobs(k)/bathy
      endif
    end do
    do k=1,N_layers
      siglay(i,k) = .5*(siglev(i,k)+siglev(i,k+1))
    end do
  end do


  !-----------------------------------------------------------
  ! netcdf header
  !-----------------------------------------------------------

  call cfcheck( nf90_create("fake_forcing.nc",nf90_clobber,ofid) ) 

  !dimensions
  call cfcheck(nf90_def_dim(ofid,"time",nf90_unlimited, time_did) )
  call cfcheck(nf90_def_dim(ofid,"nele",N_elems, nele_did) )
  call cfcheck(nf90_def_dim(ofid,"node",N_verts, node_did) )
  call cfcheck(nf90_def_dim(ofid,"siglay",N_layers, siglay_did) )
  call cfcheck(nf90_def_dim(ofid,"siglev",N_layers+1, siglev_did) )
  call cfcheck(nf90_def_dim(ofid,"three",3, three_did) )

  !time
  call cfcheck( nf90_def_var(ofid,"time",nf90_float,(/time_did/), time_vid) )
  call cfcheck( nf90_put_att(ofid, time_vid,"long_name","Time") )
  call cfcheck( nf90_put_att(ofid, time_vid,"units","days after 00:00:00") )
  call cfcheck( nf90_put_att(ofid, time_vid,"calendar","none") )


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

  !siglay shift
  call cfcheck( nf90_def_var(ofid,"siglay_shift",nf90_float,(/node_did,siglay_did/), sigshift_vid) )
  call cfcheck( nf90_put_att(ofid, sigshift_vid,"long_name","Shifted Sigma Layers") )

  !siglev
  call cfcheck( nf90_def_var(ofid,"siglev",nf90_float,(/node_did,siglev_did/), siglev_vid) )
  call cfcheck( nf90_put_att(ofid, siglev_vid,"long_name","Sigma Layers") )
  call cfcheck( nf90_put_att(ofid, siglev_vid,"standard_name","ocean_sigma/general_coordinate") )
  call cfcheck( nf90_put_att(ofid, siglev_vid,"positive","down") )
  call cfcheck( nf90_put_att(ofid, siglev_vid,"valid_min","-1") )
  call cfcheck( nf90_put_att(ofid, siglev_vid,"valid_max","0") )
  call cfcheck( nf90_put_att(ofid, siglev_vid,"formula_terms","siglev:siglev eta:zeta depth:depth") )

  !siglay
  call cfcheck( nf90_def_var(ofid,"siglay",nf90_float,(/node_did,siglay_did/), siglay_vid) )
  call cfcheck( nf90_put_att(ofid, siglay_vid,"long_name","Sigma Layers") )
  call cfcheck( nf90_put_att(ofid, siglay_vid,"standard_name","ocean_sigma/general_coordinate") )
  call cfcheck( nf90_put_att(ofid, siglay_vid,"positive","down") )
  call cfcheck( nf90_put_att(ofid, siglay_vid,"valid_min","-1") )
  call cfcheck( nf90_put_att(ofid, siglay_vid,"valid_max","0") )
  call cfcheck( nf90_put_att(ofid, siglay_vid,"formula_terms","siglay:siglay eta:zeta depth:depth") )

  !nv
  call cfcheck( nf90_def_var(ofid,"nv",nf90_int,(/nele_did,three_did/), tri_vid) )
  call cfcheck( nf90_put_att(ofid, tri_vid,"long_name","nodes surrounding element") )

  !el
  call cfcheck( nf90_def_var(ofid,"zeta",nf90_float,(/node_did,time_did/), zeta_vid) )
  call cfcheck( nf90_put_att(ofid, zeta_vid,"long_name","Water Surface Elevation") )
  call cfcheck( nf90_put_att(ofid, zeta_vid,"units","meters") )
  call cfcheck( nf90_put_att(ofid, zeta_vid,"positive","up") )
  call cfcheck( nf90_put_att(ofid, zeta_vid,"standard_name","sea_surface_elevation") )
  call cfcheck( nf90_put_att(ofid, zeta_vid,"type","data") )

  !u  
  call cfcheck( nf90_def_var(ofid,"u",nf90_float,(/nele_did,siglay_did,time_did/), u_vid) )
  call cfcheck( nf90_put_att(ofid, u_vid,"long_name","Eastward Water Velocity") )
  call cfcheck( nf90_put_att(ofid, u_vid,"units","meters s-1") )
  call cfcheck( nf90_put_att(ofid, u_vid,"grid","fvcom_grid") )
  call cfcheck( nf90_put_att(ofid, u_vid,"type","data") )

  !v  
  call cfcheck( nf90_def_var(ofid,"v",nf90_float,(/nele_did,siglay_did,time_did/), v_vid) )
  call cfcheck( nf90_put_att(ofid, v_vid,"long_name","Northward water velocity") )
  call cfcheck( nf90_put_att(ofid, v_vid,"units","meters s-1") )
  call cfcheck( nf90_put_att(ofid, v_vid,"grid","fvcom_grid") )
  call cfcheck( nf90_put_att(ofid, v_vid,"type","data") )

  !omega
  call cfcheck( nf90_def_var(ofid,"ww",nf90_float,(/node_did,siglev_did,time_did/), omega_vid) )
  call cfcheck( nf90_put_att(ofid, omega_vid,"long_name","Vertical Sigma Coordinate Velocity") )
  call cfcheck( nf90_put_att(ofid, omega_vid,"units","s-1") )
  call cfcheck( nf90_put_att(ofid, omega_vid,"grid","fvcom_grid") )
  call cfcheck( nf90_put_att(ofid, omega_vid,"type","data") )

  !kh    
  call cfcheck( nf90_def_var(ofid,"kh",nf90_float,(/node_did,siglev_did,time_did/), kh_vid) )
  call cfcheck( nf90_put_att(ofid, kh_vid,"long_name","Turbulent Eddy Viscosity For Scalars") )
  call cfcheck( nf90_put_att(ofid, kh_vid,"units","m 2 s-1") )
  call cfcheck( nf90_put_att(ofid, kh_vid,"grid","fvcom_grid") )
  call cfcheck( nf90_put_att(ofid, kh_vid,"coordinates","x y") )
  call cfcheck( nf90_put_att(ofid, kh_vid,"type","data") )

  !t  
  call cfcheck( nf90_def_var(ofid,"temp",nf90_float,(/node_did,siglay_did,time_did/), t_vid) )
  call cfcheck( nf90_put_att(ofid, t_vid,"long_name","sea_water_temperature") )
  call cfcheck( nf90_put_att(ofid, t_vid,"units","degrees_C") )
  call cfcheck( nf90_put_att(ofid, t_vid,"grid","fvcom_grid") )
  call cfcheck( nf90_put_att(ofid, t_vid,"type","data") )

  !globals
  call cfcheck( nf90_put_att(ofid,nf90_global,"title"     ,"dummy") )
  call cfcheck( nf90_put_att(ofid,nf90_global,"institution"     ,"School for Marine Science and Technology") )
  call cfcheck( nf90_put_att(ofid,nf90_global,"history"     ,"model started at: 18/07/2008   18:45") )
  call cfcheck( nf90_put_att(ofid,nf90_global,"references","http://fvcom.smast.umassd.edu")) 
  call cfcheck( nf90_put_att(ofid,nf90_global,"source"     ,"FVCOM_2.6") )
  call cfcheck( nf90_put_att(ofid,nf90_global,"Conventions"     ,"CF-1.0") )

  !end definitions
  call cfcheck( nf90_enddef(ofid) )

  !------------------------------------------------------------------------------
  ! dump static data 
  !------------------------------------------------------------------------------
  call cfcheck( nf90_put_var(ofid,x_vid,x)) 
  call cfcheck( nf90_put_var(ofid,y_vid,y)) 
  call cfcheck( nf90_put_var(ofid,h_vid,h)) 
  call cfcheck( nf90_put_var(ofid,siglay_vid,siglay)) 
  call cfcheck( nf90_put_var(ofid,sigshift_vid,siglay)) 
  call cfcheck( nf90_put_var(ofid,siglev_vid,siglev)) 
  call cfcheck( nf90_put_var(ofid,tri_vid,tri,START=(/1,1/))) 

  !------------------------------------------------------------------------------
  ! loop over time, make up data
  !------------------------------------------------------------------------------


  !oscillating solid body rotation with temperature increasing outward from core and with time
  if(case_type == 1)then
    do i=1,N_elems
      theta = atan2(yc(i),xc(i))
      u(i,:) = (-sin(theta) * sqrt(  xc(i)**2 + yc(i)**2 ))/max(xmax,ymax)  
      v(i,:) = ( cos(theta) * sqrt(  xc(i)**2 + yc(i)**2 ))/max(ymax,xmax) 
    end do
    xmax = maxval(x)
    ymax = maxval(y)
    do i=1,N_verts
      t(i,:) = 4.*sqrt(  x(i)**2 + y(i)**2 )/sqrt(xmax**2 + ymax**2) + 12.
      kh(i,:) = 0.0
    end do
  elseif(case_type == 2 .or. case_type ==3 .or. case_type == 4)then
    do i=1,N_verts
    do k=1,N_levels
      kh(i,k) = khobs(k)
    end do
    do k=1,N_layers
      t(i,k) = 30+20*siglay(i,k)
    end do
    end do
  endif



  do n=1,nframes
    call cfcheck( nf90_put_var(ofid,time_vid,float(n-1),START = (/n/))) 
    call cfcheck( nf90_put_var(ofid,zeta_vid,h*0.0,START = (/1,n/))) 
    call cfcheck( nf90_put_var(ofid,u_vid,u,START = (/1,1,n/))) 
    call cfcheck( nf90_put_var(ofid,v_vid,v,START = (/1,1,n/))) 
    call cfcheck( nf90_put_var(ofid,omega_vid,omega,START = (/1,1,n/))) 
    call cfcheck( nf90_put_var(ofid,kh_vid,kh,START = (/1,1,n/))) 
    call cfcheck( nf90_put_var(ofid,t_vid,t,START = (/1,1,n/))) 

    if(case_type == 1)then
      pi = 4*atan(1.0) 
      u = u*cos(pi*float(n-1))
      v = v*cos(pi*float(n-1))
      t = t + 1. 
    endif
  end do


  !close the file
  call cfcheck( nf90_close(ofid) )





end 

!-----------------------------------------------
! runtime errors  - netcdf
!-----------------------------------------------
subroutine cfcheck(status)
    use netcdf
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop
    end if
end subroutine cfcheck

