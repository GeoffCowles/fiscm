!=======================================================================
! Fiscm NetCDF Forcing Routines 
!
! Description
!   Setup and access NetCDF forcing 
!    - use frames (time slices) built from containers
!    - current frame interpolated from bracketing frames
!    
! Comments: 
!    - Modify to use D. Stuebe's NetCDF package once binned random-walk
!      and NURBS field reconstruction tested and working
!    - Requires FORTRAN90 NetCDF 3x libraries
!
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!
!=======================================================================
module forcing
use gparms
use mod_igroup
use netcdf
implicit none

logical :: debug = .false.

integer, parameter :: flt_type = 1
integer, parameter :: int_type = 0

integer, parameter :: CONSTANT = 0
integer, parameter :: TIME_VARYING = 1

integer, parameter :: FRAME_SETUP = 0
integer, parameter :: FRAME_STATS = 1

integer, private :: ffile_id,nDim,nVar,nAtt,uDid,fNum,uD_extant
character(len=fstr), private :: ffile(max_nf)

real(sp) :: forcing_beg_time
real(sp) :: forcing_end_time
real(sp) :: forcing_deltaT
integer  :: nframes
integer  :: lasti1,lasti2 
real(sp), allocatable :: ftimes(:),tmptime(:)
integer , allocatable :: fid_ft(:)
integer , allocatable :: rec_ft(:)
integer , allocatable :: fnr_ft(:)
type container
  integer :: dtype
  integer :: ttype
  integer :: ndims
  integer, pointer :: dims(:) 
  integer :: varid
  character(len=fstr) :: vname
  !can do in with single 1-d array if we could reshape allocatable
  real(sp), pointer :: f1(:)
  real(sp), pointer :: f2(:,:)
  real(sp), pointer :: f3(:,:,:)
  real(sp), pointer :: f4(:,:,:,:)
  integer , pointer :: i1(:)
  integer , pointer :: i2(:,:)
  integer , pointer :: i3(:,:,:)
  integer , pointer :: i4(:,:,:,:)
end type container

type dataframe 
  real(sp) :: time
  real(sp) :: ncfid
  real(sp) :: iframe
  integer  :: nvars
  integer  :: id
  logical  :: initial_read
  type(container), pointer :: fdata(:)   ![nvars]
end type dataframe

integer         :: now = 3 !frame 3 is current data
integer         :: nex = 4 !frame 4 is next data
type(dataframe) :: frame(4) 

!overload ">" to be able to compare the order of time frames in time
interface operator(>)
  module procedure frame_order
end interface

interface get_forcing
  module procedure get_forcing_f1
  module procedure get_forcing_f2
end interface

contains

!========================================================================
! update the model forcing to time (t)
!   using linear interpolation
!========================================================================
subroutine update_forcing(t,iframe0)
  implicit none
  real(sp), intent(in) :: t
  integer :: i1,i2,i1f,i2f,iframe0
  logical inew

  call bracket(t,i1,i2) 
  i1f = 0 ; i2f = 0
  
  !determine if either i1 or i2 are already loaded 
  if(frame(1)%iframe == i1) i1f = 1  
  if(frame(2)%iframe == i1) i1f = 2  
  if(frame(1)%iframe == i2) i2f = 1  
  if(frame(2)%iframe == i2) i2f = 2  

  !load i1 if necessary 
  if(i1f == 0)then
    if(i2f == 1) then
      call read_frame(frame(2),i1)
      i1f = 2
    else
      call read_frame(frame(1),i1)
      i1f = 1
    endif
  endif

  !load i2 if necessary 
  if(i2f == 0)then
    if(i1f == 1) then
      call read_frame(frame(2),i2)
    else
      call read_frame(frame(1),i2)
    endif
  endif

  call interp_two_frames(frame(iframe0),t,frame(1),frame(2))
  !call frame_info(frame(1),FRAME_STATS)
  !call frame_info(frame(2),FRAME_STATS)
 
end subroutine update_forcing

!========================================================================
! setup frames (set variable names, allocate, etc)
!   we need 3 frames, one to store two timeframes
!   from the netcdf files, and one to store data
!   interpolated to time (t)
!========================================================================
subroutine setup_forcing(nvars,varlist)
  implicit none
  integer, intent(in) :: nvars
  character(len=*)    :: varlist(nvars)
  integer i
  frame(1) = frame_(nvars,varlist,1)
  frame(2) = frame_(nvars,varlist,2)
  frame(3) = frame_(nvars,varlist,3)
  frame(4) = frame_(nvars,varlist,4)

  call frame_info(frame(1),FRAME_SETUP)

  !initialize frames 
  lasti1 = -1
  lasti2 = -1

end subroutine setup_forcing 

!========================================================================
! dump frame info to screen 
!========================================================================
subroutine frame_info(frame,itype)
  use utilities, only : drawline
  implicit none
  type(dataframe), intent(in) :: frame
  integer :: itype
  !----------------------------
  integer i,j
  real(sp) :: fmax,fmin,fave 
  integer :: ndims,dims(4)
  character(len=fstr) :: vname
  logical :: tv
 
  if(frame%nvars < 1)return

  !------------------------------------------------------------
  ! dump frame vars and basic information (dimensions and type)
  !------------------------------------------------------------
  if(itype == FRAME_SETUP)then
    call drawline("-")
    write(*,*)' forcing var   |   | dims |   dim1  dim2  dim3  dim4' 
    call drawline("-")
    do i=1,frame%nvars 
      tv = .false.
      if(frame%fdata(i)%ttype == TIME_VARYING) tv = .true.
      ndims = frame%fdata(i)%ndims
      vname = frame%fdata(i)%vname
      dims  = 0
      dims(1:ndims)  = frame%fdata(i)%dims(1:ndims)
      write(*,'(A15,A3,L1,A3,I4,A3,4I6)')vname,' | ',tv,' | ',ndims,' | ',(dims(j),j=1,ndims)
    end do
  !------------------------------------------------------------
  ! dump frame vars and stats on variables 
  !------------------------------------------------------------
  elseif(itype == FRAME_STATS) then
    call drawline("-")
    write(*,*)' forcing var   |    min |    max |    ave'  
    call drawline("-")
    do i=1,frame%nvars
      fmax = calc_stat(frame,frame%fdata(i)%vname,'max')
      fmin = calc_stat(frame,frame%fdata(i)%vname,'min')
      fave = calc_stat(frame,frame%fdata(i)%vname,'mean')
      vname = frame%fdata(i)%vname
      write(*,'(A15,A3,F6.2,A3,F6.2,A3,F6.2)')vname,' | ',fmin,' | ',fmax,' | ',fave
    end do
  else
    write(*,*)'fatal error in frame_info'
    write(*,*)'itype must be: ',FRAME_SETUP,' or ',FRAME_STATS
    stop
  endif
  call drawline("-")
  
end subroutine frame_info

!========================================================================
! frame type constructor 
!========================================================================
function frame_(nvars,varlist,id) result(frame)
  implicit none
  integer, intent(in) :: nvars
  character(len=fstr) :: varlist(nvars)
  integer, intent(in) :: id
  type(dataframe) :: frame
  !--------------------------------------
  integer :: dim1,dim2,i,i1,i2,i3,j,jj
  integer :: fid,nAtts,vid,ndims,xtype,ttype,dtype
  character(len=mstr) :: msg,msg2
  character(fstr) :: vname,dname
  integer :: dimids(NF90_MAX_VAR_DIMS)
  integer :: alldims(NF90_MAX_VAR_DIMS)
  integer :: dims(NF90_MAX_VAR_DIMS)
  integer :: varids(nvars)


  frame%id     = id
  frame%nvars  = nvars
  frame%time   = -hugenum
  frame%iframe = -1
  ffile_id=fid_ft(1)
  
  !return if nvars = 0
  if(nvars < 1)return
  
  !make sure all the variables exist
  do i=1,nvars
    if(debug)write(*,*)'making sure: ',trim(varlist(i)),' exists'
    msg  = "error: var ["//trim(varlist(i))//"]"
    msg2 = " not found in "//trim(ffile(1))
    call ncdchk( nf90_inq_varid(ffile_id,varlist(i),varids(i)),msg,msg2 ) 
  end do

  !set initial values for time and frame_num
  frame%ncfid = ffile_id 

  !allocate primary dataspace
  allocate(frame%fdata(nvars)) 

  !for each variable, set name, get dims, allocate
  do i=1,nvars
    vid = varids(i) 
    call ncdchk(nf90_inquire_variable(ffile_id,vid,vname,xtype,ndims,dimids,nAtts))
    alldims = 0
    dims = 0
    jj = 0
    !loop over dimensions
    do j=1,ndims
      jj = jj + 1
      call ncdchk(nf90_inquire_dimension(ffile_id, dimids(j), dname, alldims(j)) )
      if(dimids(j) == uDid)then
         ndims = ndims - 1
         ttype = TIME_VARYING
      else
         dims(jj) = alldims(j)
         ttype = CONSTANT
      end if
    end do

    !set vartype and check
    if(xtype == NF90_FLOAT)then
      dtype = flt_type
    elseif(xtype == NF90_INT)then
      dtype = int_type
    else
      write(*,*)'fatal error in frame_: type: ',xtype,' not defined'
      stop
    endif
       
    !set dimensions and check
    if(ndims > 4)then
      write(*,*)'fatal error in frame_: number of var dimensions exceeds 4 '
      stop
    endif
  
    !set container
    frame%fdata(i) = container_(dtype,ttype,vid,vname,ndims,dims)  

    !set initial read to false
    frame%initial_read = .false.
    
  end do

  
end function frame_

!========================================================================
! read data into a frame from file frame i 
!========================================================================
subroutine read_frame(frame,f) 
  implicit none
  type(dataframe), intent(inout) :: frame
  integer, intent(in) :: f
  integer :: i,varid,ndims,ttype,stat
  integer, allocatable, dimension(:) :: start,vsize

  !make sure frame is in bounds
  if(f > nframes .or. f < 1)then
    write(*,*)'error in read_frame'
    write(*,*)'fatal error in read_frame'
    write(*,*)'cannot read frame: ',f,' from file: ',trim(ffile(fnr_ft(f)))
    write(*,*)'max frames in netcdf file is: ',nframes
    stop 
  endif
  frame%iframe = f
  frame%time   = ftimes(f)

  if(debug)write(*,*)'============data frame: ',frame%id,' reading frame ',f

  !loop over vars, read data into frame
  do i=1,frame%nvars
    if(frame%initial_read.and. (frame%fdata(i)%ttype == CONSTANT))cycle
    if(debug)write(*,*)'reading data into: ',trim(frame%fdata(i)%vname)

    !set dimensions, stride, start
    varid = frame%fdata(i)%varid
    ndims = frame%fdata(i)%ndims
    ttype = frame%fdata(i)%ttype
    ffile_id=fid_ft(f)    

    !static var, read in entirety 
    if(ttype == CONSTANT)then

      if(ndims==1)call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f1) )
      if(ndims==2)call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f2) )
      if(ndims==3)call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f3) )
      if(ndims==4)call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f4) )

    else !time-varying variable

      allocate(start(ndims+ttype))
      allocate(vsize(ndims+ttype))
      start(1:ndims) = 1
      vsize(1:ndims) = frame%fdata(i)%dims(1:ndims)
      start(ndims+1) = rec_ft(f)
      vsize(ndims+1) = 1
      if(ndims==1) call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f1,start,vsize) )
      if(ndims==2) call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f2,start,vsize) )
      if(ndims==3) call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f3,start,vsize) )
      if(ndims==4) call ncdchk( nf90_get_var(ffile_id, varid, frame%fdata(i)%f4,start,vsize) )
      deallocate(start)
      deallocate(vsize)

    endif
    
  end do
  
  !set initial_frame to .true. will no longer read non-time-varying vars
  frame%initial_read = .true.

!  call frame_info(frame,FRAME_STATS)
  

end subroutine read_frame
  
!========================================================================
! linearly interpolate between two data frames 
!========================================================================
subroutine interp_two_frames(frame,t,frame1,frame2) !result(frame)
  implicit none
  type(dataframe), intent(inout) :: frame
  real(sp), intent(in) :: t
  type(dataframe), intent(in) :: frame1,frame2
!  type(dataframe) :: frame
  !----------------------------------------
  real(sp) :: t1,t2,c1,c2
  integer  :: ndims,dtype,i
  integer  :: d1,d2,d3,d4

  !set interpolation coefficients
  t1 = frame1%time
  t2 = frame2%time
  if(t2 > t1)then
    c2 = (t - t1)/(t2-t1)
    c1 = 1.-c2
  else
    c1 = (t - t2)/(t1-t2)
    c2 = 1.-c1
  endif
 
  !loop over vars, interpolate to intermediate frame 
  do i=1,frame1%nvars
    dtype = frame1%fdata(i)%dtype
    ndims = frame1%fdata(i)%ndims
    if(dtype == flt_type)then
      if(ndims == 1) frame%fdata(i)%f1 = c1*frame1%fdata(i)%f1 + c2*frame2%fdata(i)%f1
      if(ndims == 2) frame%fdata(i)%f2 = c1*frame1%fdata(i)%f2 + c2*frame2%fdata(i)%f2 
      if(ndims == 3) frame%fdata(i)%f3 = c1*frame1%fdata(i)%f3 + c2*frame2%fdata(i)%f3
      if(ndims == 4) frame%fdata(i)%f4 = c1*frame1%fdata(i)%f4 + c2*frame2%fdata(i)%f4
    else
      if(ndims == 1) frame%fdata(i)%i1 = c1*frame1%fdata(i)%i1 + c2*frame2%fdata(i)%i1
      if(ndims == 2) frame%fdata(i)%i2 = c1*frame1%fdata(i)%i2 + c2*frame2%fdata(i)%i2
      if(ndims == 3) frame%fdata(i)%i3 = c1*frame1%fdata(i)%i3 + c2*frame2%fdata(i)%i3
      if(ndims == 4) frame%fdata(i)%i4 = c1*frame1%fdata(i)%i4 + c2*frame2%fdata(i)%i4
    endif
  enddo

end subroutine interp_two_frames
!end function interp_two_frames

!========================================================================
! check frame order  
! returns true if frame2 is after frame1 
!========================================================================
function frame_order(frame1,frame2) result(order)
  implicit none
  type(dataframe), intent(in) :: frame1
  type(dataframe), intent(in) :: frame2
  logical :: order
  order = .false.
  if(frame2%time > frame1%time) order = .true. 
end function frame_order

!========================================================================
! container constructor
!   set name
!   set varid
!   set dimensions
!   allocate space
!========================================================================
function container_(dtype,ttype,vid,vname,ndims,dims) result(c)
  implicit none
  integer, intent(in) :: dtype 
  integer, intent(in) :: ttype 
  integer, intent(in) :: vid
  character(len=*)    :: vname
  integer, intent(in) :: ndims
  integer, intent(in) :: dims(ndims)
  type(container) :: c 
  integer :: adims(4)

  !check dimensions
  if(minval(dims) < 1)then 
    write(*,*)'fatal error in container_'
    write(*,*)'input dimensions provided < 0'
    write(*,*)dims
    stop
  endif

  !set and check the number of dimensions
  allocate(c%dims(ndims)) ; c%dims = dims 
  adims = 1
  adims(1:ndims) = dims

  c%dtype = dtype
  c%ttype = ttype
  c%vname = vname
  c%ndims = ndims
  c%varid = vid

  !allocate dataspace - float and initialize to zero
  if(dtype == flt_type)then
    if(ndims==1)then
       allocate(c%f1(adims(1)))  
       c%f1 = 0.0
    elseif(ndims==2)then
       allocate(c%f2(adims(1),adims(2)))
       c%f2 = 0.0
    elseif(ndims==3)then
       allocate(c%f3(adims(1),adims(2),adims(3)))
       c%f3 = 0.0
    elseif(ndims==4)then
       allocate(c%f4(adims(1),adims(2),adims(3),adims(4)))
       c%f4 = 0.0
    endif
  elseif(dtype == int_type)then
    if(ndims==1)then
       allocate(c%i1(adims(1)))  
       c%i1 = 0.0
    elseif(ndims==2)then
       allocate(c%i2(adims(1),adims(2)))
       c%i2 = 0.0
    elseif(ndims==3)then
       allocate(c%i3(adims(1),adims(2),adims(3)))
       c%i3 = 0.0
    elseif(ndims==4)then
       allocate(c%i4(adims(1),adims(2),adims(3),adims(4)))
       c%i4 = 0.0
    endif
  else
    write(*,*)'fatal error in container_'
    write(*,*)'data type: ',dtype,' not a valid type'
    stop
  endif

end function container_

function ncdscan(status,message,message2) result(ncd_okay)
  implicit none
  integer, intent ( in) :: status
  character(len=*), optional   :: message
  character(len=*), optional   :: message2
  logical ncd_okay  
  ncd_okay = .true. 

  if(status /= nf90_noerr) then
    ncd_okay = .false.
    if(present(message))print *, trim(message)
    if(present(message2))print *, trim(message2)
    print *, trim(nf90_strerror(status))
  end if
end function ncdscan  

subroutine ncdchk(status,message,message2)
  implicit none
  integer, intent ( in) :: status
  character(len=*), optional   :: message
  character(len=*), optional   :: message2

  if(status /= nf90_noerr) then
    if(present(message))print *, trim(message)
    if(present(message2))print *, trim(message2)
    print *, trim(nf90_strerror(status))
    stop
  end if
end subroutine ncdchk 

!========================================================================
! calculate statistic /stat_type/ from variable /vname/ in frame /frame/ 
!  options are:  min,max,mean
!========================================================================
function calc_stat(frame,vname,stat_type) result(stat)
  implicit none
  type(dataframe), intent(in) :: frame
  character(len=*) :: vname
  character(len=*) :: stat_type 
  real(sp) :: stat
  integer  :: i,ndims,dimt,j

  !determine the var
  i = valindx(frame,vname) 

  !get dimension
  ndims = frame%fdata(i)%ndims

  !-gwc need to expand to integer type
  if(stat_type == 'max')then
    if(ndims==1)stat = maxval(frame%fdata(i)%f1)
    if(ndims==2)stat = maxval(frame%fdata(i)%f2)
    if(ndims==3)stat = maxval(frame%fdata(i)%f3)
    if(ndims==4)stat = maxval(frame%fdata(i)%f4)
  elseif(stat_type == 'min')then
    if(ndims==1)stat = minval(frame%fdata(i)%f1)
    if(ndims==2)stat = minval(frame%fdata(i)%f2)
    if(ndims==3)stat = minval(frame%fdata(i)%f3)
    if(ndims==4)stat = minval(frame%fdata(i)%f4)
  elseif(stat_type == 'mean')then
    dimt = 1 
    do j=1,ndims
      dimt = dimt*frame%fdata(i)%dims(j)
    end do
    if(ndims==1)stat = sum(frame%fdata(i)%f1)/float(dimt)
    if(ndims==2)stat = sum(frame%fdata(i)%f2)/float(dimt)
    if(ndims==3)stat = sum(frame%fdata(i)%f3)/float(dimt)
    if(ndims==4)stat = sum(frame%fdata(i)%f4)/float(dimt)
  else
    write(*,*)'fatal error in calc_stat'
    write(*,*)'do not know how to compute statistic: ',trim(stat_type)
    stop
  endif

end function calc_stat

function valindx(frame,vname) result(indx)
  implicit none
  type(dataframe), intent(in) :: frame
  character(len=*) :: vname
  integer :: indx
  integer :: i

  indx = 0
  do i=1,frame%nvars
    if(frame%fdata(i)%vname == trim(vname))then
      indx = i
    endif
  end do

  if(indx == 0)then
    write(*,*)'error in valindx'
    write(*,*)'frame does not contain variable: ',trim(vname)
    stop
  endif
   
end function valindx

!========================================================================
! interface data to outside world through pointer
!  - 1D, float
!========================================================================
subroutine get_forcing_f1(vname,p,iframe) 
  implicit none
  character(len=*) :: vname
  real(sp), pointer :: p(:)
  integer :: i,ndims,iframe
  integer, allocatable :: dims(:)

  !get variable index
  i = valindx(frame(iframe),vname)

  !get dimensions
  ndims = frame(iframe)%fdata(i)%ndims
  allocate(dims(ndims)) ; dims = 0
  dims = frame(iframe)%fdata(i)%dims

  if(ndims /= 1)then
    write(*,*)'error in get_forcing_f1'
    write(*,*)'number of dimensions of variable: ',trim(vname)
    write(*,*)'is: ',ndims,' cannot put this into a 1D array'
    stop
  endif

  !set pointer
  p => frame(iframe)%fdata(i)%f1

end subroutine get_forcing_f1

!========================================================================
! interface data to outside world through pointer
!  - 2D, float
!========================================================================
subroutine get_forcing_f2(vname,p,iframe) 
  implicit none
  character(len=*) :: vname
  real(sp), pointer :: p(:,:)
  integer :: i,ndims,iframe
  integer, allocatable :: dims(:)

  !get variable index
  i = valindx(frame(iframe),vname)

  !get dimensions
  ndims = frame(iframe)%fdata(i)%ndims
  allocate(dims(ndims)) ; dims = 0
  dims = frame(iframe)%fdata(i)%dims

  if(ndims /= 2)then
    write(*,*)'error in get_forcing_f2'
    write(*,*)'number of dimensions of variable: ',trim(vname)
    write(*,*)'is: ',ndims,' cannot put this into a 2D array'
    stop
  endif

  !set pointer
  p => frame(iframe)%fdata(i)%f2
end subroutine get_forcing_f2

function get_ncfid() result(fid)
  integer :: fid
  fid = ffile_id
end function get_ncfid

subroutine bracket(t,i1,i2) 
  use utilities
  real(sp), intent(in) :: t
  integer, intent(out) :: i1
  integer, intent(out) :: i2
  !------------------------
  integer :: i

  if(nframes <= 1)then
    write(*,*)'error in bracket: netcdf forcing file only has 1 frame'
    stop
  endif

  if(t < ftimes(1))then
    write(*,*)'error in bracket'
    write(*,*)'model time is not in the range of the netcdf dataset'
    write(*,*)'model time       : ',gettime(int(t))
    write(*,*)'netcdf begin time: ',gettime(int(ftimes(1)))
    write(*,*)'netcdf end   time: ',gettime(int(ftimes(nframes)))
    stop
  endif

  do i=1,nframes-1
    if(ftimes(i) <= t .and. t <= ftimes(i+1))then
      i1 = i
      i2 = i+1
    endif
  end do

end subroutine bracket
  
subroutine exchange_forcing
  implicit none

  !----------------------------------------
  integer  :: ndims,dtype,i
  integer  :: d1,d2,d3,d4


  !loop over vars, exchange frame
  do i=1,frame(now)%nvars
    dtype = frame(now)%fdata(i)%dtype
    ndims = frame(now)%fdata(i)%ndims
    if(dtype == flt_type)then
      if(ndims == 1) frame(now)%fdata(i)%f1 = frame(nex)%fdata(i)%f1
      if(ndims == 2) frame(now)%fdata(i)%f2 = frame(nex)%fdata(i)%f2
      if(ndims == 3) frame(now)%fdata(i)%f3 = frame(nex)%fdata(i)%f3
      if(ndims == 4) frame(now)%fdata(i)%f4 = frame(nex)%fdata(i)%f4
    else
      if(ndims == 1) frame(now)%fdata(i)%i1 = frame(nex)%fdata(i)%i1
      if(ndims == 2) frame(now)%fdata(i)%i2 = frame(nex)%fdata(i)%i2
      if(ndims == 3) frame(now)%fdata(i)%i3 = frame(nex)%fdata(i)%i3
      if(ndims == 4) frame(now)%fdata(i)%i4 = frame(nex)%fdata(i)%i4
    endif
  enddo


end subroutine exchange_forcing

!========================================================================
! open the forcing file  
!   - make sure it exists
!   - read time (convert to sec if nec.)
!========================================================================
subroutine open_forcing_file(ffiles_in,nfls,fbeg,fend) 
  use utilities
  implicit none
  integer, intent(in) :: nfls
  character(len=fstr) :: ffiles_in(nfls)
  integer :: nfsfid(nfls),nfsnrec(nfls)
  character(len=fstr) :: ffile_in
  real(sp), intent(out) :: fbeg,fend
  integer :: fid
  character(len=mstr) :: msg
  character(len=fstr) :: dname,tunits
  integer :: varid,i,ierr,n
  
  nframes=0 
  do n=1,nfls
  ffile_in=ffiles_in(n) !maybe it goes wrong
write(*,*)ffile_in
  msg = "error opening forcing file: "//trim(ffile_in)
  call ncdchk( nf90_open(trim(ffile_in),nf90_nowrite,fid),msg ) 
  !set file name and id for this module
  nfsfid(n)   = fid
  ffile(n) = ffile_in

  !inquire dataset info
  msg = "reading number of dimensions from: "//trim(ffile_in)
  call ncdchk(nf90_inquire(fid, nDim,nVar,nAtt,uDid,fNum) ,msg)
  !determine number of frames in the file (size of uDid)
  uD_extant = 1
  if(uDid /= 0)then
    call ncdchk(nf90_inquire_dimension(fid, uDid, dname, uD_extant ))
  endif
  nframes = uD_extant+nframes
  nfsnrec(n) = uD_extant
  end do

  !allocate space to hold times and read in 
  !nframes = uD_extant
   
  allocate(ftimes(nframes)) ; ftimes = 0.0
  allocate(fid_ft(nframes)) ; fid_ft = 0
  allocate(rec_ft(nframes)) ; rec_ft = 0
  allocate(fnr_ft(nframes)) ; fnr_ft = 0

  msg = "error reading time variable from netcdf file"
  uD_extant=0
  do n=1,nfls
  uD_extant=uD_extant+nfsnrec(n)
  ffile_id=nfsfid(n)
  call ncdchk( nf90_inq_varid(ffile_id,'time',varid),msg )
  allocate(tmptime(nfsnrec(n))) ; tmptime = 0.0
  call ncdchk(nf90_get_var(ffile_id, varid,tmptime ),msg)
  ftimes(uD_extant-nfsnrec(n)+1 : uD_extant  )=tmptime(1:nfsnrec(n))
  deallocate(tmptime)
  !read units on time to see if conversion from days to seconds is necessary
  if ( n == nfls ) then
!     ftimes = ftimes - 43*365.0  !gwc what is this?
!  if ( nf90_get_att(ffile_id, varid, 'units', tunits) == nf90_noerr)then
!    if(index(tunits,'day') /= 0) ftimes = ftimes*day_2_sec
!    write(*,*)'%%%%%% converting netcdf time units to seconds!!!!'
!  endif
  endif
  fid_ft(uD_extant-nfsnrec(n)+1 : uD_extant  ) = ffile_id
  rec_ft(uD_extant-nfsnrec(n)+1 : uD_extant  ) = (/(i,i=1,nfsnrec(n))/)
  fnr_ft(uD_extant-nfsnrec(n)+1 : uD_extant  ) = n
  enddo
   

  !set begin/end/deltaT from forcing
  !make sure time is sequential
  forcing_beg_time = ftimes(1)
  forcing_end_time = ftimes(nframes)

  if(nframes > 1)then
    do i=2,nframes
      if(ftimes(i)-ftimes(i-1) <= 0.0)then
        write(*,*)'netcdf time is not monotonically increasing'
      endif
    end do
  endif

  !set fbeg/fend to return val
  fbeg = forcing_beg_time
  fend = forcing_end_time
    
  call drawline('-')
  write(*,*)'Opened up model forcing file: ',ffiles_in
  write(*,*)'number of frames in file: ',nfsnrec
  write(*,*)'forcing begins at:',fbeg !gettime(int(fbeg)),fbeg
  write(*,*)'forcing ends   at:',fend !gettime(int(fend)),fend
  call drawline('-')

end subroutine open_forcing_file


End Module forcing
