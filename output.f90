!=======================================================================
! Fiscm NetCDF Output Routines 
!
! Description:
!   Routines to output fiscm data to NetCDF 
!
! Comments:  Requires fortran90 NetCDF 3.x libraries
!    
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!
!=======================================================================
module output_routines

use gparms
use mod_igroup
use netcdf
implicit none


integer, private  ::  time_did
integer, private  ::  nlag_did
integer, private  ::  time_vid
integer, private  :: dynm2d(2)
integer, private  :: dynm1d(1)

public

interface add_cdfglobs
  module procedure add_cdfglobs_int
  module procedure add_cdfglobs_flt
  module procedure add_cdfglobs_str
end interface

contains

!-------------------------------------------
! output driver
!-------------------------------------------
subroutine cdf_out(ng,g,itnum,time,otype)
  integer, intent(in) :: ng
  type(igroup),intent(inout), dimension(ng) :: g
  integer , intent(in) :: itnum
  real(sp), intent(in) :: time 
  integer,  intent(in) :: OTYPE
  integer :: n
 
  select case(otype)
  !write the netcdf header file
  case(NCDO_HEADER)             
    do n=1,ng
      call write_header(g(n))    
    end do
  !add the user-defined state variables to the output header
  case(NCDO_ADD_STATES)
    do n=1,ng
      call addtocdf_states(g(n)) !
    end do
  !output the state variables 
  case(NCDO_OUTPUT)
    do n=1,ng
      if(mod(itnum-1,g(n)%intvl_out)==0)call output_group(g(n),time)
    end do
  end select

end subroutine cdf_out

!-------------------------------------------
! write netcdf header
!-------------------------------------------
subroutine write_header(g)
  type(igroup), intent(inout) :: g
  character(len=fstr)  :: fname
  !character(len=1)     :: num
  integer              :: ierr,fid
 
  g%frame_out = 1
  !write(num,'(I1)') g%id
  

  !construct the name (eventually use some sim string/directory loc)
  !fname = 'fiscm_'//num//'.nc'
  !fname = g%fname_out//'_'//num//'.nc'
  fname = trim(g%fname_out)
  !g%fname_out = fname

  !open the file - write access
  call cfcheck( nf90_create(trim(fname),nf90_clobber,fid) ) ;  g%fid_out = fid
  
  !global attributes
  call cfcheck( nf90_put_att(fid,nf90_global,"title", "FVCOM PARTICLE TRACKING")) 
  call cfcheck( nf90_put_att(fid,nf90_global,"Conventions", "CF-1.0")) 
  call cfcheck( nf90_put_att(fid,nf90_global,"source"      ,"ParticleFVCOM_3.0"))  !trim(FISCM_VERSION)) )
  call cfcheck( nf90_put_att(fid,nf90_global,"group"     ,g%id) )
  call cfcheck( nf90_put_att(fid,nf90_global,"DT_bio"    ,g%DT_bio) )

  !define the time dimension   
  call cfcheck(nf90_def_dim(fid,"time",nf90_unlimited, time_did) )
  dynm1d = (/time_did/)

  !time variable
  call cfcheck( nf90_def_var(fid,"time",nf90_float,dynm1d, time_vid) )
  call cfcheck( nf90_put_att(fid, time_vid,"long_name","time") )
  call cfcheck( nf90_put_att(fid, time_vid,"units","days since 0.0") )
  call cfcheck( nf90_put_att(fid, time_vid,"time_zone","none") )

  !lag dimension
  call cfcheck(nf90_def_dim(fid,"nlag",g%Tnind       , nlag_did) )
  dynm2d = (/nlag_did,time_did/)

  !close the file
  call cfcheck( nf90_close(g%fid_out) )
end subroutine write_header
  
!-----------------------------------------------
! add state variable definitions to netcdf file
!-----------------------------------------------
subroutine addtocdf_states(g)
type(igroup), intent(inout) :: g
integer              :: ierr

!open the file
ierr = nf90_open(g%fname_out,nf90_write,g%fid_out)
if(ierr /= nf90_noerr)then
  write(*,*)'error opening', trim(g%fname_out)
  write(*,*)trim(nf90_strerror(ierr))
endif
call cfcheck( nf90_redef(g%fid_out) )

!write the state variables slated for output
call write_cdf_header_vars(g%state,g%fid_out,dynm2d)

!close the file
call cfcheck( nf90_close(g%fid_out) )

end subroutine addtocdf_states

!-----------------------------------------------
! add additional parameters to netcdf global
! parameter list -- FLOATS
!-----------------------------------------------
subroutine add_cdfglobs_flt(g,desc,val)
type(igroup), intent(inout) :: g
character(len=*)     :: desc
real(sp), intent(in) :: val
integer              :: ierr

!open the file
ierr = nf90_open(g%fname_out,nf90_write,g%fid_out)
if(ierr /= nf90_noerr)then
  write(*,*)'error opening', trim(g%fname_out)
  write(*,*)trim(nf90_strerror(ierr))
endif
call cfcheck( nf90_redef(g%fid_out) )

!write variable to global
  call cfcheck( nf90_put_att(g%fid_out,nf90_global,trim(desc),val) )

!close the file
call cfcheck( nf90_close(g%fid_out) )

end subroutine add_cdfglobs_flt

!-----------------------------------------------
! add additional parameters to netcdf global
! parameter list -- INTEGERS
!-----------------------------------------------
subroutine add_cdfglobs_int(g,desc,val)
type(igroup), intent(inout) :: g
character(len=*)     :: desc
integer , intent(in) :: val
integer              :: ierr

!open the file
ierr = nf90_open(g%fname_out,nf90_write,g%fid_out)
if(ierr /= nf90_noerr)then
  write(*,*)'error opening', trim(g%fname_out)
  write(*,*)trim(nf90_strerror(ierr))
endif
call cfcheck( nf90_redef(g%fid_out) )

!write variable to global
  call cfcheck( nf90_put_att(g%fid_out,nf90_global,trim(desc),val) )

!close the file
call cfcheck( nf90_close(g%fid_out) )

end subroutine add_cdfglobs_int

!-----------------------------------------------
! add additional parameters to netcdf global
! parameter list -- STRINGS
!-----------------------------------------------
subroutine add_cdfglobs_str(g,desc,val)
type(igroup), intent(inout) :: g
character(len=*)     :: desc
character(len=*)     :: val
integer              :: ierr

!open the file
ierr = nf90_open(g%fname_out,nf90_write,g%fid_out)
if(ierr /= nf90_noerr)then
  write(*,*)'error opening', trim(g%fname_out)
  write(*,*)trim(nf90_strerror(ierr))
endif
call cfcheck( nf90_redef(g%fid_out) )

!write variable to global
  call cfcheck( nf90_put_att(g%fid_out,nf90_global,trim(desc),val) )

!close the file
call cfcheck( nf90_close(g%fid_out) )

end subroutine add_cdfglobs_str

!-----------------------------------------------
! output time dependent vars (states and time) 
! to netcdf files
!-----------------------------------------------
subroutine output_group(g,time)
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: time
  integer              :: ierr
  integer              :: dims(1)

  !set frame
  dims(1) = g%frame_out

  !debug
  !write(*,*)'dumping: ',g%id,' file: ',trim(g%fname_out),' frame: ',g%frame_out

  !open the file
  ierr = nf90_open(g%fname_out,nf90_write,g%fid_out)
  if(ierr /= nf90_noerr)then
    write(*,*)'error opening', trim(g%fname_out)
    write(*,*)trim(nf90_strerror(ierr))
  endif

  !dump time
  call cfcheck( nf90_put_var(g%fid_out,time_vid,time,START=dims) )

  !dump state variable data
  call write_cdf_data(g%state,g%fid_out,g%frame_out)

  !increment frame counter
  g%frame_out = g%frame_out + 1

  !close up
  call cfcheck( nf90_close(g%fid_out) )

end subroutine output_group

!-----------------------------------------------
! runtime errors 
!-----------------------------------------------
subroutine cfcheck(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 
    end if
end subroutine cfcheck

end module output_routines
