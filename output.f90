!=======================================================================
! Fiscm NetCDF Output Routines
! Copyright:    2008(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Comments:  Requires fortran90 NetCDF 3x libraries
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


contains
!-------------------------------------------
! output driver
!-------------------------------------------
subroutine output(ng,g,time)
  integer, intent(in) :: ng
  type(igroup),intent(inout), dimension(ng) :: g
  real(sp), intent(in) :: time
  integer :: n

   do n=1,ng
     if( (time-g(n)%tlast_out) < g(n)%intvl_out .or. time < g(n)%start_out) cycle
     if(g(n)%frame_out == 0) call write_header(g(n),time)
     call output_group(g(n),time)
   end do
end subroutine output

!-------------------------------------------
! write netcdf header
!-------------------------------------------
subroutine write_header(g,time)
  type(igroup), intent(inout) :: g
  real(sp), intent(in) :: time
  character(len=fstr)  :: fname
  character(len=1)     :: num
  integer              :: ierr,fid
 
  g%frame_out = 1
  write(num,'(I1)') g%id
  

  !construct the name (eventually use some sim string/directory loc)
  fname = 'fiscm_'//num//'.nc'
  g%fname_out = fname

  !open the file - write access
  call cfcheck( nf90_create(trim(fname),nf90_clobber,fid) ) ;  g%fid_out = fid
  
  !global attributes
  call cfcheck( nf90_put_att(fid,nf90_global,"code"      ,"FISCM") )
  call cfcheck( nf90_put_att(fid,nf90_global,"group"     ,g%id) )

  !write the dimensions (time, number of individuals)
  call cfcheck(nf90_def_dim(fid,"nlag",g%Tnind       , nlag_did) )
  call cfcheck(nf90_def_dim(fid,"time",nf90_unlimited, time_did) )

  !define the dimensions   
  dynm2d = (/nlag_did,time_did/)
  dynm1d = (/time_did/)

  !time variable
  call cfcheck( nf90_def_var(fid,"time",nf90_float,dynm1d, time_vid) )
  call cfcheck( nf90_put_att(fid, time_vid,"long_name","time") )
  call cfcheck( nf90_put_att(fid, time_vid,"units","days") )

  !write the state variables slated for output
  call write_cdf_header_vars(g%state,fid,dynm2d)

  !close the file
  call cfcheck( nf90_close(g%fid_out) )

end subroutine write_header
  
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
  call cfcheck( nf90_put_var(g%fid_out,time_vid,time/(3600*24),START=dims) )

  !dump state variable data
  call write_cdf_data(g%state,g%fid_out,g%frame_out)

  !increment frame counter
  g%frame_out = g%frame_out + 1

  !close up
  call cfcheck( nf90_close(g%fid_out) )

end subroutine output_group

subroutine cfcheck(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 
    end if
end subroutine cfcheck

end module output_routines
