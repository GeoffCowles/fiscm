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
     !call output_group(g(n),time)
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
  integer              :: ierr
  
  
  g%frame_out = 1
  write(num,'(I1)') g%id

  !construct the name (eventually use some sim string/directory loc)
  fname = 'fiscm_'//num//'.nc'
  write(*,*)'writing header for group: ',g%id,' file: ',trim(fname)

  !open the file - write access
  call cfcheck( nf90_create(trim(fname),nf90_clobber,g%fid_out) )
  
  !global attributes
  call cfcheck( nf90_put_att(g%fid_out,nf90_global,"code"      ,"FISCM") )

  !write the dimensions 

  !write critical scalar variables

  !write the state variables slated for output

  !close the file
  call cfcheck( nf90_close(g%fid_out) )

end subroutine write_header
  
subroutine cfcheck(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 
    end if
end subroutine cfcheck

end module output_routines
