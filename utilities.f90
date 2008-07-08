!=======================================================================
! Fiscm Utilities 
! Copyright:    2008(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
!=======================================================================
module utilities 
use gparms
implicit none


contains
!------------------------------------------------------------------------------!
!   Return a Time String Days:Hours:Minutes:Seconds from Number of Seconds     !
!------------------------------------------------------------------------------!

function gettime(insecs) result(instring)

  implicit none
  character(len=13)   :: instring
  integer, intent(in) :: insecs 
  character(len=4)    :: s0
  character(len=2)    :: s1,s2,s3
  integer :: dtcp,htcp,mtcp,stcp
 
  dtcp = insecs/(3600*24)
  htcp = mod(insecs,(3600*24))/3600
  mtcp = mod(insecs,(3600))/60
  stcp = insecs - (dtcp*3600*24 + htcp*3600 + mtcp*60)

  if(dtcp < 10000.)then
  if(dtcp >= 1000.)then
    write(s0,"(i4)")int(dtcp)
  else if(dtcp >= 100.)then
    write(s0,"('0',i3)")int(dtcp)
  else if(dtcp >= 10)then
    write(s0,"('00',i2)")int(dtcp)
  else
    write(s0,"('000',i1)")int(dtcp)
  end if
  if(htcp >= 10.)then
    write(s1,"(i2)")int(htcp)
  else
    write(s1,"('0',i1)")int(htcp)
  end if
  if(mtcp >= 10.)then
    write(s2,"(i2)")int(mtcp)
  else
    write(s2,"('0',i1)")int(mtcp)
  end if
  if(stcp >= 10.)then
    write(s3,"(i2)")int(stcp)
  else
    write(s3,"('0',i1)")int(stcp)
  end if
  instring = s0//":"//s1//":"//s2//":"//s3
  else
  instring = "> 1000 days"
  end if

  return
  end function gettime

end module utilities
