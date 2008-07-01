!=======================================================================
! Fiscm Main 
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

Program fiscm
  use gparms 
  use mod_igroup
  use bio
  type(igroup), allocatable :: igroups(:)
  integer :: ngroups,i
  real(sp) :: beg_time,end_time
  real(sp) :: deltaT

  ngroups = 2
  allocate(igroups(ngroups))
  igroups(1) = group_(21)
  igroups(2) = group_(15)
  beg_time = 0.0
  end_time   = 24*3600*20.
  deltaT     = 3600

  do i=1,ngroups
    call print_group_summary(igroups(i))
  end do

  call init_bio(igroups(1),5)
  
  t = beg_time
  do while (t <= end_time)
    call advance_bio(igroups(1),t)
    t = t + deltaT
  end do

  

  deallocate(igroups)


End Program fiscm
