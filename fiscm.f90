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
! Major Todos:
!    1.) migrate to D.S. netcdf libraries
!    2.) migrate to D.S. time type
!    3.) openmp doloops to multithread heavier sections 
!        parallel random number generator?
!=======================================================================

Program fiscm
  use gparms 
  use mod_igroup
  use mod_AD
  use bio
  use output_routines
  integer :: i

  type(igroup), allocatable :: igroups(:)
  real(sp) :: beg_time,end_time
  real(sp) :: deltaT
  integer  :: sim_direction
  integer  :: ngroups

  !---------------------------------------------------
  ! initialize groups and data
  !---------------------------------------------------
  ngroups = 2
  allocate(igroups(ngroups))
  igroups(1) = group_(21,1)
  igroups(2) = group_(15,2)
  beg_time   = 0.0
  end_time   = 24*3600*20.
  deltaT     = 3600

  do i=1,ngroups
    call print_group_summary(igroups(i))
  end do

  call init_bio(igroups(1),5)
  
  !---------------------------------------------------
  ! main loop over time 
  !---------------------------------------------------
  t = beg_time
  do while (t <= end_time)
	
	call adv_diff(ngroups,igroups,t)
	
    !call bio(igroups,t)

    call output(ngroups,igroups,t)

    t = t + deltaT
  end do

  !---------------------------------------------------
  ! cleanup data
  !---------------------------------------------------
  deallocate(igroups)


End Program fiscm
