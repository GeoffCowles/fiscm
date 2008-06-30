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
  integer :: ngroups

  ngroups = 2
  allocate(igroups(ngroups))
  igroups(1) = group_(21)
  igroups(2) = group_(15)

  call init_bio(igroups(1),100)
  call advance_bio(igroups(1))


  deallocate(igroups)


End Program fiscm
