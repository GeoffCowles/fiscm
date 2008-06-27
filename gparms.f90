!=======================================================================
! Oscar Global Parameter Class
! Copyright:    2005(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Authors:      G. Cowles 
!               School for Marine Science and Technology     
!               University of Massachusetts-Dartmouth
!
! Comments:     Oscar Global Parameters
!=======================================================================

Module Gparms  

Implicit None

!single precision kind (used for pV3)
integer, parameter :: sp = selected_real_kind(6 , 37)

!double precision kind
integer, parameter :: dp = selected_real_kind(15,307)

!double precision constants
real(dp), parameter :: an8th = 0.125_dp
real(dp), parameter :: a4th  = 0.250_dp
real(dp), parameter :: a3rd  = 1.0_dp/3.0_dp
real(dp), parameter :: ahalf = 0.500_dp
real(dp), parameter :: zero  = 0.000_dp
real(dp), parameter :: one   = 1.000_dp

!----------------------------------------------------------------
!string length
!    fstr:    filename length
!    vstr:    variable name length
!    sstr:    subroutine name length
!    tstr:    text string
!    cstr:    string variable standard length
!----------------------------------------------------------------
integer, parameter  :: fstr = 120  
integer, parameter  :: tstr = 120 
integer, parameter  :: vstr = 30    
integer, parameter  :: sstr = 30   
integer, parameter  :: cstr = 30 

!----------------------------------------------------------------
!trigonemetric
!      pi:    pi
!     d2r:    convert degrees to radians
!     r2d:    convert radians to degrees
!----------------------------------------------------------------

real(dp), parameter :: pi  = 3.14159265358979312 
real(dp), parameter :: d2r = pi/180.0_dp
real(dp), parameter :: r2d = 180.0_dp/pi

!----------------------------------------------------------------
!oceanic parameters
!    gacc        :  gravitational acceleration   [ms^-2]
!    omega_earth :  earths rotation rate         [s^-1]
!----------------------------------------------------------------
real(dp) :: gacc  = 9.8016_dp !default
real(dp), parameter :: omega_earth = 7.292e-5

!----------------------------------------------------------------
!limiting parameters
!    CFL_MAX:    max practical CFL
!----------------------------------------------------------------
real(dp), parameter :: CFL_MAX = 5.0_dp

!----------------------------------------------------------------
!large and small numbers
!    hugenum = largest float
!    tinynum = smallest float
!----------------------------------------------------------------
real(dp), parameter :: hugenum = huge(1.0_dp)
real(dp), parameter :: tinynum = tiny(1.0_dp)


End Module Gparms
