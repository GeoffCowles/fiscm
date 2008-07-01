!=======================================================================
! Fiscm Global Parameters
! Copyright:    2008(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Comments:     FISCM Global Parameters
!=======================================================================

Module Gparms  

Implicit None

!single precision kind 
integer, parameter :: sp = selected_real_kind(6 , 37)

!double precision kind
!integer, parameter :: sp = selected_real_kind(15,307)

!double precision constants
real(sp), parameter :: an8th = 0.125_sp
real(sp), parameter :: a4th  = 0.250_sp
real(sp), parameter :: a3rd  = 1.0_sp/3.0_sp
real(sp), parameter :: ahalf = 0.500_sp
real(sp), parameter :: zero  = 0.000_sp
real(sp), parameter :: one   = 1.000_sp

!----------------------------------------------------------------
!string length
!    fstr:    filename length
!    vstr:    variable name length
!    sstr:    short string length
!    tstr:    text string
!    cstr:    string variable standard length
!----------------------------------------------------------------
integer, parameter  :: fstr = 120  
integer, parameter  :: tstr = 120 
integer, parameter  :: vstr = 30    
integer, parameter  :: sstr = 15   
integer, parameter  :: cstr = 30 

!----------------------------------------------------------------
!trigonemetric
!      pi:    pi
!     d2r:    convert degrees to radians
!     r2d:    convert radians to degrees
!----------------------------------------------------------------

real(sp), parameter :: pi  = 3.14159265358979312 
real(sp), parameter :: d2r = pi/180.0_sp
real(sp), parameter :: r2d = 180.0_sp/pi

!----------------------------------------------------------------
!oceanic parameters
!    gacc        :  gravitational acceleration   [ms^-2]
!    omega_earth :  earths rotation rate         [s^-1]
!----------------------------------------------------------------
real(sp) :: gacc  = 9.8016_sp !default
real(sp), parameter :: omega_earth = 7.292e-5

!----------------------------------------------------------------
!limiting parameters
!    CFL_MAX:    max practical CFL
!----------------------------------------------------------------
real(sp), parameter :: CFL_MAX = 5.0_sp

!----------------------------------------------------------------
!large and small numbers
!    hugenum = largest float
!    tinynum = smallest float
!----------------------------------------------------------------
real(sp), parameter :: hugenum = huge(1.0_sp)
real(sp), parameter :: tinynum = tiny(1.0_sp)

!----------------------------------------------------------------
! control netcdf output of variables
!    NETCDF_YES:  yes, output
!    NETCDF_NO:   no, do not output
!----------------------------------------------------------------
integer, parameter :: NETCDF_YES = 1
integer, parameter :: NETCDF_NO  = 0

End Module Gparms
