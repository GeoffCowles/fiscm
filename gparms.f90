!=======================================================================
! Fiscm Global Parameters 
!
! Description
!   Defines global params (constants) for the fiscm code 
!    
! Comments: 
!    Originally from the OSCAR shallow water equation solver
!
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!
!=======================================================================

Module Gparms  

Implicit None

!single precision kind 
!integer, parameter :: sp = selected_real_kind(6 , 37)

!double precision kind
integer, parameter :: sp = selected_real_kind(15,307)

!double precision constants
real(sp), parameter :: an8th = 0.125_sp
real(sp), parameter :: a4th  = 0.250_sp
real(sp), parameter :: a3rd  = 1.0_sp/3.0_sp
real(sp), parameter :: ahalf = 0.500_sp
real(sp), parameter :: zero  = 0.000_sp
real(sp), parameter :: one   = 1.000_sp

!max arrays sizes
integer, parameter :: max_state_vars = 200

!----------------------------------------------------------------
!string length
!    fstr:    filename length
!    vstr:    variable name length
!    sstr:    short string length
!    tstr:    text string
!    cstr:    string variable standard length
!    mstr:    message length
!----------------------------------------------------------------
integer, parameter  :: fstr = 120  
integer, parameter  :: tstr = 120 
integer, parameter  :: vstr = 30    
integer, parameter  :: sstr = 15   
integer, parameter  :: cstr = 30 
integer, parameter  :: mstr = 120 
!----------------------------------------------------------------
!time
!     day_2_sec:    convert days to seconds 
!     sec_2_day:    convert seconds to days 
!----------------------------------------------------------------

real(sp), parameter :: day_2_sec = 86400. 
real(sp), parameter :: sec_2_day = one/day_2_sec 

!----------------------------------------------------------------
!statistical
!      rvar:  variance of a uniform random walk [-1,1]
!----------------------------------------------------------------
real(sp), parameter :: rvar = a3rd

!----------------------------------------------------------------
!trigonometric
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
integer, parameter :: NCDO_HEADER = 0
integer, parameter :: NCDO_ADD_STATES = 1
integer, parameter :: NCDO_OUTPUT = 2

!----------------------------------------------------------------
! status of individuals
!    DEAD: -3
!    SETTLED: -2
!    EXITED DOMAIN: -1
!    UNKNOWN: 0
!    ACTIVE: 1
!----------------------------------------------------------------
integer, parameter :: DEAD = -3
integer, parameter :: SETTLED = -2
integer, parameter :: EXITED = -1
integer, parameter :: UNKNOWN = 0
integer, parameter :: ACTIVE = 1

!----------------------------------------------------------------
! diffusion type
!----------------------------------------------------------------
integer, parameter :: HDIFF_NONE     = 0
integer, parameter :: HDIFF_CONSTANT = 1
integer, parameter :: HDIFF_VARIABLE = 2
integer, parameter :: VDIFF_NONE     = 0
integer, parameter :: VDIFF_VARIABLE = 1
integer, parameter :: VDIFF_SPLINED  = 2
integer, parameter :: VDIFF_BINNED   = 3
!-model setup----------------------------------------------------
integer, parameter ::  sz_cor     = 1  ! 0 - input s ;1 - input z
integer, parameter :: fix_dep     = 1  ! 0 - unfixed ;1 - fix(dep)
integer, parameter :: dvm_bio     = 0  ! 0 - nodvm   ;1 - dvm(bio)
real   , allocatable :: zpini(:),zptini(:)

!----------------------------------------------------------------
!----------------------------------------------------------------
! version & var character name 
!----------------------------------------------------------------
character(len=fstr) :: FISCM_VERSION= "fiscm1.0"

End Module Gparms
