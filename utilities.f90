!=======================================================================
! Fiscm Utilities 
!
! Description
!   Utilities for fiscm 
!    
! !REVISION HISTORY:                   
!  Original author(s): G. Cowles 
!
! Fortran90 Random Number generators from the Web
!
!=======================================================================
module utilities 
!uses
  use gparms
  implicit none

! default all is private.
  private
 
! public member functions:
  public gettime 
  public drawline 
  public isintriangle 
  public cfcheck 
  public get_unique_strings
  public ran1 
  public ran_from_range 
  public unitrand 
  public normal 
!  public spline
  public normal_circle
  public random_square

interface normal    
  module procedure normal 
  module procedure normal_zeromean 
end interface

contains

!------------------------------------------------------------------------------!
!   Return a Time String Days:Hours:Minutes:Seconds from Number of Seconds     !
!     input:                                                                   !
!       integer:  insecs [number of seconds]                                   !
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

  !---------------------------------------------
  !draw a line 72 characters long repeating c 
  !optional: dump to unit iunit
  !---------------------------------------------
  subroutine drawline(c,iunit)
    character(len=1) :: c
    integer, intent(in), optional :: iunit
    character(len=72) :: line
    integer :: i
    line(1:1) = "!"
    do i=2,72
      line(i:i) = c
    end do
    if(present(iunit))then
      write(iunit,*)line
    else
      write(*,'(A72)')line
    endif
  end subroutine drawline
    
  !------------------------------------------------------------------------------
  !  determine if point (x0,y0) is in triangle defined by nodes (xt(3),yt(3))    |
  !  using algorithm used for scene rendering in computer graphics               |
  !  algorithm works well unless particle happens to lie in a line parallel      |
  !  to the edge of a triangle.                                                  |
  !  This can cause problems if you use a regular grid, say for idealized        |
  !  modelling and you happen to see particles right on edges or parallel to     |
  !  edges.                                                                      |
  !------------------------------------------------------------------------------
   logical function isintriangle(i,x0,y0,xt,yt)
   implicit none
   integer,  intent(in) :: i
   real(sp), intent(in) :: x0,y0
   real(sp), intent(in) :: xt(3),yt(3)
   !----------------------------------
   real(sp) :: f1,f2,f3
   real(sp) :: x1(2)
   real(sp) :: x2(2)
   real(sp) :: x3(2)
   real(sp) :: p(2)
   !----------------------------------
   !revised by Xinyou Lin
    isintriangle = .true.

   f1 = (y0-yt(1))*(xt(2)-xt(1)) - (x0-xt(1))*(yt(2)-yt(1))
   f1 = f1*((yt(3)-yt(1))*(xt(2)-xt(1)) - (xt(3)-xt(1))*(yt(2)-yt(1)))  
   f2 = (y0-yt(3))*(xt(1)-xt(3)) - (x0-xt(3))*(yt(1)-yt(3))
   f2 = f2*((yt(2)-yt(3))*(xt(1)-xt(3)) - (xt(2)-xt(3))*(yt(1)-yt(3)))

   f3 = (y0-yt(2))*(xt(3)-xt(2)) - (x0-xt(2))*(yt(3)-yt(2))
   f3 =f3*((yt(1)-yt(2))*(xt(3)-xt(2)) - (xt(1)-xt(2))*(yt(3)-yt(2)))

   if(f1 <0.0_sp .or. f2 <0.0_sp .or.f3 <0.0_sp ) isintriangle = .false.
   return
   end function isintriangle

  !-----------------------------------------------
  ! runtime errors  - netcdf
  !-----------------------------------------------
  subroutine cfcheck(status)
      use netcdf
      integer, intent ( in) :: status
  
      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop
      end if
  end subroutine cfcheck

  !-----------------------------------------------
  ! search list of strings and return unique 
  !-----------------------------------------------
  subroutine get_unique_strings(ns,s,n_uniq)
    integer, intent(in)  :: ns
    character(len=*)     :: s(ns)
    integer, intent(out) :: n_uniq
    !------------------------------
    character(len=fstr)  :: stmp(ns)
    integer :: uniq(ns)
    integer :: i,j,ii

    if(ns < 1)then 
      write(*,*)'error in get_unique_strings'
      write(*,*)'number of strings to check == 0'
      stop
    endif

    uniq = 1 
    do i=1,ns    
      stmp(i) = s(i)
    end do
    
    do i=2,ns    
      do j=1,i-1
        if( s(j) == s(i)) uniq(i) = 0 
      end do
    end do

    !reinitialize s
    do i=1,ns    
      s(i) = ""
    end do

    !count unique
    n_uniq = sum(uniq)  

    !transfer uniq 
    ii = 0
    do i=1,ns    
      if(uniq(i) == 1)then
        ii = ii + 1
        s(ii) = stmp(i)
      endif
    end do

  end subroutine get_unique_strings

  !-----------------------------------------------
  ! return a random number, precision "sp" using  
  ! fortran90 intrinsic "random_number"
  !-----------------------------------------------
  real(sp) function ran1()  
    implicit none 
    real(sp) x 
    call random_number(x) 
    ran1=x 
  end function ran1

  !-----------------------------------------------
  ! return a random number in specified interval 
  !-----------------------------------------------
  function ran_from_range(fmin,fmax)  
    implicit none 
    real(sp) ran_from_range 
    real(sp) fmin,fmax 
    ran_from_range=(fmax - fmin) * ran1() + fmin 
  end function ran_from_range

  !-----------------------------------------------
  ! return unit random number (-1 or 1) 
  !-----------------------------------------------
  real(sp) function unitrand()  
    implicit none 
    real(sp) :: tmp
    tmp = ran1()
    unitrand = sign(1.0_sp,ran1()-0.5_sp)
  end function unitrand 

  !-----------------------------------------------
  ! return random number from normal distribution 
  ! with mean ->  mean and standard dev -> sigma
  ! from?
  !-----------------------------------------------
  real(sp) function normal(mean,sigma) 
    implicit none 
    real(sp), intent(in) ::  mean
    real(sp), intent(in) ::  sigma 
    
    normal=normal_zeromean()*sigma+mean !normal_zeromean() is distribution of standard normal
    return 
  end function normal


  !-----------------------------------------------
  ! return random number from normal distribution 
  ! with mean = 0.0 and dev = 1. 
  !-----------------------------------------------
  function normal_zeromean() result(mynormal)
    implicit none 
    real(sp) :: r1,r2
    real(sp) :: mynormal
    real(sp),parameter :: PI = 3.14159265358979d0
    r1 = ran1()
    r2 = ran1()
    mynormal = sqrt(DBLE(-2.)*log(r1)) * cos(DBLE(2.)*PI*r2)

    return 
  end function normal_zeromean

  !-----------------------------------------------
  ! fit a cubic splint (zero tension) to data 
  ! from numerical recipes
  ! in:  
  !   n:  dimension of data
  !   x:  independent variable 
  !   y:  dependent variable
  ! yp1:  boundary condition at i=1
  ! ypn:  boundary condition at i=n
  !  ys:  spline
  !-----------------------------------------------
!   subroutine spline(n,x,y,yp1,ypn,ysp)
! 
!     integer,  intent(in ) :: n
!     real(sp), intent(in ) :: x(n)
!     real(sp), intent(in ) :: y(n)
!     real(sp), intent(in ) :: yp1
!     real(sp), intent(in ) :: ypn
!     real(sp), intent(out) :: ysp(n)  
!     !------------------------------
!     integer, parameter :: nmax = 50
!     integer  :: i,k
!     real(sp) :: p,qn,sig,un,u(nmax)
!   
!     if (yp1.gt..99e30) then
!         ysp(1)=0.
!         u(1)=0.
!      else
!         ysp(1)=-0.5
!         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
!      endif
!      do i=2,n-1
!         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
!         p=sig*ysp(i-1)+2.
!         ysp(i)=(sig-1.)/p
!         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))--sig*u(i-1))/p
!      end do
!      if (ypn.gt..99e30) then
!         qn=0.
!         un=0.
!      else
!         qn=0.5
!         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
!      endif
!      ysp(n)=(un-qn*u(n-1))/(qn*ysp(n-1)+1.)
!      do k=n-1,1,-1
!         ysp(k)=ysp(k)*ysp(k+1)+u(k)
!      end do
!      return
!   end subroutine spline
  
  !-----------------------------------------------
  ! return x,y pseudo-normally distributed in 
  ! circle centered at (xc,yc), of radius r
  ! r will extend to two standard devs
  !-----------------------------------------------
  subroutine normal_circle(n,xc,yc,rad,x,y) 
    implicit none 
    integer,  intent(in)  :: n
    real(sp), intent(in)  :: xc
    real(sp), intent(in)  :: yc
    real(sp), intent(in)  :: rad 
    real(sp), intent(out) :: x(n) 
    real(sp), intent(out) :: y(n) 
    real(sp) :: theta,rval
    integer  :: i 

    do i=1,n
      theta = ran_from_range(0.0_sp,2*pi)
      rval  = normal(0.0_sp,rad/2.)
      x(i) = rval*cos(theta)
      y(i) = rval*sin(theta)
    end do
      
  end subroutine normal_circle 

  !-----------------------------------------------
  ! random distribution of [n] particles in a square 
  ! [xmin,xmax,ymin,ymax]
  !-----------------------------------------------
  subroutine random_square(n,xmin,xmax,ymin,ymax,x,y)
    implicit none
    integer,  intent(in)  :: n
    real(sp), intent(in)  :: xmin
    real(sp), intent(in)  :: xmax
    real(sp), intent(in)  :: ymin
    real(sp), intent(in)  :: ymax
    real(sp), intent(out) :: x(n)
    real(sp), intent(out) :: y(n)
    real(sp) :: theta,rval
    integer  :: i

    do i=1,n
      x(i) = ran_from_range(xmin,xmax)
      y(i) = ran_from_range(ymin,ymax)
    end do
    

  end subroutine random_square 


end module utilities 
