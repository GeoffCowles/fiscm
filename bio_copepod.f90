!-------------------------------------------------
!-----Calanus model, version 1.0------------------
!------Jess McNally, WHOI--------------------------
MODULE bio
     
Use gparms
Use mod_igroup

IMPLICIT NONE 
  
! REAL BEL_BETA ! Beta for Belehradek equation, common for all calanus

!------STATE VARIABLES declared in copepod.nml-----

! ---STAGE 
! Stage number keeps track of stage, gender and maturity. 

! INTEGER iStage(N_IND)
! Keeps track of stage as an integer

! REAL pasd(N_IND) 
!'Proportional Accumulative Stage Days', to keep track of progress in development

! Stages = 1 Egg
!        = 2 N1
!        = 3 N2
!        = 4 N3
!        = 5 N4
!        = 6 N5
!        = 7 N6
!        = 8 C1
!        = 9 C2
!        =10 C3
!        =11 C4
!        =12 C5
!        =13 Male
!        =14 Immature Female
!        =15 Mature Female 


!parameters
integer, parameter :: nstages = 15
real(sp) :: BEL_BETA = -2.05

!=====Passed from species specific vars module====
! For Belehradek equation Duration = a*(temperature+alpha)^beta
! will just declare them here now, for testing. These are values for C. fin
 
real(sp) :: BEL_ALPHA 
integer :: BEL_A(13)
character*100  INI_BIOF      
Namelist /NML_COPEPOD/  &
     & BEL_ALPHA,       &
     & BEL_A,           &   
     & INI_BIOF
contains

subroutine init_bio(g,Nind_start)
       use output_routines
       use utilities 
       type(igroup), intent(inout) :: g
       integer,      intent(in)    :: Nind_start
       integer :: iunit, ios, i, np,npoint
       logical :: fexist
       real(sp) , pointer :: tspawn(:)
       real(sp) , pointer :: x(:)
       real(sp) , pointer :: y(:)
       real(sp) , pointer :: s(:)
       real(sp) , pointer :: z(:)
       integer  , pointer :: istatus(:)

       real(sp)  tmp
           
      !set current problem size
       if(Nind_start < 0)then
          write(*,*)'must specify more than 0 individuals in init_bio'
       stop
       else if(Nind_start > g%Tnind)then
          write(*,*)'number of individuals in init_bio (',Nind_start,') exceeds problem size',g%Tnind
       stop
       endif
       g%Nind = Nind_start
  !set problem size
  np = g%Nind

  !read the namelist (if any) 
  if(g%paramfile /= "NONE")then
  inquire(file=trim(g%paramfile),exist=fexist)
  if(.not.fexist)then
    write(*,*)'fatal error: namelist file: ',trim(g%paramfile),' does not exist, stopping...'
    stop 
  endif
  write(*,*)'initializing namelist using paramfile: ',g%paramfile
  open(unit=iunit,file=trim(g%paramfile),form='formatted')
  read(unit=iunit,nml=nml_copepod,iostat=ios)
  if(ios /= 0)then
    write(*,*)'fatal error: could not read copepod namelist from: ',trim(g%paramfile)
    stop
  endif
  endif ! file /= NONE
  close(iunit)
 write(*, nml=nml_copepod)

 !------------------------------------
 ! set the spawning time
 !------------------------------------
  call get_state('tspawn',g,tspawn)
  tspawn(1:np) = 0.0*day_2_sec  
  nullify(tspawn)
 
  !-----------------------------------
  ! 2D,3D -> initialize x,y 
  !-----------------------------------
  if(g%space_dim > 1)then
    call get_state('x',g,x)
    call get_state('y',g,y)
    call get_state('status',g,istatus)

!    call random_square(np,-2500._sp,2500._sp,-2500._sp,2500._sp,x,y)
!    Revised by Xinyou Lin
  
  inquire(file=trim(INI_BIOF),exist=fexist)
  if(.not.fexist)then
    write(*,*)INI_BIOF ,' does not exist, stopping...'
    stop
  endif
  write(*,*)'initializing x y z using file: ',INI_BIOF
  open(unit=iunit,file=trim(INI_BIOF),form='formatted')
  read(iunit,*)npoint
  if(npoint /= np)then
   write(*,*)'The number of points are wrong in',INI_BIOF
   stop
  endif
  do i=1,np
  read(iunit,*) npoint,x(i),y(i)
  enddo 


  close(iunit)
   istatus=1;
    nullify(x)
    nullify(y)
    nullify(istatus)

  endif
  !-----------------------------------
  ! 3D -> initialize s-coordinate
  !-----------------------------------
  if(g%space_dim > 2)then
    call get_state('s',g,s)
    call get_state('z',g,z)
!    do i=1,g%Nind
!      s(i) = -float(i-1)/float(g%Nind-1)
!    end do

  open(unit=iunit,file=trim(INI_BIOF),form='formatted')
  read(iunit,*)npoint

  if(npoint /= np)then
   write(*,*)'The number of points are wrong in',INI_BIOF
   stop
  endif
  if(sz_cor == 1) then
  allocate(zpini(np))
  allocate(zptini(np))
  endif 
  do i=1,np
  read(iunit,*) npoint,tmp,tmp,zpini(i),zptini(i)
  if(sz_cor == 0)then
  s(i) = zpini(i)
  elseif(sz_cor == 1)then
  z(i) = zpini(i)
  endif
! hours to seconds
  enddo
 
  zptini=zptini*3600.0
  close(iunit)
    nullify(s)
    nullify(z)

  endif


  !add parameters to netcdf file header 
  call add_cdfglobs(g,"info","some kind of info")

end subroutine init_bio 

!----------------------------------------------------
! advance the biology (this routine is called from the main loop at
!                      a time interval of DT_bio )
!------------------------------------------------------------------
subroutine advance_bio(g,mtime)  
  type(igroup), intent(inout) :: g
  real(sp),     intent(in   ) :: mtime
  real(sp),     pointer :: PASD(:)
  real(sp),     pointer :: T(:)
  integer ,     pointer :: stage(:)
  integer ,     pointer :: istatus(:)
  integer ,     pointer :: diapause(:)
  integer               :: i,N
  real(sp)              :: deltaT, stg_duri, stg_durj, stg_durdiff

 
 !construct pointers to access and modify state variables for the group
  call get_state('PASD',g,PASD)
  call get_state('T',g,T)
  call get_state('status',g,istatus)
  call get_state('stage',g,stage)
  call get_state('diapause',g,diapause) 

 
  !set problem size
  N = g%nind
  deltaT = g%DT_bio*sec_2_day !time step in days 

!  following interface is cleaner but compiler cannot construct
!  generic interface if arguments are the same even if the
!  return value is of different type
!  PASD   => get_state('PASD',g)
!  istatus => get_state('status',g)
!  stage  => get_state('stage',g)


! Open file for output to text file:
! open(unit=21, file='caloutput.txt', form='formatted')
  
   
!============================================================

!------------MAIN LOOP-------------------

       sind_loop: DO i=1, N

!----check for diapause----
! needs to be filled in
 
             if(istatus(i) /= ACTIVE) cycle

!increment development

            develop: SELECT CASE (stage(i))

!Development from egg to C5
             CASE (1:12)
                 stg_duri = BEL_A(stage(i))*((T(i)+BEL_ALPHA)**BEL_BETA) 
                 stg_durj = BEL_A(stage(i)+1)*((T(i)+BEL_ALPHA)**BEL_BETA)
                 stg_durdiff = stg_durj - stg_duri
                 PASD(i) = PASD(i) + (deltaT/stg_durdiff)   

!If stage has surpassed 13, need to make half of those go to 14 instead (half male, half female)
!Simple solution is that 'odd individuals' stay male, and 'even' individuals become female
!Can maybe put a random generator in here later?
                  tofemale: IF (stage(i)>13.AND.MOD(i,2)==0) THEN
                           stage(i) = stage(i)+1
                  ENDIF tofemale

!Males, eventually can put some function in to model their death rate.						
               CASE (13)
                   CYCLE

!Females, immature and mature
               CASE (14:)
                   CYCLE
   
                     END SELECT develop
    


!truncate PASD to nearest whole number, to keep track of stage as int
                       stage = AINT(PASD(i))

   


       ENDDO sind_loop 
! nullify pointers 
  nullify(PASD)
  nullify(T)
  nullify(istatus)
  nullify(stage)
  nullify(diapause)

End Subroutine advance_bio


End Module  bio 
