#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/sh
            DEF_FLAGS     = -P -C -traditional 
            EXEC          = fiscm 
#--------------------------------------------------------------------------
#        NETCDF OUTPUT      NEED NETCDF 3.x 
#                           NEED FORTRAN90 INTERFACE
#                           IOLIBS - PATH TO libnetcdf.a
#                           IOINCS - PATH TO netcdf.mod
#--------------------------------------------------------------------------

#          --- following for oscar/rudy ----
	  IOLIBS       =  /usr/local/netcdf/gfortran/lib/libnetcdf.la
          IOINCS       =  -I/usr/local/netcdf/gfortran/include

#          --- following for minke ----
#             IOLIBS       =  -L/hosts/salmon01/data00/medm/lib  -lnetcdf
#             IOINCS       =  -I/hosts/salmon01/data00/medm/include


#--------------------------------------------------------------------------
#  APPLE OS X/Absoft (G5)
#--------------------------------------------------------------------------
#         CPPFLAGS = $(DEF_FLAGS)  
#         ABSOFTLIB = /Applications/Absoft/lib/f77_oldnames.o
#         CPP      = /usr/bin/cpp
#         FC       = f90 -YEXT_PFX=_ 
#         DEBFLGS  =  # -qcheck -C -g
#         OPT      =
#--------------------------------------------------------------------------
# APPLE / gfortran 
#--------------------------------------------------------------------------
         CPPFLAGS = $(DEF_FLAGS)
         CPP      = /usr/bin/cpp
         FC       = gfortran #-funderscoring # -fno-second-underscore
         FLINK    = /bin/sh /usr/local/netcdf/gfortran/bin/libtool  --mode=link gfortran
         DEBFLGS  = 
         OPT      = 
#--------------------------------------------------------------------------
# LINUX / ifort 
#--------------------------------------------------------------------------
#        CPPFLAGS = $(DEF_FLAGS)
#         CPP      = /usr/bin/cpp
#         FC       = ifort #-funderscoring # -fno-second-underscore
#         FLINK    = $(FC) 
#         DEBFLGS  =
#         OPT      =

#==========================================================================
#  END USER DEFINITION SECTION
#==========================================================================

         FFLAGS = $(DEBFLGS) $(OPT) 
         MDEPFLAGS = --cpp --fext=f90 --file=-
         RANLIB = ranlib

#--------------------------------------------------------------------------
#  Libraries / Include Files          
#--------------------------------------------------------------------------

	LIBS  =     $(IOLIBS)
	INCS  =     $(IOINCS) 

#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90 

.f90.o:
	$(FC) -c $(FFLAGS) $(FIXEDFLAGS) $(INCS) $(INCLDIR) $*.f90  

#--------------------------------------------------------------------------
#  FISCM Source Code.
#--------------------------------------------------------------------------


F95FILES=    calendar.f90 gparms.f90 mod_pvar.f90 mod_igroup.f90 bio.f90 adv_diff.f90\
			output.f90 forcing.f90 fiscm.f90


 SRCS = $(F95FILES) 

 OBJS = $(SRCS:.f90=.o)

#--------------------------------------------------------------------------
#  Linking Directives               
#--------------------------------------------------------------------------

$(EXEC):	$(OBJS)
		$(FLINK) $(FFLAGS) $(INCS) $(LDFLAGS) -o $(EXEC) $(OBJS) $(LIBS)
#		$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXEC) $(OBJS) $(LIBS)

#--------------------------------------------------------------------------
#  Target to create dependecies.
#--------------------------------------------------------------------------

depend:
		makedepf90  $(SRCS) >> makedepends


#--------------------------------------------------------------------------
#  Tar Up Code                           
#--------------------------------------------------------------------------

tarfile:
	tar cvf fiscm.tar *.f90  makefile input.txt  makedepends 

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod fiscm 

#--------------------------------------------------------------------------
#  Common rules for all Makefiles - do not edit.
#--------------------------------------------------------------------------

emptyrule::

#--------------------------------------------------------------------------
#  Empty rules for directories that do not have SUBDIRS - do not edit.
#--------------------------------------------------------------------------

install::
	@echo "install in $(CURRENT_DIR) done"

install.man::
	@echo "install.man in $(CURRENT_DIR) done"

Makefiles::

includes::
include ./makedepends
