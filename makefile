
#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/sh
            DEF_FLAGS     = -P -C -traditional 
            EXEC          = fiscm 
#==========================================================================
#  BEG USER DEFINITION SECTION
#    REQUIRED:
#       1.) Fortran90 Compiler
#       2.) NetCDF 3.x libraries with Fortran90 interface 
#           it would be best if they were built using the same compiler
#       3.) In some compilers, iargc/getarg are not standard
#           for example, Absoft 9.x requires a link to library U77 (-lU77)
#           if you are using gfortran/intel you should not have to worry
#           about this.
#==========================================================================
#         ==== REQUIRED
#          IOLIBS       =  /usr/local/netcdf/gfortran/lib/libnetcdf.la
#          IOINCS       =  -I/usr/local/netcdf/gfortran/include
#          CPPFLAGS = $(DEF_FLAGS)
#          CPP      = /usr/bin/cpp
#          FC       = gfortran #-funderscoring # -fno-second-underscore

#          ==== OPTIONAL
#          FLINK    = /bin/sh /usr/local/netcdf/gfortran/bin/libtool  --mode=link gfortran
#          DEBFLGS  =
#          OPT      =

include make.inc
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


F95FILES=    gparms.f90 utilities.f90 pvar.f90 igroup.f90 bio.f90\
		drivers.f90 output.f90 forcing.f90 fvcom_driver.f90 fiscm.f90


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
	tar cvf fiscm.tar *.f90  makefile *.txt *.nml  makedepends make.inc

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
