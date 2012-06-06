#          --- following for oscar/rudy ----
	  #IOLIBS       =  -L/hosts/mao/usr/medm/install/netcdf/3.6.3/em64t/lib -lnetcdf
          #IOINCS       =  -I/hosts/mao/usr/medm/install/netcdf/3.6.3/em64t/include
	  IOLIBS       =  -L/opt/netcdf3/ifort_gcc64/lib -lnetcdf 
          IOINCS       =  -I/opt/netcdf3/ifort_gcc64/include


#          --- following for minke ----
#             colon=:
#             empty=
#             dashI= $(empty) -I
#	     dashL= $(empty) -L
#             IOLIBS       =  -L$(subst $(colon),$(dashL),$(LIBPATH)) -lnetcdf 
#             IOINCS       =  -I$(subst $(colon),$(dashI),$(INCLUDEPATH)) 

#--------------------------------------------------------------------------
# APPLE / gfortran 
#--------------------------------------------------------------------------
#         CPPFLAGS = $(DEF_FLAGS)
#         CPP      = /usr/bin/cpp
#         FC       = gfortran #-funderscoring # -fno-second-underscore
#         FLINK    = /bin/sh /usr/local/netcdf/gfortran/bin/libtool  --mode=link gfortran
#         DEBFLGS  = 
#         OPT      = -O3
#--------------------------------------------------------------------------
# LINUX / ifort 
#--------------------------------------------------------------------------
        CPPFLAGS = $(DEF_FLAGS)
         CPP      = /usr/bin/cpp
         FC       = ifort #-funderscoring # -fno-second-underscore
         FLINK    = ifort ##/bin/sh /usr/local/netcdf/gfortran/bin/libtool  --mode=link gfortran
         DEBFLGS  = 
         OPT      = -O3 -g
