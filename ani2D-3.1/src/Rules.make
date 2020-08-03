############################################################
# DEFINE compilers, graphic viewer and make
############################################################
F77     = gfortran
CC      = gcc
PS      = gv -scalebase=2    #name of a PostScript viewer
MAKE    = make

FLINKER = $(F77)
CLINKER = $(CC)

FFLAGS  = -O2 
CFLAGS  = -O2 -fPIC
LDFLAGS = 

LIBSYS  = -lgfortran         #for linking C and Fortran

#LDFLAGS = -lg77libs              #for PGI
#LIBSYS  = -lifport -lifcore      #for Intel



############################################################
ANIAFT  = $(ANIHOME)/src/aniAFT
ANIPRJ  = $(ANIHOME)/src/aniPRJ
ANIFEM  = $(ANIHOME)/src/aniFEM
ANILMR  = $(ANIHOME)/src/aniLMR
ANILU   = $(ANIHOME)/src/aniLU 
ANIILU  = $(ANIHOME)/src/aniILU 
ANIINB  = $(ANIHOME)/src/aniINB 
ANIMBA  = $(ANIHOME)/src/aniMBA
ANIRCB  = $(ANIHOME)/src/aniRCB
ANIVIEW = $(ANIHOME)/src/aniVIEW
ANIC2F  = $(ANIHOME)/src/aniC2F

ANIBIN  = $(ANIHOME)/bin
ANILIB  = $(ANIHOME)/lib
ANIDAT  = $(ANIHOME)/data

ANIBLAS   = $(ANIHOME)/src/blas
ANILAPACK = $(ANIHOME)/src/lapack


############################################################
ANISTOKES = $(ANIHOME)/src/Tutorials/MultiPackage/Stokes


############################################################
version = 3.1


############################################################
LIBAFT  = $(ANILIB)/libaft2D-$(version).a
LIBPRJ  = $(ANILIB)/libprj2D-$(version).a 
LIBFEM  = $(ANILIB)/libfem2D-$(version).a 
LIBLMR  = $(ANILIB)/liblmr2D-$(version).a 
LIBMBA  = $(ANILIB)/libmba2D-$(version).a 
LIBLU   = $(ANILIB)/liblu-5.1.a 
LIBILU  = $(ANILIB)/libilu-$(version).a 
LIBINB  = $(ANILIB)/libinb-$(version).a 
LIBRCB  = $(ANILIB)/librcb2D-$(version).a 
LIBVIEW = $(ANILIB)/libview2D-$(version).a 
LIBC2F  = $(ANILIB)/libc2f2D-$(version).a 
LIBUTILS= $(ANILIB)/libutils2D-$(version).a 

#LIBBLAS   = -lblas
#LIBLAPACK = -llapack
LIBBLAS   = $(ANILIB)/libblas-3.2.a
LIBLAPACK = $(ANILIB)/liblapack-3.2.a $(ANILIB)/liblapack_ext-3.2.a 


############################################################
INCLUDE = -I$(ANIMBA) -I$(ANIFEM) -I$(ANIC2F)



