########################################################################
ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

LIBS	= -L$(ANILIB) -llapack-3.0 -lblas-3.0 -lm


########################################################################
LIBOBJS = template.o assemble.o fem2Dtri.o fem2Dedge.o fem2Dsub.o \
	  opIDEN.o opGRAD.o opDIV.o opCURL.o opDUDX.o \
	  opTRACE.o opDUDN.o \
          bc.o algebra.o forlibfem.o isys.o dof.o \
          utils.o error.o fem2Derr.o \
	  DG.o


########################################################################
all:	    help

lib:	    $(LIBFEM) info

clean:
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.exe $(ANIBIN)/CSRsystem

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib - create library libfem2D-2.X.a"
	    @echo " "
	    @echo "   clean - clean the package"
	    @echo "   help  - print this message"
	    @echo " "

info:
	    @echo " "
	    @echo "Summary:"
	    @echo "Executables are located in bin/"
	    @echo "Libraries   are located in lib/"
	    @echo " "


########################################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o -I$(ANIMBA)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o


########################################################################
$(LIBFEM): $(LIBOBJS)
	   @echo ''
	   @echo 'Creating library $(LIBFEM)'
	   @rm -f $(LIBFEM)
	   @ar -r $(LIBFEM)  $(LIBOBJS)


########################################################################
fem2Dtri.o : fem2Dtri.fd
fem2Dedge.o: fem2Dtri.fd
fem2Dsub.o : fem2Dtri.fd
fem2Derr.o : fem2Dtri.fd

opIDEN.o   : fem2Dtri.fd
opGRAD.o   : fem2Dtri.fd
opDIV.o    : fem2Dtri.fd
opCURL.o   : fem2Dtri.fd
opDUDX.o   : fem2Dtri.fd

opTRACE.o  : fem2Dtri.fd
opDUDN.o   : fem2Dtri.fd

bc.o       : fem2Dtri.fd

template.o : fem2Dtri.fd assemble.fd
assemble.o : fem2Dtri.fd assemble.fd

dof.o      : fem2Dtri.fd
isys.o     : assemble.fd

DG.o       : fem2Dtri.fd assemble.fd

