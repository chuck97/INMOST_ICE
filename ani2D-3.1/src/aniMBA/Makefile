ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

############################################################
# YOU DO NOT NEED TO CHANGE THE LINES BELOW
############################################################
OBJCORE     = mba_analytic.o mba_nodal.o mba_fixshape.o \
              ani2.o forlibmba.o calcrv.o \
              makM.o mesh_audit.o makQ_linear.o module_hessian.o list_new.o makS.o \
              prjcrv.o lintrp2D.o dsort.o \
              nlnfnc.o minim.o \
              colapse_edgeF1.o colapse_edgeF2.o insrt_point.o swap_edge.o \
	      move_point.o split_triangle.o \
	      auxSP.o auxSF.o tangled.o check.o status.o control.o \
	      triangle.o edge.o \
              update.o maps.o utils.o \
	      refine.o smoothing.o untangle.o \
              ZZ.o delaunay.o datacopy.o \
	      statistics.o error.o time.o isnan.o random.o \
	      ../aniIO/io.o

OBJIO       = load_mesh.o loadM_other.o saveM.o saveM_other.o

OBJDBG      = debug.o

OBJLIBMBA   = $(OBJCORE) $(OBJIO) $(OBJDBG)


############################################################
all: 	    help

lib:	    $(LIBMBA) info

clean: 
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.*ps $(ANIBIN)/*.exe 

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib   - library libmba2D-*.a"
	    @echo "   clean - clean the package"
	    @echo "   help  - print this message"
	    @echo " "

info:
	    @echo " "
	    @echo "Library files are located in ../../lib"
	    @echo " "


############################################################
$(LIBMBA): $(OBJLIBMBA)
	@echo ''
	@rm -f $(LIBMBA)
	@ar -r $(LIBMBA)  $(OBJLIBMBA)


############################################################
mba_analytic.o  : lintrp.fd  status.fd
mba_nodal.o     : lintrp.fd  status.fd
mba_fixshape.o  : lintrp.fd  status.fd

ani2.o          : makS.fd    colors.fd  status.fd  operat.fd  magic.fd
metric.o        :            magic.fd
metric_aee.o    : makS.fd    magic.fd

swap_edge.o     : makS.fd    colors.fd  status.fd  operat.fd
clpsF1.o        : makS.fd    colors.fd  status.fd  operat.fd
clpsF2.o        : makS.fd    colors.fd  status.fd  operat.fd
move_point.o    : makS.fd    colors.fd  status.fd  magic.fd
update.o 	: makS.fd
insrt_point.o 	: makS.fd    colors.fd  status.fd  operat.fd
split_triangle.o: makS.fd    colors.fd  status.fd

minim.o  	: status.fd  makS.fd
nlnfnc.o        : status.fd
makM.o   	: makS.fd    colors.fd  status.fd
mesh_audit.o   	:            colors.fd  status.fd
module_hessian.o: magic.fd 
makS.o   	: makS.fd    colors.fd  status.fd  operat.fd
lintrp2D.o 	: lintrp.fd
list.o		: list.fd

status.o        : status.fd
control.o       : status.fd

ZZ.o            : makS.fd
smoothing.o     : makS.fd

debug.o   	: colors.fd  status.fd


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC) -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

