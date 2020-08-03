############################################################
ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

############################################################
LIBRCB  = $(ANILIB)/librcb2D-$(version).a

############################################################
LIBOBJS  = auxproc.o refine.o coarse.o error.o ../aniIO/io.o

############################################################
all:	    help

lib:	    $(LIBRCB) info

clean:
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.exe $(ANIBIN)/*.ps 

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib - create library librcb2D"
	    @echo " "
	    @echo "   clean - clean this directory"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo "Summary:"
	    @echo "Libraries   are located in lib/"
	    @echo " "

############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o 


############################################################
$(LIBRCB): $(LIBOBJS)
	   @echo ''
	   @rm -f $(LIBRCB)
	   @ar -r $(LIBRCB)  $(LIBOBJS)


############################################################
