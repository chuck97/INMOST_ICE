############################################################
ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

############################################################
LIBOBJS  = bcg.o cg.o gmres.o matvecCSR.o iluoo.o ilu0.o dsort.o \
	   ../aniIO/io.o


############################################################
all:	    help

lib:	    $(LIBILU) info

clean:
	    @/bin/rm -f *.o 

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib  - create library libilu.a"
	    @echo " "
	    @echo "   clean - clean the package"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo " "
	    @echo "Summary:"
	    @echo "Libraries   are located in lib/"
	    @echo " "


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o 


############################################################
$(LIBILU): $(LIBOBJS)
	   @echo ''
	   @echo 'Creating library lib/aniILU-$(version).a'
	   @rm -f $(LIBILU)
	   @ar -r $(LIBILU)  $(LIBOBJS)



