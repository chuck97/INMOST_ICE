############################################################
ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

############################################################
LIBOBJS  = solver_inb.o jv_inb.o bt_inb.o bcg_inb.o


############################################################
all:	    help

lib:	    $(LIBINB) info

clean:
	    @/bin/rm -f *.o 

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib  - create library libinb.a"
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
solver_inb.o:	solver_inb.f inbparam.h inbdflts.h
jv_inb.o:	jv_inb.f 
bt_inb.o:	bt_inb.f inbparam.h
bcg_inb.o:	bcg_inb.f 

############################################################
$(LIBINB): $(LIBOBJS)
	   @echo ''
	   @echo 'Creating library lib/aniINB-$(version).a'
	   @rm -f $(LIBINB)
	   @ar -r $(LIBINB)  $(LIBOBJS)



