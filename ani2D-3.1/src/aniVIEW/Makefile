############################################################
ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

############################################################
OBJVIEW  = isoline.o mesh.o mesh_full.o streamlines.o matrix.o demo.o


############################################################
all: 	    help

lib:	    $(LIBVIEW) info

clean: 
	    @/bin/rm -f *.o $(ANIBIN)/*.exe *~ $(ANIBIN)/*.ps 

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib - create library libview2D.a"
	    @echo " "
	    @echo "   clean - clean this directory"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo " "
	    @echo "Library files are located in lib"
	    @echo " "


############################################################
$(LIBVIEW): $(OBJVIEW)
	@echo ''
	@rm -f $(LIBVIEW)
	@ar -r $(LIBVIEW)  $(OBJVIEW)


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

