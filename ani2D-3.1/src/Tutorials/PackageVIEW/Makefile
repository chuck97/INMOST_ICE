############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make

############################################################
LIBS    = -L$(ANILIB) $(LIBVIEW) -lm


############################################################
EXEMESH   = $(ANIBIN)/aniVIEW_mesh.exe
EXEMATRIX = $(ANIBIN)/aniVIEW_matrix.exe

OBJMESH   = main.o
OBJMATRIX = main_matrix.o


############################################################
all: 	    help

exe:        $(EXEMESH) $(EXEMATRIX) 

run:        
	    @cd $(ANIBIN) && ./aniVIEW_mesh.exe
	    @cd $(ANIBIN) && ./aniVIEW_matrix.exe

clean: 
	    @/bin/rm -f *.o $(ANIBIN)/*.exe *~ $(ANIBIN)/*.ps 

help:
	    @echo "make {exe|run|gs|clean|help}"
	    @echo " "
	    @echo "   exe - compile the regression tests"
	    @echo "   run - run the programs"
	    @echo "   gs  - draw a mesh, function isolines, and matrix"
	    @echo " "
	    @echo "   clean - clean this directory"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo " "
	    @echo "Executabel files are located in bin"
	    @echo " "

gs:
	    @echo " "
	    @echo "Postscript figures:"
	    @echo "   1. the mesh"
	    @echo "   2. isolines of solution"
	    @echo "   3. matrix"
	    @cd $(ANIBIN); $(PS) mesh.ps; $(PS) isolines.ps; $(PS) matrix.ps


############################################################
$(EXEMESH): $(OBJMESH)
	@echo ''
	@echo 'Creating executable bin/aniVIEW_mesh.exe'
	@$(FLINKER) $(LDFLAGS) -o $(EXEMESH) $(OBJMESH) $(LIBVIEW) $(LIBS)

$(EXEMATRIX): $(OBJMATRIX)
	@echo ''
	@echo 'Creating executable bin/aniVIEW_matrix.exe'
	@$(FLINKER) $(LDFLAGS) -o $(EXEMATRIX) $(OBJMATRIX) $(LIBVIEW) $(LIBS)


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

