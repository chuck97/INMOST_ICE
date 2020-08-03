########################################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make

LIBS    = -L$(ANILIB) $(LIBAFT) $(LIBFEM) $(LIBILU) $(LIBLMR) $(LIBMBA) \
	              $(LIBVIEW) $(LIBLAPACK) $(LIBBLAS) -lm


########################################################################
EXE         = $(ANIBIN)/multi_HybridMFEM.exe

OBJ         = main.o forlibfem.o 


########################################################################
all: 	    help

exe: 	    $(EXE) info

run:
	    cd $(ANIBIN) && ./multi_HybridMFEM.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps $(ANIBIN)/*.exe 

help:
	    @echo "make {exe|run|clean|help|gs}"
	    @echo " "
	    @echo "   exe - compile the code and link with libraries"
	    @echo " "
	    @echo "   run - run adaptive computation"    
	    @echo " "
	    @echo "   gs    - draw the mesh and solution "
	    @echo "   clean - clean the package"
	    @echo "   help  - print this message"
	    @echo " "

info:
	    @echo " "
	    @echo "Executable files are located in   bin"
	    @echo "Postscript figures are located in bin"
	    @echo " "

gs:
	    @echo " "
	    @echo "Postscript figures:"
	    @echo "   1. the initial mesh"
	    @echo "   2. isolines of solution interpolated from midedges to nodes"
	    @echo "   3. the final (adapted) mesh"
	    @echo "   4. isolines of solution "
	    @echo " "
	    cd $(ANIBIN); $(PS) mesh_initial.ps; $(PS) iso.ps; $(PS) mesh_final.ps; $(PS) his.ps



########################################################################
$(EXE): $(OBJ) 
	@echo ''
	@echo 'Creating executable bin/multi_HybridMFEM.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE) $(OBJ) $(LIBS)
	@echo ''


########################################################################
main.o      : main.f  
forlibfem.o : forlibfem.f


########################################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

