############################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make

LIBS	= -L$(ANILIB) $(LIBAFT) $(LIBFEM) $(LIBLU) $(LIBLMR) $(LIBMBA) \
	              $(LIBVIEW) $(LIBLAPACK) $(LIBBLAS) -lm 


############################################################
EXE         = $(ANIBIN)/multi_Stokes.exe

OBJ         =  main.o forlibaft.o forlibmba.o forlibfem.o 


############################################################
all: 	    help

exe: 	    $(EXE) info

run:
	    cd $(ANIBIN) && ./multi_Stokes.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps $(ANIBIN)/*.exe 

help:
	    @echo "make {exe|run|clean|help|gs}"
	    @echo " "
	    @echo "   exe - compile the code and link with libraries"
	    @echo " "
	    @echo "   run - run adaptive computation"
	    @echo " "
	    @echo "   gs    - draw the mesh and solution"
	    @echo "   clean - clean the package"
	    @echo "   help  - print help message"
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
	    @echo "   2. streamlines"
	    @echo "   3. isolines of the pressure"
	    @echo "   4. the final (adapted) mesh"
	    @echo "   5. streamlines"
	    @echo "   6. isolines of the pressure"
	    @echo " "
 
	    @cd $(ANIBIN); $(PS) mesh_initial.ps; $(PS) streamlines_ini.ps; $(PS) pressure_ini.ps; \
	                   $(PS) mesh_final.ps;   $(PS) streamlines_fin.ps; $(PS) pressure_fin.ps


############################################################
$(EXE): $(OBJ) 
	@echo ''
	@echo 'Creating executable bin/Stokes.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE) $(OBJ) $(LIBS) 
	@echo ''


############################################################
main.o      : main.f  
forlibaft.o : forlibaft.f
forlibfem.o : forlibfem.f
forlibmba.o : forlibmba.f


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

