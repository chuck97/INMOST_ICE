############################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make

LIBS	= -L$(ANILIB) $(LIBAFT) $(LIBFEM) $(LIBLMR) $(LIBMBA) $(LIBVIEW) \
	              $(LIBLU) $(LIBLAPACK) $(LIBBLAS) -lm


############################################################
EXE         = $(ANIBIN)/multi_ConDifAdap.exe

OBJ         = main.o forlibfem.o 


############################################################
all: 	    help

exe: 	    $(EXE) info

run:
	    cd $(ANIBIN) &&  ./multi_ConDifAdap.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps $(ANIBIN)/*.exe 

help:
	    @echo "make {exe|run|clean|help|gs}"
	    @echo " "
	    @echo "   exe - compile the code and link with libraries "
	    @echo " "
	    @echo "   run - run adaptive computation"
	    @echo " "
	    @echo "   gs    - show  mesh and solution initial and final"
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
	    @echo " "
	    @echo "Postscript figures:"
	    @echo "   1. initial mesh"
	    @echo "   2. isolines of the corresponding solution"
	    @echo "   3. final (adapted) mesh"
	    @echo "   4. isolines of the corresponding solution"
	    @echo " "
	    cd $(ANIBIN); $(PS) mesh_initial.ps; $(PS) iso_ini.ps; 
	    cd $(ANIBIN); $(PS) mesh_final.ps;   $(PS) iso_fin.ps



############################################################
$(EXE): $(OBJ) 
	@echo ''
	@echo 'Creating executable bin/ConDifAdap.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE) $(OBJ) $(LIBS)
	@echo ''


############################################################
main.o      : main.f  
forlibmba.o : forlibmba.f
forlibfem.o : forlibfem.f


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

