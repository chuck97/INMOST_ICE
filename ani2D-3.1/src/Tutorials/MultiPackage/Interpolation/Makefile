############################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make

LIBS	= -L$(ANILIB) $(LIBAFT) $(LIBLMR) $(LIBFEM) $(LIBVIEW) \
                      $(LIBMBA) $(LIBLAPACK) $(LIBBLAS) -lm 


############################################################
EXE1        = $(ANIBIN)/multi_Interpolation.exe
EXE2        = $(ANIBIN)/multi_InterpolGrad.exe

OBJ1        = main.o     forlibfem.o
OBJ2        = mainGrad.o forlibfem.o


############################################################
all: 	    help

exe: 	    $(EXE1) $(EXE2) info

run-int:
	    cd $(ANIBIN) &&  ./multi_Interpolation.exe

run-grd:
	    cd $(ANIBIN) &&  ./multi_InterpolGrad.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps $(ANIBIN)/*.exe 

help:
	    @echo "make {exe|run-int|run-grd|clean|help|gs}"
	    @echo " "
	    @echo "     exe - compile the code and link with libraries "
	    @echo " "
	    @echo "     run-int - run adaptive computation"
	    @echo "     run-grd - run adaptive computation for gradient"
	    @echo " "
	    @echo "     gs     - show  mesh and solution initial and final"
	    @echo "     clean  - clean the package"
	    @echo "     help   - print this message"
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
	    cd $(ANIBIN); $(PS) mesh_initial.ps; $(PS) iso.ps; $(PS) mesh_final.ps; $(PS) his.ps


############################################################
$(EXE1): $(OBJ1) 
	@echo ''
	@echo 'Creating executable bin/multi_Interpolation.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE1) $(OBJ1) $(LIBS)
	@echo ''

$(EXE2): $(OBJ2) 
	@echo ''
	@echo 'Creating executable bin/multi_InterpolGrad.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE2) $(OBJ2) $(LIBS)
	@echo ''


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o  $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

