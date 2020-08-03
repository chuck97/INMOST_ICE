############################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make

LIBS = -L$(ANILIB) $(LIBAFT) $(LIBILU) $(LIBFEM) $(LIBMBA) \
	           $(LIBVIEW) $(LIBLAPACK) $(LIBBLAS) -lm


############################################################
EXE1 = $(ANIBIN)/multi_DG.exe

OBJ1 = main.o forlibaft.o forlibfem.o 


############################################################
all: 	    help

exe: 	    $(EXE1) info

run:
	    cd $(ANIBIN) &&  ./multi_DG.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps $(ANIBIN)/*.exe 

help:
	    @echo "make {exe|run|clean|help|gs}"
	    @echo " "
	    @echo "   exe - compile the code and link with libraries "
	    @echo " "
	    @echo "   run - run adaptive computation with Hessian recovery"
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
	    @echo "   1. mesh"
	    @echo "   2. solution isolines"
	    @echo "   3. matrix sparsity structure"
	    @echo " "
	    cd $(ANIBIN); $(PS) mesh.ps; $(PS)iso.ps
	    cd $(ANIBIN); $(PS) matrix.ps


############################################################
$(EXE1): $(OBJ1) 
	@echo ''
	@echo 'Creating executable bin/multi_DG.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE1) $(OBJ1) $(LIBS)
	@echo ''


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

