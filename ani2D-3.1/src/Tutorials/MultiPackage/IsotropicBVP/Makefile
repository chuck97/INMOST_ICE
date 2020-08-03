############################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make

LIBS = -L$(ANILIB) $(LIBAFT) $(LIBLMR) $(LIBMBA) $(LIBFEM) \
                   $(LIBLU)  $(LIBILU) $(LIBVIEW) \
	           $(LIBLAPACK) $(LIBBLAS) -lm



############################################################
EXE1 = $(ANIBIN)/multi_IsoAdapHessian.exe
EXE2 = $(ANIBIN)/multi_IsoAdapPosteriori.exe

OBJ1 = main_hes.o forlibfem_hes.o forlibaft.o forlibmba.o 
OBJ2 = main_est.o forlibfem_est.o forlibaft.o forlibmba.o 


############################################################
all: 	    help

exe: 	    $(EXE1) $(EXE2) info

run-hes:
	    cd $(ANIBIN) && ./multi_IsoAdapHessian.exe

run-est:
	    cd $(ANIBIN) && ./multi_IsoAdapPosteriori.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps $(ANIBIN)/*.exe 

help:
	    @echo "make {exe|run-hes|run-est|clean|help|gs}"
	    @echo " "
	    @echo "   exe - compile the code and link with libraries "
	    @echo " "
	    @echo "   run-hes - run adaptive computation with Hessian recovery"
	    @echo " "
	    @echo "   run-est - run adaptive computation with aposteriori error estimate"
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
	    cd $(ANIBIN); $(PS) mesh_initial.ps; $(PS) iso.ps; $(PS) mesh_final.ps; $(PS) his.ps


############################################################
$(EXE1): $(OBJ1) 
	@echo ''
	@echo 'Creating executable bin/multi_IsoAdapHessian.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE1) $(OBJ1) $(LIBS)
	@echo ''

$(EXE2): $(OBJ2) 
	@echo ''
	@echo 'Creating executable bin/multi_IsoAdapPosteriori.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE2) $(OBJ2) $(LIBS)
	@echo ''


############################################################
main_hes.o      : main_hes.f  
forlibfem_hes.o : forlibfem_hes.f
main_est.o      : main_est.f  
forlibfem_est.o : forlibfem_est.f
forlibaft.o     : forlibaft.c
forlibmba.o     : forlibmba.f


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

