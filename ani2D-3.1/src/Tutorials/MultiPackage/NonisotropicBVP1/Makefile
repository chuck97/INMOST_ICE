############################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make

LIBS	= -L$(ANILIB) $(LIBAFT) $(LIBFEM) $(LIBLU)$(LIBLMR) $(LIBMBA) \
	              $(LIBVIEW) $(LIBLAPACK) $(LIBBLAS) -lm 


############################################################
# YOU DO NOT NEED TO CHANGE THE LINES BELOW
############################################################
EXE         = $(ANIBIN)/multi_AnisoAdap.exe

OBJ         = main.o forlibaft.o forlibmba.o forlibfem.o 


############################################################
all: 	    help

exe: 	    $(EXE) info

run:
	    cd $(ANIBIN) && ./multi_AnisoAdap.exe

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
	    @cd $(ANIBIN); $(PS) mesh_initial.ps; $(PS) iso_ini.ps;\
	                   $(PS) mesh_final.ps;   $(PS) iso_fin.ps


############################################################
$(EXE): $(OBJ) 
	@echo ''
	@echo 'Creating executable $(EXE)' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE) $(OBJ) $(LIBS)
	@echo ''


############################################################
main.o      : main.f  
forlibaft.o : forlibaft.c
forlibmba.o : forlibmba.f
forlibfem.o : forlibfem.f


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

