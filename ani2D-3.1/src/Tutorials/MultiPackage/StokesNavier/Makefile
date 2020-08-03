############################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make

LIBS	= -L$(ANILIB) $(LIBAFT) $(LIBLU) $(LIBINB) $(LIBLMR) \
	              $(LIBPRJ) $(LIBFEM) $(LIBMBA) \
	              $(LIBVIEW) $(LIBLAPACK) $(LIBBLAS) -lm 


############################################################
EXE         = $(ANIBIN)/multi_StokesNavier.exe

OBJ         =  main.o forlibfem.o forlibinb.o forlibaft.o


############################################################
all: 	    help

exe: 	    $(EXE) info

run:
	    cd $(ANIBIN) && ./multi_StokesNavier.exe

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
	    @echo "  2-3. isolines of initial velocity components (u_x and u_y)"
	    @echo "   4. velocity streamlines"
	    @echo "   5. the final (adapted) mesh"
	    @echo "  6-7. isolines of final velocity components (u_x and u_y)"
	    @echo "   8. velocity streamlines"
	    @echo " "
 
	    @cd $(ANIBIN); $(PS) mesh_initial.ps; $(PS) velx_ini.ps; $(PS) vely_ini.ps; \
                                                  $(PS) streamlines_ini.ps; \
	                   $(PS) mesh_final.ps;   $(PS) velx_fin.ps; $(PS) vely_fin.ps; \
                                                  $(PS) streamlines_fin.ps


############################################################
$(EXE): $(OBJ) 
	@echo ''
	@echo 'Creating executable bin/multi_StokesNavier.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE) $(OBJ) $(LIBS) 
	@echo ''


############################################################
main.o      : main.f  
forlibfem.o : forlibfem.f
forlibinb.o : forlibinb.f


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

