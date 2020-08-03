############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make

############################################################
LIBS	= -L$(ANILIB) $(LIBFEM) $(LIBVIEW) $(LIBLAPACK) $(LIBBLAS) -lm


############################################################
EXEFEMTRI = $(ANIBIN)/aniFEM_triangle.exe 
EXEFEMBC  = $(ANIBIN)/aniFEM_BC.exe
EXEFEMTEM = $(ANIBIN)/aniFEM_template.exe

FEMTRIOBJS = mainTriangle.o
FEMBCOBJS  = mainBC.o
FEMTEMOBJS = mainTemplate.o


############################################################
all:	    help

exe:        $(EXEFEMTRI) $(EXEFEMBC) $(EXEFEMTEM)

clean:
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.exe $(ANIBIN)/CSRsystem

help:
	    @echo " "
	    @echo "make {exe|run-triangle|run-stokes|run-bc|clean|help}"
	    @echo " "
	    @echo "   exe - compile three executables"
	    @echo " "
	    @echo "   run-triangle - generate a few elemental matrices"
	    @echo "   run-bc       - generate matrix for diffusion problem with various boundary conditions"
	    @echo "   run-stokes   - generate matrix for the Stokes problem"
	    @echo " "
	    @echo "   gs - draws sparcity structure of the generated matrix"
	    @echo " "
	    @echo "   clean - remove temporary files"
	    @echo "   help  - print this message"
	    @echo " "

info:
	    @echo " "
	    @echo "Summary:"
	    @echo "Executables are located in bin/"
	    @echo "Libraries   are located in lib/"
	    @echo " "

run-triangle:
	    cd $(ANIBIN) && ./aniFEM_triangle.exe

run-bc:
	    cd $(ANIBIN) && ./aniFEM_BC.exe

run-stokes:
	    cd $(ANIBIN) && ./aniFEM_template.exe

gs:	    
	    @cd $(ANIBIN); $(PS) matrix_fin.ps


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o


############################################################
$(EXEFEMTRI): $(FEMTRIOBJS)
	   @echo ''
	   @echo 'Creating executable bin/femTriangle.exe'
	   @$(FLINKER) $(LDFLAGS) -o $(EXEFEMTRI) $(FEMTRIOBJS) $(LIBFEM) $(LIBMBA) $(LIBS)

$(EXEFEMBC): $(FEMBCOBJS)
	   @echo ''
	   @echo 'Creating executable bin/femTriangle.exe'
	   @$(FLINKER) $(LDFLAGS) -o $(EXEFEMBC) $(FEMBCOBJS) $(LIBFEM) $(LIBMBA) $(LIBS)

$(EXEFEMTEM): $(FEMTEMOBJS)
	   @echo ''
	   @echo 'Creating executable bin/femTemplate.exe'
	   @$(FLINKER) $(LDFLAGS) -o $(EXEFEMTEM) $(FEMTEMOBJS) $(LIBFEM) $(LIBMBA) $(LIBS)


############################################################

