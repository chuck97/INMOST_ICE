############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make

LIBS	= -L$(ANILIB) $(LIBILU) $(LIBBLAS) -lm


############################################################
EXE0     = $(ANIBIN)/aniILU_bcg0.exe 
EXE1     = $(ANIBIN)/aniILU_gmr0.exe 
EXE2     = $(ANIBIN)/aniILU_bcg2.exe 


############################################################
EXEOBJ0  = main_bcg_ilu0.o 
EXEOBJ1  = main_gmres_ilu0.o 
EXEOBJ2  = main_bcg_ilu2.o 


############################################################
all:	    help

exe:        $(EXE0) $(EXE1) $(EXE2) 

clean:
	    @/bin/rm -f *.o *~ $(EXE0) $(EXE1) $(EXE2) $(ANIBIN)/CSRsystem

help:
	    @echo "make {exe|bcg-ilu0|gmres-ilu0|bcg-ilu2|clean|help}"
	    @echo " "
	    @echo "   exe  - compile the test code"
	    @echo " "
	    @echo "   bcg-ilu0   - run the test code for BCGstab-ILU0"
	    @echo "   gmres-ilu0 - run the test code for GMRES-ILU0"
	    @echo "   bcg-ilu2   - run the test code for BCGstab-ILU2"
	    @echo " "
	    @echo "   clean - clean the package"
	    @echo "   help  - print this message"
	    @echo " "

info:
	    @echo " "
	    @echo "Summary:"
	    @echo "Executables are located in bin/"
	    @echo " "

bcg-ilu0:
	    @cd $(ANIBIN) && ./aniILU_bcg0.exe

gmres-ilu0:
	    @cd $(ANIBIN) && ./aniILU_gmr0.exe

bcg-ilu2:
	    @cd $(ANIBIN) && ./aniILU_bcg2.exe


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o 


############################################################
$(EXE0): $(EXEOBJ0)
	   @echo ''
	   @echo 'Creating executable bin/aniILU_bcg0.exe'
	   @$(FLINKER) $(LDFLAGS) -o $(EXE0) $(EXEOBJ0) $(LIBS)

$(EXE1): $(EXEOBJ1)
	   @echo ''
	   @echo 'Creating executable bin/aniILU_gmr0.exe'
	   @$(FLINKER) $(LDFLAGS) -o $(EXE1) $(EXEOBJ1) $(LIBS)

$(EXE2): $(EXEOBJ2)
	   @echo ''
	   @echo 'Creating executable bin/aniILU_bcg2.exe'
	   @$(FLINKER) $(LDFLAGS) -o $(EXE2) $(EXEOBJ2) $(LIBS)


