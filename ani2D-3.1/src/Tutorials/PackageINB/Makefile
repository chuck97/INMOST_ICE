############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make

LIBS	= -L$(ANILIB) $(LIBINB) $(LIBILU) $(LIBLAPACK) $(LIBBLAS) -lm


############################################################
EXE1     = $(ANIBIN)/aniINB_simple.exe 
EXE2     = $(ANIBIN)/aniINB_bratu.exe 


############################################################
EXEOBJ1  = main_simple.o 
EXEOBJ2  = main_bratu.o 


############################################################
all:	    help

exe:        $(EXE1) $(EXE2)  

clean:
	    @/bin/rm -f *.o *~ $(EXE1) $(EXE2) 

help:
	    @echo "make {exe|run-simple|run-bratu|clean|help}"
	    @echo " "
	    @echo "   exe   - compile the test code"
	    @echo " "
	    @echo "   run-simple - run a simple problem"
	    @echo "   run-bratu  - run the test code for Bratu problem"
	    @echo " "
	    @echo "   clean - clean the package"
	    @echo "   help  - print this message"
	    @echo " "

info:
	    @echo " "
	    @echo "Summary:"
	    @echo "Executables are located in bin/"
	    @echo " "

run-simple:
	    @cd $(ANIBIN) &&  ./aniINB_simple.exe

run-bratu:
	    @cd $(ANIBIN) &&  ./aniINB_bratu.exe



############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o 


############################################################
$(EXE1): $(EXEOBJ1)
	   @echo ''
	   @echo 'Creating executable bin/aniINBsimple.exe'
	   @$(FLINKER) $(LDFLAGS) -o $(EXE1) $(EXEOBJ1) $(LIBS)

$(EXE2): $(EXEOBJ2)
	   @echo ''
	   @echo 'Creating executable bin/aniINBbratu.exe'
	   @$(FLINKER) $(LDFLAGS) -o $(EXE2) $(EXEOBJ2) $(LIBS)


