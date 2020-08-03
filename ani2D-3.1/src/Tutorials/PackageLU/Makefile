############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBS	= -L$(ANILIB) $(LIBLAPACK) $(LIBBLAS) -lm


############################################################
EXE4    = $(ANIBIN)/aniLU.exe

OBJ4    = main.o aux_umf.o


############################################################
all: 	    help

exe: 	    $(EXE4) info

run:
	    @cd $(ANIBIN) && ./aniLU.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps \
            $(ANIBIN)/*.exe $(ANIBIN)/CSRsystem

help:
	    @echo "make {exe|run|clean|help|}"
	    @echo " "
	    @echo "   exe - create executables"
	    @echo " "
	    @echo "   run - solve the finite element system with liblu.a"
	    @echo " "
	    @echo "   clean - clean this directory"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo " "
	    @echo "Executable files are located in   bin"
	    @echo "Postscript figures are located in bin"
	    @echo " "


############################################################
$(EXE4): $(OBJ4) 
	@echo ''
	@echo 'Creating executable bin/aniLU.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE4) $(OBJ4) $(LIBLU) $(LIBS)
	@echo ''


############################################################
mainSolFemSys.o: mainSolFemSys.f 
aux_umf.o:       aux_umf.f


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

