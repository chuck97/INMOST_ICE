############################################################
ANIHOME = ..
include $(ANIHOME)/src/Rules.make

LIBS	= -L$(ANILIB) $(LIBVIEW) $(LIBAFT) $(LIBFEM) $(LIBMBA) \
	              $(LIBLMR) $(LIBLU) $(LIBLAPACK) $(LIBLAPACK_EXT) $(LIBBLAS) -lm 


############################################################
EXEMBA      = $(ANIBIN)/aniPY_mba.exe
EXESTOKES   = $(ANIBIN)/aniPY_stokes.exe

OBJMBA      = mainPYmba.o  
OBJSTOKES   = mainPYstokes.o forlibaft.o 

INCLUDE     = -I$(ANIFEM)

############################################################
all: 	    help

exe: 	    $(EXEMBA) $(EXESTOKES)

run-mba:
	    cd ../bin; ./aniPY_mba.exe

run-stokes:
	    cd ../bin; ./aniPY_stokes.exe

clean: 
	    @/bin/rm -f *.o  $(ANIBIN)/*exe  *_py.f  $(ANIBIN)/*ps

help:
	    @echo "make {exe|run-mba|run-stokes|clean|help}"
	    @echo " "
	    @echo "   exe        - compile the code"
	    @echo "   run-mba    - run code from package Ani2D-MBA"
	    @echo "   run-stokes - run code from tutotial Stokes"
	    @echo " "
	    @echo "   clean  - clean the package"
	    @echo "   help   - print help message"
	    @echo " "


############################################################
$(EXEMBA): $(OBJMBA) 
	@echo ''
	@echo 'Creating executable bin/aniPY_mba.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXEMBA) $(OBJMBA) $(LIBS)
	@echo ''

$(EXESTOKES): $(OBJSTOKES) 
	@echo ''
	@echo 'Creating executable bin/aniPY_stokes.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXESTOKES) $(OBJSTOKES) $(LIBS)
	@echo ''


############################################################
.f.o:
	@echo $(F77)' -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo $(CC)'  -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o 

