############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBS     = -L$(ANILIB) $(LIBVIEW) $(LIBPRJ) $(LIBFEM) $(LIBILU) \
	               $(LIBMBA) $(LIBLAPACK) $(LIBBLAS) -lm 

############################################################
EXEPRJ	= $(ANIBIN)/aniPRJ.exe
OBJPRJ	= main.o  


############################################################
all:	    help

exe:        $(EXEPRJ) info

run:        $(EXECRB)
	    @cd $(ANIBIN) && ./aniPRJ.exe

clean:
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.exe $(ANIBIN)/*.ps $(ANIBIN)/save.ani

help:		
	    @echo "make {exe|run|gs-grd|gs-rsl|clean|help}"
	    @echo " "
	    @echo "   exe - create four executables in bin/"
	    @echo "   run - run the remap test"
	    @echo " "
	    @echo "   gs-grd - draw two meshes"
	    @echo "   gs-rsl - draw meshes intersection and remapped solution"
	    @echo " "
	    @echo "   clean  - clean auto-built files"
	    @echo "   help   - print this help message"
	    @echo " "

info:	
	    @echo " " 
	    @echo "Executables are located in ../../../bin/"
	    @echo " "

gs-grd:	    
	    @cd $(ANIBIN); $(PS) mesh_first.ps; $(PS) mesh_second.ps

gs-rsl:	    
	    @cd $(ANIBIN); $(PS) mesh_final.ps; $(PS) isolines.ps


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o


$(EXEPRJ): $(OBJPRJ)
	@echo ''
	@echo 'Creating executable $(EXEPRJ)'
	@$(FLINKER) $(LDFLAGS) -o $(EXEPRJ) $(OBJPRJ) $(LIBS)
	@echo ''


