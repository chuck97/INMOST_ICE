############################################################
# SET UP THE MODEL (crv_model) 
############################################################
OBJUSER = crv_model.o  # only if the domain is given analytically


############################################################
# YOU DO NOT NEED TO CHANGE THE LINES BELOW
############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBS	= -L$(ANILIB) $(LIBAFT) $(LIBVIEW) -lm  


############################################################
EXEFRONT1    = $(ANIBIN)/aniAFT_front1.exe
EXEFRONT2    = $(ANIBIN)/aniAFT_front2.exe
EXEBNDSQUARE = $(ANIBIN)/aniAFT_boundary_square.exe
EXEBNDWING   = $(ANIBIN)/aniAFT_boundary_wing.exe

OBJFRONT1    = main_front1.o
OBJFRONT2    = main_front2.o
OBJBNDSQUARE = main_boundary_square.o 
OBJBNDWING   = main_boundary_wing.o $(OBJUSER)


############################################################
all: 	    help

exe: 	    $(EXEFRONT1) $(EXEFRONT2) $(EXEBNDSQUARE) $(EXEBNDWING) info

run-frt:
	    cd $(ANIBIN) && ./aniAFT_front1.exe

run-fr2:
	    cd $(ANIBIN) && ./aniAFT_front2.exe

run-bnd:
	    cd $(ANIBIN) && ./aniAFT_boundary_square.exe

run-ana:
	    cd $(ANIBIN) && ./aniAFT_boundary_wing.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps \
            $(ANIBIN)/*.exe $(ANIBIN)/gmon.out 

help:
	    @echo "make {exe|run-frt|run-fr2|run-bnd|run-ana|clean|help|gs|gsgs}"
	    @echo " "
	    @echo "   exe - compile and link the code"
	    @echo " "
	    @echo "   run-frt - run the code with a given front (first way)"
	    @echo "   run-fr2 - run the code with a given front (second way)"
	    @echo "   run-bnd - run the code with a given polygonal boundary"
	    @echo "   run-ana - run the code with a given analytical boundary"
	    @echo " "
	    @echo "   gs    - draw the mesh (bin/mesh_final.ps)"
	    @echo " "
	    @echo "   clean - clean the package"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo " "
	    @echo "Executable files are located in   bin"
	    @echo "Postscript figures are located in bin"
	    @echo " "

gs:
	    @cd $(ANIBIN) && $(PS) mesh_final.ps

gsgs:
	    @cd $(ANIBIN) && $(PS) mesh_initial.ps
	    @cd $(ANIBIN) && $(PS) mesh_final.ps



############################################################
$(EXEFRONT1): $(OBJFRONT1) 
	@echo ''
	@echo 'Creating executable $(EXEFRONT1)' 
	@$(FLINKER) $(LDFLAGS) -o $(EXEFRONT1) $(OBJFRONT1) $(LIBS)
	@echo ''

$(EXEFRONT2): $(OBJFRONT2) 
	@echo ''
	@echo 'Creating executable $(EXEFRONT2)' 
	@$(FLINKER) $(LDFLAGS) -o $(EXEFRONT2) $(OBJFRONT2) $(LIBS)
	@echo ''

$(EXEBNDSQUARE): $(OBJBNDSQUARE) 
	@echo ''
	@echo 'Creating executable $(EXEBDSQUARE)' 
	@$(FLINKER) $(LDFLAGS) -o $(EXEBNDSQUARE) $(OBJBNDSQUARE) $(LIBS)
	@echo ''

$(EXEBNDWING): $(OBJBNDWING) 
	@echo ''
	@echo 'Creating executable $(EXEBDWING)' 
	@$(FLINKER) $(LDFLAGS) -o $(EXEBNDWING) $(OBJBNDWING) $(LIBS)
	@echo ''


############################################################
main_front1.o   : main_front1.f 
main_front2.o   : main_front2.f 
main_boundary.o : main_boundary.f  
crv_model.o     : crv_model.c


############################################################
.f.o:
	@echo $(F77)' -c ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo $(CC)'  -c ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

