############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBS     = -L$(ANILIB) $(LIBVIEW) $(LIBMBA) $(LIBLMR) \
	               $(LIBLAPACK) $(LIBBLAS) -lm 

############################################################
EXEANA	= $(ANIBIN)/aniMBA_analytic.exe
EXENOD	= $(ANIBIN)/aniMBA_nodal.exe
EXEFIX	= $(ANIBIN)/aniMBA_fixshape.exe
EXETNG	= $(ANIBIN)/aniMBA_untangle.exe
EXETRI	= $(ANIBIN)/aniMBA_triangle.exe


############################################################
OBJANA	= main_analytic.o 
OBJNOD	= main_nodal.o 
OBJFIX	= main_fixshape.o 
OBJTNG	= main_tangled.o 
OBJTRI  = main_triangle.o


############################################################
all:	    help

exe:        $(EXEANA) $(EXENOD) $(EXEFIX) $(EXETNG) $(EXETRI) info

run-ana:    $(EXEANA)
	    @cd $(ANIBIN) && ./aniMBA_analytic.exe

run-nod:    $(EXENOD)
	    @cd $(ANIBIN) && ./aniMBA_nodal.exe

run-fix:    $(EXEFIX)
	    @cd $(ANIBIN) && ./aniMBA_fixshape.exe

run-tng:    $(EXETNG)
	    @cd $(ANIBIN) && ./aniMBA_untangle.exe

run-tri:    $(EXETRI)
	    @cd $(ANIBIN) && ./aniMBA_triangle.exe

clean:
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.exe $(ANIBIN)/*.ps $(ANIBIN)/save.ani

help:		
	    @echo "make {exe|run-ana|run-nod|run-fix|run-tng|gs-ini|gs-fin|clean|help}"
	    @echo " "
	    @echo "   exe - create four executables in bin/"
	    @echo " "
	    @echo "   run-ana - run the code with an analytical metric"
	    @echo "   run-nod - run the code with a user-provided discrete nodal metric"
	    @echo "   run-fix - improve shape quality of a given mesh"
	    @echo "   run-tng - run the code that untangles mesh"
	    @echo "   run-tri - run the code that reads mesh in the Triangle format"
	    @echo " "
	    @echo "   gs-ini - draw the initial mesh"
	    @echo "   gs-fin - draw the final mesh"
	    @echo " "
	    @echo "   clean  - clean auto-built files"
	    @echo "   help   - print this help message"
	    @echo " "
	    @echo "2-lines: Package aniMBA generates unstructured triangular meshes that"
	    @echo "are quasi-uniform in a user defined metric, analytical or discrete."
	    @echo " "

info:	
	    @echo " " 
	    @echo "Executables are located in ../../../bin/"
	    @echo " "

gs-ini:	    
	    @cd $(ANIBIN); $(PS) mesh_initial.ps

gs-fin:	    
	    @cd $(ANIBIN); $(PS) mesh_final.ps


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o


$(EXEANA): $(OBJANA)
	@echo ''
	@echo 'Creating executable $(EXEANA)'
	@$(FLINKER) $(LDFLAGS) -o $(EXEANA) $(OBJANA) $(LIBS)
	@echo ''

$(EXENOD): $(OBJNOD)
	@echo ''
	@echo 'Creating executable $(EXENOD)'
	@$(FLINKER) $(LDFLAGS) -o $(EXENOD) $(OBJNOD) $(LIBS)
	@echo ''

$(EXEFIX): $(OBJFIX)
	@echo ''
	@echo 'Creating executable $(EXEFIX)'
	@$(FLINKER) $(LDFLAGS) -o $(EXEFIX) $(OBJFIX) $(LIBS)
	@echo ''

$(EXETNG): $(OBJTNG)
	@echo ''
	@echo 'Creating executable $(EXETNG)'
	@$(FLINKER) $(LDFLAGS) -o $(EXETNG) $(OBJTNG) $(LIBS)
	@echo ''

$(EXETRI): $(OBJTRI)
	@echo ''
	@echo 'Creating executable $(EXETRI)'
	@$(FLINKER) $(LDFLAGS) -o $(EXETRI) $(OBJTRI) $(LIBS)
	@echo ''

