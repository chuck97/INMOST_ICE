############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBS = -L$(ANILIB) $(LIBC2F) $(LIBVIEW) $(LIBMBA) \
                   $(LIBLAPACK) $(LIBBLAS) -lm  $(LIBSYS)


############################################################
EXEANA = $(ANIBIN)/aniC2F_analytic.exe
EXENOD = $(ANIBIN)/aniC2F_nodal.exe
EXEAUD = $(ANIBIN)/aniC2F_audit.exe

OBJANA = main_analytic.o
OBJNOD = main_nodal.o
OBJAUD = main_audit.o


############################################################
all: 	    help

exe: 	    $(EXEANA) $(EXENOD) $(EXEAUD) info

run-ana:
	    cd $(ANIBIN) &&  ./aniC2F_analytic.exe

run-nod:
	    cd $(ANIBIN) &&  ./aniC2F_nodal.exe

run-aud:
	    cd $(ANIBIN) &&  ./aniC2F_audit.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps \
            $(ANIBIN)/*.exe $(ANIBIN)/gmon.out $(ANIDAT)/save.ani

help:
	    @echo "make {exe|run-ana|run-nod|clean|help}"
	    @echo " "
	    @echo "   exe - compile the code"
	    @echo " "
	    @echo "   run-ana - run code with an analytic metric"	
	    @echo "   run-nod - run code with a discrete metric"
	    @echo "   run-aud - run mesh audit"
	    @echo " "
	    @echo "   gs - show PostScript files with initial and final meshes"
	    @echo " "
	    @echo "   clean - clean this directory"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo " "
	    @echo "Executable files are located in   bin"
	    @echo "Postscript figures are located in bin"
	    @echo " "

gs:
	    @cd $(ANIBIN); $(PS) initial_mesh.ps; $(PS) final_mesh.ps


############################################################
$(EXEANA): $(OBJANA) 
	@echo ''
	@echo 'Creating executable $(EXEANA)' 
	@$(CLINKER) $(LDFLAGS) -o $(EXEANA) $(OBJANA) $(LIBS)
	@echo ''

$(EXENOD): $(OBJNOD) 
	@echo ''
	@echo 'Creating executable $(EXENOD)' 
	@$(CLINKER) $(LDFLAGS) -o $(EXENOD) $(OBJNOD) $(LIBS)
	@echo ''

$(EXEAUD): $(OBJAUD) 
	@echo ''
	@echo 'Creating executable $(EXEAUD)' 
	@$(CLINKER) $(LDFLAGS) -o $(EXEAUD) $(OBJAUD) $(LIBS)
	@echo ''


############################################################
.f.o:
	@echo $(F77)' -c ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo $(CC)'  -c ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o $(INCLUDE)

