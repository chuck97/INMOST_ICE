############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBS     = -L$(ANILIB) $(LIBLMR) $(LIBMBA) \
	               $(LIBLAPACK) $(LIBBLAS) -lm


############################################################
EXELMR1	= $(ANIBIN)/aniLMR_nodal.exe
EXELMR2	= $(ANIBIN)/aniLMR_edgeest.exe
EXELMR3	= $(ANIBIN)/aniLMR_gradient.exe

OBJLMR1	= mainNodal2Metric.o 
OBJLMR2	= mainEst2Metric.o esterr.o
OBJLMR3 = mainFunc2GradMetric.o


############################################################
all:	    help

exe:        $(EXELMR1) $(EXELMR2) $(EXELMR3) info

run-nod:    $(EXELMR1)
	    @cd $(ANIBIN) && ./aniLMR_nodal.exe

run-est:    $(EXELMR2)
	    @cd $(ANIBIN) && ./aniLMR_edgeest.exe

run-grd:    $(EXELMR3)
	    @cd $(ANIBIN) && ./aniLMR_gradient.exe

clean:
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.exe $(ANIBIN)/metric

help:		
	    @echo "make {exe|run-nod|run-est|run-grd|clean|help}"
	    @echo " "
	    @echo "   exe - create three executables"
	    @echo "   run-nod - generate metric from discrete function"
	    @echo "   run-est - generate metric from error estimates"
	    @echo "   run-grd - generate gradient metric from discrete function"
	    @echo " "
	    @echo "   clean - clean the package"
	    @echo "   help  - print help message"
	    @echo " "

info:	
	    @echo " " 
	    @echo "Summary:"
	    @echo "Executables are located in ../../../bin/"
	    @echo " "


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o


$(EXELMR1): $(OBJLMR1)
	@echo ''
	@echo 'Creating executable $(EXELMR1)'
	$(FLINKER) $(LDFLAGS) -o $(EXELMR1) $(OBJLMR1) $(LIBS)
	@echo ''

$(EXELMR2): $(OBJLMR2)
	@echo ''
	@echo 'Creating executable $(EXELMR2)'
	@$(FLINKER) $(LDFLAGS) -o $(EXELMR2) $(OBJLMR2) $(LIBS)
	@echo ''

$(EXELMR3): $(OBJLMR3)
	@echo ''
	@echo 'Creating executable $(EXELMR3)'
	@$(FLINKER) $(LDFLAGS) -o $(EXELMR3) $(OBJLMR3) $(LIBS)
	@echo ''

############################################################


