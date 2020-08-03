############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBS	= -L$(ANILIB) $(LIBVIEW) $(LIBRCB) -lm 
	                

############################################################
EXERCB  = $(ANIBIN)/aniRCB.exe 


############################################################
RCBOBJS = main.o 


############################################################
all:	    help

exe:        $(EXERCB) info

clean:
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.exe $(ANIBIN)/*.ps 

help:
	    @echo "make {exe|run|gs-ini|gs-ref|gs-crs|clean|help}"
	    @echo " "
	    @echo "   exe - compile the code"
	    @echo " "
	    @echo "   run - run the code"
	    @echo " "
	    @echo "   gs-ini - draw initial mesh (bin/mesh_initial.ps)"
	    @echo "   gs-ref - draw refined mesh (bin/mesh_final.ps)"
	    @echo "   gs-crs - draw coarsened mesh (bin/mesh_coarse.ps)"
	    @echo " "
	    @echo "   clean - clean the package"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo "Summary:"
	    @echo "Executables are located in bin/"
	    @echo " "

run:
	    @cd $(ANIBIN) && ./aniRCB.exe

gs-ini:
	    @cd $(ANIBIN) && $(PS) mesh_initial.ps

gs-ref:
	    @cd $(ANIBIN) && $(PS) mesh_final.ps

gs-crs:
	    @cd $(ANIBIN) && $(PS) mesh_coarse.ps


############################################################
.f.o:
	   @echo '$(F77) -c $(FFLAGS) ' $*.f
	   @$(F77) $(FFLAGS) -c $*.f -o $*.o 


############################################################
$(EXERCB): $(RCBOBJS) 
	   @echo ''
	   @echo 'Creating executable bin/rcbBisect.exe'
	   @$(FLINKER) $(LDFLAGS) -o $(EXERCB) $(RCBOBJS) $(LIBS)


