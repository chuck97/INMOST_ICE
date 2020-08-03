############################################################
ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

############################################################
# YOU DO NOT NEED TO CHANGE THE LINES BELOW
############################################################
OBJLIBPRJ = metamesh.o rhs.o tri2tri.o


############################################################
all:	    help

lib:        $(LIBPRJ) info 

clean:
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.exe $(ANIBIN)/*.*ps

help:		
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib   - library libmba2D-*.a"
	    @echo "   clean - clean the package"
	    @echo "   help  - print this message"
	    @echo " "

info:	
	    @echo " " 
	    @echo "Library files are located in ../../lib"
	    @echo " "


############################################################
$(LIBPRJ): $(OBJLIBPRJ)
	   @echo ''
	   @rm -f $(LIBPRJ)
	   @ar -r $(LIBPRJ)  $(OBJLIBPRJ)


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o
 

