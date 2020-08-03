############################################################
# YOU DO NOT NEED TO CHANGE THE LINES BELOW
############################################################
ANIHOME = ../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBBLAS  =  $(ANILIB)/libblas-3.2.a 

OBJBLAS  =  daxpy.o dcopy.o ddot.o   dnrm2.o  \
            dgemm.o dgemv.o dger.o   dtrsv.o  \
            dtrsm.o dscal.o dsyr.o            \
            dswap.o dsymv.o dsyr2.o  dsyr2k.o dsyrk.o \
            drot.o drotg.o             \
            dtrmm.o dtrmv.o idamax.o lsame.o  \
	    xerbla.o


############################################################
all: 	    help

lib:	    $(LIBBLAS) info

clean: 
	    @/bin/rm -f *.o $(ANIBLAS)/*.o *~  

help:
	    @echo "make {lib|clean|help}"
	    @echo "   home directory: " $(ANILIB)
	    @echo " "
	    @echo "   lib   - create library libblas-3.2.a"
	    @echo " "
	    @echo "   clean - clean the package"
	    @echo "   help  - print this message"
	    @echo " "

info:
	    @echo " "
	    @echo "Executable files are located in   bin"
	    @echo "Library files are located in      lib"
	    @echo "Postscript figures are located in bin"
	    @echo " "


############################################################
$(LIBBLAS): $(OBJBLAS)
	@echo ''
	@rm -f $(ANILIB)/libblas-3.2.a
	@ar -r $(ANILIB)/libblas-3.2.a  $(OBJBLAS)
	@echo 'If this library has been installed set LIBS = -lblas -lm'


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

