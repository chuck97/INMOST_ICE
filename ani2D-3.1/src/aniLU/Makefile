############################################################
ANIHOME = ../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBLU   = $(ANILIB)/liblu-5.1.a
LIBS    = -L$(ANILIB) -llu  -lblas-3.2

AMD     = $(ANIHOME)/src/aniLU/AMD
UMFPACK = $(ANIHOME)/src/aniLU/UMFPACK
UFCONFIG= $(ANIHOME)/src/aniLU/UFconfig

############################################################
OPT = -O

CFLAGS = -c  $(OPT)  
FFLAGS = -c  $(OPT)  

LD = $(F77)
LFLAGS = $(OPT) 

############################################################
all:    help

help:
	@echo "make {lib|clean|help}"
	@echo " "
	@echo " Library option:"
	@echo "    lib   - make library liblu.a"
	@echo "    clean - clean this directory and packages AMD, UMFPACK and bin"
	@echo "    help  - print this message"
	@echo " "

clean :
	@rm -rf *.o $(EXEC) $(AMD)/Lib/*.o $(AMD)/Lib/*.a \
	       $(UMFPACK)/Lib/*.o  $(UMFPACK)/Lib/*.a \
	       UMFPACK/Demo/umf4_f77wrapper.o

lib:
	@cd $(AMD)/Lib; make
	@cd $(UMFPACK)/Lib; make
	@$(CC) $(CFLAGS) $(UMFPACK)/Demo/umf4_f77wrapper.c -o $(UMFPACK)/Demo/umf4_f77wrapper.o \
	      -I$(UMFPACK)/Include -I$(AMD)/Include -I$(UFCONFIG) 
	@rm -f $(LIBLU)
	@ar -r $(LIBLU) $(UMFPACK)/Demo/umf4_f77wrapper.o
	@mergelib $(LIBLU) $(AMD)/Lib/libamd.a 
	@mergelib $(LIBLU) $(UMFPACK)/Lib/libumfpack.a 


############################################################
.c.o :
	$(CC) $(CFLAGS) $*.c -o $*.o
.f.o :
	$(F77) $(FFLAGS) $*.f -o $*.o

