############################################################
ANIHOME = ../..
include $(ANIHOME)/src/Rules.make


############################################################
OBJCORE     = ani2D.o


############################################################
all: 	    help

lib:	    $(LIBC2F) info

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps \
            $(ANIBIN)/*.exe $(ANIBIN)/gmon.out 

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib - make the library libc2f2D.a"
	    @echo " "
	    @echo "   clean - clean this directory"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo " "
	    @echo "Library files are located in      lib"
	    @echo " "


############################################################
$(LIBC2F): $(OBJCORE) 
	@echo ''
	@rm -f $(LIBC2F)
	@ar -r $(LIBC2F) $(OBJCORE)


############################################################
ani2D.o         : ani2D.h


############################################################
.f.o:
	@echo $(F77)' -c ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo $(CC)'  -c ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

