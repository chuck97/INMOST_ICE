############################################################
ANIHOME = ../..
include $(ANIHOME)/src/Rules.make


############################################################
OBJLIBAFT = memory2.o aft2d.o refine2.o \
            region2.o struct2.o tree2.o \
	    tria2.o user2.o


############################################################
all:	    help

lib:	    $(LIBAFT) info

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps \
            $(ANIBIN)/*.exe $(ANIBIN)/gmon.out 

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib - make the library libaft2D.a"
	    @echo " "
	    @echo "   clean - clean this directory"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo " "
	    @echo "Library files are located in      lib"
	    @echo " "


############################################################
$(LIBAFT): $(OBJLIBAFT) 
	@echo ''
	@rm -f $(LIBAFT)
	@ar -r $(LIBAFT) $(OBJLIBAFT)


############################################################
.f.o:
	@echo $(F77)' -c ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo $(CC)'  -c  $(CFLAGS)' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

