############################################################
ANIHOME = ../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBLAPACK   = $(ANILIB)/liblapack-3.2.a 

OBJLAPACK   = dsyev.o  ilaenv.o dlansy.o \
              dlascl.o dsytrd.o dsterf.o dorgtr.o \
              dsteqr.o ieeeck.o dlassq.o dlatrd.o \
              dsytd2.o dlanst.o dlae2.o  dlapy2.o \
              dlasrt.o dorgql.o dorgqr.o dlaset.o \
              dlaev2.o dlasr.o  dlartg.o dlarfg.o \
              dorg2l.o dlarft.o dorg2r.o dlarfb.o \
              dlarf.o  dsysv.o  dsytrf.o dsytrs.o \
	      dsytf2.o dlasyf.o iparmq.o disnan.o \
	      dgesv.o  dgetrf.o dgetrs.o iladlr.o \
	      dgetf2.o dlaswp.o iladlr.o iladlc.o \
	      dpotrf.o dpotf2.o dposv.o  dpotrs.o \
	      dlaisnan.o


############################################################
all: 	    help

lib:	    $(LIBLAPACK) info

clean: 
	    @/bin/rm -f *.o $(ANILAPACK)/*.o *~  

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib - library liblapack-3.2.a"
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
$(LIBLAPACK): $(OBJLAPACK)
	@rm -f $(ANILIB)/liblapack-3.2.a $(ANILIB)/liblapack_ext-3.2.a
	@$(F77) -c extension/dlamch.f -o dlamch.o $(INCLUDE)
	@ar -r $(ANILIB)/liblapack-3.2.a  $(OBJLAPACK) 
	@ar -r $(ANILIB)/liblapack_ext-3.2.a  dlamch.o


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

