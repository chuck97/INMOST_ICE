ANIHOME = ./
include $(ANIHOME)/src/Rules.make


###########################################################
all: 	help


###########################################################
libs: 
	@echo "======================================="
	@echo "        Making libraries               "
	@cd $(ANIHOME)/src/blas;    make lib
	@cd $(ANIHOME)/src/lapack;  make lib
	@cd $(ANIHOME)/src/aniAFT;  make lib
	@cd $(ANIHOME)/src/aniFEM;  make lib
	@cd $(ANIHOME)/src/aniILU;  make lib
	@cd $(ANIHOME)/src/aniINB;  make lib
	@cd $(ANIHOME)/src/aniLU;   make lib
	@cd $(ANIHOME)/src/aniLMR;  make lib
	@cd $(ANIHOME)/src/aniMBA;  make lib
	@cd $(ANIHOME)/src/aniPRJ;  make lib
	@cd $(ANIHOME)/src/aniRCB;  make lib
	@cd $(ANIHOME)/src/aniC2F;  make lib
	@cd $(ANIHOME)/src/aniVIEW; make lib

packages: 
	@echo "======================================="
	@echo "        Making packages                "
	@echo "======================================="
	@cd $(ANIHOME)/src/Tutorials/PackageAFT;  make exe
	@cd $(ANIHOME)/src/Tutorials/PackageFEM;  make exe
	@cd $(ANIHOME)/src/Tutorials/PackageILU;  make exe
	@cd $(ANIHOME)/src/Tutorials/PackageINB;  make exe
	@cd $(ANIHOME)/src/Tutorials/PackageLU;   make exe
	@cd $(ANIHOME)/src/Tutorials/PackageLMR;  make exe
	@cd $(ANIHOME)/src/Tutorials/PackageMBA;  make exe
	@cd $(ANIHOME)/src/Tutorials/PackagePRJ;  make exe
	@cd $(ANIHOME)/src/Tutorials/PackageRCB;  make exe
	@cd $(ANIHOME)/src/Tutorials/PackageC2F;  make exe
	@cd $(ANIHOME)/src/Tutorials/MultiPackage/ConvectionDiffusion;  make exe
	@cd $(ANIHOME)/src/Tutorials/MultiPackage/HybridMixedFem;       make exe
	@cd $(ANIHOME)/src/Tutorials/MultiPackage/Interpolation;        make exe
	@cd $(ANIHOME)/src/Tutorials/MultiPackage/IsotropicBVP;         make exe
	@cd $(ANIHOME)/src/Tutorials/MultiPackage/NonisotropicBVP1;     make exe
	@cd $(ANIHOME)/src/Tutorials/MultiPackage/NonisotropicBVP2;     make exe
	@cd $(ANIHOME)/src/Tutorials/MultiPackage/PosterioriEstimates;  make exe
	@cd $(ANIHOME)/src/Tutorials/MultiPackage/Stokes;               make exe
	@cd $(ANIHOME)/src/Tutorials/MultiPackage/StokesNavier;         make exe
	@cd $(ANIHOME)/src/Tutorials/MultiPackage/UnsteadyConDif;       make exe
	@cd $(ANIHOME)/src/Tutorials/MultiPackage/DiscontinuousGalerkin;make exe
	@cd $(ANIHOME)/python;  make exe

test-aft: 
	@cd $(ANIHOME)/src/Tutorials/PackageAFT;  make exe run-frt gs

test-prj: 
	@cd $(ANIHOME)/src/Tutorials/PackagePRJ;  make exe run gs-grd gs-rsl

clean:
	@dirlist="src/ani* src/blas src/lapack python \
	src/Tutorials/Package* src/Tutorials/MultiPackage/[C,D,H,I,M,N,S,U,P]*"; \
	for dir in $$dirlist; do \
	   cd $$dir; make clean; cd "$(CURDIR)"; \
	done;

cleanlibs: 
	@cd $(ANIHOME)/lib;   /bin/rm -fr *.a

help:
	@echo " "
	@echo "make {libs|packages|clean|cleanlibs|help}" 
	@echo " "
	@echo "We recommend to use CMake to build the package Ani2D,"
	@echo "but will support regular Makefiles in versions 3.X"
	@echo " "
	@echo "   libs      - create all libraries"
	@echo "   packages  - create all executables"
	@echo " "
	@echo "   test-aft  - test library AFT (only 1 test)"
	@echo "   test-prj  - test library PRJ"
	@echo " "
	@echo "   clean     - remove temporary files except libraries"
	@echo "   cleanlibs - remove libraries"
	@echo " "
	@echo "   help      - print this information"
	@echo " "


