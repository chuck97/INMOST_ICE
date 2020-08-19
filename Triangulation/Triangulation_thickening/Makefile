CPPC = g++
CFLG = -c -Wall -Werror -std=c++1z
ANILIBPATH = ../../ani2D-3.1/lib
MAPXLIBPATH = ./mapx-0.10
ANILIBS = -laft2D-3.1 -lview2D-3.1 -lgfortran -lm 
MAPXLIBS = -lmapx

triangulation: netcdf_parsing main.o
	$(CPPC) -o triangulation.exe main.o -L $(ANILIBPATH) $(ANILIBS) -L $(MAPXLIBPATH) $(MAPXLIBS)
netcdf_parsing: 
	python parse_netcdf.py
main.o: main.cpp
	$(CPPC) $(CFLG) main.cpp
clean:
	rm -rf *.o triangulation.exe *.ps ./concentration_data/*.txt

