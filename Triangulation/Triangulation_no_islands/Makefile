CPPC = g++
CFLG = -c -Wall -Werror -std=c++1z
LPATH = ../ani2D-3.1/lib
LIBS = -laft2D-3.1 -lview2D-3.1 -lgfortran -lm

triangulation: main.o
	$(CPPC) -o triangulation.exe main.o -L $(LPATH) $(LIBS)
main.o: main.cpp
	$(CPPC) $(CFLG) main.cpp
clean:
	rm -rf *.o triangulation.exe *.ps

