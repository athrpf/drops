#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DROPS_ROOT = ..

DIR = geom

EXEC =
DELETE = $(EXEC) maketopo_helper.cpp

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules:

default: ClrScr

maketopo: maketopo.cpp topo.h topo.cpp.in ../misc/utils.o
	cp -f topo.cpp.in maketopo_helper.cpp
	$(CXX) $(CXXFLAGS) -c maketopo_helper.cpp
	$(CXX) $(CXXFLAGS) -o maketopo maketopo.cpp maketopo_helper.o ../misc/utils.o $(LFLAGS)
	rm -f maketopo_helper.cpp maketopo_helper.o

topo.cpp: maketopo topo.cpp.in
	./maketopo

topoclean:
	rm -f topo.cpp

# include common rules:
include ../Make.common
