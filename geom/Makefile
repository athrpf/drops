#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = geom

EXEC = maketopo

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(filter-out topo.o,$(CPP:%.cpp=%.o))


# rules:

default: ClrScr

topo.o: topo.cpp topo.h ../misc/utils.h
	$(CC) $(CFLAGS) -c topo.cpp

maketopo: maketopo.cpp topo.h topo.cpp.in ../misc/utils.o
	cp -f topo.cpp.in maketopo_helper.cpp
	$(CC) $(CFLAGS) -c maketopo_helper.cpp
	$(CC) $(CFLAGS) $(LFLAGS) -o maketopo maketopo.cpp maketopo_helper.o ../misc/utils.o
	rm -f maketopo_helper.cpp maketopo_helper.o

topo.cpp: maketopo topo.cpp.in
	./maketopo


# include common rules:
include ../Make.common
