#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = navstokes

EXEC = nsdrops

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr

nsdrops: \
    ../geom/boundary.o ../geom/builder.o nsdrops.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o \
    ../misc/problem.o ../misc/utils.o ../stokes/stokes.o ../out/output.o \
    ../num/discretize.o 
	$(CC) $(LFLAGS) -o $@ $^


# include common rules:
include ../Make.common
