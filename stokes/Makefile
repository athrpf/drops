#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = stokes

EXEC = sdrops sdropsP2

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr

sdrops: \
    ../geom/boundary.o ../geom/builder.o ../stokes/sdrops.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../num/interfacePatch.o \
    ../misc/utils.o ../num/discretize.o ../out/output.o
	$(CXX) -o $@ $^ $(LFLAGS)

sdropsP2: \
    ../stokes/stokesCoeff.o \
    ../geom/boundary.o ../geom/builder.o ../stokes/sdropsP2.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../num/interfacePatch.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../num/MGsolver.o ../misc/params.o \
    ../out/ensightOut.o ../stokes/integrTime.o \
    ../out/vtkOut.o ../geom/geomselect.o \
    ../misc/bndmap.o ../geom/bndVelFunctions.o 
	$(CXX) -o $@ $^ $(LFLAGS)

# include common rules:
include ../Make.common
