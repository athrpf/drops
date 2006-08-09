#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = navstokes

EXEC = nsdrops nsdrops_begehung insdrops insadrops

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr

nsdrops: \
    ../geom/boundary.o ../geom/builder.o ../navstokes/nsdrops.o \
    ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/fe.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

nsdrops_begehung: \
    ../geom/boundary.o ../geom/builder.o ../navstokes/nsdrops_begehung.o \
    ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/fe.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

insdrops: \
    ../geom/boundary.o ../geom/builder.o ../navstokes/insdrops.o \
    ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/fe.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

insadrops: \
    ../geom/boundary.o ../geom/builder.o ../navstokes/insadrops.o \
    ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/fe.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)


# include common rules:
include ../Make.common
