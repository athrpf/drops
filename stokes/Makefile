#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = stokes

EXEC = sdrops sdropsP2 FSdrops drcavad drivcav

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr

run: ClrScr drcavad
	drcavad 1e-14 1e-8 0 10 .01 .8 1 4

sdrops: \
    ../geom/boundary.o ../geom/builder.o sdrops.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o \
    ../misc/problem.o ../misc/utils.o ../stokes/stokes.o ../num/discretize.o \
    ../out/output.o 
	$(CC) $(LFLAGS) -o $@ $^

sdropsP2: \
    ../geom/boundary.o ../geom/builder.o sdropsP2.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o \
    ../misc/problem.o ../misc/utils.o stokes.o ../out/output.o \
    ../num/discretize.o
	$(CC) $(LFLAGS) -o $@ $^

drivcav: \
    ../geom/boundary.o ../geom/builder.o drivcav.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o \
    ../misc/problem.o ../misc/utils.o stokes.o ../out/output.o 
	$(CC) $(LFLAGS) -o $@ $^

drcavad: \
    ../geom/boundary.o ../geom/builder.o drcavad.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o \
    ../misc/problem.o ../misc/utils.o stokes.o ../out/output.o \
    ../num/discretize.o
	$(CC) $(LFLAGS) -o $@ $^

FSdrops: \
    ../geom/boundary.o ../geom/builder.o FSdrops.o ../geom/multigrid.o \
    ../out/output.o ../num/unknowns.o ../geom/topo.o \
    ../num/fe.o ../misc/problem.o ../misc/utils.o instatstokes.o stokes.o
	$(CC) $(LFLAGS) -o $@ $^


# include common rules:
include ../Make.common
