#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = stokes

EXEC = sdrops sdropsP2 FSdrops drcavad drivcav MGsdropsP2

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr

run: ClrScr drcavad
	drcavad 1e-14 1e-8 0 10 .01 .8 1 4

sdrops: \
    ../geom/boundary.o ../geom/builder.o ../stokes/sdrops.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o \
    ../misc/utils.o ../stokes/stokes.o ../num/discretize.o ../out/output.o
	$(CC) -o $@ $^ $(LFLAGS)

sdropsP2: \
    ../geom/boundary.o ../geom/builder.o ../stokes/sdropsP2.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o \
    ../misc/utils.o ../stokes/stokes.o ../out/output.o ../num/discretize.o ../num/MGsolver.o
	$(CC) -o $@ $^ $(LFLAGS)

MGsdropsP2: \
    ../geom/boundary.o ../geom/builder.o ../stokes/MGsdropsP2.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../num/MGsolver.o \
    ../misc/utils.o ../stokes/stokes.o ../out/output.o ../num/discretize.o
	$(CC) -o $@ $^ $(LFLAGS)

drivcav: \
    ../geom/boundary.o ../geom/builder.o ../stokes/drivcav.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o \
    ../misc/utils.o ../stokes/stokes.o ../out/output.o
	$(CC) -o $@ $^ $(LFLAGS)

drcavad: \
    ../geom/boundary.o ../geom/builder.o ../stokes/drcavad.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o \
    ../misc/utils.o ../stokes/stokes.o ../out/output.o ../num/discretize.o
	$(CC) -o $@ $^ $(LFLAGS)

FSdrops: \
    ../geom/boundary.o ../geom/builder.o ../stokes/FSdrops.o ../geom/multigrid.o \
    ../out/output.o ../num/unknowns.o ../geom/topo.o ../num/fe.o \
    ../misc/problem.o ../misc/utils.o ../stokes/instatstokes.o ../stokes/stokes.o
	$(CC) -o $@ $^ $(LFLAGS)


# include common rules:
include ../Make.common
