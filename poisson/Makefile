#############################################
#   D R O P S   local makefile              #
#############################################

# variables: 

DIR = poisson

EXEC = drops MGdrops testdrops

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr dropsrep

drops: \
    ../geom/boundary.o ../geom/builder.o drops.o ../geom/multigrid.o \
    ../geom/parallel.o ../num/unknowns.o ../geom/topo.o poisson.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o
	$(CC) $(LFLAGS) -o $@ $^

MGdrops: \
    ../geom/boundary.o ../geom/builder.o MGdrops.o ../geom/multigrid.o \
    ../geom/parallel.o ../out/output.o ../num/unknowns.o ../geom/topo.o \
    ../num/MGsolver.o poisson.o ../misc/problem.o ../misc/utils.o
	$(CC) $(LFLAGS) -o $@ $^

testdrops: \
    ../geom/boundary.o ../geom/builder.o testdrops.o ../geom/multigrid.o \
    ../geom/parallel.o ../num/unknowns.o ../geom/topo.o ../num/MGsolver.o \
    poisson.o ../misc/problem.o ../misc/utils.o ../out/output.o \
    ../num/discretize.o
	$(CC) $(LFLAGS) -o $@ $^


# include common rules:
include ../Make.common
