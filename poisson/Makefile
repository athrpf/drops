#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = poisson

EXEC = drops MGdrops testdrops dropsrep ipdrops ipdropsrep drops_neumann \
       MGdropsP2 ihcpfilm

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr dropsrep

drops: \
    ../geom/boundary.o ../geom/builder.o ../poisson/drops.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../poisson/poisson.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o
	$(CC) -o $@ $^ $(LFLAGS)

dropsrep: \
    ../geom/boundary.o ../geom/builder.o ../poisson/dropsrep.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../poisson/poisson.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o
	$(CC) -o $@ $^ $(LFLAGS)

drops_neumann: \
    ../geom/boundary.o ../geom/builder.o ../poisson/drops_neumann.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../poisson/poisson.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o
	$(CC) -o $@ $^ $(LFLAGS)

MGdrops: \
    ../geom/boundary.o ../geom/builder.o ../poisson/MGdrops.o ../geom/multigrid.o \
    ../out/output.o ../num/unknowns.o ../geom/topo.o \
    ../num/MGsolver.o ../poisson/poisson.o ../misc/problem.o ../misc/utils.o \
    ../num/fe.o
	$(CC) -o $@ $^ $(LFLAGS)

MGdropsP2: \
    ../geom/boundary.o ../geom/builder.o ../poisson/MGdropsP2.o ../geom/multigrid.o \
    ../out/output.o ../num/unknowns.o ../geom/topo.o \
    ../num/MGsolver.o ../poisson/poisson.o ../misc/problem.o ../misc/utils.o \
    ../num/fe.o
	$(CC) -o $@ $^ $(LFLAGS)

testdrops: \
    ../geom/boundary.o ../geom/builder.o ../poisson/testdrops.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/MGsolver.o ../poisson/poisson.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o \
    ../num/fe.o
	$(CC) -o $@ $^ $(LFLAGS)
  
ipdrops: \
    ../geom/boundary.o ../geom/builder.o ../poisson/ipdrops.o \
    ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../num/MGsolver.o ../poisson/instatpoisson.o ../misc/problem.o ../misc/utils.o \
    ../num/fe.o
	$(CC) -o $@ $^ $(LFLAGS)

ipdropsrep: \
    ../geom/boundary.o ../geom/builder.o ../poisson/ipdropsrep.o \
    ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/instatpoisson.o ../misc/problem.o ../misc/utils.o
	$(CC) -o $@ $^ $(LFLAGS)

ipfilm: \
    ../geom/boundary.o ../geom/builder.o ../poisson/ipfilm.o \
    ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/instatpoisson.o ../misc/problem.o ../misc/utils.o
	$(CC) -o $@ $^ $(LFLAGS)

matlab/ipdrops: ../geom/boundary.o ../geom/builder.o \
    ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/instatpoisson.o ../misc/problem.o ../out/output.o ../misc/utils.o
	$(CC) $(CFLAGS) matlab/ipdrops.cpp -o $@ $^ $(LFLAGS)

matlab/ipfilm: ../geom/boundary.o ../geom/builder.o \
    ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/instatpoisson.o ../misc/problem.o ../out/output.o ../misc/utils.o
	$(CC) $(CFLAGS) matlab/ipfilm.cpp -o $@ $^ $(LFLAGS)


clean: matlab/clean

matlab/clean:
	cd matlab && $(MAKE) clean

# include common rules:
include ../Make.common
