#############################################
#   D R O P S   local makefile              #
#############################################

# variables: 

DIR = poisson

EXEC = drops MGdrops testdrops dropsrep ipdrops 

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr dropsrep

drops: \
    ../geom/boundary.o ../geom/builder.o drops.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o poisson.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o
	$(CC) $(LFLAGS) -o $@ $^

dropsrep: \
    ../geom/boundary.o ../geom/builder.o dropsrep.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o poisson.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o
	$(CC) $(LFLAGS) -o $@ $^

MGdrops: \
    ../geom/boundary.o ../geom/builder.o MGdrops.o ../geom/multigrid.o \
    ../out/output.o ../num/unknowns.o ../geom/topo.o \
    ../num/MGsolver.o poisson.o ../misc/problem.o ../misc/utils.o
	$(CC) $(LFLAGS) -o $@ $^

testdrops: \
    ../geom/boundary.o ../geom/builder.o testdrops.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/MGsolver.o \
    poisson.o ../misc/problem.o ../misc/utils.o ../out/output.o \
    ../num/discretize.o
	$(CC) $(LFLAGS) -o $@ $^

ipdrops: \
    ../geom/boundary.o ../geom/builder.o ipdrops.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o instatpoisson.o \
    ../misc/problem.o ../misc/utils.o 
	$(CC) $(LFLAGS) -o $@ $^
  
matlab/ipdrops: ../geom/boundary.o ../geom/builder.o  \
    ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    instatpoisson.o ../misc/problem.o ../misc/utils.o    
	$(CC) $(CFLAGS) $(LFLAGS) matlab/ipdrops.cpp -o $@ $^


clean: matlab/clean

matlab/clean:
	cd matlab; make clean
  
# include common rules:
include ../Make.common
