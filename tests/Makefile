#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = tests

EXEC = reftest iptest mattest testfe

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr

reftest: \
    reftest.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../geom/parallel.o ../num/unknowns.o ../misc/utils.o ../geom/topo.o \
    ../out/output.o 
	$(CC) $(LFLAGS) -o $@ $^
        
iptest: \
    interp1.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o \
    ../geom/parallel.o ../num/fe.o ../geom/boundary.o ../misc/utils.o 
	$(CC) $(LFLAGS) -o $@ $^

mattest: \
    mattest.o ../misc/utils.o
	$(CC) $(LFLAGS) -o $@ $^

testfe: \
    testfe.o ../misc/utils.o ../geom/topo.o \
    ../geom/parallel.o ../geom/multigrid.o ../num/unknowns.o
	$(CC) $(LFLAGS) -o $@ $^



# include common rules:
include ../Make.common
