#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = tests

EXEC = reftest ip1test ip2test mattest testfe

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr

reftest: \
    reftest.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../misc/utils.o ../geom/topo.o \
    ../out/output.o 
	$(CC) $(LFLAGS) -o $@ $^
        
ip1test: \
    interp1.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o 
	$(CC) $(LFLAGS) -o $@ $^

ip2test: \
    interp2.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o 
	$(CC) $(LFLAGS) -o $@ $^

mattest: \
    mattest.o ../misc/utils.o
	$(CC) $(LFLAGS) -o $@ $^

testfe: \
    testfe.o ../misc/utils.o ../geom/topo.o \
    ../geom/multigrid.o ../num/unknowns.o
	$(CC) $(LFLAGS) -o $@ $^



# include common rules:
include ../Make.common
