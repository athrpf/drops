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
    ../tests/reftest.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../misc/utils.o ../geom/topo.o ../out/output.o
	$(CC) -o $@ $^ $(LFLAGS)

ip1test: \
    ../tests/interp1.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CC) -o $@ $^ $(LFLAGS)

ip2test: \
    ../tests/interp2.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CC) -o $@ $^ $(LFLAGS)

mattest: \
    ../tests/mattest.o ../misc/utils.o
	$(CC) -o $@ $^ $(LFLAGS)

testfe: \
    ../tests/testfe.o ../misc/utils.o ../geom/topo.o \
    ../geom/multigrid.o ../num/unknowns.o
	$(CC) -o $@ $^ $(LFLAGS)


# include common rules:
include ../Make.common
