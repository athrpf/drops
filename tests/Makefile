#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = tests

EXEC =  ip1test ip2test mattest prolongationp2test reftest testfe tetrabuildertest 

DELETE = $(EXEC) *.out *.diff *.off *.mg *.dat

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)

TESTNAME = drivencavity-adaptiv-schur drivencavity-adaptiv-uzawa \
    navierstokes-stationaer-uniform \
    poisson-instationaer-uniform poisson-stationaer-adaptiv \
    poisson-stationaer-multigrid poisson-stationaer-neumann \
    poisson-stationaer-uniform stokes-instationaer-uniform \
    stokes-stationaer-uniform
    
# rules

default: ClrScr

check: $(TESTNAME:%=check_%) ip1check
	@echo "--> All testcases have been executed!"

check_%:
	cd ../$$(./testinfo "makedir" $*); make $$(../${DIR}/testinfo "maketarget" $*)
	./runtest $*

reftest: \
    ../tests/reftest.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../misc/utils.o ../geom/topo.o ../out/output.o
	$(CC) -o $@ $^ $(LFLAGS)

ip1test: \
    ../tests/interp1.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CC) -o $@ $^ $(LFLAGS)

ip1check: ip1test
	@ ./ip1test > /dev/null 2>&1 || echo "WARNING: P1-repair failed."

ip2test: \
    ../tests/interp2.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CC) -o $@ $^ $(LFLAGS)

prolongationp2test: \
    ../tests/prolongationp2test.o ../geom/multigrid.o ../geom/topo.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CC) -o $@ $^ $(LFLAGS)

tetrabuildertest: \
    ../tests/tetrabuildertest.o ../geom/multigrid.o ../geom/topo.o \
    ../num/unknowns.o ../geom/builder.o ../geom/boundary.o ../misc/utils.o
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
