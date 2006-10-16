#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = tests

EXEC =  ip1test ip2test mattest prolongationp2test reftest testfe \
        tetrabuildertest sbuffer minres meshreader restrictp2 vectest \
        p2local quadbase globallist triang quadCut bicgstab gcr blockmat \
        mass quad5 downwind

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
	cd ../$$(./testinfo "makedir" $*) && $(MAKE) $$(../${DIR}/testinfo "maketarget" $*)
	./runtest $*

reftest: \
    ../tests/reftest.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../misc/utils.o ../geom/topo.o ../out/output.o
	$(CXX) -o $@ $^ $(LFLAGS)

ip1test: \
    ../tests/interp1.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

ip1check: ip1test
	@ ./ip1test > /dev/null 2>&1 || echo "WARNING: P1-repair failed."

ip2test: \
    ../tests/interp2.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

prolongationp2test: \
    ../tests/prolongationp2test.o ../geom/multigrid.o ../geom/topo.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

tetrabuildertest: \
    ../tests/tetrabuildertest.o ../geom/multigrid.o ../geom/topo.o \
    ../num/unknowns.o ../geom/builder.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

mattest: \
    ../tests/mattest.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

testfe: \
    ../tests/testfe.o ../misc/utils.o ../geom/topo.o \
    ../geom/multigrid.o ../num/unknowns.o
	$(CXX) -o $@ $^ $(LFLAGS)

sbuffer: \
    ../tests/sbuffer.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

minres: \
    ../tests/minres.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

meshreader: \
    ../tests/meshreader.o ../geom/multigrid.o ../geom/topo.o \
    ../num/unknowns.o ../geom/builder.o ../geom/boundary.o ../misc/utils.o \
    ../out/output.o
	$(CXX) -o $@ $^ $(LFLAGS)

restrictp2: \
    ../tests/restrictp2.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

vectest: \
    ../tests/vectest.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

p2local: ../tests/p2local.o ../misc/utils.o \
  ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o ../num/unknowns.o \
  ../out/output.o  ../misc/problem.o ../geom/builder.o \
  ../num/discretize.o ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

quadbase: \
    ../tests/quadbase.o ../misc/utils.o ../geom/builder.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o \
    ../num/fe.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

globallist: \
    ../tests/globallist.o  ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

triang: \
    ../tests/triang.o  ../misc/utils.o ../geom/builder.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o
	$(CXX) -o $@ $^ $(LFLAGS)

quadCut: \
    ../tests/quadCut.o  ../misc/utils.o ../geom/builder.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o \
    ../levelset/levelset.o ../levelset/fastmarch.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

bicgstab: \
    ../tests/bicgstab.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

gcr: \
    ../tests/gcr.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)
blockmat: \
    ../tests/blockmat.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

mass: \
    ../tests/mass.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

quad5: \
    ../tests/quad5.o ../misc/utils.o ../geom/builder.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o \
    ../num/fe.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

downwind: \
    ../tests/downwind.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

# include common rules:
include ../Make.common
