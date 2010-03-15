#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = tests

EXEC =  ip1test ip2test mattest prolongationp2test reftest testfe \
        tetrabuildertest sbuffer minres meshreader restrictp2 vectest \
        p2local quadbase globallist triang quadCut bicgstab gcr blockmat \
        mass quad5 downwind quad5_2D interfaceP1FE serialization xfem \
        directsolver f_Gamma neq splitboundary reparam_init reparam \
        extendP1onChild

DELETE = $(EXEC) *.out *.diff *.off *.mg *.dat

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)

# include common rules:
include ../Make.common

PAR_OBJ_ = ../parallel/parallel.o ../parallel/parmultigrid.o \
           ../parallel/partime.o ../parallel/addeddata.o ../parallel/loadbal.o \
           ../parallel/exchange.o ../parallel/memmgr_std.o ../parallel/parmgserialization.o \
           ../parallel/logger.o ../parallel/parddd.o

PAR_OBJ = $(if $(PAR_BUILD),$(PAR_OBJ_),)

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
    ../tests/reftest.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../misc/utils.o ../geom/topo.o ../out/output.o ../misc/problem.o \
    ../num/interfacePatch.o ../num/fe.o 
	$(CXX) -o $@ $^ $(LFLAGS)

ip1test: \
    ../tests/interp1.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

ip1check: ip1test
	@ ./ip1test > /dev/null 2>&1 || echo "WARNING: P1-repair failed."

ip2test: \
    ../tests/interp2.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

prolongationp2test: \
    ../tests/prolongationp2test.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

tetrabuildertest: \
    ../tests/tetrabuildertest.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o \
    ../num/unknowns.o ../geom/builder.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

mattest: \
    ../tests/mattest.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

testfe: \
    ../tests/testfe.o ../misc/utils.o ../geom/topo.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o
	$(CXX) -o $@ $^ $(LFLAGS)

sbuffer: \
    ../tests/sbuffer.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

minres: \
    ../tests/minres.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

meshreader: \
    ../tests/meshreader.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o \
    ../num/unknowns.o ../geom/builder.o ../geom/boundary.o ../misc/utils.o \
    ../out/output.o ../misc/problem.o ../num/interfacePatch.o ../num/fe.o 
	$(CXX) -o $@ $^ $(LFLAGS)

restrictp2: \
    ../tests/restrictp2.o ../geom/simplex.o ../geom/multigrid.o ../geom/topo.o ../out/output.o \
    ../num/unknowns.o ../geom/builder.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../geom/boundary.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

vectest: \
    ../tests/vectest.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

p2local: ../tests/p2local.o ../misc/utils.o \
  ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o ../num/unknowns.o \
  ../out/output.o  ../misc/problem.o ../geom/builder.o ../num/interfacePatch.o \
  ../num/discretize.o ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

quadbase: \
    ../tests/quadbase.o ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

globallist: \
    ../tests/globallist.o  ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

serialization: \
    ../tests/serialization.o  ../misc/utils.o  ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../out/output.o ../num/fe.o ../misc/problem.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

triang: \
    ../tests/triang.o  ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o
	$(CXX) -o $@ $^ $(LFLAGS)

quadCut: \
    ../tests/quadCut.o  ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o ../num/interfacePatch.o \
    ../levelset/levelset.o ../levelset/fastmarch.o ../num/discretize.o ../num/fe.o ../levelset/surfacetension.o
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
    ../tests/quad5.o ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

downwind: \
    ../tests/downwind.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

quad5_2D: \
    ../tests/quad5_2D.o ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

interfaceP1FE: \
    ../tests/interfaceP1FE.o  ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o ../num/interfacePatch.o \
    ../levelset/levelset.o ../levelset/fastmarch.o ../num/discretize.o ../num/fe.o ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

xfem: \
    ../tests/xfem.o ../misc/utils.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../geom/boundary.o ../geom/topo.o ../num/unknowns.o ../misc/problem.o ../num/interfacePatch.o \
    ../num/fe.o ../num/discretize.o ../levelset/levelset.o ../levelset/fastmarch.o \
    ../stokes/instatstokes2phase.o ../out/ensightOut.o ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

f_Gamma: \
    ../tests/f_Gamma.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/lsetparams.o ../num/interfacePatch.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o \
    ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

neq: \
    ../tests/neq.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)

extendP1onChild: \
	../tests/extendP1onChild.o ../num/discretize.o ../misc/utils.o ../geom/topo.o ../num/fe.o ../misc/problem.o \
	../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

splitboundary: \
    ../tests/splitboundary.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../misc/utils.o ../geom/topo.o ../out/output.o ../num/fe.o ../misc/problem.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

reparam_init: \
    ../tests/reparam_init.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../misc/utils.o ../geom/topo.o ../out/output.o \
    ../misc/problem.o ../levelset/levelset.o ../num/fe.o \
    ../misc/utils.o ../num/discretize.o ../num/interfacePatch.o \
    ../levelset/fastmarch.o ../levelset/surfacetension.o ../out/ensightOut.o
	$(CXX) -o $@ $^ $(LFLAGS)

reparam: \
    ../tests/reparam.o ../levelset/fastmarch.o ../levelset/levelset.o \
    ../geom/simplex.o ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../out/ensightOut.o ../num/interfacePatch.o ../levelset/surfacetension.o \
    ../levelset/params.o ../misc/params.o ../levelset/lsetparams.o ../out/vtkOut.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

../tests/directsolver.o:
ifneq ($(DSOLVER_CXXFLAGS),)
	$(CXX) $(CXXFLAGS) $(DSOLVER_CXXFLAGS) -o $*.o -c $*.cpp
else
	@echo "compilation skipped due to missing CHOLMOD/UMFPACK installation"
endif

ifneq ($(DSOLVER_LFLAGS),)
directsolver: \
    ../tests/directsolver.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS) $(DSOLVER_LFLAGS)
else
directsolver:
	@echo "compilation skipped due to missing CHOLMOD/UMFPACK installation"
endif

