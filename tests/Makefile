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
        extendP1onChild principallattice

DELETE = $(EXEC) *.out *.diff *.off *.mg *.dat

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)

# include common rules:
include ../Make.common

PAR_OBJ_ = ../parallel/parallel.o ../parallel/parmultigrid.o \
           ../parallel/partime.o ../parallel/addeddata.o ../parallel/loadbal.o \
           ../parallel/exchange.o ../parallel/memmgr_std.o ../parallel/parmgserialization.o \
           ../parallel/logger.o ../parallel/parddd.o ../parallel/partitioner.o

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
    ../tests/downwind.o ../misc/utils.o ../geom/multigrid.o ../geom/simplex.o ../geom/topo.o ../num/unknowns.o
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
    ../stokes/instatstokes2phase.o ../out/ensightOut.o ../levelset/surfacetension.o \
    ../misc/bndmap.o ../geom/bndVelFunctions.o ../geom/multigridgraph.o\
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o 
	$(CXX) -o $@ $^ $(LFLAGS)

f_Gamma: \
    ../tests/f_Gamma.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/lsetparams.o ../num/interfacePatch.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o ../geom/multigridgraph.o\
    ../levelset/surfacetension.o ../misc/bndmap.o ../geom/bndVelFunctions.o \
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o ../num/quadrature.o 
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

principallattice: \
    ../tests/principallattice.o ../misc/utils.o ../geom/principallattice.o ../num/discretize.o ../geom/topo.o \
    ../num/fe.o ../num/interfacePatch.o ../misc/problem.o ../num/unknowns.o ../geom/simplex.o \
    ../geom/multigrid.o ../geom/builder.o ../geom/boundary.o ../num/unknowns.o ../geom/reftetracut.o \
    ../geom/subtriangulation.o ../num/quadrature.o
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

# DO NOT DELETE

../tests/TestBrickflowSer.o: ../geom/multigrid.h ../geom/simplex.h
../tests/TestBrickflowSer.o: ../misc/utils.h ../geom/boundary.h
../tests/TestBrickflowSer.o: ../misc/container.h ../geom/topo.h
../tests/TestBrickflowSer.o: ../num/unknowns.h ../geom/builder.h
../tests/TestBrickflowSer.o: ../num/bndData.h ../num/solver.h ../num/spmat.h
../tests/TestBrickflowSer.o: ../num/spblockmat.h ../num/stokessolver.h
../tests/TestBrickflowSer.o: ../num/MGsolver.h ../misc/problem.h
../tests/TestBrickflowSer.o: ../num/MGsolver.tpp ../stokes/integrTime.h
../tests/TestBrickflowSer.o: ../stokes/stokes.h ../num/fe.h
../tests/TestBrickflowSer.o: ../levelset/mgobserve.h ../num/fe.tpp
../tests/TestBrickflowSer.o: ../num/discretize.h ../num/discretize.tpp
../tests/TestBrickflowSer.o: ../stokes/stokes.tpp ../levelset/adaptriang.h
../tests/TestBrickflowSer.o: ../levelset/levelset.h
../tests/TestBrickflowSer.o: ../levelset/surfacetension.h
../tests/TestBrickflowSer.o: ../levelset/levelset.tpp ../levelset/fastmarch.h
../tests/TestBrickflowSer.o: ../num/interfacePatch.h
../tests/TestBrickflowSer.o: ../num/interfacePatch.tpp ../misc/KDtree.h
../tests/TestBrickflowSer.o: ../misc/KDtree.tpp ../levelset/fastmarch.tpp
../tests/TestBrickflowSer.o: ../levelset/adaptriang.tpp ../partests/params.h
../tests/TestBrickflowSer.o: ../misc/params.h ../levelset/params.h
../tests/TestBrickflowSer.o: ../out/output.h ../out/ensightOut.h
../tests/TestBrickflowSer.o: ../navstokes/instatnavstokes2phase.h
../tests/TestBrickflowSer.o: ../stokes/instatstokes2phase.h
../tests/TestBrickflowSer.o: ../stokes/instatstokes2phase.tpp
../tests/TestBrickflowSer.o: ../navstokes/instatnavstokes2phase.tpp
../tests/TestBrickflowSer.o: ../levelset/coupling.h ../num/nssolver.h
../tests/TestBrickflowSer.o: ../levelset/coupling.tpp
../tests/TestGMRESSer.o: ../geom/multigrid.h ../geom/simplex.h
../tests/TestGMRESSer.o: ../misc/utils.h ../geom/boundary.h
../tests/TestGMRESSer.o: ../misc/container.h ../geom/topo.h ../num/unknowns.h
../tests/TestGMRESSer.o: ../geom/builder.h ../num/bndData.h ../num/solver.h
../tests/TestGMRESSer.o: ../num/spmat.h ../num/spblockmat.h ../num/fe.h
../tests/TestGMRESSer.o: ../misc/problem.h ../levelset/mgobserve.h
../tests/TestGMRESSer.o: ../num/fe.tpp ../out/output.h ../partests/params.h
../tests/TestGMRESSer.o: ../misc/params.h ../levelset/params.h
../tests/TestInterpolSer.o: ../geom/multigrid.h ../geom/simplex.h
../tests/TestInterpolSer.o: ../misc/utils.h ../geom/boundary.h
../tests/TestInterpolSer.o: ../misc/container.h ../geom/topo.h
../tests/TestInterpolSer.o: ../num/unknowns.h ../geom/builder.h
../tests/TestInterpolSer.o: ../num/bndData.h ../misc/problem.h ../num/spmat.h
../tests/TestInterpolSer.o: ../num/fe.h ../levelset/mgobserve.h ../num/fe.tpp
../tests/TestInterpolSer.o: ../out/output.h
../tests/TestSedSer.o: ../geom/multigrid.h ../geom/simplex.h ../misc/utils.h
../tests/TestSedSer.o: ../geom/boundary.h ../misc/container.h ../geom/topo.h
../tests/TestSedSer.o: ../num/unknowns.h ../geom/builder.h ../num/bndData.h
../tests/TestSedSer.o: ../partests/params.h ../misc/params.h
../tests/TestSedSer.o: ../levelset/params.h ../levelset/adaptriang.h
../tests/TestSedSer.o: ../levelset/mgobserve.h ../levelset/levelset.h
../tests/TestSedSer.o: ../num/spmat.h ../num/discretize.h ../misc/problem.h
../tests/TestSedSer.o: ../num/fe.h ../num/fe.tpp ../num/discretize.tpp
../tests/TestSedSer.o: ../num/solver.h ../num/spblockmat.h
../tests/TestSedSer.o: ../levelset/surfacetension.h ../levelset/levelset.tpp
../tests/TestSedSer.o: ../levelset/fastmarch.h ../num/interfacePatch.h
../tests/TestSedSer.o: ../num/interfacePatch.tpp ../misc/KDtree.h
../tests/TestSedSer.o: ../misc/KDtree.tpp ../levelset/fastmarch.tpp
../tests/TestSedSer.o: ../levelset/adaptriang.tpp ../out/output.h
../tests/TestSedSer.o: ../out/ensightOut.h
../tests/TestSedSer.o: ../navstokes/instatnavstokes2phase.h
../tests/TestSedSer.o: ../stokes/instatstokes2phase.h ../stokes/stokes.h
../tests/TestSedSer.o: ../stokes/stokes.tpp ../num/MGsolver.h
../tests/TestSedSer.o: ../num/MGsolver.tpp ../stokes/instatstokes2phase.tpp
../tests/TestSedSer.o: ../navstokes/instatnavstokes2phase.tpp
../tests/TestSedSer.o: ../levelset/coupling.h ../num/nssolver.h
../tests/TestSedSer.o: ../num/stokessolver.h ../stokes/integrTime.h
../tests/TestSedSer.o: ../levelset/coupling.tpp
../tests/bicgstab.o: ../num/solver.h ../misc/container.h ../misc/utils.h
../tests/bicgstab.o: ../num/spmat.h ../num/spblockmat.h
../tests/blockmat.o: ../num/spmat.h ../misc/utils.h ../misc/container.h
../tests/blockmat.o: ../num/spblockmat.h
../tests/directsolver.o: ../num/directsolver.h ../num/spmat.h ../misc/utils.h
../tests/directsolver.o: ../misc/container.h
../tests/downwind.o: ../num/spmat.h ../misc/utils.h ../misc/container.h
../tests/downwind.o: ../misc/problem.h ../geom/multigrid.h ../geom/simplex.h
../tests/downwind.o: ../geom/boundary.h ../geom/topo.h ../num/unknowns.h
../tests/downwind.o: ../geom/builder.h ../num/bndData.h
../tests/extendP1onChild.o: ../num/discretize.h ../geom/multigrid.h
../tests/extendP1onChild.o: ../geom/simplex.h ../misc/utils.h
../tests/extendP1onChild.o: ../geom/boundary.h ../misc/container.h
../tests/extendP1onChild.o: ../geom/topo.h ../num/unknowns.h
../tests/extendP1onChild.o: ../misc/problem.h ../geom/builder.h
../tests/extendP1onChild.o: ../num/bndData.h ../num/spmat.h ../num/fe.h
../tests/extendP1onChild.o: ../levelset/mgobserve.h ../num/fe.tpp
../tests/extendP1onChild.o: ../num/discretize.tpp
../tests/f_Gamma.o: ../geom/multigrid.h ../geom/simplex.h ../misc/utils.h
../tests/f_Gamma.o: ../geom/boundary.h ../misc/container.h ../geom/topo.h
../tests/f_Gamma.o: ../num/unknowns.h ../out/output.h ../misc/problem.h
../tests/f_Gamma.o: ../geom/builder.h ../num/bndData.h ../num/spmat.h
../tests/f_Gamma.o: ../stokes/instatstokes2phase.h ../stokes/stokes.h
../tests/f_Gamma.o: ../num/solver.h ../num/spblockmat.h ../num/fe.h
../tests/f_Gamma.o: ../levelset/mgobserve.h ../num/fe.tpp ../num/discretize.h
../tests/f_Gamma.o: ../num/discretize.tpp ../stokes/stokes.tpp
../tests/f_Gamma.o: ../levelset/levelset.h ../levelset/surfacetension.h
../tests/f_Gamma.o: ../levelset/levelset.tpp ../levelset/fastmarch.h
../tests/f_Gamma.o: ../num/interfacePatch.h ../num/interfacePatch.tpp
../tests/f_Gamma.o: ../misc/KDtree.h ../misc/KDtree.tpp
../tests/f_Gamma.o: ../levelset/fastmarch.tpp ../num/MGsolver.h
../tests/f_Gamma.o: ../num/MGsolver.tpp ../stokes/instatstokes2phase.tpp
../tests/f_Gamma.o: ../stokes/integrTime.h ../num/stokessolver.h
../tests/f_Gamma.o: ../levelset/coupling.h ../num/nssolver.h
../tests/f_Gamma.o: ../levelset/coupling.tpp ../levelset/params.h
../tests/f_Gamma.o: ../misc/params.h
../tests/gcr.o: ../num/solver.h ../misc/container.h ../misc/utils.h
../tests/gcr.o: ../num/spmat.h ../num/spblockmat.h
../tests/globallist.o: ../misc/container.h ../misc/utils.h
../tests/interfaceP1FE.o: ../geom/builder.h ../geom/multigrid.h
../tests/interfaceP1FE.o: ../geom/simplex.h ../misc/utils.h
../tests/interfaceP1FE.o: ../geom/boundary.h ../misc/container.h
../tests/interfaceP1FE.o: ../geom/topo.h ../num/unknowns.h ../num/bndData.h
../tests/interfaceP1FE.o: ../levelset/levelset.h ../num/spmat.h
../tests/interfaceP1FE.o: ../num/discretize.h ../misc/problem.h ../num/fe.h
../tests/interfaceP1FE.o: ../levelset/mgobserve.h ../num/fe.tpp
../tests/interfaceP1FE.o: ../num/discretize.tpp ../num/solver.h
../tests/interfaceP1FE.o: ../num/spblockmat.h ../levelset/surfacetension.h
../tests/interfaceP1FE.o: ../levelset/levelset.tpp ../levelset/fastmarch.h
../tests/interfaceP1FE.o: ../num/interfacePatch.h ../num/interfacePatch.tpp
../tests/interfaceP1FE.o: ../misc/KDtree.h ../misc/KDtree.tpp
../tests/interfaceP1FE.o: ../levelset/fastmarch.tpp
../tests/interp1.o: ../misc/utils.h ../num/spmat.h ../misc/container.h
../tests/interp1.o: ../geom/multigrid.h ../geom/simplex.h ../geom/boundary.h
../tests/interp1.o: ../geom/topo.h ../num/unknowns.h ../out/output.h
../tests/interp1.o: ../misc/problem.h ../geom/builder.h ../num/bndData.h
../tests/interp1.o: ../num/fe.h ../levelset/mgobserve.h ../num/fe.tpp
../tests/interp2.o: ../misc/utils.h ../num/spmat.h ../misc/container.h
../tests/interp2.o: ../geom/multigrid.h ../geom/simplex.h ../geom/boundary.h
../tests/interp2.o: ../geom/topo.h ../num/unknowns.h ../out/output.h
../tests/interp2.o: ../misc/problem.h ../geom/builder.h ../num/bndData.h
../tests/interp2.o: ../num/fe.h ../levelset/mgobserve.h ../num/fe.tpp
../tests/mass.o: ../num/solver.h ../misc/container.h ../misc/utils.h
../tests/mass.o: ../num/spmat.h ../num/spblockmat.h
../tests/mattest.o: ../misc/container.h ../misc/utils.h ../num/solver.h
../tests/mattest.o: ../num/spmat.h ../num/spblockmat.h
../tests/meshreader.o: ../misc/container.h ../misc/utils.h ../out/output.h
../tests/meshreader.o: ../geom/multigrid.h ../geom/simplex.h
../tests/meshreader.o: ../geom/boundary.h ../geom/topo.h ../num/unknowns.h
../tests/meshreader.o: ../misc/problem.h ../geom/builder.h ../num/bndData.h
../tests/meshreader.o: ../num/spmat.h
../tests/minres.o: ../num/solver.h ../misc/container.h ../misc/utils.h
../tests/minres.o: ../num/spmat.h ../num/spblockmat.h
../tests/neq.o: ../misc/utils.h ../num/spmat.h ../misc/container.h
../tests/neq.o: ../num/solver.h ../num/spblockmat.h
../tests/p2local.o: ../misc/utils.h ../num/spmat.h ../misc/container.h
../tests/p2local.o: ../geom/multigrid.h ../geom/simplex.h ../geom/boundary.h
../tests/p2local.o: ../geom/topo.h ../num/unknowns.h ../out/output.h
../tests/p2local.o: ../misc/problem.h ../geom/builder.h ../num/bndData.h
../tests/p2local.o: ../num/discretize.h ../num/fe.h ../levelset/mgobserve.h
../tests/p2local.o: ../num/fe.tpp ../num/discretize.tpp
../tests/prolongationp2test.o: ../misc/utils.h ../num/spmat.h
../tests/prolongationp2test.o: ../misc/container.h ../geom/multigrid.h
../tests/prolongationp2test.o: ../geom/simplex.h ../geom/boundary.h
../tests/prolongationp2test.o: ../geom/topo.h ../num/unknowns.h
../tests/prolongationp2test.o: ../out/output.h ../misc/problem.h
../tests/prolongationp2test.o: ../geom/builder.h ../num/bndData.h ../num/fe.h
../tests/prolongationp2test.o: ../levelset/mgobserve.h ../num/fe.tpp
../tests/quad5.o: ../misc/utils.h ../num/spmat.h ../misc/container.h
../tests/quad5.o: ../geom/multigrid.h ../geom/simplex.h ../geom/boundary.h
../tests/quad5.o: ../geom/topo.h ../num/unknowns.h ../out/output.h
../tests/quad5.o: ../misc/problem.h ../geom/builder.h ../num/bndData.h
../tests/quad5.o: ../num/discretize.h ../num/fe.h ../levelset/mgobserve.h
../tests/quad5.o: ../num/fe.tpp ../num/discretize.tpp
../tests/quad5_2D.o: ../misc/utils.h ../num/spmat.h ../misc/container.h
../tests/quad5_2D.o: ../geom/multigrid.h ../geom/simplex.h ../geom/boundary.h
../tests/quad5_2D.o: ../geom/topo.h ../num/unknowns.h ../geom/builder.h
../tests/quad5_2D.o: ../num/bndData.h ../num/discretize.h ../misc/problem.h
../tests/quad5_2D.o: ../num/fe.h ../levelset/mgobserve.h ../num/fe.tpp
../tests/quad5_2D.o: ../num/discretize.tpp
../tests/quadCut.o: ../geom/builder.h ../geom/multigrid.h ../geom/simplex.h
../tests/quadCut.o: ../misc/utils.h ../geom/boundary.h ../misc/container.h
../tests/quadCut.o: ../geom/topo.h ../num/unknowns.h ../num/bndData.h
../tests/quadCut.o: ../levelset/levelset.h ../num/spmat.h ../num/discretize.h
../tests/quadCut.o: ../misc/problem.h ../num/fe.h ../levelset/mgobserve.h
../tests/quadCut.o: ../num/fe.tpp ../num/discretize.tpp ../num/solver.h
../tests/quadCut.o: ../num/spblockmat.h ../levelset/surfacetension.h
../tests/quadCut.o: ../levelset/levelset.tpp ../levelset/fastmarch.h
../tests/quadCut.o: ../num/interfacePatch.h ../num/interfacePatch.tpp
../tests/quadCut.o: ../misc/KDtree.h ../misc/KDtree.tpp
../tests/quadCut.o: ../levelset/fastmarch.tpp
../tests/quadbase.o: ../misc/utils.h ../num/spmat.h ../misc/container.h
../tests/quadbase.o: ../geom/multigrid.h ../geom/simplex.h ../geom/boundary.h
../tests/quadbase.o: ../geom/topo.h ../num/unknowns.h ../out/output.h
../tests/quadbase.o: ../misc/problem.h ../geom/builder.h ../num/bndData.h
../tests/quadbase.o: ../num/discretize.h ../num/fe.h ../levelset/mgobserve.h
../tests/quadbase.o: ../num/fe.tpp ../num/discretize.tpp
../tests/reftest.o: ../geom/multigrid.h ../geom/simplex.h ../misc/utils.h
../tests/reftest.o: ../geom/boundary.h ../misc/container.h ../geom/topo.h
../tests/reftest.o: ../num/unknowns.h ../geom/builder.h ../num/bndData.h
../tests/reftest.o: ../out/output.h ../misc/problem.h ../num/spmat.h
../tests/reparam.o: ../levelset/levelset.h ../num/spmat.h ../misc/utils.h
../tests/reparam.o: ../misc/container.h ../num/discretize.h
../tests/reparam.o: ../geom/multigrid.h ../geom/simplex.h ../geom/boundary.h
../tests/reparam.o: ../geom/topo.h ../num/unknowns.h ../misc/problem.h
../tests/reparam.o: ../geom/builder.h ../num/bndData.h ../num/fe.h
../tests/reparam.o: ../levelset/mgobserve.h ../num/fe.tpp
../tests/reparam.o: ../num/discretize.tpp ../num/solver.h ../num/spblockmat.h
../tests/reparam.o: ../levelset/surfacetension.h ../levelset/levelset.tpp
../tests/reparam.o: ../levelset/fastmarch.h ../num/interfacePatch.h
../tests/reparam.o: ../num/interfacePatch.tpp ../misc/KDtree.h
../tests/reparam.o: ../misc/KDtree.tpp ../levelset/fastmarch.tpp
../tests/reparam.o: ../levelset/adaptriang.h ../levelset/adaptriang.tpp
../tests/reparam.o: ../levelset/params.h ../misc/params.h ../out/ensightOut.h
../tests/reparam.o: ../out/vtkOut.h
../tests/reparam_init.o: ../misc/utils.h ../geom/multigrid.h
../tests/reparam_init.o: ../geom/simplex.h ../geom/boundary.h
../tests/reparam_init.o: ../misc/container.h ../geom/topo.h ../num/unknowns.h
../tests/reparam_init.o: ../geom/builder.h ../num/bndData.h
../tests/reparam_init.o: ../levelset/levelset.h ../num/spmat.h
../tests/reparam_init.o: ../num/discretize.h ../misc/problem.h ../num/fe.h
../tests/reparam_init.o: ../levelset/mgobserve.h ../num/fe.tpp
../tests/reparam_init.o: ../num/discretize.tpp ../num/solver.h
../tests/reparam_init.o: ../num/spblockmat.h ../levelset/surfacetension.h
../tests/reparam_init.o: ../levelset/levelset.tpp ../levelset/fastmarch.h
../tests/reparam_init.o: ../num/interfacePatch.h ../num/interfacePatch.tpp
../tests/reparam_init.o: ../misc/KDtree.h ../misc/KDtree.tpp
../tests/reparam_init.o: ../levelset/fastmarch.tpp ../out/ensightOut.h
../tests/reparam_init.o: ../levelset/adaptriang.h ../levelset/adaptriang.tpp
../tests/restrictp2.o: ../misc/utils.h ../num/spmat.h ../misc/container.h
../tests/restrictp2.o: ../geom/multigrid.h ../geom/simplex.h
../tests/restrictp2.o: ../geom/boundary.h ../geom/topo.h ../num/unknowns.h
../tests/restrictp2.o: ../out/output.h ../misc/problem.h ../geom/builder.h
../tests/restrictp2.o: ../num/bndData.h ../num/fe.h ../levelset/mgobserve.h
../tests/restrictp2.o: ../num/fe.tpp
../tests/sbuffer.o: ../misc/container.h ../misc/utils.h
../tests/serialization.o: ../misc/container.h ../misc/utils.h
../tests/serialization.o: ../geom/multigrid.h ../geom/simplex.h
../tests/serialization.o: ../geom/boundary.h ../geom/topo.h ../num/unknowns.h
../tests/serialization.o: ../out/output.h ../misc/problem.h ../geom/builder.h
../tests/serialization.o: ../num/bndData.h ../num/spmat.h
../tests/splitboundary.o: ../geom/multigrid.h ../geom/simplex.h
../tests/splitboundary.o: ../misc/utils.h ../geom/boundary.h
../tests/splitboundary.o: ../misc/container.h ../geom/topo.h
../tests/splitboundary.o: ../num/unknowns.h ../geom/builder.h
../tests/splitboundary.o: ../num/bndData.h ../out/output.h ../misc/problem.h
../tests/splitboundary.o: ../num/spmat.h
../tests/testfe.o: ../num/fe.h ../misc/container.h ../misc/utils.h
../tests/testfe.o: ../geom/topo.h ../geom/multigrid.h ../geom/simplex.h
../tests/testfe.o: ../geom/boundary.h ../num/unknowns.h ../misc/problem.h
../tests/testfe.o: ../geom/builder.h ../num/bndData.h ../num/spmat.h
../tests/testfe.o: ../levelset/mgobserve.h ../num/fe.tpp
../tests/tetrabuildertest.o: ../geom/builder.h ../geom/multigrid.h
../tests/tetrabuildertest.o: ../geom/simplex.h ../misc/utils.h
../tests/tetrabuildertest.o: ../geom/boundary.h ../misc/container.h
../tests/tetrabuildertest.o: ../geom/topo.h ../num/unknowns.h
../tests/tetrabuildertest.o: ../num/bndData.h
../tests/triang.o: ../misc/utils.h ../geom/multigrid.h ../geom/simplex.h
../tests/triang.o: ../geom/boundary.h ../misc/container.h ../geom/topo.h
../tests/triang.o: ../num/unknowns.h ../geom/builder.h ../num/bndData.h
../tests/vectest.o: ../num/spmat.h ../misc/utils.h ../misc/container.h
../tests/xfem.o: ../misc/utils.h ../num/spmat.h ../misc/container.h
../tests/xfem.o: ../geom/multigrid.h ../geom/simplex.h ../geom/boundary.h
../tests/xfem.o: ../geom/topo.h ../num/unknowns.h ../out/output.h
../tests/xfem.o: ../misc/problem.h ../geom/builder.h ../num/bndData.h
../tests/xfem.o: ../num/discretize.h ../num/fe.h ../levelset/mgobserve.h
../tests/xfem.o: ../num/fe.tpp ../num/discretize.tpp ../levelset/levelset.h
../tests/xfem.o: ../num/solver.h ../num/spblockmat.h
../tests/xfem.o: ../levelset/surfacetension.h ../levelset/levelset.tpp
../tests/xfem.o: ../levelset/fastmarch.h ../num/interfacePatch.h
../tests/xfem.o: ../num/interfacePatch.tpp ../misc/KDtree.h
../tests/xfem.o: ../misc/KDtree.tpp ../levelset/fastmarch.tpp
../tests/xfem.o: ../stokes/instatstokes2phase.h ../stokes/stokes.h
../tests/xfem.o: ../stokes/stokes.tpp ../num/MGsolver.h ../num/MGsolver.tpp
../tests/xfem.o: ../stokes/instatstokes2phase.tpp ../out/ensightOut.h
