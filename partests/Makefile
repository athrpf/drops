#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = partests
EXEC = TestRefPar TestExchangePar TestPoissonPar TestInterpolPar TestSedPar
#       TestStokesPar TestInstatStokesPar \
#       TestMzellePar TestMzelleAdaptPar MzelleNMRParamEst TestBrickflowPar \
#       TestFilmPar
       
# compile errors; not converted to Ensight6OutCL:       TestLevelsetPar
# out of date: TestSedPar TestNavStokesPar TwoPhaseFlowPar

DELETE = $(EXEC)
CPP = $(wildcard *.cpp)
OBJ = $(CPP:%.cpp=%.o)

# include common rules:
include ../Make.common

# rules:

default:

all: $(EXEC)

PAR_OBJ_ = ../parallel/parallel.o ../parallel/parmultigrid.o \
          ../parallel/partime.o ../parallel/addeddata.o ../parallel/partitioner.o ../parallel/loadbal.o \
          ../parallel/exchange.o ../parallel/memmgr_std.o ../parallel/parmgserialization.o \
          ../parallel/logger.o ../parallel/parddd.o

PAR_OBJ = $(if $(PAR_BUILD),$(PAR_OBJ_),)

TestRefPar: \
   $(PAR_OBJ) \
   ../partests/TestRefPar.o ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/fe.o ../num/interfacePatch.o ../num/unknowns.o ../out/output.o \
    ../partests/params.o ../levelset/lsetparams.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestExchangePar: \
   $(PAR_OBJ) \
   ../partests/TestExchangePar.o ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/unknowns.o ../num/fe.o ../num/interfacePatch.o ../levelset/fastmarch.o \
   ../partests/params.o ../levelset/levelset.o ../num/discretize.o ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestPoissonPar: \
   $(PAR_OBJ) \
   ../partests/TestPoissonPar.o ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/unknowns.o ../num/discretize.o ../num/fe.o ../num/interfacePatch.o \
   ../out/output.o ../out/ensightOut.o ../partests/params.o ../levelset/lsetparams.o ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestInterpolPar: \
    $(PAR_OBJ) \
    ../partests/TestInterpolPar.o ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
    ../geom/builder.o ../misc/utils.o ../levelset/levelset.o ../levelset/lsetparams.o \
    ../misc/params.o ../misc/problem.o ../num/discretize.o ../num/unknowns.o \
    ../num/fe.o ../num/interfacePatch.o ../out/output.o ../out/ensightOut.o \
    ../partests/params.o ../levelset/fastmarch.o ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestSedPar: \
    ../partests/TestSedPar.o ../geom/boundary.o ../geom/builder.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/lsetparams.o ../geom/geomselect.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o ../num/MGsolver.o\
    ../num/fe.o ../out/ensightOut.o ../stokes/integrTime.o ../poisson/transport2phase.o ../levelset/twophaseutils.o\
    ../num/interfacePatch.o ../out/vtkOut.o ../surfactant/ifacetransp.o ../levelset/surfacetension.o ../geom/simplex.o \
    ../num/hypre.o ../levelset/coupling.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

# perform tests
###############

test_all: test_exchange test_poisson test_ref test_interpol test_sed
	echo "All tests performed"

test_exchange: TestExchangePar
	echo "Testing ExchangeCL ..."
	$(ARCH_RUN) -np 4 ./TestExchangePar param-files/Exchange.param > ex_tmp.txt
	diff ex_tmp.txt ref-out/Exchange_P004.txt ; rm -f ex_tmp.txt

test_poisson: TestPoissonPar
	echo "Testing Poisson ..."
	$(ARCH_RUN) -np 4 ./TestPoissonPar param-files/Poisson.param > poisson_tmp.txt
	diff poisson_tmp.txt ref-out/Poisson_P004.txt ; rm -f poisson_tmp.txt

test_ref: TestRefPar
	echo "Testing Refinement ... "
	$(ARCH_RUN) -np 4 ./TestRefPar param-files/Ref.param > ref_tmp.txt
	diff ref_tmp.txt ref-out/Ref_P004.txt ; rm -f ref_tmp.txt

test_interpol: TestInterpolPar
	echo "Testing Interpolation ... "
	$(ARCH_RUN) -np 4 ./TestInterpolPar > interpol_tmp.txt
	diff interpol_tmp.txt ref-out/Interpol_P004.txt  ; rm -f interpol_tmp.txt

test_sed: TestSedPar
	echo "Testing sedimenting droplet ..."
	$(ARCH_RUN) -np 4 ./TestSedPar param-files/Sed.param | sed -re "s/(.* )([0-9\.]+)( seconds| sec\.)/\1XXX\3/" > sed_tmp.txt
	diff sed_tmp.txt ref-out/Sed_P004.txt ; rm -f sed_tmp.txt


# deprecated:
#############

TestStokesPar: \
   $(PAR_OBJ) \
   ../partests/TestStokesPar.o ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/unknowns.o ../num/discretize.o ../num/fe.o ../num/interfacePatch.o \
   ../num/stokessolver.o ../out/output.o ../out/ensightOut.o ../out/vtkOut.o\
   ../partests/params.o ../levelset/lsetparams.o ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestInstatStokesPar: \
   $(PAR_OBJ) \
   ../partests/TestInstatStokesPar.o ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/unknowns.o ../num/discretize.o ../num/fe.o ../num/interfacePatch.o \
   ../num/stokessolver.o ../out/output.o ../out/ensightOut.o ../out/vtkOut.o\
   ../partests/params.o ../levelset/lsetparams.o ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestNavStokesPar: \
   $(PAR_OBJ) \
   ../partests/TestNavStokesPar.o ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../levelset/levelset.o ../misc/utils.o ../misc/problem.o ../misc/params.o\
   ../num/unknowns.o ../num/discretize.o ../num/fe.o ../num/interfacePatch.o \
   ../num/stokessolver.o ../out/output.o ../out/ensightOut.o \
   ../partests/params.o  ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestMzellePar: \
   $(PAR_OBJ) \
   ../partests/TestMzellePar.o ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/lsetparams.o \
   ../misc/utils.o ../misc/params.o ../misc/problem.o ../num/discretize.o ../num/unknowns.o \
   ../num/fe.o ../num/stokessolver.o ../num/interfacePatch.o ../out/output.o  \
   ../out/ensightOut.o ../partests/params.o ../stokes/instatstokes2phase.o ../stokes/integrTime.o \
   ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestMzelleAdaptPar: \
   $(PAR_OBJ) \
   ../partests/TestMzelleAdaptPar.o ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/lsetparams.o \
   ../misc/utils.o ../misc/params.o ../misc/problem.o ../num/discretize.o ../num/unknowns.o \
   ../num/fe.o ../num/stokessolver.o ../num/interfacePatch.o ../out/output.o \
   ../out/ensightOut.o ../out/vtkOut.o ../partests/params.o ../stokes/instatstokes2phase.o \
   ../stokes/integrTime.o ../misc/xfem.o ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

MzelleNMRParamEst: \
   $(PAR_OBJ) \
   ../partests/MzelleNMRParamEst.o ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/lsetparams.o \
   ../misc/utils.o ../misc/params.o ../misc/problem.o ../num/discretize.o ../num/unknowns.o \
   ../num/fe.o ../num/stokessolver.o ../num/interfacePatch.o ../out/output.o  \
   ../out/ensightOut.o ../out/vtkOut.o ../partests/params.o ../stokes/instatstokes2phase.o \
   ../stokes/integrTime.o ../misc/xfem.o ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestBrickflowPar: \
   $(PAR_OBJ) \
   ../partests/TestBrickflowPar.o ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/lsetparams.o \
   ../misc/utils.o ../misc/params.o ../misc/problem.o ../num/discretize.o ../num/unknowns.o \
   ../num/fe.o ../num/stokessolver.o ../num/interfacePatch.o ../out/output.o  \
   ../out/ensightOut.o ../out/vtkOut.o ../partests/params.o ../stokes/instatstokes2phase.o \
   ../stokes/integrTime.o ../misc/xfem.o ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestFilmPar: \
   $(PAR_OBJ) \
   ../partests/TestFilmPar.o ../geom/simplex.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/unknowns.o ../num/discretize.o ../num/fe.o ../num/interfacePatch.o \
   ../out/output.o ../out/ensightOut.o ../partests/params.o ../levelset/surfacetension.o
	$(CXX) -o $@ $^ $(LFLAGS)


