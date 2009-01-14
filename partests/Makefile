#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = partests
EXEC = TestRefPar TestExchangePar TestMGSerPar TestInterpolPar \
       TestSedPar TestGMRESPar TestLevelsetPar \
       TestPoissonPar TestStokesPar TestInstatStokesPar TestNavStokesPar \
       TestMzellePar TestMzelleAdaptPar MzelleNMRParamEst TestBrickflowPar \
       TestFilmPar #TestRefSedPar TwoPhaseFlowPar

DELETE = $(EXEC)
CPP = $(wildcard *.cpp)
OBJ = $(CPP:%.cpp=%.o)

# rules:

default:

all: $(EXEC)

PAR_OBJ = ../parallel/parallel.o ../parallel/parmultigrid.o \
          ../parallel/partime.o ../parallel/addeddata.o ../parallel/loadbal.o \
          ../parallel/exchange.o ../parallel/memmgr_std.o ../parallel/parmgserialization.o \
          ../parallel/logger.o

TestRefPar: \
   $(PAR_OBJ) \
   ../partests/TestRefPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/fe.o ../num/interfacePatch.o ../num/unknowns.o ../out/output.o \
    ../partests/params.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestExchangePar: \
   $(PAR_OBJ) \
   ../partests/TestExchangePar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/unknowns.o ../num/fe.o ../num/interfacePatch.o \
   ../partests/params.o ../levelset/levelset.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestMGSerPar: \
   $(PAR_OBJ) \
   ../partests/TestMGSerPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/fe.o ../num/interfacePatch.o ../num/unknowns.o ../out/output.o \
   ../out/vtkOut.o ../partests/params.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestInterpolPar: \
    $(PAR_OBJ) \
    ../partests/TestInterpolPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
    ../geom/builder.o ../misc/utils.o ../levelset/levelset.o ../levelset/params.o \
    ../misc/params.o ../misc/problem.o ../num/discretize.o ../num/unknowns.o \
    ../num/fe.o ../num/interfacePatch.o ../out/output.o ../out/ensightOut.o \
    ../partests/params.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestSedPar: \
   $(PAR_OBJ) \
   ../partests/TestSedPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
    ../geom/builder.o ../misc/utils.o ../levelset/levelset.o ../levelset/params.o \
    ../misc/params.o ../misc/problem.o ../num/discretize.o ../num/unknowns.o \
    ../num/fe.o ../num/interfacePatch.o ../out/output.o ../out/ensightOut.o \
    ../out/vtkOut.o ../partests/params.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestGMRESPar: \
   $(PAR_OBJ) \
   ../partests/TestGMRESPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/unknowns.o ../num/discretize.o ../num/fe.o ../num/interfacePatch.o \
   ../out/output.o ../out/ensightOut.o ../partests/params.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestLevelsetPar: \
   $(PAR_OBJ) \
   ../partests/TestLevelsetPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o \
   ../num/unknowns.o ../num/discretize.o ../num/fe.o ../num/interfacePatch.o \
   ../levelset/fastmarch.o ../out/ensightOut.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestPoissonPar: \
   $(PAR_OBJ) \
   ../partests/TestPoissonPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/unknowns.o ../num/discretize.o ../num/fe.o ../num/interfacePatch.o \
   ../out/output.o ../out/ensightOut.o ../partests/params.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestStokesPar: \
   $(PAR_OBJ) \
   ../partests/TestStokesPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/unknowns.o ../num/discretize.o ../num/fe.o ../num/interfacePatch.o \
   ../num/stokessolver.o ../out/output.o ../out/ensightOut.o ../out/vtkOut.o\
   ../partests/params.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestInstatStokesPar: \
   $(PAR_OBJ) \
   ../partests/TestInstatStokesPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/unknowns.o ../num/discretize.o ../num/fe.o ../num/interfacePatch.o \
   ../num/stokessolver.o ../out/output.o ../out/ensightOut.o ../out/vtkOut.o\
   ../partests/params.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestNavStokesPar: \
   $(PAR_OBJ) \
   ../partests/TestNavStokesPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../levelset/levelset.o ../misc/utils.o ../misc/problem.o ../misc/params.o\
   ../num/unknowns.o ../num/discretize.o ../num/fe.o ../num/interfacePatch.o \
   ../num/stokessolver.o ../out/output.o ../out/ensightOut.o \
   ../partests/params.o  
	$(CXX) -o $@ $^ $(LFLAGS)

TestMzellePar: \
   $(PAR_OBJ) \
   ../partests/TestMzellePar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/params.o \
   ../misc/utils.o ../misc/params.o ../misc/problem.o ../num/discretize.o ../num/unknowns.o \
   ../num/fe.o ../num/stokessolver.o ../num/interfacePatch.o ../out/output.o  \
   ../out/ensightOut.o ../partests/params.o ../stokes/instatstokes2phase.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestMzelleAdaptPar: \
   $(PAR_OBJ) \
   ../partests/TestMzelleAdaptPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/params.o \
   ../misc/utils.o ../misc/params.o ../misc/problem.o ../num/discretize.o ../num/unknowns.o \
   ../num/fe.o ../num/stokessolver.o ../num/interfacePatch.o ../out/output.o  \
   ../out/ensightOut.o ../out/vtkOut.o ../partests/params.o ../stokes/instatstokes2phase.o
	$(CXX) -o $@ $^ $(LFLAGS)

MzelleNMRParamEst: \
   $(PAR_OBJ) \
   ../partests/MzelleNMRParamEst.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/params.o \
   ../misc/utils.o ../misc/params.o ../misc/problem.o ../num/discretize.o ../num/unknowns.o \
   ../num/fe.o ../num/stokessolver.o ../num/interfacePatch.o ../out/output.o  \
   ../out/ensightOut.o ../out/vtkOut.o ../partests/params.o ../stokes/instatstokes2phase.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestBrickflowPar: \
   $(PAR_OBJ) \
   ../partests/TestBrickflowPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../levelset/levelset.o ../levelset/fastmarch.o ../levelset/params.o \
   ../misc/utils.o ../misc/params.o ../misc/problem.o ../num/discretize.o ../num/unknowns.o \
   ../num/fe.o ../num/stokessolver.o ../num/interfacePatch.o ../out/output.o  \
   ../out/ensightOut.o ../out/vtkOut.o ../partests/params.o ../stokes/instatstokes2phase.o
	$(CXX) -o $@ $^ $(LFLAGS)

TestFilmPar: \
   $(PAR_OBJ) \
   ../partests/TestFilmPar.o ../geom/multigrid.o ../geom/boundary.o ../geom/topo.o \
   ../geom/builder.o ../misc/utils.o ../misc/problem.o ../misc/params.o \
   ../num/unknowns.o ../num/discretize.o ../num/fe.o ../num/interfacePatch.o \
   ../out/output.o ../out/ensightOut.o ../partests/params.o
	$(CXX) -o $@ $^ $(LFLAGS)

#TestRefSedPar: \
#   $(PAR_OBJ) \
#   ../partests/TestRefSedPar.o ../geom/multigrid.o ../geom/boundary.o ../misc/utils.o\
#   ../geom/topo.o ../geom/builder.o ../misc/params.o ../levelset/params.o ../partests/params.o \
#   ../num/unknowns.o ../misc/problem.o ../out/output.o ../out/vtkOut.o
#	$(CXX) -o $@ $^ $(LFLAGS)

#TwoPhaseFlowPar: \
#   $(PAR_OBJ) \
#   ../partests/TwoPhaseFlowPar.o ../geom/multigrid.o ../geom/boundary.o ../misc/utils.o \
#   ../num/unknowns.o ../geom/topo.o ../geom/builder.o ../misc/problem.o ../num/discretize.o\
#   ../levelset/fastmarch.o ../levelset/params.o ../partests/params.o ../misc/params.o \
#   ../levelset/levelset.o ../num/fe.o ../stokes/instatstokes2phase.o ../out/ensightOut.o \
#   ../out/vtkOut.o
#	$(CXX) -o $@ $^ $(LFLAGS)


# include common rules:
include ../Make.common

# overrides targets that are deleted by distclean

DISTDELETE = ensight/* vtk/* output/* *.case *.info
