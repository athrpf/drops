#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = poisson

EXEC = drops_statP2 poissonP1

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)

# include common rules:
include ../Make.common

PAR_OBJ_ = ../parallel/parallel.o ../parallel/parmultigrid.o \
           ../parallel/partime.o ../parallel/addeddata.o ../parallel/partitioner.o ../parallel/loadbal.o \
           ../parallel/exchange.o ../parallel/memmgr_std.o ../parallel/parddd.o

PAR_OBJ = $(if $(PAR_BUILD),$(PAR_OBJ_),)

# rules

default: ClrScr dropsrep

ipdropsAD: \
    ../poisson/poissonCoeff.o \
    ../geom/boundary.o ../geom/builder.o ../poisson/ipdropsAD.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/interfacePatch.o \
    ../poisson/poisson.o ../misc/problem.o ../misc/utils.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

poissonP1:\
    ../poisson/poissonCoeff.o \
    ../poisson/poissonP1.o ../geom/boundary.o ../geom/builder.o ../out/ensightOut.o\
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/poisson.o ../misc/problem.o ../misc/utils.o ../out/output.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o ../geom/geomselect.o\
    ../poisson/poissonparams.o ../misc/params.o ../levelset/lsetparams.o\
    ../stokes/instatstokes2phase.o ../levelset/levelset.o ../levelset/fastmarch.o\
    ../levelset/surfacetension.o ../out/vtkOut.o ../misc/bndmap.o ../geom/bndScalarFunctions.o\
    $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)
	
drops_statP2:\
    ../poisson/poissonCoeff.o \
    ../poisson/drops_statP2.o ../geom/boundary.o ../geom/builder.o\
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/poisson.o ../misc/problem.o ../misc/utils.o ../out/output.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o ../geom/geomselect.o\
    ../poisson/poissonparams.o ../misc/params.o ../levelset/lsetparams.o \
    ../stokes/instatstokes2phase.o ../levelset/levelset.o ../levelset/fastmarch.o\
    ../levelset/surfacetension.o ../out/vtkOut.o ../misc/bndmap.o ../geom/bndScalarFunctions.o\
    ../num/hypre.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

matlab/ipdrops:\
    ../poisson/poissonCoeff.o \
    ../geom/boundary.o ../geom/builder.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/interfacePatch.o \
    ../poisson/poisson.o ../misc/problem.o ../out/output.o ../misc/utils.o
	$(CXX) $(CXXFLAGS) matlab/ipdrops.cpp -o $@ $^ $(LFLAGS)

matlab/ipfilm:\
    ../poisson/poissonCoeff.o \
    ../geom/boundary.o ../geom/builder.o ../num/interfacePatch.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/poisson.o ../misc/problem.o ../out/output.o ../misc/utils.o
	$(CXX) $(CXXFLAGS) matlab/ipfilm.cpp -o $@ $^ $(LFLAGS)


clean: matlab/clean

matlab/clean:
	cd matlab && $(MAKE) clean
