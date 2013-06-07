#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DROPS_ROOT = ..

DIR = poisson

EXEC = balapoissonP1

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)

LIBS += -lgsl -lgslcblas -lm

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
    ../poisson/poisson.o ../misc/problem.o ../misc/utils.o ../misc/bndmap.o \
    ../num/fe.o $(CXX) -o $@ $^ $(LFLAGS)

interpolate: ../poisson/interpolation/boxinterp.o\
	../poisson/interpolation/createlevelset.o\
	../poisson/interpolation/heightinterpolate.o\
	../poisson/interpolation/levelsettoheight.o\
	../poisson/interpolation/velocityinterpolation.o


balapoissonP1: ../poisson/interpolation/boxinterp.o\
	../poisson/interpolation/createlevelset.o\
	../poisson/interpolation/heightinterpolate.o\
	../poisson/interpolation/levelsettoheight.o\
	../poisson/interpolation/velocityinterpolation.o\
	../poisson/poissonCoeff.o \
	../poisson/balakotaiah.o \
	../poisson/fullDNSBala.o \
	../poisson/ale.o \
	../poisson/poissonP1.o ../geom/boundary.o ../geom/builder.o ../out/ensightOut.o\
	../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
	../poisson/poisson.o ../misc/problem.o ../misc/utils.o ../out/output.o \
	../num/fe.o ../num/discretize.o ../num/interfacePatch.o ../geom/geomselect.o\
	../misc/params.o \
	../out/vtkOut.o ../misc/bndmap.o ../geom/bndScalarFunctions.o\
	../geom/bndVelFunctions.o\
	$(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS) $(LIBS)

DNSpoissonP1:\
    ../poisson/poissonCoeff.o \
    ../poisson/balakotaiah.o \
    ../poisson/fullDNSBala.o \
    ../poisson/ale.o \
    ../poisson/poissonP1.o ../geom/boundary.o ../geom/builder.o ../out/ensightOut.o\
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/poisson.o ../misc/problem.o ../misc/utils.o ../out/output.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o ../geom/geomselect.o\
    ../misc/params.o \
    ../out/vtkOut.o ../misc/bndmap.o ../geom/bndScalarFunctions.o\
    ../geom/bndVelFunctions.o\
    $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS) $(LIBS)

poissonP2:\
    ../geom/boundary.o ../geom/builder.o ../out/ensightOut.o\
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/poisson.o ../misc/problem.o ../misc/utils.o ../out/output.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o ../geom/geomselect.o\
    ../misc/params.o \
    ../out/vtkOut.o ../misc/bndmap.o ../geom/bndScalarFunctions.o\
    ../geom/bndVelFunctions.o\
    ../poisson/poissonCoeff.o \    ../poisson/poissonP2.o \
    ../num/hypre.o $(PAR_OBJ)   $(CXX) -o $@ $^ $(LFLAGS)

matlab/ipdrops:\
    ../poisson/poissonCoeff.o \
    ../geom/boundary.o ../geom/builder.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/interfacePatch.o \
    ../poisson/poisson.o ../misc/problem.o ../out/output.o ../misc/utils.o  $(CXX) $(CXXFLAGS) matlab/ipdrops.cpp -o $@ $^ $(LFLAGS)

matlab/ipfilm:\
    ../poisson/poissonCoeff.o \
    ../geom/boundary.o ../geom/builder.o ../num/interfacePatch.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/poisson.o ../misc/problem.o ../out/output.o ../misc/utils.o  $(CXX) $(CXXFLAGS) matlab/ipfilm.cpp -o $@ $^ $(LFLAGS)


clean: matlab/clean pydrops/clean

matlab/clean:   cd matlab && $(MAKE) clean

pydrops/clean:  cd pydrops && $(MAKE) clean

check-syntax:       gcc -o nul -S ${CHK_SOURCES}
