#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = poisson

EXEC = drops MGdrops testdrops dropsrep ipdrops ipdropsrep drops_neumann \
       MGdropsP2 ipfilm drops_statP2

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)

# include common rules:
include ../Make.common

PAR_OBJ_ = ../parallel/parallel.o ../parallel/parmultigrid.o \
           ../parallel/partime.o ../parallel/addeddata.o ../parallel/loadbal.o \
           ../parallel/exchange.o ../parallel/memmgr_std.o

PAR_OBJ = $(if $(PAR_BUILD),$(PAR_OBJ_),)

# rules

default: ClrScr dropsrep

drops: \
    ../geom/boundary.o ../geom/builder.o ../poisson/drops.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../poisson/poisson.o ../num/interfacePatch.o\
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

dropsrep: \
    ../geom/boundary.o ../geom/builder.o ../poisson/dropsrep.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../poisson/poisson.o ../num/interfacePatch.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

drops_neumann: \
    ../geom/boundary.o ../geom/builder.o ../poisson/drops_neumann.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../poisson/poisson.o ../num/interfacePatch.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

MGdrops: \
    ../geom/boundary.o ../geom/builder.o ../poisson/MGdrops.o ../geom/simplex.o ../geom/multigrid.o \
    ../out/output.o ../num/unknowns.o ../geom/topo.o ../num/interfacePatch.o \
    ../num/MGsolver.o ../poisson/poisson.o ../misc/problem.o ../misc/utils.o \
    ../num/fe.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

MGdropsP2: \
    ../geom/boundary.o ../geom/builder.o ../poisson/MGdropsP2.o ../geom/simplex.o ../geom/multigrid.o \
    ../out/output.o ../num/unknowns.o ../geom/topo.o ../num/interfacePatch.o \
    ../num/MGsolver.o ../poisson/poisson.o ../misc/problem.o ../misc/utils.o \
    ../num/fe.o  ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

testdrops: \
    ../geom/boundary.o ../geom/builder.o ../poisson/testdrops.o ../geom/simplex.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/MGsolver.o ../poisson/poisson.o \
    ../misc/problem.o ../misc/utils.o ../out/output.o ../num/discretize.o \
    ../num/fe.o ../num/interfacePatch.o
	$(CXX) -o $@ $^ $(LFLAGS)

ipdrops: \
    ../geom/boundary.o ../geom/builder.o ../poisson/ipdrops.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/interfacePatch.o \
    ../num/MGsolver.o ../poisson/instatpoisson.o ../misc/problem.o ../misc/utils.o \
    ../num/fe.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

ipdropsAD: \
    ../geom/boundary.o ../geom/builder.o ../poisson/ipdropsAD.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/interfacePatch.o \
    ../poisson/instatpoisson.o ../misc/problem.o ../misc/utils.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

ipdropsrep: \
    ../geom/boundary.o ../geom/builder.o ../poisson/ipdropsrep.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/interfacePatch.o \
    ../poisson/instatpoisson.o ../misc/problem.o ../misc/utils.o ../num/fe.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

ipfilm: \
    ../geom/boundary.o ../geom/builder.o ../poisson/ipfilm.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/interfacePatch.o \
    ../poisson/instatpoisson.o ../misc/problem.o ../misc/utils.o ../out/ensightOut.o \
    ../num/fe.o ../num/discretize.o
	$(CXX) -o $@ $^ $(LFLAGS)

drops_statP2:\
    ../poisson/drops_statP2.o ../geom/boundary.o ../geom/builder.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/poisson.o ../misc/problem.o ../misc/utils.o \
    ../num/fe.o ../num/discretize.o ../num/interfacePatch.o \
    $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)


matlab/ipdrops: ../geom/boundary.o ../geom/builder.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/interfacePatch.o \
    ../poisson/instatpoisson.o ../misc/problem.o ../out/output.o ../misc/utils.o
	$(CXX) $(CXXFLAGS) matlab/ipdrops.cpp -o $@ $^ $(LFLAGS)

matlab/ipfilm: ../geom/boundary.o ../geom/builder.o ../num/interfacePatch.o \
    ../geom/simplex.o ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o \
    ../poisson/instatpoisson.o ../misc/problem.o ../out/output.o ../misc/utils.o
	$(CXX) $(CXXFLAGS) matlab/ipfilm.cpp -o $@ $^ $(LFLAGS)


clean: matlab/clean

matlab/clean:
	cd matlab && $(MAKE) clean
