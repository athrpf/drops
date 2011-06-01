#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = transport 

EXEC =  ns_transp

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)

include ../Make.common

PAR_OBJ_ = ../parallel/parallel.o ../parallel/parmultigrid.o \
           ../parallel/partime.o ../parallel/addeddata.o ../parallel/loadbal.o \
           ../parallel/exchange.o ../parallel/memmgr_std.o ../parallel/parmgserialization.o \
           ../parallel/logger.o

PAR_OBJ = $(if $(PAR_BUILD),$(PAR_OBJ_),)

# rules

default: ClrScr all

ns_transp: \
    ../transport/transportCoeff.o \
    ../transport/ns_transp.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/lsetparams.o \
    ../navstokes/instatnavstokes2phase.o ../misc/xfem.o ../geom/bndScalarFunctions.o \
    ../geom/bndVelFunctions.o ../misc/bndmap.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o ../num/MGsolver.o\
    ../num/fe.o ../out/ensightOut.o ../stokes/integrTime.o ../transport/transportNitsche.o \
    ../num/interfacePatch.o  ../misc/params.o ../levelset/fastmarch.o ../num/fe.o ../misc/xfem.o\
    ../out/vtkOut.o ../levelset/surfacetension.o ../geom/simplex.o ../surfactant/ifacetransp.o\
    ../geom/geomselect.o  ../levelset/twophaseutils.o ../num/hypre.o ../levelset/coupling.o \
    ../geom/principallattice.o ../geom/reftetracut.o ../geom/subtriangulation.o \
    ../num/quadrature.o ../geom/multigridgraph.o $(PAR_OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)


