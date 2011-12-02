#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DROPS_ROOT=../..

DIR = poisson/pydrops

EXEC = convdiff

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr

#    ../../poisson/poissonCoeff.o \

convdiff:\
    ../../geom/boundary.o ../../geom/builder.o ../../out/ensightOut.o\
    ../../geom/simplex.o ../../geom/multigrid.o ../../num/unknowns.o ../../geom/topo.o \
    ../../poisson/poissonCoeff.o\
    ../../poisson/poisson.o ../../misc/problem.o ../../misc/utils.o ../../out/output.o \
    ../../num/fe.o ../../num/discretize.o ../../num/interfacePatch.o ../../geom/geomselect.o\
    ../../misc/params.o \
    ../../out/vtkOut.o ../../misc/bndmap.o ../../geom/bndScalarFunctions.o\
    ../../geom/bndVelFunctions.o
	$(CXX) $(CXXFLAGS) -o convdiff $^ convection_diffusion.cpp

# include common rules:
include $(DROPS_ROOT)/Make.common
