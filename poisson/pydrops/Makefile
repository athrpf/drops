#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DROPS_ROOT=../..

DIR = poisson/pydrops

EXEC = convdiff test_convdiff prepy_product test_iaprob

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr

#    ../../poisson/poissonCoeff.o \

convdiffcpp:  \
    ../../geom/boundary.o ../../geom/builder.o ../../out/ensightOut.o\
    ../../geom/simplex.o ../../geom/multigrid.o ../../num/unknowns.o ../../geom/topo.o \
    ../../poisson/poissonCoeff.o\
    ../../poisson/poisson.o ../../misc/problem.o ../../misc/utils.o ../../out/output.o \
    ../../num/fe.o ../../num/discretize.o ../../num/interfacePatch.o ../../geom/geomselect.o\
    ../../misc/params.o \
    ../../out/vtkOut.o ../../misc/bndmap.o ../../geom/bndScalarFunctions.o\
    ../../geom/bndVelFunctions.o
	$(CXX) $(CXXFLAGS) -o convdiff $^ cppconvection.cpp

convection_diffusion.o: \
    ../../geom/boundary.o ../../geom/builder.o ../../out/ensightOut.o\
    ../../geom/simplex.o ../../geom/multigrid.o ../../num/unknowns.o ../../geom/topo.o \
    ../../poisson/poissonCoeff.o\
    ../../poisson/poisson.o ../../misc/problem.o ../../misc/utils.o ../../out/output.o \
    ../../num/fe.o ../../num/discretize.o ../../num/interfacePatch.o ../../geom/geomselect.o\
    ../../misc/params.o \
    ../../out/vtkOut.o ../../misc/bndmap.o ../../geom/bndScalarFunctions.o\
    ../../geom/bndVelFunctions.o
	$(CXX) $(CXXFLAGS)  -c convection_diffusion.cpp $^

test_convdiff:\
    ../../geom/boundary.o ../../geom/builder.o ../../out/ensightOut.o\
    ../../geom/simplex.o ../../geom/multigrid.o ../../num/unknowns.o ../../geom/topo.o \
    ../../poisson/poissonCoeff.o\
    ../../poisson/poisson.o ../../misc/problem.o ../../misc/utils.o ../../out/output.o \
    ../../num/fe.o ../../num/discretize.o ../../num/interfacePatch.o ../../geom/geomselect.o\
    ../../misc/params.o \
    ../../out/vtkOut.o ../../misc/bndmap.o ../../geom/bndScalarFunctions.o\
    ../../geom/bndVelFunctions.o
	$(CXX) $(CXXFLAGS) -o test_convdiff $^ test_convdiff.cpp
    
test_iaprob:\
    ../../geom/boundary.o ../../geom/builder.o ../../out/ensightOut.o\
    ../../geom/simplex.o ../../geom/multigrid.o ../../num/unknowns.o ../../geom/topo.o \
    ../../poisson/poissonCoeff.o\
    ../../poisson/poisson.o ../../misc/problem.o ../../misc/utils.o ../../out/output.o \
    ../../num/fe.o ../../num/discretize.o ../../num/interfacePatch.o ../../geom/geomselect.o\
    ../../misc/params.o \
    ../../out/vtkOut.o ../../misc/bndmap.o ../../geom/bndScalarFunctions.o\
    ../../geom/bndVelFunctions.o
	$(CXX) $(CXXFLAGS) -o test_iaprob $^ test_iaprob.cpp

prepy_product:\
    ../../geom/boundary.o ../../geom/builder.o ../../out/ensightOut.o\
    ../../geom/simplex.o ../../geom/multigrid.o ../../num/unknowns.o ../../geom/topo.o \
    ../../poisson/poissonCoeff.o\
    ../../poisson/poisson.o ../../misc/problem.o ../../misc/utils.o ../../out/output.o \
    ../../num/fe.o ../../num/discretize.o ../../num/interfacePatch.o ../../geom/geomselect.o\
    ../../misc/params.o \
    ../../out/vtkOut.o ../../misc/bndmap.o ../../geom/bndScalarFunctions.o\
    ../../geom/bndVelFunctions.o
	$(CXX) $(CXXFLAGS) -o prepy_product $^ prepy_product.cpp

# include common rules:
include $(DROPS_ROOT)/Make.common
