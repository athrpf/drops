#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = levelset

EXEC = reparam lsshear surfTens prJump film brick_transp twophasedrops

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr all

reparam: \
    ../levelset/reparam.o ../levelset/fastmarch.o ../levelset/levelset.o \
    ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o ../out/ensightOut.o
	$(CXX) -o $@ $^ $(LFLAGS)

lsshear: \
    ../levelset/lsshear.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/fastmarch.o \
    ../num/fe.o ../out/ensightOut.o ../stokes/instatstokes2phase.o ../misc/xfem.o ../num/MGsolver.o
	$(CXX) -o $@ $^ $(LFLAGS)

surfTens: \
    ../levelset/surfTens.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o \
    ../num/fe.o ../out/ensightOut.o ../stokes/integrTime.o ../misc/xfem.o ../num/MGsolver.o
	$(CXX) -o $@ $^ $(LFLAGS)

prJump: \
    ../levelset/prJump.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o ../num/MGsolver.o \
    ../num/fe.o  ../out/ensightOut.o ../stokes/integrTime.o ../misc/xfem.o
	$(CXX) -o $@ $^ $(LFLAGS)

film: \
    ../levelset/film.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o ../num/MGsolver.o\
    ../num/fe.o ../out/ensightOut.o ../stokes/integrTime.o ../misc/xfem.o 
	$(CXX) -o $@ $^ $(LFLAGS)

brick_transp: \
    ../levelset/brick_transp.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o \
    ../misc/params.o ../levelset/fastmarch.o ../num/fe.o ../misc/xfem.o\
    ../poisson/transport2phase.o ../out/ensightOut.o ../out/vtkOut.o
	$(CXX) -o $@ $^ $(LFLAGS)

twophasedrops: \
    ../levelset/twophasedrops.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o ../misc/xfem.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o ../num/MGsolver.o\
    ../num/fe.o ../out/ensightOut.o ../stokes/integrTime.o ../poisson/transport2phase.o
	$(CXX) -o $@ $^ $(LFLAGS)

# include common rules:
include ../Make.common
