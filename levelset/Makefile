#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = levelset

EXEC = lsdrops reparam lsdrivcav lsshear curvDriven risingBubble \
    risingBubbleAdap mzelle_stat mzelle_instat mzelle_instat_ns \
    surfTens prJump film brickflow mzelle_ns_adap

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr all

lsdrops: \
    ../levelset/lsdrops.o ../levelset/fastmarch.o ../levelset/levelset.o \
    ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

reparam: \
    ../levelset/reparam.o ../levelset/fastmarch.o ../levelset/levelset.o \
    ../geom/multigrid.o ../geom/builder.o ../geom/topo.o ../geom/boundary.o \
    ../num/unknowns.o ../misc/utils.o ../misc/problem.o ../num/discretize.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

lsdrivcav: \
    ../levelset/lsdrivcav.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/fastmarch.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

lsshear: \
    ../levelset/lsshear.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/fastmarch.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

curvDriven: \
    ../levelset/curvDriven.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/fastmarch.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

risingBubble: \
    ../levelset/risingBubble.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/params.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/fastmarch.o ../levelset/params.o \
    ../num/fe.o ../stokes/instatstokes2phase.o
	$(CXX) -o $@ $^ $(LFLAGS)

risingBubbleAdap: \
    ../levelset/risingBubbleAdap.o ../geom/boundary.o ../geom/builder.o ../levelset/levelset.o \
    ../geom/multigrid.o ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/params.o ../misc/problem.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/fastmarch.o ../levelset/params.o ../stokes/instatstokes2phase.o
	$(CXX) -o $@ $^ $(LFLAGS)

mzelle_stat: \
    ../levelset/mzelle_stat.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/fastmarch.o \
    ../num/fe.o ../stokes/instatstokes2phase.o
	$(CXX) -o $@ $^ $(LFLAGS)

mzelle_instat: \
    ../levelset/mzelle_instat.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o \
    ../misc/params.o ../levelset/fastmarch.o ../num/MGsolver.o \
    ../num/fe.o ../stokes/instatstokes2phase.o
	$(CXX) -o $@ $^ $(LFLAGS)

mzelle_instat_ns: \
    ../levelset/mzelle_instat_ns.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

mzelle_ns_adap: \
    ../levelset/mzelle_ns_adap.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

brick_ns_adap: \
    ../levelset/brick_ns_adap.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

brickflow: \
    ../levelset/brickflow.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

surfTens: \
    ../levelset/surfTens.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

prJump: \
    ../levelset/prJump.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)

film: \
    ../levelset/film.o ../geom/boundary.o ../geom/builder.o ../geom/multigrid.o \
    ../num/unknowns.o ../geom/topo.o ../num/fe.o ../misc/problem.o ../levelset/levelset.o \
    ../misc/utils.o ../out/output.o ../num/discretize.o ../levelset/params.o \
    ../misc/params.o ../levelset/fastmarch.o ../stokes/instatstokes2phase.o \
    ../num/fe.o
	$(CXX) -o $@ $^ $(LFLAGS)



# include common rules:
include ../Make.common
