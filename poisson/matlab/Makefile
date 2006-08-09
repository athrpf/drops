#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

MC = mex

FLAGS = -f ../poisson/matlab/cxxopts.sh

MCFLAGS=$(FLAGS) -I..


# rules

default: ipdrops

ipdrops:
	cd .. && $(MAKE) ../geom/topo.o
	cd .. && $(MAKE) matlab/ipdrops CXX=$(MC) CXXFLAGS="$(MCFLAGS)"

ipfilm:
	cd .. && $(MAKE) ../geom/topo.o
	cd .. && $(MAKE) matlab/ipfilm CXX=$(MC) CXXFLAGS="$(MCFLAGS)"

clean:
	rm -f *.mexglx


DROPS_ROOT = ../..
include ../../drops.conf
