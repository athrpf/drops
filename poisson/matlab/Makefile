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
	cd ..; make ../geom/topo.o
	cd ..; make matlab/ipdrops CC=$(MC) CFLAGS="$(MCFLAGS)"
  
clean:
	rm -f *.mexglx 
  

DROPS_ROOT = ../..
include ../../drops.conf
