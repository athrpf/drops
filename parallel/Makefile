#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DIR = parallel
EXEC =
DELETE = $(EXEC)
CPP = $(wildcard *.cpp)
OBJ = $(CPP:%.cpp=%.o)

# rules:

default:

all: $(EXEC)

PAR_OBJ = parallel.o parmultigrid.o partime.o addeddata.o partitioner.o loadbal.o exchange.o memmgr_std.o parmgserialization.o logger.o parddd.o 

# include common rules:
include ../Make.common
