#############################################
#   D R O P S   local makefile              #
#############################################

# variables:

DROPS_ROOT = ..

DIR = misc

EXEC = compensight

DELETE = $(EXEC)

CPP = $(wildcard *.cpp)

OBJ = $(CPP:%.cpp=%.o)


# rules

default: ClrScr all

compensight: \
    ../misc/compensight.o ../misc/utils.o
	$(CXX) -o $@ $^ $(LFLAGS)


# include common rules:
include ../Make.common
