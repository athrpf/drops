
CXX_FLAGS=-g -O0
LIBS += -lgsl -lgslcblas -lm
all: heightsplines
	
heightsplines: heightinterpolateLiangDataSPLINESTest.cpp 
	g++ ${CXX_FLAGS} -o heightsplines heightinterpolateLiangDataSPLINESTest.cpp -lgsl -lgslcblas -lm
