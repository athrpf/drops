CXX_FLAGS=-g -O0
LIBS += -lgsl -lgslcblas -lm
all: errorsplines error
	
errorsplines: ErrorEstimationBalakotaiahSinusTESTSPLINESsimple.cpp 
	g++ ${CXX_FLAGS} -o errorsplines ErrorEstimationBalakotaiahSinusTESTSPLINESsimple.cpp -lgsl -lgslcblas -lm
	
	
error: ErrorEstimationBalakotaiahSinusTEST.cpp 
	g++ ${CXX_FLAGS} -o error ErrorEstimationBalakotaiahSinusTEST.cpp 
