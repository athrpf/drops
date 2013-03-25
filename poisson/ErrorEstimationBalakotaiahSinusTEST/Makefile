CXX_FLAGS=-g -O0
LIBS += -lgsl -lgslcblas -lm
all: errorsplines errorOutputInput errorsplinesFileInput error fileout errorBalaTEST
	
errorsplines: ErrorEstimationBalakotaiahSinusTESTSPLINESsimple.cpp 
	g++ ${CXX_FLAGS} -o errorsplines ErrorEstimationBalakotaiahSinusTESTSPLINESsimple.cpp -lgsl -lgslcblas -lm
	
errorOutputInput: ErrorTESTSPLINESFileOutputInput.cpp 
	g++ ${CXX_FLAGS} -o errorOutputInput ErrorTESTSPLINESFileOutputInput.cpp -lgsl -lgslcblas -lm
	
errorsplinesFileInput: ErrorEstimationBalakotaiahSinusTESTSPLINESFileInput.cpp 
	g++ ${CXX_FLAGS} -o errorsplinesFileInput ErrorEstimationBalakotaiahSinusTESTSPLINESFileInput.cpp -lgsl -lgslcblas -lm
	
	
error: ErrorEstimationBalakotaiahSinusTEST.cpp 
	g++ ${CXX_FLAGS} -o error ErrorEstimationBalakotaiahSinusTEST.cpp 
	
	
fileout: enumerated_file_output_test.cpp 
	g++ ${CXX_FLAGS} -o fileout enumerated_file_output_test.cpp 
	
errorBalaTEST: ErrorEstimationBalakotaiahSPLINESFileInputTEST.cpp 
	g++ ${CXX_FLAGS} -o errorBalaTEST ErrorEstimationBalakotaiahSPLINESFileInputTEST.cpp -lgsl -lgslcblas -lm
