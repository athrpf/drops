
CXX_FLAGS=-g -O0
LIBS += -lgsl -lgslcblas -lm
all: heightsplines spacialperiodicity
	
heightsplines: heightinterpolateLiangDataSPLINESTest.cpp 
	g++ ${CXX_FLAGS} -o heightsplines heightinterpolateLiangDataSPLINESTest.cpp -lgsl -lgslcblas -lm
	
spacialperiodicity: heightinterpolateLiangDataSPLINESTestSpacialPeriodicity.cpp 
	g++ ${CXX_FLAGS} -o spacialperiodicity heightinterpolateLiangDataSPLINESTestSpacialPeriodicity.cpp -lgsl -lgslcblas -lm
