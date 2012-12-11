
CXX_FLAGS=-g -O0
LIBS += -lgsl -lgslcblas -lm
all: heightsplines spacialperiodicity temporalperiodicity testplotSPLINEFineGrid
	
heightsplines: heightinterpolateLiangDataSPLINESTest.cpp 
	g++ ${CXX_FLAGS} -o heightsplines heightinterpolateLiangDataSPLINESTest.cpp -lgsl -lgslcblas -lm
	
spacialperiodicity: heightinterpolateLiangDataSPLINESTestSpacialPeriodicity.cpp 
	g++ ${CXX_FLAGS} -o spacialperiodicity heightinterpolateLiangDataSPLINESTestSpacialPeriodicity.cpp -lgsl -lgslcblas -lm
	
temporalperiodicity: heightinterpolateLiangDataSPLINESTestTemporalPeriodicity.cpp 
	g++ ${CXX_FLAGS} -o temporalperiodicity heightinterpolateLiangDataSPLINESTestTemporalPeriodicity.cpp -lgsl -lgslcblas -lm
	
testplotSPLINEFineGrid: heightinterpolateLiangDataSPLINESTestPLOTFineGrid.cpp 
	g++ ${CXX_FLAGS} -o testplotSPLINEFineGrid heightinterpolateLiangDataSPLINESTestPLOTFineGrid.cpp -lgsl -lgslcblas -lm
