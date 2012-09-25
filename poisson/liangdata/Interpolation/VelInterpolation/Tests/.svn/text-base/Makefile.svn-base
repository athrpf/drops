
CXX_FLAGS=-g -O0

all: janTemporalPERIODICADDED janSpacialPERIODICADDED TESTValuesOnGridpointPERIODICADDED TESTFineGridPERIODICADDED
#main: main.cpp
#	g++ -o wilko main.cpp

	
janTemporalPERIODICADDED: TESTfunctionjaninterpolateLiangDataTemporalPeriodicityPERIODICADDED.cpp functionjaninterpolateLiangData boxinterpLiangData PERIODICADDEDphysicaldata GridpointCoordinatesFromGridpointNumber 
	g++ ${CXX_FLAGS} -o janTemporalPERIODICADDED TESTfunctionjaninterpolateLiangDataTemporalPeriodicityPERIODICADDED.cpp functionjaninterpolateLiangData.o boxinterpLiangData.o     GridpointCoordinatesFromGridpointNumber.o 	
	
janSpacialPERIODICADDED: TESTfunctionjaninterpolateLiangDataSpacialPeriodicityPERIODICADDED.cpp functionjaninterpolateLiangData boxinterpLiangData PERIODICADDEDphysicaldata GridpointCoordinatesFromGridpointNumber 
	g++ ${CXX_FLAGS} -o janSpacialPERIODICADDED TESTfunctionjaninterpolateLiangDataSpacialPeriodicityPERIODICADDED.cpp functionjaninterpolateLiangData.o boxinterpLiangData.o     GridpointCoordinatesFromGridpointNumber.o 		
	
TESTValuesOnGridpointPERIODICADDED: TESTfunctionjaninterpolateLiangDataValuesOnGridpointPERIODICADDED.cpp functionjaninterpolateLiangData boxinterpLiangData PERIODICADDEDphysicaldata GridpointCoordinatesFromGridpointNumber 
	g++ ${CXX_FLAGS} -o TESTValuesOnGridpointPERIODICADDED TESTfunctionjaninterpolateLiangDataValuesOnGridpointPERIODICADDED.cpp functionjaninterpolateLiangData.o boxinterpLiangData.o     GridpointCoordinatesFromGridpointNumber.o 
	
TESTFineGridPERIODICADDED: TESTfunctionjaninterpolateLiangDataFineGridWilkoOrderPERIODICADDED.cpp functionjaninterpolateLiangData boxinterpLiangData PERIODICADDEDphysicaldata  
	g++ ${CXX_FLAGS} -o TESTFineGridPERIODICADDED TESTfunctionjaninterpolateLiangDataFineGridWilkoOrderPERIODICADDED.cpp functionjaninterpolateLiangData.o boxinterpLiangData.o     GridpointCoordinatesFromGridpointNumber.o 	
	
functionjaninterpolateLiangData: ../functionjaninterpolateLiangData.hpp ../functionjaninterpolateLiangData.cpp 
	g++ ${CXX_FLAGS} -c ../functionjaninterpolateLiangData.cpp

boxinterpLiangData: ../../boxinterpLiangData.hpp ../../boxinterpLiangData.cpp
	g++ ${CXX_FLAGS} -c ../../boxinterpLiangData.cpp 

GridpointCoordinatesFromGridpointNumber: GridpointCoordinatesFromGridpointNumber.hpp GridpointCoordinatesFromGridpointNumber.cpp
	g++ ${CXX_FLAGS} -c GridpointCoordinatesFromGridpointNumber.cpp

PERIODICADDEDphysicaldata: ../../PERIODICADDEDphysicaldata.hpp
