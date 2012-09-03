
CXX_FLAGS=-g -O0

all: TESTTemporal TESTSpacial TESTTemporalSHORT TESTValuesOnGridpoint TESTFineGrid TESTfindheightUpperLowerSimple TESTfindheightUpperLowerMoreSubtle findheightviafunctionSHORTALTERNATIVE 

#main: main.cpp
#	g++ -o wilko main.cpp

TESTfindheightUpperLowerMoreSubtle: TESTfindheightUpperLowerMoreSubtle.cpp PERIODICADDEDphysicaldata  
	g++ ${CXX_FLAGS} -o TESTfindheightUpperLowerMoreSubtle TESTfindheightUpperLowerSimple.cpp 
	
TESTfindheightUpperLowerSimple: TESTfindheightUpperLowerSimple.cpp PERIODICADDEDphysicaldata  
	g++ ${CXX_FLAGS} -o TESTfindheightUpperLowerSimple TESTfindheightUpperLowerSimple.cpp 
    	
TESTFineGrid: TESTheightinterpolateLiangDataFineGrid.cpp heightinterpolateLiangData  boxinterpLiangData PERIODICADDEDphysicaldata
	g++ ${CXX_FLAGS} -o TESTFineGrid TESTheightinterpolateLiangDataFineGrid.cpp heightinterpolateLiangData.o boxinterpLiangData.o	
    
TESTValuesOnGridpoint: TESTheightinterpolateLiangDataValuesOnGridpoint.cpp heightinterpolateLiangData  boxinterpLiangData PERIODICADDEDphysicaldata
	g++ ${CXX_FLAGS} -o TESTValuesOnGridpoint TESTheightinterpolateLiangDataValuesOnGridpoint.cpp heightinterpolateLiangData.o boxinterpLiangData.o
	
TESTTemporalSHORT: TESTheightinterpolateLiangDataTemporalPeriodicitySHORT.cpp heightinterpolateLiangData  boxinterpLiangData DiscreteLevelSetToDiscreteHeight createlevelLiang PERIODICADDEDphysicaldata
	g++ ${CXX_FLAGS} -o TESTTemporalSHORT TESTheightinterpolateLiangDataTemporalPeriodicity.cpp heightinterpolateLiangData.o boxinterpLiangData.o DiscreteLevelSetToDiscreteHeight.o	createlevelLiang.o
	
TESTTemporal: TESTheightinterpolateLiangDataTemporalPeriodicity.cpp heightinterpolateLiangData boxinterpLiangData PERIODICADDEDphysicaldata
	g++ ${CXX_FLAGS} -o TESTTemporal TESTheightinterpolateLiangDataTemporalPeriodicity.cpp heightinterpolateLiangData.o boxinterpLiangData.o
	
TESTSpacial: TESTheightinterpolateLiangDataSpacialPeriodicity.cpp heightinterpolateLiangData  boxinterpLiangData PERIODICADDEDphysicaldata
	g++ ${CXX_FLAGS} -o TESTSpacial TESTheightinterpolateLiangDataSpacialPeriodicity.cpp heightinterpolateLiangData.o boxinterpLiangData.o
	
	
findheightviafunctionSHORTALTERNATIVE:  findheightviafunctionSHORTALTERNATIVE.cpp DiscreteLevelSetToDiscreteHeight createlevelLiang PERIODICADDEDphysicaldata
	g++ ${CXX_FLAGS} -o findheightviafunctionSHORTALTERNATIVE findheightviafunctionSHORTALTERNATIVE.cpp DiscreteLevelSetToDiscreteHeight.o createlevelLiang.o   
	
heightinterpolateLiangData: ../heightinterpolateLiangData.hpp ../heightinterpolateLiangData.cpp 
	g++ ${CXX_FLAGS} -c ../heightinterpolateLiangData.cpp 
		
boxinterpLiangData: ../../boxinterpLiangData.hpp ../../boxinterpLiangData.cpp
	g++ ${CXX_FLAGS} -c ../../boxinterpLiangData.cpp 
	
DiscreteLevelSetToDiscreteHeight: ../DiscreteLevelSetToDiscreteHeight.hpp ../DiscreteLevelSetToDiscreteHeight.cpp 
	g++ ${CXX_FLAGS} -c ../DiscreteLevelSetToDiscreteHeight.cpp 

createlevelLiang: ../createlevelLiang.hpp ../createlevelLiang.cpp 
	g++ ${CXX_FLAGS} -c ../createlevelLiang.cpp 	
	
PERIODICADDEDphysicaldata: ../../PERIODICADDEDphysicaldata.hpp


