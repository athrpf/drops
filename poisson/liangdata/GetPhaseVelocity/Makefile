
CXX_FLAGS=-g -O0
LIBS += -lgsl -lgslcblas -lm
all: phasevelocity phasevelocitySPLINES
	
phasevelocity: get_phase_velocity.cpp 
	g++ ${CXX_FLAGS} -o phasevelocity get_phase_velocity.cpp -lgsl -lgslcblas -lm
	
phasevelocitySPLINES: get_phase_velocitySPLINES.cpp 
	g++ ${CXX_FLAGS} -o phasevelocitySPLINES get_phase_velocitySPLINES.cpp -lgsl -lgslcblas -lm
