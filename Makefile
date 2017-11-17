#lookAtThis: main.cpp SPH_system.h	
#	g++ -o lookAtThis main.cpp -lglut -lGL -lGLU

all: main.o SPH_system.o
	g++ -o toRun main.o SPH_system.o -lglut -lGL -lGLU
main.o: main.cpp SPH_system.h
	g++ -c main.cpp -lglut -lGL -lGLU
SPH_system.o: SPH_system.cpp SPH_system.h Vector3T.h
	g++ -c SPH_system.cpp -lglut -lGL -lGLU
clean:
	rm toRun
	rm -rf *.o
