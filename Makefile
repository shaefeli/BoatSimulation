#lookAtThis: main.cpp SPH_system.h	
#	g++ -o lookAtThis main.cpp -lglut -lGL -lGLU
LIBRARIES=  -lglfw -lGLEW -lglut -lGL -lGLU
COMPILER_FLAGS = -Wall -pedantic

all: main.o Basic_SPH_System.o OpenGL_Renderer.o Maya_Interface.o
	g++ -o toRun main.o Basic_SPH_System.o OpenGL_Renderer.o Maya_Interface.o $(COMPILER_FLAGS) $(LIBRARIES)
main.o: main.cpp Basic_SPH_System.h
	g++ -c main.cpp $(COMPILER_FLAGS) $(LIBRARIES)
Basic_SPH_System.o: Basic_SPH_System.cpp Basic_SPH_System.h
	g++ -c Basic_SPH_System.cpp $(COMPILER_FLAGS) $(LIBRARIES)
OpenGL_Renderer.o: OpenGL_Renderer.cpp OpenGL_Renderer.h
	g++ -c OpenGL_Renderer.cpp $(COMPILER_FLAGS)  $(LIBRARIES)
Maya_Interface.o: Maya_Interface.cpp Maya_Interface.h
	g++ -c Maya_Interface.cpp $(COMPILER_FLAGS) $(LIBRARIES)
clean:
	rm toRun
	rm -rf *.o
	rm frameFiles/*
