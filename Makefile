LIBRARIES= -lgmp -lmpfr -lCGAL -lglfw -lGLEW -lglut -lGL -lGLU -lboost_thread -Ipoisson-disk-sampling/include -IOBJ-Loader -Itinyobjloader
#COMPILER_FLAGS = -Wall -pedantic
COMPILER_FLAGS = -O2 -std=c++11

all: main.o Basic_SPH_System.o Particle_Generator.o OpenGL_Renderer.o Uniform_Grid.o viridis.o Maya_Interface.o
	g++ -o pbs_sim main.o viridis.o Particle_Generator.o Basic_SPH_System.o OpenGL_Renderer.o Uniform_Grid.o Maya_Interface.o $(COMPILER_FLAGS) $(LIBRARIES)
main.o: main.cpp
	g++ -c main.cpp $(COMPILER_FLAGS) $(LIBRARIES)
Basic_SPH_System.o: Basic_SPH_System.cpp Basic_SPH_System.h
	g++ -c Basic_SPH_System.cpp $(COMPILER_FLAGS) $(LIBRARIES) 
OpenGL_Renderer.o: OpenGL_Renderer.cpp Uniform_Grid.h
	g++ -c OpenGL_Renderer.cpp viridis.cpp $(COMPILER_FLAGS)  $(LIBRARIES) 
Uniform_Grid.o: Uniform_Grid.cpp Uniform_Grid.h
	g++ -c Uniform_Grid.cpp $(COMPILER_FLAGS)  $(LIBRARIES)
viridis.o: viridis.cpp viridis.h
	g++ -c viridis.cpp $(COMPILER_FLAGS)  $(LIBRARIES)
Particle_Generator.o: Particle_Generator.cpp Particle_Generator.h
	g++ -c Particle_Generator.cpp $(COMPILER_FLAGS)  $(LIBRARIES)
Maya_Interface.o: Maya_Interface.cpp Maya_Interface.h
	g++ -c Maya_Interface.cpp $(COMPILER_FLAGS)  $(LIBRARIES)
clean:
	rm -rf *.o
	rm pbs_sim
	rm frameFiles/*
