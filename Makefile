LIBRARIES= -lgmp -lmpfr -lCGAL -lglfw -lGLEW -lglut -lGL -lGLU -lboost_thread -Ipoisson-disk-sampling/include -IOBJ-Loader -Itinyobjloader
#COMPILER_FLAGS = -Wall -pedantic
COMPILER_FLAGS =  -frounding-math -O2

all: main.o Basic_SPH_System.o Particle_Generator.o OpenGL_Renderer.o Uniform_Grid.o viridis.o Maya_Interface.o PCI_SPH.o
	g++ -std=c++0x -o pbs_sim main.o viridis.o Particle_Generator.o PCI_SPH.o Basic_SPH_System.o OpenGL_Renderer.o Uniform_Grid.o Maya_Interface.o $(COMPILER_FLAGS) $(LIBRARIES) -pg -no-pie
main.o: main.cpp
	g++ -std=c++0x -c main.cpp $(COMPILER_FLAGS) $(LIBRARIES)
Basic_SPH_System.o: Basic_SPH_System.cpp Basic_SPH_System.h
	g++ -std=c++0x -c Basic_SPH_System.cpp $(COMPILER_FLAGS) $(LIBRARIES) -pg  -no-pie
PCI_SPH.o: PCI_SPH.cpp PCI_SPH.h
	g++ -std=c++0x -c PCI_SPH.cpp $(COMPILER_FLAGS) $(LIBRARIES) -pg  -no-pie
OpenGL_Renderer.o: OpenGL_Renderer.cpp Uniform_Grid.h
	g++ -std=c++0x -c OpenGL_Renderer.cpp viridis.cpp $(COMPILER_FLAGS)  $(LIBRARIES) -pg  -no-pie
Uniform_Grid.o: Uniform_Grid.cpp Uniform_Grid.h
	g++ -std=c++0x -c Uniform_Grid.cpp $(COMPILER_FLAGS)  $(LIBRARIES) -pg  -no-pie
viridis.o: viridis.cpp viridis.h
	g++ -std=c++0x -c viridis.cpp $(COMPILER_FLAGS)  $(LIBRARIES) -pg  -no-pie
Particle_Generator.o: Particle_Generator.cpp Particle_Generator.h
	g++ -std=c++0x -c Particle_Generator.cpp $(COMPILER_FLAGS)  $(LIBRARIES) -pg  -no-pie
Maya_Interface.o: Maya_Interface.cpp Maya_Interface.h
	g++ -std=c++0x -c Maya_Interface.cpp $(COMPILER_FLAGS)  $(LIBRARIES) -pg  -no-pie
clean:
	rm -rf *.o
	rm pbs_sim
	rm frameFiles/*
	rm infoFile
