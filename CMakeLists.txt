cmake_minimum_required(VERSION 3.0)
#
# One need to set GLFW, GLM, GLEW paths manually to make it work for sure
#
project(PBS_project)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

FIND_PACKAGE(OpenGL 4.0 REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

find_package(GLUT)

IF(APPLE)
	INCLUDE_DIRECTORIES ( /System/Library/Frameworks )
	FIND_LIBRARY(COCOA_LIBRARY 				Cocoa)
	FIND_LIBRARY(IOKIT_LIBRARY 				IOKit)
	FIND_LIBRARY(COREVIDEO_LIBRARY 			CoreVideo)
	FIND_LIBRARY(COREGRAPHICS_LIBRARY 		CoreGraphics)
	FIND_LIBRARY(COREFOUNDATION_LIBRARY 	CoreFoundation)

	MARK_AS_ADVANCED (COCOA_LIBRARY
			GLUT_LIBRARY
			OpenGL_LIBRARY)

	SET(EXTRA_LIBS  ${COCOA_LIBRARY}
					${GLUT_LIBRARY}
					${OPENGL_LIBRARIES}
					${IOKIT_LIBRARY}
					${COREVIDEO_LIBRARY}
					${COREGRAPHICS_LIBRARY}
					${COREFOUNDATION_LIBRARY})

	include_directories("/usr/local/Cellar/glm/0.9.8.5/include")
	include_directories("/usr/local/include/GLFW/")
	find_library(GLFWlib glfw3)
ENDIF (APPLE)
include_directories(${GLUT_INCLUDE_DIRS})
find_package(GLEW REQUIRED)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/poisson-disk-sampling/include")

set( projFiles
		"main.cpp"
		"Basic_SPH_System.cpp"
		"Basic_SPH_System.h"
		"Uniform_Grid.cpp"
		"Uniform_Grid.h"
		"OpenGL_Renderer.cpp"
		"Vector3T.h"
		)

file( COPY "fragment_shader.fs"
		DESTINATION "bin")
file( COPY "vertex_shader.vs"
		DESTINATION "bin")
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(PBS_project ${projFiles})
target_link_libraries(PBS_project ${EXTRA_LIBS} ${GLFWlib} ${GLEW_LIBRARIES})
target_include_directories(PBS_project PUBLIC ${OPENGL_INCLUDE_DIR})