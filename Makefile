all:
	g++ -std=c++0x -Wno-deprecated *.cpp *.cc -DGL_GLEXT_PROTOTYPES -lglut -lGL -lGLU

osx:
	g++ -std=c++11 -pedantic -g -O3 -Wno-deprecated *.cpp *.cc -framework GLUT -framework OpenGL
icc:
	icpc -std=c++11 -O3 -openmp -pedantic -g -Wno-deprecated *.cpp *.cc -lglut -lGL -lGLU -DGL_GLEXT_PROTOTYPES
offload:
	icpc -std=c++11 -DOFFLOAD -O3 -openmp -pedantic -g -Wno-deprecated -lglut -lGL -lGLU  *.cpp *.cc -DGL_GLEXT_PROTOTYPES
clean:
	rm -f *~ *.o a3
