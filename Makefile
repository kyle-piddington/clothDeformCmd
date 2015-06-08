all: a.out
	g++ -std=c++0x -Wno-deprecated *.cpp *.cc -DGL_GLEXT_PROTOTYPES -lglut -lGL -lGLU

cuda: clean
	nvcc -std=c++11 -pedantic -g -O3 -Wno-deprecated *.cpp *.cc -framework GLUT -framework OpenGL

osx: clean
	g++ -std=c++11 -pedantic -g -O3 -Wno-deprecated *.cpp *.cc -framework GLUT -framework OpenGL
icc: clean
	icpc -std=c++11 -O3 -openmp -pedantic -g -Wno-deprecated *.cpp *.cc -lglut -lGL -lGLU -DGL_GLEXT_PROTOTYPES
offload: clean
	icpc -std=c++11 -DOFFLOAD -O3 -qopt-report-file=cloth.optrpt -qopt-report=5 -qopt-report-phase=vec -qopt-report-routine=step     -openmp -pedantic -g -Wno-deprecated -lglut -lGL -lGLU  *.cpp *.cc -DGL_GLEXT_PROTOTYPES
clean:
	rm -rf *~ *.o a3 a.out*
