all:
	g++ -std=c++11 -Wno-deprecated  *.cpp *.cc -DGL_GLEXT_PROTOTYPES -lglut -lGL -lGLU

cuda: clean
	nvcc *.cpp *.cc *.cu -Xlinker -framework,OpenGL,-framework,GLUT
cudaLinux: clean
	nvcc -std=c++11 -g  *.cpp *.cc *.cu -Xlinker -lGL -lGLU -DGL_GLEXT_PROTOTYPES -lglut 
osx: clean
	g++ -std=c++11 -pedantic -g -O3 -Wno-deprecated *.cpp *.cc -framework GLUT -framework OpenGL
icc: clean
	icpc -std=c++11 -O3 -openmp -pedantic -g -Wno-deprecated *.cpp *.cc -lglut -lGL -lGLU -DGL_GLEXT_PROTOTYPES
offload: clean
	icpc -std=c++11 -DOFFLOAD -O3 -openmp -pedantic -g -Wno-deprecated -lglut -lGL -lGLU  *.cpp *.cc -DGL_GLEXT_PROTOTYPES
clean:
	rm -rf *~ *.o a.out*
