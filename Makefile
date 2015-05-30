all:
	g++ -std=c++11 -Wno-deprecated *.cpp *.cc -DGL_GLEXT_PROTOTYPES -lglut -lGL -lGLU 

osx:
	g++ -std=c++11 -pedantic -g -Wno-deprecated *.cpp *.cc -framework GLUT -framework OpenGL 
clean:
	rm -f *~ *.o a3
