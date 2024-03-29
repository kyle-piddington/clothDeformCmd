//
// sueda
// October, 2014
//

#ifdef __APPLE__
#include <GLUT/glut.h>
#endif
#ifdef __unix__
#include <GL/glut.h>
#endif
#include <iostream>
#include "GLSL.h"
#include "Camera.h"
#include "Shape.h"
#include "MatrixStack.h"
#include "Shader.h"
#include "Material.h"
#include "Light.h"

#include "glm/gtx/rotate_vector.hpp"
#include "Program.h"
#include "Cloth.h"
#include <memory>

#define TIMESTEP 1/60.0
using namespace std;
bool keyToggles[256];
bool keyDowns[256];
double t = 0, dt = 0;
Camera camera;
Program prog_phong;
Program prog_debug;


Cloth testCloth(16,16,64);

Material defaultMaterial = Material(
		glm::vec3(0.2,0.2,0.2),
		glm::vec3(0.8,0.3,0.3),
		glm::vec3(0.0,0.0,0.0),
		10.0);
std::vector<Light> defaultLights;





void initGL()
{
	//////////////////////////////////////////////////////
	// Initialize GL for the whole scene
	//
	// Set background color
	glClearColor(0.4f, 0.4f, 0.4f, 0.4f);
	// Enable z-buffer test
	glEnable(GL_DEPTH_TEST);




	//////////////////////////////////////////////////////
	// Initialize shaders
	//

	// Create shader handles

	prog_phong.setShaderNames("shaders/passthrough_vert.glsl", "shaders/phong_frag.glsl");
	prog_phong.init();
	prog_phong.addAttribute("vertPos");
	prog_phong.addAttribute("vertNor");
	//prog_phong.addAttribute("vertTexCoords");
	prog_phong.addUniform("P");
	prog_phong.addUniform("MV");
	prog_phong.addUniform("NORM");
	prog_phong.addUniform("lightPos");
	prog_phong.addUniform("lightIntensity");
	prog_phong.addUniform("numLights");
	prog_phong.addUniform("ka");
	prog_phong.addUniform("kd");
	prog_phong.addUniform("ks");
	prog_phong.addUniform("s");


	prog_debug.setShaderNames("shaders/debug_passthrough.glsl","shaders/debug_shade.glsl");
	prog_debug.init();
	prog_debug.addAttribute("vertPos");
	prog_debug.addUniform("P");
	prog_debug.addUniform("MV");

	defaultLights.push_back(Light(glm::vec3(0.0,1.0,1.0),0.8));
	defaultLights.push_back(Light(glm::vec3(-1.0,1.0,1.0),0.2));
	GLSL::checkVersion();
	testCloth.init();
	testCloth.bind();
	
   camera.mouseClicked(0,-1200,true,false,false);
   camera.mouseClicked(600, 600,false,false,false);
}

void reshapeGL(int w, int h)
{
	// Set view size
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	// Set camera aspect ratio
	camera.setAspect((float)w/h);
}

void update(float dt)
{
	testCloth.step(dt);
}
void drawDebug()
{
	MatrixStack P, MV;
	// Apply camera transforms
	P.pushMatrix();
	camera.applyProjectionMatrix(&P);
	MV.pushMatrix();
	camera.applyCameraMatrix(&MV);
	prog_debug.bind();
	glUniformMatrix4fv(prog_debug.getUniform("P"), 1, GL_FALSE, glm::value_ptr( P.topMatrix()));
	glUniformMatrix4fv(prog_debug.getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV.topMatrix()));
	

	//Debug Draws

}
void drawGL()
{
	// Clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Create matrix stacks
	MatrixStack P, MV;
	// Apply camera transforms
	P.pushMatrix();
	camera.applyProjectionMatrix(&P);
	MV.pushMatrix();
	camera.applyCameraMatrix(&MV);


	// Bind the program
	prog_phong.bind();

	glUniformMatrix4fv(prog_phong.getUniform("P"), 1, GL_FALSE, glm::value_ptr( P.topMatrix()));
	glUniformMatrix4fv(prog_phong.getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV.topMatrix()));
	glUniformMatrix4fv(prog_phong.getUniform("NORM"),1,GL_FALSE,
		glm::value_ptr(glm::transpose(glm::inverse(MV.topMatrix()))));



	defaultMaterial.bindMaterial(prog_phong.getUniform("ka"),
										  prog_phong.getUniform("kd"),
										  prog_phong.getUniform("ks"),
										  prog_phong.getUniform("s"));

	/*lighting*/
	glm::vec3 *lPos = (glm::vec3 *)malloc(sizeof(glm::vec3)*defaultLights.size());
	float *lInt = (float *)malloc(sizeof(float) * defaultLights.size());


	for(int i = 0; i < defaultLights.size(); i++){
		lPos[i] = defaultLights[i].position;
		lInt[i] = defaultLights[i].intensity;
	}

	glUniform3fv(prog_phong.getUniform("lightPos"), defaultLights.size(), (const float *)lPos);
	glUniform1fv(prog_phong.getUniform("lightIntensity"),defaultLights.size(), (const float *)lInt);
	glUniform1i(prog_phong.getUniform("numLights"),defaultLights.size());




	free(lPos);
	free(lInt);

	// Draw shape

	testCloth.draw(prog_phong.getAttribute("vertPos"),
				  prog_phong.getAttribute("vertNor"),
				  -1);




	// Unbind the program
	prog_phong.unbind();
	// Pop stacks
	MV.popMatrix();
	P.popMatrix();

	if(keyToggles['d'])
	{
		drawDebug();
	}
	// Double buffer
	glutSwapBuffers();
}

void mouseGL(int button, int state, int x, int y)
{
	int modifier = glutGetModifiers();
	bool shift = modifier & GLUT_ACTIVE_SHIFT;
	bool ctrl  = modifier & GLUT_ACTIVE_CTRL;
	bool alt   = modifier & GLUT_ACTIVE_ALT;
	camera.mouseClicked(x, y, shift, ctrl, alt);
}

void motionGL(int x, int y)
{
	camera.mouseMoved(x, y);
	glutPostRedisplay();
}

void keyboardUpGL(unsigned char key, int x, int y)
{
   keyDowns[key] = false;
}

void keyboardGL(unsigned char key, int x, int y)
{
	keyToggles[key] = !keyToggles[key];
	keyDowns[key] = true;
   switch(key) {
		case 27:
			// ESCAPE
			exit(0);
			break;
		case 's':
			
			update(TIMESTEP);
			break;
		case 'k':
			testCloth.kickCenter();
			update(TIMESTEP);
			break;
      case 'o':

         if(keyToggles['o'])
         {
            testCloth.startOffload();
         }
         else
         {
            testCloth.endOffload();
         }
         break;
	}
   if(keyToggles['l'])
	{
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
	}
 
	else
	{
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	}
	// Refresh screen
	glutPostRedisplay();
}

void idleGL()
{
	static int frame = 0, time, timebase = 0;	// for fps
	
	frame++; // calc fps
   
   if(keyDowns['e'])
   {
      testCloth.expand(-0.1);
   }
   else if(keyDowns['w'])
   {
      testCloth.expand(0.1);
   }
   float oldTime = glutGet(GLUT_ELAPSED_TIME);
	if(keyToggles[' '])
	{
			update(TIMESTEP);
	}
	glutPostRedisplay();
	
	// calc fps
	time = glutGet(GLUT_ELAPSED_TIME);

	//printf("update: %lf\n", (time - oldTime));			
	if (time - timebase > 1000)
	{
		double fps = frame * 1000.0 / (time - timebase);
		timebase = time;
		frame = 0;
		
		printf("fps: %lf\n", fps);
	}
}

int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitWindowSize(400, 400);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Zach Arend");
	glutMouseFunc(mouseGL);
	glutMotionFunc(motionGL);
	glutKeyboardFunc(keyboardGL);
   glutKeyboardUpFunc(keyboardUpGL);
   glutReshapeFunc(reshapeGL);
	glutDisplayFunc(drawGL);
	glutIdleFunc(idleGL);
	initGL();
	glutMainLoop();
	return 0;
}
