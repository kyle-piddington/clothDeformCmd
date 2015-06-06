/**
 * A horribly ineficient shape for a triangle mesh.
 * Seriously it's awful, don't use it.
 *
 */
#ifndef __Cloth_H___
#define __Cloth_H___
#ifdef __APPLE__
#include <GLUT/glut.h>
#endif
#ifdef __unix__
#include <GL/glut.h>
#endif
#include <iostream>
#include "GLSL.h"
#include "glm/glm.hpp"
#include <vector>



class Cloth
{
public:
   Cloth(int w, int h, int res);
   ~Cloth();
   void bind();
   void init();
   void draw(GLint h_pos, GLint h_nor, GLint h_tex);
   void step(float time);
   glm::vec2 getUV(int vertIdx);
   glm::vec3 getVert(int vertIdx);
   int getNumVerts(){return verts.size();}
   int getNumTriangles(){return numTriangles;}
   std::vector<int> getInds(){return inds;}
private:
   std::vector<float> verts;
   std::vector<float> norms;
   std::vector<float> tex;
   std::vector<int> inds;
   int numTriangles;
   GLuint posBufID;
   GLuint norBufID;
   GLuint indBufID;
   GLuint texBufID;
   void rebindNorms();
   void rebindVerts();


};


#endif
