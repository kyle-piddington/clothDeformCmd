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
#include <vector>
#include "Eigen/Sparse"
#include "Eigen/Geometry"

//Prototype integrator definition.
class ClothForceIntegrator;

class Cloth
{
public:
   Cloth(float w, float h, int res);
   ~Cloth();
   void bind();
   void init();
   void draw(GLint h_pos, GLint h_nor, GLint h_tex);
   void step(double time);
   Eigen::Vector2d getUV(int vertIdx);
   Eigen::Vector3d getVert(int vertIdx);
   int getNumVerts(){return verts.size()/3;}
   int getNumTriangles(){return numTriangles;}
   std::vector<int> getInds(){return inds;}
   void kickCenter();
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
   ClothForceIntegrator * integrator;
   void rebindNorms();
   void rebindVerts();
   int res;   

};


#endif
