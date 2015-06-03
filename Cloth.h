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


typedef struct Triangle
{
   glm::vec3 vertA;
   int idxA;
   glm::vec3 vertB;
   int idxB;
   glm::vec3 vertC;
   int idxC;
}Triangle;


typedef struct Weights{
   float d;
   float ua;
   float va;
   float ub;
   float vb;
   float uc;
   float vc;

}Weights;


class Cloth
{
public:
   Cloth(int w, int h, int res, float youngMod = 1000.0f, float poissonCoeff = 0.0f, float dampening = 0.0f);
   ~Cloth();
   void bind();
   void init();
   void draw(GLint h_pos, GLint h_nor, GLint h_tex);
   void step(float time);
private:
   std::vector<float> verts;
   std::vector<float> velocities;
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
   void precalc();

   //weights calculated in precalc
   std::vector<Triangle> triList;
   std::vector<Weights> triWeights;
   //Cheaper to allocate here than every frame
   std::vector<float> forces;
   //Inline helper methods
   inline glm::vec2 getUV(int vertIdx);
   inline glm::vec3 getVert(int vertIdx);
   inline Triangle getTriangle(int triangleNumber);
   Weights precalcTriangle(Triangle t);
   void integrateForces(float dt);
   float t;

   float youngMod;
   float poissonCoeff;
   float dampening;

};


#endif
