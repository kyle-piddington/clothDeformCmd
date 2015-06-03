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


typedef struct Triangle
{
   Eigen::Vector3d vertA;
   int idxA;
   Eigen::Vector3d vertB;
   int idxB;
   Eigen::Vector3d vertC;
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
   inline Eigen::Vector2d getUV(int vertIdx);
   inline Eigen::Vector3d getVert(int vertIdx);
   inline Triangle getTriangle(int triangleNumber);
   Weights precalcTriangle(Triangle t);
   float t;

   float youngMod;
   float poissonCoeff;
   float dampening;

};


#endif
