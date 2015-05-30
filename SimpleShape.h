/**
 * A horribly ineficient shape for a triangle mesh.
 * Seriously it's awful, don't use it.
 *
 */
#ifndef __SimpleShape_H___
#define __SimpleShape_H___
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

class Triangle
{
public:
   Triangle(
      glm::vec3 v1, glm::vec3 v2, glm::vec3 v3,
      glm::vec3 n1, glm::vec3 n2, glm::vec3 n3)
   {
      verts.push_back(v1);
      verts.push_back(v2);
      verts.push_back(v3);
      norms.push_back(n1);
      norms.push_back(n2);
      norms.push_back(n3);
   }
   std::vector<glm::vec3> verts;
   std::vector<glm::vec3> norms;
};


class SimpleShape
{
public:
   SimpleShape();
   ~SimpleShape();
   void addTriangle(Triangle tris);
   void bind();
   void draw(GLint h_pos, GLint h_nor);
private:
   std::vector<float> verts;
   std::vector<float> norms;
   std::vector<int> inds;
   GLuint posBufID;
   GLuint norBufID;
   GLuint indBufID;
};


#endif
