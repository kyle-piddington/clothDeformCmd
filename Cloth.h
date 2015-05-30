#ifndef __CLOTH_H__
#define __CLOTH_H__
#include "glm/glm.hpp"
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

typedef struct Vertex
{
   glm::vec3 position;
   glm::vec3 velocity;
}Vertex;


typedef struct Triangle
{
   Vertex a;
   Vertex b;
   Vertex c;
}Triangle;


class Cloth
{
public:
   Cloth(glm::vec2 size, glm::vec2 resolution);
   ~Cloth(){}
   void precalc();
   void init();
   void step(float dt);
   void bind();
   void draw(GLint h_pos, GLint h_nor);
 
private:

   GLuint posBufID;
   GLuint norBufID;

   int numTriangles;
   std::vector<glm::vec3> vertex_position;
   std::vector<glm::vec3> vertex_velocity;
   std::vector<glm::vec2> weights;
   std::vector<float> d;
   std::vector<int> indices;

   

};
#endif
