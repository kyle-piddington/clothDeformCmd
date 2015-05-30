#include "SimpleShape.h"

SimpleShape::SimpleShape() :
posBufID(0),
norBufID(0),
indBufID(0)
{

}
SimpleShape::~SimpleShape()
{

}

void SimpleShape::addTriangle(Triangle tris)
{
   for (std::vector<glm::vec3>::iterator i = tris.verts.begin(); i != tris.verts.end(); ++i)
   {
      verts.push_back(i->x);
      verts.push_back(i->y);
      verts.push_back(i->z);
   }

   for (std::vector<glm::vec3>::iterator i = tris.norms.begin(); i != tris.norms.end(); ++i)
   {
      norms.push_back(i->x);
      norms.push_back(i->y);
      norms.push_back(i->z);
   }

   inds.push_back(inds.size());
   inds.push_back(inds.size());
   inds.push_back(inds.size());

}

void SimpleShape::bind()
{
   glGenBuffers(1, &posBufID);
   glBindBuffer(GL_ARRAY_BUFFER, posBufID);
   glBufferData(GL_ARRAY_BUFFER, verts.size()*sizeof(float), &verts[0], GL_STATIC_DRAW);

   glGenBuffers(1, &norBufID);
   glBindBuffer(GL_ARRAY_BUFFER, norBufID);
   glBufferData(GL_ARRAY_BUFFER, norms.size()*sizeof(float), &norms[0], GL_STATIC_DRAW);

   glGenBuffers(1, &indBufID);
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indBufID);
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, inds.size()*sizeof(int),&inds[0],GL_STATIC_DRAW);

   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
   assert(glGetError() == GL_NO_ERROR);
}

void SimpleShape::draw(GLint h_pos, GLint h_nor)
{
   // Enable and bind position array for drawing
   GLSL::enableVertexAttribArray(h_pos);
   glBindBuffer(GL_ARRAY_BUFFER, posBufID);
   glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, 0);

   // Enable and bind normal array (if it exists) for drawing
   if(norBufID) {
      GLSL::enableVertexAttribArray(h_nor);
      glBindBuffer(GL_ARRAY_BUFFER, norBufID);
      glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, 0);
   }

   // Bind index array for drawing
   int nIndices = inds.size();
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indBufID);

   // Draw
   glDrawElements(GL_TRIANGLES, nIndices, GL_UNSIGNED_INT, 0);

   // Disable and unbind
   if(norBufID) {
      GLSL::disableVertexAttribArray(h_nor);
   }
   GLSL::disableVertexAttribArray(h_pos);
   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

