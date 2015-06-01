#include "Cloth.h"

Cloth::Cloth(int w, int h, int res) :
posBufID(0),
norBufID(0),
indBufID(0),
texBufID(0)
{
   //init();
   for(int i = 0 ; i < res; i++)
   {
      for(int j = 0; j < res; j++)
      {
         verts.push_back(i*((float)w/(res-1)) - w/2.0);
         verts.push_back(j*((float)h/(res-1)) - h/2.0);
         verts.push_back(0);
         tex.push_back((float)i/(res-1));
         tex.push_back((float)j/(res-1));

      }
   }

   for(int i = 0; i < res-1; i++)
   {
      for(int j = 0; j < res-1; j++)
      {
         inds.push_back(i*res + j);
         inds.push_back(i*res + j + 1);
         inds.push_back(i*res + j + res + 1);
         inds.push_back(i*res + j + res + 1);
         inds.push_back(i*res + j + res);
         inds.push_back(i*res + j);
      }
   }

   for (std::vector<float>::iterator i = verts.begin(); i != verts.end(); i+=3)
   {
      std::cout << *i << " " << (*(i + 1)) << " " << (*(i + 2)) << std::endl;
   }
   std::cout << inds.size() << std::endl;

}
Cloth::~Cloth()
{


}

//precalculation
void Cloth::precalc(){
   //calculate for every triangle.
   //set constants for each triangle.  an individual point will have different constants for each triangle it's part of
}

void Cloth::init()
{
   glGenBuffers(1, &posBufID);
   //glGenBuffers(1, &texBufID);
   glGenBuffers(1, &indBufID);
   glGenBuffers(1, &norBufID);


}

void Cloth::bind()
{
   glBindBuffer(GL_ARRAY_BUFFER, posBufID);
   glBufferData(GL_ARRAY_BUFFER, verts.size()*sizeof(float), &verts[0], GL_STATIC_DRAW);


   //glBindBuffer(GL_ARRAY_BUFFER, texBufID);
   //glBufferData(GL_ARRAY_BUFFER, tex.size()*sizeof(float),&tex[0],GL_STATIC_DRAW);

   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indBufID);
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, inds.size()*sizeof(int),&inds[0],GL_STATIC_DRAW);

   assert(glGetError() == GL_NO_ERROR);
   rebindNorms();
   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
   

}

void Cloth::rebindVerts()
{
   glBindBuffer(GL_ARRAY_BUFFER, posBufID);
   glBufferData(GL_ARRAY_BUFFER, verts.size()*sizeof(float), &verts[0], GL_STATIC_DRAW);
   glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void Cloth::rebindNorms()
{
   norms.resize(verts.size());
   std::fill(norms.begin(), norms.end(), 0);
   for (std::vector<int>::iterator i = inds.begin(); i != inds.end(); i+=3)
   {
      glm::vec3 a = glm::vec3(
         verts[3*(*i)],
         verts[3*(*i) + 1],
         verts[3*(*i) + 2]
         );
      glm::vec3 b = glm::vec3(
         verts[3*(*(i + 1))],
         verts[3*(*(i + 1)) + 1],
         verts[3*(*(i + 1)) + 2]
         );
      glm::vec3 c = glm::vec3(
         verts[3*(*(i + 2))],
         verts[3*(*(i + 2)) + 1],
         verts[3*(*(i + 2)) + 2]
         );
      glm::vec3 norm = glm::cross((c-a),(b-a));
      norms[3*(*i)]    +=(norm.x);
      norms[3*(*i) + 1]+=(norm.y);
      norms[3*(*i) + 2]+=(norm.z);

      norms[3*(*(i + 1))]+=(norm.x);
      norms[3*(*(i + 1)) + 1]+=(norm.y);
      norms[3*(*(i + 1)) + 2]+=(norm.z);

      norms[3*(*(i + 2))]+=(norm.x);
      norms[3*(*(i + 2)) + 1]+=(norm.y);
      norms[3*(*(i + 2)) + 2]+=(norm.z);
   }

   glBindBuffer(GL_ARRAY_BUFFER, norBufID);
   glBufferData(GL_ARRAY_BUFFER, norms.size()*sizeof(float), &norms[0], GL_STATIC_DRAW);
   glBindBuffer(GL_ARRAY_BUFFER,0);
   assert(glGetError() == GL_NO_ERROR);
  // std::cout << "Norms: " << norms.size() << std::endl;

}
void Cloth::draw(GLint h_pos, GLint h_nor, GLint h_tex)
{
   // Enable and bind position array for drawing
   GLSL::enableVertexAttribArray(h_pos);
   glBindBuffer(GL_ARRAY_BUFFER, posBufID);
   glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, 0);
   
   // Enable and bind normal array (if it exists) for drawing
   if(norBufID && h_nor >= 0) {

      GLSL::enableVertexAttribArray(h_nor);
      glBindBuffer(GL_ARRAY_BUFFER, norBufID);
      glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, 0);
   }

   if(texBufID && h_tex >= 0){
      GLSL::enableVertexAttribArray(h_tex);
      glBindBuffer(GL_ARRAY_BUFFER,texBufID);
      glVertexAttribPointer(h_tex,2,GL_FLOAT,GL_FALSE, 0, 0);
   }
   // Bind index array for drawing
   int nIndices = inds.size();
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indBufID);

   // Draw
   glDrawElements(GL_TRIANGLES, nIndices, GL_UNSIGNED_INT, 0);

   // Disable and unbind
   if(norBufID && h_nor >= 0) {
      GLSL::disableVertexAttribArray(h_nor);
   }
   if(texBufID && h_tex >= 0)
   {
      GLSL::disableVertexAttribArray(h_tex);
   }

   GLSL::disableVertexAttribArray(h_pos);
   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}
void Cloth::step(float dt)
{
  /**
   * Physics Code Hither
   */
  rebindVerts();
  rebindNorms();
}
