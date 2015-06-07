#include "Cloth.h"
#include "RK4.h"
#include "ClothForceIntegrator.h"

Cloth::Cloth(float w, float h, int res) :
posBufID(0),
norBufID(0),
indBufID(0),
texBufID(0),
res(res)
{
   //init();
   for(int j = 0 ; j < res; j++)
   {
      for(int i = 0; i < res; i++)
      {
         verts.push_back(-i*((float)w/(res-1)) + w/2.0);
         verts.push_back(0);
         verts.push_back(j*((float)h/(res-1)) - h/2.0);
         
         tex.push_back((float)i/(res-1) * w);
         tex.push_back((float)j/(res-1) * w);

      }
   }

   for(int j = 0; j < res-1; j++)
   {
      for(int i = 0; i < res-1; i++)
      {
         inds.push_back(i*res + j);
         inds.push_back(i*res + j + 1);
         inds.push_back(i*res + j + res + 1);
         inds.push_back(i*res + j + res + 1);
         inds.push_back(i*res + j + res);
         inds.push_back(i*res + j);
      }
   }

   numTriangles = inds.size()/3;
}
Cloth::~Cloth()
{

   delete integrator;
}


/**
 * Generate the initial buffers
 */
void Cloth::init()
{
   integrator = new ClothForceIntegrator();
   integrator->init(inds,verts,tex);
   glGenBuffers(1, &posBufID);
   //glGenBuffers(1, &texBufID);
   glGenBuffers(1, &indBufID);
   glGenBuffers(1, &norBufID);


}

/**
 * Bind the initial position, index, and normal buffers
 */
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

/**
 * Update and rebind the vertices, and store on the gpu
 */
void Cloth::rebindVerts()
{
   glBindBuffer(GL_ARRAY_BUFFER, posBufID);
   glBufferData(GL_ARRAY_BUFFER, verts.size()*sizeof(float), &verts[0], GL_STATIC_DRAW);
   glBindBuffer(GL_ARRAY_BUFFER, 0);
}

/**
 * Update and rebind the normals, and store on the GPU
 */
void Cloth::rebindNorms()
{
   norms.resize(verts.size());
   std::fill(norms.begin(), norms.end(), 0);
   for (std::vector<int>::iterator i = inds.begin(); i != inds.end(); i+=3)
   {
      Eigen::Vector3d a = Eigen::Vector3d(
         verts[3*(*i)],
         verts[3*(*i) + 1],
         verts[3*(*i) + 2]
         );
      Eigen::Vector3d b = Eigen::Vector3d(
         verts[3*(*(i + 1))],
         verts[3*(*(i + 1)) + 1],
         verts[3*(*(i + 1)) + 2]
         );
      Eigen::Vector3d c = Eigen::Vector3d(
         verts[3*(*(i + 2))],
         verts[3*(*(i + 2)) + 1],
         verts[3*(*(i + 2)) + 2]
         );
      Eigen::Vector3d norm = (b-a).cross(c-a);
      norms[3*(*i)]    +=(norm.x());
      norms[3*(*i) + 1]+=(norm.y());
      norms[3*(*i) + 2]+=(norm.z());

      norms[3*(*(i + 1))]+=(norm.x());
      norms[3*(*(i + 1)) + 1]+=(norm.y());
      norms[3*(*(i + 1)) + 2]+=(norm.z());

      norms[3*(*(i + 2))]+=(norm.x());
      norms[3*(*(i + 2)) + 1]+=(norm.y());
      norms[3*(*(i + 2)) + 2]+=(norm.z());
   }

   glBindBuffer(GL_ARRAY_BUFFER, norBufID);
   glBufferData(GL_ARRAY_BUFFER, norms.size()*sizeof(float), &norms[0], GL_STATIC_DRAW);
   glBindBuffer(GL_ARRAY_BUFFER,0);
   assert(glGetError() == GL_NO_ERROR);
  // std::cout << "Norms: " << norms.size() << std::endl;

}

/**
 * Draw the cloth
 * @param h_pos position handle
 * @param h_nor normal handle
 * @param h_tex texture handle
 */
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



void Cloth::step(double dt)
{
   std::vector<int> lockedVerts;
   for(int i = 0; i < res; i++)
   {
      lockedVerts.push_back(i);
      lockedVerts.push_back(res*(res-1) + i);
      //lockedVerts.push_back(res-1 - i);
   }
   integrator->step(dt,&verts[0],lockedVerts);
   rebindVerts();
   rebindNorms();
}



/**
 * Retrieve the UV value of a vertex
 * @param  vertIdx the vertex index
 * @return         the uv value, as a vec2
 */
Eigen::Vector2d Cloth::getUV(int vertIdx)
{
     
   return Eigen::Vector2d(tex[vertIdx*2],tex[vertIdx*2+1]);
}

/**
 * Retrieve the XYZ value of a vertex
 * @param  vertIdx the indx of the vertex
 * @return         the position of the vertex
 */
Eigen::Vector3d Cloth::getVert(int vertIdx)
{
   return Eigen::Vector3d(verts[3*vertIdx],
                          verts[3*vertIdx + 1],
                          verts[3*vertIdx +2]);
}

void Cloth::kickCenter()
{
   std::vector<int> pushVerts;
   for(int i = -5; i < 5; i++)
   {
      for(int j = -5; j < 5; j++)
      {
         pushVerts.push_back(res*res/2 + res/2 + i + j*res);
      }
   }
   integrator->addForce(pushVerts,Eigen::Vector3d(0,8000,0)) ; 
}

void Cloth::expand(float amnt)
{
   for(int i = 0; i < res; i++)
   {
      verts[3*i + 2] -= amnt;
      verts[(res*(res-1) + i)*3 + 2] += amnt;
   }
   integrator->rebind(verts);
 
}
