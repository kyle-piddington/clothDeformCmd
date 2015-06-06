#include "ClothForceIntegrator.h"


/**
 * Push the force derivative structure to the Phi
 */
#pragma offload_attribute (push,target(mic))
struct Vector
{
   double x;
   double y;
   double z;
};
struct ForceDerivative
{
   double h0;
   double h1;
   double h2;
   double h3;
};
#pragma offload_attribute(pop)


//OFFLOAD METHODS
//(THESE SHOULD BE OFFLOADABLE)

inline struct Vector sumU(Vector & a, Vector & b, Vector & C, int triangle)
{

   Vector point;
   point.x = a.x * wUA[triangle] + b.x * wUB[triangle] + c.x * wUC[triangle];
   point.y = a.y * wUA[triangle] + b.y * wUB[triangle] + c.y * wUC[triangle];
   point.z = a.z * wUA[triangle] + b.z * wUB[triangle] + c.z * wUC[triangle];
   return Vector;
}
inline struct Vector sumV(Vector & a, Vector & b, Vector & C, int triangle)
{

   Vector point;
   point.x = a.x * wVA[triangle] + b.x * wVB[triangle] + c.x * wVC[triangle];
   point.y = a.y * wVA[triangle] + b.y * wVB[triangle] + c.y * wVC[triangle];
   point.z = a.z * wVA[triangle] + b.z * wVB[triangle] + c.z * wVC[triangle];
   return Vector;
}

}
double calcXForceAtTime(double * vertX, double * velX, double t, double dt)
{
   return -1;
}

double calcYForceAtTime(double * vertX, double * velX, double t, double dt)
{
   return -1;
}
double calcZForceAtTime(double * vertX, double * velX, double t, double dt)
{
   return -1;
}

/**
 * Setup triangle weights organized by indicies
 */
inline void ClothForceIntegrator::caluclateTriangleWeights(Cloth & cloth)
{

   /**
    * Iterate through each triangle and calcluate a weight grouping
    */
   for(int i = 0; i < cloth.getNumTriangles()*3; i+=3)
   {
      glm::vec2 a = cloth.getUV(indicies[i]);
      glm::vec2 b = cloth.getUV(indicies[i+1]);
      glm::vec2 c = cloth.getUV(indicies[i+2]);

      //create and set weights
      float d = a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y);
      float recrip = 1.0 /  d;

      wUA[i/3] = (b.y - c.y) * recrip;
      wVA[i/3] = (c.x - b.x) * recrip;
      wUB[i/3] = (c.y - a.y) * recrip;
      wVB[i/3] = (a.x - c.x) * recrip;
      wUC[i/3] = (a.y - b.y) * recrip;
      wVC[i/3] = (b.x - a.y) * recrip;
   }

}


/**
 * Create the initial arrays, and organize data in SOA structure, push to the Phi.
 * @param cloth [description]
 */
void ClothForceIntegrator::init(Cloth & cloth)
{
   /**
    * Replace with __mm_malloc later
    */
   indicies = new int[cloth.getNumTriangles()*3];
   wUA = new double[cloth.getNumTriangles()];
   wUB = new double[cloth.getNumTriangles()];
   wUC = new double[cloth.getNumTriangles()];
   wVA = new double[cloth.getNumTriangles()];
   wVB = new double[cloth.getNumTriangles()];
   wVC = new double[cloth.getNumTriangles()];
   wUA = new double[cloth.getNumTriangles()];
   vertsX = new double[cloth.getNumVerts()];
   vertsY = new double[cloth.getNumVerts()];
   vertsZ = new double[cloth.getNumVerts()];
   velsX = new double[cloth.getNumVerts()];
   velsY = new double[cloth.getNumVerts()];
   velsZ = new double[cloth.getNumVerts()];


   memcpy(indicies,cloth.getInds().data(), cloth.getInds().size());
   caluclateTriangleWeights(cloth);
}
void ClothForceIntegrator::step(double dt, double * outputVertices )
{



   t += dt;
   //Write to output vertices
   
}

ClothForceIntegrator::~ClothForceIntegrator()
{
   delete [] wUA ;
   delete [] wUB ;
   delete [] wUC ;
   delete [] wVA ;
   delete [] wVB ;
   delete [] wVC ;
   delete [] vertsX ;
}