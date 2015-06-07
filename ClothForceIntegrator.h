#ifndef __CLOTH_FORCE_INTEGRATOR__
#define __CLOTH_FORCE_INTEGRATOR__
#include "Cloth.h"
#include <memory>
#include <vector>
#include "Eigen/Sparse"

/**
 * An offloadable cloth simulation class to handle simulation steps;
 */

class ClothForceIntegrator
{
public:
   ClothForceIntegrator(){};
   ~ClothForceIntegrator();
   void init(Cloth & cloth);
   void step(double dt, float * outputVertices, std::vector<int> & lockedVerts);
   void addForce(std::vector<int> verts, Eigen::Vector3d force);
   void rebind(Cloth & cloth);

private:
   inline void caluclateTriangleWeights(Cloth & cloth);

   /**
    * An array of weights organized per triangle.
    */
    
   double * wUA;
   double * wUB;
   double * wUC;
   double * wVA;
   double * wVB;
   double * wVC;
   double * dArray;
   double * vertsX;
   double * vertsY;
   double * vertsZ;
   double * velsX;
   double * velsY;
   double * velsZ;
   int * indicies;

   double * forceX;
   double * forceY;
   double * forceZ;

   int numTriangles;
   int numVerts;

   double time;
   double dt;
   

};
#endif