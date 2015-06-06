#ifndef __CLOTH_FORCE_INTEGRATOR__
#define __CLOTH_FORCE_INTEGRATOR__
#include "Cloth.h"
#include <memory>
/**
 * An offloadable cloth simulation class to handle simulation steps;
 */

class ClothForceIntegrator
{
public:
   ClothForceIntegrator(Cloth & cloth);
   ~ClothForceIntegrator();
   void init(Cloth & cloth);
   void step(double dt, double * outputVertices);

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
   double * vertsX;
   double * vertsY;
   double * vertsZ;
   double * velsX;
   double * velsY;
   double * velsZ;
   int * indicies;

   double time;
   double dt;
   

};
#endif