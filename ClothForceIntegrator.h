#ifndef __CLOTH_FORCE_INTEGRATOR__
#define __CLOTH_FORCE_INTEGRATOR__
#include <memory>
#include <vector>
//#include "Eigen/Sparse"

typedef struct floatVec2D {
   float x, y;
} floatVec2d;

/**
 * An offloadable cloth simulation class to handle simulation steps;
 */
class Cloth;
class ClothForceIntegrator
{
public:
   ClothForceIntegrator(){};
   ~ClothForceIntegrator();
   void init(std::vector<int>  indices, std::vector<float>  vertices, std::vector<float>  weights);
   void step(double dt, float * outputVertices, std::vector<int> & lockedVerts);
   void startOffload();
   void endOffload();
private:
   inline void caluclateTriangleWeights( std::vector<float> & weights, std::vector<int> & inds);
   void bind();

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
   
   bool offload;
};
#endif
