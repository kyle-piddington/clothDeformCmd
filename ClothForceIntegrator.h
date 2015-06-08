#ifndef __CLOTH_FORCE_INTEGRATOR__
#define __CLOTH_FORCE_INTEGRATOR__
#include <memory>
#include <vector>
#include <cuda_runtime.h>

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
   ClothForceIntegrator(std::vector<int>  indices, std::vector<float>  vertices, std::vector<float>  weights);
   ~ClothForceIntegrator();
   void step(double dt, float * outputVertices, std::vector<int> & lockedVerts);
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

   double * d_wUA; 
   double * d_wUB;
   double * d_wUC;
   double * d_wVA;
   double * d_wVB;
   double * d_wVC;
   double * d_dArray;
   
   double * vertsX;
   double * vertsY;
   double * vertsZ;
   
   double *d_vertsX, *d_vertsY, *d_vertsZ;
   
   double *d_velsX, *d_velsY, *d_velsZ;
   
   double * velsX;
   double * velsY;
   double * velsZ;
   
   double * forceX;
   double * forceY;
   double * forceZ;
   
   int * indicies;

   int * d_indicies;
   
   unsigned int * outIdx;
   unsigned int * d_outIdx;
   
   double * expandedForceX;
   double * expandedForceY;
   double * expandedForceZ;

   double * d_expandedForceX;
   double * d_expandedForceY;
   double * d_expandedForceZ;
   
   unsigned int * counts;
   unsigned int * locs;

   const int numTriangles;
   const int numVerts;
   const int numIndicies;

   double time;
   double dt;
};
#endif
