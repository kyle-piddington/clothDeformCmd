
#include "ClothForceIntegrator.h"
#include "utils.h"
#include <algorithm>
#include <iostream>

#define YOUNG_MOD 1000.0 //N/m
#define POISSON_COEFF 0.85
#define MASS 600.0 //kg/m
#define GRAVITY -4.5
#define DAMPN  0.05
#define COLLISIONSTR 20.0

#define MIN_STEP 0.0005 //Seconds

#define MAX_THREADS 1024
#define THREADS_PER_BLOCK MAX_THREADS

const double YoungPoissonMatrixScalar = YOUNG_MOD/(1 - POISSON_COEFF * POISSON_COEFF);

struct Vector
{
   double x;
   double y;
   double z;
};

struct Vector sumPoint(Vector & a, Vector & b, Vector & c,
						  double  wUA, double  wUB, double  wUC)
{

   struct Vector point;
   point.x = a.x * wUA + b.x * wUB + c.x * wUC;
   point.y = a.y * wUA + b.y * wUB + c.y * wUC;
   point.z = a.z * wUA + b.z * wUB + c.z * wUC;
   return point;
}

Vector calculateForce(Vector U, Vector V,
							 double rUJ, double rVJ, double d)
{
   double euu = 0.5 * (U.x * U.x + U.y * U.y + U.z * U.z -1);
   double evv = 0.5 * (V.x * V.x + V.y * V.y + V.z * V.z -1);
   //EUV is 0.5* (UTV + VTU), should be 2 UtV
   double euv =  U.x * V.x + U.y * V.y + U.z * V.z;

   Vector sigmas;
   sigmas.x = euu + POISSON_COEFF * evv ;
   sigmas.y = POISSON_COEFF*evv + euu ;
   sigmas.z = euv * (1 - POISSON_COEFF) / 2 ;

   sigmas.x *= YoungPoissonMatrixScalar;
   sigmas.y *= YoungPoissonMatrixScalar;
   sigmas.z *= YoungPoissonMatrixScalar;
   Vector force;
   force.x = -fabs(d)/2 *(
			   sigmas.x * rUJ * U.x +
			   sigmas.y * rVJ * V.x +
			   sigmas.z * (rUJ * V.x + rVJ * U.x));

   force.y = -fabs(d)/2 * (
			   sigmas.x * rUJ * U.y +
			   sigmas.y * rVJ * V.y +
			   sigmas.z * (rUJ * V.y + rVJ * U.y));

   force.z = -fabs(d)/2 * (
			   sigmas.x * rUJ * U.z +
			   sigmas.y * rVJ * V.z +
			   sigmas.z * (rUJ * V.z + rVJ * U.z)) ;
   return force;
}

/**
 * Setup triangle weights organized by indicies
 */
inline void ClothForceIntegrator::caluclateTriangleWeights(std::vector<float> & weights,std::vector<int> & inds)
{
   /**
	* Iterate through each triangle and calcluate a weight grouping
	*/
   for(int i = 0; i < numTriangles*3; i+=3)
   {
     floatVec2D a, b, c;
     a.x = weights[2*inds[i]];
     a.y = weights[2*inds[i]+1];
     b.x = weights[2*inds[i+1]];
     b.y = weights[2*inds[i+1]+1];
     c.x = weights[2*inds[i+2]];
     c.y = weights[2*inds[i+2]+1];
	  //create and set weights
	  float d = a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y);
	  float recrip = 1.0 /  d;
	  dArray[i/3] = d;
	  wUA[i/3] = (b.y - c.y) * recrip;
	  wVA[i/3] = (c.x - b.x) * recrip;
	  wUB[i/3] = (c.y - a.y) * recrip;
	  wVB[i/3] = (a.x - c.x) * recrip;
	  wUC[i/3] = (a.y - b.y) * recrip;
	  wVC[i/3] = (b.x - a.x) * recrip;
   }
}

inline double lengthSq(Vector v)
{
   return v.x * v.x + v.y * v.y + v.z * v.z;
}

/**
 * Create the initial arrays, and organize data in SOA structure, push to the Phi.
 * @param cloth [description]
 */
 /**
   OrigIndicies:
   [a b c] [a b c]
   vertices:
   [x y z] [ x y z] [x y z] */
ClothForceIntegrator::ClothForceIntegrator(std::vector<int>  orig_indices, std::vector<float>  vertices, std::vector<float>  weights):
	numTriangles(orig_indices.size()/3),
	numVerts(vertices.size()/3),
	numIndicies(orig_indices.size())

{
   /**
	* Replace with __mm_malloc later
	*/
   
   std::cout << "numVerts: " << numVerts << std::endl;
   
   indicies = new int[orig_indices.size()];
   wUA = new double[numTriangles];
   wUB = new double[numTriangles];
   wUC = new double[numTriangles];
   wVA = new double[numTriangles];
   wVB = new double[numTriangles];
   wVC = new double[numTriangles];
   dArray = new double[numTriangles];

   vertsX = new double[numVerts];
   vertsY = new double[numVerts];
   vertsZ = new double[numVerts];
   velsX = new double[numVerts];
   velsY = new double[numVerts];
   velsZ = new double[numVerts];

   forceX = new double[numVerts];
   forceY = new double[numVerts];
   forceZ = new double[numVerts];
   
   expandedForceX = new double[numTriangles*3];
   expandedForceY = new double[numTriangles*3];
   expandedForceZ = new double[numTriangles*3];

   for(int i = 0; i < numVerts; i++)
   {
	  vertsX[i] = vertices[3*i];
	  vertsY[i] = vertices[3*i+1];
	  vertsZ[i] = vertices[3*i+2];
   }

   std::fill(velsX, velsX + numVerts, 0);
   std::fill(velsY, velsY + numVerts, 0);
   std::fill(velsZ, velsZ + numVerts, 0);

   std::fill(forceX, forceX + numVerts, 0);
   std::fill(forceY, forceY + numVerts, 0);
   std::fill(forceZ, forceZ + numVerts, 0);

   for(int i = 0; i < numIndicies; i++)
   {
      indicies[i] = orig_indices[i];
   }
   caluclateTriangleWeights(weights,orig_indices);
   
   
   counts = new unsigned int[numVerts];
   locs = new unsigned int[numVerts];
   outIdx = new unsigned int[numTriangles*3];
   
   for(int i = 0; i < numVerts; i++) { 	// init counts
      counts[i] = 0;
   }
   
   for(int i = 0; i < numIndicies; i++) {
      counts[indicies[i]]++;
   }
   
   std::cout << "init'd counts" << std::endl;
   
   locs[0] = 0;   	// init locs
   for(int i = 1; i < numVerts; i++) {
      locs[i] = counts[i - 1] + locs[i - 1];
   }
   
   unsigned int *tempCounts = new unsigned int[numVerts];   	// init outIdx
   for(int i = 0; i < numVerts; i++) {
      tempCounts[i] = 0;
   }
   
   for(int i = 0; i < numIndicies; i++) {
      int vertex = indicies[i];
      outIdx[i] = locs[vertex] + tempCounts[vertex];
      tempCounts[vertex]++;
   }
   
   std::cout << "dumping counts: ";
   for (int i = 0; i < numVerts; i++) {
      std::cout << counts[i] << " ";
   }
   std::cout << std::endl;
   
   std::cout << "dumping locs: ";
   for (int i = 0; i < numVerts; i++) {
      std::cout << locs[i] << " ";
   }
   std::cout << std::endl;
   
   std::cout << "dumping outIdx: ";
   for (int i = 0; i < numIndicies; i++) {
      std::cout << outIdx[i] << " ";
   }
   std::cout << std::endl;
   
    std::cout << "dumping indicies: ";
   for (int i = 0; i < numIndicies; i++) {
      std::cout << indicies[i] << " ";
   }
   std::cout << std::endl;
   
   std::cout << "ran init without dying" << std::endl;
//   exit(EXIT_FAILURE);

   const size_t vertsSize = numVerts * sizeof(double);

   checkCudaErrors(cudaMalloc(&d_vertsX, vertsSize));    // init vertsXYZ on device
   checkCudaErrors(cudaMalloc(&d_vertsY, vertsSize));
   checkCudaErrors(cudaMalloc(&d_vertsZ, vertsSize));
   
   
   
   checkCudaErrors(cudaMalloc(&d_velsX, vertsSize));  // init velsXYZ on device
   checkCudaErrors(cudaMalloc(&d_velsY, vertsSize));
   checkCudaErrors(cudaMalloc(&d_velsZ, vertsSize));

   const size_t indiciesSize = numIndicies * sizeof(int);

   // init d_indicies on device
   checkCudaErrors(cudaMalloc(&d_indicies, indiciesSize));
   checkCudaErrors(cudaMemcpy(d_indicies, indicies, indiciesSize, cudaMemcpyHostToDevice));

   // init weights on device
   const size_t weightsSize = numTriangles * sizeof(double);

   checkCudaErrors(cudaMalloc(&d_wUA, weightsSize)); 
   checkCudaErrors(cudaMalloc(&d_wUB, weightsSize));
   checkCudaErrors(cudaMalloc(&d_wUC, weightsSize));
   checkCudaErrors(cudaMalloc(&d_wVA, weightsSize));
   checkCudaErrors(cudaMalloc(&d_wVB, weightsSize));
   checkCudaErrors(cudaMalloc(&d_wVC, weightsSize));
   checkCudaErrors(cudaMalloc(&d_dArray, weightsSize));

   checkCudaErrors(cudaMemcpy(d_wUA, wUA, weightsSize, cudaMemcpyHostToDevice)); 
   checkCudaErrors(cudaMemcpy(d_wUB, wUB, weightsSize, cudaMemcpyHostToDevice));
   checkCudaErrors(cudaMemcpy(d_wUC, wUC, weightsSize, cudaMemcpyHostToDevice));
   checkCudaErrors(cudaMemcpy(d_wVA, wVA, weightsSize, cudaMemcpyHostToDevice));
   checkCudaErrors(cudaMemcpy(d_wVB, wVB, weightsSize, cudaMemcpyHostToDevice));
   checkCudaErrors(cudaMemcpy(d_wVC, wVC, weightsSize, cudaMemcpyHostToDevice));
   checkCudaErrors(cudaMemcpy(d_dArray, dArray, weightsSize, cudaMemcpyHostToDevice));

   // inid d_outIdx on device
   const size_t outIdxSize = numTriangles * 3 * sizeof(unsigned int);

   checkCudaErrors(cudaMalloc(&d_outIdx, outIdxSize));
   checkCudaErrors(cudaMemcpy(d_outIdx, outIdx, outIdxSize, cudaMemcpyHostToDevice));

   // allocate d_expandedForceXYZ
   const size_t expandedForceSize = numTriangles * 3 * sizeof(double);

   checkCudaErrors(cudaMalloc(&d_expandedForceX, expandedForceSize));
   checkCudaErrors(cudaMalloc(&d_expandedForceY, expandedForceSize));
   checkCudaErrors(cudaMalloc(&d_expandedForceZ, expandedForceSize));

   // DEBUG

   // verify init of d_indicies
   checkCudaErrors(cudaMemcpy(indicies, d_indicies, indiciesSize, cudaMemcpyDeviceToHost));

   // verify init of weights on device
   checkCudaErrors(cudaMemcpy(wUA, d_wUA, weightsSize, cudaMemcpyDeviceToHost)); 
   checkCudaErrors(cudaMemcpy(wUB, d_wUB, weightsSize, cudaMemcpyDeviceToHost));
   checkCudaErrors(cudaMemcpy(wUC, d_wUC, weightsSize, cudaMemcpyDeviceToHost));
   checkCudaErrors(cudaMemcpy(wVA, d_wVA, weightsSize, cudaMemcpyDeviceToHost));
   checkCudaErrors(cudaMemcpy(wVB, d_wVB, weightsSize, cudaMemcpyDeviceToHost));
   checkCudaErrors(cudaMemcpy(wVC, d_wVC, weightsSize, cudaMemcpyDeviceToHost));
   checkCudaErrors(cudaMemcpy(dArray, d_dArray, weightsSize, cudaMemcpyDeviceToHost));

   // verify init of d_outIdx
   checkCudaErrors(cudaMemcpy(outIdx, d_outIdx, outIdxSize, cudaMemcpyDeviceToHost));
}

#define FABS(X) ((X) < 0 ? (-(X)) : (X))

__global__
void findExpandedForce(size_t numTriangles,
                       double *vertsX, double *vertsY, double *vertsZ,
                       double *velsX, double *velsY, double *velsZ,
                       int *indicies,
                       double * wUA, double * wUB, double * wUC, double * wVA, double * wVB, double * wVC, double * dArray,
                       double * expandedForceX, double * expandedForceY, double * expandedForceZ,
                       unsigned int * outIdx) {
   size_t triNo = blockIdx.x * blockDim.x + threadIdx.x;

   if (triNo < numTriangles) {
      Vector A, B, C, vA, vB, vC;

      A.x = vertsX[indicies[triNo*3]];       // find positions and velocities of each point on this triangle
      B.x = vertsX[indicies[triNo*3 + 1]];
      C.x = vertsX[indicies[triNo*3 + 2]];

      A.y = vertsY[indicies[triNo*3]];
      B.y = vertsY[indicies[triNo*3 + 1]];
      C.y = vertsY[indicies[triNo*3 + 2]];

      A.z = vertsZ[indicies[triNo*3]];
      B.z = vertsZ[indicies[triNo*3 + 1]];
      C.z = vertsZ[indicies[triNo*3 + 2]];

      vA.x = velsX[indicies[triNo*3]];
      vB.x = velsX[indicies[triNo*3 + 1]];
      vC.x = velsX[indicies[triNo*3 + 2]];

      vA.y = velsY[indicies[triNo*3]];
      vB.y = velsY[indicies[triNo*3 + 1]];
      vC.y = velsY[indicies[triNo*3 + 2]];

      vA.z = velsZ[indicies[triNo*3]];
      vB.z = velsZ[indicies[triNo*3 + 1]];
      vC.z = velsZ[indicies[triNo*3 + 2]];


      Vector U, V;
      U.x = A.x * wUA[triNo] + B.x * wUB[triNo] + C.x * wUC[triNo];
      U.y = A.y * wUA[triNo] + B.y * wUB[triNo] + C.y * wUC[triNo];
      U.z = A.z * wUA[triNo] + B.z * wUB[triNo] + C.z * wUC[triNo];

      V.x = A.x * wVA[triNo] + B.x * wVB[triNo] + C.x * wVC[triNo];
      V.y = A.y * wVA[triNo] + B.y * wVB[triNo] + C.y * wVC[triNo];
      V.z = A.z * wVA[triNo] + B.z * wVB[triNo] + C.z * wVC[triNo];

      double d = dArray[triNo];

      Vector forceA, forceB, forceC;
      double euu, evv, euv;
      Vector sigmas;

      euu = 0.5 * (U.x * U.x + U.y * U.y + U.z * U.z -1);
      evv = 0.5 * (V.x * V.x + V.y * V.y + V.z * V.z -1);
      euv =  U.x * V.x + U.y * V.y + U.z * V.z;

      double youngPoissonMatrixScalar = YOUNG_MOD/(1 - POISSON_COEFF * POISSON_COEFF);

      sigmas.x = euu + POISSON_COEFF * evv ;
      sigmas.y = POISSON_COEFF*evv + euu ;
      sigmas.z = euv * (1 - POISSON_COEFF) / 2 ;

      sigmas.x *= youngPoissonMatrixScalar;
      sigmas.y *= youngPoissonMatrixScalar;
      sigmas.z *= youngPoissonMatrixScalar;

      if (d < 0) {
         d = -d;
      }

      forceA.x = -(d)/2 *(
               sigmas.x * wUA[triNo] * U.x +
               sigmas.y * wVA[triNo] * V.x +
               sigmas.z * (wUA[triNo] * V.x + wVA[triNo] * U.x));

      forceA.y = -(d)/2 * (
               sigmas.x * wUA[triNo] * U.y +
               sigmas.y * wVA[triNo] * V.y +
               sigmas.z * (wUA[triNo] * V.y + wVA[triNo] * U.y));

      forceA.z = -(d)/2 * (
               sigmas.x * wUA[triNo] * U.z +
               sigmas.y * wVA[triNo] * V.z +
               sigmas.z * (wUA[triNo] * V.z + wVA[triNo] * U.z));

      expandedForceX[outIdx[triNo*3]] = forceA.x - DAMPN * vA.x;
      expandedForceY[outIdx[triNo*3]] = forceA.y - DAMPN * vA.y;
      expandedForceZ[outIdx[triNo*3]] = forceA.z - DAMPN * vA.z;

      forceB.x = -(d)/2 *(
               sigmas.x * wUB[triNo] * U.x +
               sigmas.y * wVB[triNo] * V.x +
               sigmas.z * (wUB[triNo] * V.x + wVB[triNo] * U.x));

      forceB.y = -(d)/2 * (
               sigmas.x * wUB[triNo] * U.y +
               sigmas.y * wVB[triNo] * V.y +
               sigmas.z * (wUB[triNo] * V.y + wVB[triNo] * U.y));

      forceB.z = -(d)/2 * (
               sigmas.x * wUB[triNo] * U.z +
               sigmas.y * wVB[triNo] * V.z +
               sigmas.z * (wUB[triNo] * V.z + wVB[triNo] * U.z));

      forceC.x = -(d)/2 *(
               sigmas.x * wUC[triNo] * U.x +
               sigmas.y * wVC[triNo] * V.x +
               sigmas.z * (wUC[triNo] * V.x + wVC[triNo] * U.x));

      forceC.y = -(d)/2 * (
               sigmas.x * wUC[triNo] * U.y +
               sigmas.y * wVC[triNo] * V.y +
               sigmas.z * (wUC[triNo] * V.y + wVC[triNo] * U.y));

      forceC.z = -(d)/2 * (
               sigmas.x * wUC[triNo] * U.z +
               sigmas.y * wVC[triNo] * V.z +
               sigmas.z * (wUC[triNo] * V.z + wVC[triNo] * U.z));

      

      //Update second vertex
      expandedForceX[outIdx[triNo*3+1]] = forceB.x - DAMPN * vB.x;
      expandedForceY[outIdx[triNo*3+1]] = forceB.y - DAMPN * vB.y;
      expandedForceZ[outIdx[triNo*3+1]] = forceB.z - DAMPN * vB.z;

      //update third vertex
      expandedForceX[outIdx[triNo*3 + 2]] = forceC.x - DAMPN * vC.x;
      expandedForceY[outIdx[triNo*3 + 2]] = forceC.y - DAMPN * vC.y;
      expandedForceZ[outIdx[triNo*3 + 2]] = forceC.z - DAMPN * vC.z;
   }
}

void ClothForceIntegrator::step(double stepAmnt, float * outputVertices, std::vector<int> & theLockedVerts )
{
   
   const size_t vertsSize = numVerts * sizeof(double);
   
   
   
   //Calculate a force vector
   int numLockedVerts = theLockedVerts.size();
   int *lockedVerts = new int[numLockedVerts];  	  // initialize locked verts (map)
   for(int i = 0; i < theLockedVerts.size(); i++)
   {
      lockedVerts[i] = theLockedVerts[i];
   }

   for(int steps = 0; steps < (int)(stepAmnt/MIN_STEP); steps++)  	// for each step in a frame
   {
      const double dt = MIN_STEP;
      const double recipMass = 1/(MASS/numVerts);

      /* launch a thread for each triandle and write results to expanded force
      shit we will need initialized:
         d_vertsXYZ   -
         d_velsXYZ    -
         d_indicies   -
         d_(weights)  -
         d_outIdx     -
      shit we will need allocated:
         d_expandedForceX -
      */

      checkCudaErrors(cudaMemcpy(d_vertsX, vertsX, vertsSize, cudaMemcpyHostToDevice));
      checkCudaErrors(cudaMemcpy(d_vertsY, vertsY, vertsSize, cudaMemcpyHostToDevice));
      checkCudaErrors(cudaMemcpy(d_vertsZ, vertsZ, vertsSize, cudaMemcpyHostToDevice));

      checkCudaErrors(cudaMemcpy(d_velsX, velsX, vertsSize, cudaMemcpyHostToDevice));
      checkCudaErrors(cudaMemcpy(d_velsY, velsY, vertsSize, cudaMemcpyHostToDevice));
      checkCudaErrors(cudaMemcpy(d_velsZ, velsZ, vertsSize, cudaMemcpyHostToDevice));

      dim3 blocks, threads;

      blocks = dim3(numTriangles / THREADS_PER_BLOCK + 1);
      threads = dim3(THREADS_PER_BLOCK);

      findExpandedForce<<<blocks, threads>>>(numTriangles, d_vertsX, d_vertsY, d_vertsZ, d_velsX, d_velsY, d_velsZ, d_indicies, d_wUA, d_wUB, d_wUC, d_wVA, d_wVB, d_wVC, d_dArray, d_expandedForceX, d_expandedForceY, d_expandedForceZ, d_outIdx);
      cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());

	   // for(int i = 0; i < numTriangles; i++)
	   // {
		  // Vector A, B, C, vA, vB, vC;

		  // A.x = vertsX[indicies[i*3]]; 	     // find positions and velocities of each point on this triangle
		  // B.x = vertsX[indicies[i*3 + 1]];
		  // C.x = vertsX[indicies[i*3 + 2]];

		  // A.y = vertsY[indicies[i*3]];
		  // B.y = vertsY[indicies[i*3 + 1]];
		  // C.y = vertsY[indicies[i*3 + 2]];

		  // A.z = vertsZ[indicies[i*3]];
		  // B.z = vertsZ[indicies[i*3 + 1]];
		  // C.z = vertsZ[indicies[i*3 + 2]];

		  // vA.x = velsX[indicies[i*3]];
		  // vB.x = velsX[indicies[i*3 + 1]];
		  // vC.x = velsX[indicies[i*3 + 2]];

		  // vA.y = velsY[indicies[i*3]];
		  // vB.y = velsY[indicies[i*3 + 1]];
		  // vC.y = velsY[indicies[i*3 + 2]];

		  // vA.z = velsZ[indicies[i*3]];
		  // vB.z = velsZ[indicies[i*3 + 1]];
		  // vC.z = velsZ[indicies[i*3 + 2]];

		  // Vector U = sumPoint(A,B,C,wUA[i],wUB[i],wUC[i]);  	// find U,V for this triangle
		  // Vector V = sumPoint(A,B,C,wVA[i],wVB[i],wVC[i]);

		  // Vector forceA = calculateForce(U,V,wUA[i],wVA[i],dArray[i]);  	// calculate force on each point in triandle
		  // Vector forceB = calculateForce(U,V,wUB[i],wVB[i],dArray[i]);
		  // Vector forceC = calculateForce(U,V,wUC[i],wVC[i],dArray[i]);

		  // //Update first vertex
		  // expandedForceX[outIdx[i*3]] = forceA.x - DAMPN * vA.x;
		  // expandedForceY[outIdx[i*3]] = forceA.y - DAMPN * vA.y;
		  // expandedForceZ[outIdx[i*3]] = forceA.z - DAMPN * vA.z;

		  // //Update second vertex
		  // expandedForceX[outIdx[i*3+1]] = forceB.x - DAMPN * vB.x;
		  // expandedForceY[outIdx[i*3+1]] = forceB.y - DAMPN * vB.y;
		  // expandedForceZ[outIdx[i*3+1]] = forceB.z - DAMPN * vB.z;

		  // //update third vertex
		  // expandedForceX[outIdx[i*3 + 2]] = forceC.x - DAMPN * vC.x;
		  // expandedForceY[outIdx[i*3 + 2]] = forceC.y - DAMPN * vC.y;
		  // expandedForceZ[outIdx[i*3 + 2]] = forceC.z - DAMPN * vC.z;
	   // }

      const size_t expandedForceSize = numTriangles * 3 * sizeof(double); // FIXME: this should be const class member

      checkCudaErrors(cudaMemcpy(expandedForceX, d_expandedForceX, expandedForceSize, cudaMemcpyDeviceToHost));
      checkCudaErrors(cudaMemcpy(expandedForceY, d_expandedForceY, expandedForceSize, cudaMemcpyDeviceToHost));
      checkCudaErrors(cudaMemcpy(expandedForceZ, d_expandedForceZ, expandedForceSize, cudaMemcpyDeviceToHost));

      for(int i = 0; i < numVerts; i++) { // calc forceXYZ with expandedForceXYZ
         int locMax = (i+1 == numVerts)? numVerts :  locs[i+1];
         for(int j = locs[i]; j < locMax; j++  ){
            forceX[i] += expandedForceX[j];
            forceY[i] += expandedForceY[j];
            forceZ[i] += expandedForceZ[j];
         }
      }

	   for(int i = 0; i < numVerts; i++)
	   {
		  forceY[i] += GRAVITY;
	   }
      
	   for(int i = 0; i < numLockedVerts; ++i)
	   {
		  forceX[lockedVerts[i]]=0;
		  forceY[lockedVerts[i]]=0;
		  forceZ[lockedVerts[i]]=0;
	   }

	   for(int i = 0; i < numVerts; i++)
	   {
		  velsX[i]  += forceX[i] * dt * recipMass;
		  vertsX[i] += velsX[i] * dt;

		  velsY[i]  += forceY[i] * dt * recipMass;
		  vertsY[i] += velsY[i] * dt;

		  velsZ[i]  += forceZ[i] * dt * recipMass;
		  vertsZ[i] += velsZ[i] * dt;
	   }

      for(int i = 0; i < numVerts; i++) {
         forceX[i] = forceY[i] = forceZ[i] = 0;
      }
   }
   delete [] lockedVerts;
   /**
	* Set final vertex positions
	*/
   for(int i = 0; i < numVerts; i++)
   {
	  outputVertices[i*3] = (float)vertsX[i];
	  outputVertices[i*3+1] = (float)vertsY[i];
	  outputVertices[i*3+2] = (float)vertsZ[i];  
   }
}

ClothForceIntegrator::~ClothForceIntegrator()
{
   delete [] wUA ;
   delete [] wUB ;
   delete [] wUC ;
   delete [] wVA ;
   delete [] wVB ;
   delete [] wVC ;
   delete [] vertsX;
   delete [] vertsY;
   delete [] vertsZ;
   delete [] velsX;
   delete [] velsY;
   delete [] velsZ;
   delete [] indicies;
   delete [] forceX;
   delete [] forceY;
   delete [] forceZ;
   delete [] dArray;
}
