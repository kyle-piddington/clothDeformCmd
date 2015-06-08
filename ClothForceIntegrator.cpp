
#include "ClothForceIntegrator.h"
#include <algorithm>
#include <iostream>

#define YOUNG_MOD 1000.0 //N/m
#define POISSON_COEFF 0.85
#define MASS 800.0 //kg/m
#define GRAVITY -4.5
#define DAMPN  0.05
#define COLLISIONSTR 20.0

#define ALLOC alloc_if(1)
#define FREE free_if(1)
#define RETAIN free_if(0)
#define REUSE alloc_if(0)
#define MIN_STEP 0.001 //Seconds


//#include <omp.h>

void ClothForceIntegrator::rebind(std::vector<float> vertices)
{
   for(int i = 0; i < vertices.size()/3; i++)
   {

	  vertsX[i] = vertices[3*i];
	  vertsY[i] = vertices[3*i+1];
	  vertsZ[i] = vertices[3*i+2];
   }

   #ifdef OFFLOAD
   #pragma offload_transfer if(__offload) target(mic:0)\
			in(vertsX: length(numVerts) REUSE RETAIN)\
			in(vertsY: length(numVerts) REUSE RETAIN)\
			in(vertsZ: length(numVerts) REUSE RETAIN)
   #endif
}

__attribute__((target(mic:0))) const double YoungPoissonMatrixScalar = YOUNG_MOD/(1 - POISSON_COEFF * POISSON_COEFF);


struct Vector
{
   double x;
   double y;
   double z;
};


//fMETHODS
//(THESE SHOULD BE OFFLOADABLE)



__attribute__((target(mic:0))) struct Vector sumPoint(Vector & a, Vector & b, Vector & c,
					  double  wUA, double  wUB, double  wUC)
{

   struct Vector point;
point.x = a.x * wUA + b.x * wUB + c.x * wUC;
point.y = a.y * wUA + b.y * wUB + c.y * wUC;
point.z = a.z * wUA + b.z * wUB + c.z * wUC;
   return point;
}



__attribute__((target(mic:0))) Vector calculateForce(Vector U, Vector V,
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
	  std::cout << i << std::endl;
	  Eigen::Vector2d a = Eigen::Vector2d(weights[2*inds[i]],weights[2*inds[i]+1]);
	  Eigen::Vector2d b = Eigen::Vector2d(weights[2*inds[i+1]],weights[2*inds[i+1]+1]);
	  Eigen::Vector2d c = Eigen::Vector2d(weights[2*inds[i+2]],weights[2*inds[i+2]+1]);
	  std::cout << std::endl;
	  //create and set weights
	  float d = a.x() * (b.y() - c.y()) + b.x() * (c.y() - a.y()) + c.x() * (a.y() - b.y());
	  float recrip = 1.0 /  d;
	  dArray[i/3] = d;
	  wUA[i/3] = (b.y() - c.y()) * recrip;
	  wVA[i/3] = (c.x() - b.x()) * recrip;
	  wUB[i/3] = (c.y() - a.y()) * recrip;
	  wVB[i/3] = (a.x() - c.x()) * recrip;
	  wUC[i/3] = (a.y() - b.y()) * recrip;
	  wVC[i/3] = (b.x() - a.x()) * recrip;
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
void ClothForceIntegrator::init(std::vector<int>  orig_indices, std::vector<float>  vertices, std::vector<float>  weights)
{

   /**
	* Replace with __mm_malloc later
	*/
   numTriangles = orig_indices.size()/3;
   numVerts = vertices.size()/3;
   indicies = (int *) _mm_malloc(sizeof(int)*numTriangles*3,64);
   wUA = (double *) _mm_malloc(sizeof(double)*numTriangles,64);
   wUB = (double *) _mm_malloc(sizeof(double)*numTriangles,64);
   wUC = (double *) _mm_malloc(sizeof(double)*numTriangles,64);
   wVA = (double *) _mm_malloc(sizeof(double)*numTriangles,64);
   wVB = (double *) _mm_malloc(sizeof(double)*numTriangles,64);
   wVC = (double *) _mm_malloc(sizeof(double)*numTriangles,64);
   dArray = (double *) _mm_malloc(sizeof(double)*numTriangles,64);

   vertsX = (double *) _mm_malloc(sizeof(double)*numVerts,64);
   vertsY = (double *) _mm_malloc(sizeof(double)*numVerts,64);
   vertsZ = (double *) _mm_malloc(sizeof(double)*numVerts,64);
   velsX = (double *) _mm_malloc(sizeof(double)*numVerts,64);
   velsY = (double *) _mm_malloc(sizeof(double)*numVerts,64);
   velsZ = (double *) _mm_malloc(sizeof(double)*numVerts,64);

   forceX = (double *) _mm_malloc(sizeof(double)*numVerts,64);
   forceY = (double *) _mm_malloc(sizeof(double)*numVerts,64);
   forceZ = (double *) _mm_malloc(sizeof(double)*numVerts,64);

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


   memcpy(indicies,orig_indices.data(), sizeof(int)*orig_indices.size());
   caluclateTriangleWeights(weights,orig_indices);

   std::cout << "numVerts: " << numVerts << std::endl;
   std::cout << "vertsX: " << vertsX << std::endl;
   std::cout << "vertsY: " << vertsY << std::endl;
   std::cout << "vertsZ: " << vertsZ << std::endl;
   #ifdef OFFLOAD
	#pragma offload_transfer target(mic:0)\
		in(wUA: length(numTriangles) ALLOC RETAIN)\
		in(wUB: length(numTriangles) ALLOC RETAIN)\
		in(wUC: length(numTriangles) ALLOC RETAIN)\
		in(wVA: length(numTriangles) ALLOC RETAIN)\
		in(wVB: length(numTriangles) ALLOC RETAIN)\
		in(wVC: length(numTriangles) ALLOC RETAIN)\
		in(dArray: length(numTriangles) ALLOC RETAIN)\
		in(vertsX: length(numVerts) ALLOC RETAIN)\
		in(vertsY: length(numVerts) ALLOC RETAIN)\
		in(vertsZ: length(numVerts) ALLOC RETAIN)\
		in(velsX: length(numVerts) ALLOC RETAIN)\
		in(velsY: length(numVerts) ALLOC RETAIN)\
		in(velsZ: length(numVerts) ALLOC RETAIN)\
		in(forceX: length(numVerts) ALLOC RETAIN)\
		in(forceY: length(numVerts) ALLOC RETAIN)\
		in(forceZ: length(numVerts) ALLOC RETAIN)\
		in(indicies: length(numTriangles*3) ALLOC RETAIN)
	#endif
  
   std::cout << "numVerts: " << numVerts << std::endl;
   std::cout << "vertsX: " << vertsX << std::endl;
   std::cout << "vertsY: " << vertsY << std::endl;
   std::cout << "vertsZ: " << vertsZ << std::endl;
}

void ClothForceIntegrator::step(double stepAmnt, float * outputVertices, std::vector<int> & theLockedVerts )
{
   //Calculate a force vector

   int numLockedVerts = theLockedVerts.size();
   int *lockedVerts = new int[numLockedVerts];
   memcpy(lockedVerts, theLockedVerts.data(), numLockedVerts * sizeof(int));
   double * l_wUA = wUA;
   double * l_wUB = wUB;
   double * l_wUC = wUC;
   double * l_wVA = wVA;
   double * l_wVB = wVB;
   double * l_wVC = wVC;
   double * l_dArray = dArray;
   double * l_vertsX = vertsX;
   double * l_vertsY = vertsY;
   double * l_vertsZ = vertsZ;
   double * l_velsX = velsX;
   double * l_velsY = velsY;
   double * l_velsZ = velsZ;
   double * l_forceX = forceX;
   double * l_forceY = forceY;
   double * l_forceZ = forceZ;
   int * l_indicies = indicies;

   __assume_aligned(indicies,64);
   __assume_aligned(wUA,64);
   __assume_aligned(wUB,64);
   __assume_aligned(wUC,64);
   __assume_aligned(wVA,64);
   __assume_aligned(wVB,64);
   __assume_aligned(wVC,64);
   __assume_aligned(dArray,64);

   __assume_aligned(vertsX,64);
   __assume_aligned(vertsY,64);
   __assume_aligned(vertsZ,64);
   __assume_aligned(velsX,64);
   __assume_aligned(velsY,64);
   __assume_aligned(velsZ,64);

   __assume_aligned(forceX,64);
   __assume_aligned(forceY,64);
   __assume_aligned(forceZ,64);

   #ifdef OFFLOAD
   #pragma offload if(__offload) target(mic:0) \
      in(wUA: length(0) REUSE RETAIN)\
      in(wUB: length(0) REUSE RETAIN)\
      in(wUC: length(0) REUSE RETAIN)\
      in(wVA: length(0) REUSE RETAIN)\
      in(wVB: length(0) REUSE RETAIN)\
      in(wVC: length(0) REUSE RETAIN)\
      in(dArray: length(0) REUSE RETAIN)\
      in(vertsX: length(0) REUSE RETAIN)\
      in(vertsY: length(0) REUSE RETAIN)\
      in(vertsZ: length(0) REUSE RETAIN)\
      in(velsX: length(0) REUSE RETAIN)\
      in(velsY: length(0) REUSE RETAIN)\
      in(velsZ: length(0) REUSE RETAIN)\
      in(forceX: length(0) REUSE RETAIN)\
      in(forceY: length(0) REUSE RETAIN)\
      in(forceZ: length(0) REUSE RETAIN)\
      in(indicies: length(0) REUSE RETAIN)\
      inout(lockedVerts: length(numLockedVerts))
   #endif
   {
       for(int steps = 0; steps < (int)(stepAmnt/MIN_STEP); steps++)
	   {
	      // std::cout << "offloaded some stuff" << std::endl;
	
	      const double dt = MIN_STEP;
	      const double recipMass = 1/(MASS/numVerts);

		   #pragma omp for simd
		   for(int i = 0; i < numTriangles; i++)
		   {
			  Vector A, B, C, vA, vB, vC;

			  A.x = vertsX[indicies[i*3]];
			  B.x = vertsX[indicies[i*3 + 1]];
			  C.x = vertsX[indicies[i*3 + 2]];

			  A.y = vertsY[indicies[i*3]];
			  B.y = vertsY[indicies[i*3 + 1]];
			  C.y = vertsY[indicies[i*3 + 2]];

			  A.z = vertsZ[indicies[i*3]];
			  B.z = vertsZ[indicies[i*3 + 1]];
			  C.z = vertsZ[indicies[i*3 + 2]];

			  vA.x = velsX[indicies[i*3]];
			  vB.x = velsX[indicies[i*3 + 1]];
			  vC.x = velsX[indicies[i*3 + 2]];

			  vA.y = velsY[indicies[i*3]];
			  vB.y = velsY[indicies[i*3 + 1]];
			  vC.y = velsY[indicies[i*3 + 2]];

			  vA.z = velsZ[indicies[i*3]];
			  vB.z = velsZ[indicies[i*3 + 1]];
			  vC.z = velsZ[indicies[i*3 + 2]];

			  Vector U = sumPoint(A,B,C,wUA[i],wUB[i],wUC[i]);
			  Vector V = sumPoint(A,B,C,wVA[i],wVB[i],wVC[i]);

			  Vector forceA = calculateForce(U,V,wUA[i],wVA[i],dArray[i]);
			  Vector forceB = calculateForce(U,V,wUB[i],wVB[i],dArray[i]);
			  Vector forceC = calculateForce(U,V,wUC[i],wVC[i],dArray[i]);

			  //Update first vertex
//			  #pragma omp atomic
			  forceX[indicies[i*3]] += forceA.x - DAMPN * vA.x;
//			  #pragma omp atomic
			  forceY[indicies[i*3]] += forceA.y - DAMPN * vA.y;
//			  #pragma omp atomic
			  forceZ[indicies[i*3]] += forceA.z - DAMPN * vA.z;

			  //Update second vertex
//			  #pragma omp atomic
			  forceX[indicies[i*3+1]] += forceB.x - DAMPN * vB.x;
//			  #pragma omp atomic     
			  forceY[indicies[i*3+1]] += forceB.y - DAMPN * vB.y;
//			  #pragma omp atomic
			  forceZ[indicies[i*3+1]] += forceB.z - DAMPN * vB.z;

			  //update third vertex
//			  #pragma omp atomic
			  forceX[indicies[i*3 + 2]] += forceC.x - DAMPN * vC.x;
//			  #pragma omp atomic
			  forceY[indicies[i*3 + 2]] += forceC.y - DAMPN * vC.y;
//			  #pragma omp atomic
			  forceZ[indicies[i*3 + 2]] += forceC.z - DAMPN * vC.z;
		   }


		   {
  	   		   #pragma simd
			   for(int i = 0; i < numVerts; i++)
			   {
				  forceY[i] += GRAVITY;
			   }
			   #pragma simd
			   for(int i = 0; i < numLockedVerts; ++i)
			   {
				  forceX[lockedVerts[i]]=0;
				  forceY[lockedVerts[i]]=0;
				  forceZ[lockedVerts[i]]=0;

			   }
			   #pragma simd
			   for(int i = 0; i < numVerts; i++)
			   {
				  velsX[i]  += forceX[i] * dt * recipMass;
				  vertsX[i] += velsX[i] * dt;

				  velsY[i]  += forceY[i] * dt * recipMass;
				  vertsY[i] += velsY[i] * dt;

				  velsZ[i]  += forceZ[i] * dt * recipMass;
				  vertsZ[i] += velsZ[i] * dt;
			    }
	   		   #pragma simd
			    for(int i = 0; i < numVerts; i++) {
			        forceX[i] = forceY[i] = forceZ[i] = 0;
			    }
			}
		   
		   // memset(forceX, numVerts * sizeof(double), 0);
		   // memset(forceY, numVerts * sizeof(double), 0);
		   // memset(forceZ, numVerts * sizeof(double), 0);
	   }
   }
   delete [] lockedVerts;

   wUA = l_wUA;
   wUB = l_wUB;
   wUC = l_wUC;
   wVA = l_wVA;
   wVB = l_wVB;
   wVC = l_wVC;
   dArray = l_dArray;
   vertsX = l_vertsX;
   vertsY = l_vertsY;
   vertsZ = l_vertsZ;
   velsX = l_velsX;
   velsY = l_velsY;
   velsZ = l_velsZ;
   forceX = l_forceX;
   forceY = l_forceY;
   forceZ = l_forceZ;
   indicies = l_indicies;

   #ifdef OFFLOAD
   #pragma offload_transfer if(__offload) target(mic:0)\
			out(vertsX: length(numVerts) REUSE RETAIN)\
			out(vertsY: length(numVerts) REUSE RETAIN)\
			out(vertsZ: length(numVerts) REUSE RETAIN)
   #endif
   /**
	* Set final vertex positions
	*/
   
   for(int i = 0; i < numVerts; i++)
   {
	  outputVertices[i*3] = (float)vertsX[i];
	  outputVertices[i*3+1] = (float)vertsY[i];
	  outputVertices[i*3+2] = (float)vertsZ[i];  
   }
 	
   
   //Print data for debugging.  But FP error...
   /*
   if(time < 2.0){
	  //printf("%f ", time);
	  for(int i = 0; i < numVerts; i++)
	  {
		 printf("(x:%f y:%f z:%f),", outputVertices[i*3], outputVertices[i*3+1], outputVertices[i*3+2]);
	  }
	  printf("\n\n");
   }

   */



   //Write to output vertices

   
}
void ClothForceIntegrator::addForce(std::vector<int> verts, Eigen::Vector3d force)
{
   for (std::vector<int>::iterator i = verts.begin(); i != verts.end(); ++i)
   {
	  forceX[*i] = force.x();
	  forceY[*i] = force.y();
	  forceZ[*i] = force.z();

   }

   #ifdef OFFLOAD
   #pragma offload_transfer if(__offload) target(mic:0)\
			in(forceX: length(numVerts) REUSE RETAIN)\
			in(forceY: length(numVerts) REUSE RETAIN)\
			in(forceZ: length(numVerts) REUSE RETAIN)
   #endif
}
void ClothForceIntegrator::startOffload()
{
   
   #ifdef OFFLOAD
   #pragma offload_transfer target(mic:0)\
		in(vertsX: length(numVerts) REUSE RETAIN)\
		in(vertsY: length(numVerts) REUSE RETAIN)\
		in(vertsZ: length(numVerts) REUSE RETAIN)\
		in(velsX: length(numVerts) REUSE RETAIN)\
		in(velsY: length(numVerts) REUSE RETAIN)\
		in(velsZ: length(numVerts) REUSE RETAIN)
	#endif
	#ifndef OFFLOAD
	  std::cout << "startOffload: Not compiled with Offload, please recompile" << std::endl;
	#endif
   __offload = true;
}
void ClothForceIntegrator::endOffload()
{

	#ifdef OFFLOAD
	#pragma offload_transfer target(mic:0)\
		out(vertsX: length(numVerts) REUSE RETAIN)\
		out(vertsY: length(numVerts) REUSE RETAIN)\
		out(vertsZ: length(numVerts) REUSE RETAIN)\
		out(velsX: length(numVerts) REUSE RETAIN)\
		out(velsY: length(numVerts) REUSE RETAIN)\
		out(velsZ: length(numVerts) REUSE RETAIN)
	#endif
	#ifndef OFFLOAD
	  std::cout << "endOffload: Not compiled with Offload, please recompile" << std::endl;
	#endif
   __offload = false;
}
ClothForceIntegrator::~ClothForceIntegrator()
{
   
	#ifdef OFFLOAD
	#pragma offload_transfer target(mic:0)\
		nocopy(wUA: length(numTriangles) FREE)\
		nocopy(wUB: length(numTriangles) FREE)\
		nocopy(wUC: length(numTriangles) FREE)\
		nocopy(wVA: length(numTriangles) FREE)\
		nocopy(wVB: length(numTriangles) FREE)\
		nocopy(wVC: length(numTriangles) FREE)\
		nocopy(dArray: length(numTriangles) FREE)\
		nocopy(vertsX: length(numVerts) FREE)\
		nocopy(vertsY: length(numVerts) FREE)\
		nocopy(vertsZ: length(numVerts) FREE)\
		nocopy(velsX: length(numVerts) FREE)\
		nocopy(velsY: length(numVerts) FREE)\
		nocopy(velsZ: length(numVerts) FREE)\
		nocopy(forceX: length(numVerts) FREE)\
		nocopy(forceY: length(numVerts) FREE)\
		nocopy(forceZ: length(numVerts) FREE)\
		nocopy(indicies: length(numTriangles*3) FREE)
	#endif
   _mm_free(wUA );
   _mm_free(wUB );
   _mm_free(wUC );
   _mm_free(wVA );
   _mm_free(wVB );
   _mm_free(wVC );
   _mm_free(vertsX);
   _mm_free(vertsY);
   _mm_free(vertsZ);
   _mm_free(velsX);
   _mm_free(velsY);
   _mm_free(velsZ);
   _mm_free(indicies);
   _mm_free(forceX);
   _mm_free(forceY);
   _mm_free(forceZ);
   _mm_free(dArray);
}
