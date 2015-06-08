
#include "ClothForceIntegrator.h"
#include <algorithm>
#include <iostream>

#define YOUNG_MOD 1000.0 //N/m
#define POISSON_COEFF 0.85
#define MASS 600.0 //kg/m
#define GRAVITY -4.5
#define DAMPN  0.05
#define COLLISIONSTR 20.0

#define ALLOC alloc_if(1)
#define FREE free_if(1)
#define RETAIN free_if(0)
#define REUSE alloc_if(0)
#define MIN_STEP 0.0005 //Seconds


//#include <omp.h>

 const double YoungPoissonMatrixScalar = YOUNG_MOD/(1 - POISSON_COEFF * POISSON_COEFF);

/**
 * Push the force derivative structure to the Phi
 */
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
	  std::cout << i << std::endl;
     floatVec2D a, b, c;
     a.x = weights[2*inds[i]];
     a.y = weights[2*inds[i]+1];
     b.x = weights[2*inds[i+1]];
     b.y = weights[2*inds[i+1]+1];
     c.x = weights[2*inds[i+2]];
     c.y = weights[2*inds[i+2]+1];
	  std::cout << std::endl;
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
void ClothForceIntegrator::init(std::vector<int>  orig_indices, std::vector<float>  vertices, std::vector<float>  weights)
{

   /**
	* Replace with __mm_malloc later
	*/
   numTriangles = orig_indices.size()/3;
   numVerts = vertices.size()/3;
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
	#ifdef OFFLOAD
	#pragma offload_transfer target(mic)\
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

}

void ClothForceIntegrator::step(double stepAmnt, float * outputVertices, std::vector<int> & theLockedVerts )
{
   //Calculate a force vector

   int numLockedVerts = theLockedVerts.size();
   int *lockedVerts = new int[numLockedVerts];
   memcpy(lockedVerts, theLockedVerts.data(), numLockedVerts * sizeof(int));

   std::cout << "sup, bitches" << std::endl;
   std::cout << "numVerts: " << numVerts << std::endl;
   std::cout << "outputVertices: " << outputVertices << std::endl;
   std::cout << "vertsX: " << vertsX << std::endl;
   std::cout << "vertsY: " << vertsY << std::endl;
   std::cout << "vertsZ: " << vertsZ << std::endl;
   
   #ifdef OFFLOAD
   #pragma offload target(mic)\
      nocopy(wUA: length(numTriangles) REUSE RETAIN)\
      nocopy(wUB: length(numTriangles) REUSE RETAIN)\
      nocopy(wUC: length(numTriangles) REUSE RETAIN)\
      nocopy(wVA: length(numTriangles) REUSE RETAIN)\
      nocopy(wVB: length(numTriangles) REUSE RETAIN)\
      nocopy(wVC: length(numTriangles) REUSE RETAIN)\
      nocopy(dArray: length(numTriangles) REUSE RETAIN)\
      nocopy(vertsX: length(numVerts) REUSE RETAIN)\
      nocopy(vertsY: length(numVerts) REUSE RETAIN)\
      nocopy(vertsZ: length(numVerts) REUSE RETAIN)\
      nocopy(velsX: length(numVerts) REUSE RETAIN)\
      nocopy(velsY: length(numVerts) REUSE RETAIN)\
      nocopy(velsZ: length(numVerts) REUSE RETAIN)\
      nocopy(forceX: length(numVerts) REUSE RETAIN)\
      nocopy(forceY: length(numVerts) REUSE RETAIN)\
      nocopy(forceZ: length(numVerts) REUSE RETAIN)\
      nocopy(indicies: length(numTriangles*3) REUSE RETAIN)\
      in(stepAmnt)\
      in(numTriangles)\
      in(numVerts)\
      in(lockedVerts: length(numLockedVerts))\
      in(numLockedVerts)
   #endif
  
   for(int steps = 0; steps < (int)(stepAmnt/MIN_STEP); steps++)
   {
      // std::cout << "offloaded some stuff" << std::endl;
      printf(" offloaded some stuff\n");

      const double dt = MIN_STEP;
      const double recipMass = 1/(MASS/numVerts);

	   #pragma omp parallel for
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
		  #pragma omp atomic
		  forceX[indicies[i*3]] += forceA.x - DAMPN * vA.x;
		  #pragma omp atomic
		  forceY[indicies[i*3]] += forceA.y - DAMPN * vA.y;
		  #pragma omp atomic
		  forceZ[indicies[i*3]] += forceA.z - DAMPN * vA.z;

		  //Update second vertex
		  #pragma omp atomic
		  forceX[indicies[i*3+1]] += forceB.x - DAMPN * vB.x;
		  #pragma omp atomic     
		  forceY[indicies[i*3+1]] += forceB.y - DAMPN * vB.y;
		  #pragma omp atomic
		  forceZ[indicies[i*3+1]] += forceB.z - DAMPN * vB.z;

		  //update third vertex
		  #pragma omp atomic
		  forceX[indicies[i*3 + 2]] += forceC.x - DAMPN * vC.x;
		  #pragma omp atomic
		  forceY[indicies[i*3 + 2]] += forceC.y - DAMPN * vC.y;
		  #pragma omp atomic
		  forceZ[indicies[i*3 + 2]] += forceC.z - DAMPN * vC.z;
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
	   
	   // memset(forceX, numVerts * sizeof(double), 0);
	   // memset(forceY, numVerts * sizeof(double), 0);
	   // memset(forceZ, numVerts * sizeof(double), 0);
   }

   std::cout << "didn't die" << std::endl;
   std::cout << "numVerts: " << numVerts << std::endl;
   std::cout << "outputVertices: " << outputVertices << std::endl;
   std::cout << "vertsX: " << vertsX << std::endl;
   std::cout << "vertsY: " << vertsY << std::endl;
   std::cout << "vertsZ: " << vertsZ << std::endl;
   /**
	* Set final vertex positions
	*/
  
   for(int i = 0; i < numVerts; i++)
   {
	  outputVertices[i*3] = (float)vertsX[i];
	  outputVertices[i*3+1] = (float)vertsY[i];
	  outputVertices[i*3+2] = (float)vertsZ[i];  
   }
 
   std::cout << "didn't die the sequel" << std::endl;  

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
void ClothForceIntegrator::startOffload()
{
   
   #ifdef OFFLOAD
   #pragma offload_transfer target(mic)\
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
}
void ClothForceIntegrator::endOffload()
{

	#ifdef OFFLOAD
	#pragma offload_transfer target(mic)\
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
}
ClothForceIntegrator::~ClothForceIntegrator()
{
   
	#ifdef OFFLOAD
	#pragma offload_transfer target(mic)\
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
