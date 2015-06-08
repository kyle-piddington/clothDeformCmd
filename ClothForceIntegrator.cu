
#include "ClothForceIntegrator.h"
#include <algorithm>
#include <iostream>

#define YOUNG_MOD 1000.0 //N/m
#define POISSON_COEFF 0.85
#define MASS 600.0 //kg/m
#define GRAVITY -4.5
#define DAMPN  0.05
#define COLLISIONSTR 20.0

#define MIN_STEP 0.0005 //Seconds

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
}

void ClothForceIntegrator::step(double stepAmnt, float * outputVertices, std::vector<int> & theLockedVerts )
{
   //Calculate a force vector
   int numLockedVerts = theLockedVerts.size();
   int *lockedVerts = new int[numLockedVerts];
   memcpy(lockedVerts, theLockedVerts.data(), numLockedVerts * sizeof(int));
  
   for(int steps = 0; steps < (int)(stepAmnt/MIN_STEP); steps++)
   {
      const double dt = MIN_STEP;
      const double recipMass = 1/(MASS/numVerts);

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
		  forceX[indicies[i*3]] += forceA.x - DAMPN * vA.x;
		  forceY[indicies[i*3]] += forceA.y - DAMPN * vA.y;
		  forceZ[indicies[i*3]] += forceA.z - DAMPN * vA.z;

		  //Update second vertex
		  forceX[indicies[i*3+1]] += forceB.x - DAMPN * vB.x;
		  forceY[indicies[i*3+1]] += forceB.y - DAMPN * vB.y;
		  forceZ[indicies[i*3+1]] += forceB.z - DAMPN * vB.z;

		  //update third vertex
		  forceX[indicies[i*3 + 2]] += forceC.x - DAMPN * vC.x;
		  forceY[indicies[i*3 + 2]] += forceC.y - DAMPN * vC.y;
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
   }
   
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
