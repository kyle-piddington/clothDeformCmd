#include "ClothForceIntegrator.h"
#define YOUNG_MOD 5000.0 //N/m
#define POISSON_COEFF 0.65
#define MASS 1.0 //kg/m
#define GRAVITY -0.5
#define DAMPN  0.35
#define COLLISIONSTR 20.0
#include <algorithm>
#include <iostream>
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


//OFFLOAD METHODS
//(THESE SHOULD BE OFFLOADABLE)



inline struct Vector sumPoint(Vector & a, Vector & b, Vector & c,
                          double  wUA, double  wUB, double  wUC)
{

   struct Vector point;
   point.x = a.x * wUA + b.x * wUB + c.x * wUC;
   point.y = a.y * wUA + b.y * wUB + c.y * wUC;
   point.z = a.z * wUA + b.z * wUB + c.z * wUC;
   return point;
}



inline Vector calculateForce(Vector U, Vector V,
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
inline void ClothForceIntegrator::caluclateTriangleWeights(Cloth & cloth)
{

   /**
    * Iterate through each triangle and calcluate a weight grouping
    */
   for(int i = 0; i < cloth.getNumTriangles()*3; i+=3)
   {
      std::cout << i << std::endl;
      Eigen::Vector2d a = cloth.getUV(indicies[i]);
      Eigen::Vector2d b = cloth.getUV(indicies[i+1]);
      Eigen::Vector2d c = cloth.getUV(indicies[i+2]);
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

inline void caluclatePositionsAtTime(double * vertsX, double * vertsY, double * vertsZ,
                              double * velsX, double * velsY, double * velsZ,
                              double dt, int numVerts,
                              double * outX, double * outY, double * outZ    )
{
   for(int i = 0; i < numVerts; i++)
   {
      outX[i] = vertsX[i] + velsX[i] * dt;
      outY[i] = vertsY[i] + velsY[i] * dt;
      outZ[i] = vertsZ[i] + velsZ[i] * dt;
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
void ClothForceIntegrator::init(Cloth & cloth)
{
   /**
    * Replace with __mm_malloc later
    */
   numTriangles = cloth.getNumTriangles();
   indicies = new int[cloth.getNumTriangles()*3];
   wUA = new double[cloth.getNumTriangles()];
   wUB = new double[cloth.getNumTriangles()];
   wUC = new double[cloth.getNumTriangles()];
   wVA = new double[cloth.getNumTriangles()];
   wVB = new double[cloth.getNumTriangles()];
   wVC = new double[cloth.getNumTriangles()];
   dArray = new double[cloth.getNumTriangles()];
   numVerts = cloth.getNumVerts();
   vertsX = new double[cloth.getNumVerts()];
   vertsY = new double[cloth.getNumVerts()];
   vertsZ = new double[cloth.getNumVerts()];
   velsX = new double[cloth.getNumVerts()];
   velsY = new double[cloth.getNumVerts()];
   velsZ = new double[cloth.getNumVerts()];

   forceX = new double[cloth.getNumVerts()];
   forceY = new double[cloth.getNumVerts()];
   forceZ = new double[cloth.getNumVerts()];

   for(int i = 0; i < cloth.getNumVerts(); i++)
   {
      Eigen::Vector3d vert = cloth.getVert(i);
      vertsX[i] = vert.x();
      vertsY[i] = vert.y();
      vertsZ[i] = vert.z();
   }

   std::fill(velsX, velsX + numVerts, 0);
   std::fill(velsY, velsY + numVerts, 0);
   std::fill(velsZ, velsZ + numVerts, 0);

   std::fill(forceX, forceX + numVerts, 0);
   std::fill(forceY, forceY + numVerts, 0);
   std::fill(forceZ, forceZ + numVerts, 0);


   memcpy(indicies,cloth.getInds().data(), sizeof(int)*cloth.getInds().size());
   caluclateTriangleWeights(cloth);
}
void ClothForceIntegrator::step(double dt, float * outputVertices, std::vector<int> & lockedVerts )
{




   time += dt;
   //Calculate a force vector
   //
   //EXPLICIT EULER FOR NOW, JESUS THIS IS GOING TO SUCK
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

   /**
    * Cloth-self-intersection test
    */
   for(int i = 0; i < numVerts; i++)
   {



      forceY[i] += GRAVITY;

   }
   for (std::vector<int>::iterator i = lockedVerts.begin(); i != lockedVerts.end(); ++i)
   {
      forceX[*i]=0;
      forceY[*i]=0;
      forceZ[*i]=0;

   }
   const double recipMass = 1/MASS;

   for(int i = 0; i < numVerts; i++)
   {
      velsX[i]  += forceX[i] * dt * recipMass;
      vertsX[i] += velsX[i] * dt;

      velsY[i]  += forceY[i] * dt * recipMass;
      vertsY[i] += velsY[i] * dt;

      velsZ[i]  += forceZ[i] * dt * recipMass;
      vertsZ[i] += velsZ[i] * dt;


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
   std::fill(forceX, forceX + numVerts, 0);
   std::fill(forceY, forceY + numVerts, 0);
   std::fill(forceZ, forceZ + numVerts, 0);

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