#ifndef __VOLUME_H__
#define __VOLUME_H__
#include <unordered_map>
#include <functional>
#include "glm/glm.hpp"
using namespace std;

struct Coord
{
   Coord(int _x,int _y,int _z):x(_x),y(_y),z(_z){}
   int x;
   int y;
   int z;
};

struct VolBlock
{
   
   VolBlock(){};
   float values[8];
   glm::vec3 normals[8];
   glm::vec3 positions[8];

};


class Volume
{
public:
   /**
    * Declare a new volume
    */
   Volume(glm::vec3 _origin, glm::vec3 _lhw, int _resolution);
   /**
    * Evaluate a functor at some coordinate.
    * This function also preloads nearby volume points.
    */
   float evaluate(Coord coord,float (*fp)(glm::vec3));
   
   glm::vec3 position(Coord coord)
      {
         glm::vec3 ret = step*glm::vec3(coord.x,coord.y,coord.z) + origin;
         //printf("%d %d %d: %f %f %f\n", coord.x, coord.y, coord.z, ret.x, ret.y, ret.z );
         return ret;
      }
   
   glm::vec3 evaluateNormal(Coord coord, float (*fp)(glm::vec3));
   /**
    * Evaluate a block for marching cubes
    */
   VolBlock evaluateBlock(Coord coord, float(*fp)(glm::vec3));
   ~Volume();

   int getResolution()
   {
      return resolution;
   }

private:
   glm::vec3 origin;
   glm::vec3 lhw;
   glm::vec3 step;
   int resolution;



};




#endif
