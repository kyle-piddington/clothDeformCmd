#include "CubeMarcher.h"
#include "glm/glm.hpp"
#include <cassert>
#include <cstdint>
float lerpEdgeVals(float v1Val, float v2Val, float isoTarg)
{
   if(abs(isoTarg - v1Val)<0.00001)
   {
      return 0;
   }
   else if(abs(isoTarg - v2Val)<0.00001)
   {
      return 1;
   }
   else if(abs(v1Val-v2Val) < 0.00001)
   {
      return 0;
   }
   return (isoTarg - v1Val) / (v2Val-v1Val);

}
uint8_t evaluateBlkVerts(const VolBlock & blk , float isoLevel)
{
   uint8_t res = 0;
   for(int i = 0; i < 8; i++)
   {
      res |= (blk.values[i] < isoLevel) << i;
   }
   return res;
}

void buildTriangles(SimpleShape & buffer,VolBlock & block, float isoLevel)
{
   //Get verts
   uint8_t verts = evaluateBlkVerts(block, isoLevel);
   //Now get the set of triangles
   std::vector<glm::vec3> trisVerts;
   std::vector<glm::vec3> trisNorms;
   for(int i = 0; i < 16; i++)
   {
      int edgeIdx =  triTable[verts][i];
      if(edgeIdx != -1)
      {
         std::pair<int,int> verts = vertEdgeTable[edgeIdx];
         //Get vertices
         float lerpVal =
         lerpEdgeVals(block.values[verts.first],block.values[verts.second],isoLevel);

         glm::vec3 vert = block.positions[verts.first] + (block.positions[verts.second]-block.positions[verts.first])*lerpVal;
         glm::vec3 norm = block.normals[verts.first] + (block.normals[verts.second]-block.normals[verts.first])*lerpVal;
         trisVerts.push_back(vert);
         trisNorms.push_back(-norm);

      //Lerp vals
      //Lerp norms
      //store, continue
      }

   }
   assert(trisVerts.size()%3 == 0);
   for(int i = 0; i < trisVerts.size(); i+=3)
   {

      buffer.addTriangle(Triangle(
         trisVerts[i+0],trisVerts[i+1],trisVerts[i+2],
         trisNorms[i+0],trisNorms[i+1],trisNorms[i+2]
         ));
   }

}
void CubeMarcher::GenerateTriangles(SimpleShape & buffer, Volume & volume, float isoLevel, float(*fp)(glm::vec3) )
{
   int vResol = volume.getResolution();
   for(int i = 0; i < vResol-1; i++)
   {
      for(int j = 0; j < vResol-1; j++)
      {
         for(int k = 0; k < vResol-1; k++)
         {
            VolBlock blk = volume.evaluateBlock(Coord(i,j,k),fp);
            buildTriangles(buffer,blk,isoLevel);
         }
      }
   }

}