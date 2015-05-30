#include "Volume.h"
#include <stdio.h>
#include <iostream>
Volume::Volume(glm::vec3 _origin, glm::vec3 _lhw, int _resolution):
origin(_origin),
lhw(_lhw),
resolution(_resolution)
{
   step = lhw/(float)resolution;
}

float Volume::evaluate(Coord coord,float (*fp)(glm::vec3))
{
   glm::vec3 point =  step*glm::vec3(coord.x,coord.y,coord.z)+origin;
   //std::cout << point.x << ' ' << point.y << ' ' << point.z << std::endl;
   float val =  fp(point);

   return val;
}

glm::vec3 Volume::evaluateNormal(Coord coord, float(*fp)(glm::vec3))
{

   glm::vec3 center = step*glm::vec3(coord.x,coord.y,coord.z)+origin;
   float dx = (fp(center+glm::vec3(step.x,0,0)) - fp(center-glm::vec3(step.x,0,0)))/(2*step.x);
   float dy = (fp(center+glm::vec3(0,step.y,0)) - fp(center-glm::vec3(0,step.y,0)))/(2*step.y);
   float dz = (fp(center+glm::vec3(0,0,step.z)) - fp(center-glm::vec3(0,0,step.z)))/(2*step.z);

   //std::cout << dx << ' ' << dy << ' ' << dz << ' ' << std::endl;
   return glm::normalize(glm::vec3(dx,dy,dz));
}
Volume::~Volume(){}




VolBlock Volume::evaluateBlock(Coord coord, float(*fp)(glm::vec3))
{
   VolBlock block = VolBlock();
   block.values[0] = this->evaluate(Coord(coord.x,coord.y,coord.z+1), fp);
   block.values[1] = this->evaluate(Coord(coord.x+1,coord.y,coord.z+1), fp);
   block.values[2] = this->evaluate(Coord(coord.x+1,coord.y,coord.z), fp);
   block.values[3] = this->evaluate(Coord(coord.x,coord.y,coord.z), fp);
   block.values[4] = this->evaluate(Coord(coord.x,coord.y+1,coord.z+1), fp);
   block.values[5] = this->evaluate(Coord(coord.x+1,coord.y+1,coord.z+1), fp);
   block.values[6] = this->evaluate(Coord(coord.x+1,coord.y+1,coord.z), fp);
   block.values[7] = this->evaluate(Coord(coord.x,coord.y+1,coord.z), fp);

   block.normals[0] = this->evaluateNormal(Coord(coord.x,coord.y,coord.z+1), fp);
   block.normals[1] = this->evaluateNormal(Coord(coord.x+1,coord.y,coord.z+1), fp);
   block.normals[2] = this->evaluateNormal(Coord(coord.x+1,coord.y,coord.z), fp);
   block.normals[3] = this->evaluateNormal(Coord(coord.x,coord.y,coord.z), fp);
   block.normals[4] = this->evaluateNormal(Coord(coord.x,coord.y+1,coord.z+1), fp);
   block.normals[5] = this->evaluateNormal(Coord(coord.x+1,coord.y+1,coord.z+1), fp);
   block.normals[6] = this->evaluateNormal(Coord(coord.x+1,coord.y+1,coord.z), fp);
   block.normals[7] = this->evaluateNormal(Coord(coord.x,coord.y+1,coord.z), fp);
   
   block.positions[0] = this->position(Coord(coord.x,coord.y,coord.z+1));
   block.positions[1] = this->position(Coord(coord.x+1,coord.y,coord.z+1));
   block.positions[2] = this->position(Coord(coord.x+1,coord.y,coord.z));
   block.positions[3] = this->position(Coord(coord.x,coord.y,coord.z));
   block.positions[4] = this->position(Coord(coord.x,coord.y+1,coord.z+1));
   block.positions[5] = this->position(Coord(coord.x+1,coord.y+1,coord.z+1));
   block.positions[6] = this->position(Coord(coord.x+1,coord.y+1,coord.z));
   block.positions[7] = this->position(Coord(coord.x,coord.y+1,coord.z));
   return block;
}


