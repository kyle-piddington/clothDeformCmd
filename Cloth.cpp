#include "Cloth.h"


Cloth::Cloth(glm::vec2 size, glm::vec2 resolution)
{
   /*Cloth Generation*/
}

//precalculation
void Cloth::precalc(){
   //calculate for every triangle.
   //set constants for each triangle.  an individual point will have different constants for each triangle it's part of
}

void Cloth::init()
{
   /*Buffer generation*/
}


/*Cloth simulation*/
void Cloth::step(float dt)
{
   //Loop through every triangle. Put parallel loop here?
   
   //Calc U, V, U', V'

   //Calc forces

   //Update acc, vel, pos
}

void Cloth::bind()
{
   /*update vertices on GPU*/
}
