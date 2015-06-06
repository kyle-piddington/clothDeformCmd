#include "glm/glm.hpp"
#ifndef __RK4_SOLVER__
#define __RK4_SOLVER__ value
class RK4Solver
{
public:

   struct State
   {
       glm::vec3 pos;      // position
       glm::vec3 vel;      // velocity
       glm::vec3 nextForce;
       float mass;
    };
    struct Derivative
    {
       glm::vec3 dx;      // dx/dt = velocity
       glm::vec3 dv;      // dv/dt = acceleration
    };

   static void integrate(State & state,
      float t,
      float dt);
private:
   static  Derivative evaluate(const State &init,
      float t,
      float dt,
      const Derivative &d);
   inline static glm::vec3 acceleration(const State &state, float dt);

};
#endif