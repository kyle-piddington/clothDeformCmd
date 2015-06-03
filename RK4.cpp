
#include "RK4.h"


glm::vec3 RK4Solver::acceleration( const RK4Solver::State &state, float dt )
{
 return state.nextForce/state.mass;
}


RK4Solver::Derivative RK4Solver::evaluate( const RK4Solver::State &initial, 
   float t, 
   float dt, 
   const RK4Solver::Derivative &d )
{
 State state;
 state.pos = initial.pos + d.dx*dt;
 state.vel = initial.vel + d.dv*dt;

 RK4Solver::Derivative output;
 output.dx = state.vel;
 output.dv = acceleration(state, dt);
 return output;
}



void RK4Solver::integrate( RK4Solver::State &state, 
   float t, 
   float dt )
{
 RK4Solver::Derivative a,b,c,d;

 a = evaluate( state, t, 0.0f, Derivative() );
 b = evaluate( state, t, dt*0.5f, a );
 c = evaluate( state, t, dt*0.5f, b );
 d = evaluate( state, t, dt, c );

 glm::vec3 dxdt = 1.0f / 6.0f * 
 ( a.dx + 2.0f*(b.dx + c.dx) + d.dx );

 glm::vec3 dvdt = 1.0f / 6.0f * 
 ( a.dv + 2.0f*(b.dv + c.dv) + d.dv );

 state.pos = state.pos + dxdt * dt;
 state.vel = state.vel + dvdt * dt;
}
