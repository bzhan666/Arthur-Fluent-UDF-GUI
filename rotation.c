#include "udf.h"
#include "math.h"
#include "dynamesh_tools.h"


DEFINE_CG_MOTION(object, dt, vel, omega, time, dtime)

{

  /* reset velocities */
  NV_S (vel, =, 0.0);
  NV_S (omega, =, 0.0);

  
  /* compute velocity formula */
  vel[0]=0.47058*cos(2*M_PI*time);
  vel[1]=0.47058*sin(2*M_PI*time);
  omega[2]=15;
  printf("\n");
  printf("\n x_velocity = %g \n",vel[0]);
  
}