#include "udf.h"
#define PI 3.141592654

DEFINE_CG_MOTION(joker,dt,vel,omega,time,dtime)

{
	vel[0]=0.47058*cos(2*PI*time);
	vel[1]=0.47058*cos(2*PI*time);
	omega[2]=15;
}
