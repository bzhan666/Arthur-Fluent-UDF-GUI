#include "udf.h"
DEFINE_ZONE_MOTION(zonemove,omega,axis,origin,vel,time,dtime)
{
	vel[0]=0.05;
	vel[1]=0.06*sin(6.28*time);
	return;
}