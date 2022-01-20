#include "udf.h"
DEFINE_ZONE_MOTION(Zonemotion,omega,axis,origin,vel,time,dtime)
{
	if (time<0.5)
	{
		vel[0]=2.0*(0.5-time);
	}
	else
	{
		vel[0]=0;
	}
	return;
}
