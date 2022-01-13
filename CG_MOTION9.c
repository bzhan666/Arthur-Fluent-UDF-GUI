#include "udf.h"
DEFINE_CG_MOTION(rot,dt,vel,omega,time,dtime)
{
	omega[2]=40*3.14*2/60;
}