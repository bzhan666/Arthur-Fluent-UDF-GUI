#include "udf.h"
DEFINE_CG_MOTION(dage, dt, vel, omega, time, dtime)
{
	vel[2]=1.5-0.5*time;
}
