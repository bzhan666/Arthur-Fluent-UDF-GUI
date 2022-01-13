#include "udf.h"
#define PI 3.141592654
DEFINE_PROFILE(temperature,thread,position)
{
	real r[ND_ND];
	real x;
	face_t f;
	begin_f_loop(f,thread)
	{
		F_CENTROID(r,f,thread);
		x=r[0];
		F_PROFILE(f,thread,position)=300.+100.*sin(PI*x/0.005);
	}
	end_f_loop(f,thread)
}