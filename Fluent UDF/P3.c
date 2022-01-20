#include "udf.h"
DEFINE_PROFILE(unsteady_pressure,thread,position)
{
	face_t f;
	real t=CURRENT_TIME;
	begin_f_loop(f,thread)
	{
		F_PROFILE(f,thread,position)=101325.0+5.0*sin(10.*t);
	}
	end_f_loop(f,thread)
}
