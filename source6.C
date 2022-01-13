#include "udf.h"
DEFINE_SOURCE(x_mon,cell,thread,ds,eqn)
{
	real source,time;
	time=CURRENT_TIME;
	if(time<1)
	{
		source=1*C_R(cell,thread);
		ds[eqn]=0;
	}
	else
	{
		source=0;
		ds[eqn]=0;
	}
	return source;
}