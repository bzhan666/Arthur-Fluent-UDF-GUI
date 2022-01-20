#include "udf.h"
 DEFINE_INIT(my_init_func,d)
 {
 	cell_t c;
 	Thread *t;
 	real xc[ND_ND];
	thread_loop_c(t,d)
		{
			begin_c_loop_all(c,t)
 			{
 				C_CENTROID(xc,c,t);
 				if (xc[0]*xc[0]/0.04+xc[1]*xc[1]/0.0225 < 1)
 					C_T(c,t) = 400.;
 				else
 					C_T(c,t) = 300.;
 			}
 			end_c_loop_all(c,t)
 		}
 } 
