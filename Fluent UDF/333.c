#include "udf.h"
DEFINE_PROFILE(temp_wall_profile, t, i)
{
	real x[ND_ND];
	face_t f;
	begin_f_loop(f, t)
	{
		F_CENTROID(x, f, t);
		F_PROFILE(f, t, i) = -3200*x[0]*x[0]+1600*x[0]+300;
	}
	end_f_loop(f, t)
}

DEFINE_SOURCE(heat_source, c, t, dS, eqn)
{
	const real C1 = 1.0;
	const real C2 = 0.03;
	real source;
	real cst_aux = C1 / (C2-C_U(c, t));
	source = cst_aux * C_T(c, t);
	dS[eqn] = cst_aux;
	return source;
}
/*DEFINE_ADJUST(post_adjust, d)*/
DEFINE_ON_DEMAND(post_demand)
/*DEFINE_EXECUTE_AT_END(post_end)*/
{
	real T_min = REAL_MAX, T_max = 0.0;
	Thread *t;
	cell_t c;
	Domain *d = Get_Domain(1);

	if(N_UDM < 1)
	{
		Message0("\n Error: No UDM defined! Abort UDF execution.\n");
		return;
	}

	thread_loop_c(t, d)
	{
		begin_c_loop(c, t)
		{
			C_UDMI(c, t, 0) = C_T(c, t) - 273.15;
			if (C_T(c,t) < T_min)
				T_min = C_T(c, t);
			if (C_T(c, t) > T_max)
				T_max = C_T(c, t);
		}
		end_c_loop(c, t)
	}

	T_min =T_min- 273.15;
	T_max =T_max- 273.15;
	Message0(" Minimum temperature = %.1f degC\tMaximum temperature = %.1f degC\n", T_min, T_max);
}

