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

DEFINE_ADJUST(post_adjust, d)
{
	real T_min = REAL_MAX, T_max = 0.0;

#if !PR_HOST
	Thread *t;
	cell_t c;
#endif

	if(N_UDM < 1)
	{
		Message0("\n Error: No UDM defined! Abort UDF execution.\n");
		return;
	}

#if !PR_HOST
	thread_loop_c(t, d)
	{
		begin_c_loop_int(c, t)
		{
			C_UDMI(c, t, 0) = C_T(c, t) - 273.15;
			if (C_T(c,t) < T_min)
				T_min = C_T(c, t);
			if (C_T(c, t) > T_max)
				T_max = C_T(c, t);
		}
		end_c_loop(c, t)
	}
#endif

	T_min = PRF_GRLOW1(T_min);
	T_max = PRF_GRHIGH1(T_max);

	node_to_host_real_2(T_min, T_max);
	T_min =T_min- 273.15;
	T_max =T_max- 273.15;

#if !RP_NODE
	Message0(" Minimum temperature = %.1f degC\tMaximum temperature = %.1f degC\n", T_min, T_max);
#endif
}


