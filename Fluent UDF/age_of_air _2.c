#include "udf.h"
DEFINE_DIFFUSIVITY(age_diff, c, t, i)
{
	return C_MU_L(c, t) + C_MU_T(c, t);
}

DEFINE_SOURCE(age_sou, c, t, dS, eqn)
{
	dS[eqn] = 0.0;
	
	return C_R(c, t);
}

