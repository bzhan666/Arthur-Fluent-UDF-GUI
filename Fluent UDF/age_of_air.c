#include "udf.h"

DEFINE_DIFFUSIVITY(age_diff, c, t, i)
{
	return C_R(c,t) * 2.88e-05 + C_MU_EFF(c,t) / 0.7;
}

DEFINE_SOURCE(age_sou, c, t, dS, eqn)
{
	dS[eqn] = 0.0;
	
	return C_R(c, t);
}
