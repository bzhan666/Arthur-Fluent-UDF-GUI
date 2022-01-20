#include "udf.h"
DEFINE_SOURCE(xmom_source,c,t,dS,eqn)
{
	real C2=100;
	real source;
	real x[ND_ND];
	C_CENTROID(x,c,t);
	source=-C2*0.5*C_R(c,t)*x[1]*fabs(C_U(c,t))*C_U(c,t);
	dS[eqn]=-1*C2*C_R(c,t)*x[1]*fabs(C_U(c,t));
	return source;
}
